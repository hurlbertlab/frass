## comparing frass and caterpillar occurrence (#traps with frass >0.1mg.d and #of surveys with caterpillars present)
#Savannah Carter
#4.4.2

#libraries
library(gsheet)
library(dplyr)
library(tidyr)
library(purrr)
library(rstatix)
library(tidyverse)
library(jsonlite)
library(daymetr)
library(corrplot)


# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#   reading in CC and altering it per julian week :
# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+

#Loading in CC data 
options(timeout = 300)  
api_url <- "https://api.github.com/repos/hurlbertlab/caterpillars-analysis-public/contents/data"
files <- fromJSON(api_url)
dataset_file <- files$name[grepl("fullDataset", files$name, ignore.case = TRUE)]
# pick the latest one
latest_file <- dataset_file[1]
github_raw <- "https://raw.githubusercontent.com/hurlbertlab/caterpillars-analysis-public/master/data/"
fullDataset <- read.csv(paste0(github_raw, latest_file))

#-------------------------------------------------------------------------------
#Meandensitybyweek function for caterpillar data
# Function for calculating the mode of a series of values
# --in this particular use case, if there multiple modes, we want the largest value
Mode = function(x){ 
  if (!is.numeric(x)) {
    stop("values must be numeric for mode calculation")
  }
  ta = table(x)
  tam = max(ta)
  mod = as.numeric(names(ta)[ta == tam])
  return(max(mod))
}

# Function for substituting values based on a condition using dplyr::mutate
# Modification of dplyr's mutate function that only acts on the rows meeting a condition
mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}

# Function for calculating and displaying arthropod phenology by week (density, mean biomas)
meanDensityByWeek = function(surveyData, # merged dataframe of Survey and arthropodSighting tables for a single site
                             ordersToInclude = 'All',       # which arthropod orders to calculate density for (codes)
                             
                             minLength = 0,         # minimum arthropod size to include 
                             jdRange = c(1,365),
                             outlierCount = 10000,
                             plot = FALSE,
                             plotVar = 'fracSurveys', # 'meanDensity' or 'fracSurveys' or 'meanBiomass'
                             minSurveyCoverage = 0.8, # minimum proportion of unique survey branches examined per week in order to include the week as a data point
                             allDates = TRUE,
                             new = TRUE,
                             color = 'black',
                             allCats = TRUE,
                             ...)                  

{
  
  if(length(ordersToInclude)==1 & ordersToInclude[1]=='All') {
    ordersToInclude = unique(surveyData$Group)
  }
  
  numUniqueBranches = length(unique(surveyData$PlantFK))
  
  firstFilter = surveyData %>%
    filter(julianday >= jdRange[1], julianday <= jdRange[2]) %>%
    mutate(julianweek = 7*floor(julianday/7) + 4)
  
  effortByWeek = firstFilter %>%
    group_by(julianweek) %>%
    summarize(nSurveyBranches = n_distinct(PlantFK),
              nSurveys = n_distinct(ID)) %>%
    mutate(modalBranchesSurveyed = Mode(5*ceiling(nSurveyBranches/5)),
           nSurveySets = nSurveys/modalBranchesSurveyed,
           modalSurveySets = Mode(round(nSurveySets)),
           okWeek = ifelse(nSurveySets/modalSurveySets >= minSurveyCoverage, 1, 0))
  
  if (allDates) {
    effortByWeek$okWeek = 1
  }
  
  if (!allCats) {
    secondFilter = firstFilter %>%
      filter(Hairy != 1, Tented != 1, Rolled != 1)
  } else {
    secondFilter = firstFilter
  }
  
  arthCount = secondFilter %>%
    filter(Length >= minLength, 
           Group %in% ordersToInclude) %>%
    mutate(Quantity2 = ifelse(Quantity > outlierCount, 1, Quantity)) %>% #outlier counts replaced with 1
    group_by(julianweek) %>%
    summarize(totalCount = sum(Quantity2, na.rm = TRUE),
              numSurveysGTzero = length(unique(ID[Quantity > 0])),
              totalBiomass = sum(Biomass_mg, na.rm = TRUE)) %>% 
    right_join(effortByWeek, by = 'julianweek') %>%
    filter(okWeek == 1) %>%
    #next line replaces 3 fields with 0 if the totalCount is NA
    mutate_cond(is.na(totalCount), totalCount = 0, numSurveysGTzero = 0, totalBiomass = 0) %>%
    mutate(meanDensity = totalCount/nSurveys,
           fracSurveys = 100*numSurveysGTzero/nSurveys,
           meanBiomass = totalBiomass/nSurveys) %>%
    arrange(julianweek) %>%
    data.frame()
  
  if (plot & new) {
    plot(arthCount$julianweek, arthCount[, plotVar], type = 'l', 
         col = color, las = 1, ...)
    points(arthCount$julianweek, arthCount[, plotVar], pch = 16, col = color, ...)
  } else if (plot & new==F) {
    points(arthCount$julianweek, arthCount[, plotVar], type = 'l', col = color, ...)
    points(arthCount$julianweek, arthCount[, plotVar], pch = 16, col = color, ...)
  }
  return(arthCount)
}
#-------------------------------------------------------------------------------
#site filter fulldataset for all years
cat_data_all_years <-fullDataset %>%
  filter(Name %in% c("NC Botanical Garden", "Prairie Ridge Ecostation"),
         Year %in% 2015:2026) 

#have meandensitybyweek aggregate caterpillar stuff by week 
cat_data_byweek <- cat_data_all_years %>%
  group_by(Year, Name) %>%
  group_split() %>%                 # split into a list, one dataframe per year
  map_dfr(~ {
    out <- meanDensityByWeek(
      surveyData = .x,
      ordersToInclude = "caterpillar",
      allDates = TRUE
    )
    out$Year <- unique(.x$Year)      # add Year back
    out$Name <- unique(.x$Name)
    out
  })%>%
  rename(Site=Name)%>%
  mutate(Site = case_when(
    Site == "NC Botanical Garden" ~ "117",
    Site == "Prairie Ridge Ecostation" ~ "8892356"  )) 
#with cat_data_byweek want to use fracSurvey and meanDensity columns


#filter by cutoff days and correct years for both sites
cat_data_byweek <- cat_data_byweek %>%
  filter(
    (Site == 117 & Year %in% c(2015, 2018, 2019, 2021, 2022) & julianweek %in% 142:200) |
      (Site == 8892356 & Year %in% c(2015:2019, 2021:2026) & julianweek %in% 154:198))


# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#   reading in frass and altering it per julian week :
# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+

# Function for reading in frass data from GoogleDoc
# *if aim is to backup GoogleDoc and write to disk only, then open =F and write = T
# *if aim is to use data without writing to disk, then open = T and write = F
frassData = function(open = T, write = F) {
  require(gsheet)
  url = "https://docs.google.com/spreadsheets/d/1RwXzwhHUbP0m5gKSOVhnKZbS1C_NrbdfHLglIVCzyFc/edit#gid=1479231778"
  data = gsheet2tbl(url)
  
  if (write) {
    # Write a copy
    write.csv(data, paste('data/frass_', Sys.Date(), '.csv', sep = ''),
              row.names = F)
  }
  if (open) { return (data) }
}
#----------------------------------------------------------------------------------
#Function for fixing time format and downloading corrected csv
TimeCleaning = function() {
  read_in_data <- gsheet2tbl('https://docs.google.com/spreadsheets/d/1RwXzwhHUbP0m5gKSOVhnKZbS1C_NrbdfHLglIVCzyFc/edit#gid=1479231778')
  
  remove_NAs <- read_in_data %>%
    filter(!is.na(Time.Set) & !is.na(Time.Collected))
  
  write.csv(remove_NAs %>% 
              mutate(Time.Set = ifelse(test = grepl(":", remove_NAs$Time.Set), 
                                       yes = remove_NAs$Time.Set, 
                                       no = paste(substr(remove_NAs$Time.Set, 1, nchar(remove_NAs$Time.Set)-2), ":", substr(remove_NAs$Time.Set, 3, 4), sep = "")), 
                     Time.Collected = ifelse(test = grepl(":", remove_NAs$Time.Collected), 
                                             yes = remove_NAs$Time.Collected, 
                                             no = paste(substr(remove_NAs$Time.Collected, 1, nchar(remove_NAs$Time.Collected)-2), ":", substr(remove_NAs$Time.Collected, 3, 4), sep = ""))), 
            paste('data/frass_', Sys.Date(), '.csv', sep = ''), row.names = F)
}
#----------------------------------------------------------------------------------
# Function that takes a date field (formatted as %m/%d/%Y) and a time field (hh:mm in 24h time), converts the date to julian day and adds the fractional
julianDayTime = function(date, hour_min) {
  require(lubridate)
  jday = yday(date)
  temp = sapply(strsplit(hour_min, ":"), function(x) { #(day represented by the hours and minutes)
    x = as.numeric(x)
    x[1] + x[2]/60
  })
  output = jday + temp/24
  return(output)
}
#-----------------------------------------------------------------------------------------------
# altering frassdata so that times and days are corrected
data = frassData(open = T) %>%
  filter(!is.na(Time.Set) & !is.na(Time.Collected)) %>%
  mutate(Date.Set = as.Date(Date.Set, format = "%m/%d/%Y"),
         Time.Set = as.character(Time.Set),
         Time.Collected = as.character(Time.Collected),
         Date.Collected = as.Date(Date.Collected, format = "%m/%d/%Y"),
         Year = format(Date.Collected, "%Y"),
         jday.Set = julianDayTime(Date.Set, Time.Set),
         jday.Collected = julianDayTime(Date.Collected, Time.Collected),
         frass.mg.d = Frass.mass..mg./(jday.Collected - jday.Set),
         frass.no.d = Frass.number/(jday.Collected - jday.Set),
         jday = (floor(jday.Collected) + floor(jday.Set))/2)

#-----------------------------------------------------------------------------------------------
#Filtering Data for Frass Occurence -> use occurance_frass_combined_weeks DF
#-----------------------------------------------------------------------------------------------
#filter data so only reliable rows are left then filter frass.mg.d so that only traps with total mass >4mg are left
filtered_mass <- data %>%
  filter(OK==1)%>% #only days deemed reliable left
  mutate(julianweek = 7 * floor(jday / 7) + 4)%>%
  mutate(included_in_trap_count = if_else(Frass.mass..mg. > 4,1,0))  #CHANGE THRESHOLD HERE
  

#-----------------------------------------------------------------------------------------------
#now need to count number of traps with same julian day  
occurance_frass <- filtered_mass %>%
  group_by(Site, Year, julianweek)%>%
  mutate(trap_occurance_percent = mean(included_in_trap_count == 1, na.rm=TRUE)) #mean here calculates the proportion along the length of the groupby columns


#issue that in 2015-2023 sites sampled twice a week, so still percentage works since dividing by total seen???
occurance_frass_combined_weeks <- occurance_frass %>%
  group_by(Site, Year, julianweek) %>%
  summarise(
    # average frass measurements
    trap_occurance_percent = mean(trap_occurance_percent, na.rm = TRUE),
    # keep representative values for the rest
    jday = min(jday, na.rm = TRUE),
    .groups = "drop")%>%
  mutate(Year=as.integer(Year)) %>%
  mutate(Site = case_when(
    Site == "Botanical Garden" ~ "8892356",
    Site == "Prairie Ridge" ~ "117"  )) %>%
  filter(
    (Site == 117 & Year %in% c(2015, 2018, 2019, 2021, 2022) & julianweek %in% 142:200) |
      (Site == 8892356 & Year %in% c(2015:2019, 2021:2026) & julianweek %in% 154:198))



#-----------------------------------------------------------------------------------------------
#Filtering Data for Frass Mass and cateprillar biomass -> IMPUTATION use imputation_data DF
#-----------------------------------------------------------------------------------------------
# using data to find mean frass per day for reliable frass only 
#read in proper url, change dates and label events for below meanfrass
url = "https://docs.google.com/spreadsheets/d/1RwXzwhHUbP0m5gKSOVhnKZbS1C_NrbdfHLglIVCzyFc/edit#gid=1611171427"
events = gsheet2tbl(url)
events$date = as.Date(events$date, format = "%m/%d/%Y")

#fiter out OK days and creating by trap area for mass values
meanfrass = data %>%
  filter(!is.na(Frass.mass..mg.)) %>%
  filter(OK == 1) %>% #keeps reliable frass row
  mutate(site = as.character(ifelse(Site=="Botanical Garden", 8892356, 117))) %>%
  group_by(site, Date.Collected, Year, jday) %>%
  summarize(
    mass = mean(frass.mg.d, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate( #altering the density based on size of the trap
    trap_area_cm2 = ifelse(Year <= 2018, 309.74, 197.71), #before 2018 use 309cm^2 after use 197cm^2 
    density_mg_cm2 = mass / trap_area_cm2,
    density_mg_m2 = density_mg_cm2 * 10000   # optional but recommended
  ) %>%
  left_join(events[, c('date', 'site', 'reliability')],
            by = c('Date.Collected' = 'date', 'site' = 'site')) %>%
  rename(date = Date.Collected)

#make sure frass has correct columns--------------------------------------
#having mean Frass data sorted by julian week
meanfrass <- meanfrass %>%
  mutate(julianweek = 7 * floor(jday / 7) + 4)%>% #this is using jday calculated from 'data' where it uses day and time set to alter jday
  mutate(Year = as.integer(Year)) #make sure integer
#combine and average meanfrass data that comes from same week (2015-2023 this happened), is mean frass per day
meanfrass_combinedweeks <- meanfrass %>%
  group_by(site, Year, julianweek) %>%
  summarise(
    # average frass measurements
    mass = mean(mass, na.rm = TRUE),
    density = mean(density_mg_cm2, na.rm = TRUE),
    # keep representative values for the rest
    date = min(date, na.rm = TRUE),   # or first(date)
    jday = mean(jday, na.rm = TRUE),
    reliability = first(reliability),
    
    .groups = "drop")

#Caterpillar Count data---------------------------------------------------
#NCBG site filter fulldataset for all years
NCBG <- fullDataset %>%
  filter(Name %in% c("NC Botanical Garden"),
         Year %in% 2015:2026)
#PR site filter fulldataset for all years
PR <- fullDataset %>%
  filter(Name %in% c("Prairie Ridge Ecostation"),
         Year %in% 2015:2026)
#have meandensitybyweek aggregate caterpillar stuff by week for NCBG
cats_NCBG <- NCBG %>%
  group_by(Year) %>%
  group_split() %>%                 # split into a list, one dataframe per year
  map_dfr(~ {
    out <- meanDensityByWeek(
      surveyData = .x,
      ordersToInclude = "caterpillar",
      allDates = TRUE
    )
    out$Year <- unique(.x$Year)      # add Year back
    out
  })%>%
  mutate(site=8892356)
#have meandensitybyweek aggregate caterpillar stuff by week for PR
cats_PR <- PR %>%
  group_by(Year) %>%
  group_split() %>%                 # split into a list, one dataframe per year
  map_dfr(~ {
    out <- meanDensityByWeek(
      surveyData = .x,
      ordersToInclude = "caterpillar",
      allDates = TRUE
    )
    out$Year <- unique(.x$Year)      # add Year back
    out
  })%>%
  mutate(site=117)
#all caterpillar data together
cats_all <- rbind(cats_NCBG, cats_PR)

#combine into one df
all_data <- cats_all %>%
  mutate(site = as.character(site)) %>% #make sure same type to join
  full_join(meanfrass_combinedweeks, by = c("julianweek", "Year", "site")) #shows where missing frass data is? is full join appropriate here?

#filter by cutoff days and correct years for both sites
all_data <- all_data %>%
  filter(
    (site == 117 & Year %in% c(2015, 2018, 2019, 2021, 2022) & julianweek %in% 142:200) |
      (site == 8892356 & Year %in% c(2015:2019, 2021:2026) & julianweek %in% 154:198)
  ) #ok kinda shows weeks where no frass data compared to CC but doesnt address issues of individual days where no data

#clean all data so only have columns I want and divide by trap area
all_data_clean <- all_data %>%
  select(site, Year, jday, julianweek, meanBiomass, date, mass) %>% #make sure no duplicate weeks
  group_by(site, Year, julianweek) %>%
  summarise(
    meanBiomass = mean(meanBiomass, na.rm = TRUE),
    mass        = mean(mass, na.rm = TRUE),
    .groups = "drop") %>%
  mutate(biomass_density = meanBiomass/(ifelse(Year <= 2018, 309.74, 197.71))) %>% #dividing by 209 for years 2018 and before, and 197 for years after
  mutate(frass_density = mass/ (ifelse(Year <= 2018, 309.74, 197.71)))


#bring in imputation data from last semester:
#imputation function
impute_biomass_data <- function(data, site_mass_defaults) {
  
  # --- Precompute first/last julianweek reference values across years ---
  first_jweek_refs <- data %>%
    group_by(site, Year) %>%
    slice_min(julianweek, n = 1) %>%
    ungroup() %>%
    group_by(site) %>%
    summarise(
      ref_first_meanBiomass = mean(meanBiomass, na.rm = TRUE),
      ref_first_mass        = mean(mass, na.rm = TRUE),
      .groups = "drop"
    )
  
  last_jweek_refs <- data %>%
    group_by(site, Year) %>%
    slice_max(julianweek, n = 1) %>%
    ungroup() %>%
    group_by(site) %>%
    summarise(
      ref_last_meanBiomass = mean(meanBiomass, na.rm = TRUE),
      ref_last_mass        = mean(mass, na.rm = TRUE),
      .groups = "drop"
    )
  
  data %>%
    left_join(first_jweek_refs, by = "site") %>%
    left_join(last_jweek_refs,  by = "site") %>%
    group_by(site, Year) %>%
    arrange(julianweek, .by_group = TRUE) %>%
    mutate(
      # --- save originals before imputation ---
      orig_meanBiomass = meanBiomass,
      orig_mass        = mass,
      
      # --- fill first julianweek NA with cross-year average for that site ---
      meanBiomass = if_else(
        is.na(meanBiomass) & julianweek == first(julianweek),
        ref_first_meanBiomass,
        meanBiomass
      ),
      mass = if_else(
        is.na(mass) & julianweek == first(julianweek),
        coalesce(site_mass_defaults[as.character(site)], ref_first_mass),
        mass
      ),
      
      # --- fill last julianweek NA with cross-year average for that site ---
      meanBiomass = if_else(
        is.na(meanBiomass) & julianweek == last(julianweek),
        ref_last_meanBiomass,
        meanBiomass
      ),
      mass = if_else(
        is.na(mass) & julianweek == last(julianweek),
        ref_last_mass,
        mass
      ),
      
      # --- interpolate interior NAs ---
      meanBiomass = na.approx(meanBiomass, x = julianweek, na.rm = FALSE, rule = 2, maxgap = 2),
      mass        = na.approx(mass,        x = julianweek, na.rm = FALSE, rule = 2, maxgap = 2),
      
      # --- flag imputed rows ---
      imputed_biomass = as.integer(is.na(orig_meanBiomass) & !is.na(meanBiomass)),
      imputed_mass    = as.integer(is.na(orig_mass)        & !is.na(mass)),
      
      # --- drop helper columns ---
      orig_meanBiomass      = NULL,
      orig_mass             = NULL,
      ref_first_meanBiomass = NULL,
      ref_first_mass        = NULL,
      ref_last_meanBiomass  = NULL,
      ref_last_mass         = NULL
    ) %>%
    ungroup()
}
#run function
imputation_data <- impute_biomass_data(all_data_clean, site_mass_defaults = c() )
#make site become Site so can combine later
imputation_data <- imputation_data %>%
  rename(Site = site)



# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#   combining frass mass/occurrence and caterpillar occurrence,biomass, and density into one dataframe: 
# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
all_five_variables_dataframe <- cat_data_byweek %>%
  select(Site, Year, julianweek, fracSurveys, meanDensity)%>%
  left_join(occurance_frass_combined_weeks, by=c("Site", "Year", "julianweek")) %>%
    left_join(imputation_data %>% select(Site, Year, julianweek, meanBiomass, mass), by=c("Site", "Year", "julianweek")) %>%
  rename(frass_mass = mass) #all years previous standardized (same jday range before joining)


# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#   running correlations between site/year/julianweek for all 5 variables to see which ones best correlated
# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
vars_of_interest <- c("fracSurveys", "meanBiomass", "meanDensity", "trap_occurance_percent", "frass_mass")

#nest data by Site x Year combo
nested_data_pearson <- all_five_variables_dataframe %>%
  select(Site, Year, all_of(vars_of_interest)) %>%
  group_by(Site, Year) %>%
  nest()%>%
  mutate(cor_matrix = map(data, ~ cor_mat(.x, vars = vars_of_interest, method ="pearson", use ="pairwise.complete.obs"))) #run correlation, does pearson feel right?

#nest data by Site x Year combo
nested_data_spearmans <- all_five_variables_dataframe %>%
  select(Site, Year, all_of(vars_of_interest)) %>%
  group_by(Site, Year) %>%
  nest()%>%
  mutate(cor_matrix = map(data, ~ cor_mat(.x, vars = vars_of_interest, method ="spearman", use ="pairwise.complete.obs"))) #run correlation, does pearson feel right?

#look at a group
nested_data_spearmans$cor_matrix[[1]]

#---------------------------------------------------------------------------
#visualizations- may have to do for loop
correlation_plotting <- function(data, year_choice, site_choice) {
  df <- data %>%
    filter(Year == year_choice, Site == site_choice)
  #pull the cor_matrix out of the list-column
  cor_df <- df$cor_matrix[[1]]
  
  #convert rstatix cor_mat() output (has a rowname/var column) into a real matrix
  cor_matrix <- cor_df %>%
    column_to_rownames(var = colnames(cor_df)[1]) %>%
    as.matrix()
  ## plot:
  corrplot(cor_matrix, 
           type = "upper", 
           title = paste(site_choice, year_choice),
           mar = c(0, 0, 2, 0),
           method = "shade", 
           order = "original", #original,hclust, alphabet
           tl.col = "black", 
           tl.srt = 45,
           cl.align.text="l",
           cl.offset = .5,
           addCoef.col = "black",   
           number.cex = 0.8,
           addgrid.col = "black",
           col = colorRampPalette(c("firebrick2", "white", "dodgerblue3"))(200))
  
  invisible(cor_matrix)
}
correlation_plotting(nested_data_spearmans, 2025, 8892356)  

##saving as a pdf------------------------ ^^^^^
# Years for each site
years_PR   <- c(2015, 2018, 2019, 2021, 2022)
years_NCBG <- setdiff(2015:2026, 2020)   
setwd("C:/Z_School/school/HurlbertLab/graphs")
#set up pdf
pdf(
  file = "correlation_spearman_pearson_OG_upper.pdf",
  width = 8,
  height = 8)
#layout for pdf
par(
  mfrow = c(3, 2),
  mar = c(4, 4, 3, 6),  
  oma = c(0, 0, 2, 0))
#loops over each sites
for (yr in years_NCBG) {
  try(
    correlation_plotting(
      data = nested_data_spearmans,
      year_choice = yr,
      site_choice = 8892356  
    ),
    silent = TRUE)}
for (yr in years_PR) {
  try(
    correlation_plotting(
      data = nested_data_spearmans,
      year_choice = yr,
      site_choice = 117   
    ),
    silent = TRUE)}
dev.off()

#---------------------------------------------------------------------------
#make new column where we keep only significant values from correlation matrix
nested_data_pearson <- nested_data_pearson %>%
  mutate(cor_signif = map(cor_matrix, cor_mark_significant))
#filter out significant pairs each matrix
signif_pairs <- nested_data_pearson %>%
  mutate(signif_long = map(cor_matrix, ~ cor_gather(.x) %>% filter(p < 0.05))) %>%
  select(Year, Site, signif_long) %>%
  unnest(signif_long)

signif_pairs
  
# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#     Do mean of correlation values when we stack them for every square 
# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#Create an array that have every one of these as a different layer and then apply a function so 
array_cormatrix <- abind(nested_data_spearmans$cor_matrix, along = 3)

#grab the row labels from the first slice (assumes same order across slices)
row_labels <- array_cormatrix[, "rowname", 1]

#drop the "rowname" column, keep only the numeric columns
num_cols <- setdiff(colnames(array_cormatrix), "rowname")
array_cormatrix_num <- array(
  as.numeric(array_cormatrix[, num_cols, ]),
  dim = c(nrow(array_cormatrix), length(num_cols), dim(array_cormatrix)[3]),
  dimnames = list(row_labels, num_cols, dimnames(array_cormatrix)[[3]])
)
#calculate the mean for each square:
mean_cormatrix <- apply(array_cormatrix_num, c(1,2), mean, na.rm= TRUE)

# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#     visualize all of the 5 variables on one line chart graph 
# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
all_variables_plotting <- function(data, year_choice, site_choice) {
  
  df <- data %>%
    filter(Year == year_choice, Site == site_choice)
  
  ## ---- Plot ----
  par(mar = c(5, 6, 4, 6))  # space for one right axis
  
  # Caterpillar fracSurveys
  plot(
    df$julianweek, df$fracSurveys,
    type = "l",
    col = "forestgreen",
    lty = "solid",
    lwd = 2,
    xlab = "Julian week",
    ylab = "",
    ylim = range(df$fracSurveys, na.rm = TRUE),
    main = paste(site_choice, year_choice)
  )
  
  #Caterpillar meanDensity
  par(new = TRUE)
  
  plot(
    df$julianweek, df$meanDensity,
    type = "l",
    col = "forestgreen",
    lty="dashed",
    lwd = 2,
    axes = FALSE,
    xlab = "",
    ylab = "",
    ylim = range(df$meanDensity, na.rm = TRUE)
  )
  #Caterpillar meanBiomass
  par(new = TRUE)
  
  plot(
    df$julianweek, df$meanBiomass,
    type = "l",
    col = "forestgreen",
    lty="dotted",
    lwd = 2,
    axes = FALSE,
    xlab = "",
    ylab = "",
    ylim = range(df$meanBiomass, na.rm = TRUE)
  ) 
  #Frass Occurence 
  par(new = TRUE)
  
  plot(
    df$julianweek, df$trap_occurance_percent,
    type = "l",
    col = "sienna",
    lty="solid",
    lwd = 2,
    axes = FALSE,
    xlab = "",
    ylab = "",
    ylim = range(df$trap_occurance_percent, na.rm = TRUE)
  )
  #Frass Occurence 
  par(new = TRUE)
  
  plot(
    df$julianweek, df$frass_mass,
    type = "l",
    col = "sienna",
    lty="dashed",
    lwd = 2,
    axes = FALSE,
    xlab = "",
    ylab = "",
    ylim = range(df$frass_mass, na.rm = TRUE)
  )
  
  ## ---- Legend ----
  legend(
    "topleft",
    legend = expression(
      paste("Cat Occurrence"),
      paste("Cat Density"),
      paste("Cat Biomass"),
      paste("Frass Occurrence"),
      paste("Frass Mass")
    ),
    col = c(
      "forestgreen", "forestgreen",
      "forestgreen", "sienna", "sienna"
    ),
    lwd = 2,
    lty = c(1, 2, 3, 1, 2),
    bty = "n",
    cex = 0.8
  )
  
  invisible(df)
}
all_variables_plotting(all_five_variables_dataframe, 2026, 8892356)

##saving as a pdf------------------------ ^^^^^
# Years for each site
years_PR   <- c(2015, 2018, 2019, 2021, 2022)
years_NCBG <- setdiff(2015:2026, 2020)   
setwd("C:/Z_School/school/HurlbertLab/graphs")
#set up pdf
pdf(
  file = "all_5_variables_plot.pdf",
  width = 8,
  height = 8)
#layout for pdf
par(
  mfrow = c(3, 2),
  mar = c(4, 4, 3, 6),  
  oma = c(0, 0, 2, 0))
#loops over each sites
for (yr in years_NCBG) {
  try(
    all_variables_plotting(
      data = all_five_variables_dataframe,
      year_choice = yr,
      site_choice = 8892356  
    ),
    silent = TRUE)}
for (yr in years_PR) {
  try(
    all_variables_plotting(
      data = all_five_variables_dataframe,
      year_choice = yr,
      site_choice = 117   
    ),
    silent = TRUE)}
dev.off()

#---------------------------------------------------
#figure out a way to make correlation graphs and line graphs on same document to view 

##saving as a pdf------------------------ ^^^^^
# Years for each site
years_PR   <- c(2015, 2018, 2019, 2021, 2022)
years_NCBG <- setdiff(2015:2026, 2020)   
datasets <- c(all_five_variables_dataframe, nested_data_spearmans)
num_datasets <- length(datasets)
setwd("C:/Z_School/school/HurlbertLab/graphs")
#set up pdf
pdf(
  file = "correlation_and_linecharts.pdf",
  width = 8,
  height = 8)
#layout for pdf
par(
  mfrow = c(3, 2),
  mar = c(4, 4, 3, 6),  
  oma = c(0, 0, 2, 0))
#loops over each sites
for (yr in seq_along(years_NCBG)) {
  x <- years_NCBG[yr] #loop and alternate
  dataset_index <- ((yr -1)%% num_datasets) +1 #determine what dataset to use
  current_dataet <- datasets[[dataset_index]]
  try(
    all_variables_plotting(
      data = all_five_variables_dataframe,
      year_choice = yr,
      site_choice = 8892356  
    ),
    silent = TRUE)}
for (yr in seq_along(years_PR)) {
  try(
    all_variables_plotting(
      data = all_five_variables_dataframe,
      year_choice = yr,
      site_choice = 117   
    ),
    silent = TRUE)}
dev.off()

# 3. Loop and alternate
for (i in seq_along(master_vector)) {
  x <- master_vector[i]
  
  # Determine which dataset to use (1, 2, 1, 2...)
  dataset_index <- ((i - 1) %% num_datasets) + 1
  current_dataset <- datasets[[dataset_index]]
  
  # Your logic here
  print(paste("Processing:", x, "using Dataset", dataset_index))
  # print(head(current_dataset)) 
}






  
  
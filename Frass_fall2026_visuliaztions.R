

#Load Libraries
library(gsheet)
library(dplyr)
library(tidyr)
library(purrr)
library(rstatix)
library(tidyverse)
library(jsonlite)
library(daymetr)
library(corrplot)
library(zoo)
library(abind)
library(pastecs)
library(StatMatch)
library(MASS)
library(ggplot2)
library(patchwork)

#loading CC data
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
#Filtering Data for Frass OCCURRENCE -> use occurrence_frass_combined_weeks DF
#-----------------------------------------------------------------------------------------------
#filter data so only reliable rows are left then filter frass.mg.d so that only traps with total mass >4mg are left
filtered_mass <- data %>%
  filter(OK==1)%>% #only days deemed reliable left
  mutate(julianweek = 7 * floor(jday / 7) + 4)%>%
  mutate(included_in_trap_count = if_else(Frass.mass..mg. > 4,1,0))  #CHANGE THRESHOLD HERE

#---------------------------------------------------------
#now need to count number of traps with same julian day  
occurrence_frass <- filtered_mass %>%
  group_by(Site, Year, julianweek)%>%
  mutate(trap_occurrence_percent = mean(included_in_trap_count == 1, na.rm=TRUE)) #mean here calculates the proportion along the length of the groupby columns


#issue that in 2015-2023 sites sampled twice a week, so still percentage works since dividing by total seen???
occurrence_frass_combined_weeks <- occurrence_frass %>%
  group_by(Site, Year, julianweek) %>%
  summarise(
    # average frass measurements
    trap_occurrence_percent = mean(trap_occurrence_percent, na.rm = TRUE),
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

#imputation function -> fills NA values
impute_data <- function(data, var) {
  
  # Turn the bare column name into a string for building new column names
  var_chr <- rlang::as_name(rlang::ensym(var))
  
  imputed_flag_name <- paste0("imputed_", var_chr)
  
  # --- Cross-year reference values for the first julianweek of each Site/Year ---
  first_jweek_refs <- data %>%
    group_by(Site, Year) %>%
    slice_min(julianweek, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    group_by(Site) %>%
    summarise(ref_first = mean({{ var }}, na.rm = TRUE), .groups = "drop")
  
  # --- Cross-year reference values for the last julianweek of each Site/Year ---
  last_jweek_refs <- data %>%
    group_by(Site, Year) %>%
    slice_max(julianweek, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    group_by(Site) %>%
    summarise(ref_last = mean({{ var }}, na.rm = TRUE), .groups = "drop")
  
  data %>%
    left_join(first_jweek_refs, by = "Site") %>%
    left_join(last_jweek_refs,  by = "Site") %>%
    group_by(Site, Year) %>%
    arrange(julianweek, .by_group = TRUE) %>%
    mutate(
      .orig = {{ var }},
      
      # fill first julianweek NA with cross-year site average
      {{ var }} := if_else(
        is.na({{ var }}) & julianweek == first(julianweek),
        ref_first,
        {{ var }}
      ),
      
      # fill last julianweek NA with cross-year site average
      {{ var }} := if_else(
        is.na({{ var }}) & julianweek == last(julianweek),
        ref_last,
        {{ var }}
      ),
      
      # interpolate interior NAs
      {{ var }} := na.approx({{ var }}, x = julianweek, na.rm = FALSE,
                             rule = 2, maxgap = 2),
      
      # flag imputed rows
      !!imputed_flag_name := as.integer(is.na(.orig) & !is.na({{ var }}))
    ) %>%
    dplyr::select(-.orig, -ref_first, -ref_last) %>%
    ungroup()
}

#imputation function -> fills zero values
impute_zeros <- function(data, var) {
  
  var_chr   <- rlang::as_name(rlang::ensym(var))
  flag_name <- paste0("imputed_", var_chr)
  
  # Keep the original values, then treat 0 as missing
  data <- data %>%
    mutate(
      .orig = {{ var }},
      {{ var }} := na_if({{ var }}, 0)
    )
  
  # Cross-year average of the FIRST week per Site/Year (zeros/NAs excluded)
  first_jweek_refs <- data %>%
    group_by(Site, Year) %>%
    slice_min(julianweek, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    group_by(Site) %>%
    summarise(ref_first = mean({{ var }}, na.rm = TRUE), .groups = "drop") %>%
    mutate(ref_first = if_else(is.nan(ref_first), NA_real_, ref_first))
  
  # Cross-year average of the LAST week per Site/Year (zeros/NAs excluded)
  last_jweek_refs <- data %>%
    group_by(Site, Year) %>%
    slice_max(julianweek, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    group_by(Site) %>%
    summarise(ref_last = mean({{ var }}, na.rm = TRUE), .groups = "drop") %>%
    mutate(ref_last = if_else(is.nan(ref_last), NA_real_, ref_last))
  
  data %>%
    left_join(first_jweek_refs, by = "Site") %>%
    left_join(last_jweek_refs,  by = "Site") %>%
    group_by(Site, Year) %>%
    arrange(julianweek, .by_group = TRUE) %>%
    mutate(
      .was_missing = is.na({{ var }}),
      
      # first week missing -> site average for first weeks in other years
      {{ var }} := if_else(
        is.na({{ var }}) & julianweek == first(julianweek),
        ref_first, {{ var }}
      ),
      
      # last week missing -> site average for last weeks in other years
      {{ var }} := if_else(
        is.na({{ var }}) & julianweek == last(julianweek),
        ref_last, {{ var }}
      ),
      
      # interior missing -> interpolate between values above and below
      {{ var }} := na.approx({{ var }}, x = julianweek, na.rm = FALSE,
                             rule = 2, maxgap = 2),
      
      # flag rows that were zero/NA and now have an imputed value
      !!flag_name := as.integer(.was_missing & !is.na({{ var }})),
      
      # anything that couldn't be filled goes back to its original value
      {{ var }} := coalesce({{ var }}, .orig)
    ) %>%
    ungroup() %>%
    dplyr::select(-.orig, -.was_missing, -ref_first, -ref_last)
}

#use the zero one
occurrence_frass_combined_weeks_imputation_data <-
  impute_zeros(occurrence_frass_combined_weeks, trap_occurrence_percent)



#-----------------------------------------------------------------------------------------------
#Filtering Data for Frass NUMBER of PELLETS 
#-----------------------------------------------------------------------------------------------
frass_pellet_number <- data %>%
  filter(OK==1)%>% #only days deemed reliable left
  mutate(julianweek = 7 * floor(jday / 7) + 4)%>%
  group_by(julianweek, Site, Year) %>%
  summarize(number_pellets = sum(Frass.number), .groups = "drop") %>%
  mutate(Site = dplyr::recode(Site, "Botanical Garden" = "8892356"))%>%
  mutate(Site = dplyr::recode(Site, "Prairie Ridge" = "117"))
frass_pellet_number$Year <- as.integer(frass_pellet_number$Year)

#impute using the na one
frass_pellet_number_imputation_data <-
  impute_data(frass_pellet_number, number_pellets)

#-----------------------------------------------------------------------------------------------
#Filtering Data for Frass and Cat MASS
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
  dplyr::select(site, Year, jday, julianweek, meanBiomass, date, mass) %>% #make sure no duplicate weeks
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


#-----------------------------------------------------------------------------------------------
#Combing all 6 variables into one DF
#-----------------------------------------------------------------------------------------------
#combining frass mass/occurrence and caterpillar occurrence, biomass, and density into one dataframe
all_six_variables_dataframe <- cat_data_byweek %>%
  dplyr::select(Site, Year, julianweek, fracSurveys, meanDensity)%>%
  left_join(occurrence_frass_combined_weeks, by=c("Site", "Year", "julianweek")) %>%
  left_join(imputation_data %>% dplyr::select(Site, Year, julianweek, meanBiomass, mass), by=c("Site", "Year", "julianweek")) %>%
  left_join(frass_pellet_number, dplyr::select(Site, Year, julianweek, number_pellets), by =c("Site", "Year", "julianweek")) %>%
  rename(frass_mass = mass) #all years previous standardized (same jday range before joining)

#-----------------------------------------------------------------------------------------------
 # Loading in temp data
#-----------------------------------------------------------------------------------------------
#loading in temp data (from Nosa)------------------------------------------------------
tmp_file = tempfile(fileext = ".csv")

AnomalySites = fullDataset[, c("Name","Year", "Latitude", "Longitude")] %>% 
  filter(Name %in% c("NC Botanical Garden", "Prairie Ridge Ecostation")) %>% 
  group_by(Name, Latitude, Longitude) %>% 
  summarise(n = n()) %>% 
  dplyr::select(-n) %>% as.data.frame()   

AnomalyDaymetr = AnomalySites %>% 
  rename(
    site = Name,
    lat = Latitude,
    lon = Longitude
  ) %>%
  write.csv(tmp_file, row.names = FALSE)

# pass temp CSV to function
TempAnomaly = download_daymet_batch(
  file_location = tmp_file,
  start = 2015,
  end = 2025, # this is the most recent available in daymetr
  internal = TRUE
)
# remove temporary file 
unlink(tmp_file)

TempAnomalyData_clean <- lapply(TempAnomaly, function(x) {
  x$data %>% mutate(site = x$site,
                    Latidue = x$latitude,
                    Longitude = x$longitude)})

AllTemp= bind_rows(TempAnomalyData_clean) #this is the file we want (has all NCBG and PR data)
#clean up globals
rm(AnomalySites, TempAnomaly, files, TempAnomalyData_clean)



# ---------------------------------------------------------------------------------------------
#   altering temperature to be what sites I want and have jday and jweek columns
# -----------------------------------------------------------------------------------------------
#rename site names so match other dataframes, also rename so can be joined properly with meanfrass data
AllTemp <- AllTemp %>%
  mutate(site = case_when(
    site == "NC Botanical Garden" ~ "117",
    site == "Prairie Ridge Ecostation" ~ "8892356"  )) 
AllTemp <- rename(AllTemp, jday=yday)
AllTemp<- rename(AllTemp, Year=year)
AllTemp <- rename (AllTemp, Site=site)


#adding average weekly temp column to alltemp data, and averaging temp for all days with same jweek value
AllTemp_weekly <- AllTemp %>%
  mutate(julianweek = 7 * floor(jday / 7) + 4) %>%
  mutate(avgtemp = (tmax..deg.c. + tmin..deg.c.) / 2)%>%
  group_by(Year, Site, julianweek) %>%
  summarise(weeklytemp = mean(avgtemp, na.rm = TRUE),
            .groups = "drop")

# ---------------------------------------------------------------------------------------------
#   combining weekly temp and 6 variables dataframe
# -----------------------------------------------------------------------------------------------
visualization_data <- all_six_variables_dataframe %>%
  left_join(AllTemp_weekly, by=c("Year", "Site", "julianweek")) %>% #2026 data not available
  group_by(Site, julianweek) %>%
  mutate(fracSurveys_avg = mean(fracSurveys, na.rm = TRUE),
         meanDensity_avg = mean(meanDensity, na.rm = TRUE),
         meanBiomass_avg = mean(meanBiomass, na.rm = TRUE),
         trap_occurrence_percent_avg = mean(trap_occurrence_percent, na.rm = TRUE),
         frass_mass_avg = mean(frass_mass, na.rm = TRUE),
         number_pellets_avg = mean(number_pellets, na.rm = TRUE),
         weeklytemp_avg = mean(weeklytemp, na.rm = TRUE)) #maybe change to summarize

#also create average columns for all 6 variables and temp
visualization_data_average <- visualization_data %>%
  group_by(Site, julianweek) %>% #should I be grouping by site??????
  summarise(fracSurveys_avg = mean(fracSurveys, na.rm = TRUE),
         meanDensity_avg = mean(meanDensity, na.rm = TRUE),
         meanBiomass_avg = mean(meanBiomass, na.rm = TRUE),
         trap_occurrence_percent_avg = mean(trap_occurrence_percent, na.rm = TRUE),
         frass_mass_avg = mean(frass_mass, na.rm = TRUE),
         number_pellets_avg = mean(number_pellets, na.rm = TRUE),
         weeklytemp_avg = mean(weeklytemp, na.rm = TRUE)) #maybe change to summarize

# ---------------------------------------------------------------------------------------------
#   creating visualization for basic anomaly plots
# -----------------------------------------------------------------------------------------------
#make a function for plotting each variable for all years and then also plotting that average as a bold line
anomaly_plotting <- function(data, data2, var_choice, avg_var_choice,
                             site_choice,
                             x_var = julianweek, year_var = Year,
                             site_var = Site) {
  
  # keep only the chosen site in both data frames
  df1 <- data  %>% filter({{ site_var }} == .env$site_choice)
  df2 <- data2 %>% filter({{ site_var }} == .env$site_choice)
  
  ggplot(df1) +
    # one thin line per year
    geom_line(aes(x = {{ x_var }},
                  y = {{ var_choice }},
                  group = {{ year_var }},
                  color = factor({{ year_var }})),
              alpha = 0.7, linewidth = 0.5) +
    # bold average line
    geom_line(data = df2,
              aes(x = {{ x_var }},
                  y = {{ avg_var_choice }}),
              color = "black", linewidth = 1.2) +
    labs(x = rlang::as_label(enquo(x_var)),
         y = rlang::as_label(enquo(var_choice)),
         color = "Year",
         title = paste("Site", site_choice, "-", rlang::as_label(enquo(var_choice)))) +
    theme_classic()
}

#use it
anomaly_plotting(visualization_data, visualization_data_average,
                 trap_occurrence_percent, trap_occurrence_percent_avg, 8892356)

#saving as a pdf------------------------ ^^^^^
vars     <- c("fracSurveys", "meanDensity", "trap_occurrence_percent",
              "meanBiomass", "number_pellets", "frass_mass", "weeklytemp")
avg_vars <- paste0(vars, "_avg")
sites    <- c(8892356, 117)

setwd("C:/Z_School/school/HurlbertLab/graphs")
pdf("6vars_linecharts.pdf", width = 8, height = 10)

for (s in sites) {
  plots <- map2(vars, avg_vars,
                ~ anomaly_plotting(visualization_data,
                                   visualization_data_average,
                                   !!rlang::sym(.x),
                                   !!rlang::sym(.y),
                                   site_choice = s))
  
  print(wrap_plots(plots, ncol = 2, nrow = 4, guides = "collect") +
          plot_annotation(title = paste("Site", s)))
}

dev.off()


# ---------------------------------------------------------------------------------------------
#   adding anomaly columns, also add two week temperature bins
# -----------------------------------------------------------------------------------------------
visualization_data_anomaly <- visualization_data%>%
  group_by(Site, julianweek) %>%
  mutate(fracSurveys_anomaly = fracSurveys - fracSurveys_avg,
         meanDensity_anomaly = meanDensity - meanDensity_avg,
         meanBiomass_anomaly = meanBiomass - meanBiomass_avg,
         trap_occurrence_percent_anomaly = trap_occurrence_percent - trap_occurrence_percent_avg,
         frass_mass_anomaly = frass_mass - frass_mass_avg,
         number_pellets_anomaly = number_pellets - number_pellets_avg,
         weeklytemp_anomaly = weeklytemp - weeklytemp_avg)


#code from frass_early_spring_investigations document
#have variables in df: all_six_variables_dataframe

x <- data.frame (all_six_variables_dataframe$Year, all_six_variables_dataframe$Site, all_six_variables_dataframe$julianweek)

#create 2 week avg temp blocks (start with alltemp_weekly)

# caterpillar key: one row per site/year/julianweek
variable_weeks <- all_six_variables_dataframe %>%
  distinct(Site, Year, julianweek) %>%
  mutate(win_start = julianweek - 14,   # two weeks before
         win_end   = julianweek - 1)    # up to, but not including, the caterpillar week

AllTemp_weekly_2week <- variable_weeks %>%
  left_join(
    AllTemp_weekly %>% dplyr::select(Site, Year, temp_week = julianweek, weeklytemp),
    by = join_by(Site, Year, win_start <= temp_week, win_end >= temp_week)
  ) %>%
  group_by(Site, Year, julianweek) %>%      # julianweek comes from the caterpillar data
  summarise(spring_avg_2week = mean(weeklytemp, na.rm = TRUE),
            n_temp_weeks = sum(!is.na(weeklytemp)),
            .groups = "drop")

#adding two week temp values
visualization_data_anomaly <- visualization_data_anomaly %>%
  left_join(AllTemp_weekly_2week, by = c("Site", "Year", "julianweek"))%>%
  group_by(Site, julianweek) %>%
  mutate(spring_avg_2week_avg = mean(spring_avg_2week, na.rm = TRUE),
         spring_avg_2week_anomaly = spring_avg_2week - spring_avg_2week_avg)


# ---------------------------------------------------------------------------------------------
#   re-running mixed effects models & interaction model
# -----------------------------------------------------------------------------------------------
#MIXED EFFECTS MODEL
library(lme4)
#all variables
allvars <- c("fracSurveys_anomaly", "meanDensity_anomaly", "trap_occurrence_percent_anomaly", "meanBiomass_anomaly", "frass_mass_anomaly", "number_pellets_anomaly")
library(lmerTest)
#lmer(response(all6)_anomaly ~ spring_temp_anomaly + site + (1 | year), data=df)
models <- lapply(allvars, function(v) {
  lmer(as.formula(paste(v, "~ spring_avg_2week_anomaly + Site + (1| Year)")), data = visualization_data_anomaly)
})
names(models) <- allvars

lapply(models, summary)

#---------------------------
#INTERACTION TERMS
#lm(response_var_anomaly ~ spring_avg_2weeks_anomaly + site + spring_avg_2weeks_anomaly * site):
models <- lapply(allvars, function(v) {
  lm(as.formula(paste(v, "~ spring_avg_2week_anomaly + Site + Site * spring_avg_2week_anomaly")), data = visualization_data_anomaly)
})
names(models) <- allvars

#to see all models at once
lapply(models, summary)

#-----------------------------------------
#lm(response_var_anomaly ~ spring_avg_2weeks_anomaly + site + Year + Site*year):
models <- lapply(allvars, function(v) {
  lm(as.formula(paste(v, "~ spring_avg_2week_anomaly + Site + Year + Site*Year")), data = visualization_data_anomaly)
})
names(models) <- allvars

#to see all models at once
lapply(models, summary) #does the spring 2-week temperature anomaly predict the response, after accounting for differences among sites and years?


# ---------------------------------------------------------------------------------------------
#   creating visualization for centered ANOMALY plots & z score ones
# -----------------------------------------------------------------------------------------------
anomaly_centered <- all_six_variables_dataframe %>%
  left_join(AllTemp_weekly, by=c("Year", "Site", "julianweek")) %>%
  group_by(Site, Year) %>%
  mutate(fracSurveys_avg = mean(fracSurveys, na.rm = TRUE), #MEAN START HERE
            meanDensity_avg = mean(meanDensity, na.rm = TRUE),
            meanBiomass_avg = mean(meanBiomass, na.rm = TRUE),
            trap_occurrence_percent_avg = mean(trap_occurrence_percent, na.rm = TRUE),
            frass_mass_avg = mean(frass_mass, na.rm = TRUE),
            number_pellets_avg = mean(number_pellets, na.rm = TRUE),
            weeklytemp_avg = mean(weeklytemp, na.rm = TRUE), 
         fracSurveys_anomaly = fracSurveys - fracSurveys_avg, #ANAOMALY START HERE
         meanDensity_anomaly = meanDensity - meanDensity_avg,
         meanBiomass_anomaly = meanBiomass - meanBiomass_avg,
         trap_occurrence_percent_anomaly = trap_occurrence_percent - trap_occurrence_percent_avg,
         frass_mass_anomaly = frass_mass - frass_mass_avg,
         number_pellets_anomaly = number_pellets - number_pellets_avg,
         weeklytemp_anomaly = weeklytemp - weeklytemp_avg,
         fracSurveys_sd = sd(fracSurveys, na.rm = TRUE), #STANDARD DEVIATION START HERE
         meanDensity_sd = sd(meanDensity, na.rm = TRUE),
         meanBiomass_sd = sd(meanBiomass, na.rm = TRUE),
         trap_occurrence_percent_sd = sd(trap_occurrence_percent, na.rm = TRUE),
         frass_mass_sd = sd(frass_mass, na.rm = TRUE),
         number_pellets_sd = sd(number_pellets, na.rm = TRUE),
         weeklytemp_sd = sd(weeklytemp, na.rm = TRUE),
         fracSurveys_zscore = (fracSurveys_anomaly/fracSurveys_sd), #Z SCORE START HERE
         meanDensity_zscore = (meanDensity_anomaly/meanDensity_sd),
         meanBiomass_zscore = (meanBiomass_anomaly/meanBiomass_sd),
         trap_occurrence_percent_zscore = (trap_occurrence_percent_anomaly/trap_occurrence_percent_sd),
         frass_mass_zscore = (frass_mass_anomaly/frass_mass_sd),
         number_pellets_zscore = (number_pellets_anomaly/number_pellets_sd),
         weeklytemp_zscore = (weeklytemp_anomaly/weeklytemp_sd))

#make a function for plotting each variable for all years and where they are all centered at 0
anomaly_centered_plotting <- function(data, var_choice, site_choice, x_var = julianweek, year_var = Year, site_var = Site) {
  
  # keep only the chosen site
  df1 <- data  %>% filter({{ site_var }} == .env$site_choice)

  ggplot(df1) +
    # one thin line per year
    geom_line(aes(x = {{ x_var }},
                  y = {{ var_choice }},
                  group = {{ year_var }},
                  color = factor({{ year_var }})),
              alpha = 0.7, linewidth = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    theme_classic()
}

#use it
anomaly_centered_plotting(anomaly_centered, trap_occurrence_percent_anomaly, 8892356)
anomaly_centered_plotting(anomaly_centered, trap_occurrence_percent_zscore, 8892356)


#saving as a pdf------------------------ ^^^^^
vars  <- c("fracSurveys_anomaly", "meanDensity_anomaly", "trap_occurrence_percent_anomaly",
           "meanBiomass_anomaly", "number_pellets_anomaly", "frass_mass_anomaly", "weeklytemp_anomaly")
vars_z  <- c("fracSurveys_zscore", "meanDensity_zscore", "trap_occurrence_percent_zscore",
           "meanBiomass_zscore", "number_pellets_zscore", "frass_mass_zscore", "weeklytemp_zscore")
sites <- c(8892356, 117)

setwd("C:/Z_School/school/HurlbertLab/graphs")
pdf("anomalytest.pdf", width = 8, height = 10)

for (s in sites) {
  plots <- map(vars_z, ~ anomaly_centered_plotting(
    anomaly_centered,
    !!rlang::sym(.x),     # turns the string into a column name
    site_choice = s       # must be named, since the 3rd positional arg is site_choice
  ))
  
  print(wrap_plots(plots, ncol = 2, nrow = 4, guides = "collect") +
          plot_annotation(title = paste("Site", s)))
}

dev.off()   # closes the PDF, otherwise it may be corrupt or empty

#--------------------------
#linear model for zscore anomaly
vars_z  <- c("fracSurveys_zscore", "meanDensity_zscore", "trap_occurrence_percent_zscore",
             "meanBiomass_zscore", "number_pellets_zscore", "frass_mass_zscore", "weeklytemp_zscore")

#lmer(response(all6)_anomaly ~ spring_temp_anomaly, data=df)
models <- lapply(vars_z, function(v) {
  lm(as.formula(paste(v, "~ weeklytemp_zscore")), data = anomaly_centered)
})
names(models) <- vars_z

lapply(models, summary)


# ---------------------------------------------------------------------------------------------
#   creating visualization for centered IMPUTED ANOMALY plots & z score ones
# -----------------------------------------------------------------------------------------------
#combining frass mass/occurrence and caterpillar occurrence, biomass, and density into one dataframe
all_six_variables_dataframe_impute <- cat_data_byweek %>%
  dplyr::select(Site, Year, julianweek, fracSurveys, meanDensity)%>%
  left_join(occurrence_frass_combined_weeks_imputation_data %>% dplyr::select(Site, Year, julianweek, trap_occurrence_percent),by=c("Site", "Year", "julianweek")) %>%
  left_join(imputation_data %>% dplyr::select(Site, Year, julianweek, meanBiomass, mass), by=c("Site", "Year", "julianweek")) %>%
  left_join(frass_pellet_number_imputation_data, dplyr::select(Site, Year, julianweek, number_pellets), by =c("Site", "Year", "julianweek")) %>%
  rename(frass_mass = mass) #all years previous standardized (same jday range before joining)

#ok lowkey use imputation function again on the combined dataframe to fill in missing trap_occurrence_percent and number pellets columns
all_six_variables_dataframe_impute2 <- impute_data(all_six_variables_dataframe_impute, trap_occurrence_percent) 

all_six_variables_dataframe_impute3 <- impute_data(all_six_variables_dataframe_impute2, number_pellets)

#make all impute columns
anomaly_centered_impute <- all_six_variables_dataframe_impute3 %>%
  left_join(AllTemp_weekly, by=c("Year", "Site", "julianweek")) %>%
  group_by(Site, Year) %>%
  mutate(fracSurveys_avg = mean(fracSurveys, na.rm = TRUE), #MEAN START HERE
         meanDensity_avg = mean(meanDensity, na.rm = TRUE),
         meanBiomass_avg = mean(meanBiomass, na.rm = TRUE),
         trap_occurrence_percent_avg = mean(trap_occurrence_percent, na.rm = TRUE),
         frass_mass_avg = mean(frass_mass, na.rm = TRUE),
         number_pellets_avg = mean(number_pellets, na.rm = TRUE),
         weeklytemp_avg = mean(weeklytemp, na.rm = TRUE), 
         fracSurveys_anomaly = fracSurveys - fracSurveys_avg, #ANAOMALY START HERE
         meanDensity_anomaly = meanDensity - meanDensity_avg,
         meanBiomass_anomaly = meanBiomass - meanBiomass_avg,
         trap_occurrence_percent_anomaly = trap_occurrence_percent - trap_occurrence_percent_avg,
         frass_mass_anomaly = frass_mass - frass_mass_avg,
         number_pellets_anomaly = number_pellets - number_pellets_avg,
         weeklytemp_anomaly = weeklytemp - weeklytemp_avg,
         fracSurveys_sd = sd(fracSurveys, na.rm = TRUE), #STANDARD DEVIATION START HERE
         meanDensity_sd = sd(meanDensity, na.rm = TRUE),
         meanBiomass_sd = sd(meanBiomass, na.rm = TRUE),
         trap_occurrence_percent_sd = sd(trap_occurrence_percent, na.rm = TRUE),
         frass_mass_sd = sd(frass_mass, na.rm = TRUE),
         number_pellets_sd = sd(number_pellets, na.rm = TRUE),
         weeklytemp_sd = sd(weeklytemp, na.rm = TRUE),
         fracSurveys_zscore = (fracSurveys_anomaly/fracSurveys_sd), #Z SCORE START HERE
         meanDensity_zscore = (meanDensity_anomaly/meanDensity_sd),
         meanBiomass_zscore = (meanBiomass_anomaly/meanBiomass_sd),
         trap_occurrence_percent_zscore = (trap_occurrence_percent_anomaly/trap_occurrence_percent_sd),
         frass_mass_zscore = (frass_mass_anomaly/frass_mass_sd),
         number_pellets_zscore = (number_pellets_anomaly/number_pellets_sd),
         weeklytemp_zscore = (weeklytemp_anomaly/weeklytemp_sd))

#make a function for plotting each variable for all years and where they are all centered at 0
anomaly_centered_plotting <- function(data, var_choice, site_choice, x_var = julianweek, year_var = Year, site_var = Site) {
  
  # keep only the chosen site
  df1 <- data  %>% filter({{ site_var }} == .env$site_choice)
  
  ggplot(df1) +
    # one thin line per year
    geom_line(aes(x = {{ x_var }},
                  y = {{ var_choice }},
                  group = {{ year_var }},
                  color = factor({{ year_var }})),
              alpha = 0.7, linewidth = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    theme_classic()
}

#use it
anomaly_centered_plotting(anomaly_centered_impute, trap_occurrence_percent_anomaly, 8892356)
anomaly_centered_plotting(anomaly_centered_impute, trap_occurrence_percent_zscore, 8892356)


#saving as a pdf------------------------ ^^^^^
vars  <- c("fracSurveys_anomaly", "meanDensity_anomaly", "trap_occurrence_percent_anomaly",
           "meanBiomass_anomaly", "number_pellets_anomaly", "frass_mass_anomaly", "weeklytemp_anomaly")
vars_z  <- c("fracSurveys_zscore", "meanDensity_zscore", "trap_occurrence_percent_zscore",
             "meanBiomass_zscore", "number_pellets_zscore", "frass_mass_zscore", "weeklytemp_zscore")
sites <- c(8892356, 117)

setwd("C:/Z_School/school/HurlbertLab/graphs")
pdf("anomalytest_impute_zscore.pdf", width = 8, height = 10)

for (s in sites) {
  plots <- map(vars_z, ~ anomaly_centered_plotting(
    anomaly_centered_impute,
    !!rlang::sym(.x),     # turns the string into a column name
    site_choice = s       # must be named, since the 3rd positional arg is site_choice
  ))
  
  print(wrap_plots(plots, ncol = 2, nrow = 4, guides = "collect") +
          plot_annotation(title = paste("Site", s)))
}

dev.off()   

#--------------------------
#linear model for zscore anomaly
vars_z  <- c("fracSurveys_zscore", "meanDensity_zscore", "trap_occurrence_percent_zscore",
             "meanBiomass_zscore", "number_pellets_zscore", "frass_mass_zscore", "weeklytemp_zscore")

#lmer(response(all6)_anomaly ~ spring_temp_anomaly, data=df)
models <- lapply(vars_z, function(v) {
  lm(as.formula(paste(v, "~ weeklytemp_zscore")), data = anomaly_centered_impute)
})
names(models) <- vars_z

lapply(models, summary)



#SCATTERPLOT OF THINGS WITH ABLINE INCLUDED BY COLOR




















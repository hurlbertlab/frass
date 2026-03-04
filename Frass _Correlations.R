## looking into correlations for Frass data
#Savannah Carter
#R 4.4.2

#libraries:
library(stringr)
library(tidyverse)
library(readr)
library(data.table)
library(dplyr)
library(lubridate)
library(readxl)
library(gsheet)
library(tidyr)
library(purrr)
library(rstatix)
library(tidyverse)
library(jsonlite)
library(daymetr)
library(lme4)


# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#   Datasets needed:
# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+

#Loading in CC data 
options(timeout = 300)  
api_url <- "https://api.github.com/repos/hurlbertlab/caterpillars-analysis-public/contents/data"
files <- fromJSON(api_url)
dataset_file <- files$name[grepl("fullDataset", files$name, ignore.case = TRUE)]
# pick the latest one
latest_file <- dataset_file[1]
github_raw <- "https://raw.githubusercontent.com/hurlbertlab/caterpillars-analysis-public/master/data/"
fullDataset <- read.csv(paste0(github_raw, latest_file))

#----------------------------------------------------------------------------------
#loading in temp data
tmp_file = tempfile(fileext = ".csv")

AnomalySites = fullDataset[, c("Name","Year", "Latitude", "Longitude")] %>% 
  filter(Name %in% c("NC Botanical Garden", "Prairie Ridge Ecostation")) %>% 
  group_by(Name, Latitude, Longitude) %>% 
  summarise(n = n()) %>% 
  select(-n) %>% as.data.frame()   

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
  end = 2024, # this is the most recent available in daymetr
  internal = TRUE
)

# remove temporary file 
unlink(tmp_file)

TempAnomalyData_clean <- lapply(TempAnomaly, function(x) {
  x$data %>% mutate(site = x$site,
                    Latidue = x$latitude,
                    Longitude = x$longitude)})

AllTemp= bind_rows(TempAnomalyData_clean) #this is the file we want (has all NCBG and PR data)

#----------------------------------------------------------------------------------
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
#--------------------------------------------------------------------------------------------------------------------------------
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
#--------------------------------------------------------------------------------------------------------------------------------
# using data to find mean frass per day for reliable frass only 
#read in proper url, change dates and label events for below meanfrass
url = "https://docs.google.com/spreadsheets/d/1RwXzwhHUbP0m5gKSOVhnKZbS1C_NrbdfHLglIVCzyFc/edit#gid=1611171427"
events = gsheet2tbl(url)
events$date = as.Date(events$date, format = "%m/%d/%Y")

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

#--------------------------------------------------------------------------------------------------------------------------------
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

# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#   Combing meanfrass, temp, and fulldataset into one data sheet:
# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#correcting functions so can combine properly----------------------------
#having mean Frass data sorted by julian week
meanfrass <- meanfrass %>%
  mutate(julianweek = 7 * floor(jday / 7) + 4)
#make sure year is a integar
meanfrass <- meanfrass %>%
  mutate(Year = as.integer(Year))
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
    
    .groups = "drop"
  )
#filter meanfrass by NCBG
meanfrass_NCBG <- meanfrass_combinedweeks %>%
  filter(site %in% c("8892356"))
#filter meanfrass by PR
meanfrass_PR <- meanfrass_combinedweeks %>%
  filter(site %in% c("117"))
#NCBG site filter fulldataset for all years
NCBG <- fullDataset %>%
  filter(Name %in% c("NC Botanical Garden"),
         Year %in% 2015:2025)
#PR site filter fulldataset for all years
PR <- fullDataset %>%
  filter(Name %in% c("Prairie Ridge Ecostation"),
         Year %in% 2015:2025)
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
  })
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
  })
#left join fulldataset and meanfrass by julianweek
cats_NCBG <- meanfrass_NCBG %>%
  left_join(cats_NCBG, by = c("julianweek", "Year"))
cats_PR <- meanfrass_PR %>%
  left_join(cats_PR, by = c("julianweek", "Year"))  
#rename columns to make sense:
cats_NCBG <- rename(cats_NCBG, frass_mass=mass)
cats_NCBG <- rename(cats_NCBG, frass_density=density)
cats_NCBG <- rename(cats_NCBG, frass_reliability=reliability)
cats_PR <- rename(cats_PR, frass_mass=mass)
cats_PR <- rename(cats_PR, frass_density=density)
cats_PR <- rename(cats_PR, frass_reliability=reliability)

# Temp stuff--------------------------------------
#rename site names so match other dataframes, also rename to can be joined properly with meanfrass data
AllTemp <- AllTemp %>%
  mutate(site = case_when(
    site == "NC Botanical Garden" ~ "117",
    site == "Prairie Ridge Ecostation" ~ "8892356"  )) 
AllTemp <- rename(AllTemp, jday=yday)
AllTemp<- rename(AllTemp, Year=year)
#add average temp column to alltemp data
AllTemp <- AllTemp %>%
  mutate(avgtemp = (tmax..deg.c. + tmin..deg.c.) / 2)  #avg max and min temps
#adding average weekly temp column to alltemp data, and averaging temp for all days with same jweek value
AllTemp_weekly <- AllTemp %>%
  mutate(julianweek = 7 * floor(jday / 7) + 4) %>%
  group_by(Year, site, julianweek) %>%
  summarise(weeklytemp = mean(avgtemp, na.rm = TRUE),
            .groups = "drop")
#combining temp data and frass data into one dataset
Temp_with_frass <- meanfrass_combinedweeks %>%
  left_join(AllTemp_weekly, by = c("julianweek", "Year", "site"))
#first biomass estimate (using tinbergen2024 HV)
Tinbergen_biomass <- Temp_with_frass %>%
  mutate(Tin_biomass = mass * exp(3.8 - 0.10 * weeklytemp))

# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#   combining caterpillar data with tinbergen data set so all in one df
# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+ 
#create nontemp corrected biomass dataframe:
cats_all <- bind_rows(cats_NCBG, cats_PR)
#combine tinbergen biomass with nontemp corrected dataframe by jweek, site, year
cats_tinbergen_biomass <- Tinbergen_biomass %>%
  left_join(
    cats_all %>%
      select(julianweek, site, Year, meanBiomass, nSurveys), #meanBiomass is NOT temp corrected
    by = c("julianweek", "site", "Year")
  )
#cutoff dates that are not shared between all years
cats_tinbergen_biomass_cutoff <- cats_tinbergen_biomass %>%
  filter(
    (site == 117 & between(jday, 142, 200)) |
      (site == 8892356 & between(jday, 154, 198)) |
      !(site %in% c(117, 8892356))
  )
#density values
cats_tinbergen_biomass_density <- cats_tinbergen_biomass_cutoff %>%
  mutate(Tin_biomass_density = Tin_biomass/(ifelse(Year <= 2018, 309.74, 197.71))) %>% #dividing by 209 for years 2018 and before, and 197 for years after
  mutate(biomass_density = meanBiomass/ (ifelse(Year <= 2018, 309.74, 197.71)))



# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#   testing if width / height of tinbergen biomass density AND non temp corrected peaks different significantly from year to year
# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+ 
#make sure all_centroid is this one:
all_centroids <- cats_tinbergen_biomass_density %>%
  group_by(Year, site) %>%
  summarise(centroid_frass_density = weighted.mean(jday, density, na.rm = TRUE),
            centroid_TinbergenBiomass_density =weighted.mean(jday, Tin_biomass_density, na.rm = TRUE),
            centroid_ActualBiomass_density =weighted.mean(jday, biomass_density, na.rm = TRUE))%>%
  mutate(diff_centroid = centroid_TinbergenBiomass_density - centroid_ActualBiomass_density) %>%
  mutate(start_day = case_when(
    site == 8892356 ~ 154,
    site == 117 ~ 142)) %>%
  mutate(bin_frass_density = floor((centroid_frass_density - start_day) / 3) + 1) %>%
  mutate(bin_tinbergenbiomass_density = floor((centroid_TinbergenBiomass_density - start_day) / 3) + 1) %>%
  mutate(bin_Actualbiomass = floor((centroid_ActualBiomass_density - start_day) / 3) + 1) 
#Is centroid timing shifting earlier or later over time? Now the Year coefficient tells you: Positive slope → peak happening later. Negative slope → peak happening earlier
#tinbergen biomass-------------------------
all_centroids$Year_c <- all_centroids$Year - mean(all_centroids$Year) #center year so easier to interpret data
model1 <- lm(centroid_ActualBiomass_density ~ Year_c * factor(site),
             data = all_centroids)
summary(model1) #centroid∼Yearc​×site
#frass density-------------------------------
all_centroids$Year_c <- all_centroids$Year - mean(all_centroids$Year) #center year so easier to interpret data
model2 <- lm(centroid_frass_density ~ Year_c * factor(site),
             data = all_centroids)
summary(model2) #centroid∼Yearc​×site

# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#   Plotting within sites: correlation between frass density and actual biomass density at each site (looking at cutoffs)
# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+ 
#calculating spearmans
spearmans <- function(dataset) {
  dataset %>%
    group_by(site, Year) %>%
    filter(
      n() >= 2,
      !all(is.na(density)),
      !all(is.na(biomass_density))) %>%
    summarise(
      rho = cor(density, biomass_density, method = "spearman"),
      p_value = cor.test(density, biomass_density, method = "spearman")$p.value,
      n = n(),
      .groups = "drop")}
spearmans_results <- spearmans(cats_tinbergen_biomass_density)

#plot results, Here filled circles = p < 0.05, open circles = not significant
plot(
  NA,
  xlim = range(spearmans_results$Year),
  ylim = c(-1, 1),
  xlab = "Year",
  ylab = "Spearman's rho",
  main = "Spearman Correlation Through Time"
)
abline(h = 0, lty = 2, col = "gray50")
# 2. Define sites and colors
sites <- unique(spearmans_results$site)
cols <- c(
  "117" = "lightblue4",
  "8892356" = "darkred"
)# 3. Loop over sites
for (i in seq_along(sites)) {
  
  sub <- spearmans_results[spearmans_results$site == sites[i], ]
  sub <- sub[order(sub$Year), ]  # important!
  
  lines(sub$Year, sub$rho, col = cols[i], lwd = 2)
  
  points(
    sub$Year,
    sub$rho,
    col = cols[i],
    pch = ifelse(sub$p_value < 0.05, 16, 1),
    cex = 1.2)}
# 4. Legend
legend(
  "topleft",
  legend = sites,
  col = cols,
  lwd = 2,
  pch = 16,
  title = "Site",
  bty = "n"
)
#all years liner model data
model1 <- lm(biomass_density ~ density + factor(Year) + factor(site),
             data = cats_tinbergen_biomass_density)

summary(model1)



# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#   Plotting between sites: correlation between frass density and actual biomass between years that line up at NCBG and PR
# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+ 
#shared years
shared_years <- c(2015, 2018, 2019, 2021, 2022)
#single correlation across both sites & all shared years: one overall Spearman for those two sites during shared years:
shared_data <- cats_tinbergen_biomass_density %>%
  filter(
    site %in% c(117, 8892356),
    Year %in% shared_years)

cor.test(
  shared_data$density,
  shared_data$biomass_density,
  method = "spearman")

##add year as fixed effect:
shared_data <- cats_tinbergen_biomass_density %>%
  filter(
    site %in% c(117, 8892356),
    Year %in% shared_years
  )

model1 <- lm(biomass_density ~ density + factor(Year) + factor(site),
             data = shared_data)

summary(model1)


# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#   Plotting between sites: adding in temp and precipitation data
# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+ 
#modify AllTemp to include precipitation average up for the week 
AllTemp_weekly_prcp <- AllTemp %>%
  mutate(julianweek = 7 * floor(jday / 7) + 4) %>%
  group_by(Year, site, julianweek) %>%
  summarise(
    weeklytemp = mean(avgtemp, na.rm = TRUE),
    avg_prcp   = mean(prcp..mm.day., na.rm = TRUE),
    .groups = "drop"
  )
#combined new precipation data to shared_data so has all weather stuff
shared_data_new <- shared_data %>%
  left_join(
    AllTemp_weekly_prcp %>%
      select(Year, site, julianweek, avg_prcp),
    by = c("Year", "site", "julianweek"))
#center data
shared_data_new <- shared_data_new %>%
  mutate(
    density_c = scale(density, scale = FALSE),
    prcp_c = scale(avg_prcp, scale = FALSE),
    temp_c = scale(weeklytemp, scale = FALSE)
  )
##Add an actual climate variable---------------------------------------
model_mixed <- lmer(biomass_density ~ density_c + prcp_c + temp_c +
                      (1 | site) + (1 | Year),
                    data = shared_data_new)
summary(model_mixed)
#plotting
plot(fitted(model_mixed), resid(model_mixed))
abline(h = 0, col = "red")
#checking qqplots
qqnorm(resid(model_mixed))
qqline(resid(model_mixed)) #indicating right skewed as some large frass measurements are not captured by line
ranef(model_mixed)

#adding log transformation to account for skew ---------------------------------
model_log <- lmer(
  log(biomass_density) ~ density_c + prcp_c + temp_c +
    (1 | site) + (1 | Year),
  data = shared_data_new)
summary(model_log)
#plot
plot(model_log)
qqnorm(resid(model_log))
qqline(resid(model_log))












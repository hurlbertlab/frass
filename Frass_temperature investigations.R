## Investigating temperature as explanatory variable for varince in caterillar and frass patterns over time
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
library(gridExtra)

# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#   Datasets needed:
# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#Loading in CC data ------------------------------------------------------
options(timeout = 300)  
api_url <- "https://api.github.com/repos/hurlbertlab/caterpillars-analysis-public/contents/data"
files <- fromJSON(api_url)
dataset_file <- files$name[grepl("fullDataset", files$name, ignore.case = TRUE)]
# pick the latest one
latest_file <- dataset_file[1]
github_raw <- "https://raw.githubusercontent.com/hurlbertlab/caterpillars-analysis-public/master/data/"
fullDataset <- read.csv(paste0(github_raw, latest_file))

#loading in temp data (from Nosa)------------------------------------------------------
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



#loading in temp data (from Nosa) 26 YEARS------------------------------------------------------
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
TempAnomaly_long = download_daymet_batch(
  file_location = tmp_file,
  start = 2000, #ALTER TO INCLUDE DATA FROM LAST 26 YEARS
  end = 2025, # this is the most recent available in daymetr
  internal = TRUE
)
# remove temporary file 
unlink(tmp_file)

TempAnomalyData_clean <- lapply(TempAnomaly_long, function(x) {
  x$data %>% mutate(site = x$site,
                    Latidue = x$latitude,
                    Longitude = x$longitude)})

AllTemp_long= bind_rows(TempAnomalyData_clean) #this is the file we want (has all NCBG and PR data)
#clean up globals
rm(AnomalySites, TempAnomaly_long, files, TempAnomalyData_clean)

# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#   altering temperature to be what sites I want and have jday and jweek columns
# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+ 
#rename site names so match other dataframes, also rename so can be joined properly with meanfrass data
AllTemp <- AllTemp %>%
  mutate(site = case_when(
    site == "NC Botanical Garden" ~ "117",
    site == "Prairie Ridge Ecostation" ~ "8892356"  )) 
AllTemp <- rename(AllTemp, jday=yday)
AllTemp<- rename(AllTemp, Year=year)
#add average temp column to alltemp data and nonoptimal column
AllTemp <- AllTemp %>%
  mutate(avgtemp = (tmax..deg.c. + tmin..deg.c.) / 2)%>%#avg max and min temps 
  mutate(nonoptimal = ifelse(tmax..deg.c.>= 40, 1, 0 )) %>% #add binary column for nonoptimal temp (0 if no 1 if >=40)
  mutate(optimal = ifelse(tmax..deg.c.>=32, 1, 0))



#adding average weekly temp column to alltemp data, and averaging temp for all days with same jweek value
AllTemp_weekly <- AllTemp %>%
  mutate(julianweek = 7 * floor(jday / 7) + 4) %>%
  group_by(Year, site, julianweek) %>%
  summarise(weeklytemp = mean(avgtemp, na.rm = TRUE),
            .groups = "drop")

# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#Do years with more optimal days have worse correlations between frass and cat variables?
# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+ 
#create one df with all years summed optimal days and column for mean optimal day
AllTemp_long <- AllTemp_long %>%
  mutate(site = case_when(
    site == "NC Botanical Garden" ~ "117",
    site == "Prairie Ridge Ecostation" ~ "8892356"  )) 
AllTemp_long <- rename(AllTemp_long, jday=yday)
AllTemp_long<- rename(AllTemp_long, Year=year)
#sum and make mean column
summed_optimal_long <- AllTemp_long %>%
  mutate(optimal = ifelse(tmax..deg.c.>=32, 1, 0))%>%
  group_by(site, Year) %>%
  summarize(total_optimal = sum(optimal)) %>%
  mutate(mean_optimal_days = mean(total_optimal))
#based on total data from last 26ish years is the data point for an individual year significantly different from site mean?
years_PR   <- c(2015, 2018, 2019, 2021, 2022)
years_NCBG <- setdiff(2015:2026, 2020) 

#PR
for (yr in years_PR) {
  dat <- summed_optimal_long[summed_optimal_long$Year == yr & summed_optimal_long$site == 117, ]
  
  comparison <- t.test(dat$total_optimal, dat$mean_optimal_days)
  
  cat("\nYear:", yr, "\n")
  print(comparison)}

#NCBG
for (yr in years_NCBG) {
  dat <- summed_optimal_long[summed_optimal_long$Year == yr, ]
  
  comparison <- t.test(dat$total_optimal, dat$mean_optimal_days)
  
  cat("\nYear:", yr, "\n")
  print(comparison)}













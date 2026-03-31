##Looking at weather in region for Frass
#savannah Carter
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

# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#   altering temperature to be what sites I want and have jday and jweek columns
# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+ 
#rename site names so match other dataframes, also rename to can be joined properly with meanfrass data
AllTemp <- AllTemp %>%
  mutate(site = case_when(
    site == "NC Botanical Garden" ~ "117",
    site == "Prairie Ridge Ecostation" ~ "8892356"  )) 
AllTemp <- rename(AllTemp, jday=yday)
AllTemp<- rename(AllTemp, Year=year)
#add average temp column to alltemp data
AllTemp <- AllTemp %>%
  mutate(avgtemp = (tmax..deg.c. + tmin..deg.c.) / 2)  #avg max and min temps#adding average weekly temp column to alltemp data, and averaging temp for all days with same jweek value
AllTemp_weekly <- AllTemp %>%
  mutate(julianweek = 7 * floor(jday / 7) + 4) %>%
  group_by(Year, site, julianweek) %>%
  summarise(weeklytemp = mean(avgtemp, na.rm = TRUE),
            .groups = "drop")

# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#   making temperature graphs for each year of weekly average temp to help explain drops in total frass/biomass abundance
# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+ 
plotting_temp_yearly <- function(data, year_choice, site_choice) { 
  # ---- Filter by year and site ----
  df <- data %>%
    filter(Year == year_choice, site == site_choice)
  
  if(nrow(df) == 0) {
    stop("No data for the specified year and site.")
  }
  
  # ---- Plot ----
  par(mar = c(5, 4, 4, 4))
  
  plot(
    df$julianweek,
    df$weeklytemp,
    type = "l",
    col = "cornflowerblue",
    lwd = 2,
    xlab = "Julian week",
    ylab = "Weekly average temp",
    main = paste(site_choice, year_choice))
  #add lines at start and end of smapling period for each site
  if(site_choice == 117) {
    abline(v = c(142, 200), lty = 5, col="darkred")
  }
  if(site_choice == 8892356){
    abline(v = c(154, 198), lty = 5, col="darkred")
  }
  
  legend(
    "topleft",
    legend = c("weekly temp"),
    col = c("cornflowerblue"),
    lwd = 2,
    lty = 1,
    bty = "n"
  )
  invisible(df)}
plotting_temp_yearly(AllTemp_weekly, 2016, 8892356)
#SAVING TO A PDF-------------------------------
#Years for each site
years_NCBG <- c(2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022,2023, 2024)
years_PR   <- c(2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022)
setwd("C:/Z_School/HurlbertLab/graphs")
#set pdf size
pdf(file = "plotting_temp.pdf",
    width = 11,
    height = 12)
#set up how many plots per page
par(mfrow = c(3, 2),     
    mar = c(3, 4, 3, 4), 
    oma = c(0, 0, 2, 0)) 
#loop through all plots and years
# NCBG plots
# NCBG plots
for (yr in years_NCBG) {
  try(
    plotting_temp_yearly(
      data = AllTemp_weekly,
      year_choice = yr,        # was years_NCBG, should be yr
      site_choice = 8892356),
    silent = TRUE)}
# PR plots
for (yr in years_PR) {
  try(
    plotting_temp_yearly(
      data = AllTemp_weekly,
      year_choice = yr,        # was years_PR, should be yr
      site_choice = 117),
    silent = TRUE)}
dev.off()

# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#   average temp per week across all years for each site
# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+ 
AllTemp_byweek <- AllTemp_weekly %>% #avg jweek temp from years 2015-2024 for each site
  group_by(site, julianweek) %>%
  summarise(week_avg = mean(weeklytemp, na.rm = TRUE))

#do temp anomalies (then plot them)
temp_anomolies <- AllTemp_weekly %>%
  left_join(AllTemp_byweek, by =c("julianweek", 'site')) %>%
  mutate(diff = weeklytemp-week_avg) #how different year weeklytemp is from average across all years during 2015-2024

#plot anomoliesby site per year
# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+ 
plotting_temp_anomaly <- function(data, year_choice, site_choice) { 
  # ---- Filter by year and site ----
  df <- data %>%
    filter(Year == year_choice, site == site_choice)
  
  if(nrow(df) == 0) {
    stop("No data for the specified year and site.")
  }
  
  # ---- Plot ----
  par(mar = c(5, 4, 4, 4))
  
  plot(
    df$julianweek,
    df$diff,
    type = "l",
    col = "cornflowerblue",
    lwd = 2,
    xlab = "Julian week",
    ylab = "Temp anomaly",
    main = paste(site_choice, year_choice))
  #add lines at start and end of sampling period for each site
  if(site_choice == 117) {
    abline(v = c(142, 200), lty = 5, col="darkred")
  }
  if(site_choice == 8892356){
    abline(v = c(154, 198), lty = 5, col="darkred")
  }
  #adding y=0 line
  abline(h=0, lty=1, col="grey21")
  
  legend(
    "topleft",
    legend = c("Temp anomaly"),
    col = c("cornflowerblue"),
    lwd = 2,
    lty = 1,
    bty = "n"
  )
  invisible(df)}
plotting_temp_anomaly(temp_anomolies, 2020, 117)
#SAVING TO A PDF-------------------------------
#Years for each site
years_NCBG <- c(2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022,2023, 2024)
years_PR   <- c(2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022)
setwd("C:/Z_School/HurlbertLab/graphs")
#set pdf size
pdf(file = "plotting_temp_anomolies.pdf",
    width = 11,
    height = 12)
#set up how many plots per page
par(mfrow = c(3, 2),     
    mar = c(3, 4, 3, 4), 
    oma = c(0, 0, 2, 0)) 
#loop through all plots and years
# NCBG plots
# NCBG plots
for (yr in years_NCBG) {
  try(
    plotting_temp_anomaly(
      data = temp_anomolies,
      year_choice = yr,        # was years_NCBG, should be yr
      site_choice = 8892356),
    silent = TRUE)}
# PR plots
for (yr in years_PR) {
  try(
    plotting_temp_anomaly(
      data = temp_anomolies,
      year_choice = yr,        # was years_PR, should be yr
      site_choice = 117),
    silent = TRUE)}
dev.off()




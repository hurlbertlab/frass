#frass imputation stuff in own file
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
library(zoo)

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

##same thing but by day because then can have individual Jday calulations to do imputation on
meanDensityByDay = function(surveyData,
                            ordersToInclude = 'All',
                            minLength = 0,
                            jdRange = c(1,365),
                            outlierCount = 10000,
                            plot = FALSE,
                            plotVar = 'fracSurveys',
                            minSurveyCoverage = 0.8,
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
    filter(julianday >= jdRange[1], julianday <= jdRange[2])
  # removed the julianweek mutate line
  
  effortByDay = firstFilter %>%
    group_by(julianday) %>%                          # group by julianday instead of julianweek
    summarize(nSurveyBranches = n_distinct(PlantFK),
              nSurveys = n_distinct(ID)) %>%
    mutate(modalBranchesSurveyed = Mode(5*ceiling(nSurveyBranches/5)),
           nSurveySets = nSurveys/modalBranchesSurveyed,
           modalSurveySets = Mode(round(nSurveySets)),
           okWeek = ifelse(nSurveySets/modalSurveySets >= minSurveyCoverage, 1, 0))
  
  if (allDates) {
    effortByDay$okWeek = 1
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
    mutate(Quantity2 = ifelse(Quantity > outlierCount, 1, Quantity)) %>%
    group_by(julianday) %>%                          # group by julianday instead of julianweek
    summarize(totalCount = sum(Quantity2, na.rm = TRUE),
              numSurveysGTzero = length(unique(ID[Quantity > 0])),
              totalBiomass = sum(Biomass_mg, na.rm = TRUE)) %>% 
    right_join(effortByDay, by = 'julianday') %>%    # join by julianday
    filter(okWeek == 1) %>%
    mutate_cond(is.na(totalCount), totalCount = 0, numSurveysGTzero = 0, totalBiomass = 0) %>%
    mutate(meanDensity = totalCount/nSurveys,
           fracSurveys = 100*numSurveysGTzero/nSurveys,
           meanBiomass = totalBiomass/nSurveys) %>%
    arrange(julianday) %>%
    data.frame()
  
  if (plot & new) {
    plot(arthCount$julianday, arthCount[, plotVar], type = 'l', 
         col = color, las = 1, ...)
    points(arthCount$julianday, arthCount[, plotVar], pch = 16, col = color, ...)
  } else if (plot & new==F) {
    points(arthCount$julianday, arthCount[, plotVar], type = 'l', col = color, ...)
    points(arthCount$julianday, arthCount[, plotVar], pch = 16, col = color, ...)
  }
  return(arthCount)
}

# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#   Creating dataframe with caterpillar meanbiomass and frass values:
# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
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
      (site == 8892356 & Year %in% c(2015:2019, 2021:2025) & julianweek %in% 154:198)
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


# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#   Total season estimations for frass and biomass for NON imputated data
# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+ 
#use number of surveys from imputation data
check_date <- all_data_clean %>%
  group_by(site, Year) %>%
  summarise(
    number_surveys_bugs = sum(!is.na(meanBiomass)),  # counts non-NA rows for biomass
    number_surveys_frass = sum(!is.na(mass)),         # counts non-NA rows for frass
    .groups = "drop"
  )
#add up mass column and mean biomass in the same year
season_totals <- all_data_clean %>%
  group_by(site, Year) %>%
  mutate(
    total_frass = sum(frass_density, na.rm = TRUE), #doing trap/area counts
    total_actualbiomass = sum(biomass_density, na.rm = TRUE)
  ) %>%
  select(site, Year, total_frass, total_actualbiomass) %>%
  distinct() %>%
  left_join(check_date, by = c("Year", "site")) %>% #left join with check data so can divide by the number of surveys to standardize across years 
  mutate(
    total_frass_divided = total_frass / number_surveys_frass,
    total_actualbiomass_divided = total_actualbiomass / number_surveys_bugs
  )
#plot total estimations for frass and biomass across years
par(mar = c(5, 4, 4, 5))   # increase right margin
plot_seasonal_estimates <- function(data, site_name){
  # Filter data for the chosen site
  site_data <- data[data$site == site_name, ]
  # First plot (left axis)
  plot(site_data$Year,
       site_data$total_frass_divided,
       type="b",
       col="sienna",
       ylab="Frass Density Totals",
       xlab="Year",
       main=paste("Site:", site_name))
  
  # Allow second plot on same figure
  par(new=TRUE)
  # Second plot (right axis)
  plot(site_data$Year,
       site_data$total_actualbiomass_divided,
       type="b",
       col="forestgreen",
       axes=FALSE,
       xlab="",
       ylab="")
  axis(4)
  mtext("Biomass Density Totals", side=4, line=3)
  
  legend("topleft",
         legend=c("Frass","Biomass"),
         col=c("sienna","forestgreen"),
         lty=1,
         pch=1, cex=0.8)}
plot_seasonal_estimates(season_totals, 8892356) #totals from frass/area and biomass/area

# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#   Not imputed basic plotting charts (RAW DATA) 
# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+ 
#plot on linecharts with centroid dates
density_plotting <- function(data, year_choice, site_choice) {
  
  df <- data %>%
    filter(Year == year_choice, site == site_choice)
  
  if (nrow(df) == 0 ||
      all(is.na(df$julianweek)) ||
      all(is.na(df$frass_density))) {
    
    plot.new()
    text(0.5, 0.5,
         paste("No data for", site_choice, year_choice),
         cex = 1.2)
    return(invisible(NULL))
  }
  
  ## ---- Centroids ----
  frass_density_centroid <- weighted.mean(
    df$julianweek, df$frass_density, na.rm = TRUE
  )
  
  biomass_centroid <- weighted.mean(
    df$julianweek, df$biomass_density, na.rm = TRUE
  )
  
  ## ---- Plot ----
  par(mar = c(5, 4, 4, 6))  # space for one right axis
  
  # LEFT AXIS — FRASS
  plot(
    df$julianweek, df$frass_density,
    type = "l",
    col = "sienna",
    lwd = 2,
    xlab = "Julian week",
    ylab = "Frass density",
    main = paste(site_choice, year_choice)
  )
  
  if (is.finite(frass_density_centroid)) {
    abline(v = frass_density_centroid, col = "sienna", lty = 2, lwd = 2)
  }
  
  ## ---- RIGHT AXIS — BIOMASS (NOT TEMP CORRECTED) ----
  par(new = TRUE)
  
  plot(
    df$julianweek, df$biomass_density,
    type = "l",
    col = "forestgreen",
    lwd = 2,
    axes = FALSE,
    xlab = "",
    ylab = "",
    ylim = range(df$biomass_density, na.rm = TRUE)
  )
  
  axis(side = 4, col.axis = "forestgreen")
  mtext("Actual Caterpillar Biomass Density",
        side = 4, line = 2.0, col = "forestgreen", cex=0.8)
  
  if (is.finite(biomass_centroid)) {
    abline(v = biomass_centroid, col = "forestgreen", lty = 2, lwd = 2)
  }
  
  ## ---- Legend ----
  legend(
    "topleft",
    legend = c(
      "Frass Mass density",
      "Frass centroid",
      "Actual Caterpillar Biomass density",
      "Biomass centroid"
    ),
    col = c(
      "sienna", "sienna",
      "forestgreen", "forestgreen"
    ),
    lwd = 2,
    lty = c(1, 2, 1, 2),
    bty = "n",
    cex = 0.8
  )
  
  invisible(df)
}
density_plotting(all_data_clean, 2016, 8892356)
#saving as a pdf------------------------ ^^^^^
# Years for each site
years_PR   <- c(2015, 2018, 2019, 2021, 2022)
years_NCBG <- setdiff(2015:2025, 2020)   
setwd("C:/Z_School/HurlbertLab/graphs")
#set up pdf
pdf(
  file = "rawdata_plotting.pdf",
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
    density_plotting(
      data = all_data_clean,
      year_choice = yr,
      site_choice = 8892356  
    ),
    silent = TRUE)}
for (yr in years_PR) {
  try(
    density_plotting(
      data = all_data_clean,
      year_choice = yr,
      site_choice = 117   
    ),
    silent = TRUE)}
dev.off()


#######################################################################################################
#######################################################################################################


# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#   Doing imputation by julianweek, only if whole week missing will we impute data
# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+ 
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

# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#   Total season estimations for frass and biomass for imputated data
# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+ 
#count number of surveys done per site per year
check_date <- imputation_data %>%
  group_by(site, Year) %>%
  summarise(
    number_surveys_bugs = sum(!is.na(meanBiomass)),  # counts non-NA rows for biomass
    number_surveys_frass = sum(!is.na(mass)),         # counts non-NA rows for frass
    .groups = "drop"
  )

#collapse duplicate julian weeks into one week with one mean value, compute per area calculations
imputation_totals <- imputation_data %>%
  group_by(site, Year, julianweek) %>%
  summarise(
    meanBiomass = mean(meanBiomass, na.rm = TRUE),
    mass        = mean(mass, na.rm = TRUE),
    imputed_frass     = max(imputed_mass, na.rm = TRUE),
    imputed_biomass     = max(imputed_biomass, na.rm = TRUE),
    .groups = "drop") %>%
  mutate(biomass_density = meanBiomass/(ifelse(Year <= 2018, 309.74, 197.71))) %>% #dividing by 209 for years 2018 and before, and 197 for years after
  mutate(frass_density = mass/ (ifelse(Year <= 2018, 309.74, 197.71)))

#calculate seasonal totals by site year
season_totals_impute <- imputation_totals %>%
  group_by(site, Year) %>%
  mutate(total_frass = sum(frass_density), #sums from density
         total_actualbiomass = sum(biomass_density))%>% #total biomass seems super low since its been divided by survey.... idk what to do about that
  select(site, Year, total_frass, total_actualbiomass)%>%
  distinct()%>%
  left_join(check_date, by = c("Year", "site"))%>% #divide the totals by the correct number of surveys and stuff from check_date 
  mutate(
    total_frass_divided = total_frass/number_surveys_frass,
    total_actualbiomass_divided = total_actualbiomass/number_surveys_bugs)
#plot new estimations
plot_seasonal_estimates <- function(data, site_name){
  # Filter data for the chosen site
  site_data <- data[data$site == site_name, ]
  # First plot (left axis)
  plot(site_data$Year,
       site_data$total_frass_divided,
       type="b",
       col="sienna",
       ylab="Frass Density Totals",
       xlab="Year",
       main=paste("Site:", site_name))
  
  # Allow second plot on same figure
  par(new=TRUE)
  # Second plot (right axis)
  plot(site_data$Year,
       site_data$total_actualbiomass_divided,
       type="b",
       col="forestgreen",
       axes=FALSE,
       xlab="",
       ylab="")
  axis(4)
  mtext("Biomass Density Totals", side=4, line=3)
  
  legend("topleft",
         legend=c("Frass","Biomass"),
         col=c("sienna","forestgreen"),
         lty=1,
         pch=1, cex=0.8)}
plot_seasonal_estimates(season_totals_impute, 8892356)

#----------------------------------------------------------------------------------------
#plot on linecharts with centroid dates
density_plotting_impute <- function(data, year_choice, site_choice) {
  
  df <- data %>%
    filter(Year == year_choice, site == site_choice)
  
  if (nrow(df) == 0 ||
      all(is.na(df$julianweek)) ||
      all(is.na(df$frass_density))) {
    
    plot.new()
    text(0.5, 0.5,
         paste("No data for", site_choice, year_choice),
         cex = 1.2)
    return(invisible(NULL))
  }
  
  ## ---- Centroids ----
  frass_density_centroid <- weighted.mean(
    df$julianweek, df$frass_density, na.rm = TRUE
  )
  
  biomass_centroid <- weighted.mean(
    df$julianweek, df$biomass_density, na.rm = TRUE
  )
  
  ## ---- subset imputed rows ----
  imputed_frass   <- df %>% filter(imputed_frass == 1)
  imputed_biomass <- df %>% filter(imputed_biomass == 1)
  
  ## ---- Plot ----
  par(mar = c(5, 4, 4, 6))
  
  # LEFT AXIS — FRASS
  plot(
    df$julianweek, df$frass_density,
    type = "l",
    col = "sienna",
    lwd = 2,
    xlab = "Julian week",
    ylab = "Frass density",
    main = paste(site_choice, year_choice)
  )
  
  # red points where frass was imputed
  if (nrow(imputed_frass) > 0) {
    points(imputed_frass$julianweek, imputed_frass$frass_density,
           col = "red", pch = 19, cex = 1.2)
  }
  
  if (is.finite(frass_density_centroid)) {
    abline(v = frass_density_centroid, col = "sienna", lty = 2, lwd = 2)
  }
  
  ## ---- RIGHT AXIS — BIOMASS ----
  par(new = TRUE)
  
  plot(
    df$julianweek, df$biomass_density,
    type = "l",
    col = "forestgreen",
    lwd = 2,
    axes = FALSE,
    xlab = "",
    ylab = "",
    ylim = range(df$biomass_density, na.rm = TRUE)
  )
  
  # red points where biomass was imputed
  if (nrow(imputed_biomass) > 0) {
    points(imputed_biomass$julianweek, imputed_biomass$biomass_density,
           col = "red", pch = 19, cex = 1.2)
  }
  
  axis(side = 4, col.axis = "forestgreen")
  mtext("Actual Caterpillar Biomass Density",
        side = 4, line = 2.0, col = "forestgreen", cex = 0.8)
  
  if (is.finite(biomass_centroid)) {
    abline(v = biomass_centroid, col = "forestgreen", lty = 2, lwd = 2)
  }
  
  ## ---- Legend ----
  legend(
    "topleft",
    legend = c(
      "Frass density",
      "Frass centroid",
      "Caterpillar biomass density",
      "Biomass centroid",
      "Imputed value"
    ),
    col = c("sienna", "sienna", "forestgreen", "forestgreen", "red"),
    lwd = c(2, 2, 2, 2, NA),
    lty = c(1, 2, 1, 2, NA),
    pch = c(NA, NA, NA, NA, 19),
    bty = "n",
    cex = 0.8
  )
  
  invisible(df)
}
density_plotting_impute(imputation_totals, 2025, 8892356)
#saving as a pdf------------------------ ^^^^^
# Years for each site
years_PR   <- c(2015, 2018, 2019, 2021, 2022)
years_NCBG <- setdiff(2015:2025, 2020)   
setwd("C:/Z_School/HurlbertLab/graphs")
#set up pdf
pdf(
  file = "imputation_plotting.pdf",
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
    density_plotting_impute(
      data = imputation_totals,
      year_choice = yr,
      site_choice = 8892356  
    ),
    silent = TRUE)}
for (yr in years_PR) {
  try(
    density_plotting_impute(
      data = imputation_totals,
      year_choice = yr,
      site_choice = 117   
    ),
    silent = TRUE)}
dev.off()



#getting centroid dates for imputated values
imputation_centroid <- imputation_totals_week %>%
  group_by(site, Year) %>%
  summarise(
    frass_centroid_D =
      sum(julianweek * frass_density, na.rm = TRUE) /
      sum(frass_density, na.rm = TRUE),
    
    biomass_centroid_D =
      sum(julianweek * biomass_density, na.rm = TRUE) /
      sum(biomass_density, na.rm = TRUE),
    
    .groups = "drop")

# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#   Re-Doing centroid anomolies on new imputation data
# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+ 
#calculate mean of frass and actual caterpillar biomass centroids
mean_imputation_centroids <- imputation_centroid %>%
  group_by(site)%>%
  mutate(mean_frass =mean(frass_centroid_D),
         mean_actualbiomass = mean(biomass_centroid_D))%>%
  ungroup()%>%
  mutate(
    actualBiomass_diff = biomass_centroid_D - mean_actualbiomass,
    frass_diff = frass_centroid_D - mean_frass) %>%
  select(Year, site, frass_centroid_D, mean_frass, frass_diff, biomass_centroid_D, mean_actualbiomass, actualBiomass_diff)
##plot on 1:1 line so we can see if same years come before and after---------------------------------
shared_years <- c(2015, 2018, 2019, 2021, 2022)
plot_data <- mean_imputation_centroids[mean_imputation_centroids$Year %in% shared_years, ]
years <- sort(unique(plot_data$Year))
cols <- colorRampPalette(c("red","goldenrod1", "springgreen3","dodgerblue", "purple"))(length(years)) #color ramp
year_cols <- cols[match(plot_data$Year, years)]
site_pch <- ifelse(plot_data$site == 117, 8, 17) #symbols
#plot
plot(plot_data$frass_diff,
     plot_data$actualBiomass_diff,
     col = year_cols,
     pch = site_pch,
     cex = 1.5,
     xlab = "Frass centroid difference",
     ylab = "Actual biomass centroid difference")
abline(h = 0, lty = 2)  # horizontal line at y = 0
abline(v = 0, lty = 2)  # vertical line at x = 0
#legend
legend("topleft",
       legend = c("117", "8892356"),
       pch = c(17,8),
       title = "Site",
       bty = "n")
legend("bottomright",
       legend = years,
       col = cols,
       pch = 16,
       title = "Year",
       bty = "n")
#Plot anomaly of cat centroid at given year and site and anomaly of frass centroid at given year and given site and find r^2 and p value
plot_site_correlation <- function(data, site_name){
  
  # Filter site
  site_data <- data[data$site == site_name, ]
  
  # Spearman correlation
  cor_test <- cor.test(site_data$frass_centroid_D,
                       site_data$biomass_centroid_D,
                       method = "spearman")
  
  rho <- cor_test$estimate
  r2 <- rho^2
  
  # Scatter plot
  plot(site_data$frass_centroid_D,
       site_data$biomass_centroid_D,
       pch = 19,
       col = "black",
       xlab = "Frass Density",
       ylab = "Actual Biomass Density",
       main = paste("Site:", site_name))
  
  # Add regression line
  model <- lm(biomass_centroid_D ~ frass_centroid_D, data = site_data)
  abline(model, col = "blue", lwd = 2)
  
  # Add Spearman R² text
  legend("topleft",
         legend = paste("Spearman rho =", round(rho, 3)),
         bty = "n")}

plot_site_correlation(mean_imputation_centroids, "117")

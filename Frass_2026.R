### Explaining seasonal and annual variation in caterpillar occurrence using frass
#R 4.4.2
#Savannah Carter

#Loading in appropriate libraries 
library(gsheet)
library(dplyr)
library(tidyr)

#--------------------------------------------------------------------------------------------------------------------------------
# reading in frassdata and then functions for correcting time and date
#--------------------------------------------------------------------------------------------------------------------------------

# Function for reading in frass data from GoogleDoc
# *if aim is to backup GoogleDoc and write to disk only, then open =F and write = T
# *if aim is to use data without writing to disk, then open = T and write = F
frassData = function(open = F, write = F) {
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

# Function that takes a date field (formatted as %m/%d/%Y) and a time field (hh:mm in 24h time), converts the date to julian day and adds the fractional
#(day represented by the hours and minutes)
julianDayTime = function(date, hour_min) {
  require(lubridate)
  jday = yday(date)
  temp = sapply(strsplit(hour_min, ":"), function(x) {
    x = as.numeric(x)
    x[1] + x[2]/60
  })
  output = jday + temp/24
  return(output)
}
#--------------------------------------------------------------------------------------------------------------------------------
# altering frassdata so that times and days are corrected
#--------------------------------------------------------------------------------------------------------------------------------
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
# using data to find mean frass 
#--------------------------------------------------------------------------------------------------------------------------------
#read in proper url, change dates and label events for below meanfrass
url = "https://docs.google.com/spreadsheets/d/1RwXzwhHUbP0m5gKSOVhnKZbS1C_NrbdfHLglIVCzyFc/edit#gid=1611171427"
events = gsheet2tbl(url)
events$date = as.Date(events$date, format = "%m/%d/%Y")

meanfrass = data %>%
  filter(!is.na(Frass.mass..mg.)) %>%
  mutate(site = as.character(ifelse(Site=="Botanical Garden", 8892356, 117))) %>%
  group_by(site, Date.Collected, Year, jday) %>%
  summarize(mass = mean(frass.mg.d, na.rm=T),
            density = mean(frass.no.d, na.rm=T)) %>%
  left_join(events[, c('date', 'site', 'reliability')], by = c('Date.Collected' = 'date', 
                                                               'site' = 'site')) %>%
  rename(date = Date.Collected)

#--------------------------------------------------------------------------------------------------------------------------------
# function for plotting frass phenology
#--------------------------------------------------------------------------------------------------------------------------------
# minReliability is the minimum reliability score for including in the analysis.
#    3 - reliable, no obvious problems
#    2 - frass traps wet, or potential minor issues
#    1 - major problems, unreliable frass data
frassplot = function(frassdata, inputSite, year, color = 'black', new = T, 
                     var = 'mass', minReliability = 0, xlab = 'Julian day', ylab = '', 
                     jds = c(136, 167, 197), # May 15, Jun 15, Jul 15
                     ...) {
  
  temp = filter(frassdata, site == inputSite, Year == year, reliability >= minReliability) %>%
    data.frame()
  
  if (new) {
    plot(temp$jday, temp[, var], xlab = xlab, ylab = ylab,
         type = 'l', col = color, xaxt = 'n',...)
    points(temp$jday, temp[, var], pch = 16, col = color,...)
    mtext(jds, 1, at = jds, line = 1)
    axis(1, at = c(jds, jds+14), tck = -.02, labels = FALSE)
  } else {
    points(temp$jday, temp[, var], type = 'l', col = color, ...)
    points(temp$jday, temp[, var], pch = 16, col = color, ...)
  }
}
#general frassplot example to alter
frassplot(meanfrass, inputSite = 117, 2022, 'red', new = T, var = 'mass', xlim = c(138,205),
          ylim = c(0, 10.5), lwd = 2, minReliability = 2, xlab = "Julian Day", ylab = "Frass (mg./day)", lty = 'solid', main = 'NCBG Frass')

#+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
#--------------------------------------------------------------------------------------------------------------------------------
# bringing in caterpillar count data
#--------------------------------------------------------------------------------------------------------------------------------
# Reading in Caterpillar Count dataset:
fullDataset = read.csv('data/fullDataset_2025-06-17.csv') #**need to update this to at least end of summer
#filtering caterpillar count dataset to just show NCBG
NCBG = fullDataset %>%
  filter(Name == "NC Botanical Garden", Year == 2025)
#filtering caterpillar count dataset to just show PR
PR = fullDataset %>%
  filter(Name =="Prairie Ridge Ecostation", Year == 2025)

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

#as of now showing 2025 caterpillar count vs density??? need to look into this
# Make sure to establish beatvis.bg to use meanDensityByDay function
beatvis.bg = meanDensityByWeek(NCBG, ordersToInclude = 'caterpillar', plot = TRUE) #just showing 2025 data
beatvis.pr = meanDensityByWeek(PR, ordersToInclude = 'caterpillar', plot = TRUE, new = TRUE)

#+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
#--------------------------------------------------------------------------------------------------------------------------------
#Plots looking at frass vs caterpillar density and biomass 
#--------------------------------------------------------------------------------------------------------------------------------
#linear model of frass mass vs caterpillar density -> compared by jweek
frass_density_lm <- function(fullDataset, data, events,
                             site_name, year) {
  
  # 1. Filter full dataset for caterpillar density
  sitefilter_fulldataset <- fullDataset %>%
    dplyr::filter(Name == site_name, Year == year)
  
  catsfiltered <- meanDensityByWeek(
    sitefilter_fulldataset,
    ordersToInclude = "caterpillar"
  )
  
  # 2. Build mean frass per week
  meanfrass_week <- data %>%
    dplyr::filter(!is.na(Frass.mass..mg.)) %>%
    dplyr::mutate(
      site = as.character(ifelse(Site == "Botanical Garden", 8892356, 117)),
      julianweek = 7 * floor(jday / 7) + 4
    ) %>%
    dplyr::group_by(site, Date.Collected, Year, julianweek) %>%
    dplyr::summarize(
      mass = mean(frass.mg.d, na.rm = TRUE),
      density = mean(frass.no.d, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::left_join(
      events[, c("date", "site", "reliability")],
      by = c("Date.Collected" = "date", "site" = "site")
    ) %>%
    dplyr::rename(date = Date.Collected)
  
  sitefilter_meanfrass <- meanfrass_week %>%
    dplyr::filter(site == 8892356, Year == year)
  
  # 3. Join frass + caterpillar density
  joined_frasscat <- dplyr::left_join(
    catsfiltered,
    sitefilter_meanfrass,
    by = "julianweek"
  )
  
  # 4. Linear regression
  linear_model <- lm(mass ~ meanDensity, data = joined_frasscat) #mean density comes from average density per survey (n surveys)
  
  # 5. Plot
  plot(joined_frasscat$meanDensity, joined_frasscat$mass,
       xlab = "Mean Caterpillar Density",
       ylab = "Mean Frass Mass (mg)",
       main = paste(site_name, year, "Linear Regression: Frass ~ Density"),
       pch = 16, col = "darkgray")
  
  abline(linear_model, col = "red", lwd = 2)
  
  # 6. Legend with equation + R²
  coefs <- coef(linear_model)
  intercept <- round(coefs[1], 3)
  slope <- round(coefs[2], 3)
  r2 <- round(summary(linear_model)$r.squared, 3)
  
  eq <- paste0("y = ", slope, "x + ", intercept, "\nR² = ", r2)
  
  legend("topleft", legend = eq, bty = "n", text.col = "blue", cex = 1)
  
  # 7. Return model + joined data
  return(list(
    model = linear_model,
    data = joined_frasscat))
}

#example of use:
frass_density_lm( fullDataset = fullDataset, data = data, events = events, site_name = "Prairie Ridge Ecostation", year = 2021 )

#linear model of frass mass vs caterpillar biomass -> compared by jweek
frass_biomass_lm <- function(fullDataset, data, events,
                             site_name, year) {
  
  # 1. Filter full dataset for caterpillar biomass
  sitefilter_fulldataset <- fullDataset %>%
    dplyr::filter(Name == site_name, Year == year)
  
  catsfiltered <- meanDensityByWeek(
    sitefilter_fulldataset,
    ordersToInclude = "caterpillar"
  )
  
  # 2. Build mean frass per week
  meanfrass_week <- data %>%
    dplyr::filter(!is.na(Frass.mass..mg.)) %>%
    dplyr::mutate(
      site = as.character(ifelse(Site == "Botanical Garden", 8892356, 117)),
      julianweek = 7 * floor(jday / 7) + 4
    ) %>%
    dplyr::group_by(site, Date.Collected, Year, julianweek) %>%
    dplyr::summarize(
      mass = mean(frass.mg.d, na.rm = TRUE),
      density = mean(frass.no.d, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::left_join(
      events[, c("date", "site", "reliability")],
      by = c("Date.Collected" = "date", "site" = "site")
    ) %>%
    dplyr::rename(date = Date.Collected)
  
  sitefilter_meanfrass <- meanfrass_week %>%
    dplyr::filter(site == 8892356, Year == year)
  
  # 3. Join frass + caterpillar biomass
  joined_frasscat <- dplyr::left_join(
    catsfiltered,
    sitefilter_meanfrass,
    by = "julianweek"
  )
  
  # 4. Linear regression: frass mass ~ mean biomass
  linear_model <- lm(mass ~ meanBiomass, data = joined_frasscat) #mean biomass comes from average biomass per survey (n surveys)
  
  # 5. Plot
  plot(joined_frasscat$meanBiomass, joined_frasscat$mass,
       xlab = "Mean Caterpillar Biomass",
       ylab = "Frass Mass (mg)",
       main = paste(site_name, year, "Linear Regression: Frass ~ Biomass"),
       pch = 16, col = "darkgray")
  
  abline(linear_model, col = "red", lwd = 2)
  
  # 6. Legend with equation + R²
  coefs <- coef(linear_model)
  intercept <- round(coefs[1], 3)
  slope <- round(coefs[2], 3)
  r2 <- round(summary(linear_model)$r.squared, 3)
  
  eq <- paste0("y = ", slope, "x + ", intercept, "\nR² = ", r2)
  
  legend("topleft", legend = eq, bty = "n", text.col = "blue", cex = 1)
  
  # 7. Return model + joined data
  return(list(
    model = linear_model,
    data = joined_frasscat))
}
#example of use
frass_biomass_lm( fullDataset = fullDataset, data = data, events = events, site_name = "Prairie Ridge Ecostation", year = 2022 )















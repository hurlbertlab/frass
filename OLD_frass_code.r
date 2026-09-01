#code dating before 2026: 

# Plotting frass density or mass over time
library(gsheet)
library(dplyr)
library(tidyr)


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



# Function that takes a date field (formatted as %m/%d/%Y) and a time field
# (hh:mm in 24h time), converts the date to julian day and adds the fractional
# day represented by the hours and minutes

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


# Function for plotting frass phenology
#   minReliability is the minimum reliability score for including in the analysis.
#    3 - reliable, no obvious problems
#    2 - frass traps wet, or potential minor issues
#    1 - major problems, unreliable frass data
# 'jds' causes error (object not found) - what is this supposed to be?


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
  
  


#################################################################

# CC dataset
fullDataset = read.csv('data/fullDataset_2025-06-17.csv')

NCBG = fullDataset %>%
  filter(Name == "NC Botanical Garden", Year == 2025)

PR = fullDataset %>%
  filter(Name =="Prairie Ridge Ecostation", Year == 2025)

#################################################################
# Function for substituting values based on a condition using dplyr::mutate
# Modification of dplyr's mutate function that only acts on the rows meeting a condition
mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}

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

# Function for calculating and displaying arthropod phenology by week
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


# Make sure to establish beatvis.bg to use meanDensityByDay function
beatvis.bg = meanDensityByWeek(NCBG, ordersToInclude = 'caterpillar', plot = TRUE)
beatvis.pr = meanDensityByWeek(PR, ordersToInclude = 'caterpillar', plot = TRUE, new = TRUE)



# Get frass data and then get julian days and times
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

# Sampling event data that specify reliability of data on any given date
# (due to storms, etc that may affect frass recovery)

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

write.csv(meanfrass, "data/frass_by_day_2015-2021.csv", row.names = F)


#########################################################################
#plot compiling Bot Garden frass from 2015 through 2025
frassplot(meanfrass, inputSite = 8892356, 2015, 'red', new = T, var = 'mass', xlim = c(138,205),
          ylim = c(0, 10.14), lwd = 2, minReliability = 2, xlab = "Julian Day", ylab = "Frass (mg./day)", lty = 'solid', main = 'NCBG Frass')
frassplot(meanfrass, inputSite = 8892356, 2016, 'green', new = F, var = 'mass', xlim = c(138,205),
          ylim = c(0, 10.14), lwd = 2, minReliability = 2, lty = 'twodash', main = 'NCBG Frass')
frassplot(meanfrass, inputSite = 8892356, 2017, 'orange', new = F, var = 'mass', xlim = c(138,205),
          ylim = c(0, 10.14), lwd = 2, minReliability = 2, lty = 'dotted', main = 'NCBG Frass')
frassplot(meanfrass, inputSite = 8892356, 2018, 'blue', new = F, var = 'mass', xlim = c(138,205),
          ylim = c(0, 10.14), lwd = 2, minReliability = 2, lty = 'dashed', main = 'NCBG Frass')
frassplot(meanfrass, inputSite = 8892356, 2018, 'blueviolet', new = F, var = 'mass', xlim = c(138,205),
          ylim = c(0, 10.14), lwd = 2, minReliability = 2, lty = 'longdash', main = 'NCBG Frass')
frassplot(meanfrass, inputSite = 8892356, 2019, 'darkgreen', new = F, var = 'mass', xlim = c(138,205),
          ylim = c(0, 10.14), lwd = 2, minReliability = 2, lty = 'dotdash', main = 'NCBG Frass')
frassplot(meanfrass, inputSite = 8892356, 2021, 'violet', new = F, var = 'mass', xlim = c(138,205),
          ylim = c(0, 10.14), lwd = 2, minReliability = 2, lty = 'solid', main = 'NCBG Frass')
frassplot(meanfrass, inputSite = 8892356, 2022, 'yellow', new = F, var = 'mass', xlim = c(138,205),
          ylim = c(0, 10.14), lwd = 2, minReliability = 2, lty = 'twodash', main = 'NCBG Frass')
frassplot(meanfrass, inputSite = 8892356, 2023, 'grey', new = F, var = 'mass', xlim = c(138,205),
          ylim = c(0, 10.14), lwd = 2, minReliability = 2, lty = 'dotted', main = 'NCBG Frass')
frassplot(meanfrass, inputSite = 8892356, 2024, 'navy', new = F, var = 'mass', xlim = c(138,205),
          ylim = c(0, 10.14), lwd = 2, minReliability = 2, lty = 'dashed', main = 'NCBG Frass')
frassplot(meanfrass, inputSite = 8892356, 2025, 'black', new = F, var = 'mass', xlim = c(138,205),
          ylim = c(0, 10.14), lwd = 2, minReliability = 2, lty = 'longdash', main = 'NCBG Frass')

#legend to decode graphic
legend(136, 10.2, title = "Survey Year", c("2015", "2016", "2017", "2018", "2019", "2021", "2022", "2023", "2024", "2025"), cex = .7, bty = "n", y.intersp = .8,
       lty=c("solid", "twodash", "dotted", "dashed", "longdash", "dotdash", "solid", "twodash", "dotted", "dashed", "longdash"), col=c("red", "green", "orange", "blue", "blueviolet", "darkgreen", "violet", "yellow", "grey", "navy", "black"), lwd = 2)

#plot compiling Prairie Ridge frass from 2015 through 2022. 2016 and 2017 weren't sampled so not included in figure
frassplot(meanfrass, inputSite = 117, 2015, 'red', new = T, var = 'mass', xlim = c(138,205),
          ylim = c(0, 11.5), lwd = 2, minReliability = 2, xlab = "Julian Day", ylab = "Frass (mg./day)", lty = 'solid', main = 'Prairie Ridge Frass')
frassplot(meanfrass, inputSite = 117, 2018, 'blue', new = F, var = 'mass', xlim = c(138,205),
          ylim = c(0, 11.5), lwd = 2, minReliability = 2, lty = 'dashed', main = 'Prairie Ridge Frass')
frassplot(meanfrass, inputSite = 117, 2019, 'purple', new = F, var = 'mass', xlim = c(138,205),
          ylim = c(0, 11.5), lwd = 2, minReliability = 2, lty = 'longdash', main = 'Prairie Ridge Frass')
frassplot(meanfrass, inputSite = 117, 2021, 'yellow', new = F, var = 'mass', xlim = c(138,205),
          ylim = c(0, 11.5), lwd = 2, minReliability = 2, lty = 'solid', main = 'Prairie Ridge Frass')
frassplot(meanfrass, inputSite = 117, 2022, 'pink', new = F, var = 'mass', xlim = c(138,205),
          ylim = c(0, 11.5), lwd = 2, minReliability = 2, lty = 'twodash', main = 'Prairie Ridge Frass')

#legend to decode graphic
legend(137, 11.6, title = "Survey Year", c("2015", "2016", "2017", "2018", "2019", "2021", "2022"), cex = .7, bty = "y", y.intersp = .8,
       lty=c("solid", "twodash", "dotted", "dashed", "longdash", "solid", "twodash"), col=c("red", "green", "orange", "blue", "purple", "yellow", "pink"))

#plot compiling Prairie Ridge and Bot Garden frass from 2015 through 2022 .Not showing 2016 & 2017 prairie ridge only data due to an error - needs trouble shooting.
frassplot(meanfrass, inputSite = 8892356, 2015, 'violet', new = T, var = 'mass', xlim = c(138,205),
          ylim = c(0, 11.5), lwd = 2, minReliability = 2, xlab = "Julian Day", ylab = "Frass (mg./day)", lty = 'solid', main = 'Prairie Ridge vs. Botanical Garden Frass')
frassplot(meanfrass, inputSite = 8892356, 2016, 'green', new = F, var = 'mass', xlim = c(138,205),
          ylim = c(0, 11.5), lwd = 2, minReliability = 2, lty = 'twodash', main = '')
frassplot(meanfrass, inputSite = 8892356, 2017, 'yellow', new = F, var = 'mass', xlim = c(138,205),
          ylim = c(0, 11.5), lwd = 2, minReliability = 2, lty = 'dotted', main = '')
frassplot(meanfrass, inputSite = 8892356, 2018, 'cyan', new = F, var = 'mass', xlim = c(138,205),
          ylim = c(0, 11.5), lwd = 2, minReliability = 2, lty = 'dashed', main = '')
frassplot(meanfrass, inputSite = 8892356, 2019, 'bisque4', new = F, var = 'mass', xlim = c(138,205),
          ylim = c(0, 11.5), lwd = 2, minReliability = 2, lty = 'longdash', main = '')
frassplot(meanfrass, inputSite = 8892356, 2021, 'orange', new = F, var = 'mass', xlim = c(138,205),
          ylim = c(0, 11.5), lwd = 2, minReliability = 2, lty = 'dotdash', main = '')
frassplot(meanfrass, inputSite = 8892356, 2022, 'red', new = F, var = 'mass', xlim = c(138,205),
          ylim = c(0, 11.5), lwd = 2, minReliability = 2, lty = 'solid', main = '')
frassplot(meanfrass, inputSite = 117, 2015, 'blueviolet', new = F, var = 'mass', xlim = c(138,205),
          ylim = c(0, 11.5), lwd = 2, minReliability = 2, xlab = "Julian Day", ylab = "Frass (mg./day)", lty = 'solid', main = '')
frassplot(meanfrass, inputSite = 117, 2016, 'darkgreen', new = F, var = 'mass', xlim = c(138,205),
          ylim = c(0, 11.5), lwd = 2, minReliability = 2, lty = 'twodash', main = '')
frassplot(meanfrass, inputSite = 117, 2017, 'darkgoldenrod', new = F, var = 'mass', xlim = c(138,205),
          ylim = c(0, 11.5), lwd = 2, minReliability = 2, lty = 'dotted', main = '')
frassplot(meanfrass, inputSite = 117, 2018, 'deepskyblue4', new = F, var = 'mass', xlim = c(138,205),
          ylim = c(0, 11.5), lwd = 2, minReliability = 2, lty = 'dashed', main = '')
frassplot(meanfrass, inputSite = 117, 2019, 'black', new = F, var = 'mass', xlim = c(138,205),
          ylim = c(0, 11.5), lwd = 2, minReliability = 2, lty = 'longdash', main = '')
frassplot(meanfrass, inputSite = 117, 2021, 'darkorange3', new = F, var = 'mass', xlim = c(138,205),
          ylim = c(0, 11.5), lwd = 2, minReliability = 2, lty = 'dotdash', main = '')
frassplot(meanfrass, inputSite = 117, 2022, 'darkred', new = F, var = 'mass', xlim = c(138,205),
          ylim = c(0, 11.5), lwd = 2, minReliability = 2, lty = 'solid', main = '')
#legend to decode graphic
legend("topleft", cex = .58 , title = "Survey Site & Year", c("BG 2015", "BG 2016", "BG 2017", "BG 2018", "BG 2019", "BG 2021", "BG 2022", "PR 2015", "PR 2016", "PR 2017", "PR 2018", "PR 2019", "PR 2021", "PR 2022"), lwd = 2, bty = "n",
       lty=c("solid", "twodash", "dotted", "dashed", "longdash", "dotdash", "solid", "solid", "twodash", "dotted", "dashed", "longdash", "dotdash", "solid"), col=c("violet", "green", "yellow", "cyan", "bisque4", "orange", "red", "blueviolet", "darkgreen", "darkgoldenrod", "deepskyblue4", "black", "darkorange3", "darkred"))


##### savannah's attempt at connecting data to output from frass_imagedataRemastered so that they can be compared
#first need to convert area to volume size by doing to 1.5 power then adding up all particles with same date, then we can join the output df with datafiltered which has mass of frass, then compare caluclated area vs actual mass measured 



DataFiltered = filter(data, Year %in% c(2021:2025))
#combining data on frass pieces mass and particle number from 2021-2025 and the area data 
combined <- left_join(DataFiltered, output, by = c("Trap", "Date.Collected", "Year"))

areafrass = combined
#adding reliability score from events df
events <- events %>%
  rename(Date.Collected = date)  
areafrass <- left_join(
  areafrass,
  events %>% select(Date.Collected, reliability),
  by = c("Date.Collected")
)
#change back events column to be just date
events <- events %>%
  rename(date = Date.Collected)  
#taking inspo from meanfrass and doing the same thing with the area
meanarea <- areafrass %>%
  filter(!is.na(Area)) %>%
  mutate(Site.x = as.character(ifelse(Site.x == "Botanical Garden", 8892356, 117))) %>%
  group_by(Site.x, Date.Collected, Year.x, jday) %>%
  summarize(
    Area = mean(Area, na.rm = TRUE),
    density = mean(frass.no.d, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(events[, c("date", "site", "reliability")], by = c("Date.Collected" = "date", "Site.x" = "site")) %>%
  rename(date = Date.Collected)



  
  
#################################################################################################
#plotting both frass and arthrocount, Future goal: *make this a actual function so that can just plug in year/site
                          

#### linear model of frass mass vs caterpillar density, change name and year to get data ####
#filter the fulldataset so that cat count info can be found 
sitefilter_fulldataset = filter(fullDataset, Name== "Prairie Ridge Ecostation", Year== 2015)
catsfiltered = meanDensityByWeek(sitefilter_fulldataset, ordersToInclude = "caterpillar")
#meanfrass filtered by year, site, and convert jday to julianweek
meanfrass_week <- data %>%
  filter(!is.na(Frass.mass..mg.)) %>%
  mutate(
    site = as.character(ifelse(Site == "Botanical Garden", 8892356, 117)),
    julianweek = 7 * floor(jday / 7) + 4   # midpoint of each 7‑day block
  ) %>%
  group_by(site, Date.Collected, Year, julianweek) %>%
  summarize(
    mass = mean(frass.mg.d, na.rm = TRUE),
    density = mean(frass.no.d, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    events[, c("date", "site", "reliability")],
    by = c("Date.Collected" = "date", "site" = "site")
  ) %>%
  rename(date = Date.Collected)

sitefilter_meanfrass <- meanfrass_week %>%
  filter(site == 8892356, Year == 2015)  

#do a table join with two data sets so LM works
joined_frasscat <- left_join(catsfiltered, sitefilter_meanfrass, by = "julianweek")
#linear regression using catsdensity as indep var and frass as depen var
linear_regeression_frass_catdensity <- lm(mass ~ meanDensity, data = joined_frasscat)
summary(linear_regeression_frass_catdensity)
#Scatterplot of the data
plot(joined_frasscat$meanDensity, joined_frasscat$mass,
     xlab = "Mean Caterpillar Density",
     ylab = "Mean Frass Mass (mg)",
     main = "PR 2015 Linear Regression: Frass ~ Density",
     pch = 16, col = "darkgray")
#Add regression line
abline(linear_regeression_frass_catdensity, col = "red", lwd = 2)
#legend
# Extract coefficients
coefs <- coef(linear_regeression_frass_catdensity)
intercept <- round(coefs[1], 3)
slope <- round(coefs[2], 3)
# Get R-squared
r2 <- summary(linear_regeression_frass_catdensity)$r.squared
r2 <- round(r2, 3)
# Build equation text
eq <- paste0("y = ", slope, "x + ", intercept, "\nR² = ", r2)
# Add legend to plot
legend("topleft", legend = eq, bty = "n", text.col = "blue", cex = 1)

#### same thing but caterpillar meanbiomass vs meanfrass ####
#filter the fulldataset so that cat count info can be found 
sitefilter_fulldataset = filter(fullDataset, Name== "Prairie Ridge Ecostation", Year== 2015)
catsfiltered = meanDensityByWeek(sitefilter_fulldataset, ordersToInclude = "caterpillar")
#meanfrass filtered by year, site, and convert jday to julianweek
meanfrass_week <- data %>%
  filter(!is.na(Frass.mass..mg.)) %>%
  mutate(
    site = as.character(ifelse(Site == "Botanical Garden", 8892356, 117)),
    julianweek = 7 * floor(jday / 7) + 4   # midpoint of each 7‑day block
  ) %>%
  group_by(site, Date.Collected, Year, julianweek) %>%
  summarize(
    mass = mean(frass.mg.d, na.rm = TRUE),
    density = mean(frass.no.d, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    events[, c("date", "site", "reliability")],
    by = c("Date.Collected" = "date", "site" = "site")
  ) %>%
  rename(date = Date.Collected)

sitefilter_meanfrass <- meanfrass_week %>%
  filter(site == 8892356, Year == 2015)  

#do a table join with two data sets so LM works
joined_frasscat <- left_join(catsfiltered, sitefilter_meanfrass, by = "julianweek")
#linear regression using catsdensity as indep var and frass as depen var
linear_regeression_frass_catdensity <- lm(mass ~ meanBiomass, data = joined_frasscat)
summary(linear_regeression_frass_catdensity)
#Scatterplot of data
plot(joined_frasscat$meanBiomass, joined_frasscat$mass,
     xlab = "Mean Caterpillar Biomass",
     ylab = "Frass Mass (mg)",
     main = "PR 2015 Linear Regression: Frass ~ Biomass",
     pch = 16, col = "darkgray")
#Add regression line
abline(linear_regeression_frass_catdensity, col = "red", lwd = 2)
#legend
# Extract coefficients
coefs <- coef(linear_regeression_frass_catdensity)
intercept <- round(coefs[1], 3)
slope <- round(coefs[2], 3)
# Get R-squared
r2 <- summary(linear_regeression_frass_catdensity)$r.squared
r2 <- round(r2, 3)
# Build equation text
eq <- paste0("y = ", slope, "x + ", intercept, "\nR² = ", r2)
# Add legend to plot
legend("topleft", legend = eq, bty = "n", text.col = "blue", cex = 1)


##############################################################################################################
#code about comparing two different methodologies of frass collection:
#########################################################
## Data Setup ##

# Milk jug frass data
# Function for reading in frass data from GoogleDoc
# *if aim is to backup GoogleDoc and write to disk only, then open = F and write = T
# *if aim is to use data without writing to disk, then open = T and write = F

frassLoad = function(open = T, write = F) {
  require(gsheet)
  url = "https://docs.google.com/spreadsheets/d/1RwXzwhHUbP0m5gKSOVhnKZbS1C_NrbdfHLglIVCzyFc/edit#gid=806965256"
  data = gsheet2tbl(url)
  
  if (write) {
    # Write a copy
    write.csv(data, paste('data/frass_', Sys.Date(), '.csv', sep = ''),
              row.names = F)
  }
  if (open) { return (data) }
}

# Filter paper frass data
# Function for reading in frass data from GoogleDoc
# *if aim is to backup GoogleDoc and write to disk only, then open = F and write = T
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

# renaming data sets
data = frassLoad(open = T)
NCBG_PR_frassdata = frassData(open = T) 

# removing outliers in frassLoad
dataWO = data[data$Weight_Raw < 50,]
data_rawpcsWO = data[data$Pieces_Raw < 60,]
data_srtdpcsWO = data[data$Pieces_Sorted < 50,]
data_img_exlc_outlier = data[data$Img_Sorted < 20, ]

library(gsheet)
library(dplyr)
library(tidyr)
library(data.table)

##########################################################

## Plotting and Analysis: Milk Jug Frass Collection ## (7 plots) 

# Weight of Raw Frass / Weight of Sorted Frass
raw_sort_outlier_excl = lm(Weight_Sorted ~ Weight_Raw, data = dataWO)
plot(data$Weight_Raw[data$Weight_Raw<50], data$Weight_Sorted[data$Weight_Raw<50],
     main = "Milk Jug Method:\nFrass Mass Comparison (mg.)", 
     xlab = "Raw", ylab = "Sorted", 
     pch = 20, cex = 1, col = 'deepskyblue2')
abline(raw_sort_outlier_excl, col = "deepskyblue2")
sortraw_sum = summary(raw_sort_outlier_excl)
sortraw_sum_r2 = sortraw_sum$adj.r.squared
mylabel = bquote(italic(R)^2 == .(format(sortraw_sum_r2, digits = 3)))
text(x = 2.5, y = 30.9, labels = mylabel)
abline(a = 0, b = 1, col = "lightgrey")
legend("topleft", cex = .7, title = " ", c("Regression", "1:1"), lwd = 2, bty = "n",
       lty=c("dashed", "solid"), col=c("deepskyblue2", "lightgrey"))

# Sorted Frass Comparison: Pieces vs. Weight (mg.)
plot(data$Pieces_Sorted[data$Pieces_Sorted<100], data$Weight_Sorted[data$Pieces_Sorted<100], 
     main = "Sorted Frass Comparison:\n Pieces vs. Weight", 
     xlab = "Pieces", ylab ="Weight (mg.)", col = 'orange', pch = 20)
sorted_lm = lm(data$Weight_Sorted[data$Pieces_Sorted<100] ~ data$Pieces_Sorted[data$Pieces_Sorted<100], data = data_srtdpcsWO)
abline(sorted_lm)
sorted_lm_sum = summary(sorted_lm)
sorted_lm_r2 = sorted_lm_sum$adj.r.squared
mylabel = bquote(italic(R)^2 == .(format(sorted_lm_r2, digits = 3)))
text(x = 12, y = 30, labels = mylabel)

# % of Area (Unsorted) / Sorted Weight
plot(dataWO$Img_Raw, dataWO$Weight_Sorted, 
     main = "Frass Comparison:\n Image-based Estimate vs. Mass", 
     xlab = "% of Area (Unsorted)",  ylab = "Sorted Weight (mg.)", 
     ylim = c(-.5,20), col = 'orange' , pch = 18)
raw_img = lm(dataWO$Weight_Sorted ~ dataWO$Img_Raw, data = dataWO)
abline(raw_img)
imgwght_sum = summary(raw_img)
imgwght_sum_r2 = imgwght_sum$adj.r.squared
mylabel = bquote(italic(R)^2 == .(format(imgwght_sum_r2, digits = 3)))
text(x = 1.9, y = 19.7, labels = mylabel)

# % of Area Comparison: Unsorted Frass / Sorted Frass
plot(data$Img_Raw[data$Img_Sorted<20], data$Img_Sorted[data$Img_Sorted<20], 
     main = "% of Area Comparison", 
     xlab = "Unsorted Frass", ylab ="Sorted Frass", 
     ylim = c(-.5, 10), col = 'deepskyblue2', pch = 20)
rawsort_img = lm(data$Img_Sorted ~ data$Img_Raw, data = data)
abline(rawsort_img)
rawsort_img_sum = summary(rawsort_img)
img_rawsrt_r2 = rawsort_img_sum$adj.r.squared
mylabel = bquote(italic(R)^2 == .(format(img_rawsrt_r2, digits = 3)))
text(x = 2, y = 9.7, labels = mylabel)

# Sorted Pieces / % of Area (Unsorted) -- to evaluate an image based estimate for mass
plot(data$Pieces_Sorted[data$Pieces_Raw<60], data$Img_Raw[data$Pieces_Raw<60], 
     main = "Frass Comparison:\n # of Pieces vs. % of Area", 
     xlab = "Sorted Frass", ylab = "% of Area (Unsorted)", 
     xlim = c(0,30.009), pch = 20, cex = 1, col = 'orange')
raw_pcs = lm(Img_Raw ~ Pieces_Sorted, data = data)
raw_pcs_outlier_excl = lm(Img_Raw ~ Pieces_Sorted, data = data_rawpcsWO)
abline(raw_pcs_outlier_excl)
sortrawimg_sum = summary(raw_pcs_outlier_excl)
sortrawimg_sum_r2 = sortrawimg_sum$adj.r.squared
mylabel = bquote(italic(R)^2 == .(format(sortrawimg_sum_r2, digits = 3)))
text(x = 4.9, y = 6.1, labels = mylabel)

# Sorted Pieces / % of Area (Sorted) -- to evaluate an image based estimate for mass
plot(data$Pieces_Sorted[data$Pieces_Sorted<50], data$Img_Sorted[data$Pieces_Sorted<50], 
     main = "Frass Comparison:\n # of Pieces vs. % of Area", 
     xlab = "Sorted Frass", ylab = "% of Area (Sorted)",
     pch = 20, cex = 1, col = 'deepskyblue2')
sort_pcs_outlier_excl = lm(Img_Sorted ~ Pieces_Sorted, data = data_srtdpcsWO)
abline(sort_pcs_outlier_excl)
sortdpcs_sum = summary(sort_pcs_outlier_excl)
sortdpcs_sum_r2 = sortdpcs_sum$adj.r.squared
mylabel = bquote(italic(R)^2 == .(format(sortdpcs_sum_r2, digits = 3)))
text(x = 7, y = 7.3, labels = mylabel)

# Volume vs. Mass -- to evaluate an image-based estimate for weight
plot(data$Volume_raw, data$Weight_Sorted, 
     main = "Estimated Volume As A Proxy For Frass Weight", 
     xlab = "Unsorted Frass Volume (mm^3)", 
     ylab = "Sorted Frass Weight (mg.)", ylim = c(.5, 35), pch = 16, 
     col = "orange")
volwght.lm = lm(data$Weight_Sorted ~ data$Volume_raw, data = data)
abline(volwght.lm)
volwght.lm_sum = summary(volwght.lm)
volwght.lm_sum_r2 = volwght.lm_sum $adj.r.squared
mylabel = bquote(italic(R)^2 == .(format(volwght.lm_sum_r2, digits = 3)))
text(x = 37, y = 34, labels = mylabel)

####################################################################

## Methods Comparison: Filter paper vs. Milk jug collection ## 

# Data setup #

# Traps were set on Jul 3rd, Filter paper collected on the 6th & 10th; 
# but milk jug collected on 10th - must sum filter paper frass by circle 
# in order to make an accurate comparison

# Create new data table for filter paper, select for all NCBG rows starting 7/6 
filterfrass_all = NCBG_PR_frassdata[c(1154:1169,1182:1197),]
# Isolate by frass trap site (remove non-milk jug adjacent code)
srtd_filterpaper = filterfrass_all[ ! filterfrass_all$Survey %in% c("1DBD","2DBS", "3DBV","4DCE","5DCI","6DCM","7DCQ","8DCV"), ]
# Sum values of identical survey codes (to account for the additional collection days)
filtermass = aggregate(Frass.mass..mg. ~ Survey, data = srtd_filterpaper, sum)
filterpcs = aggregate(Frass.number ~ Survey, data = srtd_filterpaper, sum)
# Create new data table with these summed values, merging pieces & mass data
filter_sum = merge(filterpcs, filtermass, by = "Survey")
filter_sum$Date.Collected <- "7/10/2018"

# Create data table for normal values and collection dates
filterdates_nonsum = NCBG_PR_frassdata[c(1210:1225, 1238:1253), c("Survey","Frass.mass..mg.","Frass.number", "Date.Collected")]
# Isolate by frass trap site (remove non-milk jug adjacent code)
filter_normal = filterdates_nonsum[ ! filterdates_nonsum$Survey %in% c("1DBD","2DBS", "3DBV","4DCE","5DCI","6DCM","7DCQ","8DCV"), ]
# Combine this "normal" data table with the adjusted table 
filterpaper = rbind(filter_sum, filter_normal)

# Create new table for milk jug data, isolating by frass trap site, get rid of NAs
milkjugs = data[c(73:96),c("Survey","Weight_Sorted", "Pieces_Sorted", "Date.Collected")]
na.omit(milkjugs)

# Merge both data tables to compare milk jug and filter paper by mass and pieces
frasstrapscomp <- filterpaper %>% 
  left_join(milkjugs, by = c("Survey", "Date.Collected"))
setnames(frasstrapscomp, old=c("Weight_Sorted","Pieces_Sorted", "Frass.number","Frass.mass..mg."), new=c("FrassNumber_milkjug", "FrassMass_milkjug","FrassNumber_filterpaper","FrassMass_filterpaper"))
frasstrapscomp = frasstrapscomp[-c(4),]

# Add new columns with mass and pieces data normalized for area
# Area of milk jug = 171.9; Area of filter paper = 433.6 cm^2 
frasstrapscomp = transform(frasstrapscomp, FrassNumber.adj_filterpaper = FrassNumber_filterpaper / 433.6)
frasstrapscomp = transform(frasstrapscomp, FrassMass.adj_filterpaper = FrassMass_filterpaper / 433.6)
frasstrapscomp = transform(frasstrapscomp, FrassNumber.adj_milkjug = FrassNumber_milkjug / 171.9)
frasstrapscomp = transform(frasstrapscomp, FrassMass.adj_milkjug = FrassMass_milkjug / 171.9)


## Plotting and Analysis: Filter Paper vs. Milk Jug  ## (4 plots)

# Filter Paper vs. Milk Jug  - both mass & pieces (non-normalized)
plot(frasstrapscomp$FrassNumber_filterpaper, frasstrapscomp$FrassMass_filterpaper, 
     main = "Frass Collection Method:\nFilter Paper vs. Milk Jug (non-normalized)", 
     xlab = "Pieces", ylab ="Mass (mg.)",  
     col = 'orange', pch = 20, xlim=c(-5, 90), ylim=c(5.5, 95))
points(frasstrapscomp$FrassNumber_milkjug, frasstrapscomp$FrassMass_milkjug, 
       col = 'deepskyblue2', pch = 20)

# Filter Paper vs. Milk Jug  - both mass & pieces (normalized)
par(mar=c(4, 5, 5, 3)) # Bottom, Left, Top, Right
plot(frasstrapscomp$FrassNumber.adj_filterpaper, frasstrapscomp$FrassMass.adj_filterpaper, 
     main = "Frass Collection Method Comparison", 
     xlab = expression(paste("Pieces per ", cm^2)), 
     ylab = expression(paste("Mg. per ", cm^2)),  
     col = 'orange', pch = 20, cex = 1, xlim=c(0, .5), ylim=c(0, .5))
points(frasstrapscomp$FrassNumber.adj_milkjug, frasstrapscomp$FrassMass.adj_milkjug, 
       col = 'deepskyblue2', pch = 20)
# Add regression through both, add 1:1 ref line
filterpaper.lm = lm(FrassMass.adj_filterpaper ~ FrassNumber.adj_filterpaper, data = frasstrapscomp)
milkjug.lm = lm(FrassMass.adj_milkjug ~ FrassNumber.adj_milkjug, data = frasstrapscomp)
abline(filterpaper.lm, col = "orange", lty = "dashed")
abline(milkjug.lm, col = "deepskyblue2", lty = "dashed")
legend("topright", cex = .7, title = " ", c("Filter Paper", "Milk Jug"), lwd = 2, bty = "n",
       lty=c("dashed", "dashed"), col=c("orange", "deepskyblue2"))


# Filter Paper vs. Milk Jug by mass collected 
par(mar=c(4, 5, 5, 3)) # Bottom, Left, Top, Right
plot(frasstrapscomp$FrassMass.adj_filterpaper, frasstrapscomp$FrassMass.adj_milkjug, 
     main = expression(paste("Frass Collection Method Comparison:" ~ "Mass (mg.) per" ~ cm^{2})),
     xlab = expression(paste("Filter Paper ")), 
     ylab = expression(paste("Milk Jug ")),  
     col = 'deepskyblue2', pch = 19, cex = .8, 
     xlim = c(0, .4), ylim=c(0, .5))
methodcompare.lm = lm(FrassMass.adj_milkjug ~ FrassMass.adj_filterpaper, data = frasstrapscomp )
abline(methodcompare.lm, lty = "dashed", col = "deepskyblue2")
methodcompare_sum = summary(methodcompare.lm)
methodcompare_sum_r2 = methodcompare_sum$adj.r.squared
mylabel = bquote(italic(R)^2 == .(format(methodcompare_sum_r2, digits = 3)))
text(x = .36, y = .06, labels = mylabel)
abline(a = 0, b = 1, col = "lightgrey")
legend("bottomright", cex = .7, title = " ", c("Regression", "1:1"), lwd = 2, bty = "n",
       lty=c("dashed", "solid"), col=c("deepskyblue2", "lightgrey"))


# Filter Paper vs. Milk Jug by pieces collected
par(mar=c(4, 5, 5, 3)) # Bottom, Left, Top, Right
plot(frasstrapscomp$FrassNumber.adj_filterpaper, frasstrapscomp$FrassNumber.adj_milkjug, 
     main = expression(paste("Frass Collection Method Comparison:" ~ "Pieces per" ~ cm^{2})), 
     xlab = expression(paste("Filter Paper ")), 
     ylab = expression(paste("Milk Jug ")),  
     col = 'green3', pch = 19, cex = .8, 
     xlim = c(0, .6), ylim=c(0, .21))
methodcompare_num.lm = lm(frasstrapscomp$FrassNumber.adj_milkjug ~ frasstrapscomp$FrassNumber.adj_filterpaper, data = frasstrapscomp )
abline(methodcompare_num.lm, lty = "dashed", col = "green3")
methodcompare_numsum = summary(methodcompare_num.lm)
methodcompare_numsum_r2 = methodcompare_numsum$adj.r.squared
mylabel = bquote(italic(R)^2 == .(format(methodcompare_numsum_r2, digits = 3)))
text(x = .57, y = .02, labels = mylabel)
abline(a = 0, b = 1, col = "lightgrey")
legend("bottomright", cex = .7, title = " ", c("Regression", "1:1"), lwd = 2, bty = "n",
       lty=c("dashed", "solid"), col=c("green3", "lightgrey"))




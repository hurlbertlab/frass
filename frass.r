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
frassplot2 = function(frassdata, inputSite, year, color = 'black', new = T, #julian week
                      var = 'mass', minReliability = 0, xlab = 'Julian Week', ylab = '', 
                      ...) {
  
  temp = filter(frassdata, site == inputSite, Year == year, reliability >= minReliability) %>%
    mutate(julianweek = 7 * floor(jday / 7) + 4) %>%
    data.frame()
  
  if (new) {
    plot(temp$julianweek, temp[, var], xlab = xlab, ylab = ylab,
         type = 'l', col = color, xaxt = 'n', ...)
    points(temp$julianweek, temp[, var], pch = 16, col = color, ...)
    mtext(jds, 1, at = jds, line = 1)
    axis(1, at = c(jds, jds + 14), tck = -0.02, labels = FALSE)
  } else {
    points(temp$julianweek, temp[, var], type = 'l', col = color, ...)
    points(temp$julianweek, temp[, var], pch = 16, col = color, ...)
  }
} 
frassplot2 = function(frassdata, inputSite, year, color = 'black', new = T, #julian day
                      var = 'mass', minReliability = 0, xlab = 'Julian day', ylab = '', 
                      ...) {
  
  temp = filter(frassdata, site == inputSite, Year == year, reliability >= minReliability) %>%
    mutate(julianweek = 7 * floor(jday / 7) + 4) %>%
    data.frame()
  
  if (new) {
    plot(temp$jday, temp[, var], xlab = xlab, ylab = ylab,
         type = 'l', col = color, xaxt = 'n', ...)
    points(temp$jday, temp[, var], pch = 16, col = color, ...)
    mtext(jds, 1, at = jds, line = 1)
    axis(1, at = c(jds, jds + 14), tck = -0.02, labels = FALSE)
  } else {
    points(temp$jday, temp[, var], type = 'l', col = color, ...)
    points(temp$jday, temp[, var], pch = 16, col = color, ...)
  }
}
#wider margins
par(mar = c(4, 4, 4, 4))
#2025NCBG
sitefilter_fulldataset = filter(fullDataset, Name== "NC Botanical Garden", Year== 2025)
catsfiltered = meanDensityByWeek(sitefilter_fulldataset, ordersToInclude = "caterpillar", plot = TRUE, xlim = c(138, 205), main = 'NCBG 2025', xlab = 'Julian Week', ylab = 'Cat Count')
par(new = TRUE)
par(yaxt = "n")
frassplot2(meanfrass, inputSite = 8892356, 2025, 'red', new = T, var = 'mass',  xlim = c(138, 205),
          ylim = c(0, .5), lwd = 2, minReliability = 2, lty = 'dotted', xlab = '', xaxt = "n")
par(yaxt = "s")
par(xaxt = "s")
axis(side = 4, at = seq(0, 0.5, by = 0.1), col.axis = "red")
mtext("Frass mass (mg)", side = 4, line = 2, col = "red", cex = 1, las = 0) #not working for some reason??
#2024NCBG
sitefilter_fulldataset = filter(fullDataset, Name== "NC Botanical Garden", Year== 2024)
catsfiltered = meanDensityByWeek(sitefilter_fulldataset, ordersToInclude = "caterpillar", plot = TRUE, xlim = c(138, 205), main = 'NCBG 2024', xlab = 'Julian Week', ylab = 'Cat Count')
par(new = TRUE)
par(yaxt = "n")
frassplot2(meanfrass, inputSite = 8892356, 2024, 'red', new = T, var = 'mass',  xlim = c(138, 205),
           ylim = c(0, 1), lwd = 2, minReliability = 2, lty = 'dotted', xlab = '', xaxt = "n")
par(yaxt = "s")
par(xaxt = "s")
axis(side = 4, at = seq(0, 1, by = 0.1), col.axis = "red")
mtext("Frass mass (mg)", side = 4, line = 2, col = "red", cex = 1, las = 0) 
#2023NCBG
sitefilter_fulldataset = filter(fullDataset, Name== "NC Botanical Garden", Year== 2023)
catsfiltered = meanDensityByWeek(sitefilter_fulldataset, ordersToInclude = "caterpillar", plot = TRUE, xlim = c(138, 205), main = 'NCBG 2023', xlab = 'Julian Week', ylab = 'Cat Count')
par(new = TRUE)
par(yaxt = "n")
frassplot2(meanfrass, inputSite = 8892356, 2023, 'red', new = T, var = 'mass',  xlim = c(138, 205),
           ylim = c(0, 3), lwd = 2, minReliability = 2, lty = 'dotted', xlab = '', xaxt = "n")
par(yaxt = "s")
par(xaxt = "s")
axis(side = 4, at = seq(0, 3, by = 0.3), col.axis = "red")
mtext("Frass mass (mg)", side = 4, line = 3, col = "darkgreen", cex = 1.2, las = 0) 
#2022NCBG
sitefilter_fulldataset = filter(fullDataset, Name== "NC Botanical Garden", Year== 2022)
catsfiltered = meanDensityByWeek(sitefilter_fulldataset, ordersToInclude = "caterpillar", plot = TRUE, xlim = c(138, 205), main = 'NCBG 2022', xlab = 'Julian Week', ylab = 'Cat Count')
par(new = TRUE)
par(yaxt = "n")
frassplot2(meanfrass, inputSite = 8892356, 2022, 'red', new = T, var = 'mass',  xlim = c(138, 205),
           ylim = c(0, 1.5), lwd = 2, minReliability = 2, lty = 'dotted', xlab = '', xaxt = "n")
par(yaxt = "s")
par(xaxt = "s")
axis(side = 4, at = seq(0, 1.5, by = 0.1), col.axis = "red")
mtext("Frass mass (mg)", side = 4, line = 3, col = "darkgreen", cex = 1.2, las = 0)
#2022PR
sitefilter_fulldataset = filter(fullDataset, Name== "Prairie Ridge Ecostation", Year== 2022)
catsfiltered = meanDensityByWeek(sitefilter_fulldataset, ordersToInclude = "caterpillar", plot = TRUE, xlim = c(138, 205), main = 'PR 2022', xlab = 'Julian Week', ylab = 'Cat Count')
par(new = TRUE)
par(yaxt = "n")
frassplot2(meanfrass, inputSite = 117, 2022, 'red', new = T, var = 'mass',  xlim = c(138, 205),
           ylim = c(0, 3), lwd = 2, minReliability = 2, lty = 'dotted', xlab = '', xaxt = "n")
par(yaxt = "s")
par(xaxt = "s")
axis(side = 4, at = seq(0, 3, by = 0.1), col.axis = "red")
mtext("Frass mass (mg)", side = 4, line = 3, col = "darkgreen", cex = 1.2, las = 0)
#2021NCBG
sitefilter_fulldataset = filter(fullDataset, Name== "NC Botanical Garden", Year== 2021)
catsfiltered = meanDensityByWeek(sitefilter_fulldataset, ordersToInclude = "caterpillar", plot = TRUE, xlim = c(138, 205), main = 'NCBG 2021', xlab = 'Julian Week', ylab = 'Cat Count')
par(new = TRUE)
par(yaxt = "n")
frassplot2(meanfrass, inputSite = 8892356, 2021, 'red', new = T, var = 'mass',  xlim = c(138, 205),
           ylim = c(0, 1.5), lwd = 2, minReliability = 2, lty = 'dotted', xlab = '', xaxt = "n")
par(yaxt = "s")
par(xaxt = "s")
axis(side = 4, at = seq(0, 1.5, by = 0.1), col.axis = "red")
mtext("Frass mass (mg)", side = 4, line = 3, col = "darkgreen", cex = 1.2, las = 0)
#2021PR
sitefilter_fulldataset = filter(fullDataset, Name== "Prairie Ridge Ecostation", Year== 2021)
catsfiltered = meanDensityByWeek(sitefilter_fulldataset, ordersToInclude = "caterpillar", plot = TRUE, xlim = c(138, 205), main = 'PR 2021', xlab = 'Julian Week', ylab = 'Cat Count')
par(new = TRUE)
par(yaxt = "n")
frassplot2(meanfrass, inputSite = 117, 2021, 'red', new = T, var = 'mass',  xlim = c(138, 205),
           ylim = c(0, 16), lwd = 2, minReliability = 2, lty = 'dotted', xlab = '', xaxt = "n")
par(yaxt = "s")
par(xaxt = "s")
axis(side = 4, at = seq(0, 16, by = 1), col.axis = "red")
mtext("Frass mass (mg)", side = 4, line = 3, col = "darkgreen", cex = 1.2, las = 0)
#2019NCBG
sitefilter_fulldataset = filter(fullDataset, Name== "NC Botanical Garden", Year== 2019)
catsfiltered = meanDensityByWeek(sitefilter_fulldataset, ordersToInclude = "caterpillar", plot = TRUE, xlim = c(138, 205), main = 'NCBG 2019', xlab = 'Julian Week', ylab = 'Cat Count')
par(new = TRUE)
par(yaxt = "n")
frassplot2(meanfrass, inputSite = 8892356, 2019, 'red', new = T, var = 'mass',  xlim = c(138, 205),
           ylim = c(0, 3), lwd = 2, minReliability = 2, lty = 'dotted', xlab = '', xaxt = "n")
par(yaxt = "s")
par(xaxt = "s")
axis(side = 4, at = seq(0, 3, by = 0.2), col.axis = "red")
mtext("Frass mass (mg)", side = 4, line = 3, col = "darkgreen", cex = 1.2, las = 0)
#2019PR
sitefilter_fulldataset = filter(fullDataset, Name== "Prairie Ridge Ecostation", Year== 2019)
catsfiltered = meanDensityByWeek(sitefilter_fulldataset, ordersToInclude = "caterpillar", plot = TRUE, xlim = c(138, 205), main = 'PR 2019', xlab = 'Julian Week', ylab = 'Cat Count')
par(new = TRUE)
par(yaxt = "n")
frassplot2(meanfrass, inputSite = 117, 2019, 'red', new = T, var = 'mass',  xlim = c(138, 205),
           ylim = c(0, 5), lwd = 2, minReliability = 2, lty = 'dotted', xlab = '', xaxt = "n")
par(yaxt = "s")
par(xaxt = "s")
axis(side = 4, at = seq(0, 5, by = 0.4), col.axis = "red")
mtext("Frass mass (mg)", side = 4, line = 3, col = "darkgreen", cex = 1.2, las = 0)
#2018NCBG
sitefilter_fulldataset = filter(fullDataset, Name== "NC Botanical Garden", Year== 2018)
catsfiltered = meanDensityByWeek(sitefilter_fulldataset, ordersToInclude = "caterpillar", plot = TRUE, xlim = c(138, 205), main = 'NCBG 2018', xlab = 'Julian Week', ylab = 'Cat Count')
par(new = TRUE)
par(yaxt = "n")
frassplot2(meanfrass, inputSite = 8892356, 2018, 'red', new = T, var = 'mass',  xlim = c(138, 205),
           ylim = c(0, 9), lwd = 2, minReliability = 2, lty = 'dotted', xlab = '', xaxt = "n")
par(yaxt = "s")
par(xaxt = "s")
axis(side = 4, at = seq(0, 9, by = 1), col.axis = "red")
mtext("Frass mass (mg)", side = 4, line = 3, col = "darkgreen", cex = 1.2, las = 0)
#2018PR
sitefilter_fulldataset = filter(fullDataset, Name== "Prairie Ridge Ecostation", Year== 2018)
catsfiltered = meanDensityByWeek(sitefilter_fulldataset, ordersToInclude = "caterpillar", plot = TRUE, xlim = c(138, 205), main = 'PR 2018', xlab = 'Julian Week', ylab = 'Cat Count')
par(new = TRUE)
par(yaxt = "n")
frassplot2(meanfrass, inputSite = 117, 2018, 'red', new = T, var = 'mass',  xlim = c(138, 205),
           ylim = c(0, 12), lwd = 2, minReliability = 2, lty = 'dotted', xlab = '', xaxt = "n")
par(yaxt = "s")
par(xaxt = "s")
axis(side = 4, at = seq(0, 12, by = 1), col.axis = "red") #skips 10 and 12 on axis??
mtext("Frass mass (mg)", side = 4, line = 3, col = "darkgreen", cex = 1.2, las = 0)
#2017NCBG
sitefilter_fulldataset = filter(fullDataset, Name== "NC Botanical Garden", Year== 2017)
catsfiltered = meanDensityByWeek(sitefilter_fulldataset, ordersToInclude = "caterpillar", plot = TRUE, xlim = c(138, 205), main = 'NCBG 2017', xlab = 'Julian Week', ylab = 'Cat Count')
par(new = TRUE)
par(yaxt = "n")
frassplot2(meanfrass, inputSite = 8892356, 2017, 'red', new = T, var = 'mass',  xlim = c(138, 205),
           ylim = c(0, 10.5), lwd = 2, minReliability = 2, lty = 'dotted', xlab = '', xaxt = "n")
par(yaxt = "s")
par(xaxt = "s")
axis(side = 4, at = seq(0, 10.5, by = 1), col.axis = "red")
mtext("Frass mass (mg)", side = 4, line = 3, col = "darkgreen", cex = 1.2, las = 0)
#2016NCBG
sitefilter_fulldataset = filter(fullDataset, Name== "NC Botanical Garden", Year== 2016)
catsfiltered = meanDensityByWeek(sitefilter_fulldataset, ordersToInclude = "caterpillar", plot = TRUE, xlim = c(138, 205), main = 'NCBG 2016', xlab = 'Julian Week', ylab = 'Cat Count')
par(new = TRUE)
par(yaxt = "n")
frassplot2(meanfrass, inputSite = 8892356, 2016, 'red', new = T, var = 'mass',  xlim = c(138, 205),
           ylim = c(0, 3), lwd = 2, minReliability = 2, lty = 'dotted', xlab = '', xaxt = "n")
par(yaxt = "s")
par(xaxt = "s")
axis(side = 4, at = seq(0, 3, by = .2), col.axis = "red")
mtext("Frass mass (mg)", side = 4, line = 3, col = "darkgreen", cex = 1.2, las = 0)
#2015NCBG
sitefilter_fulldataset = filter(fullDataset, Name== "NC Botanical Garden", Year== 2015)
catsfiltered = meanDensityByWeek(sitefilter_fulldataset, ordersToInclude = "caterpillar", plot = TRUE, xlim = c(138, 205), main = 'NCBG 2015', xlab = 'Julian Week', ylab = 'Cat Count')
par(new = TRUE)
par(yaxt = "n")
frassplot2(meanfrass, inputSite = 8892356, 2015, 'red', new = T, var = 'mass',  xlim = c(138, 205),
           ylim = c(0, 3), lwd = 2, minReliability = 2, lty = 'dotted', xlab = '', xaxt = "n")
par(yaxt = "s")
par(xaxt = "s")
axis(side = 4, at = seq(0, 3, by = .2), col.axis = "red")
mtext("Frass mass (mg)", side = 4, line = 3, col = "darkgreen", cex = 1.2, las = 0)
#2015PR
sitefilter_fulldataset = filter(fullDataset, Name== "Prairie Ridge Ecostation", Year== 2015)
catsfiltered = meanDensityByWeek(sitefilter_fulldataset, ordersToInclude = "caterpillar", plot = TRUE, xlim = c(138, 205), main = 'PR 2015', xlab = 'Julian Week', ylab = 'Cat Count')
par(new = TRUE)
par(yaxt = "n")
frassplot2(meanfrass, inputSite = 117, 2015, 'red', new = T, var = 'mass',  xlim = c(138, 205),
           ylim = c(0, 9), lwd = 2, minReliability = 2, lty = 'dotted', xlab = '', xaxt = "n")
par(yaxt = "s")
par(xaxt = "s")
axis(side = 4, at = seq(0, 9, by = 1), col.axis = "red") 
mtext("Frass mass (mg)", side = 4, line = 3, col = "darkgreen", cex = 1.2, las = 0)
#legend to add to graphs
legend("topleft",                            
       legend = c("Cat count", "Frass mass (mg)"),  
       col = c("black", "red"),              
       lty = c(1, 3),                 
       lwd = 2,                               
       bty = "n",
       cex = .75)                             

#### linear model of frass mass vs caterpillar density, change name and year to get data ####
#filter the fulldataset so that cat count info can be found 
sitefilter_fulldataset = filter(fullDataset, Name== "Prairie Ridge Ecostation", Year== 2015)
catsfiltered = meanDensityByWeek(sitefilter_fulldataset, ordersToInclude = "caterpillar")
#meanfrass filtered by year, site, and convert jday to julianweek
sitefilter_meanfrass <- meanfrass %>%
  filter(jday >= jdRange[1], jday <= jdRange[2]) %>%
  mutate(julianweek = 7 * floor(jday / 7) + 4) %>%
  filter(site == 8892356, Year == 2015)  
#do a table join with two data sets so LM works
joined_frasscat <- left_join(catsfiltered, sitefilter_meanfrass, by = "julianweek")
#linear regression using catsdensity as indep var and frass as depen var
linear_regeression_frass_catdensity <- lm(mass ~ meanDensity, data = joined_frasscat)
summary(linear_regeression_frass_catdensity)
#Scatterplot of the data
plot(joined_frasscat$meanDensity, joined_frasscat$mass,
     xlab = "Mean Caterpillar Density",
     ylab = "Frass Mass (mg)",
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
sitefilter_meanfrass <- meanfrass %>%
  filter(jday >= jdRange[1], jday <= jdRange[2]) %>%
  mutate(julianweek = 7 * floor(jday / 7) + 4) %>%
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


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
                     var = 'mass', minReliability = 0, xlab = 'Julian day', ylab = '', ...) {
  
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
  mutate(site = ifelse(Site=="Botanical Garden", 8892356, 117)) %>%
  group_by(site, Date.Collected, Year, jday) %>%
  summarize(mass = mean(frass.mg.d, na.rm=T),
            density = mean(frass.no.d, na.rm=T)) %>%
  left_join(events[, c('date', 'site', 'reliability')], by = c('Date.Collected' = 'date', 
                                                              'site' = 'site')) %>%
  rename(date = Date.Collected)

write.csv(meanfrass, "data/frass_by_day_2015-2021.csv", row.names = F)



beatvis.pr = rbind(beatsheet.pr, amsurvey.pr)
beatvis.bg = rbind(beatsheet.bg, amsurvey.bg)


# Frass plotting

par(mfcol = c(2,1), mar = c(4,4,1,1), mgp = c(2.25, .75, 0))

## Frass Mass
# Bot Garden
frassplot(meanfrass, inputSite = 8892356, 2019, 'red', new = T, var = 'mass', xlim = c(138,205),
          ylim = c(0, 4), lwd = 2, minReliability = 1, lty = 'dotted', main = 'NCBG, 2019')
frassplot(meanfrass, inputSite = 8892356, 2019, 'red', new = F, var = 'mass', 
          lwd = 3, minReliability = 2, lty = 'dashed')
frassplot(meanfrass, inputSite = 8892356, 2019, 'red', new = F, var = 'mass', 
          lwd = 4, minReliability = 3, lty = 'solid')
par(new = T)
bglep15.mass = meanDensityByDay(beatvis.bg, ordersToInclude = "LEPL", inputYear = 2015,
                                 inputSite = 8892356, jdRange = c(138,205), outlierCount = 30,
                                 plot = T, new = T, plotVar = 'meanBiomass',  xlim = c(138,205),
                                 lwd = 4, col = 'blueviolet', yaxt = 'n', ylab = '')
legend("topleft", c('frass', 'LEPL mass'), lwd = 4, col = c('red', 'blueviolet'))


frassplot(meanfrass, inputSite = 8892356, 2019, 'red', new = T, var = 'mass', xlim = c(138,205), 
          ylim = c(0, 4), lwd = 2, minReliability = 1, lty = 'dotted', main = 'NCBG, 2019')
frassplot(meanfrass, inputSite = 8892356, 2019, 'red', new = F, var = 'mass', 
          lwd = 3, minReliability = 2, lty = 'dashed')
frassplot(meanfrass, inputSite = 8892356, 2019, 'red', new = F, var = 'mass', 
          lwd = 4, minReliability = 3, lty = 'solid')
par(new = T)
bglep16.mass = meanDensityByDay(beatvis.bg, ordersToInclude = "LEPL", inputYear = 2016,
                                 inputSite = 8892356, jdRange = c(138,205), outlierCount = 30,
                                 plot = T, new = T, plotVar = 'meanBiomass', xlim = c(138, 205),
                                 lwd = 4, col = 'blueviolet', yaxt = 'n', ylab = '')


frassplot(meanfrass, inputSite = 8892356, 2017, 'red', new = T, var = 'mass', xlim = c(138,205),
          ylim = c(0, 12), lwd = 2, minReliability = 1, lty = 'dotted', main = 'NCBG, 2017')
frassplot(meanfrass, inputSite = 8892356, 2017, 'red', new = F, var = 'mass', 
          lwd = 3, minReliability = 2, lty = 'dashed')
frassplot(meanfrass, inputSite = 8892356, 2017, 'red', new = F, var = 'mass', 
          lwd = 4, minReliability = 3, lty = 'solid')
par(new = T)
bglep17.mass = meanDensityByDay(beatvis.bg, ordersToInclude = "LEPL", inputYear = 2017,
                                 inputSite = 8892356, jdRange = c(138,205), outlierCount = 30,
                                 plot = T, new = T, plotVar = 'meanBiomass', xlim = c(138, 205),
                                 lwd = 4, col = 'blueviolet', yaxt = 'n', ylab = '')
#plot compiling Bot Garden frass from 2015 through 2018
frassplot(meanfrass, inputSite = 8892356, 2015, 'red', new = T, var = 'mass', xlim = c(138,205),
          ylim = c(0, 10.14), lwd = 2, minReliability = 2, xlab = "Julian Day", ylab = "Frass (mg./day)", lty = 'solid', main = 'NCBG Frass')
frassplot(meanfrass, inputSite = 8892356, 2016, 'green', new = F, var = 'mass', xlim = c(138,205),
          ylim = c(0, 10.14), lwd = 2, minReliability = 2, lty = 'twodash', main = 'NCBG Frass')
frassplot(meanfrass, inputSite = 8892356, 2017, 'orange', new = F, var = 'mass', xlim = c(138,205),
          ylim = c(0, 10.14), lwd = 2, minReliability = 2, lty = 'dotted', main = 'NCBG Frass')
frassplot(meanfrass, inputSite = 8892356, 2018, 'blue', new = F, var = 'mass', xlim = c(138,205),
          ylim = c(0, 10.14), lwd = 2, minReliability = 2, lty = 'dashed', main = 'NCBG Frass')
#legend to decode graphic
legend(136, 10.2, title = "Survey Year", c("2015", "2016", "2017", "2018"), cex = .7, bty = "n", y.intersp = .8,
       lty=c("solid", "twodash", "dotted", "dashed"), col=c("red", "green", "orange", "blue"), lwd = 2)

# Prairie Ridge
frassplot(meanfrass, inputSite = 117, 2019, 'red', new = T, var = 'mass', xlim = c(138, 205),
          ylim = c(0, 6), lwd = 2, minReliability = 1, lty = 'dotted', main = 'PR, 2019')
frassplot(meanfrass, inputSite = 117, 2019, 'red', new = F, var = 'mass', 
          lwd = 3, minReliability = 2, lty = 'dashed')
frassplot(meanfrass, inputSite = 117, 2019, 'red', new = F, var = 'mass', 
          lwd = 4, minReliability = 3, lty = 'solid')
par(new=T)
prlep15.mass = meanDensityByDay(beatvis.pr, ordersToInclude = "LEPL", inputYear = 2015,
                                inputSite = 117, jdRange = c(138,205), outlierCount = 30,
                                plot = T, plotVar = 'meanBiomass', xlim = c(138, 205),
                                lwd = 4, col = 'blueviolet', yaxt = 'n', ylab = '')
#plot compiling Prairie Ridge frass from 2015 through 2018. Not showing 2016 & 2017 data due to an error - needs trouble shooting.
frassplot(meanfrass, inputSite = 117, 2015, 'red', new = T, var = 'mass', xlim = c(138,205),
          ylim = c(0, 11.5), lwd = 2, minReliability = 2, xlab = "Julian Day", ylab = "Frass (mg./day)", lty = 'solid', main = 'Prairie Ridge Frass')
frassplot(meanfrass, inputSite = 117, 2016, 'green', new = F, var = 'mass', xlim = c(138,205),
          ylim = c(0, 11.5), lwd = 2, minReliability = 2, lty = 'twodash', main = 'Prairie Ridge Frass')
frassplot(meanfrass, inputSite = 117, 2017, 'orange', new = F, var = 'mass', xlim = c(138,205),
          ylim = c(0, 11.5), lwd = 2, minReliability = 2, lty = 'dotted', main = 'Prairie Ridge Frass')
frassplot(meanfrass, inputSite = 117, 2018, 'blue', new = F, var = 'mass', xlim = c(138,205),
          ylim = c(0, 11.5), lwd = 2, minReliability = 2, lty = 'dashed', main = 'Prairie Ridge Frass')
#legend to decode graphic
legend(137, 11.6, title = "Survey Year", c("2015", "2016", "2017", "2018"), cex = .7, bty = "y", y.intersp = .8,
       lty=c("solid", "twodash", "dotted", "dashed"), col=c("red", "green", "orange", "blue"))


#plot compiling Prairie Ridge and Bot Garden frass from 2015 & 2018  .Not showing 2016 & 2017 data due to an error - needs trouble shooting.
frassplot(meanfrass, inputSite = 8892356, 2015, 'violet', new = T, var = 'mass', xlim = c(138,205),
          ylim = c(0, 11.5), lwd = 2, minReliability = 2, xlab = "Julian Day", ylab = "Frass (mg./day)", lty = 'solid', main = 'Prairie Ridge vs. Botanical Garden Frass')
frassplot(meanfrass, inputSite = 8892356, 2018, 'blue', new = F, var = 'mass', xlim = c(138,205),
          ylim = c(0, 11.5), lwd = 2, minReliability = 2, lty = 'dashed', main = '')
frassplot(meanfrass, inputSite = 117, 2015, 'green', new = F, var = 'mass', xlim = c(138,205),
          ylim = c(0, 11.5), lwd = 2, minReliability = 2, xlab = "Julian Day", ylab = "Frass (mg./day)", lty = 'twodash', main = '')
frassplot(meanfrass, inputSite = 117, 2018, 'orange', new = F, var = 'mass', xlim = c(138,205),
          ylim = c(0, 11.5), lwd = 2, minReliability = 2, lty = 'dotted', main = '')
#legend to decode graphic
legend("topleft", cex = .58, title = "Survey Site & Year", c("BG 2015", "BG 2018", "PR 2015", "PR 2018"), lwd = 2, bty = "n",
       lty=c("solid", "dashed", "twodash", "dotted"), col=c("violet", "blue", "green", "orange"))


## Frass Density
# Bot Garden
frassplot(meanfrass, inputSite = 8892356, 2015, 'red', new = T, var = 'density', xlim = c(138, 205),
          ylim = c(0, 8), lwd = 2, minReliability = 1, lty = 'dotted', main = 'NCBG, 2015')
frassplot(meanfrass, inputSite = 8892356, 2015, 'red', new = F, var = 'density', 
          lwd = 3, minReliability = 2, lty = 'dashed')
frassplot(meanfrass, inputSite = 8892356, 2015, 'red', new = F, var = 'density', 
          lwd = 4, minReliability = 3, lty = 'solid')
par(new = T)
bglep15.den = meanDensityByDay(beatvis.bg, ordersToInclude = "LEPL", inputYear = 2015,
                                inputSite = 8892356, jdRange = c(138,205), outlierCount = 30,
                                plot = T, new = T, plotVar = 'meanDensity', xlim = c(138, 205),
                                lwd = 4, col = 'blueviolet', yaxt = 'n', ylab = '')
legend("topleft", c('frass', 'LEPL density'), lwd = 4, col = c('red', 'blueviolet'))


frassplot(meanfrass, inputSite = 8892356, 2016, 'red', new = T, var = 'density', xlim = c(138, 205), 
          ylim = c(0, 8), lwd = 2, minReliability = 1, lty = 'dotted', main = 'NCBG, 2016')
frassplot(meanfrass, inputSite = 8892356, 2016, 'red', new = F, var = 'density', 
          lwd = 3, minReliability = 2, lty = 'dashed')
frassplot(meanfrass, inputSite = 8892356, 2016, 'red', new = F, var = 'density', 
          lwd = 4, minReliability = 3, lty = 'solid')
par(new = T)
bglep16.den = meanDensityByDay(beatvis.bg, ordersToInclude = "LEPL", inputYear = 2016,
                               inputSite = 8892356, jdRange = c(138,205), outlierCount = 30,
                               plot = T, new = T, plotVar = 'meanDensity', xlim = c(138, 205),
                               lwd = 4, col = 'blueviolet', yaxt = 'n', ylab = '')

frassplot(meanfrass, inputSite = 8892356, 2017, 'red', new = T, var = 'density', xlim = c(138, 205), 
          ylim = c(0, 10), lwd = 2, minReliability = 1, lty = 'dotted', main = 'NCBG, 2017')
frassplot(meanfrass, inputSite = 8892356, 2017, 'red', new = F, var = 'density', 
          lwd = 3, minReliability = 2, lty = 'dashed')
frassplot(meanfrass, inputSite = 8892356, 2017, 'red', new = F, var = 'density', 
          lwd = 4, minReliability = 3, lty = 'solid')
par(new = T)
bglep17.den = meanDensityByDay(beatvis.bg, ordersToInclude = "LEPL", inputYear = 2017,
                               inputSite = 8892356, jdRange = c(138,205), outlierCount = 30,
                               plot = T, new = T, plotVar = 'meanBiomass', xlim = c(138, 205),
                               lwd = 4, col = 'blueviolet', yaxt = 'n', ylab = '')

# Prairie Ridge
frassplot(meanfrass, inputSite = 117, 2015, 'red', new = T, var = 'mass',  xlim = c(138, 205),
          ylim = c(0, 8), lwd = 2, minReliability = 1, lty = 'dotted', main = 'PR, 2015')
frassplot(meanfrass, inputSite = 117, 2015, 'red', new = F, var = 'mass', 
          lwd = 3, minReliability = 2, lty = 'dashed')
frassplot(meanfrass, inputSite = 117, 2015, 'red', new = F, var = 'mass', 
          lwd = 4, minReliability = 3, lty = 'solid')
par(new = T)
prlep15.den = meanDensityByDay(beatvis.pr, ordersToInclude = "LEPL", inputYear = 2015,
                               inputSite = 117, jdRange = c(138,205), outlierCount = 30,
                               plot = T, new = T, plotVar = 'meanDensity', xlim = c(138, 205),
                               lwd = 4, col = 'blueviolet', yaxt = 'n', ylab = '')



# Plot comparing frass to overall (BS + Vis) lep occurrence at Prairie Ridge, 2015
pdf('output/plots/paper_plots/frass_v_caterpillars.pdf', height = 5, width = 10)
par(mfrow = c(1,2), mar = c(4, 4, 2.5, 4), mgp = c(2.5, .75, 0))

# Prairie Ridge
frassplot(meanfrass, inputSite = 117, 2015, 'red', new = T, var = 'mass', lwd = 4, 
          minReliability = 3, xlim = c(138, 205), ylim = c(1, 7), las = 1, 
          cex.lab = 1.5, ylab = 'Frass (mg / trap / day')

prlep15.fr = filter(meanfrass, site == 117, Year == 2015, jday >= 138, jday <=197)

par(new = T)
prlep15.occ = meanDensityByDay(beatvis.pr, ordersToInclude = "LEPL", inputYear = 2015,
                               inputSite = 117, jdRange = c(138,205), outlierCount = 30,
                               plot = T, new = T, plotVar = 'fracSurveys', xlim = c(138, 205), ylim = c(0, 18),
                               lwd = 4, col = 'blueviolet', yaxt = 'n', ylab = '', xaxt = 'n', xlab = '')
axis(4, at = seq(0, 16, by = 4), las = 1)
axis(4, at = seq(2, 18, by = 4), tcl = -.25, labels = F)
mtext("A", 3, adj = -.2, line = .75, cex = 2)

legend("topleft", c('frass', 'caterpillars'), lwd = 4, col = c('red', 'blueviolet'))


# Bot Garden
frassplot(meanfrass, inputSite = 8892356, 2015, 'red', new = T, var = 'mass', lwd = 4, 
          minReliability = 3, xlim = c(138, 205), ylim = c(0, 3), las = 1, cex.lab = 1.5, ylab = "")

par(new = T)
bglep15.occ = meanDensityByDay(beatvis.bg, ordersToInclude = "LEPL", inputYear = 2015,
                               inputSite = 8892356, jdRange = c(138,205), outlierCount = 30,
                               plot = T, new = T, plotVar = 'fracSurveys', xlim = c(138, 205), ylim = c(0, 24),
                               lwd = 4, col = 'blueviolet', yaxt = 'n', ylab = '', xaxt = 'n', xlab = '')
axis(4, at = seq(0, 24, by = 4), las = 1)
mtext("Caterpillars (% of surveys)", 4, line = 2.5, cex = 1.5)

mtext("B", 3, adj = -.2, line = .75, cex = 2)

dev.off()


# eliminate last 2 survey dates which are part of a late summer peak
prlep15.occ2 = filter(prlep15.occ, julianday <= 197)

occGfit = fitG(prlep15.occ2$julianday, prlep15.occ2$fracSurveys, 
              weighted.mean(prlep15.occ2$julianday, prlep15.occ2$fracSurveys),
              14, 200)  

#lines(138:200, occGfit$par[3]*dnorm(138:200, occGfit$par[1], occGfit$par[2]), col = 'blueviolet', lwd = 2, lty = 'dotted')



par(new=T)
frassplot(meanfrass, "Prairie Ridge", 2015, 'green', new = T, var = 'density', lwd = 7, 
          xlim = c(138, 200), yaxt = "n", minReliability = )

par(new=T)
frassplot(meanfrass, "Prairie Ridge", 2015, 'green', new = T, var = 'mass', lwd = 5, 
          xlim = c(138, 200), yaxt = "n", lty = 'dashed')


prlep15.frden = filter(meanfrass, Site == "Prairie Ridge", Year == 2015, jday >= 138, jday <=197)

frGfit.m = fitG(prlep15.frden$jday, prlep15.frden$density, 
              weighted.mean(prlep15.frden$jday, prlep15.frden$density),
              14, 200)  

lines(138:200, frGfit$par[3]*dnorm(138:200, frGfit$par[1], frGfit$par[2]), col = 'green', lwd = 2)


#convert frass image analysis results file name to a date
titleToDate = function(results) { #results is a file address for results of image analysis
  paste(substr(basename(results), 1, 4), substr(basename(results), 5, 6), substr(basename(results), 7, 8), sep = "-")
}

#convert frass image analysis results file name to a site name
titleToSite = function(results) { #results is a file address for results of image analysis
  paste(ifelse(test = grepl('ncbg', basename(results)), yes = 'Botanical Garden', no = 'Prairie Ridge'), sep = "")
}

#convert frass image analysis results file name to a trap number
titleToTrap = function(results) { #results is a file address for results of image analysis
  paste(ifelse(test = grepl('ncbg', basename(results)), yes = substr(basename(results), 15, 16), no = 
          ifelse(test = grepl(11, basename(results)) | grepl(12, basename(results)), yes = substr(basename(results), 13, 14), no = substr(basename(results), 13, 13))))
}

#extract total volume from an image analysis results file
extractVolTotal = function(frassFile) {
  temp <- read.csv(frassFile)
  sum(temp$Area)
}

#compile total frass volumes to single file with date, site, and trap info
frassVolumes = function(address) { #address is a file address for the folder containing results files
  frassVol <<- tibble(
    sapply(list.files(address), titleToSite),
    sapply(list.files(address), titleToDate),
    sapply(list.files(address), titleToTrap),
    sapply(list.files(address, full.names = TRUE), extractVolTotal))
  names(frassVol) <<- c('Site', 'Date.Collected', 'Trap', 'Frass.Volume.sqmm')
}
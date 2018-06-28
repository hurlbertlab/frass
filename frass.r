# Plotting frass density or mass over time
library(gsheet)

# Get frass data and then get julian days and times
data = frassData(open = T) %>%
  filter(!is.na(Time.Set) & !is.na(Time.Collected)) %>%
  mutate(Date.Set = as.Date(Date.Set, format = "%m/%d/%Y"),
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

write.csv(meanfrass, "data/arthropods/frass_by_day_2015-2017.csv", row.names = F)



beatvis.pr = rbind(beatsheet.pr, amsurvey.pr)
beatvis.bg = rbind(beatsheet.bg, amsurvey.bg)


# Frass plotting

par(mfcol = c(4,2), mar = c(4,4,1,1), mgp = c(2.25, .75, 0))

## Frass Mass
# Bot Garden
frassplot(meanfrass, inputSite = 8892356, 2015, 'red', new = T, var = 'mass', xlim = c(138,205),
          ylim = c(0, 4), lwd = 2, minReliability = 1, lty = 'dotted', main = 'NCBG, 2015')
frassplot(meanfrass, inputSite = 8892356, 2015, 'red', new = F, var = 'mass', 
          lwd = 3, minReliability = 2, lty = 'dashed')
frassplot(meanfrass, inputSite = 8892356, 2015, 'red', new = F, var = 'mass', 
          lwd = 4, minReliability = 3, lty = 'solid')
par(new = T)
bglep15.mass = meanDensityByDay(beatvis.bg, ordersToInclude = "LEPL", inputYear = 2015,
                                 inputSite = 8892356, jdRange = c(138,205), outlierCount = 30,
                                 plot = T, new = T, plotVar = 'meanBiomass',  xlim = c(138,205),
                                 lwd = 4, col = 'blueviolet', yaxt = 'n', ylab = '')
legend("topleft", c('frass', 'LEPL mass'), lwd = 4, col = c('red', 'blueviolet'))


frassplot(meanfrass, inputSite = 8892356, 2016, 'red', new = T, var = 'mass', xlim = c(138,205), 
          ylim = c(0, 4), lwd = 2, minReliability = 1, lty = 'dotted', main = 'NCBG, 2016')
frassplot(meanfrass, inputSite = 8892356, 2016, 'red', new = F, var = 'mass', 
          lwd = 3, minReliability = 2, lty = 'dashed')
frassplot(meanfrass, inputSite = 8892356, 2016, 'red', new = F, var = 'mass', 
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


# Prairie Ridge
frassplot(meanfrass, inputSite = 117, 2015, 'red', new = T, var = 'mass', xlim = c(138, 205),
          ylim = c(0, 6), lwd = 2, minReliability = 1, lty = 'dotted', main = 'PR, 2015')
frassplot(meanfrass, inputSite = 117, 2015, 'red', new = F, var = 'mass', 
          lwd = 3, minReliability = 2, lty = 'dashed')
frassplot(meanfrass, inputSite = 117, 2015, 'red', new = F, var = 'mass', 
          lwd = 4, minReliability = 3, lty = 'solid')
par(new=T)
prlep15.mass = meanDensityByDay(beatvis.pr, ordersToInclude = "LEPL", inputYear = 2015,
                                inputSite = 117, jdRange = c(138,205), outlierCount = 30,
                                plot = T, plotVar = 'meanBiomass', xlim = c(138, 205),
                                lwd = 4, col = 'blueviolet', yaxt = 'n', ylab = '')


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



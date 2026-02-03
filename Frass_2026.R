### Explaining seasonal and annual variation in caterpillar occurrence using frass
#R 4.4.2
#Savannah Carter

#Loading in appropriate libraries 
library(gsheet)
library(dplyr)
library(tidyr)
library(purrr)
library(rstatix)

#--------------------------------------------------------------------------------------------------------------------------------
# reading in frassdata and then functions for correcting time and date
#--------------------------------------------------------------------------------------------------------------------------------

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
#filtering caterpillar count dataset to just show NCBG and Prairie Ridge for selected years
NCBG_PR <- fullDataset %>%
  filter(Name %in% c("NC Botanical Garden", "Prairie Ridge Ecostation"),
    Year %in% 2015:2025)
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
beatvis.bg = meanDensityByWeek(NCBG_PR, ordersToInclude = 'caterpillar', plot = TRUE) #just showing 2025 data
beatvis.pr = meanDensityByWeek(PR, ordersToInclude = 'caterpillar', plot = TRUE, new = TRUE)

#+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
#--------------------------------------------------------------------------------------------------------------------------------
#Plots looking at frass mass vs caterpillar ABUNDANCE (line chart)
#--------------------------------------------------------------------------------------------------------------------------------
#having mean Frass data sorted by julian week
meanfrass <- meanfrass %>%
  mutate(julianweek = 7 * floor(jday / 7) + 4)
#make sure year is a integar
meanfrass <- meanfrass %>%
  mutate(Year = as.integer(Year))
#combine and average meanfrass data that comes from same week (2015-2023 this happened)
meanfrass_combinedweeks <- meanfrass %>%
  group_by(site, Year, julianweek) %>%
  summarise(
    # average frass measurements
    mass = mean(mass, na.rm = TRUE),
    density = mean(density, na.rm = TRUE),
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

#combining code above to plot cat abundance vs frassmass
plot_frass_vs_cats <- function(site, year,
                               frass_data = meanfrass_combinedweeks,
                               survey_data = fullDataset) {
  
  ## ---- Site mapping ----
  if (site == "NCBG") {
    frass_site <- "8892356"
    survey_site <- "NC Botanical Garden"
  } else if (site == "PR") {
    frass_site <- "117"
    survey_site <- "Prairie Ridge Ecostation"
  } else {
    stop("site must be 'NCBG' or 'PR'")
  }
  
  ## ---- Frass data ----
  frass <- frass_data %>%
    mutate(
      julianweek = 7 * floor(jday / 7) + 4,
      Year = as.integer(Year)
    ) %>%
    filter(
      site == frass_site,
      Year == year
    ) %>%
    rename(
      frass_mass = mass,
      frass_density = density,
      frass_reliability = reliability
    )
  
  ## ---- Caterpillar data (add Year back explicitly) ----
  cats <- survey_data %>%
    filter(
      Name == survey_site,
      Year == year
    ) %>%
    group_by(Year) %>%
    group_split() %>%
    purrr::map_dfr(~ {
      out <- meanDensityByWeek(
        surveyData = .x,
        ordersToInclude = "caterpillar",
        allDates = TRUE
      )
      out$Year <- unique(.x$Year)
      out
    })
  
  ## ---- Join by julianweek + Year ----
  dat <- frass %>%
    left_join(cats, by = c("julianweek", "Year"))
  
  ## ---- Plot ----
  par(mar = c(5, 4, 4, 4) + 0.1)
  
  plot(
    dat$julianweek,
    dat$frass_mass,
    type = "l",
    col = "forestgreen",
    lwd = 2,
    xlab = "Julian week",
    ylab = "Frass mass"
  )
  
  par(new = TRUE)
  
  plot(
    dat$julianweek,
    dat$totalCount,
    type = "l",
    col = "sienna",
    lwd = 2,
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
  
  axis(side = 4, cex.axis = 0.8)
  mtext("Caterpillar abundance", side = 4, line = 2, cex = 1.0)
  
  title(
    main = paste(year, survey_site)
  )
  
  legend(
    "topleft",
    legend = c("Frass mass", "Caterpillar abundance"),
    col = c("forestgreen", "sienna"),
    lwd = 1.5,
    cex =0.8,
    bty = "n"
  )
  
  invisible(dat)
}

plot_frass_vs_cats("PR", 2022)


#+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
#--------------------------------------------------------------------------------------------------------------------------------
#Plots looking at frass mass vs caterpillar BIOMASS (total vs mean) (line charts)
#--------------------------------------------------------------------------------------------------------------------------------
## builds off work done before this so make sure to run all the stuff that comes before frass vs cat abundace graphs to make this work

#combining code above to plot cat totalBIOMASS vs frassmass
plot_frass_vs_catstotBiomass <- function(site, year,
                               frass_data = meanfrass_combinedweeks,
                               survey_data = fullDataset) {
  
  ## ---- Site mapping ----
  if (site == "NCBG") {
    frass_site <- "8892356"
    survey_site <- "NC Botanical Garden"
  } else if (site == "PR") {
    frass_site <- "117"
    survey_site <- "Prairie Ridge Ecostation"
  } else {
    stop("site must be 'NCBG' or 'PR'")
  }
  
  ## ---- Frass data ----
  frass <- frass_data %>%
    mutate(
      julianweek = 7 * floor(jday / 7) + 4,
      Year = as.integer(Year)
    ) %>%
    filter(
      site == frass_site,
      Year == year
    ) %>%
    rename(
      frass_mass = mass,
      frass_density = density,
      frass_reliability = reliability
    )
  
  ## ---- Caterpillar data (add Year back explicitly) ----
  cats <- survey_data %>%
    filter(
      Name == survey_site,
      Year == year
    ) %>%
    group_by(Year) %>%
    group_split() %>%
    purrr::map_dfr(~ {
      out <- meanDensityByWeek(
        surveyData = .x,
        ordersToInclude = "caterpillar",
        allDates = TRUE
      )
      out$Year <- unique(.x$Year)
      out
    })
  
  ## ---- Join by julianweek + Year ----
  dat <- frass %>%
    left_join(cats, by = c("julianweek", "Year"))
  
  ## ---- Plot ----
  par(mar = c(5, 4, 4, 4) + 0.1)
  
  plot(
    dat$julianweek,
    dat$frass_mass,
    type = "l",
    col = "forestgreen",
    lwd = 2,
    xlab = "Julian week",
    ylab = "Frass mass"
  )
  
  par(new = TRUE)
  
  plot(
    dat$julianweek,
    dat$totalBiomass,
    type = "l",
    col = "sienna",
    lwd = 2,
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
  
  axis(side = 4, cex.axis = 0.8)
  mtext("Total Caterpillar Biomass", side = 4, line = 2, cex = 1.0)
  
  title(
    main = paste(year, survey_site)
  )
  
  legend(
    "topleft",
    legend = c("Frass mass", "Total cat biomass"),
    col = c("forestgreen", "sienna"),
    lwd = 1.5,
    cex =0.8,
    bty = "n"
  )
  
  invisible(dat)
}

plot_frass_vs_catstotBiomass("PR", 2022)
plot_frass_vs_catsmeanBiomass("PR", 2022)

####same thing but for meanbiomass (means divided by #surveys).... very similar to total biomass graphs
#combining code above to plot mean cat BIOMASS vs frassmass
plot_frass_vs_catsmeanBiomass <- function(site, year,
                                      frass_data = meanfrass_combinedweeks,
                                      survey_data = fullDataset) {
  
  ## ---- Site mapping ----
  if (site == "NCBG") {
    frass_site <- "8892356"
    survey_site <- "NC Botanical Garden"
  } else if (site == "PR") {
    frass_site <- "117"
    survey_site <- "Prairie Ridge Ecostation"
  } else {
    stop("site must be 'NCBG' or 'PR'")
  }
  
  ## ---- Frass data ----
  frass <- frass_data %>%
    mutate(
      julianweek = 7 * floor(jday / 7) + 4,
      Year = as.integer(Year)
    ) %>%
    filter(
      site == frass_site,
      Year == year
    ) %>%
    rename(
      frass_mass = mass,
      frass_density = density,
      frass_reliability = reliability
    )
  
  ## ---- Caterpillar data (add Year back explicitly) ----
  cats <- survey_data %>%
    filter(
      Name == survey_site,
      Year == year
    ) %>%
    group_by(Year) %>%
    group_split() %>%
    purrr::map_dfr(~ {
      out <- meanDensityByWeek(
        surveyData = .x,
        ordersToInclude = "caterpillar",
        allDates = TRUE
      )
      out$Year <- unique(.x$Year)
      out
    })
  
  ## ---- Join by julianweek + Year ----
  dat <- frass %>%
    left_join(cats, by = c("julianweek", "Year"))
  
  ## ---- Plot ----
  par(mar = c(5, 4, 4, 4) + 0.1)
  
  plot(
    dat$julianweek,
    dat$frass_mass,
    type = "l",
    col = "forestgreen",
    lwd = 2,
    xlab = "Julian week",
    ylab = "Frass mass"
  )
  
  par(new = TRUE)
  
  plot(
    dat$julianweek,
    dat$meanBiomass,
    type = "l",
    col = "sienna",
    lwd = 2,
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
  
  axis(side = 4, cex.axis = 0.8)
  mtext("Mean Caterpillar Biomass", side = 4, line = 2, cex = 1.0)
  
  title(
    main = paste(year, survey_site)
  )
  
  legend(
    "topleft",
    legend = c("Frass mass", "Mean cat biomass"),
    col = c("forestgreen", "sienna"),
    lwd = 1.5,
    cex =0.8,
    bty = "n"
  )
  
  invisible(dat)
}
plot_frass_vs_catsmeanBiomass("NCBG", 2015)


#+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
#--------------------------------------------------------------------------------------------------------------------------------
#Determining Normality and variance 
#--------------------------------------------------------------------------------------------------------------------------------
#Shapiro Wilkes test NCBG
normality_tests_NCBG <- cats_NCBG %>%
  select(site, Year, meanBiomass, meanDensity, frass_mass) %>%
  pivot_longer(
    cols = c(meanBiomass, meanDensity, frass_mass),
    names_to = "measure",   
    values_to = "value"
  ) %>%
  group_by(site, Year, measure) %>%
  shapiro_test(value)
#shapiro wilkes test PR
normality_tests_PR <- cats_PR %>%
  select(site, Year, meanBiomass, meanDensity, frass_mass) %>%
  pivot_longer(
    cols = c(meanBiomass, meanDensity, frass_mass),
    names_to = "measure",   
    values_to = "value"
  ) %>%
  group_by(site, Year, measure) %>%
  shapiro_test(value)


#+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
#--------------------------------------------------------------------------------------------------------------------------------
#Plots looking at frass vs caterpillar biomass and density (spearman correlation)
#--------------------------------------------------------------------------------------------------------------------------------
#linear model of frass mass vs caterpillar biomass -> compared by jweek
plot_spearman_frass_catbiomass <- function(site, year,
                                           frass_data = meanfrass_combinedweeks,
                                           survey_data = fullDataset) {
  
  ## ---- Site mapping ----
  if (site == "NCBG") {
    frass_site <- "8892356"
    survey_site <- "NC Botanical Garden"
  } else if (site == "PR") {
    frass_site <- "117"
    survey_site <- "Prairie Ridge Ecostation"
  } else {
    stop("site must be 'NCBG' or 'PR'")
  }
  
  ## ---- Frass data ----
  frass <- frass_data %>%
    mutate(
      julianweek = 7 * floor(jday / 7) + 4,
      Year = as.integer(Year)
    ) %>%
    filter(
      site == frass_site,
      Year == year
    ) %>%
    rename(
      frass_mass = mass
    )
  
  ## ---- Caterpillar data ----
  cats <- survey_data %>%
    filter(
      Name == survey_site,
      Year == year
    ) %>%
    group_by(Year) %>%
    group_split() %>%
    purrr::map_dfr(~ {
      out <- meanDensityByWeek(
        surveyData = .x,
        ordersToInclude = "caterpillar",
        allDates = TRUE
      )
      out$Year <- unique(.x$Year)
      out
    })
  
  ## ---- Join by julianweek + Year ----
  dat <- frass %>%
    left_join(cats, by = c("julianweek", "Year")) %>%
    filter(
      !is.na(frass_mass),
      !is.na(meanBiomass)
    )
  
  ## ---- Spearman correlation ----
  cor_test <- cor.test(
    dat$frass_mass,
    dat$meanBiomass,
    method = "spearman",
    exact = FALSE   # safer for ties
  )
  rho_val <- round(cor_test$estimate, 2)
  p_val <- signif(cor_test$p.value, 2)
  
  ## ---- Plot ----
  par(mar = c(5, 5, 4, 2) + 0.1)
  
  plot(
    dat$frass_mass,
    dat$meanBiomass,
    pch = 16,
    col = "sienna",
    xlab = "Mean weekly frass mass",
    ylab = "Mean caterpillar biomass",
    main = paste(year, survey_site)
  )
  ## Optional: monotonic trend (loess is okay for visualization)
  lines(
    lowess(dat$frass_mass, dat$meanBiomass),
    col = "forestgreen",
    lwd = 2
  )
  ## correlation text
  legend(
    "topleft",
    legend = c(
      paste0("\u03C1 = ", rho_val),   # rho symbol
      paste0("p = ", p_val)
    ),
    bty = "n",
    cex = 0.9
  )
  invisible(dat)
}

#example of use:
plot_spearman_frass_catbiomass("NCBG", 2015)

#linear model of frass mass vs caterpillar density -> compared by jweek
plot_spearman_frass_catdensity <- function(site, year,
                                           frass_data = meanfrass_combinedweeks,
                                           survey_data = fullDataset) {
  
  ## ---- Site mapping ----
  if (site == "NCBG") {
    frass_site <- "8892356"
    survey_site <- "NC Botanical Garden"
  } else if (site == "PR") {
    frass_site <- "117"
    survey_site <- "Prairie Ridge Ecostation"
  } else {
    stop("site must be 'NCBG' or 'PR'")
  }
  
  ## ---- Frass data ----
  frass <- frass_data %>%
    mutate(
      julianweek = 7 * floor(jday / 7) + 4,
      Year = as.integer(Year)
    ) %>%
    filter(
      site == frass_site,
      Year == year
    ) %>%
    rename(
      frass_mass = mass
    )
  
  ## ---- Caterpillar data ----
  cats <- survey_data %>%
    filter(
      Name == survey_site,
      Year == year
    ) %>%
    group_by(Year) %>%
    group_split() %>%
    purrr::map_dfr(~ {
      out <- meanDensityByWeek(
        surveyData = .x,
        ordersToInclude = "caterpillar",
        allDates = TRUE
      )
      out$Year <- unique(.x$Year)
      out
    })
  
  ## ---- Join by julianweek + Year ----
  dat <- frass %>%
    left_join(cats, by = c("julianweek", "Year")) %>%
    filter(
      !is.na(frass_mass),
      !is.na(meanDensity)
    )
  
  ## ---- Spearman correlation ----
  cor_test <- cor.test(
    dat$frass_mass,
    dat$meanDensity,
    method = "spearman",
    exact = FALSE   # safer for ties
  )
  
  rho_val <- round(cor_test$estimate, 2)
  p_val <- signif(cor_test$p.value, 2)
  
  ## ---- Plot ----
  par(mar = c(5, 5, 4, 2) + 0.1)
  
  plot(
    dat$frass_mass,
    dat$meanBiomass,
    pch = 16,
    col = "sienna",
    xlab = "Mean weekly frass mass",
    ylab = "Mean caterpillar Density",
    ylim = c(0, 2.5),
    main = paste(year, survey_site)
  )
  
  ## Optional: monotonic trend (loess is okay for visualization)
  lines(
    lowess(dat$frass_mass, dat$meanDensity),
    col = "forestgreen",
    lwd = 2
  )
  
  ## correlation text
  legend(
    "topleft",
    legend = c(
      paste0("\u03C1 = ", rho_val),   # rho symbol
      paste0("p = ", p_val)
    ),
    bty = "n",
    cex = 0.9
  )
  
  invisible(dat)
}

#example of use:
plot_spearman_frass_catdensity("PR", 2022)





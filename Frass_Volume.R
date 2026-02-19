#Reading in Frass pellet size data, calculating biomass based on size
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

# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#   Datasets needed:
# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#reading in all folders with frass particle data and combining into one dataframe
yearsWithData = 2021:2025

frassPath = "//ad.unc.edu/bio/HurlbertLab/Databases/CaterpillarsCount/Frass"


output = data.frame(Year = NULL, Site = NULL, Trap = NULL, Date = NULL, Particle = NULL, Area = NULL)


for (year in yearsWithData) {
  
  tmpPath = paste0(frassPath, "/", year, "/Results")
  
  filelist <- list.files(path = tmpPath, recursive = TRUE,
                         pattern = "\\.(txt|csv)$", 
                         full.names = TRUE)
  
  # For loop, read in each file
  
  
  for (file in filelist) {
    
    if(year == 2021) {
      
      tmpfile = read.csv(file, header = T)
      
    } else {
      
      tmpfile = read.table(file, sep = '\t', header = T)
      
    }
    
    
    # extracting site, date, and trap from the filename using word()
    # example: "//ad.unc.edu/bio/HurlbertLab/Databases/CaterpillarsCount/Frass/2021/Results/20210730_pr_12.txt"  
    #fix date (8th slash leads to date)
    
    file2 = file %>% 
      str_extract("Results/(.*)$") %>%
      str_remove("Results/") 
    
    Datestring <- word(file2, sep = "_", 1) 
    Site <- word(file2, sep = "_", 2) 
    Trap <- word(file2, sep = "_", 3) %>% word(sep = "\\.", 1)
    
    Date = paste(substr(Datestring, 1, 4), substr(Datestring, 5, 6), substr(Datestring, 7, 8), sep = "-")
    
    
    tmpdf = tmpfile[, 1:2]
    tmpdf$Area = tmpfile$Area
    tmpdf$Year = rep(year, nrow(tmpdf))
    tmpdf$Date = rep(Date, nrow(tmpdf))
    tmpdf$Site = rep(Site, nrow(tmpdf))
    tmpdf$Trap = rep(Trap, nrow(tmpdf))
    
    names(tmpdf)[1] = "Particle"
    
    output = rbind(output, tmpdf)
    
  } # end loop
  
  
}

output = output[, c("Year", "Site", "Trap", "Date", "Particle", "Area")]
output$Date = as.Date(output$Date, format = "%Y-%m-%d")
names(output)[names(output) == "Date"] <- "Date.Collected"

#----------------------------------------------------------------------------------------
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
#   getting Tinbergen biomass dataframe so I can combine with particle stuff
# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#correcting functions so can combine properly----------------------------
#having mean Frass data sorted by julian week
meanfrass <- meanfrass %>%
  mutate(julianweek = 7 * floor(jday / 7) + 4) %>%
  mutate(Year = as.integer(Year)) #make sure year is a integar
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
#   Creating one df that has tinbergen biomass, non corrected biomass, and weekly temps
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


# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#   combining output area df with cats_tinbergen_biomass so that all data in one df
# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#thinking about combining with frass mass OK column to make sure only reliable frass in there...
#alter df so they can be combined
data <- data %>%
  mutate(Trap = tolower(Trap)) %>%
  mutate(Year = as.integer(Year))
area_data <- output %>%
  mutate(Trap = tolower(Trap)) %>%
  mutate(Site = tolower(Site)) %>%
  mutate(Site = case_when(
    Site == "ncbg" ~ "Botanical Garden",
    Site == "pr" ~ "Prairie Ridge"  )) 
#combined dataframes: only reliable measurments in there now
combined_area_frass = area_data %>%
  left_join(
    data %>%
      select(Trap, Site, Date.Collected, Year, OK,),
    by = c("Date.Collected", "Trap", "Site", "Year")) %>%
  filter(OK == 1)

#use the area of the particle to find volume (spheres assumed)
combined_area_frass <- combined_area_frass %>%
  mutate(volume_particles = (4/3)*(Area^(3/2))/sqrt(pi)) %>%
  mutate(jday = as.numeric(format(Date.Collected, "%j"))) %>%
  mutate(julianweek = 7 * floor(jday / 7) + 4)
#find weekly volume
weekly_volume <- combined_area_frass %>%
  group_by(Year, Site, julianweek) %>%
  summarise(
    avg_volume_per_particle = mean(volume_particles, na.rm = TRUE),
    n_particles = n(),
    .groups = "drop") %>%
  rename(site = Site)%>%
  mutate(site = case_when(
    site == "Botanical Garden" ~ 8892356,
    site == "Prairie Ridge" ~ 117 )) %>%
  mutate(site = as.character(site))
#left join with cats_tinbergen_biomass so all data in one table
volume_all_data <- weekly_volume %>%
  left_join(cats_tinbergen_biomass, by = c("julianweek", "Year", "site"))
 
#plotting volume and frass mass and volume on same plot
plot_volume <- function(data, year, site) {
  
  df <- data %>%
    dplyr::filter(Year == year, .data$site == site)
  
  if(nrow(df) == 0) {
    stop("No data for the specified year and site.")
  }
  
  df <- df[order(df$julianweek), ]
  
  mass_x_centroid <- weighted.mean(df$julianweek, df$mass, na.rm = TRUE)
  volume_x_centroid <- weighted.mean(df$julianweek, df$avg_volume_per_particle, na.rm = TRUE)
  
  par(mar = c(5, 4, 4, 4) + 0.1)
  
  plot(
    df$julianweek,
    df$mass,
    type = "l",
    col = "forestgreen",
    lwd = 2,
    xlab = "Julian week",
    ylab = "frass Mass",
    main = paste("Site:", site, "Year:", year)
  )
  
  abline(v = mass_x_centroid, col = "forestgreen", lty = 2, lwd = 2)
  
  par(new = TRUE)
  plot(
    df$julianweek,
    df$avg_volume_per_particle,
    type = "l",
    col = "sienna",
    lwd = 2,
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
  
  axis(side = 4)
  mtext("Volume", side = 4, line = 2)
  
  abline(v = volume_x_centroid, col = "sienna", lty = 2, lwd = 2)
  
  legend(
    "topleft",
    legend = c(
      "frass Mass",
      "frass centroid",
      "volume",
      "volume centroid"
    ),
    col = c("forestgreen", "forestgreen", "sienna", "sienna"),
    lwd = 2,
    lty = c(1, 2, 1, 2),
    bty = "n"
  )
  
  invisible(df)
}

plot_volume(volume_all_data, year = 2022, site = "8892356")









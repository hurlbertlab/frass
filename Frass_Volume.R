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

# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#   Combining date frames
# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#thinking about combining with frass mass OK column to make sure only reliable frass in there...
#lower trap ID in both dataframes
data <- data %>%
  mutate(Trap = tolower(Trap))
output <- output %>%
  mutate(Trap = tolower(Trap))
#make data same type Year as output
data <- data %>%
  mutate(Year = as.integer(Year))
#Lower output site names
output <- output %>%
  mutate(Site = tolower(Site))
#transform PR into Prairie Ridge and NCBG to Botanical Garden
output <- output %>%
  mutate(Site = case_when(
    Site == "ncbg" ~ "Botanical Garden",
    Site == "pr" ~ "Prairie Ridge"  )) 
#combined dataframes: only reliable measurments in there now
combined_area_frass = output %>%
  left_join(
    data %>%
      select(Trap, Site, Date.Collected, Year, OK, Frass.mass..mg., Frass.number, frass.mg.d, frass.no.d),
    by = c("Date.Collected", "Trap", "Site", "Year")) %>%
  filter(OK == 1)


# final data frame has columns for Site, Trap, Date, Particle, Area, and Estimated Volume
FrassAreaVolume = cbind(tmp2,cleaned_output) %>%
  select(-'Filename')
FrassAreaVolume$new <- FrassAreaVolume$Area^1.5
colnames(FrassAreaVolume) = c('Date','Site','Trap', 'Particle', 'Area','estVolume')


























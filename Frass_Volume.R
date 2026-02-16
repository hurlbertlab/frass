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

#

























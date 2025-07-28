# Script for collating all of the image analysis data files into a single table
# with columns for site, date (collection), trap#, particle size
# different years have different ways of being stored in the Frass folder so will be broken down by year to address these differences
# years 2015-2017 didn't do Image analysis, years 2018 and 2019 are confusing so will come back to

#2015-2017: no data with frass
#2018-2019: Frass data accumulated into one results table but no sorting of particles was done before so unable to tell what particles are frass vs other material
#2021: results in .csv form and lowercase lettering used 20210713_ncbg_1a, prairie ridge also surveyed
#2022: results in .txt form and uppercase form 20220612_NCBG_1a, prairie ridge also surveyed
#2023-2025: results in .txt form, uppercase used, only NCBG surveyed 


## want table of combine data to look like:
#Site Date      Trap Particle Area
#NCBG 20250605  1a    1        .8
#NCBG 20250605  1a    2        1.23
#PR   20190605  2b    1        .23
#PR   20190605  2b    2        1.34

library(stringr)
library(tidyverse)
library(readr)
library(data.table)
library(dplyr)
library(lubridate)
library(readxl)

# works with every year! :) ########################################

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


# write as a file


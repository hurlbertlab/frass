# Script for collating all of the image analysis data files into a single table
# with columns for site, date (collection), trap#, particle size
# different years have different ways of being stored in the Frass folder so will be broken down by year to address these differences
# years 2015-2017 didn't do Image analysis, years 2018 and 2019 are confusing so will come back to

#2015-2017: no data with frass
#2018-2019: Frass data accumulated into one results table but no sorting of particles was done before so unable to tell what particles are frass vs other material
#2021: results in .csv form and lowercase lettering used 20210713_ncbg_1a, prairie ridge also surveyed
#2022: results in .txt form and uppercase form 20220612_NCBG_1A, prairie ridge also surveyed
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

## stuff below will work with years 2023-2025 because formatting of names is the same
# error for 2021/2022:: breaks down in line "Error in `[.data.frame`(tmpfile, , c("X", "Area")) : undefined columns selected" because formating different for 2021 and 2022

yearsWithData = c(2023:2025)

frassPath = "//ad.unc.edu/bio/HurlbertLab/Databases/CaterpillarsCount/Frass"

output = data.frame(Year = NULL, Site = NULL, Trap = NULL, Date = NULL, X = NULL, Area = NULL)


for (year in yearsWithData) {
  
  tmpPath = paste0(frassPath, "/", year, "/Results")
  
  filelist <- list.files(path = tmpPath, recursive = TRUE,
                              pattern = "\\.(txt|csv)$", 
                              full.names = TRUE)
  
  # For loop, read in each file
  
  
  for (file in filelist) {
    
    tmpfile = read.table(file, sep = '\t', header = T)
    
    # extracting site, date, and trap from the filename using word()
    # example:
    
    Datestring <- word(file, sep = "_", 1) 
    Site <- word(file, sep = "_", 2) 
    Trap <- word(file, sep = "_", 3) %>% word(sep = "\\.", 1)
    
    Date = paste(substr(Datestring, 1, 4), substr(Datestring, 5, 6), substr(Datestring, 7, 8), sep = "-")
    
    tmpdf = data.frame(Year = NULL, Site = NULL, Trap = NULL, Date = NULL, X = NULL, Area = NULL)
    
    tmpdf = tmpfile[, c("X", "Area")]
    tmpdf$Year = rep(year, nrow(tmpdf))
    tmpdf$Date = rep(Date, nrow(tmpdf))
    tmpdf$Site = rep(Site, nrow(tmpdf))
    tmpdf$Trap = rep(Trap, nrow(tmpdf))
    
    output = rbind(output, tmpdf)
    
  } # end loop
  

}

output = output[, c("Year", "Site", "Trap", "Date", "X", "Area")]
output$Date = as.Date(output$Date, format = "%Y-%m-%d")
names(output)[4] = "Particle"



## 2021/2022 ##


for (file in filelist){
  
  # Add a Volume Column from Area being raised to the 3/2
  Volume = output[5] ^ (3/2)
}
output = cbind(output, new_col = Volume)
names(output)[6] = "Volume"

df2023 <- rbindlist(sapply(list_of_files, fread, simplify = FALSE),
                 use.names = TRUE, idcol = "FileName", fill = TRUE)

# delete columns
df2023 = subset(df2023, select = -c(7:11) )

# rename columns 
colnames(df2023)[2] ="Number of Frass"

# create site and date columns by slicing name
# new dataframe for only the file name
list_of_names <- data.table(df2023$FileName)

# extracting site, date, and trap from the filename using word()
# example:

Date1 <- word(list_of_names, sep = "_", 1) 
Site1 <- word(list_of_names, sep = "_", 2) 
Trap1 <- word(list_of_names, sep = "_", 3)

# subset Trap and Date by removing unnecessary characters
# loop it for all files
length <- row(list_of_names)


for(i in 1:ncol(list_of_names)) {                                             # Loop over character vector
  Date <- word(Date1, sep = "/", 2)
}  

for(i in 1:ncol(list_of_names)) {                                             # Loop over character vector
  Trap <- word(Trap1, sep = ".txt", 1)
}    


## 2024 ##





###########################################################################

# 2021 # 
setwd("//ad.unc.edu/bio/hurlbertlab/Databases/CaterpillarsCount/Frass/2021/Results")


filelist = list.files()

list_of_files <- list.files(path = ".", recursive = TRUE,
                            pattern = "\\.txt$", 
                            full.names = TRUE)

  


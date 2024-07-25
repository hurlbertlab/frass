# Script for reading in ImageJ results files and getting them together in one dataframe

# 1 Create an empty dataframe to store results 

# 2 For loop to read in each file

# 3 add the results of that file to the dataframe

# 4 make sure to add information about site, date, trap, etc. (use strsplit())
library(stringr)
library(dplyr)
options(scipen=999) #disables scientific notation --> aesthetic preference 

# Nested For loop to read in files and add them to a data frame
#initate a blank output dataframe for use in forloop
output = data.frame(x = NULL)

frassyearsfolder= "z:/Databases/CaterpillarsCount/Frass"
formattedfrassyears = c("2022", "2023", "2024")

for (i in 1:length(formattedfrassyears)) {
  frassfolder = paste(frassyearsfolder, "/", formattedfrassyears[i], "/", "Results", sep="")

  for (f in list.files(frassfolder)) {
    tmp = read.table(paste(frassfolder, "/", f, sep  = ""), sep = "\t", header = T)
    tmp$Filename = f
    output = rbind(output, tmp)
  }
}
# remove columns we don't care about (Mean, Min, Max)

cleaned_output = output %>% 
  select(-one_of('Mean', 'Max', 'Min'))

# create new columns for Date, Site, Trap and fill values using strsplit

tmp2 = do.call(rbind, strsplit(as.character(cleaned_output$Filename), "_")) %>%
  as.data.frame()
colnames(tmp2) <- c("Date", "Site", "Trap")
tmp2$Trap = str_extract(tmp2$Trap, "[^.]+")
tmp2$Date = as.Date(tmp2$Date, format = "%Y%m%d")

# final data frame has columns for Site, Trap, Date, Particle, Area, and Estimated Volume
FrassAreaVolume = cbind(tmp2,cleaned_output) %>%
  select(-'Filename')
FrassAreaVolume$new <- FrassAreaVolume$Area^1.5
colnames(FrassAreaVolume) = c('Date','Site','Trap', 'Particle', 'Area','estVolume')

# join with frass gSheet
read_in_data <- gsheet2tbl('https://docs.google.com/spreadsheets/d/1RwXzwhHUbP0m5gKSOVhnKZbS1C_NrbdfHLglIVCzyFc/edit#gid=1479231778')
read_in_data$Site <- ifelse(read_in_data$Site == "Botanical Garden", "NCBG", "PR")
read_in_data$Trap <- tolower(read_in_data$Trap)

summarizedDf <- group_by(FrassAreaVolume, Date, Site, Trap) %>%
  summarize(Area = sum(Area),
            Volume = sum(estVolume)) %>%
  mutate(Trap = tolower(Trap))
read_in_data$Date.Collected = as.Date(read_in_data$Date.Collected, format = "%m/%d/%Y")
FrassAreaVolumeMass = left_join(summarizedDf, read_in_data, by = c("Site", "Date"="Date.Collected", "Trap"))


  
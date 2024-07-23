# Script for reading in ImageJ results files and getting them together in one dataframe

# 1 Create an empty dataframe to store results 

# 2 For loop to read in each file

# 3 add the results of that file to the dataframe

# 4 make sure to add information about site, date, trap, etc. (use strsplit())
library(stringr)
library(dplyr)

# example for loop (1-dimensional output)

output2 = c()

blah = c(3, 5, 2, 0, 12)

for (i in 1:10) {
  
  tmp = 3*i
  
  output2 = c(output2, tmp)

}


# Nested for loop to read in files and add them to a data frame
output = data.frame(Filename = NULL,
                    Particle = NULL,
                    Area = NULL)

frassfolder = "z:/Databases/CaterpillarsCount/Frass/2024/Results"
for (f in list.files(frassfolder)) {
  tmp = read.table(paste(frassfolder, "/", f, sep  = ""), sep = "\t", header = T)
  tmp$Filename = f
  output = rbind(output, tmp)
  }

# remove columns we don't care about (Mean, Min, Max)

cleaned_output = output %>% 
  select(-one_of('Mean', 'Max', 'Min'))

# create new columns for Date, Site, Trap and fill values using strsplit

tmp2 = do.call(rbind, strsplit(as.character(cleaned_output$Filename), "_")) %>%
  as.data.frame()
colnames(tmp2) <- c("Date", "Site", "Trap")
tmp2$Trap = str_extract(tmp2$Trap, "[^.]+")

# final dataframe has columns for Site, Trap, Date, Particle, Area
finaldf = cbind(tmp2,cleaned_output) %>%
  select(-'Filename')
 colnames(finaldf) <- c('Date','Site','Trap', 'Particle', 'Area')
 
 
 
 
 
 
  
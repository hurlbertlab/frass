# Script for analyzing milk jug frass data
frassLoad = function(open = T, write = F) {
  require(gsheet)
  url = "https://docs.google.com/spreadsheets/d/1RwXzwhHUbP0m5gKSOVhnKZbS1C_NrbdfHLglIVCzyFc/edit#gid=806965256"
  data = gsheet2tbl(url)
  
  if (write) {
    # Write a copy
    write.csv(data, paste('data/frass_', Sys.Date(), '.csv', sep = ''),
              row.names = F)
  }
  if (open) { return (data) }
}

# Script for analyzing filter paper data
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

NCBG_PR_frassdata = frassData(open = T)

#removing outliers in frassLoad
data = frassLoad(open = T)
dataWO = data[data$Weight_Raw < 50,]
data_rawpcsWO = data[data$Pieces_Raw < 60,]
data_srtdpcsWO = data[data$Pieces_Sorted < 50,]
data_img_exlc_outlier = data[data$Img_Sorted < 20, ]

# Linear model & plot describing weight_sorted vs weight_raw
raw_sort = lm(Weight_Sorted ~ Weight_Raw, data = data)
raw_sort_outlier_excl = lm(Weight_Sorted ~ Weight_Raw, data = dataWO)
sort_img = lm(Img_Sorted ~ Weight_Sorted, data = dataWO)
plot(data$Weight_Raw, data$Weight_Sorted, main = "Frass Weight Comparison (mg.)", xlab = "Weight Raw", ylab = "Weight Sorted", pch = 17, cex = 1, col = 'red')
plot(data$Weight_Raw[data$Weight_Raw<50], data$Weight_Sorted[data$Weight_Raw<50],main = "Frass Weight Comparison (mg.)", xlab = "Weight Raw", ylab = "Weight Sorted", pch = 17, cex = 1, col = 'red')
abline(raw_sort_outlier_excl)

## Linear models & plots showing raw/sorted img against raw/sorted pcs to describe how much sorting changes % of area 
# Raw 
plot(data$Pieces_Sorted[data$Pieces_Raw<60], data$Img_Raw[data$Pieces_Raw<60], main = "Raw Frass Comparison: # of Pieces vs. % of Area", 
     xlab = "Total Pieces", ylab = "% of Area (unsorted)", pch = 20, cex = 1, col = 'orange')
raw_pcs = lm(Img_Raw ~ Pieces_Sorted, data = data)
raw_pcs_outlier_excl = lm(Img_Raw ~ Pieces_Sorted, data = data_rawpcsWO)
abline(raw_pcs_outlier_excl)
summary(raw_pcs_outlier_excl)

# Sorted
plot(data$Pieces_Sorted[data$Pieces_Sorted<50], data$Img_Sorted[data$Pieces_Sorted<50], main = "Sorted Frass Comparison: # of Pieces vs. % of Area", 
     xlab = "Pieces Sorted", ylab = "% of Area", pch = 20, cex = 1, col = 'blue')
sort_pcs_outlier_excl = lm(Img_Sorted ~ Pieces_Sorted, data = data_srtdpcsWO)
abline(sort_pcs_outlier_excl)
summary(sort_pcs_outlier_excl)

plot(data$Pieces_Sorted[data$Pieces_Sorted<50], data$Img_Sorted[data$Pieces_Sorted<50], main = "Sorted Frass Comparison: # of Pieces vs. % of Area", 
     xlab = "Pieces Sorted", ylab = "% of Area", pch = 20, cex = 1, col = 'blue')
sort_pcs_outlier_excl = lm(Img_Sorted ~ Pieces_Sorted, data = data_srtdpcsWO)
abline(sort_pcs_outlier_excl)
summary(sort_pcs_outlier_excl)

# Raw image versus weight sorted
# Excluded outlier
plot(dataWO$Img_Raw, dataWO$Weight_Sorted, main = "Raw Image v. Weight of sorted Frass", xlab = "Raw Image Percent Coverage", 
     ylab = "Weight Sorted (mg)", col = 'orange' , pch = 18)
raw_img = lm(dataWO$Weight_Sorted ~ dataWO$Img_Raw, data = dataWO)
abline(raw_img)
summary(raw_img)

# Img_raw vs.Img_sort
plot(data$Img_Raw[data$Img_Sorted<20], data$Img_Sorted[data$Img_Sorted<20], main = "Comparison Img_Raw vs. Img_Sorted (% of area estimate)", xlab = "Raw Img.", ylab ="Sorted Img.", col = 'violet', pch = 20)
rawsort_img = lm(data$Img_Sorted ~ data$Img_Raw, data = data)
abline(rawsort_img)
rawsort_img_sum = summary(rawsort_img)
r2 = rawsort_img_sum$adj.r.squared
rawsort_img_pval = rawsort_img_sum$coefficients[2,4]
mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
text(x = 11.5, y = 16, labels = mylabel)

#NCBG Comparison: Filter paper vs. Milk jug collection. Traps set on 3rd.
#Filter paper collected on the 6th & 10th
#milk jug collected on 10th
#must sum filter paper frass per circle to make accurate comparison
filterpaper = NCBG_PR_frassdata[c(1154:1169,1182:1197),]
#subset weight & pcs by unique circles, then sum (to account for additional days collected in milkjug)
filtermass = aggregate(Frass.mass..mg. ~ Survey, data = filterpaper, sum)
filterpcs = aggregate(Frass.number ~ Survey, data = filterpaper, sum)
#next step is to identify which milk jug trap is near which filter trap, then combine all into one data set  
#plot comparing sorted pcs
t = merge(filterpcs, filtermass, by = "Survey")
# look at dplyr left_merge

#plot comparing sorted weight

# COMPARISONS TO DO

# Raw Img to Sorted Img
# Sorted Pieces to Sorted Weight
# Filter paper to Milk jug (sorted/sorted)
# Before and after "rain"
# Sorted weight to Img_raw  







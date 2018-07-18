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

#renaming data sets
data = frassLoad(open = T)
NCBG_PR_frassdata = frassData(open = T)

#removing outliers in frassLoad
dataWO = data[data$Weight_Raw < 50,]
data_rawpcsWO = data[data$Pieces_Raw < 60,]
data_srtdpcsWO = data[data$Pieces_Sorted < 50,]
data_img_exlc_outlier = data[data$Img_Sorted < 20, ]

# Linear model & plot describing weight_sorted vs weight_raw
raw_sort_outlier_excl = lm(Weight_Sorted ~ Weight_Raw, data = dataWO)
plot(data$Weight_Raw[data$Weight_Raw<50], data$Weight_Sorted[data$Weight_Raw<50],main = "Frass Weight Comparison (mg.)", xlab = "Weight Raw", ylab = "Weight Sorted", pch = 17, cex = 1, col = 'red')
abline(raw_sort_outlier_excl)
sortraw_sum = summary(raw_sort_outlier_excl)
sortraw_sum_r2 = sortraw_sum$adj.r.squared
mylabel = bquote(italic(R)^2 == .(format(sortraw_sum_r2, digits = 3)))
text(x = 3.3, y = 16.2, labels = mylabel)

## Linear models & plots showing raw/sorted img against raw/sorted pcs to describe how much sorting changes % of area 
# Raw/Sorted img
plot(data$Pieces_Sorted[data$Pieces_Raw<60], data$Img_Raw[data$Pieces_Raw<60], main = "Raw Frass Comparison: # of Pieces vs. % of Area", 
     xlab = "Total Pieces", ylab = "% of Area (unsorted)", pch = 20, cex = 1, col = 'orange')
raw_pcs = lm(Img_Raw ~ Pieces_Sorted, data = data)
raw_pcs_outlier_excl = lm(Img_Raw ~ Pieces_Sorted, data = data_rawpcsWO)
abline(raw_pcs_outlier_excl)
sortrawimg_sum = summary(raw_pcs_outlier_excl)
sortrawimg_sum_r2 = sortrawimg_sum$adj.r.squared
mylabel = bquote(italic(R)^2 == .(format(sortrawimg_sum_r2, digits = 3)))
text(x = 4.9, y = 6.1, labels = mylabel)

# Raw/Sorted pieces
plot(data$Pieces_Sorted[data$Pieces_Sorted<50], data$Img_Sorted[data$Pieces_Sorted<50], main = "Sorted Frass Comparison: # of Pieces vs. % of Area", 
     xlab = "Pieces Sorted", ylab = "% of Area", pch = 20, cex = 1, col = 'blue')
sort_pcs_outlier_excl = lm(Img_Sorted ~ Pieces_Sorted, data = data_srtdpcsWO)
abline(sort_pcs_outlier_excl)
sortdpcs_sum = summary(sort_pcs_outlier_excl)
sortdpcs_sum_r2 = sortdpcs_sum$adj.r.squared
mylabel = bquote(italic(R)^2 == .(format(sortdpcs_sum_r2, digits = 3)))
text(x = 7, y = 7.3, labels = mylabel)

# Raw image versus weight sorted
# Excluded outlier
plot(dataWO$Img_Raw, dataWO$Weight_Sorted, main = "Raw Image v. Weight of sorted Frass", xlab = "Raw Image Percent Coverage", 
     ylab = "Weight Sorted (mg)", col = 'orange' , pch = 18)
raw_img = lm(dataWO$Weight_Sorted ~ dataWO$Img_Raw, data = dataWO)
abline(raw_img)
imgwght_sum = summary(raw_img)
imgwght_sum_r2 = imgwght_sum$adj.r.squared
mylabel = bquote(italic(R)^2 == .(format(imgwght_sum_r2, digits = 3)))
text(x = 1.7, y = 16.2, labels = mylabel)

# Img_raw vs.Img_sort
plot(data$Img_Raw[data$Img_Sorted<20], data$Img_Sorted[data$Img_Sorted<20], main = "Comparison Img_Raw vs. Img_Sorted (% of area estimate)", xlab = "Raw Img.", ylab ="Sorted Img.", col = 'violet', pch = 20)
rawsort_img = lm(data$Img_Sorted ~ data$Img_Raw, data = data)
abline(rawsort_img)
rawsort_img_sum = summary(rawsort_img)
img_rawsrt_r2 = rawsort_img_sum$adj.r.squared
mylabel = bquote(italic(R)^2 == .(format(img_rawsrt_r2, digits = 3)))
text(x = 2, y = 16, labels = mylabel)

# plot comparing sorted pcs to sorted weight
plot(data$Pieces_Sorted[data$Pieces_Sorted<100], data$Weight_Sorted[data$Pieces_Sorted<100], main = "Comparison: Sorted Pieces vs. Sorted Weight (% of area)", xlab = "Sorted Pieces", ylab ="Sorted Weight", col = 'orange', pch = 20)
sorted_lm = lm(data$Weight_Sorted[data$Pieces_Sorted<100] ~ data$Pieces_Sorted[data$Pieces_Sorted<100], data = data_srtdpcsWO)
abline(sorted_lm)
sorted_lm_sum = summary(sorted_lm)
sorted_lm_r2 = sorted_lm_sum$adj.r.squared
mylabel = bquote(italic(R)^2 == .(format(sorted_lm_r2, digits = 3)))
text(x = 9, y = 30, labels = mylabel)


##below plot in progress
#NCBG Comparison: Filter paper vs. Milk jug collection. Traps set on 3rd.
#Filter paper collected on the 6th & 10th
#milk jug collected on 10th
#must sum filter paper frass by circle to make accurate comparison
#create new data table for filter paper from 7/6
filterfrass_all = NCBG_PR_frassdata[c(1154:1169,1182:1197),]
#isolate by frass trap site
srtd_filterpaper = filterfrass_all[ ! filterfrass_all$Survey %in% c("1DBN","2DBS", "3DBV","4DCE","5DCI","6DCM","7DCQ","8DCV"), ]
#Sum values of same filter paper frass traps (to account for the additional collection days)
filtermass = aggregate(Frass.mass..mg. ~ Survey, data = srtd_filterpaper, sum)
filterpcs = aggregate(Frass.number ~ Survey, data = srtd_filterpaper, sum)
#create data set for summed values of filter paper traps, with pieces and mass merged
filter_sum = merge(filterpcs, filtermass, by = "Survey")
filter_sum$Date.Collected <- "7/10/2018"

#create data set with normal values and collection dates
filterdates_nonsum = NCBG_PR_frassdata[c(1210:1225, 1238:1253), c("Survey","Frass.mass..mg.","Frass.number", "Date.Collected")]
#isolate by frass trap site
filter_normal = filterdates_nonsum[ ! filterdates_nonsum$Survey %in% c("1DBN","2DBS", "3DBV","4DCE","5DCI","6DCM","7DCQ","8DCV"), ]
#combine adjusted filter paper frass data set with normal data set
filterpaper = rbind(filter_sum, filter_normal)
  
#create new data table for milk jug isolating by frass trap site
milkjugs = data[c(73:96),c("Survey","Weight_Sorted", "Pieces_Sorted", "Date.Collected")]

#merge both data sets to compare milk jug and filter paper mass and peices
frasstrapscomp <- filterpaper %>% 
  left_join(milkjugs, by = c("Survey", "Date.Collected"))
setnames(frasstrapscomp, old=c("Weight_Sorted","Pieces_Sorted", "Frass.number","Frass.mass..mg."), new=c("FrassNumber_milkjug", "FrassMass_milkjug","FrassNumber_filterpaper","FrassMass_filterpaper"))

#plotting filter paper vs. milk jug mass & pieces
#next step is to normalize the data
plot(frasstrapscomp$FrassNumber_filterpaper, frasstrapscomp$FrassMass_filterpaper, main = "Frass Collection Method: Filter Paper vs. Milk Jug (un-normalized)", xlab = "Number of Pieces", ylab ="Mass /mg.",  
     col = 'orange', pch = 20, xlim=c(-5, 92), ylim=c(0, 100))
points(frasstrapscomp$FrassNumber_milkjug, frasstrapscomp$FrassMass_filterpaper, col = 'blue', pch = 20)

# COMPARISONS TO DO 

# Raw Img to Sorted Img - AD complete
# Sorted Pieces to Sorted Weight - AD complete
# Filter paper to Milk jug (sorted/sorted) - AD  in progress
# Before and after "rain"
# Sorted weight to Img_raw - AZ







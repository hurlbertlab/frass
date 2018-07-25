#########################################################
## Data Setup ##

# Milk jug frass data
# Function for reading in frass data from GoogleDoc
# *if aim is to backup GoogleDoc and write to disk only, then open = F and write = T
# *if aim is to use data without writing to disk, then open = T and write = F

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

# Filter paper frass data
# Function for reading in frass data from GoogleDoc
# *if aim is to backup GoogleDoc and write to disk only, then open = F and write = T
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

# renaming data sets
data = frassLoad(open = T)
NCBG_PR_frassdata = frassData(open = T) 

# removing outliers in frassLoad
dataWO = data[data$Weight_Raw < 50,]
data_rawpcsWO = data[data$Pieces_Raw < 60,]
data_srtdpcsWO = data[data$Pieces_Sorted < 50,]
data_img_exlc_outlier = data[data$Img_Sorted < 20, ]

library(gsheet)
library(dplyr)
library(tidyr)
library(data.table)

##########################################################

## Plotting and Analysis: Milk Jug Frass Collection ## (7 plots) 

# Weight of Raw Frass / Weight of Sorted Frass
raw_sort_outlier_excl = lm(Weight_Sorted ~ Weight_Raw, data = dataWO)
plot(data$Weight_Raw[data$Weight_Raw<50], data$Weight_Sorted[data$Weight_Raw<50],
     main = "Frass Weight Comparison (mg.)", 
     xlab = "Weight Raw", ylab = "Weight Sorted", 
     pch = 17, cex = 1, col = 'goldenrod2')
abline(raw_sort_outlier_excl)
sortraw_sum = summary(raw_sort_outlier_excl)
sortraw_sum_r2 = sortraw_sum$adj.r.squared
mylabel = bquote(italic(R)^2 == .(format(sortraw_sum_r2, digits = 3)))
text(x = 3.3, y = 30, labels = mylabel)

# Sorted Frass Comparison: Pieces vs. Weight (mg.)
plot(data$Pieces_Sorted[data$Pieces_Sorted<100], data$Weight_Sorted[data$Pieces_Sorted<100], 
     main = "Sorted Frass Comparison:\n Pieces vs. Weight", 
     xlab = "Pieces", ylab ="Weight (mg.)", col = 'orange', pch = 20)
sorted_lm = lm(data$Weight_Sorted[data$Pieces_Sorted<100] ~ data$Pieces_Sorted[data$Pieces_Sorted<100], data = data_srtdpcsWO)
abline(sorted_lm)
sorted_lm_sum = summary(sorted_lm)
sorted_lm_r2 = sorted_lm_sum$adj.r.squared
mylabel = bquote(italic(R)^2 == .(format(sorted_lm_r2, digits = 3)))
text(x = 12, y = 30, labels = mylabel)

# % of Area (Unsorted) / Sorted Weight
plot(dataWO$Img_Raw, dataWO$Weight_Sorted, 
     main = "Frass Comparison:\n Image-based Estimate vs. Mass", 
     xlab = "% of Area (Unsorted)",  ylab = "Sorted Weight (mg.)", 
     ylim = c(-.5,20), col = 'orange' , pch = 18)
raw_img = lm(dataWO$Weight_Sorted ~ dataWO$Img_Raw, data = dataWO)
abline(raw_img)
imgwght_sum = summary(raw_img)
imgwght_sum_r2 = imgwght_sum$adj.r.squared
mylabel = bquote(italic(R)^2 == .(format(imgwght_sum_r2, digits = 3)))
text(x = 1.9, y = 19.7, labels = mylabel)

# % of Area Comparison: Unsorted Frass / Sorted Frass
plot(data$Img_Raw[data$Img_Sorted<20], data$Img_Sorted[data$Img_Sorted<20], 
     main = "% of Area Comparison", 
     xlab = "Unsorted Frass", ylab ="Sorted Frass", 
     ylim = c(-.5, 10), col = 'deepskyblue2', pch = 20)
rawsort_img = lm(data$Img_Sorted ~ data$Img_Raw, data = data)
abline(rawsort_img)
rawsort_img_sum = summary(rawsort_img)
img_rawsrt_r2 = rawsort_img_sum$adj.r.squared
mylabel = bquote(italic(R)^2 == .(format(img_rawsrt_r2, digits = 3)))
text(x = 2, y = 9.7, labels = mylabel)

# Sorted Pieces / % of Area (Unsorted) -- to evaluate an image based estimate for mass
plot(data$Pieces_Sorted[data$Pieces_Raw<60], data$Img_Raw[data$Pieces_Raw<60], 
     main = "Frass Comparison:\n # of Pieces vs. % of Area", 
     xlab = "Sorted Frass", ylab = "% of Area (Unsorted)", 
     xlim = c(0,30.009), pch = 20, cex = 1, col = 'orange')
raw_pcs = lm(Img_Raw ~ Pieces_Sorted, data = data)
raw_pcs_outlier_excl = lm(Img_Raw ~ Pieces_Sorted, data = data_rawpcsWO)
abline(raw_pcs_outlier_excl)
sortrawimg_sum = summary(raw_pcs_outlier_excl)
sortrawimg_sum_r2 = sortrawimg_sum$adj.r.squared
mylabel = bquote(italic(R)^2 == .(format(sortrawimg_sum_r2, digits = 3)))
text(x = 4.9, y = 6.1, labels = mylabel)

# Sorted Pieces / % of Area (Sorted) -- to evaluate an image based estimate for mass
plot(data$Pieces_Sorted[data$Pieces_Sorted<50], data$Img_Sorted[data$Pieces_Sorted<50], 
     main = "Frass Comparison:\n # of Pieces vs. % of Area", 
     xlab = "Sorted Frass", ylab = "% of Area (Sorted)",
     pch = 20, cex = 1, col = 'deepskyblue2')
sort_pcs_outlier_excl = lm(Img_Sorted ~ Pieces_Sorted, data = data_srtdpcsWO)
abline(sort_pcs_outlier_excl)
sortdpcs_sum = summary(sort_pcs_outlier_excl)
sortdpcs_sum_r2 = sortdpcs_sum$adj.r.squared
mylabel = bquote(italic(R)^2 == .(format(sortdpcs_sum_r2, digits = 3)))
text(x = 7, y = 7.3, labels = mylabel)

# Volume vs. Mass -- to evaluate an image-based estimate for weight
plot(data$Volume_raw, data$Weight_Sorted, 
     main = "Estimated Volume As A Proxy For Frass Weight", 
     xlab = "Unsorted Frass Volume (mm^3)", 
     ylab = "Sorted Frass Weight (mg.)", ylim = c(.5, 35), pch = 16, 
     col = "orange")
volwght.lm = lm(data$Weight_Sorted ~ data$Volume_raw, data = data)
abline(volwght.lm)
volwght.lm_sum = summary(volwght.lm)
volwght.lm_sum_r2 = volwght.lm_sum $adj.r.squared
mylabel = bquote(italic(R)^2 == .(format(volwght.lm_sum_r2, digits = 3)))
text(x = 37, y = 34, labels = mylabel)

####################################################################

## Methods Comparison: Filter paper vs. Milk jug collection ## 

# Data setup #

# Traps were set on Jul 3rd, Filter paper collected on the 6th & 10th; 
# but milk jug collected on 10th - must sum filter paper frass by circle 
# in order to make an accurate comparison

# Create new data table for filter paper, select for all NCBG rows starting 7/6 
filterfrass_all = NCBG_PR_frassdata[c(1154:1169,1182:1197),]
# Isolate by frass trap site (remove non-milk jug adjacent code)
srtd_filterpaper = filterfrass_all[ ! filterfrass_all$Survey %in% c("1DBD","2DBS", "3DBV","4DCE","5DCI","6DCM","7DCQ","8DCV"), ]
# Sum values of identical survey codes (to account for the additional collection days)
filtermass = aggregate(Frass.mass..mg. ~ Survey, data = srtd_filterpaper, sum)
filterpcs = aggregate(Frass.number ~ Survey, data = srtd_filterpaper, sum)
# Create new data table with these summed values, merging pieces & mass data
filter_sum = merge(filterpcs, filtermass, by = "Survey")
filter_sum$Date.Collected <- "7/10/2018"

# Create data table for normal values and collection dates
filterdates_nonsum = NCBG_PR_frassdata[c(1210:1225, 1238:1253), c("Survey","Frass.mass..mg.","Frass.number", "Date.Collected")]
# Isolate by frass trap site (remove non-milk jug adjacent code)
filter_normal = filterdates_nonsum[ ! filterdates_nonsum$Survey %in% c("1DBD","2DBS", "3DBV","4DCE","5DCI","6DCM","7DCQ","8DCV"), ]
# Combine this "normal" data table with the adjusted table 
filterpaper = rbind(filter_sum, filter_normal)
  
# Create new table for milk jug data, isolating by frass trap site, get rid of NAs
milkjugs = data[c(73:96),c("Survey","Weight_Sorted", "Pieces_Sorted", "Date.Collected")]
na.omit(milkjugs)

# Merge both data tables to compare milk jug and filter paper by mass and pieces
frasstrapscomp <- filterpaper %>% 
  left_join(milkjugs, by = c("Survey", "Date.Collected"))
setnames(frasstrapscomp, old=c("Weight_Sorted","Pieces_Sorted", "Frass.number","Frass.mass..mg."), new=c("FrassNumber_milkjug", "FrassMass_milkjug","FrassNumber_filterpaper","FrassMass_filterpaper"))
frasstrapscomp = frasstrapscomp[-c(4),]

# Add new columns with mass and pieces data normalized for area
# Area of milk jug = 171.9; Area of filter paper = 433.6 cm^2 
frasstrapscomp = transform(frasstrapscomp, FrassNumber.adj_filterpaper = FrassNumber_filterpaper / 433.6)
frasstrapscomp = transform(frasstrapscomp, FrassMass.adj_filterpaper = FrassMass_filterpaper / 433.6)
frasstrapscomp = transform(frasstrapscomp, FrassNumber.adj_milkjug = FrassNumber_milkjug / 171.9)
frasstrapscomp = transform(frasstrapscomp, FrassMass.adj_milkjug = FrassMass_milkjug / 171.9)


## Plotting and Analysis: Filter Paper vs. Milk Jug  ## (4 plots)

# Filter Paper vs. Milk Jug  - both mass & pieces (non-normalized)
plot(frasstrapscomp$FrassNumber_filterpaper, frasstrapscomp$FrassMass_filterpaper, 
     main = "Frass Collection Method:\nFilter Paper vs. Milk Jug (non-normalized)", 
     xlab = "Pieces", ylab ="Mass (mg.)",  
     col = 'orange', pch = 20, xlim=c(-5, 90), ylim=c(5.5, 95))
points(frasstrapscomp$FrassNumber_milkjug, frasstrapscomp$FrassMass_milkjug, 
       col = 'deepskyblue2', pch = 20)

# Filter Paper vs. Milk Jug  - both mass & pieces (normalized)
par(mar=c(4, 5, 5, 3)) # Bottom, Left, Top, Right
plot(frasstrapscomp$FrassNumber.adj_filterpaper, frasstrapscomp$FrassMass.adj_filterpaper, 
     main = "Frass Collection Method:\nFilter Paper vs. Milk Jug (normalized)", 
     xlab = expression(paste("Pieces per ", cm^2)), 
     ylab = expression(paste("Mg. per ", cm^2)),  
     col = 'orange', pch = 20, cex = 1, xlim=c(-.01, .2), ylim=c(.015, .21))
points(frasstrapscomp$FrassNumber.adj_milkjug, frasstrapscomp$FrassMass.adj_milkjug, 
     col = 'deepskyblue2', pch = 20)

# Filter Paper vs. Milk Jug by mass collected
par(mar=c(4, 5, 5, 3)) # Bottom, Left, Top, Right
plot(frasstrapscomp$FrassMass.adj_filterpaper, frasstrapscomp$FrassMass.adj_milkjug, 
     main = expression(paste("Frass Collection Method Comparison:" ~ "Mass (mg.) per" ~ cm^{2})),
     xlab = expression(paste("Filter Paper ")), 
     ylab = expression(paste("Milk Jug ")),  
     col = 'deepskyblue2', pch = 19, cex = .8, ylim=c(-.01, .52))
methodcompare.lm = lm(frasstrapscomp$FrassMass.adj_milkjug ~ frasstrapscomp$FrassMass.adj_filterpaper, data = frasstrapscomp )
abline(methodcompare.lm)
methodcompare_sum = summary(methodcompare.lm)
methodcompare_sum_r2 = methodcompare_sum$adj.r.squared
mylabel = bquote(italic(R)^2 == .(format(methodcompare_sum_r2, digits = 3)))
text(x = .035, y = .52, labels = mylabel)

# Filter Paper vs. Milk Jug by pieces collected
par(mar=c(4, 5, 5, 3)) # Bottom, Left, Top, Right
plot(frasstrapscomp$FrassNumber.adj_filterpaper, frasstrapscomp$FrassNumber.adj_milkjug, 
     main = expression(paste("Frass Collection Method Comparison:" ~ "Pieces per" ~ cm^{2})), 
     xlab = expression(paste("Filter Paper ")), 
     ylab = expression(paste("Milk Jug ")),  
     col = 'deepskyblue2', pch = 19, cex = .8, 
     xlim = c(.03, .505), ylim=c(-.01, .2))
methodcompare_num.lm = lm(frasstrapscomp$FrassNumber.adj_milkjug ~ frasstrapscomp$FrassNumber.adj_filterpaper, data = frasstrapscomp )
abline(methodcompare_num.lm)
methodcompare_numsum = summary(methodcompare_num.lm)
methodcompare_numsum_r2 = methodcompare_numsum$adj.r.squared
mylabel = bquote(italic(R)^2 == .(format(methodcompare_numsum_r2, digits = 3)))
text(x = .085, y = .2, labels = mylabel)
# Script for analyzing frass data

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


data = frassLoad(open = T)
dataWO = data[data$Weight_Raw < 50,]
data_rawpcsWO = data[data$Pieces_Raw < 60,]
data_srtdpcsWO = data[data$Pieces_Sorted < 50,]
data_img_exlc_outlier = data[data$Img_Sorted < 20, ]


# Linear model 
raw_sort = lm(Weight_Sorted ~ Weight_Raw, data = data)
raw_sort_outlier_excl = lm(Weight_Sorted ~ Weight_Raw, data = dataWO)
sort_img = lm(Img_Sorted ~ Weight_Sorted, data = dataWO)

#Plotting
plot(data$Weight_Raw, data$Weight_Sorted, main = "Frass Weight Comparison (mg.)", xlab = "Weight Raw", ylab = "Weight Sorted", pch = 17, cex = 1, col = 'red')
plot(data$Weight_Raw[data$Weight_Raw<50], data$Weight_Sorted[data$Weight_Raw<50],main = "Frass Weight Comparison (mg.)", xlab = "Weight Raw", ylab = "Weight Sorted", pch = 17, cex = 1, col = 'red')
abline(raw_sort_outlier_excl)

## Linear models & plots showing raw/sorted img against raw/sorted pcs to describe how much sorting changes % of area 
# Raw 
plot(data$Pieces_Sorted[data$Pieces_Raw<60], data$Img_Raw[data$Pieces_Raw<60], main = "Raw Frass Comparison: # of Pieces vs. % of Area", 
     xlab = "Total Pieces", ylab = "% of Area (unsorted)", pch = 20, cex = 1, col = 'orange')
raw_pcs = lm(Img_Raw ~ Pieces_Sorted, data = data)
raw_pcs_outlier_excl = lm(Img_Raw ~ Pieces_Raw, data = data_rawpcsWO)
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
#excluded the outlier
plot(dataWO$Img_Raw, dataWO$Weight_Sorted, main = "Raw Image v. Weight of sorted Frass", xlab = "Raw Image Percent Coverage", 
     ylab = "Weight Sorted (mg)", col = 'orange' , pch = 18)
raw_img = lm(dataWO$Weight_Sorted ~ dataWO$Img_Raw, data = dataWO)
abline(raw_img)
summary(raw_img)

# Img_raw vs. Img_sort
plot(data$Img_Raw[data$Img_Sorted<20], data$Img_Sorted[data$Img_Sorted<20], main = "Comparison Img_Raw vs. Img_Sorted (% of area estimate)", xlab = "Raw Img.", ylab ="Sorted Img.", col = 'violet', pch = 20)
rawsort_img = lm(data_img_exlc_outlier$Img_Sorted ~ data_img_exlc_outlier$Img_Raw, data = data_img_exlc_outlier)
abline(rawsort_img)
summary(rawsort_img)

# COMPARISONS TO DO

# Raw Img to Sorted Img
# Sorted Pieces to Sorted Weight
# Filter paper to Milk jug (sorted/sorted)
# Before and after "rain"
# Sorted weight to Img_raw  







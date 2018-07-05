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
plot(data$Pieces_Raw, data$Img_Raw, main = "Raw Frass Comparison: # of Pieces vs. % of Area", 
     xlab = "Total Pieces (frass w/ debris)", ylab = "% of Area", pch = 20, cex = 1, col = 'orange')
raw_pcs = lm(Img_Raw ~ Pieces_Raw, data = data)
abline(raw_pcs)
summary(raw_pcs)

# Sorted
plot(data$Pieces_Sorted, data$Img_Sorted, main = "Sorted Frass Comparison: # of Pieces vs. % of Area", 
     xlab = "Pieces Sorted", ylab = "% of Area", pch = 20, cex = 1, col = 'blue')
sort_pcs = lm(Img_Sorted ~ Pieces_Sorted, data = data)
abline(sort_pcs)
summary(sort_pcs)

# Raw image versus weight sorted
#excluded the outlier
plot(dataWO$Img_Raw, dataWO$Weight_Sorted, main = "Raw Image v. Weight of sorted Frass", xlab = "Raw Image Percent Coverage", 
     ylab = "Weight Sorted (mg)", col = 'orange' , pch = 18)
raw_img = lm(dataWO$Img_Raw ~ dataWO$Weight_Sorted)
abline(raw_img)
summary(raw_img)




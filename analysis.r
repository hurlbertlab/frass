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

#Next plot should be comparing img data to wgt. raw then wgt. sorted to see how much sorting changes weight estimates
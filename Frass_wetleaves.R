# Looking at wet leaf data and impact on bug abundance
#Savannah Carter
#4.4.2

#libraries
library(gsheet)
library(dplyr)
library(tidyr)
library(purrr)
library(rstatix)
library(tidyverse)
library(jsonlite)
library(daymetr)
library(gridExtra)
library(zoo)

# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#   Datasets needed:
# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+

#Loading in CC data 
options(timeout = 300)  
api_url <- "https://api.github.com/repos/hurlbertlab/caterpillars-analysis-public/contents/data"
files <- fromJSON(api_url)
dataset_file <- files$name[grepl("fullDataset", files$name, ignore.case = TRUE)]
# pick the latest one
latest_file <- dataset_file[1]
github_raw <- "https://raw.githubusercontent.com/hurlbertlab/caterpillars-analysis-public/master/data/"
fullDataset <- read.csv(paste0(github_raw, latest_file))

#----------------------------------------------------------------------------------
#loading in temp data
tmp_file = tempfile(fileext = ".csv")

AnomalySites = fullDataset[, c("Name","Year", "Latitude", "Longitude")] %>% 
  filter(Name %in% c("NC Botanical Garden", "Prairie Ridge Ecostation")) %>% 
  group_by(Name, Latitude, Longitude) %>% 
  summarise(n = n()) %>% 
  select(-n) %>% as.data.frame()   

AnomalyDaymetr = AnomalySites %>% 
  rename(
    site = Name,
    lat = Latitude,
    lon = Longitude
  ) %>%
  write.csv(tmp_file, row.names = FALSE)

# pass temp CSV to function
TempAnomaly = download_daymet_batch(
  file_location = tmp_file,
  start = 2015,
  end = 2024, # this is the most recent available in daymetr
  internal = TRUE)
# remove temporary file 
unlink(tmp_file)

TempAnomalyData_clean <- lapply(TempAnomaly, function(x) {
  x$data %>% mutate(site = x$site,
                    Latidue = x$latitude,
                    Longitude = x$longitude)})

AllTemp= bind_rows(TempAnomalyData_clean) #this is the file we want (has all NCBG and PR data)

# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#   adding up wet leaf vs dry leaf totals for each arthropod group (all sites merged and all years merged)
# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
raw_totals <- fullDataset %>% 
  distinct(ID, .keep_all = TRUE)%>% #only keeps first instance of ID, make sure all unique
  group_by(Group) %>% #fixed group with inaturalists?
  summarise(
    sum_all_dry = sum(WetLeaves == 0),
    sum_all_wet = sum(WetLeaves == 1),
    sum_all = sum(sum_all_dry, sum_all_wet)) %>%
  mutate(dry_ratio = sum_all_dry/sum_all,
         wet_ratio = sum_all_wet/sum_all)

raw_totals2 <- fullDataset %>%
  group_by(OriginalGroup) %>% #original bug group
  summarise(
    sum_all_dry2 = sum(WetLeaves == 0),
    sum_all_wet2 = sum(WetLeaves == 1))
raw_totals_both <- left_join(raw_totals, raw_totals2, by = c("Group" = "OriginalGroup"))

# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#   adding up wet leaf vs dry leaf totals for each arthropod group by site (NCBG and PR) and year
# *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#botanical garden
raw_totals_ncbg <- fullDataset %>%
  filter(Name == "NC Botanical Garden") %>%
  distinct(ID, .keep_all = TRUE)%>% #only keeps first instance of ID, make sure all unique
  group_by(Group, Year) %>% #fixed group with inaturalists?
  summarise(
    sum_all_dry = sum(WetLeaves == 0),
    sum_all_wet = sum(WetLeaves == 1),
    sum_all = sum(sum_all_dry, sum_all_wet)) %>%
  mutate(dry_ratio = sum_all_dry/sum_all,
         wet_ratio = sum_all_wet/sum_all) %>%
  mutate(site = "8892356")
#prairie ridge
raw_totals_pr <- fullDataset %>%
  filter(Name == "Prairie Ridge Ecostation") %>%
  distinct(ID, .keep_all = TRUE)%>% #only keeps first instance of ID, make sure all unique
  group_by(Group, Year) %>% #fixed group with inaturalists?
  summarise(
    sum_all_dry = sum(WetLeaves == 0),
    sum_all_wet = sum(WetLeaves == 1),
    sum_all = sum(sum_all_dry, sum_all_wet)) %>%
  mutate(dry_ratio = sum_all_dry/sum_all,
         wet_ratio = sum_all_wet/sum_all) %>%
  mutate(site = "117")
#combine by rows 
raw_totals_frass_sites <- rbind(raw_totals_ncbg, raw_totals_pr)















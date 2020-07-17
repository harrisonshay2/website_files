knitr::opts_chunk$set(echo=F)
knitr::opts_chunk$set(warning=F)
knitr::opts_chunk$set(message=F)

setwd("C:/Users/Harri/OneDrive/Documents/centauri/conjunction_analysis")
library(tidyverse)
library(readr)
library(kableExtra)
library(knitr)
library(lubridate)
library(DT)
library(RColorBrewer)
library(stringr)
library(scales)
Sys.setenv(TZ='EST')
mcma_objs = readRDS("~/centauri/conjunction_analysis/RDSfiles/mcma_objs")
all_conjs = readRDS("~/centauri/conjunction_analysis/RDSfiles/all_conjs")
all_conjs_expanded = readRDS("~/centauri/conjunction_analysis/RDSfiles/all_conjs_expanded")
derelicts = readRDS("~/centauri/conjunction_analysis/RDSfiles/derelicts")
derelictDat = readRDS("~/centauri/conjunction_analysis/RDSfiles/derelictDatNew")
alt_bins = readRDS("~/centauri/conjunction_analysis/RDSfiles/alt_bins")
file_list = readRDS("~/centauri/conjunction_analysis/RDSfiles/file_list")
all_conjs_2016 = readRDS("~/centauri/conjunction_analysis/RDSfiles/all_conjs_2016")
today = toupper(strftime(Sys.Date(), format="%d%b%Y")) # current day
path = "conj_data/"





# add new conjunction files to all_conjs dataframe

# read in new conjunction files
file_list_new = list.files("C:/Users/Harri/OneDrive/Documents/centauri/conjunction_analysis/conj_data")
file_list_new = file_list_new[!(file_list_new %in% file_list)] # only the new conjunctions

colnames = c("PrimarySatellite","SecondarySatellite","TCA_EpDay",
             "TCA_UTCG","Range","RangeX","RangeY","RangeZ","Velocity",
             "VelocityX","VelocityY","VelocityZ","Latitude","Longitude",
             "Altitude","PrimaryAge","SecondaryAge","PrimaryCluster",
             "SecondaryCluster","DateGenerated","del")

all_conjs_new = data.frame()
for (i in 1:length(file_list_new)) {
  file = paste0("C:/Users/Harri/OneDrive/Documents/centauri/conjunction_analysis/conj_data/", file_list_new[i])
  
  firstLine = readLines(file, n=2)[2]
  
  if (str_count(firstLine, ',') == 20) { # if file has trailing commas
    temp_data = read_csv(file, skip=1, col_names = colnames, 
                         col_types = "ccncnnnnnnnnnnncccccc") %>%
      select(-del)
  } else {
    temp_data = read_csv(file, skip=1, 
                         col_names = colnames[-length(colnames)], 
                         col_types = "ccncnnnnnnnnnnnccccc")
  }
  
  all_conjs_new = rbind(all_conjs_new, temp_data) #for each iteration, bind the new data to the building dataset
}

mycols <- '(PrimaryCluster, SecondaryCluster)'
minf <- paste0('min',mycols)
maxf <- paste0('max',mycols)

all_conjs_new = all_conjs_new %>%
  mutate(DateGenerated = parse_date_time(DateGenerated, tz="EST", 
                                         orders=c("%Y-%m-%d %H:%M:%S", "%m/%d/%y %H:%M")),
         date = DateGenerated - 24*60*60,
         utcg = if_else(nchar(TCA_UTCG) > 7,
                        as.POSIXct(TCA_UTCG, format="%Y-%m-%d %H:%M:%S"),
                        date + TCA_EpDay*24*60*60),
         TCA_UTCG = utcg) %>% 
  select(-c(date, utcg)) %>%
  rowwise() %>% 
  mutate(firstClust = eval(parse(text=minf)),
         secondClust = eval(parse(text=maxf)),
         clusters = paste(firstClust, secondClust, sep="-")) %>% 
  ungroup() %>%
  mutate(clusterLab = if_else(firstClust=="LEO" & secondClust=="LEO", "LEO",
                              if_else((firstClust=="LEO" & secondClust!="LEO") |
                                        (firstClust!="LEO" & secondClust=="LEO"), "LEO-other",
                                      if_else(firstClust=="HIGH" & secondClust=="HIGH", "HIGH",
                                              if_else((firstClust=="HIGH" & secondClust!="HIGH") | 
                                                        (firstClust!="HIGH" & secondClust=="HIGH"), "HIGH-other",
                                                      firstClust)))),
         clusterLab = factor(clusterLab,
                             levels = c("615", "775", "850", "975", "1200", "1500", "LEO","LEO-other","HIGH","HIGH-other"),
                             ordered = T))

# update file list
file_list = append(file_list, file_list_new)
saveRDS(file_list, "~/centauri/conjunction_analysis/RDSfiles/file_list")

#########
# WORST OFFENDER alg for new conjunctions
# persistence
alts = c(615,775,850,975, 1200,1500)
pers = c(25, 90, 150,1000,1600,1800)
lw1 <- loess(pers ~ alts)

# get operational satellites
opSats = derelictDat %>% filter(avgAlt < 2000 & operational)

combinedMass_v = vector()
persistence_v = vector()
#numOpSats_v = vector()
for (i in 1:nrow(all_conjs_new)) {
  conj = all_conjs_new[i, ]
  noradId1 = gsub("--.*", "", conj$PrimarySatellite)
  noradId2 = gsub("--.*", "", conj$SecondarySatellite)
  obj1 = filter(mcma_objs, noradId == noradId1)
  obj2 = filter(mcma_objs, noradId == noradId2)
  
  combinedMass = as.numeric(obj1$mass) + as.numeric(obj2$mass)
  persistence = if_else(conj$Altitude <= 615, 25,
                        if_else(conj$Altitude > 615 & conj$Altitude <= 1500, 
                                predict(lw1, conj$Altitude), 1000)) 
  
  combinedMass_v = append(combinedMass_v, toString(combinedMass))
  persistence_v = append(persistence_v, persistence)
}
all_conjs_new$combinedMass = combinedMass_v
all_conjs_new$persistence = persistence_v

# replace missing mass values
all_conjs_new = all_conjs_new %>%
  mutate(combinedMass = if_else(grepl(",", combinedMass, fixed = T), # if it contains a comma
                                as.numeric(gsub(",.*", "", combinedMass)), # make substring up to the comma
                                as.numeric(combinedMass) )) # otherwise don't change

## get SD op sats per conj
alt_bins = readRDS("~/centauri/conjunction_analysis/RDSfiles/alt_bins")
roundDown <- function(x) 10*floor(x/10)
library(zoo)
alt_bins = derelictDat %>% 
  filter(avgAlt < 2000 & operational) %>%
  mutate(altitude = roundDown((as.numeric(apogee) + as.numeric(perigee))/2)) %>% 
  group_by(altitude) %>% 
  summarise(numOpSats = n()) %>% 
  right_join(alt_bins, by="altitude") %>%
  mutate(numOpSats = replace_na(numOpSats, 0)) %>% 
  mutate(spatDensOpSats_1 = numOpSats / volume * (10^10)) %>%
  mutate(SD = rollmean(spatDensOpSats_1, k=5, na.pad=T)) %>% na.omit()

all_conjs_new = all_conjs_new %>% 
  mutate(altitude = roundDown(Altitude)) %>% 
  left_join(select(alt_bins, c(altitude, SD)), by="altitude") 

all_conjs_new = all_conjs_new %>%
  mutate(risk = (combinedMass + persistence + SD) / Range,
         conseq = combinedMass + persistence + SD,
         # initialize scaled variables
         combinedMass_s = 1, persistence_s = 1,
         SD_s = 1, Range_s = 1, risk_s = 1, conseq_s=1)

# append new conjunctions to previous
all_conjs = rbind(all_conjs, all_conjs_new)
saveRDS(all_conjs, "~/centauri/conjunction_analysis/RDSfiles/all_conjs") # save to RDS file

# adjust scaling 
all_conjs = all_conjs %>%
  mutate(combinedMass_s = scale(combinedMass),#rescale(combinedMass, to=c(0,1000)),
         persistence_s = scale(persistence), #rescale(persistence, to=c(0,1000)),
         SD_s = scale(SD), #rescale(SD, to=c(0,1000)),
         Range_s = if_else(Range >= 1, 10^(Range), 10), #rescale(Range, to=c(1e-4, 100)),
         conseq = combinedMass + persistence + SD,
         conseq_s = combinedMass_s + persistence_s + SD_s,
         risk_s = (combinedMass_s + persistence_s + SD_s) / (Range_s),
         risk = (combinedMass + persistence + SD) / Range_s) 

min_cm = min(all_conjs$combinedMass_s, na.rm=T)
min_pers = min(all_conjs$persistence_s)
min_SD = min(all_conjs$SD_s)

all_conjs = all_conjs %>%
  mutate(combinedMass_s = combinedMass_s + abs(min_cm),
         persistence_s = persistence_s + abs(min_pers),
         SD_s = SD_s + abs(min_SD),
         risk_s = (combinedMass_s + persistence_s + SD_s) / (Range_s))

# sum risk for each object:
# list all conjunctions by first sat, then by second sat, then bind by rows
firstSet = all_conjs %>%
  mutate(noradId = as.numeric(gsub("--.*", "", PrimarySatellite))) %>%
  dplyr::select(-c(PrimarySatellite, SecondarySatellite))

secondSet = all_conjs %>%
  mutate(noradId = as.numeric(gsub("--.*", "", SecondarySatellite))) %>%
  dplyr::select(-c(PrimarySatellite, SecondarySatellite))

# append new conjunctions to previous
all_conjs_expanded = rbind(firstSet, secondSet)
saveRDS(all_conjs_expanded, "~/centauri/conjunction_analysis/RDSfiles/all_conjs_expanded") # save to RDS file

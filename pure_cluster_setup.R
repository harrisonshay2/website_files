
knitr::opts_chunk$set(echo=F)
knitr::opts_chunk$set(warning=F)
knitr::opts_chunk$set(message=F)

setwd("~/centauri/")
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
all_conjs_pure = readRDS("~/centauri/conjunction_analysis/pureRDSfiles/all_conjs_pure")
derelicts = readRDS("~/centauri/conjunction_analysis/RDSfiles/derelicts")
derelictDat = readRDS("~/centauri/conjunction_analysis/RDSfiles/derelictDatNew")
alt_bins = readRDS("~/centauri/conjunction_analysis/RDSfiles/alt_bins")
file_list_pure = readRDS("~/centauri/conjunction_analysis/pureRDSfiles/file_list_pure")
today = toupper(strftime(Sys.Date(), format="%d%b%Y")) # current day
path = "~conjunction_analysis/conj_data/"
all_conjs_expanded_pure = readRDS("~/centauri/conjunction_analysis/pureRDSfiles/all_conjs_expanded_pure")



 # add new conjunction files to all_conjs dataframe
 
 # read in new conjunction files
 file_list_new_pure = list.files("C:/Users/Harri/OneDrive/Documents/centauri/conjunction_analysis/conj_data")
 file_list_new_pure = file_list_new_pure[!(file_list_new_pure %in% file_list_pure)] # only the new conjunctions
 
 colnames = c("PrimarySatellite","SecondarySatellite","TCA_EpDay",
              "TCA_UTCG","Range","RangeX","RangeY","RangeZ","Velocity",
              "VelocityX","VelocityY","VelocityZ","Latitude","Longitude",
              "Altitude","PrimaryAge","SecondaryAge","PrimaryCluster",
              "SecondaryCluster","DateGenerated","del")
 
 all_conjs_pure_new = data.frame()
 for (i in 1:length(file_list_new_pure)) {
   file = paste0("C:/Users/Harri/OneDrive/Documents/centauri/conjunction_analysis/conj_data/", file_list_new_pure[i])
   
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
   
   all_conjs_pure_new = rbind(all_conjs_pure_new, temp_data) #for each iteration, bind the new data to the building dataset
 }
 
 mycols <- '(PrimaryCluster, SecondaryCluster)'
 minf <- paste0('min',mycols)
 
 maxf <- paste0('max',mycols)
 
 all_conjs_pure_new = all_conjs_pure_new %>%
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
          clusters = paste(firstClust, secondClust, sep="-")) #%>% 
 
 # update file list
 file_list_pure = append(file_list_pure, file_list_new_pure)
 saveRDS(file_list_pure, "~/centauri/conjunction_analysis/pureRDSfiles/file_list_pure")
 
 # adding new clustelab
 all_conjs_pure_new = all_conjs_pure_new %>%
   mutate(noradId_1 = as.numeric(gsub("--.*", "", PrimarySatellite)),
          noradId_2 = as.numeric(gsub("--.*", "", SecondarySatellite))) %>%
   left_join(dplyr::select(derelictDat, c(noradId, cluster, mass)), by=c("noradId_1" = "noradId")) %>%
   rename_at(vars(c(cluster, mass)), function(x) paste0(x, "_1")) %>%
   left_join(dplyr::select(derelictDat, c(noradId, cluster, mass)), by=c("noradId_2" = "noradId")) %>%
   rename_at(vars(c(cluster, mass)), function(x) paste0(x, "_2")) %>%
   dplyr::select(-c(noradId_1, noradId_2))
 
 #getting rid of 'N' in cluster_1 
  all_conjs_pure_new = all_conjs_pure_new %>%
   mutate(cluster_1 = if_else(cluster_1 == "c775N", "elsewhere", if_else(cluster_1 == "c850N","elsewhere",if_else(cluster_1 == "c975N", "elsewhere",if_else(cluster_1 == "c1200N","elsewhere",if_else(cluster_1 == "c1500N","elsewhere",if_else(
     cluster_1 == "cc615","elsewhere", if_else(cluster_1 == "cc775","elsewhere",if_else(cluster_1 == "cc850","elsewhere",if_else(cluster_1 == "cc975","elsewhere", if_else(cluster_1 == "cc1200",
     "elsewhere",if_else(cluster_1 == "cc1500","elsewhere",as.character(cluster_1)))))))))))))
 
 #getting rid of 'N' in cluster_2 
  all_conjs_pure_new = all_conjs_pure_new %>%
    mutate(cluster_2 = if_else(cluster_2 == "c775N", "elsewhere", if_else(cluster_2 == "c850N","elsewhere",if_else(cluster_2 == "c975N", "elsewhere",if_else(cluster_2 == "c1200N","elsewhere",if_else(cluster_2 == "c1500N","elsewhere",if_else(
      cluster_1 == "cc615","elsewhere", if_else(cluster_1 == "cc775","elsewhere",if_else(cluster_1 == "cc850","elsewhere",if_else(cluster_1 == "cc975","elsewhere", if_else(cluster_1 == "cc1200",
       "elsewhere",if_else(cluster_1 == "cc1500","elsewhere",as.character(cluster_2)))))))))))))
 
 # pure clusterlab
 all_conjs_pure_new = all_conjs_pure_new%>%
 mutate(clusterLab_pure = if_else(cluster_1 == "elsewhere" & cluster_2 == "elsewhere", "no collision", if_else(cluster_1 == "elsewhere" | cluster_2 == "elsewhere", "no collision", if_else(cluster_1 == "cleo" | cluster_2 == "cleo", "no collision", if_else(cluster_1 == "cleo" & cluster_2 == "cleo", "no collision",if_else(cluster_1 == "chigh" | cluster_2 == "chigh","no collision",if_else(cluster_1 == "chigh" 
          & cluster_2 == "chigh", "no collision",as.character(cluster_1))))))),
                       clusterLab_pure = factor(clusterLab_pure,
                            levels = c("c775", "c850", "c975", "c1500","no collision"),
                           ordered = T))
 all_conjs_pure_new<-all_conjs_pure_new[!(all_conjs_pure_new$clusterLab_pure=="no collision"),]
 #all_conjs_pure = rbind(all_conjs_pure,all_conjs_pure_new)
 #saveRDS(all_conjs_pure,"~/centauri/all_conjs_pure")
 
 
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
  numOpSats_v = vector()
  for (i in 1:nrow(all_conjs_pure_new)) {
    conj = all_conjs_pure_new[i, ]
    noradId1 = gsub("--.*", "", conj$PrimarySatellite)
    noradId2 = gsub("--.*", "", conj$SecondarySatellite)
    obj1 = filter(derelictDat, noradId == noradId1)
    obj2 = filter(derelictDat, noradId == noradId2)
    
    combinedMass = as.numeric(obj1$mass) + as.numeric(obj2$mass)
    persistence = if_else(conj$Altitude <= 615, 25,
                          if_else(conj$Altitude > 615 & conj$Altitude <= 1500, 
                                  predict(lw1, conj$Altitude), 1000)) 
    
     combinedMass_v = append(combinedMass_v, toString(combinedMass))
    persistence_v = append(persistence_v, persistence)
  }
  all_conjs_pure_new$combinedMass = as.numeric(all_conjs_pure_new$mass_1) + as.numeric(all_conjs_pure_new$mass_2)
  all_conjs_pure_new$persistence = persistence_v
 
 # replace missing mass values
  all_conjs_pure_new = all_conjs_pure_new %>%
   mutate(combinedMass = if_else(grepl(",", combinedMass, fixed = T), # if it contains a comma
                                   as.numeric(gsub(",.*", "", combinedMass)), # make substring up to the comma
                                   as.numeric(combinedMass) )) # otherwise don't change
 
 # get SD op sats per conj
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
  
  all_conjs_pure_new = all_conjs_pure_new %>% 
    mutate(altitude = roundDown(Altitude)) %>% 
    left_join(select(alt_bins, c(altitude, SD)), by="altitude") 
  
  all_conjs_pure_new = all_conjs_pure_new %>%
    mutate(risk = (combinedMass + persistence + SD) / Range,
           conseq = combinedMass + persistence + SD,
          # initialize scaled variables
           combinedMass_s = 1, persistence_s = 1,
           SD_s = 1, Range_s = 1, risk_s = 1, conseq_s=1)
 
 # append new conjunctions to previous
 all_conjs_pure = rbind(all_conjs_pure, all_conjs_pure_new)
 #saveRDS(all_conjs, "RDSfiles/all_conjs") # save to RDS file
 
 # adjust scaling 
 all_conjs_pure$combinedMass_s = scale(all_conjs_pure$combinedMass)
 all_conjs_pure$persistence_s = scale(all_conjs_pure$persistence)
 all_conjs_pure$SD_s = scale(all_conjs_pure$SD)
 all_conjs_pure = all_conjs_pure %>%
   mutate(#combinedMass_s = scale(combinedMass),#rescale(combinedMass, to=c(0,1000)),
            #persistence_s = scale(persistence))#, #rescale(persistence, to=c(0,1000)),
            #SD_s = scale(SD), #rescale(SD, to=c(0,1000)),
            Range_s = if_else(Range >= 1, 10^(Range), 10), #rescale(Range, to=c(1e-4, 100)),
            conseq = combinedMass + persistence + SD,
            conseq_s = combinedMass_s + persistence_s + SD_s,
            risk_s = (combinedMass_s + persistence_s + SD_s) / (Range_s),
            risk = (combinedMass + persistence + SD) / Range_s) 
 
 min_cm = min(all_conjs_pure$combinedMass_s, na.rm=T)
 min_pers = min(all_conjs_pure$persistence_s)
 min_SD = min(all_conjs_pure$SD_s)
 
 all_conjs_pure = all_conjs_pure %>%
   mutate(combinedMass_s = combinedMass_s + abs(min_cm),
           persistence_s = persistence_s + abs(min_pers),
           SD_s = SD_s + abs(min_SD),
           risk_s = (combinedMass_s + persistence_s + SD_s) / (Range_s))
  
 #all_conjs_pure = rbind(all_conjs_pure,all_conjs_pure_new)
 saveRDS(all_conjs_pure,"~/centauri/conjunction_analysis/pureRDSfiles/all_conjs_pure")
  
 # pure_cluster_df = all_conjs_pure %>% 
 #   mutate(Range = Range) %>%
 #   group_by(clusterLab_pure) %>%
 #   arrange(Range) %>%
 #   mutate(rowid = 1, cumnum = cumsum(rowid))
 #all_conjs_expanded_pure = data.frame()

 firstSet = all_conjs_pure %>%
   mutate(noradId = as.numeric(gsub("--.*", "", PrimarySatellite))) %>%
   dplyr::select(-c(PrimarySatellite, SecondarySatellite))
 
 secondSet = all_conjs_pure %>%
   mutate(noradId = as.numeric(gsub("--.*", "", SecondarySatellite))) %>%
   dplyr::select(-c(PrimarySatellite, SecondarySatellite))
 
 # append new conjunctions to previous
 all_conjs_expanded_pure = rbind(firstSet, secondSet)
 saveRDS(all_conjs_expanded_pure, "~/centauri/conjunction_analysis/pureRDSfiles/all_conjs_expanded_pure") # save to RDS file

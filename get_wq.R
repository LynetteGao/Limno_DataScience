rm(list = ls())

library(rLakeAnalyzer)
library(tidyverse)
library(lubridate)
library(glmtools)
library(viridis)
library(patchwork)
library(RPostgreSQL)
library(odbc)
library(berryFunctions)
library(broom)
library(zoo)
library(signal)
library(Hmisc)


setwd('/Users/robertladwig/Documents/DSI/Limno_DataScience')
lks <- list.dirs(path = 'inst/extdata/', full.names = TRUE, recursive = F)
a=1
for (ii in lks){
  print(paste0('Running ',ii))
  
  data <- read.csv(paste0(ii,'/', list.files(ii, pattern = 'pball', include.dirs = T)))
  meteo <- read.csv(paste0(ii,'/', list.files(ii, pattern = 'NLDAS', include.dirs = T)))
  eg_nml <- read_nml(paste0(ii,'/', list.files(ii, pattern = 'nml', include.dirs = T)))
  H <- abs(eg_nml$morphometry$H - max(eg_nml$morphometry$H))
  A <- eg_nml$morphometry$A
  
  if (length( list.files(ii, pattern = 'wq_data', include.dirs = T)) > 0){
    wq_data<- paste0(ii,'/', list.files(ii, pattern = 'wq_data', include.dirs = T))
    obs <- NULL
    
    for(jj in wq_data){
      raw_obs <- read.csv(jj)
      secchi <- raw_obs %>%
        dplyr::filter(CharacteristicName== "Depth, Secchi disk depth") %>%
        dplyr::select(c('ActivityStartDate', 'ResultMeasureValue'))
    }
  }
  secchi.df <- data.frame('time' = secchi$ActivityStartDate, 'secchi' = secchi$ResultMeasureValue)
  
  secchi.df$year = year(secchi$ActivityStartDate)
  secchi.df$doy = yday(secchi$ActivityStartDate)
  annual.secchi.df <- secchi.df %>%
    dplyr::filter(doy >= 120 & doy <= 275) %>%
    group_by(year) %>%
    summarise(mean.secchi = mean(secchi, na.rm = TRUE))
  
  idx.years <- rep(NA, length(unique(year(meteo$time))))
  idx.years[which(unique(year(meteo$time)) %in% annual.secchi.df$year)] = annual.secchi.df$mean.secchi
  if (is.na(idx.years[1])){
    idx.years[1] = idx.years[!is.na(idx.years)][1]
  }
  if (is.na(idx.years[length(idx.years)])){
    idx.years[length(idx.years)] = rev(idx.years[!is.na(idx.years)])[1]
  }
  secchiapp <- na.approx(idx.years)
  
  if (length( list.files(ii, pattern = 'wq_data', include.dirs = T)) > 0){
    wq_data<- paste0(ii,'/', list.files(ii, pattern = 'wq_data', include.dirs = T))
    obs <- NULL
    
    for(jj in wq_data){
      raw_obs <- read.csv(jj)
      doc <- raw_obs %>%
        dplyr::filter(CharacteristicName== "Dissolved organic carbon") %>%
        dplyr::select(c('ActivityStartDate', 'ResultMeasureValue'))
    }
  }
  doc.df <- data.frame('time' = doc$ActivityStartDate, 'doc' = doc$ResultMeasureValue)
  doc.df$doc[which(doc.df$doc<= 0)] = 0
  
  doc.df$year = year(doc$ActivityStartDate)
  doc.df$doy = yday(doc$ActivityStartDate)
  annual.doc.df <- doc.df %>%
    dplyr::filter(doy >= 120 & doy <= 275) %>%
    group_by(year) %>%
    summarise(mean.doc = mean(doc, na.rm = TRUE))
  
  idx.years <- rep(NA, length(unique(year(meteo$time))))
  idx.years[which(unique(year(meteo$time)) %in% annual.doc.df$year)] = annual.doc.df$mean.doc
  if (is.na(idx.years[1])){
    idx.years[1] = idx.years[!is.na(idx.years)][1]
  }
  if (is.na(idx.years[length(idx.years)])){
    idx.years[length(idx.years)] = rev(idx.years[!is.na(idx.years)])[1]
  }
  docapp <- na.approx(idx.years)
  
  if (length( list.files(ii, pattern = 'wq_data', include.dirs = T)) > 0){
    wq_data<- paste0(ii,'/', list.files(ii, pattern = 'wq_data', include.dirs = T))
    obs <- NULL
    
    for(jj in wq_data){
      raw_obs <- read.csv(jj)
      phos <- raw_obs %>%
        dplyr::filter(CharacteristicName== "Total Phosphorous, filtered") %>%
        dplyr::select(c('ActivityStartDate', 'ResultMeasureValue'))
    }
  }
  phos.df <- data.frame('time' = phos$ActivityStartDate, 'phos' = phos$ResultMeasureValue)
  
  phos.df$year = year(phos$ActivityStartDate)
  phos.df$doy = yday(phos$ActivityStartDate)
  annual.phos.df <- phos.df %>%
    dplyr::filter(doy >= 120 & doy <= 275) %>%
    group_by(year) %>%
    summarise(mean.phos = mean(phos, na.rm = TRUE))
  
  idx.years <- rep(NA, length(unique(year(meteo$time))))
  idx.years[which(unique(year(meteo$time)) %in% annual.phos.df$year)] = annual.phos.df$mean.phos
  if (is.na(idx.years[1])){
    idx.years[1] = idx.years[!is.na(idx.years)][1]
  }
  if (is.na(idx.years[length(idx.years)])){
    idx.years[length(idx.years)] = rev(idx.years[!is.na(idx.years)])[1]
  }
  phosapp <- na.approx(idx.years)
  
  # secchiapp <-  Hmisc::approxExtrap(annual.wq.df$year, annual.wq.df$mean.secchi, unique(year(meteo$time)), 
  #                                   method = 'constant')$y
  # 
  
  
  if (a==1){
    df <- data.frame('year' = unique(year(meteo$time)),
                         'secchidep' = secchiapp,
                          'doc' = docapp,
                     'tp' = phosapp,
                         'area' = max(A),
                     'meandoc' = mean(docapp),
                     'id' = sub("\\).*", "", sub(".*\\(", "", ii))
)
    
    # all.dat <- data.frame('date' = wind.ts$datetime[match(data$date[1], meteo$time) : (length(wind.ts$datetime)-1)],
    #                       'jb' = (207 * 10^(-6) * 9.81)/(4180 *1000) * wind.ts$sw[match(data$date[1], meteo$time) : (length(wind.ts$datetime)-1)] / 1e3 * 1000,
    #                       'st' = E.schmidt$St,
    #                       'id' = sub("\\).*", "", sub(".*\\(", "", ii)),
    #                       'z' = rep(max(H)))
    # all.dat$lmo <- (sqrt(wind.ts$ux[match(data$date[1], meteo$time) : (length(wind.ts$datetime)-1)])^3 /all.dat$jb) / zv
    #
    # lmo.all <- all.dat
  } else {
    df <- rbind(df, data.frame('year' = unique(year(meteo$time)),
                                   'secchidep' = secchiapp,
                                   'doc' = docapp,
                                   'tp' = phosapp,
                                   'area' = max(A),
                               'meandoc' = mean(docapp),
                               'id' = sub("\\).*", "", sub(".*\\(", "", ii))))
    
    # all.dat <- data.frame('date' = wind.ts$datetime[match(data$date[1], meteo$time) : (length(wind.ts$datetime)-1)],
    #                       'jb' = (207 * 10^(-6) * 9.81)/(4180 *1000) * wind.ts$sw[match(data$date[1], meteo$time) : (length(wind.ts$datetime)-1)] / 1e3 * 1000,
    #                       'st' = E.schmidt$St,
    #                       'id' = sub("\\).*", "", sub(".*\\(", "", ii)),
    #                       'z' = rep(max(H)))
    # all.dat$lmo <- (sqrt(wind.ts$ux[match(data$date[1], meteo$time) : (length(wind.ts$datetime)-1)])^3 /all.dat$jb) / zv
    #
    # lmo.all <- rbind(lmo.all,all.dat )
  }
  a=a+1
  
}

ggplot(df, aes(tp, doc, col = id)) + geom_point() 

unique(df$id)
unique(df$meandoc)

df.lake <- data.frame('id' = unique(df$id), 'ratio' = unique(df$meandoc)/5.681246)

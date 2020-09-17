rm(list = ls())

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

## packages
library(devtools)
library(glmtools) 
library(dplyr)
library(ggplot2)
library(lubridate)
library(pracma)
library(readr)
library(LakeMetabolizer)
library(adagio)
library(zoo)
library(GenSA)
library(patchwork)

library(simpleAnoxia)

## load example data
lks <- list.dirs(path = 'inst/extdata/', full.names = TRUE, recursive = F)

for (ii in lks){
  print(paste0('Running ',ii))
  data <- read.csv(paste0(ii,'/', list.files(ii, pattern = 'pball', include.dirs = T)))
  meteo <- read.csv(paste0(ii,'/', list.files(ii, pattern = 'NLDAS', include.dirs = T)))
  
  chidx <- match(as.POSIXct(data$date),as.POSIXct(meteo$time))
  wind <- meteo$WindSpeed[chidx]
  airtemp <- meteo$AirTemp[chidx]
  
  if (length( list.files(ii, pattern = 'wq_data', include.dirs = T)) > 0){
    wq_data<- paste0(ii,'/', list.files(ii, pattern = 'wq_data', include.dirs = T))
    obs <- NULL
    
    for(jj in wq_data){
      raw_obs <- read.csv(jj)
      if ('ActivityDepthHeighMeasure.MeasureValue' %in% colnames(raw_obs)){
        wq <- raw_obs %>%
          dplyr::filter(CharacteristicName== "Dissolved oxygen (DO)") %>%
          dplyr::select(c('ActivityStartDate', 'ActivityDepthHeighMeasure.MeasureValue', 'ResultMeasureValue'))
        wq <- rename(wq, 'ActivityDepthHeightMeasure.MeasureValue' = 'ActivityDepthHeighMeasure.MeasureValue')
      } else {
        wq <- raw_obs %>%
          dplyr::filter(CharacteristicName== "Dissolved oxygen (DO)") %>%
          dplyr::select(c('ActivityStartDate', 'ActivityDepthHeightMeasure.MeasureValue', 'ResultMeasureValue'))
      }
      
      if (length(wq_data) == 1 | is.null(obs)){
        obs <- wq
      } else {
        obs <- rbind(obs, wq)
      }
    }
    obs$ActivityStartDate<-as.POSIXct(obs$ActivityStartDate)
  }
  
  if (is.factor(obs$ActivityDepthHeightMeasure.MeasureValue)){
    obs$ActivityDepthHeightMeasure.MeasureValue <-  as.numeric(as.character(obs$ActivityDepthHeightMeasure.MeasureValue))
  }
  if (is.factor(obs$ResultMeasureValue)){
    obs$ResultMeasureValue <-  as.numeric(as.character(obs$ResultMeasureValue))
  }
  
  # outlier detection
  outlier_values <- boxplot.stats(obs$ResultMeasureValue)$out 
  uvx <- match(outlier_values, obs$ResultMeasureValue)
  obs$ResultMeasureValue[uvx] <- NA
  
  eg_nml <- read_nml(paste0(ii,'/', list.files(ii, pattern = 'nml', include.dirs = T)))
  H <- abs(eg_nml$morphometry$H - max(eg_nml$morphometry$H))
  A <- eg_nml$morphometry$A
  
  input.values <- input(wtemp = data, H = H, A = A)
  
  input.values$year <- year(input.values$datetime)
  input.values$doy <- yday(input.values$datetime)
  input.values$max.d <- max(H)
  
  input.values$wind = wind
  input.values$airtemp = airtemp
  write.table(input.values, paste0('input_',sub("\\).*", "", sub(".*\\(", "", ii)) ,'.txt'), append = FALSE, sep = ",", dec = ".",
              row.names = FALSE, col.names = FALSE)
  
  # proc.obs <- preprocess_obs(obs,input.values = input.values, H, A)
  w.obs <- weigh_obs(obs,input.values = input.values, H, A)
  obs_long <- w.obs[[1]]
  obs_weigh <- w.obs[[2]]

  write.table(obs_weigh,paste0('oxy_observed_',sub("\\).*", "", sub(".*\\(", "", ii)) ,'.txt'), append = FALSE, sep = " ", dec = ".",
              row.names = FALSE, col.names = FALSE)
  write.table(obs_long,paste0('oxy_observed_long_',sub("\\).*", "", sub(".*\\(", "", ii)) ,'.txt'), append = FALSE, sep = " ", dec = ".",
              row.names = FALSE, col.names = FALSE)
  
  w.obs <- weigh_obs(obs,input.values = input.values, H, A)
  obs_long <- w.obs[[1]]
  obs_weigh <- w.obs[[2]]
  
  get_wq_data <- function(kword){
    if (length( list.files(ii, pattern = 'wq_data', include.dirs = T)) > 0){
      wq_data<- paste0(ii,'/', list.files(ii, pattern = 'wq_data', include.dirs = T))
      obs <- NULL
      
      for(jj in wq_data){
        raw_obs <- read.csv(jj)
        if ('ActivityDepthHeighMeasure.MeasureValue' %in% colnames(raw_obs)){
          wq <- raw_obs %>%
            dplyr::filter(CharacteristicName== kword) %>%
            dplyr::select(c('ActivityStartDate', 'ActivityDepthHeighMeasure.MeasureValue', 'ResultMeasureValue'))
          wq <- rename(wq, 'ActivityDepthHeightMeasure.MeasureValue' = 'ActivityDepthHeighMeasure.MeasureValue')
        } else {
          wq <- raw_obs %>%
            dplyr::filter(CharacteristicName== kword) %>%
            dplyr::select(c('ActivityStartDate', 'ActivityDepthHeightMeasure.MeasureValue', 'ResultMeasureValue'))
        }
        
        if (length(wq_data) == 1 | is.null(obs)){
          obs <- wq
        } else {
          obs <- rbind(obs, wq)
        }
      }
      obs$ActivityStartDate<-as.POSIXct(obs$ActivityStartDate)
    }
    
    if (is.factor(obs$ActivityDepthHeightMeasure.MeasureValue)){
      obs$ActivityDepthHeightMeasure.MeasureValue <-  as.numeric(as.character(obs$ActivityDepthHeightMeasure.MeasureValue))
    }
    if (is.factor(obs$ResultMeasureValue)){
      obs$ResultMeasureValue <-  as.numeric(as.character(obs$ResultMeasureValue))
    }
    
    # outlier detection
    outlier_values <- boxplot.stats(obs$ResultMeasureValue)$out 
    uvx <- match(outlier_values, obs$ResultMeasureValue)
    obs$ResultMeasureValue[uvx] <- NA
    
    w.obs <- weigh_obs(obs,input.values = input.values, H, A)
    obs_long <- w.obs[[1]]
    obs_weigh <- w.obs[[2]]
    
    write.table(obs_weigh,paste0(kword,'_observed_',sub("\\).*", "", sub(".*\\(", "", ii)) ,'.txt'), append = FALSE, sep = " ", dec = ".",
                row.names = FALSE, col.names = FALSE)
    write.table(obs_long,paste0(kword,'_observed_long_',sub("\\).*", "", sub(".*\\(", "", ii)) ,'.txt'), append = FALSE, sep = " ", dec = ".",
                row.names = FALSE, col.names = FALSE)
    return(NULL)
  }
  unique(raw_obs$CharacteristicName)
  get_wq_data(kword = "Dissolved organic carbon")
  get_wq_data(kword = "Dissolved inorganic carbon")
  get_wq_data(kword = "Total organic carbon")
  get_wq_data(kword = "Total Nitrogen, filtererd")
  get_wq_data(kword = "Total Phosphorous, filtered")
  get_wq_data(kword = "Dissolved Reactive Phosphorus, lab")
  get_wq_data(kword = "Chlorophyll a, free of pheophytin")
  get_wq_data(kword = "Depth, Secchi disk depth")
  get_wq_data(kword = "Inorganic nitrogen (nitrate and nitrite) as N, lab")
  

  print('Got that data, moving on.')
}

for (ii in lks){
  print(paste0('Running ',ii))
  data <- read.csv(paste0(ii,'/', list.files(ii, pattern = 'pball', include.dirs = T))) # GLM sim water temp
  meteo <- read.csv(paste0(ii,'/', list.files(ii, pattern = 'NLDAS', include.dirs = T))) # meteo NLDAS
  
  chidx <- match(as.POSIXct(data$date),as.POSIXct(meteo$time))
  wind <- meteo$WindSpeed[chidx]
  
  if (length( list.files(ii, pattern = 'wq_data', include.dirs = T)) > 0){ # wq LTER
  wq_data<- paste0(ii,'/', list.files(ii, pattern = 'wq_data', include.dirs = T))
  obs <- NULL

  for(jj in wq_data){
      raw_obs <- read.csv(jj)
      if ('ActivityDepthHeighMeasure.MeasureValue' %in% colnames(raw_obs)){
        wq <- raw_obs %>%
          dplyr::filter(CharacteristicName== "Dissolved oxygen (DO)") %>%
          dplyr::select(c('ActivityStartDate', 'ActivityDepthHeighMeasure.MeasureValue', 'ResultMeasureValue'))
       wq <- rename(wq, 'ActivityDepthHeightMeasure.MeasureValue' = 'ActivityDepthHeighMeasure.MeasureValue')
      } else {
        wq <- raw_obs %>%
          dplyr::filter(CharacteristicName== "Dissolved oxygen (DO)") %>%
          dplyr::select(c('ActivityStartDate', 'ActivityDepthHeightMeasure.MeasureValue', 'ResultMeasureValue'))
      }
      
      if (length(wq_data) == 1 | is.null(obs)){
        obs <- wq
      } else {
        obs <- rbind(obs, wq)
      }
  }
  obs$ActivityStartDate<-as.POSIXct(obs$ActivityStartDate)
  }

  if (is.factor(obs$ActivityDepthHeightMeasure.MeasureValue)){
    obs$ActivityDepthHeightMeasure.MeasureValue <-  as.numeric(as.character(obs$ActivityDepthHeightMeasure.MeasureValue))
  }
  if (is.factor(obs$ResultMeasureValue)){
    obs$ResultMeasureValue <-  as.numeric(as.character(obs$ResultMeasureValue))
  }
  
  # outlier detection
  outlier_values <- boxplot.stats(obs$ResultMeasureValue)$out 
  uvx <- match(outlier_values, obs$ResultMeasureValue)
  obs$ResultMeasureValue[uvx] <- NA
  
  eg_nml <- read_nml(paste0(ii,'/', list.files(ii, pattern = 'nml', include.dirs = T))) # NML GLM
  H <- abs(eg_nml$morphometry$H - max(eg_nml$morphometry$H)) # DEPTH
  A <- eg_nml$morphometry$A # AREA

  input.values <- input(wtemp = data, H = H, A = A)
  
  input.values$year <- year(input.values$datetime)
  input.values$doy <- yday(input.values$datetime)
  
  input.values$wind = wind
  write.table(input.values, paste0('input.txt'), append = FALSE, sep = ",", dec = ".",
              row.names = FALSE, col.names = FALSE) # input data --> B-ODEM
  
  ggplot(input.values)+
    geom_line(aes(doy, td.depth, col='td.depth')) +
    geom_line(aes(doy,  upper.metalim, col =  'upper.metalim')) +
    geom_line(aes(doy,  lower.metalim, col =  'lower.metalim')) +
    facet_wrap(~year)
  # ggplot(input.values)+
  #   geom_line(aes(doy, vol_epil, col='t.epil')) +
  #   facet_wrap(~year)


  # data[,1] <- as.POSIXct(as.character(data[,1]))
  # filled.contour(x = data[,1]
  #                , y = seq(0, 22.5, 0.5)
  #                , as.matrix(data[,2:ncol(data)]))
  
  # proc.obs <- preprocess_obs(obs,input.values = input.values, H, A)
  w.obs <- weigh_obs(obs,input.values = input.values, H, A)
  obs_long <- w.obs[[1]]
  obs_weigh <- w.obs[[2]]
  
  ggplot(obs_long) +
    geom_freqpoly(aes(ResultMeasureValue , col = Layer),binwidth = 1)
  write.table(obs_weigh, paste0('observed.txt'), append = FALSE, sep = " ", dec = ".",
              row.names = FALSE, col.names = FALSE)

  nep_stratified = 100 # mg/m3/d
  min_stratified = -50 #-200 # mg/m3/d
  min_not_stratified  =  0 # mg/m2/d
  fsed_stratified_epi = -1500 # mg/m3/d
  nep_not_stratified = 0 # mg/m3/d
  fsed_stratified_hypo = -4500 # mg/m3/d
  fsed_not_stratified = 0 # mg/m3/d
  khalf <- 3000 # mg/m3
  # 
  # startdate = 1
  # enddate = nrow(input.values)
  # o2<- calc_do(input.values = input.values,
  #              fsed_stratified_epi = fsed_stratified_epi,
  #              fsed_stratified_hypo = fsed_stratified_hypo,
  #              fsed_not_stratified = fsed_not_stratified,
  #              nep_stratified = nep_stratified,
  #              nep_not_stratified = nep_not_stratified,
  #              min_stratified = min_stratified,
  #              min_not_stratified =min_not_stratified,
  #              wind = wind, 
  #              khalf = khalf,
  #              startdate = startdate,
  #              enddate = enddate)

  
  #mass balance code
  # o2 <- o2[startdate: enddate,]
  # o2$time <- input.values$datetime[startdate:enddate]
  # print(paste0('mass balance kept by ',round((100*sum(abs(o2$massbal) <1, na.rm = TRUE))/ sum(!is.na(o2$massbal)), 0),'% of the data'))
  # plot(o2$o2_epil)
  # plot(o2$o2_hypo)
  # t1 <- ggplot(o2, aes(time, Fsed_epi, col = Fsed_epi)) + geom_point()
  # t2 <- ggplot(o2, aes(time, NEP_epi, col = NEP_epi)) + geom_point()
  # t3 <- ggplot(o2, aes(time, Fatm_epi, col = Fatm_epi)) + geom_point()
  # t4 <- ggplot(o2, aes(time, Entrain_epi, col = Entrain_epi)) + geom_point()
  # t5 <- ggplot(o2, aes(time, Fsed_hypo, col = Fsed_hypo)) + geom_point()
  # t6 <- ggplot(o2, aes(time, Mineral_hypo, col = Mineral_hypo)) + geom_point()
  # t7 <- ggplot(o2, aes(time, Entrain_hypo, col = Entrain_hypo)) + geom_point()
  # t8 <- ggplot(o2, aes(time, o2_epil, col = o2_epil)) + geom_point()
  # t9 <- ggplot(o2, aes(time, o2_hypo, col = o2_hypo)) + geom_point()
  # t1 + t2 + t3 + t4 +t5 + t6+ t7 + t8 + t9
  # 
  # plot(input.values$datetime[startdate: enddate],o2$massbal[startdate: enddate]/input.values$vol_total[startdate: enddate]/1000,
  #      ylab = 'Mass balance O2 difference (mg/L)',
  #      xlab = '')
  # plot(input.values$datetime[startdate: enddate],o2$massbal[startdate: enddate],
  #      ylab = 'Mass balance O2 difference (mg)',
  #      xlab = '')
  # boxplot(o2$massbal, ylab = 'Mass balance O2 difference (mg/L')
  # sum(o2$massbal/input.values$vol_total[startdate: enddate]/1000, na.rm = TRUE)
  # 
  # all_fluxes <- cbind((o2$NEP_total + o2$Fatm_total + o2$Mineral_total + o2$Fsed_total) * input.values$vol_total[startdate: enddate],
  #                     (o2$NEP_epi + o2$Fatm_epi + o2$Entrain_epi + o2$Fsed_epi) * input.values$vol_epil[startdate: enddate],
  #                     (o2$Mineral_hypo + o2$Entrain_hypo + o2$Fsed_hypo) * input.values$vol_hypo[startdate: enddate])
  # sum.all_fluxes <- apply(all_fluxes, 1, sum, na.rm= TRUE)
  # 
  # all_concfluxes <- cbind((o2$NEP_total + o2$Fatm_total + o2$Mineral_total + o2$Fsed_total),
  #                         (o2$NEP_epi + o2$Fatm_epi + o2$Entrain_epi + o2$Fsed_epi),
  #                         (o2$Mineral_hypo + o2$Entrain_hypo + o2$Fsed_hypo) )
  # sum.all_concfluxes <- apply(all_concfluxes, 1, sum, na.rm= TRUE)
  # 
  # df.msblc <- data.frame('time' = input.values$datetime[startdate: enddate], 'cumsum_do' = cumsum(c(0,diff(o2$o2_total * input.values$vol_total[startdate: enddate]))),
  #                        'cumsum_fluxes' = cumsum(sum.all_fluxes))
  # ggplot(df.msblc, aes((cumsum_do), (cumsum_fluxes), col = time)) + geom_point()
  # 
  # plot((diff(o2$o2_total * input.values$vol_total[startdate: enddate])), ylim = c(0,5e11))
  # points((sum.all_fluxes), col = 'red')
  # 
  # 
  # plot(c(0,diff(o2$o2_total * input.values$vol_total[startdate: enddate])) - sum.all_fluxes)
  # plot(c(0,diff(o2$o2_total * input.values$vol_total[startdate: enddate])) - sum.all_fluxes, ylim = c(-1e13, 1e13))
  # 
  # sum(c(0,diff(o2$o2_total * input.values$vol_total[startdate: enddate]))) - sum(sum.all_fluxes)
  # sum(abs(c(0,diff(o2$o2_total * input.values$vol_total[startdate: enddate])))) - sum(abs(sum.all_fluxes))
  # sum(abs(c(0,diff(o2$o2_total)))) - sum(abs(sum.all_concfluxes))

  
  # fsed_stratified_epi fsed_stratified_hypo nep_stratified min_stratified, khalf
  init.val = c(5, 5, 5, 5, 5)
  target.iter = 60
  lb <<- c(-100, -100, 10, -1000, 1000)
  ub <<- c(-6000, -6000, 1000, 1000, 7000) 

  # # calibration-validation
  val.ratio <- 1/3
  cal.ratio <- 1 - val.ratio

  val_data <- obs_weigh[, 1 : round(ncol(obs_weigh) * val.ratio)]
  cal_data <- obs_weigh[, ((round(ncol(obs_weigh) * val.ratio))) : ncol(obs_weigh)]
  # 
  datval = data.frame('id' = ii, 'start' =input.values$datetime[obs_weigh[1, (1 + (round(ncol(obs_weigh) * val.ratio)))]])
  write.table(datval, file ='calval.csv', append=TRUE, quote = FALSE, col.names = FALSE,
              row.names = FALSE)
  # 
  modelopt <- pureCMAES(par = init.val, fun = optim_do, lower = rep(0,5),
                        upper = rep(10,5), sigma = 0.5,
                        stopfitness = -Inf,
                        stopeval = target.iter,
                        input.values = input.values,
                        fsed_not_stratified = fsed_not_stratified,
                        nep_not_stratified = nep_not_stratified, min_not_stratified = min_not_stratified,
                        wind = wind, proc.obs = cal_data, # caL_data
                        verbose = verbose, startdate = NULL, enddate = NULL)
  
  print(modelopt$fmin)
  print(lb+(ub - lb)/(10)*(modelopt$xmin))
  modelopt$xmin <- lb+(ub - lb)/(10)*(modelopt$xmin)
  
  o2<- calc_do(input.values = input.values,fsed_stratified_epi = modelopt$xmin[1],
               fsed_stratified_hypo = modelopt$xmin[2],
               fsed_not_stratified,
               nep_stratified = modelopt$xmin[3],
               nep_not_stratified,
               min_stratified = modelopt$xmin[4],
               min_not_stratified, wind,
               startdate = NULL, enddate = NULL)
  
  input.values$o2_epil <- o2[,"o2_epil"]
  input.values$o2_hypo <- o2[,"o2_hypo"]
  input.values$o2_total <- o2[,"o2_total"]
  input.values$sat_o2_epil <- o2[,"sat_o2_epil"]
  input.values$sat_o2_hypo <- o2[,"sat_o2_hypo"]
  input.values$sat_o2_total <- o2[,"sat_o2_total"]
  #'Fsed_total', 'NEP_total', 'Fatm_total', 'Mineral_total',
  #'Fsed_epi', 'NEP_epi', 'Fatm_epi', 'Entrain_epi', 'Fsed_hypo', 'Mineral_hypo', 'Entrain_hypo')
  input.values$Fsed_total <- o2[,"Fsed_total"]
  input.values$NEP_total <- o2[,"NEP_total"]
  input.values$Fatm_total <- o2[,"Fatm_total"]
  input.values$Mineral_total <- o2[,"Mineral_total"]
  input.values$Fsed_epi <- o2[,"Fsed_epi"]
  input.values$NEP_epi <- o2[,"NEP_epi"]
  input.values$Fatm_epi <- o2[,"Fatm_epi"]
  input.values$Entrain_epi <- o2[,"Entrain_epi"]
  input.values$Fsed_hypo <- o2[,"Fsed_hypo"]
  input.values$Mineral_hypo <- o2[,"Mineral_hypo"]
  input.values$Entrain_hypo <- o2[,"Entrain_hypo"]
  
  all_params <- c()
  all_o2 <- c()
  target.iter = 60
  
  # for loop for single-year cal
  # for (ikx in unique(year(obs_long$ActivityStartDate ))){
  #   kind <- which(ikx == input.values$year)
  #   datfm <- input.values$datetime[kind]
  #   if (sum(!is.na(input.values$td.depth[kind])) < 7 | (ikx %in% year(input.values$datetime[obs_weigh[1,]]) == FALSE) |
  #       !any(datfm[which(!is.na(input.values$td.depth[kind]))] %in% input.values$datetime[obs_weigh[1,]]) ){
  #     next
  #   } 
  #   startdate <- min(kind)
  #   enddate <- max(kind)
  #   
  #   modelopt <- pureCMAES(par = init.val, fun = optim_do, lower = rep(0,5),
  #                         upper = rep(10,5), sigma = 0.5,
  #                         stopfitness = -Inf,
  #                         stopeval = target.iter,
  #                         input.values = input.values,
  #                         fsed_not_stratified = fsed_not_stratified,
  #                         nep_not_stratified = nep_not_stratified, min_not_stratified = min_not_stratified,
  #                         wind = wind, proc.obs = obs_weigh, # caL_data
  #                         verbose = verbose, startdate = startdate, enddate = enddate)
  #   
  #   modelopt$xmin <- lb+(ub - lb)/(10)*(modelopt$xmin)
  #   
  #   calibrated_param <- c(ikx, modelopt$xmin, modelopt$fmin)
  #   all_params <- rbind(all_params, calibrated_param)
  #   
  #   
  #   o2<- calc_do(input.values = input.values,fsed_stratified_epi = modelopt$xmin[1],
  #                fsed_stratified_hypo = modelopt$xmin[2],
  #                fsed_not_stratified,
  #                nep_stratified = modelopt$xmin[3],
  #                nep_not_stratified,
  #                min_stratified = modelopt$xmin[4],
  #                min_not_stratified, wind,
  #                startdate = startdate, enddate = enddate)
  #   
  #   input.values$o2_epil[startdate:enddate] <- o2[startdate:enddate,"o2_epil"]
  #   input.values$o2_hypo[startdate:enddate] <- o2[startdate:enddate,"o2_hypo"]
  #   input.values$o2_total[startdate:enddate] <- o2[startdate:enddate,"o2_total"]
  #   input.values$sat_o2_epil[startdate:enddate] <- o2[startdate:enddate,"sat_o2_epil"]
  #   input.values$sat_o2_hypo[startdate:enddate] <- o2[startdate:enddate,"sat_o2_hypo"]
  #   input.values$sat_o2_total[startdate:enddate] <- o2[startdate:enddate,"sat_o2_total"]
  #   input.values$Fsed_total[startdate:enddate] <- o2[startdate:enddate,"Fsed_total"]
  #   input.values$NEP_total[startdate:enddate] <- o2[startdate:enddate,"NEP_total"]
  #   input.values$Fatm_total[startdate:enddate] <- o2[startdate:enddate,"Fatm_total"]
  #   input.values$Mineral_total[startdate:enddate] <- o2[startdate:enddate,"Mineral_total"]
  #   input.values$Fsed_epi[startdate:enddate] <- o2[startdate:enddate,"Fsed_epi"]
  #   input.values$NEP_epi[startdate:enddate] <- o2[startdate:enddate,"NEP_epi"]
  #   input.values$Fatm_epi[startdate:enddate] <- o2[startdate:enddate,"Fatm_epi"]
  #   input.values$Entrain_epi[startdate:enddate] <- o2[startdate:enddate,"Entrain_epi"]
  #   input.values$Fsed_hypo[startdate:enddate] <- o2[startdate:enddate,"Fsed_hypo"]
  #   input.values$Mineral_hypo[startdate:enddate] <- o2[startdate:enddate,"Mineral_hypo"]
  #   input.values$Entrain_hypo[startdate:enddate] <- o2[startdate:enddate,"Entrain_hypo"]
  #   
  #   all_o2 <- rbind(all_o2, o2[startdate:enddate,])
  # }
  # 
  # m.all_params <- as.data.frame(all_params)
  # colnames(m.all_params) <- c('year', 'Fsed_epi', 'Fsed_hypo', 'NEP', 'Mineral',
  #                           'khalf', 'RMSE')
  # i1 <- ggplot(m.all_params, aes(year, RMSE)) +
  #   geom_line() + geom_hline(yintercept=2.12)
  # i2 <- ggplot(m.all_params, aes(year, Fsed_epi)) +
  #   geom_line() + geom_hline(yintercept=4536.41)
  # i3 <- ggplot(m.all_params, aes(year, Fsed_hypo)) +
  #   geom_line() + geom_hline(yintercept=154.85)
  # i4 <- ggplot(m.all_params, aes(year, NEP)) +
  #   geom_line() + geom_hline(yintercept=10)
  # i5 <- ggplot(m.all_params, aes(year, Mineral)) +
  #   geom_line() + geom_hline(yintercept=-1000)
  # i6 <- ggplot(m.all_params, aes(year, khalf)) +
  #   geom_line() + geom_hline(yintercept=6591.92)
  # 
  # p <- i1 + i2 + i3 + i4 + i5 + i6 & theme_minimal()
  # p 
  # ggsave(file = paste0(ii,'/comparison_indyearcal.png'), p, dpi=300,width=316,height=190,units='mm')
  # 
  
  # input.values$o2_epil <- o2[,"o2_epil"]
  # input.values$o2_hypo <- o2[,"o2_hypo"]
  # input.values$o2_total <- o2[,"o2_total"]
  # input.values$sat_o2_epil <- o2[,"sat_o2_epil"]
  # input.values$sat_o2_hypo <- o2[,"sat_o2_hypo"]
  # input.values$sat_o2_total <- o2[,"sat_o2_total"]

  
  test_data<-compare_predict_versus_observed(obs,input.values) 
  test_data$year <- year(test_data$day)
  test_data$doy <- yday(test_data$day)
  
  fit_val = calc_fit(input.values , proc.obs = val_data )
  print(paste0('validation error: ', round(fit_val,2)))
  fit_cal = calc_fit(input.values , proc.obs = cal_data )
  print(paste0('calibration error: ', round(fit_cal,2)))
  fit_total = calc_fit(input.values , proc.obs = obs_weigh )
  print(paste0('total error: ', round(fit_total,2)))
  
  pgm <- input.values %>%
    dplyr::select(datetime, o2_epil, o2_hypo, o2_total, vol_epil, vol_hypo, 'vol_total' = vol_total)
  input.pgm <- input.values %>%
    dplyr::select(datetime, "depth_td" = td.depth, "area_td" = td_area, "area_surf" = surf_area,
                  "wtemp_total" = t.total,
                  'wtemp_epil' = t.epil, "wtemp_hypo" = t.hypo,
                  'vol_total' = vol_total, vol_epil, vol_hypo)
  o2$datetime <- input.values$datetime
  fluxes <- input.values %>%
    dplyr::select('datetime', 'Fsed_total', 'NEP_total', 'Fatm_total', 'Mineral_total',
                  'Fsed_epi', 'NEP_epi', 'Fatm_epi', 'Entrain_epi', 'Fsed_hypo', 'Mineral_hypo', 'Entrain_hypo')
  fluxes[,-c(1)] <- fluxes[,-c(1)]#/1000/input.pgm$area_surf
  calibrated_param <- as.matrix(modelopt$xmin, byrow = TRUE)
  rownames(calibrated_param) <- c('fsed_stratified_epi', 'fsed_stratified_hypo', 'nep_stratified', 'min_stratified', 'khalf')
  
  write_delim(input.pgm, path = paste0(ii,'/',sub("\\).*", "","ODEM_", sub(".*\\(", "", ii)) ,'_',round(fit_total,1),'_alldata.txt'), delim = '\t')
  write_delim(pgm, path = paste0(ii,'/',sub("\\).*", "","ODEM_",  sub(".*\\(", "", ii)) ,'_',round(fit_total,1),'_oxymodel.txt'), delim = '\t')
  write_delim(obs_long, path = paste0(ii,'/',sub("\\).*", "","ODEM_",  sub(".*\\(", "", ii)) ,'_obs.txt'), delim = '\t')
  write_delim(fluxes, path = paste0(ii,'/',sub("\\).*", "","ODEM_",  sub(".*\\(", "", ii)) ,'_fluxes.txt'), delim = '\t')
  write_delim(as.data.frame(calibrated_param), path = paste0(ii,'/',sub("\\).*", "","ODEM_",  sub(".*\\(", "", ii)) ,'_calibparams.txt'), delim = '\t')
  
  observed <- data.frame('time' = input.values$datetime[obs_weigh[1,]],
                            'year' = input.values$year[obs_weigh[1,]],
                            'doy' = input.values$doy[obs_weigh[1,]],
                            'epi' =obs_weigh[3,],
                            'hypo' = obs_weigh[4,],
                         'epi_sim' = input.values$o2_epil[obs_weigh[1,]]/input.values$vol_epil[obs_weigh[1,]]/1000,
                         'hypo_sim' = input.values$o2_hypo[obs_weigh[1,]]/input.values$vol_hypo[obs_weigh[1,]]/1000)
  
  g1 <- ggplot(input.values) +
    geom_point(aes(doy, (o2_total/1000), col = 'Total')) +
    geom_point(aes(doy, (o2_epil/1000), col = 'Epi')) +
    geom_point(aes(doy, (o2_hypo/1000), col = 'Hypo')) +
    # ylim(0,20)+
    facet_wrap(~year) +
    theme_bw() +
    geom_point(data = observed, aes(doy, epi, col = 'Obs_Epi'), size =2, alpha = 0.5) +
    geom_point(data = observed, aes(doy, hypo, col = 'Obs_Hypo'), size =2, alpha = 0.5);g1
  ggsave(file = paste0(ii,'/ODEM_oxymodel.png'), g1, dpi=300,width=316,height=190,units='mm');g1


  
  fluxes$doy <- yday(fluxes$datetime)
  fluxes$year <- year(fluxes$datetime)
  

  g23 <- ggplot(fluxes) +
    geom_point(aes(doy, Fsed_epi, col = 'Fsed_epi')) +
    geom_point(aes(doy, NEP_epi, col = 'NEP_epi')) +
    geom_point(aes(doy, Fatm_epi, col = 'Fatm_epi')) +
    # geom_point(aes(doy, Entrain_epi, col = 'Entrain_epi')) +
    geom_point(aes(doy, Fsed_hypo, col = 'Fsed_hypo')) +
    geom_point(aes(doy, Mineral_hypo, col = 'Mineral_hypo')) +
    # geom_point(aes(doy, Entrain_hypo, col = 'Entrain_hypo')) +
    # ylim(0,20)+
    facet_wrap(~year) +
    theme_bw();g23
  ggsave(file = paste0(ii,'/ODEM_fluxes.png'), g23, dpi=300,width=316,height=190,units='mm')

  

  g3 <- ggplot(observed, aes(time, epi - epi_sim, col = 'epilimnion')) +
    geom_point(size = 2) +
    geom_line() +
    geom_point(aes(time, hypo - hypo_sim, col = 'hypolimnion'), size =2) +
    geom_line(aes(time, hypo - hypo_sim, col = 'hypolimnion'))+
    ylab('Residuals (obs-mod)')+
    theme_bw();g3
  ggsave(file = paste0(ii,'/ODEM_predicted_ag_observed.png'), g3, dpi=300, width=316,height=190,units='mm')
  
  
 g10 <- ggplot(input.values) +
   geom_point(aes(doy, (sat_o2_total), col = 'Total')) +
   geom_point(aes(doy, (sat_o2_epil), col = 'Epi')) +
   geom_point(aes(doy, (sat_o2_hypo), col = 'Hypo')) +
   ylim(0,150)+
   facet_wrap(~year) +
   theme_bw();g10
 ggsave(file = paste0(ii,'/ODEM_sat_o2.png'), g10, dpi=300,width=316,height=190,units='mm')
 

  eval.info <- data.frame('id' = sub("\\).*", "", sub(".*\\(", "", ii)) ,
                          'time' = Sys.time(),
             'A' = max(A),
             'z' = max(H),
             'obs' = ncol(obs_weigh),
             'RMSE_cal' = fit_cal,
             'RMSE_val' = fit_val,
             'RMSE_total' = fit_total)
    write.table(eval.info, file ='eval.csv', append=TRUE, quote = FALSE, col.names = FALSE,
                row.names = FALSE)
  print('Nothing got broken, good job!')
}

library(viridis)
library(patchwork)
eval.df <- read.table('eval.csv', header = TRUE)

g1<- ggplot(subset(eval.df,gen=='threequarter'), aes(MaxZ, RMSE_total, col = (Asurf), label = ID)) + 
  geom_point() +   
  scale_color_viridis(option="viridis") +
  xlab('Depth') +
  ylab('RMSE in mg DO/L') + 
  ggtitle("total") +
  ylim(0,3)+ xlim(0,35)+
  geom_text(check_overlap = FALSE,hjust = 1.1, nudge_x = 0.09) + 
  theme_bw();g1
g2<- ggplot(subset(eval.df,gen=='threequarter'), aes(MaxZ, RMSE_cal, col = Asurf, label = ID)) +
  geom_point() +
  scale_color_viridis(option="viridis") +
  xlab('Depth') +
  ylab('RMSE in mg DO/L') +
  ggtitle("calibration") +
  ylim(0,3)+ xlim(0,35)+
  geom_text(check_overlap = FALSE,hjust = 1.1, nudge_x = 0.05) +
  theme_bw();g2
g3<- ggplot(subset(eval.df,gen=='threequarter'), aes(MaxZ, RMSE_val, col = Asurf, label = ID)) +
  geom_point() + #aes(size = nobs)
  scale_color_viridis(option="viridis") +
  xlab('Depth') +
  ylab('RMSE in mg DO/L') +
  ggtitle("validation") +
  ylim(0,3)+ xlim(0,35)+
  geom_text(check_overlap = FALSE,hjust = 1.1, nudge_x = 0.05) +
  theme_bw();g3
g <- g1 / g2 / g3 + plot_annotation(tag_levels = 'A',
  caption = paste0('ODEM ',Sys.time())
);g
ggsave(file = paste0('lake_results.png'), g, dpi=300, width=175,height=275,units='mm')

h<- g1 + g2 + g3 + plot_annotation(tag_levels = 'A',
                                    caption = paste0('ODEM ',Sys.time())
);h
ggsave(file = paste0('lake_results_talk.png'), h, dpi=300, width=400,height=150,units='mm')

g4<- ggplot(subset(eval.df,gen=='RNN'), aes(MaxZ, RMSE_cal, col = Asurf, label = ID)) +
  geom_point() +
  scale_color_viridis(option="viridis") +
  xlab('Depth') +
  ylab('RMSE in mg DO/L') +
  ggtitle("RNN calibration") +
  ylim(0,3)+ xlim(0,35)+
  geom_text(check_overlap = FALSE,hjust = 1.1, nudge_x = 0.05) +
  theme_bw();g4

h2 <- g2 + g4 + plot_annotation(tag_levels = 'A',
                                caption = paste0('ODEM ',Sys.time())
);h2
ggsave(file = paste0('lake_results_RNN.png'), h2, dpi=300, width=400,height=150,units='mm')


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

library(simpleAnoxia)

## load example data
lks <- list.dirs(path = 'inst/extdata/', full.names = TRUE, recursive = F)

for (ii in lks[c(9, 10, 11, 13, 14, 15, 16, 17, 8)]){
  print(paste0('Running ',ii))
  data <- read.csv(paste0(ii,'/', list.files(ii, pattern = 'temperatures.csv', include.dirs = T)))
  meteo <- read.csv(paste0(ii,'/', list.files(ii, pattern = 'NLDAS', include.dirs = T)))
  
  chidx <- match(as.POSIXct(data$date),as.POSIXct(meteo$time))
  wind <- meteo$WindSpeed[chidx]
  
  if (length( list.files(ii, pattern = 'wq_data', include.dirs = T)) > 0){
  wq_data<- paste0(ii,'/', list.files(ii, pattern = 'wq_data', include.dirs = T))
  obs <- NULL
  # obs <- data.frame() 
  for(jj in wq_data){
      raw_obs <- read.csv(jj)
      # filter_obs<-filter(raw_obs,CharacteristicName== "Dissolved oxygen (DO)")
      # more_obs <- filter_obs %>%
        # dplyr::select(c('ActivityStartDate', 'ActivityStartTime.Time', 'ActivityDepthHeightMeasure.MeasureValue', 'ResultMeasureValue'))
      
      if ('ActivityDepthHeighMeasure.MeasureValue' %in% colnames(raw_obs)){
        wq <- raw_obs %>%
          dplyr::filter(CharacteristicName== "Dissolved oxygen (DO)") %>%
          # dplyr::select(c('ActivityStartDate', 'ActivityStartTime.Time', 'ActivityDepthHeightMeasure.MeasureValue', 'ResultMeasureValue'))
          dplyr::select(c('ActivityStartDate', 'ActivityDepthHeighMeasure.MeasureValue', 'ResultMeasureValue'))
       wq <- rename(wq, 'ActivityDepthHeightMeasure.MeasureValue' = 'ActivityDepthHeighMeasure.MeasureValue')
      } else {
        wq <- raw_obs %>%
          dplyr::filter(CharacteristicName== "Dissolved oxygen (DO)") %>%
          # dplyr::select(c('ActivityStartDate', 'ActivityStartTime.Time', 'ActivityDepthHeightMeasure.MeasureValue', 'ResultMeasureValue'))
          dplyr::select(c('ActivityStartDate', 'ActivityDepthHeightMeasure.MeasureValue', 'ResultMeasureValue'))
        
      }
      
      
      if (length(wq_data) == 1 | is.null(obs)){
        obs <- wq
      } else {
        obs <- rbind(obs, wq)
      }
      # obs<-rbind(obs,more_obs)
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

  # here's the promised input function (see R/helper.R), you can add the volume and
  # temperature conversion there
  input.values <- input(wtemp = data, H = H, A = A)
  
  proc.obs <- preprocess_obs(obs,input.values = input.values, H, A)
  w.obs <- weigh_obs(obs,input.values = input.values, H, A)
  obs_long <- w.obs[[1]]
  obs_weigh <- w.obs[[2]]
  

  fsed_stratified_epi = 0.01 *100
  fsed_stratified_hypo = 0.01 *100
  fsed_not_stratified  =  0.0002
  nep_stratified = 0.1
  nep_not_stratified = 0
  min_stratified = 0.1
  min_not_stratified = 0
  
  init.val = c(0.5, 0.5, 0.1, 0.01)
  target.iter = 30
  
  modelopt <- pureCMAES(par = init.val, fun = optim_do, lower = c(0.0, 0.0, -0.5, -0.1),
                    upper = c(1.0, 1.0, 0.5, 0.1), sigma = 0.5,
                    stopfitness = -Inf, 
                    stopeval = target.iter,
                    input.values = input.values,
                    fsed_not_stratified = fsed_not_stratified,
                    nep_not_stratified = nep_not_stratified, min_not_stratified = min_not_stratified, 
                    wind, proc.obs = obs_weigh,
                    verbose = verbose)

  o2<- calc_do(input.values = input.values,fsed_stratified_epi = modelopt$xmin[1],
               fsed_stratified_hypo = modelopt$xmin[2],
               fsed_not_stratified,
               nep_stratified = modelopt$xmin[3],
               nep_not_stratified,
               min_stratified = modelopt$xmin[4],
               min_not_stratified, wind)
  
  input.values$o2_epil <- o2[,"o2_epil"]
  input.values$o2_hypo <- o2[,"o2_hypo"]
  input.values$o2_total <- o2[,"o2_total"]

  input.values$year <- year(input.values$datetime)
  input.values$doy <- yday(input.values$datetime)

  
  test_data<-compare_predict_versus_observed(obs,input.values) 
  test_data$year <- year(test_data$day)
  test_data$doy <- yday(test_data$day)
  
  # fit = calc_rmse(test_data)
  fit = calc_fit(input.values , proc.obs = obs_weigh )
  
  pgm <- input.values %>%
    dplyr::select(datetime, o2_epil, o2_hypo, o2_total, vol_epil, vol_hypo, 'vol_total' = total_vol)
  input.pgm <- input.values %>%
    dplyr::select(datetime, "depth_td" = td.depth, "area_td" = td_area, "area_surf" = surf_area,
                  "wtemp_total" = t.total,
                  'wtemp_epil' = t.epil, "wtemp_hypo" = t.hypo,
                  'vol_total' = total_vol, vol_epil, vol_hypo)
  
  write_delim(input.pgm, path = paste0(ii,'/',sub("\\).*", "", sub(".*\\(", "", ii)) ,'_',round(fit,1),'_alldata.txt'), delim = '\t')
  write_delim(pgm, path = paste0(ii,'/',sub("\\).*", "", sub(".*\\(", "", ii)) ,'_',round(fit,1),'_oxymodel.txt'), delim = '\t')
  write_delim(obs_long, path = paste0(ii,'/',sub("\\).*", "", sub(".*\\(", "", ii)) ,'_obs.txt'), delim = '\t')

  observed <- data.frame('time' = input.values$datetime[obs_weigh[1,]],
                            'year' = input.values$year[obs_weigh[1,]],
                            'doy' = input.values$doy[obs_weigh[1,]],
                            'epi' =obs_weigh[3,],
                            'hypo' = obs_weigh[4,],
                         'epi_sim' = input.values$o2_epil[obs_weigh[1,]]/input.values$vol_epil[obs_weigh[1,]]/1000,
                         'hypo_sim' = input.values$o2_hypo[obs_weigh[1,]]/input.values$vol_hypo[obs_weigh[1,]]/1000)
  
  g1 <- ggplot(input.values) +
    geom_point(aes(doy, (o2_total/total_vol/1000), col = 'Total')) +
    geom_point(aes(doy, (o2_epil/vol_epil/1000), col = 'Epi')) +
    geom_point(aes(doy, (o2_hypo/vol_hypo/1000), col = 'Hypo')) +
    ylim(0,20)+
    facet_wrap(~year) +
    theme_bw() +
    geom_point(data = observed, aes(doy, epi, col = 'Obs_Epi'), size =2, alpha = 0.5) +
    geom_point(data = observed, aes(doy, hypo, col = 'Obs_Hypo'), size =2, alpha = 0.5)
  ggsave(file = paste0(ii,'/oxymodel.png'), g1, dpi=300,width=316,height=190,units='mm')


  g3 <- ggplot(observed, aes(time, epi - epi_sim, col = 'epilimnion')) +
    geom_point(size = 2) +
    geom_line() +
    geom_point(aes(time, hypo - hypo_sim, col = 'hypolimnion'), size =2) +
    geom_line(aes(time, hypo - hypo_sim, col = 'hypolimnion'))+
    ylab('Residuals (obs-mod)')+
    theme_bw()
  ggsave(file = paste0(ii,'/predicted_ag_observed.png'), g3, dpi=300, width=316,height=190,units='mm')
  
  anno = 2000
  
  g4 <- ggplot(subset(input.values, year == anno)) +
    geom_point(aes(doy, (o2_total/total_vol/1000), col = 'Total')) +
    geom_point(aes(doy, (o2_epil/vol_epil/1000), col = 'Epi')) +
    geom_point(aes(doy, (o2_hypo/vol_hypo/1000), col = 'Hypo')) +
    ylim(0,20)+
    #facet_wrap(~year) +
    theme_bw() +
    geom_point(data = subset(observed, year == anno), aes(doy, epi, col = 'Obs_Epi'), size =2, alpha = 0.5) +
    geom_point(data = subset(observed, year == anno), aes(doy, hypo, col = 'Obs_Hypo'), size =2, alpha = 0.5);g4
  ggsave(file = paste0(ii,'/oxymodel_2000.png'), g4, dpi=300,width=316,height=190,units='mm')
  
  
  eval.info <- data.frame('id' = sub("\\).*", "", sub(".*\\(", "", ii)) ,
                          'time' = Sys.time(),
             'A' = max(A),
             'z' = max(H),
             'obs' = ncol(proc.obs),
             'RMSE' = fit)
    write.table(eval.info, file ='eval.csv', append=TRUE, quote = FALSE, col.names = FALSE,
                row.names = FALSE)
  print('Nothing got broken, good job!')
}

library(viridis)
eval.df <- read.table('eval.csv', header = TRUE)

g<- ggplot(eval.df, aes(Asurf, RMSE, col = MaxZ, label = ID)) + 
  geom_point(aes(size = nobs)) +   
  scale_color_viridis(option="viridis") +
  xlab('Surface Area') +
  ylab('RMSE in mg DO/L') +
  geom_text(check_overlap = TRUE,hjust = 0.05, nudge_x = 0.05) + 
  theme_bw();g
ggsave(file = paste0('lake_results.png'), g, dpi=300, width=316,height=190,units='mm')


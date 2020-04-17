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
  
  ggplot(input.values)+
    geom_line(aes(doy, td.depth, col='td.depth')) +
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
  

  fsed_stratified_epi = 1500 # mg/m2/d
  fsed_stratified_hypo = 4500 # mg/m2/d
  fsed_not_stratified  =  0.0002 # mg/m2/d
  nep_stratified = 0.1 # mg/m3/d
  nep_not_stratified = 0 # mg/m3/d
  min_stratified = -5 # mg/m3/d
  min_not_stratified = 0 # mg/m3/d
  khalf <- 3000 # mg/m3

  o2<- calc_do(input.values = input.values,
               fsed_stratified_epi = fsed_stratified_epi,
               fsed_stratified_hypo = fsed_stratified_hypo,
               fsed_not_stratified = fsed_not_stratified,
               nep_stratified = nep_stratified,
               nep_not_stratified = nep_not_stratified,
               min_stratified = min_stratified,
               min_not_stratified =min_not_stratified,
               wind = wind,
               khalf = khalf,
               startdate = NULL,
               enddate = NULL)
  
  # fsed_stratified_epi fsed_stratified_hypo nep_stratified min_stratified, khalf
  init.val = c(5, 5, 5, 5, 5)
  target.iter = 30
  lb <<- c(100, 100, -5, -5, 1000)
  ub <<- c(6000, 6000, 5, 5, 6000)

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
                        wind = wind, proc.obs = obs_weigh, # caL_data
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
  for (ikx in unique(year(obs_long$ActivityStartDate ))){
    kind <- which(ikx == input.values$year)
    startdate <- min(kind)
    enddate <- max(kind)
    
    modelopt <- pureCMAES(par = init.val, fun = optim_do, lower = rep(0,5),
                          upper = rep(10,5), sigma = 0.5,
                          stopfitness = -Inf,
                          stopeval = target.iter,
                          input.values = input.values,
                          fsed_not_stratified = fsed_not_stratified,
                          nep_not_stratified = nep_not_stratified, min_not_stratified = min_not_stratified,
                          wind = wind, proc.obs = obs_weigh, # caL_data
                          verbose = verbose, startdate = startdate, enddate = enddate)
    
    modelopt$xmin <- lb+(ub - lb)/(10)*(modelopt$xmin)
    
    calibrated_param <- c(ikx, modelopt$xmin, modelopt$fmin)
    all_params <- rbind(all_params, calibrated_param)
    
    
    o2<- calc_do(input.values = input.values,fsed_stratified_epi = modelopt$xmin[1],
                 fsed_stratified_hypo = modelopt$xmin[2],
                 fsed_not_stratified,
                 nep_stratified = modelopt$xmin[3],
                 nep_not_stratified,
                 min_stratified = modelopt$xmin[4],
                 min_not_stratified, wind,
                 startdate = startdate, enddate = enddate)
    
    input.values$o2_epil[startdate:enddate] <- o2[startdate:enddate,"o2_epil"]
    input.values$o2_hypo[startdate:enddate] <- o2[startdate:enddate,"o2_hypo"]
    input.values$o2_total[startdate:enddate] <- o2[startdate:enddate,"o2_total"]
    input.values$sat_o2_epil[startdate:enddate] <- o2[startdate:enddate,"sat_o2_epil"]
    input.values$sat_o2_hypo[startdate:enddate] <- o2[startdate:enddate,"sat_o2_hypo"]
    input.values$sat_o2_total[startdate:enddate] <- o2[startdate:enddate,"sat_o2_total"]
    input.values$Fsed_total[startdate:enddate] <- o2[startdate:enddate,"Fsed_total"]
    input.values$NEP_total[startdate:enddate] <- o2[startdate:enddate,"NEP_total"]
    input.values$Fatm_total[startdate:enddate] <- o2[startdate:enddate,"Fatm_total"]
    input.values$Mineral_total[startdate:enddate] <- o2[startdate:enddate,"Mineral_total"]
    input.values$Fsed_epi[startdate:enddate] <- o2[startdate:enddate,"Fsed_epi"]
    input.values$NEP_epi[startdate:enddate] <- o2[startdate:enddate,"NEP_epi"]
    input.values$Fatm_epi[startdate:enddate] <- o2[startdate:enddate,"Fatm_epi"]
    input.values$Entrain_epi[startdate:enddate] <- o2[startdate:enddate,"Entrain_epi"]
    input.values$Fsed_hypo[startdate:enddate] <- o2[startdate:enddate,"Fsed_hypo"]
    input.values$Mineral_hypo[startdate:enddate] <- o2[startdate:enddate,"Mineral_hypo"]
    input.values$Entrain_hypo[startdate:enddate] <- o2[startdate:enddate,"Entrain_hypo"]
    
    all_o2 <- rbind(all_o2, o2[startdate:enddate,])
  }
  
  m.all_params <- as.data.frame(all_params)
  colnames(m.all_params) <- c('year', 'Fsed_epi', 'Fsed_hypo', 'NEP', 'Mineral',
                            'khalf', 'RMSE')
  i1 <- ggplot(m.all_params, aes(year, RMSE)) +
    geom_line() + geom_hline(yintercept=1.97)
  i2 <- ggplot(m.all_params, aes(year, Fsed_epi)) +
    geom_line() + geom_hline(yintercept=1921.07)
  i3 <- ggplot(m.all_params, aes(year, Fsed_hypo)) +
    geom_line() + geom_hline(yintercept=6000)
  i4 <- ggplot(m.all_params, aes(year, NEP)) +
    geom_line() + geom_hline(yintercept=3.92)
  i5 <- ggplot(m.all_params, aes(year, Mineral)) +
    geom_line() + geom_hline(yintercept=-2.52)
  i6 <- ggplot(m.all_params, aes(year, khalf)) +
    geom_line() + geom_hline(yintercept=3430.83)
  
  p <- i1 + i2 + i3 + i4 + i5 + i6 & theme_minimal()
  p 
  ggsave(file = paste0(ii,'/comparison_indyearcal.png'), p, dpi=300,width=316,height=190,units='mm')
  
  
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
    geom_point(aes(doy, Entrain_epi, col = 'Entrain_epi')) +
    geom_point(aes(doy, Fsed_hypo, col = 'Fsed_hypo')) +
    geom_point(aes(doy, Mineral_hypo, col = 'Mineral_hypo')) +
    geom_point(aes(doy, Entrain_hypo, col = 'Entrain_hypo')) +
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

g1<- ggplot(eval.df, aes(Asurf, RMSE_cal, col = MaxZ, label = ID)) + 
  geom_point(aes(size = nobs)) +   
  scale_color_viridis(option="viridis") +
  xlab('Surface Area') +
  ylab('RMSE in mg DO/L') + 
  ggtitle("calibration") +
  ylim(0,6)+
  geom_text(check_overlap = TRUE,hjust = 0.05, nudge_x = 0.05) + 
  theme_bw();g1
g2<- ggplot(eval.df, aes(Asurf, RMSE_val, col = MaxZ, label = ID)) + 
  geom_point(aes(size = nobs)) +   
  scale_color_viridis(option="viridis") +
  xlab('Surface Area') +
  ylab('RMSE in mg DO/L') +
  ggtitle("validation") +
  ylim(0,6)+
  geom_text(check_overlap = TRUE,hjust = 0.05, nudge_x = 0.05) + 
  theme_bw();g2
g <- g1 + g2 + plot_annotation(tag_levels = 'A',
  caption = paste0('ODEM ',Sys.time()) 
)
ggsave(file = paste0('lake_results.png'), g, dpi=300, width=316,height=190,units='mm')


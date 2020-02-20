rm(list = ls())

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

## packages
# get glmtools by first installing devtools, e.g. install.packages('devtools), then
# devtools::install_github('USGS-R/glmtools', ref='ggplot_overhaul')
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

#devtools::install_github('LynetteGao/Limno_DataScience')
library(simpleAnoxia)

## source all functions
# source('R/helper.R')

# # load shapefiles
# library(sf)
# lakes = st_read(file.path('inst/extdata/study_lakes.shp'))
# ggplot(lakes) + geom_sf()
# poly.distance <- function(polygon){
#   polygon <- st_sfc(polygon, crs = 4269)
#   return(max(st_distance(st_cast(polygon, "POINT"))))
# }
# str(lakes)
# poly.distance(st_geometry(lakes[4,]))
# df.poly <- matrix(NA,ncol=2,nrow=nrow(lakes))
# for (jj in 2001:2500){
#   df.poly[jj,1] <- poly.distance(st_geometry(lakes[jj,]))
#   df.poly[jj,2] <- 10^(0.336 * log10(df.poly[jj,1]) - 0.245) # https://www.nrcresearchpress.com/doi/pdfplus/10.1139/f90-108
#   print(paste0(jj,' of ', nrow(df.poly),' and we use ',object.size(df.poly)))
# }


# sapply(st_geometry(lakes), poly.distance)
# purrr::map(st_geometry(lakes), poly.distance)

## load example data
lks <- list.dirs(path = 'inst/extdata/', full.names = TRUE, recursive = F)

for (ii in lks[-9]){
  print(paste0('Running ',ii))
  data <- read.csv(paste0(ii,'/', list.files(ii, pattern = 'temperatures.csv', include.dirs = T)))
  
  
  if (length( list.files(ii, pattern = 'wq_data', include.dirs = T)) > 0){
  wq_data<- paste0(ii,'/', list.files(ii, pattern = 'wq_data', include.dirs = T))
  obs <- NULL
  # obs <- data.frame() 
  for(jj in wq_data){
      raw_obs <- read.csv(jj)
      # filter_obs<-filter(raw_obs,CharacteristicName== "Dissolved oxygen (DO)")
      # more_obs <- filter_obs %>%
        # dplyr::select(c('ActivityStartDate', 'ActivityStartTime.Time', 'ActivityDepthHeightMeasure.MeasureValue', 'ResultMeasureValue'))
      
      wq <- raw_obs %>%
        dplyr::filter(CharacteristicName== "Dissolved oxygen (DO)") %>%
        dplyr::select(c('ActivityStartDate', 'ActivityStartTime.Time', 'ActivityDepthHeightMeasure.MeasureValue', 'ResultMeasureValue'))
        
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
  
  
  eg_nml <- read_nml(paste0(ii,'/', list.files(ii, pattern = 'nml', include.dirs = T)))
  H <- abs(eg_nml$morphometry$H - max(eg_nml$morphometry$H))
  A <- eg_nml$morphometry$A

  # here's the promised input function (see R/helper.R), you can add the volume and
  # temperature conversion there
  input.values <- input(wtemp = data, H = H, A = A)

  fsed_stratified = 0.01 *100
  fsed_not_stratified  =  0.0002
  nep_stratified = 0.1
  nep_not_stratified = 0
  
  init.val = c(0.5, 0.1)
  target.iter = 10
  modelopt <- neldermeadb(fn = optim_do, init.val, lower = c(0., -0.5),
                          upper = c(1.0, 0.5), adapt = TRUE, tol = 1e-2,
                          maxfeval = target.iter, input.values = input.values,
                          fsed_not_stratified = fsed_not_stratified, 
                          nep_not_stratified = nep_not_stratified, verbose = verbose)

  o2<- calc_do(input.values = input.values,fsed_stratified = modelopt$xmin[1],
               fsed_not_stratified,
               nep_stratified = modelopt$xmin[2],
               nep_not_stratified)
  
  input.values$o2_epil <- o2[,"o2_epil"]
  input.values$o2_hypo <- o2[,"o2_hypo"]
  input.values$o2_total <- o2[,"o2_total"]

  input.values$year <- year(input.values$datetime)
  input.values$doy <- yday(input.values$datetime)

  
  test_data<-compare_predict_versus_observed(obs,input.values) 
  
  fit = calc_rmse(test_data)
  
  pgm <- input.values %>%
    dplyr::select(datetime, o2_epil, o2_hypo, o2_total, vol_epil, vol_hypo, 'vol_total' = total_vol)
  
  write_delim(input.values, path = paste0(ii,'/',sub("\\).*", "", sub(".*\\(", "", ii)) ,'_',round(fit,1),'_alldata.txt'), delim = '\t')
  write_delim(pgm, path = paste0(ii,'/',sub("\\).*", "", sub(".*\\(", "", ii)) ,'_',round(fit,1),'_oxymodel.txt'), delim = '\t')
  write_delim(obs, path = paste0(ii,'/',sub("\\).*", "", sub(".*\\(", "", ii)) ,'_obs.txt'), delim = '\t')

  g1 <- ggplot(input.values) +
    geom_point(aes(doy, (o2_total/total_vol/1000), col = 'Total')) +
    geom_point(aes(doy, (o2_epil/vol_epil/1000), col = 'Epi')) +
    geom_point(aes(doy, (o2_hypo/vol_hypo/1000), col = 'Hypo')) +
    ylim(0,25)+
    facet_wrap(~year) +
    theme_bw()
  ggsave(file = paste0(ii,'/oxymodel.png'), g1, dpi=300, width=216,height=150,units='mm')
  
  g2 <- ggplot(input.values, aes(doy, td.depth)) +
    geom_line() +
    facet_wrap(~year) +
    theme_bw()
  ggsave(file = paste0(ii,'/predicted_td.png'), g2, dpi=300, width=216,height=150,units='mm')
  

  g3 <- ggplot(test_data, aes(observed_hypo_do, predict_hypo_do, col = 'hypolimnion')) +
    geom_point() +
    geom_point(aes(observed_epil_do, predict_epil_do, col = 'epilimnion')) +
    ylim(0,25)+
    xlim(0,25)+
    theme_bw()
  ggsave(file = paste0(ii,'/predicted_ag_observed.png'), g3, dpi=300, width=216,height=216,units='mm')
  
    
  print('Nothing got broken, good job!')
}

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
library(LakeMetabolizer)
library(rgdal)


devtools::install_github('LynetteGao/Limno_DataScience')
library(simpleAnoxia)

## source all functions
#source('R/helper.R')

## load example data
lks <- list.dirs(path = 'inst/extdata/', full.names = TRUE, recursive = F)
shp <- readOGR(dsn = file.path('inst/extdata/study_lakes.shp'), stringsAsFactors = F)
map <- ggplot() + geom_polygon(data = shp, aes(x = long, y = lat, group = group), colour = "black", fill = NA)
map <- map + theme_void()
ggsave(file=paste0('lakes_location.png'), map, dpi = 300,width = 450,height = 400, units = 'mm')


for (ii in lks[5]){
  data <- read.csv(paste0(ii,'/', list.files(ii, pattern = 'temperatures.csv', include.dirs = T)))
  raw_obs <- read.csv(paste0(ii,'/', list.files(ii, pattern = 'observed', include.dirs = T)))
  obs <- raw_obs %>%
    dplyr::select(c('ActivityStartDate', 'ActivityStartTime.Time', 'ActivityDepthHeightMeasure.MeasureValue', 'ResultMeasureValue'))
  
  eg_nml <- read_nml(paste0(ii,'/', list.files(ii, pattern = 'nml', include.dirs = T)))
  H <- abs(eg_nml$morphometry$H - max(eg_nml$morphometry$H))
  A <- eg_nml$morphometry$A
  
  # here's the promised input function (see R/helper.R), you can add the volume and 
  # temperature conversion there
  input.values <- input(wtemp = data, H = H, A = A)
  
  fsed_stratified = 0.01 *100
  fsed_not_stratified  =  0.0002
  nep_stratified = 0.1#max(A)*2
  nep_not_stratified = 0
  
  o2<- calc_do(input.values = input.values,fsed_stratified,fsed_not_stratified,nep_stratified,nep_not_stratified)
  input.values$o2_epil <- o2[,"o2_epil"]
  input.values$o2_hypo <- o2[,"o2_hypo"]
  input.values$o2_total <- o2[,"o2_total"]
}

input.values$year <- year(input.values$datetime)
input.values$doy <- yday(input.values$datetime)

ggplot(subset(input.values, year == '1979')) +
  geom_line(aes(doy, td.depth, col = 'Total')) 

ggplot(subset(input.values, year == '1995')) +
 # geom_line(aes(doy, (o2_total/total_vol/1000), col = 'Total')) +
  geom_line(aes(doy, (o2_epil/vol_epil/1000), col = 'Epi')) +
  geom_line(aes(doy, (o2_hypo/vol_hypo/1000), col = 'Hypo')) +
  # ylim(-10, 1e30)+
  facet_wrap(~year) +
  theme_bw()

ggplot(input.values) +
  geom_point(aes(doy, (o2_total/total_vol/1000), col = 'Total')) +
  geom_point(aes(doy, (o2_epil/vol_epil/1000), col = 'Epi')) +
  geom_point(aes(doy, (o2_hypo/vol_hypo/1000), col = 'Hypo')) +
  ylim(0,25)+
  facet_wrap(~year) +
  theme_bw()


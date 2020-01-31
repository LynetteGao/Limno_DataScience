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

devtools::install_github('LynetteGao/Limno_DataScience')
library(simpleAnoxia)

## source all functions
#source('R/helper.R')

## load example data
lks <- list.dirs(path = 'inst/extdata/', full.names = TRUE, recursive = F)



for (ii in lks[1]){
  data <- read.csv(paste0(ii,'/', list.files(ii, pattern = 'csv', include.dirs = T)))
  
  eg_nml <- read_nml(paste0(ii,'/', list.files(ii, pattern = 'nml', include.dirs = T)))
  H <- abs(eg_nml$morphometry$H - max(eg_nml$morphometry$H))
  A <- eg_nml$morphometry$A
  
  # here's the promised input function (see R/helper.R), you can add the volume and 
  # temperature conversion there
  input.values <- input(wtemp = data, H = H, A = A)
  
  fsed_stratifed = 0.01
  fsed_not_stratified  =  0.0002
  nep_stratified = max(A)*2
  nep_not_stratified = 0
  
  o2<- calc_do(input.values = input.values,fsed_stratifed,fsed_not_stratified,nep_stratified,nep_not_stratified)
  input.values$o2_epil <- o2[,"o2_epil"]
  input.values$o2_hypo <- o2[,"o2_hypo"]
  input.values$o2_total <- o2[,"o2_total"]
}
space_time <- extract_time_space(data)

input.values$timedate <- space_time$datetime
input.values$year <- year(space_time$datetime)
input.values$doy <- yday(space_time$datetime)


ggplot(input.values) +
  geom_line(aes(doy, o2_epil, col = 'Epi')) +
  geom_line(aes(doy, o2_hypo, col = 'Hypo')) +
  ylim(...)+
  facet_wrap(~year) +
  theme_bw()

ggplot(input.values) +
  geom_line(aes(doy, td.depth)) +
  facet_wrap(~year) +
  theme_bw()


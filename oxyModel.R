rm(list = ls())

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

## packages
# get glmtools by first installing devtools, e.g. install.packages('devtools), then
# devtools::install_github('USGS-R/glmtools', ref='ggplot_overhaul')
library(glmtools) 
library(dplyr)

## source all functions
source('R/helper.R')

## load example data
lks <- list.dirs(path = 'inst/extdata/', full.names = TRUE, recursive = F)

for (ii in lks){
  data <- read.csv(paste0(lks,'/', list.files(lks, pattern = 'csv', include.dirs = T)))
  
  eg_nml <- read_nml(paste0(lks,'/', list.files(lks, pattern = 'nml', include.dirs = T)))
  H <- abs(eg_nml$morphometry$H - max(eg_nml$morphometry$H))
  A <- eg_nml$morphometry$A
  
  # here's the promised input function (see R/helper.R), you can add the volume and 
  # temperature conversion there
  input.values <- input(wtemp = data, H = H, A = A)
  
  
  
  
  
  
  
  
}


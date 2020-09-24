# ii <- "analysis/data/143249470_Mendota"
args <- commandArgs(trailingOnly=TRUE)
if (length(args) > 0){
  ii   <- args[1]
}else{
  if (!exists("ii")){
    stop("Must set the variable 'ii' to a data folder path.")
  }
}

# helper functions used:
#   input
#     extract_time_space
#     calc_td_depth
#     calc_metalimdepth
#     calc_epil_hypo_temp
#     calc_vol_total
#     calc_epil_hypo_vol

source("analysis/scripts/99_packages.R")

# ii <- "analysis/data//141288680_(Fish)"
print(paste0('Running ',ii))

data  <- read.csv(paste0(ii,'/', 
                         list.files(ii, pattern = 'pball.*_temperatures.csv$', 
                                    include.dirs = T))) # GLM sim water temp
meteo <- read.csv(paste0(ii,'/', 
                         list.files(ii, pattern = 'NLDAS',
                                    include.dirs = T))) # meteo NLDAS

chidx   <- match(as.POSIXct(data$date),as.POSIXct(meteo$time))
wind    <- meteo$WindSpeed[chidx]
airtemp <- meteo$AirTemp[chidx]

if (length(list.files(ii, pattern = 'wq_data', include.dirs = T)) > 0){ # wq LTER
  wq_data <- paste0(ii,'/', 
                    list.files(ii, pattern = 'wq_data', include.dirs = T))
  obs     <- NULL
  
  for(jj in wq_data){
    # jj <- wq_data[1]
    raw_obs <- read.csv(jj)
    if ('ActivityDepthHeighMeasure.MeasureValue' %in% colnames(raw_obs)){
      wq <- raw_obs %>%
        dplyr::filter(CharacteristicName== "Dissolved oxygen (DO)") %>%
        dplyr::select(
          c('ActivityStartDate', 'ActivityDepthHeighMeasure.MeasureValue',
            'ResultMeasureValue'))
      wq <- rename(wq, 
                   'ActivityDepthHeightMeasure.MeasureValue' = 
                     'ActivityDepthHeighMeasure.MeasureValue')
    } else {
      wq <- raw_obs %>%
        dplyr::filter(CharacteristicName== "Dissolved oxygen (DO)") %>%
        dplyr::select(
          c('ActivityStartDate', 'ActivityDepthHeightMeasure.MeasureValue',
            'ResultMeasureValue'))
    }
    
    if (length(wq_data) == 1 | is.null(obs)){
      obs <- wq
    } else {
      obs <- rbind(obs, wq)
    }
  }
  obs$ActivityStartDate <- as.POSIXct(obs$ActivityStartDate)
}

if (is.factor(obs$ActivityDepthHeightMeasure.MeasureValue)){
  obs$ActivityDepthHeightMeasure.MeasureValue <- 
    as.numeric(as.character(obs$ActivityDepthHeightMeasure.MeasureValue))
}
if (is.factor(obs$ResultMeasureValue)){
  obs$ResultMeasureValue <- 
    as.numeric(as.character(obs$ResultMeasureValue))
}

# outlier detection
outlier_values              <- boxplot.stats(obs$ResultMeasureValue)$out 
uvx                         <- match(outlier_values, obs$ResultMeasureValue)
obs$ResultMeasureValue[uvx] <- NA

eg_nml <- read_nml(paste0(ii,'/',  # NML GLM 
                          list.files(ii, pattern = 'nml', include.dirs = T)))

H <- abs(eg_nml$morphometry$H - max(eg_nml$morphometry$H)) # DEPTH
A <- eg_nml$morphometry$A # AREA

input.values         <- input(wtemp = data, H = H, A = A)
input.values$year    <- year(input.values$datetime)
input.values$doy     <- yday(input.values$datetime)
input.values$wind    <- wind
input.values$airtemp <- airtemp

input.values <-
  dplyr::rename(input.values,
    datetime = datetime,
    thermocline_depth = td.depth,
    temperature_epi = t.epil,
    temperature_hypo = t.hypo,
    temperature_total = t.total,
    volume_total = vol_total,
    volume_epi = vol_epil,
    volume_hypo = vol_hypo,
    area_thermocline = td_area,
    area_surface = surf_area,
    upper_meta = upper.metalim,
    lower_meta = lower.metalim,
    year = year,
    day_of_year = doy,
    max.d = max.d,
    wind = wind,
    airtemp = airtemp)

write.table(input.values, paste0(ii, '/input.txt'),
            append = FALSE, sep = ",", dec = ".",
            row.names = FALSE, col.names = TRUE) # input data --> B-ODEM

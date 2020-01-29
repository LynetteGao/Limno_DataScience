## Helper functions

#' Calculate water density from temperature
#'
#' Calculate water density from water temperature using the formula from (Millero & Poisson, 1981).
#'
#' @param wtemp vector or matrix; Water temperatures
#' @return vector or matrix; Water densities in kg/m3
#' @export
calc_dens <- function(wtemp){
  dens = 999.842594 + (6.793952 * 10^-2 * wtemp) - (9.095290 * 10^-3 *wtemp^2) + (1.001685 * 10^-4 * wtemp^3) - (1.120083 * 10^-6* wtemp^4) + (6.536336 * 10^-9 * wtemp^5)
  return(dens)
}

#' Extract time and space information
#'
#' Extracts time (from date column) and space (aka depth) information
#'
#' @param wtemp matrix; Water temperatures (rows correspond to time, cols to depth)
#' @return list of datetimes and depths
#' @export
extract_time_space <- function(wtemp){
  time <- as.Date(as.character(wtemp$date))
  depth <- sub("^[^_]*_", "", colnames(wtemp[2:ncol(wtemp)]))
  return(list('datetime' = time, 'depth' = depth))
}

#' Calculate thermocline depth
#'
#' Calculate planar thermocline depth by checking the highest density gradient over depth.
#'
#' @param wtemp matrix; Water temperatures (rows correspond to time, cols to depth)
#' @return vector of thermocline depths in m
#' @export
calc_td_depth <- function(wtemp){
  
  grd.info <- extract_time_space(wtemp)
  temp <- as.matrix(wtemp[,-c(1)])
  dens <- calc_dens(temp)
  td.depth <- rep(NA, length(grd.info$datetime))
  
  for (ii in 1:length(td.depth)){
    # simplify data
    idx = !is.na(temp[ii,])
    temp_data = temp[ii,idx]
    dens_data = dens[ii,idx]
    depth_data = as.double(grd.info$depth[idx])
    
    td.diff <- rep(NA, length(depth_data))
    # forward differencing d rho / d z
    td.diff <- (lead(dens_data) - dens_data)/
      (lead(depth_data) - depth_data) 
    # backward differencing d rho / d z
    td.diff[length(td.diff)] <- (dens_data[length(td.diff)] - dens_data[length(td.diff)-1])/
      (depth_data[length(td.diff)] - depth_data[length(td.diff)-1]) 
    
    # find density difference between epilimnion and hypolimnion
    dens.diff = rev(dens_data)[1] - dens_data[1]
    
    if (min(temp_data) > 4 && abs(dens.diff) > 0.1){
      d.idx = which(td.diff == max(td.diff))[1]
      td.depth[ii] = depth_data[d.idx]
    } else {td.depth[ii] = NA}
  }
  return(td.depth)
}

#'
#' Calculate mean epilimnion and hypolimnion surface temperature,using the thermocline depth
#'
#' @param wtemp matrix; Water temperatures (rows correspond to time, cols to depth)
#' @return vector of thermocline depths in m
calc_epil_hypo_temp<-function(wtemp,td.depth){
  grd.info <- extract_time_space(wtemp)
  temp <- as.matrix(wtemp[,-c(1)])
  depth_data = as.double(grd.info$depth)
  
  epil_temp <-  rep(NA, length(td.depth))
  hypo_temp <- rep(NA, length(td.depth))
  total_temp<- rep(NA, length(td.depth))
  
  for (ii in 1:length(td.depth)){
    idx = !is.na(temp[ii,])
    temp_data = as.numeric(temp[ii,idx])
    if(is.na(td.depth[ii])){
      total_temp[ii]<-sum(temp_data)/length(temp_data)
    }else{
      td_idx <- min(which(td.depth[ii]>depth_data))
      epil_temp[ii] <- mean(temp_data[1:td_idx])
      hypo_temp[ii] <- mean(temp_data[(td_idx+1):length(temp_data)])
    }
  }
  return(list(epil_temp,hypo_temp,total_temp))
  
}



#' Create temperature and volume input values for oxygen model
#'
#' Calculate mean epilimnion and hypolimnion surface temperature, as well as volumes.
#'
#' @param wtemp matrix; Water temperatures (rows correspond to time, cols to depth)
#' @return vector of thermocline depths in m
#' @export
input <- function(wtemp, H, A){
  td.depth <- calc_td_depth(wtemp)
    
  return(data.frame(td.depth))
}









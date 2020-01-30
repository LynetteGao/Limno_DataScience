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
#' @param td.depth matrix; thermocline depth
#' @param H the depth info of the lake
#' @return list of temperatures vector
calc_epil_hypo_temp<-function(wtemp,td.depth,H){
  grd.info <- extract_time_space(wtemp)
  temp <- as.matrix(wtemp[,-c(1)])
  depth_data = as.double(grd.info$depth)
  
  epil_temp <-  rep(NA, length(td.depth))
  hypo_temp <- rep(NA, length(td.depth))
  total_temp<- rep(NA, length(td.depth))
  
  for (ii in 1:length(td.depth)){
    idx = !is.na(temp[ii,])
    temp_data = as.numeric(temp[ii,idx])
    total_temp[ii]<-sum(temp_data)/length(temp_data)
    if(!is.na(td.depth[ii])){
      td_idx <- max(which(td.depth[ii]>=depth_data))
      epil_temp[ii] <- mean(temp_data[1:td_idx])
      hypo_temp[ii] <- mean(temp_data[(td_idx+1):length(temp_data)])
    }
  }
  return(list('t_epil' = epil_temp,'t_hypo' = hypo_temp,'t_total' = total_temp))
  
}

#'
#' Calculate water total volume,using the thermocline depth
#'
#' @param wtemp matrix; Water temperatures (rows correspond to time, cols to depth)
#' @param td.depth matrix; thermocline depth
#' @return list of temperatures
#' @importFrom pracma trapz 
calc_total_vol<-function(H,A){
  if (length(H)==1){
    vol_total <- 1/3.0 * A * H
  }else{
    H.diff <- rep(NA, length(A))
    # integration of A over dh
    # H.diff <-  lag(H)-H
    # vol_total<- sum(A * (lag(H)-H),na.rm=TRUE)
    vol_total <- trapz(rev(H),rev(A))
  }
  return (vol_total)
}

#'
#' Calculate water epilimnion and hypolimnine volume,using the thermocline depth
#'
#' @param H vector; the depth info of the lake
#' @param A vector; the area info of the lake for each depth
#' @param td.depth matrix; thermocline depth
#' @param vol_total number;the total volume of the lake
#' @return matrix of epil and hypo volume
calc_epil_hypo_vol <- function(H,A,td.depth,vol_total){
  vol_data <- matrix(NA, nrow = length(td.depth), ncol = 2)
  colnames(vol_data) <- c("vol_epil","vol_hypo")
  for (ii in 1:length(td.depth)){
    if(!is.na(td.depth[ii])){
      h_idx <- min(which(td.depth[ii]>H))
      vol_data[ii,2] <- trapz(rev(H[h_idx:length(H)]),rev(A[h_idx:length(H)]))
      print(rev(H[h_idx:length(H)]))
      print(rev(A[h_idx:length(H)]))
      vol_data[ii,1] <- vol_total - vol_data[ii,2]
    }
  }
  return (vol_data)
}
  

#' Create temperature and volume input values for oxygen model
#'
#' Calculate mean epilimnion and hypolimnion surface temperature, as well as volumes.
#'
#' @param wtemp matrix; Water temperatures (rows correspond to time, cols to depth)
#' @return vector of thermocline depths in m
#' @export
input <- function(wtemp, H, A){
  grd.info <- extract_time_space(wtemp)
  td.depth <- calc_td_depth(wtemp)
  temp_out<-calc_epil_hypo_temp(wtemp,td.depth,H)
  total_vol<- calc_total_vol(H,A)
  vol<-calc_epil_hypo_vol(H,A,td.depth,total_vol)
  return(data.frame(cbind(datetime = grd.info$datetime,td.depth,t.epil = temp_out$t_epil,t.hypo=temp_out$t_hypo,
                          t.total = temp_out$t_total,total_vol,vol)))
}






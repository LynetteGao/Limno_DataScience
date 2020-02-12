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
    total_temp[ii]<-max(sum(temp_data)/length(temp_data),4)
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
      h_idx <- min(which(td.depth[ii]>=H))
      approx_td.area<-approx(H, A, c(0, td.depth[ii]))$y[2]
      H_with_td<-c(td.depth[ii],H[h_idx:length(H)])
      A_with_td<-c(approx_td.area,A[h_idx:length(A)])
      vol_data[ii,1] <- trapz(rev(H_with_td),rev(A_with_td))##epil
      vol_data[ii,2] <- vol_total - vol_data[ii,1]##hypo
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
  return(data.frame(datetime = as.POSIXct(grd.info$datetime),td.depth,t.epil = temp_out$t_epil,t.hypo=temp_out$t_hypo,
                          t.total = temp_out$t_total,total_vol,vol))
}



#' Calculate the mass of dissolved oxgen in epil and hypolimnion(during stratified period)
#' and the total dissovled oxygen in the lake
#' @param wtemp matrix; Water temperatures (rows correspond to time, cols to depth)
#' @return list of datetimes and depths
#' @export
#' 
calc_do<-function(input.values,fsed_stratified,fsed_not_stratified,nep_stratified,nep_not_stratified){
  ##initialize the o2(hypo),o2(epil),o2(total)
  o2_data <- matrix(NA, nrow = length(input.values$td.depth), ncol = 3)
  colnames(o2_data) <- c("o2_epil","o2_hypo","o2_total")
  o2_data <- as.data.frame(o2_data)
  
  init_o2sat <- o2.at.sat.base(temp=input.values$t.total[1],altitude = 300)*1000 # returns mg O2/L 
  o2_data$o2_total[1] <- init_o2sat*input.values$total_vol[1] # returns mg O2 (m3 = 1000 L)
  
  K600 <- k.cole.base(2) # returns m/day, assuming low wind conditions
  theta<-1.08
  for(day in 2:length(input.values$td.depth)){
    ## not stratified period, only consider the o2(total)
    if(is.na(input.values$td.depth[day])){
      theta_total <- theta^(input.values$t.total[day]-20)
      
      NEP <- nep_not_stratified * input.values$total_vol[day] * theta_total
      
      kO2 <- k600.2.kGAS.base(k600=K600,temperature=input.values$t.total[day],gas='O2') # velocity value m/d?
      o2sat<-o2.at.sat.base(temp=input.values$t.total[day],altitude = 300)*1000 # mg O2/L 
      Fatm <- kO2*(o2sat*input.values$total_vol[day] - o2_data[day-1,"o2_total"])/max(H)  # mg/m * m/d = mg/d
  
      Fsed <- fsed_not_stratified * o2_data[day-1,"o2_total"]/max(H) * theta_total  # mg/m * m/d = mg/d
      
      o2_data[day,"o2_total"] <- o2_data[day-1,"o2_total"] - Fsed + NEP + Fatm # units make sense bc every term is actually multiplied
      # print((o2_data[day,"o2_total"]/input.values$total_vol[day])/1000)
      # with delta t
    }
    # the day it turns to stratified, need to reassign the o2 to hypo and epil
    else if(is.na(input.values$td.depth[day-1])){
      o2_data[day,"o2_epil"] <- (o2_data[day-1,"o2_total"]*input.values$vol_epil[day])/input.values$total_vol[day]#/input.values$t.total[day]*input.values$vol_epil[day]
      o2_data[day,"o2_hypo"] <- (o2_data[day-1,"o2_total"]*input.values$vol_hypo[day])/input.values$total_vol[day]#/input.values$t.total[day]*input.values$vol_hypo[day]
      o2_data[day,"o2_total"]<- o2_data[day,"o2_epil"]+o2_data[day,"o2_hypo"]
    }else{
      theta_epil <- theta^(input.values$t.epil[day]-20)
      theta_hypo <-  theta^(input.values$t.hypo[day]-20)
      
      ##epil = epil.O2[i-1]+NEP[i]+Fatm[i]
      NEP_epil <- nep_stratified * input.values$vol_epil[day] * theta_epil # has to return mg/d, mg/m3/d * m3
      
      kO2_epil <- k600.2.kGAS.base(k600=K600,temperature=input.values$t.epil[day],gas='O2')
      o2sat_epil<-o2.at.sat.base(temp=input.values$t.epil[day],altitude = 300)*1000
      Fatm_epil <- kO2_epil*(o2sat_epil*input.values$vol_epil[day]-o2_data[day-1,"o2_epil"] )/input.values$td.depth[day] #possible wrong z
      
      o2_data[day,"o2_epil"] <- o2_data[day-1,"o2_epil"]+NEP_epil+Fatm_epil
      # print((o2_data[day,"o2_epil"]/input.values$vol_epil[day])/1000)
      
      ##hypo = hypo_o2[i-1]+Fhypo[i] - Fsed[i] + NEP[i]+Fatm[i]
      volumechange_hypo = input.values$vol_hypo[day]-input.values$vol_hypo[day-1]  #in m^3
      volumechange_hypo_proportion =  volumechange_hypo/input.values$vol_hypo[day-1] 
      Fhypo <- volumechange_hypo_proportion* o2_data[day-1,"o2_hypo"]
      
      Fsed <- fsed_stratified *o2_data[day-1,"o2_hypo"]/(max(H) - input.values$td.depth[day] ) * theta_hypo
      
      # NEP_hypo <- nep_stratified * input.values$vol_hypo[day]
      
      #kO2_hypo <- k600.2.kGAS.base(k600=K600,temperature=input.values$t.hypo[day],gas='O2')
      #o2sat_hypo<-o2.at.sat.base(temp=input.values$t.hypo[day],altitude = 300)*1000
      # Fatm_hypo <- kO2_hypo*(o2sat_hypo*input.values$vol_hypo[day]-o2_data[day-1,"o2_hypo"] )/(max(H)-input.values$td.depth[day])
      
      o2_data[day,"o2_hypo"] <- o2_data[day-1,"o2_hypo"] + Fhypo - Fsed #+ NEP_hypo +Fatm_hypo
      # print((o2_data[day,"o2_hypo"]/input.values$vol_hypo[day])/1000)
      ## total = hypo+epil
      o2_data[day,"o2_total"] <- o2_data[day,"o2_hypo"] + o2_data[day,"o2_epil"]
      # print((o2_data[day,"o2_total"]/input.values$total_vol[day])/1000)
    }
  }
  return (o2_data)
}











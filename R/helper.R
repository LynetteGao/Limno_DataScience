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
#' 
calc_td_depth <- function(wtemp){
  
  grd.info <- extract_time_space(wtemp)
  temp <- as.matrix(wtemp[,-c(1)])
  dens <- calc_dens(temp)
  
  cbuoy.depth <- rep(NA, length(grd.info$datetime))
  td.depth <- rep(NA, length(grd.info$datetime))
  
  condition<- apply(temp, 1, FUN=min,na.rm=TRUE) > 4
  
  for (ii in 1:length(cbuoy.depth)){
    idx = !is.na(temp[ii,])
    dens_data = dens[ii,idx]
    dens.diff = rev(dens_data)[1] - dens_data[1]
    
    if (condition[ii] && abs(dens.diff) > 0.05){
    cbuoy.depth[ii] <- center.buoyancy(temp[ii,], as.numeric(grd.info$depth))
    td.depth[ii] <- thermo.depth(temp[ii,], as.numeric(grd.info$depth))
    }
  }
  
  zdeps <- as.numeric(grd.info$depth)
  wlm.depth <- rep(NA, length(grd.info$datetime))
  
  for (ii in 1:length(wlm.depth)){
    idx = !is.na(temp[ii,])
    dens_data = dens[ii,idx]
    dens.diff = rev(dens_data)[1] - dens_data[1]
    
    if (condition[ii] && abs(dens.diff) > 0.05){
      
      Ch <- rep(NA, length = length(dens_data))
      for (jj in 1:(length(dens_data)-1)){
        Ah = 1/(zdeps[jj+1]) * sum(dens_data[1:jj])
        Bh = 1/(zdeps[length(zdeps)] -zdeps[jj]) * sum(dens_data[(jj + 1): length(dens_data)]) 
        
        diffAh = sum( (dens_data[jj:(jj+1)] - Ah)^2 )
        diffBh = sum( (dens_data[jj:(jj+1)] - Bh)^2 )
        
        Ch[jj] = diffAh + diffBh
      }
      clineDep = zdeps[which.min(na.omit(Ch))]
      wlm.depth[ii] <- clineDep
    }
  }
  
  # plot(cbuoy.depth)
  # points(td.depth, col='red')
  # 
  # 
  test <- data.frame('year' = year(grd.info$datetime), 'doy' = yday(grd.info$datetime),
                     'depth' = cbuoy.depth) #wlm.depth
  
  for (kk in unique(test$year)){
    idx <- which(kk == test$year)
    
    dx <- test[idx,3]
    dx[which(dx == (max(zdeps-1)))] = NA
    dx[which(dx == 0)] = NA
    
    NonNAindex <- which(!is.na(dx))
    if (length(na.omit(dx)) != 0){
      firstNonNA <- min(NonNAindex)
      lastNonNA <- max(NonNAindex)
      dx[firstNonNA:lastNonNA] =na.approx(dx)
    }
    test[idx,3] <-  dx
  }
  # 

  return(test$depth)#return(cbuoy.depth)#return(test$depth)
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
  
  total<- rep(NA, length(td.depth))
  hypo<- rep(NA, length(td.depth))
  total<- rep(NA, length(td.depth))
  
  td_not_exist <- is.na(td.depth)
  
 
  for (ii in (1:length(td.depth))){
    idx = !is.na(temp[ii,])
    temp_data = as.numeric(temp[ii,idx])
    total_temp[ii]<-max(sum(temp_data)/length(temp_data),4)
    if(!td_not_exist[ii]){
      td_idx <- max(which(td.depth[ii]>=depth_data))
      epil_temp[ii] <- mean(temp_data[1:td_idx])
      if (td_idx >= length(temp_data)) {td_idx = length(temp_data)-1}
      hypo_temp[ii] <- mean(temp_data[(td_idx+1):length(temp_data)])
      # if (is.na(hypo_temp[ii]) && !is.na(epil_temp[ii])){
      #   break
      # }
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
calc_vol_total<-function(H,A){
  if (length(H)==1){
    vol_total <- 1/3.0 * A * H
  }else{
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
  
  td_not_exist<-is.na(td.depth)
  for (ii in 1:length(td.depth)){
    if(!td_not_exist[ii]){
      h_idx <- min(which(td.depth[ii]>=H))
      approx_td.area<-approx(H, A, c(0, td.depth[ii]))$y[2]
      if (is.na(approx_td.area)){approx_td.area = min(A)}
      H_with_td<-c(td.depth[ii],H[h_idx:length(H)])
      A_with_td<-c(approx_td.area,A[h_idx:length(A)])
      vol_data[ii,1] <- trapz(rev(H_with_td),rev(A_with_td))##epil
      vol_data[ii,2] <- vol_total - vol_data[ii,1]##hypo
      if(vol_data[ii,2] <= 0) {
        vol_data[ii,2]= min(A) * 0.5
        vol_data[ii,1] = vol_total - vol_data[ii,2]
      }
      # if (is.na( vol_data[ii,1]) && !is.na( vol_data[ii,2])){
      #   break
      # }
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
  td_area <- approx(H, A, td.depth)$y
  surf_area <- rep(max(A), length(grd.info$datetime))
  temp_out<-calc_epil_hypo_temp(wtemp,td.depth,H)
  vol_total<- calc_vol_total(H,A)
  vol<-calc_epil_hypo_vol(H,A,td.depth,vol_total)
  return(data.frame(datetime = as.POSIXct(grd.info$datetime),td.depth,t.epil = temp_out$t_epil,t.hypo=temp_out$t_hypo,
                          t.total = temp_out$t_total,vol_total,vol,
                    td_area, surf_area))
}



#' Calculate the mass of dissolved oxgen in epil and hypolimnion(during stratified period)
#' and the total dissovled oxygen in the lake
#' @param wtemp matrix; Water temperatures (rows correspond to time, cols to depth)
#' @return list of datetimes and depths
#' @export
#' 
calc_do<-function(input.values,fsed_stratified_epi,fsed_stratified_hypo,fsed_not_stratified,nep_stratified,nep_not_stratified,
                  min_stratified, min_not_stratified, wind = NULL, khalf = NULL, startdate = NULL, enddate = NULL){
  
  ##initialize matrix
  o2_data <- matrix(NA, nrow = length(input.values$td.depth), ncol = 18) 
  colnames(o2_data) <- c("o2_epil","o2_hypo","o2_total",
                         'Fsed_total', "NEP_total", "Fatm_total", "Mineral_total",
                         'Fsed_epi', "NEP_epi", 'Fatm_epi', 'Entrain_epi',
                         'Fsed_hypo', 'Mineral_hypo','Entrain_hypo',
                         'sat_o2_epil', 'sat_o2_hypo', 'sat_o2_total', 'massbal')
  o2_data <- as.data.frame(o2_data)
  
  init_o2sat <- o2.at.sat.base(temp=input.values$t.total[1],altitude = 300)*1000 # returns mg O2/L 
  
  theta<-1.08
  
   if (is.null(khalf)){
     khalf <- 4800.
   }
  
  td_not_exist <- is.na(input.values$td.depth)
  
  if (is.null(startdate)){
    startdate = 1
  }
  if (is.null(enddate)){
    enddate =  length(input.values$td.depth)
  }
  
  o2_data$o2_total[startdate] <- init_o2sat # returns mg O2 (m3 = 1000 L)
  
  for(day in (startdate + 1):enddate){
    
    K600<-ifelse(is.null(wind), k.cole.base(2),k.cole.base(wind[day]))
    
    ## not stratified period, only consider o2_total dynamics
    if(td_not_exist[day]){
      theta_total <- theta^(input.values$t.total[day-1]-20)
      
      NEP <- (nep_not_stratified * theta_total *(o2_data[day-1,"o2_total"])/(khalf + o2_data[day-1,"o2_total"]))
      
      MINER <- (min_not_stratified * theta_total*(o2_data[day-1,"o2_total"])/(khalf + o2_data[day-1,"o2_total"]))
      
      kO2 <- k600.2.kGAS.base(k600=K600,temperature=input.values$t.total[day],gas='O2') # velocity value m/d?
      o2sat<-o2.at.sat.base(temp=input.values$t.total[day],altitude = 300)*1000 # mg O2/L -> mg/m3

      Fatm <- (kO2*(o2sat - o2_data[day-1,"o2_total"])/max(H))
      
      Fsed <- (fsed_not_stratified * (o2_data[day-1,"o2_total"])/(khalf + o2_data[day-1,"o2_total"]) * theta_total * input.values$surf_area[day-1] / input.values$vol_total[day-1])
      
      o2_data[day,"o2_total"] <-( o2_data[day-1,"o2_total"]+ valid((Fsed + NEP + Fatm + MINER),o2_data[day-1,"o2_total"]) ) * ( input.values$vol_total[day-1]/ input.values$vol_total[day])
      
      if (valid((Fsed + NEP + Fatm + MINER), o2_data[day-1,"o2_total"]) != (Fsed + NEP + Fatm + MINER)){
        Fsed <- Fsed/(Fsed + NEP + Fatm + MINER) * o2_data[day-1,"o2_total"]
        NEP <- NEP/(Fsed + NEP + Fatm + MINER) * o2_data[day-1,"o2_total"]
        Fatm <- Fatm/(Fsed + NEP + Fatm + MINER)* o2_data[day-1,"o2_total"]
        MINER <- MINER/(Fsed + NEP + Fatm + MINER) * o2_data[day-1,"o2_total"]
      }
      
      o2_data[day,"Fsed_total"] <- Fsed 
      o2_data[day,"NEP_total"] <- NEP
      o2_data[day,"Fatm_total"] <- Fatm
      o2_data[day,"Mineral_total"] <- MINER
      o2_data[day,"sat_o2_total"] <- (100. * o2_data[day,"o2_total"] )/ o2sat

    }
    # the day it turns to stratified, need to reassign the o2 to hypo and epil
    else if(is.na(input.values$td.depth[day-1])){
      o2_data[day,"o2_epil"] <- (o2_data[day-1,"o2_total"]*input.values$vol_epil[day])/input.values$vol_total[day]
      o2_data[day,"o2_hypo"] <- (o2_data[day-1,"o2_total"]*input.values$vol_hypo[day])/input.values$vol_total[day]
      o2_data[day,"o2_total"] <- (o2_data[day,"o2_hypo"]* input.values$vol_hypo[day] + o2_data[day,"o2_epil"] * input.values$vol_epil[day])/input.values$vol_total[day]
      # o2_data[day,"o2_epil"] <- o2_data[day-1,"o2_total"]#((o2_data[day-1,"o2_total"]*input.values$vol_epil[day])/input.values$vol_total[day]) #/ input.values$vol_epil[day] #/input.values$t.total[day]*input.values$vol_epil[day]
      # o2_data[day,"o2_hypo"] <- o2_data[day-1,"o2_total"]#((o2_data[day-1,"o2_total"]*input.values$vol_hypo[day])/input.values$vol_total[day])  #/ input.values$vol_hypo[day] #/input.values$t.total[day]*input.values$vol_hypo[day]
      # o2_data[day,"o2_total"]<- (o2_data[day,"o2_epil"]+o2_data[day,"o2_hypo"]) /2 #/ input.values$vol_total[day] 
      
      o2sat<-o2.at.sat.base(temp=input.values$t.epil[day],altitude = 300)*1000 
      o2_data[day,"sat_o2_epil"] <- (100. * o2_data[day,"o2_epil"])/ o2sat
      o2sat<-o2.at.sat.base(temp=input.values$t.hypo[day],altitude = 300)*1000 
      o2_data[day,"sat_o2_hypo"] <- (100. * o2_data[day,"o2_hypo"]  )/ o2sat
      o2sat<-o2.at.sat.base(temp=input.values$t.total[day],altitude = 300)*1000 
      o2_data[day,"sat_o2_total"] <- (100. * o2_data[day,"o2_total"] )/ o2sat
    }else{
      theta_epil <- theta^(input.values$t.epil[day-1]-20)
      theta_hypo <-  theta^(input.values$t.hypo[day-1]-20)
      
      NEP_epil <- (nep_stratified * theta_epil *(o2_data[day-1,"o2_epil"])/(khalf + o2_data[day-1,"o2_epil"]))
              
      kO2_epil <- k600.2.kGAS.base(k600=K600,temperature=input.values$t.epil[day],gas='O2')
      o2sat_epil<-o2.at.sat.base(temp=input.values$t.epil[day],altitude = 300)*1000

      Fatm_epil <- (kO2_epil*(o2sat_epil-o2_data[day-1,"o2_epil"] )/input.values$td.depth[day-1])
     
      
      volumechange_epi = input.values$vol_epil[day]-input.values$vol_epil[day-1]  
      volumechange_epi_proportion =  volumechange_epi/input.values$vol_epil[day-1] 
      if (volumechange_epi_proportion >= 0){
        x_do <- o2_data[day - 1,"o2_hypo"] #( o2_data[day - 1,"o2_hypo"] * abs(volumechange_epi)) / input.values$vol_hypo[day-1]
      } else {
        x_do <- o2_data[day - 1,"o2_epil"] #( o2_data[day - 1,"o2_epil"] * abs(volumechange_epi)) / input.values$vol_epil[day-1]
      }
      
      Fepi <-  (volumechange_epi_proportion* x_do)

      Fsed_epi <- (fsed_stratified_epi * (o2_data[day-1,"o2_epil"])/(khalf + o2_data[day-1,"o2_epil"]) * theta_epil * input.values$surf_area[day-1] /input.values$vol_epil[day-1])

      o2_data[day,"o2_epil"] <- (o2_data[day-1,"o2_epil"] + valid((NEP_epil+Fatm_epil + Fepi + Fsed_epi) , o2_data[day-1,"o2_epil"]) ) * (input.values$vol_epil[day-1]/input.values$vol_epil[day])
      
      if (valid((NEP_epil+Fatm_epil + Fepi + Fsed_epi), o2_data[day-1,"o2_epil"]) != (NEP_epil+Fatm_epil + Fepi + Fsed_epi)){
        Fsed_epi <- Fsed_epi/(NEP_epil+Fatm_epil + Fepi + Fsed_epi) * o2_data[day-1,"o2_epil"]
        NEP_epil <- NEP_epil/(NEP_epil+Fatm_epil + Fepi + Fsed_epi) * o2_data[day-1,"o2_epil"]
        Fatm_epil <- Fatm_epil/(NEP_epil+Fatm_epil + Fepi + Fsed_epi) * o2_data[day-1,"o2_epil"]
        Fepi <- Fepi/(NEP_epil+Fatm_epil + Fepi + Fsed_epi) * o2_data[day-1,"o2_epil"]
      }
      
      o2_data[day,"Fsed_epi"] <- Fsed_epi 
      o2_data[day,"NEP_epi"] <- NEP_epil
      o2_data[day,"Fatm_epi"] <- Fatm_epil
      o2_data[day,"Entrain_epi"] <- Fepi

      volumechange_hypo = input.values$vol_hypo[day]-input.values$vol_hypo[day-1]  #in m^3
      volumechange_hypo_proportion =  volumechange_hypo/input.values$vol_hypo[day-1] 
      
      if (volumechange_hypo_proportion >= 0){
        x_do <-  o2_data[day - 1,"o2_epil"] #( o2_data[day - 1,"o2_epil"] * abs(volumechange_hypo)) / input.values$vol_epil[day-1]
      } else {
        x_do <- o2_data[day - 1,"o2_hypo"] #( o2_data[day - 1,"o2_hypo"] * abs(volumechange_hypo)) / input.values$vol_hypo[day-1]
      }

      Fhypo <- (volumechange_hypo_proportion* x_do ) 

      MINER_hypo <- (min_stratified * theta_hypo  * (o2_data[day-1,"o2_hypo"])/(khalf + o2_data[day-1,"o2_hypo"]))
      
      
      Fsed_hypo <- (fsed_stratified_hypo * (o2_data[day-1,"o2_hypo"])/(khalf + o2_data[day-1,"o2_hypo"]) * theta_hypo * input.values$td_area[day-1] / input.values$vol_hypo[day-1])
      

      o2_data[day,"o2_hypo"] <- (o2_data[day-1,"o2_hypo"] + valid((Fhypo + Fsed_hypo + MINER_hypo), o2_data[day-1,"o2_hypo"]) ) * (input.values$vol_hypo[day-1]/input.values$vol_hypo[day])

      if (valid((Fhypo + Fsed_hypo + MINER_hypo), o2_data[day-1,"o2_hypo"]) != (Fhypo + Fsed_hypo + MINER_hypo)){
        Fsed_hypo <- Fsed_hypo/(Fhypo + Fsed_hypo + MINER_hypo) * o2_data[day-1,"o2_hypo"]
        Fhypo <- Fhypo/(Fhypo + Fsed_hypo + MINER_hypo) * o2_data[day-1,"o2_hypo"]
        MINER_hypo <- MINER_hypo/(Fhypo + Fsed_hypo + MINER_hypo) * o2_data[day-1,"o2_hypo"]
      }


      o2_data[day,"o2_total"] <- (o2_data[day,"o2_hypo"]* input.values$vol_hypo[day] + o2_data[day,"o2_epil"] * input.values$vol_epil[day])/input.values$vol_total[day]
      o2_data[day,"Fsed_hypo"] <- Fsed_hypo 
      o2_data[day,"Mineral_hypo"] <- MINER_hypo
      o2_data[day,"Entrain_hypo"] <- Fhypo

      
      o2sat<-o2.at.sat.base(temp=input.values$t.epil[day],altitude = 300)*1000 
      o2_data[day,"sat_o2_epil"] <- (100. * o2_data[day,"o2_epil"] )/ o2sat
      o2sat<-o2.at.sat.base(temp=input.values$t.hypo[day],altitude = 300)*1000 
      o2_data[day,"sat_o2_hypo"] <- (100. * o2_data[day,"o2_hypo"]  )/ o2sat
      o2sat<-o2.at.sat.base(temp=input.values$t.total[day],altitude = 300)*1000 
      o2_data[day,"sat_o2_total"] <- (100. * o2_data[day,"o2_total"] )/ o2sat
      
      # mass balance
      mass_thr <- sum(c(o2_data$o2_epil[day-1] * input.values$vol_epil[day-1], o2_data$o2_hypo[day-1] * input.values$vol_hypo[day-1])) +
       ( (NEP_epil + Fatm_epil + Fepi + Fsed_epi) * (input.values$vol_epil[day-1]) + (MINER_hypo + Fhypo + Fsed_hypo ) * (input.values$vol_hypo[day-1]))
      mass_balance <-  sum(c(o2_data$o2_epil[day]*input.values$vol_epil[day], o2_data$o2_hypo[day]* input.values$vol_hypo[day])) - mass_thr
      o2_data[day,"massbal"] <- mass_balance
      # print(mass_balance)
    }
    if (is.na(o2_data[day, "o2_total"])) {
      break
      print('RED ALERT!')
    }
  }

  return (o2_data) 
}

#' Compares simulated against measured variables
#' @param obs matrix; observational data
#' @param input.values data frame, simulaedt data
#' @return matched data frame of observed and simulated data
#' @export
#' 
compare_predict_versus_observed<-function(obs,input.values){
  ndate <- c()
  ndate <- unique(obs$ActivityStartDate)
  
  test_data<- matrix(NA, nrow = length(ndate), ncol = 7)  
  colnames(test_data) <- c("day","observed_epil_do","predict_epil_do","observed_hypo_do","predict_hypo_do",
                           "observed_total_do","predict_total_do")
  test_data <- as.data.frame(test_data)
  test_data$day<-as.POSIXct(test_data$day)
  
  td_not_exist<- is.na(input.values$td.depth)
  
  for(jj in 1:length(ndate)){
    test_data$day[jj]<-ndate[jj]
    ##identify the index of the day from input values tested on the output
    kk = which(year(test_data$day[jj]) == year(input.values$datetime) & 
                 yday(test_data$day[jj])== yday(input.values$datetime))
    if(length(kk)!=1){
      next
    }
    if(td_not_exist[kk]){
      test_data$predict_total_do[jj] <- input.values$o2_total[kk]
      
      obs_data<-  obs %>% filter(ActivityStartDate == ndate[jj]) %>% 
        slice(which.max(ActivityDepthHeightMeasure.MeasureValue)) 
      test_data$observed_total_do[jj] <- obs_data$ResultMeasureValue
      
      next
    }else{
      td <- input.values$td.depth[kk]
      test_data$predict_epil_do[jj] <- input.values$o2_epil[kk]/input.values$vol_epil[kk]/1000
      test_data$predict_hypo_do[jj] <- input.values$o2_hypo[kk]/input.values$vol_hypo[kk]/1000
      # here the mean functions will weigh the result, which will give different values in contrast to our integrated approach
      
      #shallowest layer in epi;
      obs_data_epil<- obs %>% filter(ActivityStartDate == ndate[jj]) %>% 
        filter(ActivityDepthHeightMeasure.MeasureValue<=td) %>%  slice(which.min(ActivityDepthHeightMeasure.MeasureValue)) 
      if(!nrow(obs_data_epil)==0){
      test_data$observed_epil_do[jj]<- obs_data_epil$ResultMeasureValue
      }
      
      #deepset layer in hypo
      obs_data_hypo<- obs %>% filter(ActivityStartDate == ndate[jj]) %>% 
        filter(ActivityDepthHeightMeasure.MeasureValue>td) %>%  slice(which.max(ActivityDepthHeightMeasure.MeasureValue)) 
      if(!nrow(obs_data_hypo)==0){
      test_data$observed_hypo_do[jj]<-obs_data_hypo$ResultMeasureValue
      }
    }
  }
  return(test_data)
}

#' Calculates RMSE for epilimnion
#' @param test_data matrix; Matched observed to simulated data
#' @return double value of RMSE
#' @export
#' 
calc_rmse_epil<-function(test_data){
  predicted <- test_data$predict_epil_do 
  actual <- test_data$observed_epil_do
  return (sqrt(mean((predicted-actual)**2,na.rm = TRUE))) # RMSE
}


#' Calculates RMSE for hypolimnion
#' @param test_data matrix; Matched observed to simulated data
#' @return double value of RMSE
#' @export
#' 
calc_rmse_hypo<-function(test_data){
  predicted <- test_data$predict_hypo_do
  actual <- test_data$observed_hypo_do
  return (sqrt(mean((predicted-actual)**2,na.rm = TRUE))) # RMSE
}

#' Calculates RMSE for total lake
#' @param test_data matrix; Matched observed to simulated data
#' @return double value of RMSE
#' @export
#' 
calc_rmse <- function(test_data){
  predicted <- rbind(test_data$predict_hypo_do, test_data$predict_epil_do)
  actual <- rbind(test_data$observed_hypo_do,test_data$observed_epil_do)
  return (sqrt(mean((predicted-actual)**2,na.rm = TRUE))) # RMSE
}


#' Calculates RMSE for total lake
#' @param test_data matrix; Matched observed to simulated data
#' @return double value of RMSE
#' @export
#' 
calc_fit <- function(input.values, proc.obs){
  obs <- cbind(proc.obs[3,], proc.obs[4,])
  mod <- cbind(input.values$o2_epil[proc.obs[1,]]/1000,
                  input.values$o2_hypo[proc.obs[1,]]/1000)
  return (sqrt(mean((obs-mod)**2,na.rm = TRUE))) # RMSE
}

#' Calculates RMSE for total lake
#' @param test_data matrix; Matched observed to simulated data
#' @return double value of RMSE
#' @export
#' 
optim_do <- function(p, input.values, fsed_not_stratified = 0.0002, nep_not_stratified = 0.0, min_not_stratified = 0.0,
                     wind = NULL, proc.obs,verbose,  startdate = NULL, enddate = NULL){
  p <- lb+(ub - lb)/(10)*(p)
  
  o2<- calc_do(input.values = input.values,fsed_stratified_epi = p[1],
               fsed_stratified_hypo = p[2],
               fsed_not_stratified =fsed_not_stratified,
               nep_stratified = p[3],
               nep_not_stratified = nep_not_stratified,
               min_stratified = p[4],
               min_not_stratified =min_not_stratified, wind =wind, khalf = p[5],
               startdate = startdate, enddate = enddate)
  
  input.values$o2_epil <- o2[,"o2_epil"]
  input.values$o2_hypo <- o2[,"o2_hypo"]
  input.values$o2_total <- o2[,"o2_total"]
  
  input.values$year <- year(input.values$datetime)
  input.values$doy <- yday(input.values$datetime)
  
  
  # test_data<-compare_predict_versus_observed(obs,input.values) 
  # 
  # fit = calc_rmse(test_data)
  
  fit = calc_fit(input.values = input.values, proc.obs)
  
  print(paste(round(p[1],5),round(p[2],5),round(p[3],2),round(p[4],2),round(p[5],2),'with RMSE: ',round(fit,4)))
  
  return(fit)
}

#' check whether the flux is valid and return a valid flux
#' @param flux to be checked
#' @param pool cuurent existing o2 in the lake
#' @return a valid flux
#' @export
#' 
valid<-function(flux,pool){
  if(abs(flux)>abs(pool)){
    if(flux<0){
      return (-pool)
    }else{
      return (flux)
    }
  
  }else{
    return (flux)
  }
}

#' preprocesses observed data and area-weighs them
#' @param obs observed data
#' @param pool input matrix of for instance thermocline depth
#' @param H depths
#' @param A areas
#' @return matched and weighted-averaged data
#' @export
#' 
preprocess_obs <- function(obs, input.values, H, A){
  deps <- seq(round(max(H),4), round(min(H),4), by = -0.5)
  
  if (max(H) > max(deps)){
    deps <- c(max(H), deps)
  }
  if (min(H) < min(deps)){
    deps <- c(deps, min(H))
  }
  
  areas <- approx(round(H,4), round(A,4), round(deps,4))$y
  
  apprObs <- matrix(NA, nrow= length(deps), ncol =length(unique(zoo::as.Date(obs$ActivityStartDate))))
  ts.apprObs <- matrix(NA, nrow= 3, ncol =length(unique(zoo::as.Date(obs$ActivityStartDate))))
  idx <- c()
  for (jj in unique(zoo::as.Date(obs$ActivityStartDate))){

    idy =  (match(zoo::as.Date(obs$ActivityStartDate),zoo::as.Date(jj)))
    idy <- which(!is.na(idy))
    dat <- obs[idy,]
  
    
    if (sd(dat$ActivityDepthHeightMeasure.MeasureValue) == 0 | length(dat$ActivityDepthHeightMeasure.MeasureValue) <= 1 |
        length(na.omit( round(dat$ResultMeasureValue,2))) <= 1){
      next} else {
           
        if (max(deps) > max(round(dat$ActivityDepthHeightMeasure.MeasureValue,2))) {
          dat <- rbind(dat, data.frame('ActivityStartDate' = zoo::as.Date(jj),
                                'ActivityDepthHeightMeasure.MeasureValue' = max(deps),
                                'ResultMeasureValue' = dat$ResultMeasureValue[nrow(dat)]))
        }
        
        if (any(is.na(approx(round(dat$ActivityDepthHeightMeasure.MeasureValue,2), round(dat$ResultMeasureValue,2),
                             deps)$y))){
          intvec <- (approx(round(dat$ActivityDepthHeightMeasure.MeasureValue,2), round(dat$ResultMeasureValue,2),
                                     deps)$y)
          intvec[ which(is.na(intvec))] <- intvec[( which(is.na(intvec)))+1]
          } else {
           intvec <- approx(round(dat$ActivityDepthHeightMeasure.MeasureValue,2), round(dat$ResultMeasureValue,2),
                                                        deps)$y
                                     }
    apprObs[,match(jj, unique(zoo::as.Date(obs$ActivityStartDate)))] <-intvec
    
    if (zoo::as.Date(jj) < min(zoo::as.Date(input.values$datetime)) | zoo::as.Date(jj) > max(zoo::as.Date(input.values$datetime))){
      next 
    } else {
      idx <- append(idx,  match(zoo::as.Date(jj), zoo::as.Date(input.values$datetime)))
    idz <-  which(zoo::as.Date(jj) == zoo::as.Date(input.values$datetime))
    if (is.na(input.values$td.depth[abs(idz)])){
      dz.areas <- (1*areas)/sum(areas, na.rm= TRUE)#(areas - min(areas)) / (max(areas) - min(areas))
      ts.apprObs[1, match(zoo::as.Date(jj), unique(zoo::as.Date(obs$ActivityStartDate)))] <- weighted.mean(apprObs[,match(zoo::as.Date(jj), unique(as.Date(obs$ActivityStartDate)))],
                                                                                 dz.areas, na.rm = TRUE)
    } else{
      z.td <- which(abs(input.values$td.depth[abs(idz)] - deps) == (min(abs(input.values$td.depth[abs(idz)] - deps))[1]))
      dz.hypo <- (1*areas[1:z.td])/sum(areas[1:z.td])
      dz.epi <- (1*areas[(z.td+1):length(areas)])/sum(areas[(z.td+1):length(areas)])
      
      ts.apprObs[2, match(jj, unique(zoo::as.Date(obs$ActivityStartDate)))] <- weighted.mean(apprObs[(z.td+1):length(areas),match(zoo::as.Date(jj), unique(as.Date(obs$ActivityStartDate)))],
                                                                                        dz.epi, na.rm = TRUE)
      
      ts.apprObs[3, match(jj, unique(zoo::as.Date(obs$ActivityStartDate)))] <- weighted.mean(apprObs[1:z.td,match(jj, unique(zoo::as.Date(obs$ActivityStartDate)))],
                                                                                        dz.hypo, na.rm = TRUE)
    }
    }
  }
  }
  check.na <- c()
  for (p in 1:ncol(ts.apprObs)){
    if (all(is.na(ts.apprObs[,p]))){
      check.na <- append(check.na, p)
    }
  }
  # idx <- match(unique(as.Date(obs$ActivityStartDate)), as.Date(input.values$datetime))
  if (is.null(check.na)){
    return(rbind(idx, ts.apprObs))
  } else {
    return(rbind(idx, ts.apprObs[,-(check.na)]))
  }
 
}

#' preprocesses observed data and area-weighs them w/o interpolation
#' @param obs observed data
#' @param input.values input matrix of for instance thermocline depth
#' @param H depths
#' @param A areas
#' @return matched and weighted-averaged data
#' @export
#' 
weigh_obs <- function(obs, input.values, H, A){
  
  data_long <- obs %>% arrange(ActivityStartDate)
  data_long$Area <- approx(H, A, data_long$ActivityDepthHeightMeasure.MeasureValue)$y
  
  idx <- match(zoo::as.Date(data_long$ActivityStartDate), zoo::as.Date(input.values$datetime))
  data_long$Layer <- data_long$ActivityDepthHeightMeasure.MeasureValue <= input.values$td.depth[idx]
  
  data_long$Layer[which(data_long$Layer == TRUE)] = 'EPILIMNION'
  data_long$Layer[which(data_long$Layer == FALSE)] = 'HYPOLIMNION'
  data_long$Layer[which(is.na(data_long$Layer))] = 'TOTAL'
  
  data_long$WeightValue <- rep(NA, nrow(data_long))
  weight_obs <- matrix(NA, nrow = 4, ncol = length(unique(zoo::as.Date(data_long$ActivityStartDate))))
  
  for (ii in unique(zoo::as.Date(data_long$ActivityStartDate))){
    # print(zoo::as.Date(ii))
    idx <- which(zoo::as.Date(ii) == zoo::as.Date(data_long$ActivityStartDate))
    idz <- match(zoo::as.Date(ii), zoo::as.Date(input.values$datetime))
    thdepth <- input.values$td.depth[idz]
    data <- data_long[idx, ]
    
    weight_obs[1, match(ii, unique(zoo::as.Date(data_long$ActivityStartDate)))] <- idz
    
    if (all(data$Layer == 'TOTAL')){
      total_areas <- approx(H, A, seq(from = max(H), to = 0, by = -0.5))$y
      perc <- (1 * data$Area) / max(total_areas)
      data_long$WeightValue[idx] <- data$ResultMeasureValue * perc
      data$WeightValue<- data$ResultMeasureValue * perc
      
      # print(paste(
      #   mean(data$WeightValue[which(data$Layer == 'TOTAL')])))
      
      weight_obs[2, match(ii, unique(zoo::as.Date(data_long$ActivityStartDate)))] <- mean(data$WeightValue[which(data$Layer == 'TOTAL')])
    } else {
      idy = which(data$Layer == 'EPILIMNION')
      epi_areas <- approx(H, A, seq(from = round(thdepth,1), to = 0, by = -0.5))$y
      epi_perc <- (1 * data$Area[idy]) / max(epi_areas)
      
      idt = which(data$Layer == 'HYPOLIMNION')
      hypo_areas <- approx(H, A, seq(from = max(H), to = round(thdepth,1), by = -0.5))$y
      hypo_perc <- (1 * data$Area[idt]) / max(hypo_areas)
      
      data_long$WeightValue[idx] <- c(data$ResultMeasureValue[idy] * epi_perc,
                                      data$ResultMeasureValue[idt] * hypo_perc)
      data$WeightValue <- c(data$ResultMeasureValue[idy] * epi_perc,
                            data$ResultMeasureValue[idt] * hypo_perc)
      
      # print(paste(
      #   mean(data$WeightValue[which(data$Layer == 'EPILIMNION')]),
      #   mean(data$WeightValue[which(data$Layer == 'HYPOLIMNION')])
      # ))
      weight_obs[3, match(ii, unique(zoo::as.Date(data_long$ActivityStartDate)))] <- mean(data$WeightValue[which(data$Layer == 'EPILIMNION')])
      weight_obs[4, match(ii, unique(zoo::as.Date(data_long$ActivityStartDate)))] <- mean(data$WeightValue[which(data$Layer == 'HYPOLIMNION')])
      
    }
  }
  return(list(data_long, weight_obs))
}

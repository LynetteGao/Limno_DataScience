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
  
  for (ii in 1:length(cbuoy.depth)){
    idx = !is.na(temp[ii,])
    dens_data = dens[ii,idx]
    dens.diff = rev(dens_data)[1] - dens_data[1]
    
    if (min(temp[ii,],na.rm=TRUE) > 4 && abs(dens.diff) > 0.1){
    cbuoy.depth[ii] <- center.buoyancy(temp[ii,], as.numeric(grd.info$depth))
    }
  }
  
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
    
    if (min(temp_data) > 4 && abs(dens.diff) > 0.05){
      d.idx = which(td.diff == max(td.diff)[1])
      td.depth[ii] = depth_data[d.idx]
    } else {td.depth[ii] = NA}
  }
  years = year(grd.info$datetime)
  
  for (ii in unique(years)){
    idx = which( years %in% ii)
    ydata = td.depth[idx]
    # ydata[ydata[1:200] > 10] = NA

    
   

    
     if (all(is.na(ydata)) == FALSE){
       ydata = na.approx(ydata)
       testdata = data.frame('depth' = (ydata), 'time' = seq(1,length(ydata)))
    #   # exponential.model <- lm(log(depth)~ time, data = testdata)
       poly.model <- lm(depth ~ poly(time,3), data = testdata)
    #   # expdata = predict(exponential.model)
       polydata = predict(poly.model)
       td.depth[idx[min(which(!is.na(td.depth[idx])))]:idx[max(which(!is.na(td.depth[idx])))]] =as.numeric(polydata)
    #
    # # movdata = movavg(na.approx(ydata), 7, type = 's')
    # # td.depth[idx] = NA
    # # td.depth[idx[min(which(!is.na(ydata)))]:idx[max(which(!is.na(ydata)))]] =movdata
     }

    
  }
  td.depth[td.depth < 0] = NA
  
  return(cbuoy.depth)
}


track_td_depth <- function(wtemp){

  grd.info <- extract_time_space(wtemp)
  temp <- as.matrix(wtemp[,-c(1)])
  dens <- calc_dens(temp)
  td.depth <- rep(NA, length(grd.info$datetime))

  td.diff <- matrix(NA, nrow=length(grd.info$depth), ncol=length(td.depth))
  sign.td.diff <- matrix(NA, nrow=length(grd.info$depth), ncol=length(td.depth))
  check.td <- rep(NA, length(td.depth),1)
  for (ii in 1:length(td.depth)){
    # simplify data
    idx = which(!is.na(temp[ii,]))
    temp_data = temp[ii,idx]
    dens_data = dens[ii,idx]
    depth_data = as.double(grd.info$depth[idx])

    # td.diff <- rep(NA, length(depth_data))
    # forward differencing d rho / d z
    td.diff[idx,ii] <- (lead(dens_data) - dens_data)/
      (lead(depth_data) - depth_data)
    # backward differencing d rho / d z
    td.diff[length(td.diff[idx,ii]),ii] <- (dens_data[nrow(td.diff)] - dens_data[nrow(td.diff)-1])/
      (depth_data[nrow(td.diff)] - depth_data[nrow(td.diff)-1])

    tz <- which(td.diff[idx,ii] > 0.1)
    sign.td.diff[tz,ii] <- td.diff[tz,ii]

    # find density difference between epilimnion and hypolimnion
    dens.diff = rev(dens_data)[1] - dens_data[1]

    if (min(temp_data) > 4 && abs(dens.diff) > 0.1){
      check.td[ii] <- 1
      # d.idx = which(td.diff == max(td.diff))[1]
      # td.depth[ii] = depth_data[d.idx]
    } else {
      check.td[ii] <- NA
      # td.depth[ii] = NA
      }
  }
  
  years <- year(grd.info$datetime)
  idy <- which(years %in% '1992') 
  
  new.df <- reshape2::melt(t(td.diff[,idy]))
  dens.df <- reshape2::melt(dens[idy,])
  both.df <- data.frame('Time' = new.df$Var1, 'Depth' = new.df$Var2,
                        'Gradient' = new.df$value, 'Density' = dens.df$value)
  
  contour(x = as.numeric(grd.info$datetime[idy]), y = as.numeric(grd.info$depth), z=t(td.diff[,idy]))
  contour(x = as.numeric(grd.info$datetime[idy]), y =  as.numeric(grd.info$depth), z=(dens[idy,]))
  
  library(gganimate)
  for (ii in 1:length(idy)){
  
  p <- ggplot(subset(both.df, Time == ii), aes(x = Density, y = Depth, col = Gradient, size = Gradient)) +
    #geom_line() +
    geom_point(show.legend = FALSE,alpha = 0.7) +
    scale_colour_gradientn(colours=topo.colors(10)) +
    scale_y_reverse()+
    xlim(996.5,1000)+
     annotate('text', x = 997.1, y = 10, label = paste('Day',ii), size = 3) +
    theme_bw();p
  ggsave(file = paste0('td_',ii,'.png'), p, dpi=300, width = 200, height = 200, units = 'mm')
  }
  p <- ggplot(both.df, aes(x = Density, y = Depth, col = Gradient, size = Gradient)) +
    #geom_line() +
    geom_point(show.legend = FALSE,alpha = 0.7) +
    scale_colour_gradientn(colours=topo.colors(10)) +
    scale_y_reverse()+
    theme_bw();p
 p <- p + transition_time(Time) +
    labs(title = 'Day: {frame_time}')
  anim_save('td_find.gif',p)
  
  


  

  for (jj in unique(years)){
    idy <- which(years %in% jj)
    data = td.diff[,idy]
    depth_data = as.double(grd.info$depth)
    library(rlist)
    td.pks <- c()#rep( list(list()), (ncol(data)) )
    for (kk in 1:ncol(data)){
     # td.pks[[kk]] <- list.append(td.pks[[kk]], findpeaks(data[,kk]))
      if (!is.null(findpeaks(data[,kk]))){
        td.pks <- rbind(td.pks, cbind(findpeaks(data[,kk]), rep(kk, nrow(findpeaks(data[,kk])))))
      } else {
        td.pks <- rbind(td.pks, matrix(NA,ncol=5,nrow=1))
      }

    }
    periodint <- as.double(na.contiguous(td.pks[,5]))
    idx = which((td.pks[,5]) %in% periodint)
    td.pks <- as.data.frame(td.pks)
    colnames(td.pks) <- c('grad','pos','min','max','doy')
    td.df <- td.pks %>%
      filter(doy >= min(periodint) & doy <= max(periodint))

    td.depths <- matrix(NA, ncol=10e5,nrow=nrow(data))
    test = NULL
    timtst = NULL
    for (kk in (unique(td.df$doy))){
      dat = subset(td.df, doy == kk)
        for (oo in 1:nrow(dat)){

          ds = TRUE
          hh=kk

          bub = c(dat[oo,2])
          ub = c(kk)
          if (bub < 45 & !is.na(bub)){
          while (ds){

            p1 <- subset(td.df, doy == hh+1)
            minimapks <-(p1[,2]-dat[oo,2])
            locpks <- which(abs(minimapks) < 5)
            if (length(locpks) > 1){
              locpks <- which.max(p1[locpks,1])
            }

            if (length(locpks) == 0) {
              ds = FALSE
            } else {

            dat = p1
            hh=hh+1
            oo = locpks

            bub = append(bub, p1[locpks,2])
            ub = append(ub, hh)}
          }
          pop = rep(NA, 500)
          pop[1:length(bub)] = bub
          bob = rep(NA, 500)
          bob[1:length(ub)] = ub
          if (is.null(test)){
            test = pop
            timtst = bob
          } else {
            test = cbind(test, pop)
            timtst = cbind(timtst, bob)
          }
}

        }
    }
    l = c()
    for (ii in 1:ncol(test)){
      l = append(l, sum(!is.na(test[,ii])))
    }
    sort(l)
    test[,which.max(l)]
    timtst[,which.max(l)]

    plot(na.omit(timtst[,which.max(l)]), depth_data[na.omit(test[,which.max(l)])])


    td.depth <- rep( list(list()), (ncol(data)-1) )#matrix(NA,ncol= ncol(data),nrow= (1000))
    for (kk in 1:(ncol(data)-1)){
        currentpks <- findpeaks(data[,kk])
        nextpks <- findpeaks(data[,kk+1])
        if (!is.null(nrow(currentpks))){
        for (pp in 1:nrow(currentpks)){
          minimapks <-(nextpks[,2]-currentpks[pp,2])
          locpks <- which(abs(minimapks) < 2)
          if (length(locpks) > 0){
          if (length(locpks) > 1){
            locpks <- which.max(nextpks[1,locpks])
          }
          if (( nextpks[locpks,1] - currentpks[pp,1]) > -0.001){
            td.depth[[kk]] <- list.append(td.depth[[kk]], depth_data[currentpks[pp,2]])
          }
          }
        }
        }
    }
    test = td.depth[!is.na(td.depth)]
    while (any(diff(test[!is.na(test)]) < 0)){
      test <- test[!is.na(test)]
      td.depth.diff <- (diff(test))
      test[td.depth.diff <0] <- NA
      test[td.depth.diff > 10] <- NA
    }
    td.depth.diff <- (diff(test))
    test[td.depth.diff <0] <- NA


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
calc_epil_hypo_vol <- function(H,A,td.depth,total_vol){
  vol_data <- matrix(NA, nrow = length(td.depth), ncol = 2)
  colnames(vol_data) <- c("vol_epil","vol_hypo")
  for (ii in 1:length(td.depth)){
    if(!is.na(td.depth[ii])){
      h_idx <- min(which(td.depth[ii]>=H))
      approx_td.area<-approx(H, A, c(0, td.depth[ii]))$y[2]
      if (is.na(approx_td.area)){approx_td.area = min(A)}
      H_with_td<-c(td.depth[ii],H[h_idx:length(H)])
      A_with_td<-c(approx_td.area,A[h_idx:length(A)])
      vol_data[ii,1] <- trapz(rev(H_with_td),rev(A_with_td))##epil
      vol_data[ii,2] <- total_vol - vol_data[ii,1]##hypo
      if(vol_data[ii,2] <= 0) {
        vol_data[ii,2]= min(A) * 0.5
        vol_data[ii,1] = total_vol - vol_data[ii,2]
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
calc_do<-function(input.values,fsed_stratified,fsed_not_stratified,nep_stratified,nep_not_stratified,
                  min_stratified, min_not_stratified, wind = NULL){
  ##initialize the o2(hypo),o2(epil),o2(total)
  o2_data <- matrix(NA, nrow = length(input.values$td.depth), ncol = 3)
  colnames(o2_data) <- c("o2_epil","o2_hypo","o2_total")
  o2_data <- as.data.frame(o2_data)
  
  init_o2sat <- o2.at.sat.base(temp=input.values$t.total[1],altitude = 300)*1000 # returns mg O2/L 
  o2_data$o2_total[1] <- init_o2sat*input.values$total_vol[1] # returns mg O2 (m3 = 1000 L)
  
  theta<-1.08
  for(day in 2:length(input.values$td.depth)){
    
    if (is.null(wind)){
      K600 <- k.cole.base(2) # returns m/day, assuming low wind conditions --> change to dynamic
    } else {
      K600 <- k.cole.base(wind[day])
    }
    
    ## not stratified period, only consider the o2(total)
    if(is.na(input.values$td.depth[day])){
      theta_total <- theta^(input.values$t.total[day]-20)
      
      NEP <- valid(nep_not_stratified * input.values$total_vol[day] * theta_total,
                   o2_data[day-1,"o2_total"])
      
      MINER <- valid(min_not_stratified * input.values$total_vol[day] * theta_total,
                   o2_data[day-1,"o2_total"])
      
      kO2 <- k600.2.kGAS.base(k600=K600,temperature=input.values$t.total[day],gas='O2') # velocity value m/d?
      o2sat<-o2.at.sat.base(temp=input.values$t.total[day],altitude = 300)*1000 # mg O2/L 
      Fatm <- valid(kO2*(o2sat*input.values$total_vol[day] - o2_data[day-1,"o2_total"])/max(H),
                    o2_data[day-1,"o2_total"])# mg/m * m/d = mg/d
  
      Fsed <- valid(fsed_not_stratified * o2_data[day-1,"o2_total"]/max(H) * theta_total,
                    o2_data[day-1,"o2_total"])# mg/m * m/d = mg/d
      
      o2_data[day,"o2_total"] <- o2_data[day-1,"o2_total"] - Fsed + NEP + Fatm + MINER# units make sense bc every term is actually multiplied
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
      NEP_epil <- valid(nep_stratified * input.values$vol_epil[day] * theta_epil,
                        o2_data[day-1,"o2_epil"])# has to return mg/d, mg/m3/d * m3
              
      
      kO2_epil <- k600.2.kGAS.base(k600=K600,temperature=input.values$t.epil[day],gas='O2')
      o2sat_epil<-o2.at.sat.base(temp=input.values$t.epil[day],altitude = 300)*1000
      Fatm_epil <- valid(kO2_epil*(o2sat_epil*input.values$vol_epil[day]-o2_data[day-1,"o2_epil"] )/input.values$td.depth[day],
                         o2_data[day-1,"o2_epil"])#possible wrong z
      
                     
      
      volumechange_epi = input.values$vol_epil[day]-input.values$vol_epil[day-1]  #in m^3
      volumechange_epi_proportion =  volumechange_epi/input.values$vol_epil[day-1] 
      Fepi <-  valid(volumechange_epi_proportion* o2_data[day-1,"o2_epil"],
                     o2_data[day-1,"o2_epil"])
      
      o2_data[day,"o2_epil"] <- o2_data[day-1,"o2_epil"]+NEP_epil+Fatm_epil + Fepi
      # print((o2_data[day,"o2_epil"]/input.values$vol_epil[day])/1000)
      
      ##hypo = hypo_o2[i-1]+Fhypo[i] - Fsed[i] 
      volumechange_hypo = input.values$vol_hypo[day]-input.values$vol_hypo[day-1]  #in m^3
      volumechange_hypo_proportion =  volumechange_hypo/input.values$vol_hypo[day-1] 
      Fhypo <- valid(volumechange_hypo_proportion* o2_data[day-1,"o2_hypo"],
                     o2_data[day-1,"o2_hypo"])
      
      MINER_hypo <- valid(min_stratified * input.values$vol_hypo[day] * theta_hypo,
                        o2_data[day-1,"o2_hypo"])# has to return mg/d, mg/m3/d * m3
      
      Fsed <- valid(fsed_stratified *o2_data[day-1,"o2_hypo"]/(max(H) - input.values$td.depth[day] ) * theta_hypo,
                    o2_data[day-1,"o2_hypo"])
      
      o2_data[day,"o2_hypo"] <- o2_data[day-1,"o2_hypo"] + Fhypo - Fsed + MINER_hypo#+ NEP_hypo +Fatm_hypo
# <<<<<<< HEAD
#       if(o2_data[day,"o2_hypo"]  < 0){
#         o2_data[day,"o2_hypo"]  = 0
#       }
      # print((o2_data[day,"o2_hypo"]/input.values$vol_hypo[day])/1000)
      ## total = hypo+epil

      #print((o2_data[day,"o2_hypo"]/input.values$vol_hypo[day])/1000)
     
       ## total = hypo+epil

      o2_data[day,"o2_total"] <- o2_data[day,"o2_hypo"] + o2_data[day,"o2_epil"]
      # print((o2_data[day,"o2_total"]/input.values$total_vol[day])/1000)
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
  
  for(jj in 1:length(ndate)){
    test_data$day[jj]<-ndate[jj]
    ##identify the index of the day from input values tested on the output
    kk = which(year(test_data$day[jj]) == year(input.values$datetime) & 
                 yday(test_data$day[jj])== yday(input.values$datetime))
    if(length(kk)!=1){
      next
    }
    if(is.na(input.values$td.depth[kk])){
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
optim_do <- function(p, input.values, fsed_not_stratified = 0.0002, nep_not_stratified = 0.0, min_not_stratified = 0.0,
                     wind = NULL, 
                     verbose){

  o2<- calc_do(input.values = input.values,fsed_stratified = p[1],
               fsed_not_stratified,
               nep_stratified = p[2],
               nep_not_stratified,
               min_stratified = p[3],
               min_not_stratified, wind)
  
  input.values$o2_epil <- o2[,"o2_epil"]
  input.values$o2_hypo <- o2[,"o2_hypo"]
  input.values$o2_total <- o2[,"o2_total"]
  
  input.values$year <- year(input.values$datetime)
  input.values$doy <- yday(input.values$datetime)
  
  
  test_data<-compare_predict_versus_observed(obs,input.values) 
  
  fit = calc_rmse(test_data)
  
  print(paste(round(p[1],2),round(p[2],2),round(p[3],2),'with RMSE: ',round(fit,3)))
  
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
      return (pool)
    }
    
  }else{
    return (flux)
  }
}





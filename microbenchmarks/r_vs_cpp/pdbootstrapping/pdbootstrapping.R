# original source: https://analyticsrusers.blog/2017/08/08/turbo-charge-your-r-code-with-rcpp/

source("data.R")

lagpad <- function(x, k) {
  if (!is.vector(x))
    stop('x must be a vector')
  if (!is.numeric(x))
    stop('x must be numeric')
  if (!is.numeric(k))
    stop('k must be numeric')
  if (1 != length(k))
    stop('k must be a single number')
  c(rep(NA, k), x)[1 : length(x)]
}

PDSensitivity_r <- function(list_dat1, LMResid, coef1, coef2, coef3, simNum) {
  # initialize variables
  avg_PD <- list()
  modeledTimeSeries <- list() # n is an array of 40 integers
  incrmnt <- 0
  
  for (s in 1:simNum){
    #initialize variable for each simulation
    ar1 <- 0
    ar2 <- 0
    sumPD <- 0
    avgPD <- 0
    getResid <- 0
    getRand  <- 0
    
    #first select the id
    id1 <- sample(1:nrow(dat1), 1, replace = TRUE)
    id2 <- id1 + 1
    incrmnt <- incrmnt + 1
    
    for (l in 1:nrow(dat1)) {
      if(l == 1) {
        modeledTimeSeries[l] <- as.numeric(list_dat1[id1])
      } else if (l == 2) {
        modeledTimeSeries[l] <- as.numeric(list_dat1[id2])
      } else {
        ar1  <- as.numeric(coef2) * as.numeric(modeledTimeSeries[l - 1]) #use the first coeff with the 1 periods value
        ar2  <- as.numeric(coef3) * as.numeric(modeledTimeSeries[l - 2]) #use the first coeff with the 2 periods value
        getResid  <-  LMResid[sample(1:nrow(dat1), 1, replace = TRUE)]
        modeledTimeSeries[l]  <- as.numeric(ar1) + as.numeric(ar2) + getResid + as.numeric(coef1)
      }
    }
    avg_PD[s]  <- mean(as.numeric(modeledTimeSeries))
  }
  return (avg_PD)
}



## Initial variable
avgPD_ALL_r <- list()
avgPD_ALL_cpp <- list()

## get dataset column names
gNames <- names(dat)

getData <- function(dat, x) {
  ## Create the dataset for the auto-regressive model of 2 lags
  dat1 <- dat[gNames[x]]
  dat1 <- as.data.frame(rep(dat1, 3))
  dat1 <- data.frame(x = as.numeric(unlist(dat1[1])), 
                     y = lagpad(as.numeric(unlist(dat1[2])), 1), 
                     z = lagpad(as.numeric(unlist(dat1[3])), 2))
  dat1 <- tail(dat1, -2) # remove first 2 record
  dat1
}

GNUR <- is.null(R.Version()$engine)

if(GNUR) {
  library(Renjin)
  
  t1 <- system.time( 
    for (i in 1:length(gNames)) {
      simNum <- 1000
      dat1 <- getData(dat, i)
      
      ## Run the model and get residuals
      LMfit          <- lm(x ~ y + z, data = dat1)
      LMResid        <- as.numeric(LMfit$residuals)
      
      ## Run bootstrap procedure
      avgPD <- PDSensitivity_r(
        as.numeric(unlist(dat1[1])), 
        as.numeric(LMfit$residuals), 
        LMfit$coef[1], 
        LMfit$coef[2], 
        LMfit$coef[3], simNum)
      
      ## store the 50th percentile of the average PDs
      avgPD_ALL_r[i] <- quantile(unlist(avgPD), c(.50), na.rm = TRUE)
    }
  print(avgPD_ALL_r)
  )
  
  t2 <- system.time( 
    for (i in 1:length(gNames)) {
      simNum <- 1000
      dat1 <- getData(dat, i)
      
      ## Run the model and get residuals
      LMfit          <- lm(x ~ y + z, data = dat1)
      LMResid        <- as.numeric(LMfit$residuals)
      
      ## Run bootstrap procedure
      avgPD <- PDSensitivity_cpp(
        as.numeric(unlist(dat1[1])), 
        as.numeric(LMfit$residuals), 
        LMfit$coef[1], 
        LMfit$coef[2], 
        LMfit$coef[3], simNum)
      
      ## store the 50th percentile of the average PDs
      avgPD_ALL_cpp[i] <- quantile(unlist(avgPD), c(.50), na.rm = TRUE)
    }
  print(avgPD_ALL_cpp)
  )
  
  t3 <- system.time(renjin( 
    avgPD_ALL_rj <- list()
    for (i in 1:length(gNames)) {
      simNum <- 1000
      dat1 <- getData(dat, i)
      
      ## Run the model and get residuals
      LMfit          <- lm(x ~ y + z, data = dat1)
      LMResid        <- as.numeric(LMfit$residuals)
      
      ## Run bootstrap procedure
      avgPD <- PDSensitivity_r(
        as.numeric(unlist(dat1[1])), 
        as.numeric(LMfit$residuals), 
        LMfit$coef[1], 
        LMfit$coef[2], 
        LMfit$coef[3], simNum)
      
      ## store the 50th percentile of the average PDs
      avgPD_ALL_rj[i] <- quantile(unlist(avgPD), c(.50), na.rm = TRUE)
    }
    print(avgPD_ALL_rj)
  ))
  
  timings <- rbind(t1, t2, t3)
  print(timings)
  
} else {
  
  t1 <- system.time( 
    for (i in 1:length(gNames)) {
      simNum <- 1000
      dat1 <- getData(dat, i)
      
      ## Run the model and get residuals
      LMfit          <- lm(x ~ y + z, data = dat1)
      LMResid        <- as.numeric(LMfit$residuals)
      
      ## Run bootstrap procedure
      avgPD <- PDSensitivity_r(
        as.numeric(unlist(dat1[1])), 
        as.numeric(LMfit$residuals), 
        LMfit$coef[1], 
        LMfit$coef[2], 
        LMfit$coef[3], simNum)
      
      ## store the 50th percentile of the average PDs
      avgPD_ALL_r[i] <- quantile(unlist(avgPD), c(.50), na.rm = TRUE)
    }
    print(avgPD_ALL_r)
  )
  
  t2 <- system.time( 
    for (i in 1:length(gNames)) {
      simNum <- 1000
      dat1 <- getData(dat, i)
      
      ## Run the model and get residuals
      LMfit          <- lm(x ~ y + z, data = dat1)
      LMResid        <- as.numeric(LMfit$residuals)
      
      ## Run bootstrap procedure
      avgPD <- PDSensitivity_cpp(
        as.numeric(unlist(dat1[1])), 
        as.numeric(LMfit$residuals), 
        LMfit$coef[1], 
        LMfit$coef[2], 
        LMfit$coef[3], simNum)
      
      ## store the 50th percentile of the average PDs
      avgPD_ALL_cpp[i] <- quantile(unlist(avgPD), c(.50), na.rm = TRUE)
    }
    print(avgPD_ALL_cpp)
  )
  
  timings <- rbind(t1, t2)
  print(timings)
}


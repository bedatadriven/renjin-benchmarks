
riginal source: https://analyticsrusers.blog/2017/08/08/turbo-charge-your-r-code-with-rcpp/

set.seed(1) #set random number seed
library(Rcpp)
## Create sample random data set to use in code
sampNum <- 40          # sample number
dat    <-  data.frame(
  cbind(
    runif(sampNum, min = 0.00235, max = 0.0053),
    runif(sampNum, min = 0.00902, max = 0.0242),
    runif(sampNum, min = 0.03046, max = 0.0660),
    runif(sampNum, min = 0.05843, max = 0.1168),
    runif(sampNum, min = 0.10427, max = 0.2055),
    runif(sampNum, min = 0.13828, max = 0.2600)
  )
)

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

PDSensitivity_cpp <- cppFunction("
                                 //[[Rcpp::export]]
                                 NumericVector PDSensitivity_cpp(NumericVector list_dat1, 
                                 NumericVector LMResid, int numRows, 
                                 double coef1, double coef2, 
                                 double coef3, double simNum) {
                                 
                                 NumericVector avg_PD(simNum); double sumPD = 0.0; double avgPD = 0.0; 
                                 int getRand; int id1; int id2; double modeledTimeSeries[numRows]; 
                                 double ar1; double ar2; double getResid; long incrmnt = 0;
                                 
                                 for (int s = 0.0; s < simNum ; s++) {
                                   getRand = 0; srand(42); id1 = rand() %numRows + 1; id2 = id1 + 1; 
                                   incrmnt = incrmnt + 1; ar1 = 0; ar2 = 0; getResid = 0; sumPD = 0.0; 
                                   avgPD = 0.0;
                                   
                                   for (int l = 0; l < numRows; l++) {
                                     if(l == 1) {
                                        modeledTimeSeries[l] = list_dat1[id1];
                                     } else if (l == 2) {
                                        modeledTimeSeries[l] = list_dat1[id2];
                                     } else {
                                       ar1 = coef2 * modeledTimeSeries[l - 1];
                                       ar2 = coef3 * modeledTimeSeries[l - 2];
                                       srand(42);
                                       getRand = rand() %numRows + 1;
                                       getResid = LMResid[getRand];
                                       modeledTimeSeries[l] = ar1 + ar2 + getResid + coef1;
                                     }
                                   }
                                   
                                   for (int k = 0; k < numRows; k++) {
                                    sumPD = sumPD + modeledTimeSeries[k];
                                   }
                                   
                                   avgPD = sumPD / numRows;
                                   avg_PD[s] = avgPD;
                                 }
                                 return avg_PD;
                                 }
                                 ")


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


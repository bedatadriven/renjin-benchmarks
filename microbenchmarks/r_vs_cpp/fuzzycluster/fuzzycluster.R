# original source: http://blog.revolutionanalytics.com/2017/08/kmeans-r-rcpp.html

source("data.R")

fuzzyClustering_r = function(data, centers, m) {
  c <- nrow(centers)
  rows <- nrow(data)
  cols <- ncol(data)
  result <- matrix(0, nrow = rows,ncol=c)  #defining membership matrix
  denominator <- 0
  
  for(i in 1:rows){
    for(j in 1:c){
      tempDist <- sqrt(sum((centers[j,]-data[i,])^2))  #euclidean distance, nominator inside a sum operator 
      for(k in 1:c){
        Dist <- sqrt(sum((centers[k,]-data[i,])^2))    #euclidean distance, denominator inside a sum operator 
        denominator <- denominator +((tempDist/Dist)^(2/(m-1))) #denominator of an equation 
      }    
      result[i,j] <- 1/denominator    #inserting value into membership matrix
      denominator <- 0
    }
  }
  return(result);
}


run <- function() {
  fuzzyClustering_r(mat, cxy, m)
}


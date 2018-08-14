# original source: http://blog.revolutionanalytics.com/2017/08/kmeans-r-rcpp.html

library(Rcpp)

GNUR <- is.null(R.Version()$engine)

if(GNUR) {
  library(Renjin)
  
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
  
  fuzzyClustering_cpp <- cppFunction("
  #include <Rcpp.h>
  #include <math.h>
  using namespace Rcpp;
  // [[Rcpp::export]]
  NumericMatrix fuzzyClustering_cpp(NumericMatrix data, NumericMatrix centers, int m) {  
    /* data is a matrix with observations(rows) and variables, 
    centers is a matrix with cluster centers coordinates, 
    m is a parameter of equation, c is a number of clusters
    */
      int c = centers.rows();
      int rows = data.rows();
      int cols = data.cols(); /* number of columns equals number of variables, the same as is in centers matrix */
      double tempDist = 0;        /* dist and tempDist are variables storing temporary euclidean distances */
      double dist = 0;      
      double denominator = 0;    // denominator of “main” equation
      NumericMatrix result(rows, c);    // declaration of matrix of results
      
      for(int i = 0; i < rows; i++) {
        for(int j = 0; j < c; j++) {
          for(int k = 0; k < c ; k++) {
            for(int p = 0; p < cols; p++) {
              tempDist = tempDist + pow(centers(j, p) - data(i, p), 2);      
              //in innermost loop an euclidean distance is calculated.
              dist = dist + pow(centers(k, p) - data(i, p), 2);              
              /* tempDist is nominator inside the sum operator in the equation, dist is the denominator inside the sum operator in the equation */
            }
            tempDist = sqrt(tempDist);
            dist = sqrt(dist);
            denominator = denominator + pow((tempDist / dist), (2 / (m - 1)));
            tempDist = 0;
            dist = 0;
          }
          result(i, j) = 1 / denominator;    
          // nominator / denominator in the  main equation
          denominator = 0;                           
        }
      }
      return result;
  }
  ")

  
  mat <- matrix(rnorm(300000), nrow = 1e5, ncol = 3)
  cxy <- matrix(c(1:10, rnorm(20)), nrow = 10, ncol = 3)
  m   <- 2
  
  
  t1 <- system.time(fuzzyClustering_r(mat, cxy, m))
  t2 <- system.time(fuzzyClustering_cpp(mat, cxy, m))
  t3 <- system.time(renjin(fuzzyClustering_r(mat, cxy, m)))
  timings <- rbind(t1, t2, t3)
  print(timings)
} else {
  
  fuzzyClustering_r = function(data, centers, m) {
    c <- nrow(centers)
    rows <- nrow(data)
    cols <- ncol(data)
    result <- matrix(0, nrow = rows, ncol = c)  #defining membership matrix
    denominator <- 0
    
    for(i in 1:rows){
      for(j in 1:c){
        tempDist <- sqrt(sum((centers[j, ] - data[i, ])^2))  #euclidean distance, nominator inside a sum operator 
        for(k in 1:c){
          Dist <- sqrt(sum((centers[k, ] - data[i, ])^2))    #euclidean distance, denominator inside a sum operator 
          denominator <- denominator + ((tempDist / Dist)^(2 / (m - 1))) #denominator of an equation 
        }    
        result[i, j] <- 1/denominator    #inserting value into membership matrix
        denominator <- 0           
      }
    }
    return(result);
  }
  
  fuzzyClustering_cpp <- cppFunction("
                                     #include <Rcpp.h>
                                     #include <math.h>
                                     using namespace Rcpp;
                                     // [[Rcpp::export]]
                                     NumericMatrix fuzzyClustering_cpp(NumericMatrix data, NumericMatrix centers, int m) {  
                                     /* data is a matrix with observations(rows) and variables, 
                                     centers is a matrix with cluster centers coordinates, 
                                     m is a parameter of equation, c is a number of clusters
                                     */
                                     int c = centers.rows();
                                     int rows = data.rows();
                                     int cols = data.cols(); /* number of columns equals number of variables, the same as is in centers matrix */
                                     double tempDist = 0;        /* dist and tempDist are variables storing temporary euclidean distances */
                                     double dist = 0;      
                                     double denominator = 0;    // denominator of “main” equation
                                     NumericMatrix result(rows, c);    // declaration of matrix of results
                                     
                                     for(int i = 0; i < rows; i++) {
                                     for(int j = 0; j < c; j++) {
                                     for(int k = 0; k < c ; k++) {
                                     for(int p = 0; p < cols; p++) {
                                     tempDist = tempDist + pow(centers(j, p) - data(i, p), 2);      
                                     //in innermost loop an euclidean distance is calculated.
                                     dist = dist + pow(centers(k, p) - data(i, p), 2);              
                                     /* tempDist is nominator inside the sum operator in the equation, dist is the denominator inside the sum operator in the equation */
                                     }
                                     tempDist = sqrt(tempDist);
                                     dist = sqrt(dist);
                                     denominator = denominator + pow((tempDist / dist), (2 / (m - 1)));
                                     tempDist = 0;
                                     dist = 0;
                                     }
                                     result(i, j) = 1 / denominator;    
                                     // nominator / denominator in the  main equation
                                     denominator = 0;                           
                                     }
                                     }
                                     return result;
                                     }
                                     ")

  
  mat <- matrix(rnorm(300000), nrow = 1e5, ncol = 3)
  cxy <- matrix(c(1:10, rnorm(20)), nrow = 10, ncol = 3)
  m   <- 2
  
  
  t1 <- system.time(fuzzyClustering_r(mat, cxy, m))
  t2 <- system.time(fuzzyClustering_cpp(mat, cxy, m))
  timings <- rbind(t1, t2)
  print(timings)
}


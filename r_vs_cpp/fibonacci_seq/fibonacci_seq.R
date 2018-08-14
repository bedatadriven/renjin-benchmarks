# original source: http://gallery.rcpp.org/articles/fibonacci-sequence/
library(Rcpp)

fibR <- function(n) {
  if ((n == 0) | (n == 1)) 
    return(1)
  else
    return(fibR(n - 1) + fibR(n - 2))
}

fibCpp <- cppFunction("
  #include <Rcpp.h>
  // [[Rcpp::export]]
  int fibCpp(int n) {
    if ((n == 0) | (n == 1)) 
      return 1;
    else
      return fibCpp(n - 1) + fibCpp(n - 2);
  }
")

GNUR <- is.null(R.Version()$engine)

if(GNUR) {
  library(Renjin)
  
  t1 <- system.time(fibR(20))
  t2 <- system.time(fibCpp(20))
  t3 <- system.time(renjin(fibR(20)))
  
  timings <- rbind(t1, t2, t3)
  print(timings)
} else {
  t1 <- system.time(fibR(20))
  t2 <- system.time(fibCpp(20))
  
  timings <- rbind(t1, t2)
  print(timings)
}

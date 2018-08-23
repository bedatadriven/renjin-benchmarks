# original source: http://gallery.rcpp.org/articles/fibonacci-sequence/
library(Rcpp)

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

run <- function(n = 20) {
  fibCpp(n)
}


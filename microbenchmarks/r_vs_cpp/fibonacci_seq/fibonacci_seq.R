# original source: http://gallery.rcpp.org/articles/fibonacci-sequence/
library(Rcpp)

fibR <- function(n) {
  if ((n == 0) | (n == 1)) 
    return(1)
  else
    return(fibR(n - 1) + fibR(n - 2))
}

run <- function(n = 20) {
  fibR(n)
}

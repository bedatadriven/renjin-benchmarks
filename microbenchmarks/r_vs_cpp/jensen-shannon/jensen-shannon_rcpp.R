# original source: http://gallery.rcpp.org/articles/parallel-distance-matrix/

library(Rcpp)
library(RcppArmadillo)

sourceCpp("jensen-shannon.cpp")

run <- function() {
  rcpp_js_distance(m)
  #rcpp_parallel_js_distance(m)
}


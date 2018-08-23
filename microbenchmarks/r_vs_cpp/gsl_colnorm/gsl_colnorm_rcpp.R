# original source: http://gallery.rcpp.org/articles/gsl-colnorm-example/
# Using GSL functions from R

library(Rcpp)
library(RcppGSL)

sourceCpp("functions.cpp")

M <- outer(sin(0:9), rep(1, 10), "*") + outer(rep(1, 10), cos(0:9), "*")

run <- function() {
  colNorm(M)
}


# original source: http://gallery.rcpp.org/articles/gsl-colnorm-example/
# Using GSL functions from R

library(Rcpp)
library(RcppGSL)

colNormR <- function(x) apply(x, 2, function(y) sqrt(sum(y^2)))

sourceCpp("functions.cpp")

GNUR <- is.null(R.Version()$engine)

if(GNUR) {
  library(Renjin)
  
  M <- outer(sin(0:9), rep(1, 10), "*") + outer(rep(1, 10), cos(0:9), "*")
  
  t1 <- system.time(colNorm(M))
  t2 <- system.time(colNormR(M))
  t3 <- system.time(renjin(colNormR(M)))
  timings <- rbind(t1, t2, t3)
  print(timings)
} else {
  M <- outer(sin(0:9), rep(1, 10), "*") + outer(rep(1, 10), cos(0:9), "*")
  
  t1 <- system.time(colNorm(M))
  t2 <- system.time(colNormR(M))
  timings <- rbind(t1, t2)
  print(timings)
}


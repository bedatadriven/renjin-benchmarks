# original source: http://gallery.rcpp.org/articles/parallel-matrix-transform/
# Transforming a Matrix in Parallel using RcppParallel

library(Rcpp)

m <- matrix(as.numeric(c(1:100000000)), nrow = 10000, ncol = 10000)

sourceCpp("functions.cpp")

GNUR <- is.null(R.Version()$engine)

if(GNUR) {
  library(Renjin)
  
  stopifnot(identical(matrixSqrt(m[1:5, 1:5]), parallelMatrixSqrt(m[1:5, 1:5])))
  stopifnot(identical(sqrt(m[1:5, 1:5]), parallelMatrixSqrt(m[1:5, 1:5])))
  
  t1 <- system.time(sqrt(m))
  t2 <- system.time(matrixSqrt(m))
  t3 <- system.time(parallelMatrixSqrt(m))
  t4 <- system.time(renjin(sqrt(m)))
  
  timings <- rbind(t1, t2, t3, t4)
  print(timings)
} else {
  
  stopifnot(identical(matrixSqrt(m[1:5, 1:5]), parallelMatrixSqrt(m[1:5, 1:5])))
  stopifnot(identical(sqrt(m[1:5, 1:5]), parallelMatrixSqrt(m[1:5, 1:5])))
  
  t1 <- system.time(sqrt(m))
  t2 <- system.time(matrixSqrt(m))
  t3 <- system.time(parallelMatrixSqrt(m))
  
  timings <- rbind(t1, t2, t3)
  print(timings)
}


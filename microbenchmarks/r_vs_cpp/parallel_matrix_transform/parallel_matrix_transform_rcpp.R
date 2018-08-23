# original source: http://gallery.rcpp.org/articles/parallel-matrix-transform/
# Transforming a Matrix in Parallel using RcppParallel

library(Rcpp)

source("data.R")
sourceCpp("parallel_matrix_transform.cpp")

run <- function() {
  parallelMatrixSqrt(m)
}


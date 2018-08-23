# original source: http://blog.revolutionanalytics.com/2017/08/kmeans-r-rcpp.html

library(Rcpp)

source("data.R")
sourceCpp("fuzzycluster.cpp")

run <- function() {
  fuzzyClustering_cpp(mat, cxy, m)
}
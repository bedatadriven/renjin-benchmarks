# original source: http://gallery.rcpp.org/articles/vector-minimum/
# Finding the minimum of a vector
library(Rcpp)

sourceCpp("functions.cpp")

vecminInd_R(x) which(x == min(x))[1]

set.seed(5)
x <- sample(1:1000, 1e9, replace = TRUE)  # ten out 100

GNUR <- is.null(R.Version()$engine)

if(GNUR) {
  library(Renjin)
  
  t1 <- system.time(vecmin(x))
  t2 <- system.time(min(x))
  t3 <- system.time(renjin(min(x)))
  t4 <- system.time(vecminInd(x) + 1)
  t5 <- system.time(vecminInd_R(x))
  t6 <- system.time(renjin(vecminInd_R(x)))
  
  timings <- rbind(t1, t2, t3, t4, t5, t6)
  print(timings)
} else {
  
  t1 <- system.time(vecmin(x))
  t2 <- system.time(min(x))
  t4 <- system.time(vecminInd(x) + 1)
  t5 <- system.time(vecminInd_R(x))
  
  timings <- rbind(t1, t2, t3, t4, t5, t6)
  print(timings)
}

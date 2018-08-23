# original source: http://gallery.rcpp.org/articles/parallel-distance-matrix/

library(Rcpp)
library(RcppArmadillo)

js_distance <- function(mat) {
  kld = function(p, q) sum(ifelse(p == 0 | q == 0, 0, log(p / q) * p))
  res = matrix(0, nrow(mat), nrow(mat))
  for (i in 1:(nrow(mat) - 1)) {
    for (j in (i + 1):nrow(mat)) {
      m = (mat[i, ] + mat[j, ]) / 2
      d1 = kld(mat[i, ], m)
      d2 = kld(mat[j, ], m)
      res[j, i] = sqrt(.5 * (d1 + d2))
    }
  }
  res
}

sourceCpp("functions.cpp")

# create a matrix
n  = 1000
m = matrix(runif(n * 10), ncol = 10)
m = m / rowSums(m)

# ensure that serial and parallel versions give the same result
r_res <- js_distance(m)
rcpp_res <- rcpp_js_distance(m)
rcpp_parallel_res <- rcpp_parallel_js_distance(m)
stopifnot(all(rcpp_res == rcpp_parallel_res))
stopifnot(all(rcpp_parallel_res - r_res < 1e-10)) ## precision differences

GNUR <- is.null(R.Version()$engine)

if(GNUR) {
  library(Renjin)

  t1 <- system.time(js_distance(m))
  t2 <- system.time(rcpp_js_distance(m))
  t3 <- system.time(rcpp_parallel_js_distance(m))
  t4 <- system.time(renjin(js_distance(m)))
  
  timings <- rbind(t1, t2, t3, t4)
  print(timings)
} else {
  
  t1 <- system.time(js_distance(m))
  t2 <- system.time(rcpp_js_distance(m))
  t3 <- system.time(rcpp_parallel_js_distance(m))
  
  timings <- rbind(t1, t2, t3)
  print(timings)
}


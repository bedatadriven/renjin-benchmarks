
# original source: http://gallery.rcpp.org/articles/dmvnorm-deriv-arma/

library('RcppArmadillo')
library('mvtnorm')
library('fields')

set.seed(123456789)
s <- rWishart(1, 2, diag(2))[ , , 1]
m <- rnorm(2)
X <- rmvnorm(10000, m, s)

dmvnorm_deriv1 <- function(X, mu = rep(0, ncol(X)), sigma = diag(ncol(X))) {
  fn <- function(x) -1 * c((1 / sqrt(det(2 * pi * sigma))) * exp(-0.5 * t(x-mu) %*% solve(sigma) %*% (x-mu))) * solve(sigma, (x-mu))
  out <- t(apply(X, 1, fn))
  return(out)
}

dmvnorm_deriv2 <- function(X, mean, sigma) {
  if (is.vector(X)) X <- matrix(X, ncol = length(X))
  if (missing(mean)) mean <- rep(0, length = ncol(X))
  if (missing(sigma)) sigma <- diag(ncol(X))
  
  n <- nrow(X)
  mvnorm <- dmvnorm(X, mean = mean, sigma = sigma)
  deriv <- array(NA, c(n, ncol(X)))
  
  for (i in 1:n) {
    deriv[i, ] <- -mvnorm[i] * solve(sigma, (X[i, ] - mean))
  }
  
  return(deriv)
}

sourceCpp("functions.cpp")

GNUR <- is.null(R.Version()$engine)

if(GNUR) {
  library(Renjin)
  
  t1 <- system.time(dmvnorm_deriv_arma(X, m, s))
  t2 <- system.time(dmvnorm_deriv1(X, mu = m, sigma = s))
  t3 <- system.time(dmvnorm_deriv2(X, mean = m, sigma = s))
  t4 <- system.time(renjin(dmvnorm_deriv1(X, mu = m, sigma = s)))
  t5 <- system.time(renjin(dmvnorm_deriv2(X, mean = m, sigma = s)))
  
  timings <- rbind(t1, t2, t3, t4, t5)
  print(timings)
} else {
  
  t1 <- system.time(dmvnorm_deriv_arma(X, m, s))
  t2 <- system.time(dmvnorm_deriv1(X, mu = m, sigma = s))
  t3 <- system.time(dmvnorm_deriv2(X, mean = m, sigma = s))
  
  timings <- rbind(t1, t2, t3)
  print(timings)
}


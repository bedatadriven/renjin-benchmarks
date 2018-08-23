

# original source: http://gallery.rcpp.org/articles/dmvnorm-deriv-arma/

library('RcppArmadillo')

load("simulated.RData")

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

run <- function() {
  dmvnorm_deriv2(X, mean = m, sigma = s)
}





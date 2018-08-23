

# original source: http://gallery.rcpp.org/articles/dmvnorm-deriv-arma/

load("simulated.RData")

dmvnorm_deriv <- function(X, mu = rep(0, ncol(X)), sigma = diag(ncol(X))) {
  fn <- function(x) -1 * c((1 / sqrt(det(2 * pi * sigma))) * exp(-0.5 * t(x-mu) %*% solve(sigma) %*% (x-mu))) * solve(sigma, (x-mu))
  out <- t(apply(X, 1, fn))
  return(out)
}

run <- function() {
  dmvnorm_deriv(X,mu=m,sigma=s)
}





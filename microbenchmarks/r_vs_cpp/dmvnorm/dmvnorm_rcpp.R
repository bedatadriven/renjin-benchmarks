# original source: http://gallery.rcpp.org/articles/dmvnorm-deriv-arma/

library('RcppArmadillo')

load("simulated.RData")

sourceCpp("functions.cpp")

run <- function() {
  dmvnorm_deriv_arma(X, m, s)
}
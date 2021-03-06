# original source: http://gallery.rcpp.org/articles/gsl-colnorm-example/
# Using GSL functions from R


M <- outer(sin(0:9), rep(1, 10), "*") + outer(rep(1, 10), cos(0:9), "*")

colNorm <- function(x) apply(x, 2, function(y) sqrt(sum(y^2)))

run <- function() {
  colNorm(M)
}
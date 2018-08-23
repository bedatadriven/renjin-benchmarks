# original source: http://rpubs.com/WolfgangHuber/404211
library(Rcpp)


compute_pi <- function(m) {
  s = 0
  sign = 1
  for (n in 0:m) {
    s = s + sign / (2 * n + 1)
    sign = -sign
  }
  4 * s
}

run <- function(n) {
  compute_pi(n)
}


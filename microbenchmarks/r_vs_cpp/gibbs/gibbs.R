# original source: http://adv-r.had.co.nz/Rcpp.html


run <- function(n = 1000, thin = 10) {
  mat <- matrix(nrow = n, ncol = 2)
  x <- y <- 0
  
  for (i in 1:n) {
    for (j in 1:thin) {
      x <- rgamma(1, 3, y * y + 4)
      y <- rnorm(1, 1 / (x + 1), 1 / sqrt(2 * (x + 1)))
    }
    mat[i, ] <- c(x, y)
  }
  mat
}

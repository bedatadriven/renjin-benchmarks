

a <- abs(matrix(rnorm(2500*2500)/2, ncol=2500, nrow=2500));

run <- function() {
  b <- a^1000 
  colMeans(b)
}

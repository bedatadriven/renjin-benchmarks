

a <- matrix(rnorm(2500*2500)/10, ncol=2500, nrow=2500);

run <- function() {
  b <- t(a);
  dim(b) <- c(1250, 5000);
  a <- t(b)
  colMeans(a)
}

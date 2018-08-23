
a <- rnorm(900*900)
dim(a) <- c(900, 900)

run <- function() {
  a <- crossprod(a, a)
  b <- chol(a)
  b
}

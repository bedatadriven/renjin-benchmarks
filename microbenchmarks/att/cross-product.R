
a <- rnorm(2800*2800)
dim(a) <- c(2800, 2800)

run <- function() {
  b <- crossprod(a)  	# equivalent to: b <- t(a) %*% a
  b
}
a <- rnorm(900*900)
dim(a) <- c(900, 900)
a <- crossprod(a, a)
b <- chol(a)

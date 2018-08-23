
a <- array(rnorm(600*600), dim = c(600, 600))

run <- function() {
  b <- eigen(a, symmetric=FALSE, only.values=TRUE)$Value
  b
}



a <- rnorm(2500*2500)
dim(a) <- c(2500, 2500)

run <- function() {
  b <- det(a)
  b
}
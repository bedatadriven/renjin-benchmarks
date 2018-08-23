
a <- 0
b <- 0
phi <- 1.6180339887498949
a <- floor(runif(3500000)*1000)

run <- function() {
  (phi^a - (-phi)^(-a))/sqrt(5)
}

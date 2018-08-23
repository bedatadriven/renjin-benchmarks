

a <- rnorm(2400000)

run <- function() {
  b <- fft(a)
  b  
}

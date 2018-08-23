# Function that creates a gaussian vector and compute its summary
load_test <- function() {
  ptm <- proc.time()["elapsed"]
  tata <- rnorm(n=10^7,mean = 0, sd = 1)
  s <- summary(tata)
  d <- proc.time()["elapsed"] - ptm
  return(d)
}

# Run the function 10 times and compute the mean value of execution
toto <- replicate(10,load_test())
print(mean(toto))

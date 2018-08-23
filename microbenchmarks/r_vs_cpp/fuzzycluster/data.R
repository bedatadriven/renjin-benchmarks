
# Generates test data 

set.seed(25352)
mat <- matrix(rnorm(300000), nrow = 1e5, ncol = 3)
cxy <- matrix(c(1:10, rnorm(20)), nrow = 10, ncol = 3)
m   <- 2

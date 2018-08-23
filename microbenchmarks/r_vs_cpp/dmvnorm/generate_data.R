
# This script generates simulated data for the 
# benchmark. It is not run during the actual benchmark execution

library(fields)
library(mvtnorm)

set.seed(123456789)
s <- rWishart(1, 2, diag(2))[,,1]
m <- rnorm(2)
X <- rmvnorm(10000, m, s)

save(s, m, X, file = "simulated.RData")
a <- abs(matrix(rnorm(2500*2500)/2, ncol=2500, nrow=2500));
b <- a^1000 

# Do something with result to ensure that
# Renjin doesn't optimize away the operation
print(sum(colMeans(b))


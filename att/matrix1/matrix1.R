a <- matrix(rnorm(2500*2500)/10, ncol=2500, nrow=2500);
b <- t(a);
dim(b) <- c(1250, 5000);
a <- t(b)

# Ensure that the result is actually used to ensure
# that Renjin doesn't optimize away the result
print(sum(colMeans(a)))



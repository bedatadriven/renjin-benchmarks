
# create a matrix
n  = 1000
m = matrix(runif(n * 10), ncol = 10)
m = m / rowSums(m)


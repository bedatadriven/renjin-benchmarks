a <- rnorm(7000000)
b <- sort(a, method="quick")

# Consume the output to ensure the operation is not 
# optimized away
stopifnot(which.min(b) == 1)
a <- 3000
b <- 0

b <- rep(1:a, a)
dim(b) <- c(a, a);

b <- 1 / (t(b) + 0:(a-1))

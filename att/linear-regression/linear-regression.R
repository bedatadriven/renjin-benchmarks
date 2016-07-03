
a <- rnorm(600*600)
dim(a) <- c(600,600)
b <- 1:600
qra <- qr(a, tol = 1e-7);
c <- qr.coef(qra, b)

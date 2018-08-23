a <- array(rnorm(600*600), dim = c(600, 600))
b <- eigen(a, symmetric=FALSE, only.values=TRUE)$Value

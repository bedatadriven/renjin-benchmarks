# original source: http://gallery.rcpp.org/articles/EM-algorithm-example/
# Implementing an EM Algorithm for Probit Regressions

library(Rcpp)
library(RcppArmadillo)

Sys.setenv("PKG_CXXFLAGS" = "-fopenmp")
Sys.setenv("PKG_LIBS" = "-fopenmp")

sourceCpp("em3.cpp")
sourceCpp("em4.cpp")

turnout <- read.table("turnout.csv", sep = ";", stringsAsFactors = TRUE, col.names = TRUE, row.names = FALSE)
mY <- matrix(turnout$vote)
mX <- cbind(1, turnout$income, turnout$educate, turnout$age)

GNUR <- is.null(R.Version()$engine)

if(GNUR) {
  library(Renjin)
  
  t1 <- system.time(fit0 <- glm(vote ~ income + educate + age, data = turnout, family = binomial(link = "probit")))
  t2 <- system.time(fit3 <- em3(y = mY, X = mX, maxit = 100))
  t3 <- system.time(fit4 <- em4(y = mY, X = mX, maxit = 100, nthr = 4))
  t4 <- system.time(renjin( glm(vote ~ income + educate + age, data = turnout, family = binomial(link = "probit")) ))
  
  stopifnot(identical(fit4$beta, fit3$beta))
  stopifnot(identical(fit0$beta, fit3$beta))
  
  timings <- rbind(t1, t2, t3, t4)
  print(timings)
} else {
  
  t1 <- system.time(fit0 <- glm(vote ~ income + educate + age, data = turnout, family = binomial(link = "probit")))
  t2 <- system.time(fit3 <- em3(y = mY, X = mX, maxit = 100))
  t3 <- system.time(fit4 <- em4(y = mY, X = mX, maxit = 100, nthr = 4))
  
  stopifnot(identical(fit4$beta, fit3$beta))
  stopifnot(identical(fit0$beta, fit3$beta))
  
  timings <- rbind(t1, t2, t3)
  print(timings)
}


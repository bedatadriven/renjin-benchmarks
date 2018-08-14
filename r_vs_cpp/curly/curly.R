# original source:http://dirk.eddelbuettel.com/blog/2010/09/07/
library(Rcpp)
library(inline)

  f <- function(n, x = 1) for (i in 1:n) x = 1 / (1 + x)
  g <- function(n, x = 1) for (i in 1:n) x = (1 / (1 + x))
  h <- function(n, x = 1) for (i in 1:n) x = (1 + x)^(-1)
  j <- function(n, x = 1) for (i in 1:n) x = {1 / {1 + x}}
  k <- function(n, x = 1) for (i in 1:n) x = 1 / {1 + x}
  l <- cxxfunction(signature(ns = "integer", xs = "numeric"),
                   'int n = as<int>(ns); double x=as<double>(xs);
                  for (int i = 0; i < n; i++) x = 1 / (1 + x);
                  return wrap(x); ',
                   plugin="Rcpp")
  N <- 1e6

GNUR <- is.null(R.Version()$engine)

if(GNUR) {
  library(Renjin)
  
  t1  <- system.time(f(N, 1))
  t2  <- system.time(g(N, 1))
  t3  <- system.time(h(N, 1))
  t4  <- system.time(j(N, 1))
  t5  <- system.time(k(N, 1))
  t6  <- system.time(l(N, 1))
  
  t1r <- system.time(renjin(f(N, 1)))
  t2r <- system.time(renjin(g(N, 1)))
  t3r <- system.time(renjin(h(N, 1)))
  t4r <- system.time(renjin(j(N, 1)))
  t5r <- system.time(renjin(k(N, 1)))
  t6r <- system.time(renjin(l(N, 1)))
  
  timings <- rbind(t1, t1r, t2, t2r, t3, t3r, t4, t4r, t5, t5r, t6, t6r)
  print(timings)
} else {
  
  t1  <- system.time(f(N, 1))
  t2  <- system.time(g(N, 1))
  t3  <- system.time(h(N, 1))
  t4  <- system.time(j(N, 1))
  t5  <- system.time(k(N, 1))
  t6  <- system.time(l(N, 1))
  
  timings <- rbind(t1, t2, t3, t4, t5, t6)
  print(timings)
}


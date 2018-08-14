# original source: http://rpubs.com/WolfgangHuber/404211
library(Rcpp)

  
  compute_pi <- function(m) {
    s = 0
    sign = 1
    for (n in 0:m) {
      s = s + sign / (2 * n + 1)
      sign = -sign
    }
    4 * s
  }
  
  compute_pi_vec_1 <- function(m) {
    n <- 0:m
    4 * sum((-1)^n / (2 * n + 1))
  }
  
  compute_pi_vec_2 <- function(m) {
    n <- seq(0, (m - 1) / 2, by = 2)
    den <- 2 * n + 1
    4 * (sum(1 / den) - sum(1 / (den + 2)))
  }
  
  compute_pi_cpp <- cppFunction("
    double compute_pi_cpp(int m) {
        double s = 0, sign = 1;
        for (int n = 0; n <= m; ++n) {
            s = s + sign / (2 * n + 1);
            sign = -sign;
        }
        return 4*s;
    }
  ")

  compute_pi_bcc <- compiler::cmpfun(compute_pi) 

GNUR <- is.null(R.Version()$engine)

if(GNUR) {
  library("Rcpp")
  
  pis <- rep(NA_real_, 6)
  m <- 5e8
  
  t1 = system.time(pis[1] <- compute_pi(m))
  t2 = system.time(pis[2] <- compute_pi_bcc(m))
  t3 = system.time(pis[3] <- renjin(compute_pi(m)))
  t4 = system.time(pis[4] <- compute_pi_vec_1(m))
  t5 = system.time(pis[5] <- compute_pi_vec_2(m))
  t6 = system.time(pis[6] <- compute_pi_cpp(m))
  
  timings <- rbind(t1, t2, t3, t4, t5, t6)
  
  print(pis)
  print(timings)
} else {
  
  pis <- rep(NA_real_, 6)
  m <- 5e8
  
  t1 = system.time(pis[1] <- compute_pi(m))
  t2 = system.time(pis[2] <- compute_pi_bcc(m)) 
  t3 = system.time(pis[4] <- compute_pi_vec_1(m))
  t4 = system.time(pis[5] <- compute_pi_vec_2(m))
  t5 = system.time(pis[6] <- compute_pi_cpp(m))
  
  timings <- rbind(t1, t2, t3, t4, t5)
  
  print(timings)
  print(pis)
}



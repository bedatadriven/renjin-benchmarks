# original source: http://rpubs.com/WolfgangHuber/404211

compute_pi_vec_1 <- function(m) {
  n <- 0:m
  4 * sum((-1)^n / (2 * n + 1))
}

compute_pi_vec_2 <- function(m) {
  n <- seq(0, (m - 1) / 2, by = 2)
  den <- 2 * n + 1
  4 * (sum(1 / den) - sum(1 / (den + 2)))
}
  
run <- function(n = 5e8) {
  compute_pi_vec_2(n)
}




# original source: http://rpubs.com/WolfgangHuber/404211
library(Rcpp)

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

run <- function(n = 5e8) {
  compute_pi_cpp(n)
}


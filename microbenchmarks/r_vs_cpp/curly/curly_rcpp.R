# original source:http://dirk.eddelbuettel.com/blog/2010/09/07/

library(Rcpp)

l <- cppFunction("
    double compute_curly(int n) {
        double x= 0;
        for (int i = 0; i < n; i++) x = 1 / (1 + x);
        return x;
    }
  ")

run <- function(n = 1e6) {
  l(n)
}
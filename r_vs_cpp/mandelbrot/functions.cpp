#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace RcppArmadillo;

int mandelbrot_cpp(double creal, double cimag, int maxiter) {
  double real = creal;
  double imag = cimag;
  for(int i = 0; i < maxiter; i++) {
    double real2 = real * real;
    double imag2 = imag * imag;
    if ((real2 + imag2) > 4) {
      return i + 1;
    }
    imag = 2 * real * imag + cimag;
    real = real2 - imag2 + creal;
  }
  return 0;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericMatrix mandelbrot_set_cpp(NumericVector r1, NumericVector r2, int maxiter) {
  arma::mat n3(r1.length(), r2.length());
  for (int i = 0; i < r1.length(); i++) {
    for (int j = 0; j < r2.length(); j++) {
      n3(i, j) = mandelbrot_cpp(r1[i], r2[j], maxiter);
    }
  }
  return wrap(n3);
}

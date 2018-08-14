
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec Mahalanobis(arma::mat x, arma::rowvec center, arma::mat cov){
  int n = x.n_rows;
  arma::mat x_cen;
  x_cen.copy_size(x);
  for (int i = 0; i < n; i++) {
    x_cen.row(i) = x.row(i) - center;
  }
  return sum((x_cen * cov.i()) % x_cen, 1);
}

// [[Rcpp::export]]
arma::mat dmvnorm_deriv_arma(arma::mat x, arma::rowvec mean, arma::mat sigma) {
  // get result for mv normal
  arma::vec distval = Mahalanobis(x,  mean, sigma);
  double logdet = sum(arma::log(arma::eig_sym(sigma)));
  double log2pi = std::log(2.0 * M_PI);
  arma::vec mvnorm = exp(-( (x.n_cols * log2pi + logdet + distval)/2));
  
  // get derivative of multivariate normal
  int n = x.n_rows;
  arma::mat deriv;
  deriv.copy_size(x);
  for (int i = 0; i < n; i++) {
    deriv.row(i) = -1 * mvnorm(i) * trans(solve(sigma, trans(x.row(i) - mean)));
  }
  return(deriv);
}

# include <RcppArmadillo.h>
# include <omp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

double f (double mu) {
  double val = ((R::dnorm(-mu, 0, 1, false)) /
                  (1 - R::pnorm(-mu, 0, 1, true, false))
  ) ;
  return(val) ;
}

double g (double mu) {
  double val = ((R::dnorm(-mu, 0, 1, false)) /
                  (R::pnorm(-mu, 0, 1, true, false))
  ) ;
  return(val) ;
}

// [[Rcpp::export()]]
List em4 (const arma::mat y,
          const arma::mat X,
          const int maxit = 10,
          const int nthr = 1
) {
  // inputs
  const int N = y.n_rows ;
  const int K = X.n_cols ;
  omp_set_num_threads(nthr) ;
  
  // containers
  arma::mat beta(K, 1) ;
  beta.fill(0.0) ; // initialize betas to 0
  arma::mat eystar(N, 1) ;
  eystar.fill(0) ;
  
  // algorithm
  for (int it = 0 ; it < maxit ; it++) {
    arma::mat mu = X * beta ;
    // augmentation step
    #pragma omp parallel for
    for (int n = 0 ; n < N ; n++) {
      if (y(n, 0) == 1) { // y = 1
      eystar(n, 0) = mu(n, 0) + f(mu(n, 0)) ;
      }
      if (y(n, 0) == 0) { // y = 0
      eystar(n, 0) = mu(n, 0) - g(mu(n, 0)) ;
      }
    }
    // maximization step
    beta = (X.t() * X).i() * X.t() * eystar ;
  }
  
  // returns
  List ret ;
  ret["N"] = N ;
  ret["K"] = K ;
  ret["beta"] = beta ;
  ret["eystar"] = eystar ;
  return(ret) ;
}


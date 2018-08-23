// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
// [[Rcpp::export]]
Rcpp::NumericVector colNorm(const RcppGSL::Matrix & M) {
  
  int k = M.ncol();
  Rcpp::NumericVector n(k); 		// to store results 
  
  for (int j = 0; j < k; j++) {
    RcppGSL::VectorView colview = gsl_matrix_const_column (M, j);
    n[j] = gsl_blas_dnrm2(colview);
  }
  
  return n;				// return vector  
}


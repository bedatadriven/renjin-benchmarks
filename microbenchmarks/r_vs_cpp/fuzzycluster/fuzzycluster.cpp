#include <Rcpp.h>
using namespace Rcpp;

#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix fuzzyClustering_cpp(NumericMatrix data, NumericMatrix centers, int m) {  
  /* data is a matrix with observations(rows) and variables, 
   centers is a matrix with cluster centers coordinates, 
   m is a parameter of equation, c is a number of clusters
   */
  int c = centers.rows();
  int rows = data.rows();
  int cols = data.cols(); /* number of columns equals number of variables, the same as is in centers matrix */
  double tempDist = 0;        /* dist and tempDist are variables storing temporary euclidean distances */
  double dist = 0;      
  double denominator = 0;    // denominator of “main” equation
  NumericMatrix result(rows, c);    // declaration of matrix of results
  
  for(int i = 0; i < rows; i++) {
    for(int j = 0; j < c; j++) {
      for(int k = 0; k < c ; k++) {
        for(int p = 0; p < cols; p++) {
          tempDist = tempDist + pow(centers(j, p) - data(i, p), 2);      
          //in innermost loop an euclidean distance is calculated.
          dist = dist + pow(centers(k, p) - data(i, p), 2);              
          /* tempDist is nominator inside the sum operator in the equation, dist is the denominator inside the sum operator in the equation */
        }
        tempDist = sqrt(tempDist);
        dist = sqrt(dist);
        denominator = denominator + pow((tempDist / dist), (2 / (m - 1)));
        tempDist = 0;
        dist = 0;
      }
      result(i, j) = 1 / denominator;    
      // nominator / denominator in the  main equation
      denominator = 0;                           
    }
  }
  return result;
}
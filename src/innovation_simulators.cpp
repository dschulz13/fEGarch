#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec rnorm_fEGarch(const int n) {

  arma::vec probs = Rcpp::runif(n, 0.0, 1.0);



  return probs;

}

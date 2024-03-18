#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector pdf_norm(NumericVector x) {

  return Rcpp::dnorm(x);

}

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Author: Dominik Schulz
// Version: February 20, 2025

// [[Rcpp::export]]
arma::vec autocorCpp(const arma::vec& xt, const int lag_max) {

  const double xmean = arma::mean(xt);
  const arma::vec x_dm = xt - xmean;

  const double var = arma::accu(arma::pow(x_dm, 2.0));

  arma::vec autocov = arma::zeros<arma::vec>(lag_max);
  const int n = x_dm.size();
  const int nm1 = n - 1;

  for (int i = 1; i <= lag_max; ++i) {
    autocov(i - 1) = arma::accu(x_dm.subvec(0, nm1 - i) % x_dm.subvec(i, nm1));
  }

  return autocov / var;

}

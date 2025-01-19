#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec sigt_lmgarch_Cpp(const arma::vec& x, const double mu,
                       const double omega,
                       const arma::rowvec& coef_inf,
                       const double E_e2) {

  const int n = x.size();
  const int l = n + 1000;
  arma::vec sig2 = arma::vec(n).fill(0);

  const arma::vec x_cent = x - mu;

  //Initialization based on sample
  arma::vec e2 = arma::vec(l).fill(E_e2);
  e2.subvec(1000, l - 1) = arma::pow(x_cent, 2.0);

  const arma::rowvec coef_r = arma::reverse(coef_inf);

  for (int i = 0; i < n; ++i) {

    sig2(i) = omega +
      arma::as_scalar(coef_r.subvec(l - 1000 - i, l - 1) * e2.subvec(0, 999 + i));

  }

  return arma::sqrt(sig2);

}

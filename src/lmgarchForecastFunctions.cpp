#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
arma::vec sigt_lmgarch_forecast(const arma::vec& rt,
                       const arma::vec& coef_inf,
                       const double omega,
                       const int horizon) {

  const int n = rt.size();
  arma::vec sig2 = arma::vec(horizon);
  arma::vec g_fun = arma::zeros<arma::vec>(n + horizon);
  g_fun.subvec(0, n - 1) = arma::pow(rt, 2.0);

  const arma::vec coef_inf_r = arma::reverse(coef_inf);

  const int nm1 = n - 1;
  const int nHm2 = nm1 + horizon - 1;
  const int np1 = - n + 1;
  const int z = nHm2 + np1;

  for (int i = 0; i < horizon; ++i) {

    sig2(i) = omega +
      arma::as_scalar(arma::dot(coef_inf_r.subvec(z - i, nHm2), g_fun.subvec(0, nm1 + i)));

    g_fun(n + i) = sig2(i);

  }

  return arma::sqrt(sig2);

}

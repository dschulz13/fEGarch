#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
arma::vec sigt_loggarch_forecast_long(const arma::vec& et,
                       const arma::rowvec& coef_inf,
                       const double Elneta2,
                       const double Elnsig2,
                       const int horizon) {

  const int n = et.size();
  arma::vec lnsig2 = arma::vec(horizon);
  arma::vec g_fun = arma::zeros<arma::vec>(n + horizon);
  g_fun.subvec(0, n - 1) = arma::log(arma::pow(et, 2.0)) - Elneta2;

  const arma::rowvec coef_inf_r = arma::reverse(coef_inf);

  int nm1 = n - 1;
  const int nHm2 = nm1 + horizon - 1;
  const int np1 = - n + 1;
  const int z = nHm2 + np1;

  for (int i = 0; i < horizon; ++i) {

    lnsig2(i) = Elnsig2 +
      arma::as_scalar(coef_inf_r.subvec(z - i, nHm2) * g_fun.subvec(0, nm1 + i));

  }

  return arma::exp(0.5 * lnsig2);

}

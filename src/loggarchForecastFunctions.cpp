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

// [[Rcpp::export]]
arma::vec sigt_loggarch_forecast_short(const arma::vec& et,
                        const arma::vec& sigt,
                       const double omega,
                       const arma::rowvec& phi,
                       const arma::rowvec& psi,
                       const double Elneta2,
                       const int horizon) {

  const int n = et.size();
  const int p = phi.size();
  const int q = psi.size();
  const int l = std::max(p, q);
  arma::vec lnsig2 = arma::zeros<arma::vec>(l + horizon);
  lnsig2.subvec(0, l - 1) = arma::log(arma::pow(sigt.subvec(n - l, n - 1), 2.0));
  arma::vec g_fun = arma::vec(l + horizon).fill(0);
  g_fun.subvec(0, l - 1) = arma::log(arma::pow(et.subvec(n - l, n - 1), 2.0));

  const arma::rowvec psi_r = arma::reverse(psi);
  const arma::rowvec phi_r = arma::reverse(phi);

  int il;
  int il1;
  int ilp;
  int ilq;

  for (int i = 0; i < horizon; ++i) {

    il = i + l;
    il1 = il - 1;
    ilp = il - p;
    ilq = il - q;
    lnsig2(l + i) = omega +
      arma::as_scalar(psi_r * g_fun.subvec(ilq, il1)) +
      arma::as_scalar(phi_r * lnsig2.subvec(ilp, il1));

  }

  return arma::exp(0.5 * lnsig2.subvec(l, l + horizon - 1));

}

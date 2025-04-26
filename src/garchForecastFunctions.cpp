#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec sigt_garch_forecast_short(
                               const arma::vec& et,
                               const arma::vec& sigt,
                               const double omega,
                               const arma::vec& phi,
                               const arma::vec& beta,
                               const int horizon) {

  const int n = sigt.size();
  const int p = phi.size();
  const int q = beta.size();
  const int l = std::max(p, q);
  const int Hl = horizon + l;

  int il, il1, ilp, ilq;

  arma::vec et2= arma::ones<arma::vec>(Hl);
  et2.subvec(0, l - 1) = arma::pow(et.subvec(n - l, n - 1), 2.0);

  arma::vec sig2(Hl);
  sig2.subvec(0, l - 1) = arma::pow(sigt.subvec(n - l, n - 1), 2.0);

  const arma::vec phi_r = arma::reverse(phi);
  const arma::vec beta_r = arma::reverse(beta);

  for (int i = 0; i < horizon; ++i) {
    il = i + l;
    il1 = il - 1;
    ilp = il - p;
    ilq = il - q;

    sig2(il) = omega +
      arma::as_scalar(arma::dot(phi_r, (sig2.subvec(ilp, il1) % et2.subvec(ilp, il1)))) +
      arma::as_scalar(arma::dot(beta_r, sig2.subvec(ilq, il1)));

  }

  return arma::sqrt(sig2.subvec(l, Hl - 1));
}

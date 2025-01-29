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

// [[Rcpp::export]]
Rcpp::List lmgarch_sim_Cpp(const arma::vec& innov, const double omega,
                          const arma::rowvec& coef_inf,
                          const int n_out,
                          const double mu) {

  const int n = innov.size();

  const int n_pre = 1;
  const int n_all = n + n_pre;

  double sig2;
  arma::vec sigt = arma::vec(n);

  const arma::rowvec coef_rev = arma::reverse(coef_inf);

  arma::vec rt = arma::zeros<arma::vec>(n_all);
  arma::vec eps2 = arma::zeros<arma::vec>(n_all);

  const int n_coef_m1 = coef_inf.size() - 1;
  int npp1;

  for (int i = 0; i < n; ++i) {

    npp1 = 1 + i;

    sig2 = omega +
      arma::as_scalar(coef_rev.subvec(n_coef_m1 - i, n_coef_m1) * eps2.subvec(0, i));

    sigt(i) = std::sqrt(sig2);

    rt(npp1) = innov(i) * sigt(i);
    eps2(npp1) = std::pow(rt(npp1), 2.0);

  }

  const arma::vec sigt_out = sigt.subvec(n - n_out, n - 1);
  const arma::vec innov_out = innov.subvec(n - n_out, n - 1);

  const arma::vec rt_out = rt.subvec(n_all - n_out, n_all - 1) + mu;

  return Rcpp::List::create(_["rt"] = rt_out, _["sigt"] = sigt_out, _["etat"] = innov_out);
}

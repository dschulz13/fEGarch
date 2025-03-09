#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
arma::vec sigt_fiaparch_Cpp(const arma::vec& x, const double mu,
                       const arma::vec& coef_inf,
                       const double gamma,
                       const double delta,
                       const double omega_sig,
                       const double etransf_init) {

  const int coef_len = coef_inf.size();
  const int n = x.size();
  const int l = n + coef_len;
  arma::vec sig_delta = arma::vec(n).fill(0);

  const arma::vec x_cent = x - mu;
  const arma::vec transf = arma::pow(arma::abs(x_cent) - gamma * x_cent, delta);

  arma::vec g_fun = arma::vec(l).fill(etransf_init);
  g_fun.subvec(coef_len, l - 1) = transf;

  const arma::vec coef_r = arma::reverse(coef_inf);

  const int cm1 = coef_len - 1;
  int cm1i;

  for (int i = 0; i < n; ++i) {

    cm1i = cm1 + i;

    sig_delta(i) = omega_sig +
      arma::as_scalar(arma::dot(coef_r, g_fun.subvec(i, cm1i)));

  }

  return arma::pow(sig_delta, 1.0 / delta);

}

// [[Rcpp::export]]
Rcpp::List fiaparch_sim_Cpp(const arma::vec& innov, const double omega,
                          const arma::vec& coef_inf,
                          const double gamma, const double delta,
                          const int n_out,
                          const double mu) {

  const int n = innov.size();

  const int n_pre = 1;
  const int n_all = n + n_pre;

  double sig_d;
  arma::vec sigt = arma::vec(n);

  const arma::vec coef_rev = arma::reverse(coef_inf);

  arma::vec rt = arma::zeros<arma::vec>(n_all);
  arma::vec g_fun = arma::zeros<arma::vec>(n_all);

  const int n_coef_m1 = coef_inf.size() - 1;
  int npp1;

  for (int i = 0; i < n; ++i) {

    npp1 = 1 + i;

    sig_d = omega +
      arma::as_scalar(arma::dot(coef_rev.subvec(n_coef_m1 - i, n_coef_m1), g_fun.subvec(0, i)));

    sigt(i) = std::pow(sig_d, 1.0 / delta);

    rt(npp1) = innov(i) * sigt(i);
    g_fun(npp1) = std::pow(std::abs(rt(npp1)) - gamma * rt(npp1), delta);

  }

  const arma::vec sigt_out = sigt.subvec(n - n_out, n - 1);
  const arma::vec innov_out = innov.subvec(n - n_out, n - 1);

  const arma::vec rt_out = rt.subvec(n_all - n_out, n_all - 1) + mu;

  return Rcpp::List::create(_["rt"] = rt_out, _["sigt"] = sigt_out, _["etat"] = innov_out);
}

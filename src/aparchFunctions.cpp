#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
arma::vec sigt_aparch_Cpp(const arma::vec& x, const double mu,
                       const arma::rowvec& phi, const arma::rowvec& beta,
                       const arma::vec& gamma,
                       const double delta,
                       const double omega,
                       const arma::vec E_g) {

  const int n = x.size();
  const int p = phi.size();
  const int q = beta.size();
  const int l = std::max(p, q);
  const int nl = n + l;

  const arma::vec x_dm = x - mu;

  const double E_sig_delta = omega / (1.0 - arma::as_scalar(phi * E_g) - arma::accu(beta));
  const arma::vec E_init = E_g * E_sig_delta;
  arma::vec sig_delta = arma::vec(nl).fill(E_sig_delta);

  const arma::rowvec phi_r = arma::reverse(phi);
  const arma::vec gamma_r = arma::reverse(gamma);
  const arma::rowvec beta_r = arma::reverse(beta);

  arma::vec g_fun = E_init;

  const int pm1 = p - 1;

  int jl;
  int jm1;

  for (int j = 0; j < p; ++j) {
    if (j > 0) {
      jm1 = j - 1;
      g_fun.subvec(p - j, pm1) = arma::pow(arma::abs(x_dm.subvec(0, jm1)) - gamma_r.subvec(pm1 + 1 - j, pm1) % x_dm.subvec(0, jm1), delta);
    }
    jl = j + l;
    sig_delta(jl) = omega +
      arma::as_scalar(phi_r * g_fun) +
      arma::as_scalar(beta_r * sig_delta.subvec(jl - q, jl - 1));
  }

  arma::vec xs(p);

  int lp = l + p;

  for (int i = l + p; i < nl; ++i) {

    xs = x_dm.subvec(i - lp, i - l - 1);
    g_fun = arma::pow(arma::abs(xs) - gamma_r % xs, delta);

    sig_delta(i) = omega +
      arma::as_scalar(phi_r * g_fun) +
      arma::as_scalar(beta_r * sig_delta.subvec(i - q, i - 1));

  }

  return arma::pow(sig_delta.subvec(l, nl - 1), 1.0 / delta);

}

// [[Rcpp::export]]
Rcpp::List aparch_sim_Cpp(const arma::vec& innov, const double omega,
                          const arma::rowvec& phi, const arma::rowvec& beta,
                          const arma::vec& gamma, const double delta,
                          const int n_out,
                          const double mu,
                          const arma::vec& E_g,
                          const double E_sigd) {

  const int n = innov.size();

  const int p = phi.size();
  const int q = beta.size();
  const int l = std::max(p, q);

  const int n_pre = l;
  const int n_all = n + n_pre;

  arma::vec sig_d(n_all);
  sig_d.fill(E_sigd);
  arma::vec g_fun = E_g;

  const arma::rowvec phi_rev = arma::reverse(phi);
  const arma::rowvec beta_rev = arma::reverse(beta);
  const arma::vec gamma_r = arma::reverse(gamma);

  int lj;

  for (int j = 0; j < p; ++j) {
    if (j > 0) {
      g_fun.subvec(p - j, p - 1) = arma::pow(arma::abs(innov.subvec(0, j - 1)) - gamma_r.subvec(p - j, p - 1) % innov.subvec(0, j - 1), delta);
    }
    lj = l + j;
    sig_d(lj) = omega +
      arma::as_scalar(phi_rev * (sig_d.subvec(lj - p, lj - 1) % g_fun)) +
      arma::as_scalar(beta_rev * sig_d.subvec(lj - q, lj - 1));
  }

  int im1;
  const int lp = l + p;

  for (int i = l + p; i < n_all; ++i) {

    im1 = i - 1;

    g_fun = arma::pow(arma::abs(innov.subvec(i - lp, im1 - l)) - gamma * innov.subvec(i - lp, im1 - l), delta);

    sig_d(i) = omega +
      arma::as_scalar(phi_rev * (sig_d.subvec(i - p, im1) % g_fun)) +
      arma::as_scalar(beta_rev * sig_d.subvec(i - q, im1));

  }

  const arma::vec sigt_out = arma::pow(sig_d.subvec(n_all - n_out, n_all - 1), 1.0 / delta);
  const arma::vec innov_out = innov.subvec(n_all - n_out, n_all - 1);

  const arma::vec rt_out = (sigt_out % innov_out) + mu;

  return Rcpp::List::create(_["rt"] = rt_out, _["sigt"] = sigt_out, _["etat"] = innov_out);
}

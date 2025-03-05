#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
arma::vec sigt_aparch_Cpp(const arma::vec& x, const double mu,
                       const arma::rowvec& phi, const arma::rowvec& beta,
                       const double gamma,
                       const double delta,
                       const double omega) {

  const int n = x.size();
  const int p = phi.size();
  const int q = beta.size();
  const int l = std::max(p, q);
  const int nl = n + l;

  const arma::vec x_dm = x - mu;

  const arma::vec transf = arma::pow(arma::abs(x_dm) - gamma * x_dm, delta);
  //Initialization based on sample
  const double E_g_delta = arma::mean(transf);
  //const double E_g_delta = 0.02;
  const double E_sig_delta = (omega + arma::accu(phi) * E_g_delta) / (1.0 - arma::accu(beta));
  arma::vec g_fun = arma::vec(nl).fill(E_g_delta);
  arma::vec sig_delta = arma::vec(nl).fill(E_sig_delta);
  g_fun.subvec(l, nl - 1) = transf;

  const arma::rowvec phi_r = arma::reverse(phi);
  const arma::rowvec beta_r = arma::reverse(beta);

  for (int i = l; i < nl; ++i) {

    sig_delta(i) = omega +
      arma::as_scalar(phi_r * g_fun.subvec(i - p, i - 1)) +
      arma::as_scalar(beta_r * sig_delta.subvec(i - q, i - 1));

  }

  return arma::pow(sig_delta.subvec(l, nl - 1), 1.0 / delta);

}

// [[Rcpp::export]]
Rcpp::List aparch_sim_Cpp(const arma::vec& innov, const double omega,
                          const arma::rowvec& phi, const arma::rowvec& beta,
                          const double gamma, const double delta,
                          const int n_out,
                          const double mu,
                          const double E_g,
                          const double E_sigd) {

  const int n = innov.size();

  const int p = phi.size();
  const int q = beta.size();
  const int l = std::max(p, q);

  const int n_pre = l;
  const int n_all = n + n_pre;

  arma::vec sig_d(n_all);
  sig_d.fill(E_sigd);
  arma::vec g_fun(n_all);
  g_fun.fill(E_g);
  g_fun.subvec(l, n_all - 1) = arma::pow(arma::abs(innov) - gamma * innov, delta);

  const arma::rowvec phi_rev = arma::reverse(phi);
  const arma::rowvec beta_rev = arma::reverse(beta);

  int im1;

  for (int i = l; i < n_all; ++i) {

    im1 = i - 1;

    sig_d(i) = omega +
      arma::as_scalar(phi_rev * (sig_d.subvec(i - p, im1) % g_fun.subvec(i - p, im1))) +
      arma::as_scalar(beta_rev * sig_d.subvec(i - q, im1));

  }

  const arma::vec sigt_out = arma::pow(sig_d.subvec(n_all - n_out, n_all - 1), 1.0 / delta);
  const arma::vec innov_out = innov.subvec(n_all - n_out, n_all - 1);

  const arma::vec rt_out = (sigt_out % innov_out) + mu;

  return Rcpp::List::create(_["rt"] = rt_out, _["sigt"] = sigt_out, _["etat"] = innov_out);
}

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
arma::vec sigt_garch_Cpp(const arma::vec& x, const double mu,
                       const arma::rowvec& phi, const arma::rowvec& beta,
                       const double omega) {

  const int n = x.size();
  const int p = phi.size();
  const int q = beta.size();
  const int l = std::max(p, q);
  const int nl = n + l;

  const arma::vec x_dm = x - mu;
  const arma::vec Y2 = arma::pow(x_dm, 2.0);

  //Initialization based on sample
  const double E_Y2 = omega / (1 - arma::accu(phi) - arma::accu(beta));
  const double E_sig2 = E_Y2;

  arma::vec Y_fun = arma::vec(nl).fill(E_Y2);
  arma::vec sig_2 = arma::vec(nl).fill(E_sig2);
  Y_fun.subvec(l, nl - 1) = Y2;

  const arma::rowvec phi_r = arma::reverse(phi);
  const arma::rowvec beta_r = arma::reverse(beta);

  for (int i = l; i < nl; ++i) {

    sig_2(i) = omega +
      arma::as_scalar(phi_r * Y_fun.subvec(i - p, i - 1)) +
      arma::as_scalar(beta_r * sig_2.subvec(i - q, i - 1));

  }

  return arma::sqrt(sig_2.subvec(l, nl - 1));

}

// [[Rcpp::export]]
Rcpp::List garch_sim_Cpp(const arma::vec& innov, const double omega,
                          const arma::rowvec& phi, const arma::rowvec& beta,
                          const int n_out,
                          const double mu) {

  const int n = innov.size();

  const int p = phi.size();
  const int q = beta.size();
  const int l = std::max(p, q);
  const int n_pre = l;
  const int n_all = n + n_pre;

  const double EY2 = omega / (1 - arma::accu(phi) - arma::accu(beta));

  arma::vec sig2 = arma::vec(n_all).fill(EY2);

  const arma::rowvec phi_r = arma::reverse(phi);
  const arma::rowvec beta_r = arma::reverse(beta);

  arma::vec innov2 = arma::ones(n_all);
  innov2.subvec(l, n_all - 1) = arma::pow(innov, 2.0);

  for (int i = l; i < n_all; ++i) {

    sig2(i) = omega +
      arma::as_scalar(phi_r * (sig2.subvec(i - p, i - 1) % innov2.subvec(i - p, i - 1))) +
      arma::as_scalar(beta_r * sig2.subvec(i - q, i - 1));

  }

  const arma::vec sigt_out = sig2.subvec(n_all - n_out, n_all - 1);
  const arma::vec innov_out = innov.subvec(n - n_out, n - 1);

  const arma::vec rt_out = sigt_out % innov_out + mu;

  return Rcpp::List::create(_["rt"] = rt_out, _["sigt"] = sigt_out, _["etat"] = innov_out);
}

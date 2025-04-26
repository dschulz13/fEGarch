#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
arma::vec sigt_garch_Cpp(const arma::vec& x, const double mu,
                       const arma::vec& phi, const arma::vec& beta,
                       const double omega, const double sig2_init) {

  const int n = x.size();
  const int p = phi.size();
  const int q = beta.size();
  const int l = std::max(p, q);
  const int nl = n + l;

  const arma::vec x_dm = x - mu;
  const arma::vec Y2 = arma::pow(x_dm, 2.0);

  arma::vec Y_fun = arma::vec(nl).fill(sig2_init);
  arma::vec sig_2 = arma::vec(nl).fill(sig2_init);
  Y_fun.subvec(l, nl - 1) = Y2;

  const arma::vec phi_r = arma::reverse(phi);
  const arma::vec beta_r = arma::reverse(beta);

  for (int i = l; i < nl; ++i) {

    sig_2(i) = omega +
      arma::as_scalar(arma::dot(phi_r, Y_fun.subvec(i - p, i - 1))) +
      arma::as_scalar(arma::dot(beta_r, sig_2.subvec(i - q, i - 1)));

  }

  return arma::sqrt(sig_2.subvec(l, nl - 1));

}

// [[Rcpp::export]]
Rcpp::List garch_sim_Cpp(const arma::vec& innov, const double omega,
                          const arma::vec& phi, const arma::vec& beta,
                          const int n_out,
                          const double mu, const double EY2) {

  const int n = innov.size();

  const int p = phi.size();
  const int q = beta.size();
  const int l = std::max(p, q);
  const int n_pre = l;
  const int n_all = n + n_pre;

  arma::vec sig2 = arma::vec(n_all).fill(EY2);
  arma::vec sig = arma::zeros<arma::vec>(n_all);

  const arma::vec phi_r = arma::reverse(phi);
  const arma::vec beta_r = arma::reverse(beta);

  arma::vec Yt2 = arma::vec(n_all).fill(EY2);
  arma::vec Yt = arma::vec(n);

  int im1;
  int imp;
  int imq;
  int iml;

  double sig2i;

  for (int i = l; i < n_all; ++i) {

    im1 = i - 1;
    imp = i - p;
    imq = i - q;
    iml = i - l;

    sig2i = omega +
      arma::as_scalar(arma::dot(phi_r, Yt2.subvec(imp, im1))) +
      arma::as_scalar(arma::dot(beta_r, sig2.subvec(imq, im1)));

    sig2(i) = sig2i;
    sig(iml) = std::sqrt(sig2i);

    Yt(iml) = sig(iml) * innov(iml);
    Yt2(i) = std::pow(Yt(iml), 2.0);

  }

  const arma::vec sigt_out = sig.subvec(n - n_out, n - 1);
  const arma::vec innov_out = innov.subvec(n - n_out, n - 1);

  const arma::vec rt_out = Yt.subvec(n - n_out, n - 1) + mu;

  return Rcpp::List::create(_["rt"] = rt_out, _["sigt"] = sigt_out, _["etat"] = innov_out);
}

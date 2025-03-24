#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List figarch_sim_Cpp(const arma::vec& innov, const double omega,
                          const arma::vec& coef_inf,
                          const int n_out,
                          const double mu) {

  const int n = innov.size();

  double sig2;
  arma::vec sigt = arma::vec(n);

  const arma::vec coef_rev = arma::reverse(coef_inf);

  arma::vec rt = arma::zeros<arma::vec>(n);
  arma::vec eps2 = arma::zeros<arma::vec>(n);

  const int lc = coef_inf.size();

  const int n_coef_m1 = lc - 1;


  for (int j = 0; j < n_coef_m1; ++j) {


    sig2 = omega +
      arma::as_scalar(arma::dot(coef_rev.subvec(n_coef_m1 - j, n_coef_m1), eps2.subvec(0, j)));

    sigt(j) = std::sqrt(sig2);

    rt(j) = innov(j) * sigt(j);
    eps2(j) = std::pow(rt(j), 2.0);

  }

  for (int i = n_coef_m1; i < n; ++i) {

    sig2 = omega +
      arma::as_scalar(arma::dot(coef_rev, eps2.subvec(0, i)));

    sigt(i) = std::sqrt(sig2);

    rt(i) = innov(i) * sigt(i);
    eps2(i) = std::pow(rt(i), 2.0);

  }

  const arma::vec sigt_out = sigt.subvec(n - n_out, n - 1);
  const arma::vec innov_out = innov.subvec(n - n_out, n - 1);

  const arma::vec rt_out = rt.subvec(n - n_out, n - 1) + mu;

  return Rcpp::List::create(_["rt"] = rt_out, _["sigt"] = sigt_out, _["etat"] = innov_out);
}

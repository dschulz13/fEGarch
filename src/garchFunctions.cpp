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
  const double E_Y2 = arma::mean(Y2);
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

// // [[Rcpp::export]]
// Rcpp::List garch_sim_Cpp(const arma::vec& innov, const double omega,
//                           const arma::rowvec& phi, const arma::rowvec& beta,
//                           const double gamma, const double delta,
//                           const int n_out,
//                           const double mu) {
//
//   const int n = innov.size();
//
//   const int n_pre = 1;
//   const int n_all = n + n_pre;
//
//   double sig_d;
//   arma::vec sigt = arma::vec(n);
//
//   const arma::rowvec coef_rev = arma::reverse(coef_inf);
//
//   arma::vec rt = arma::zeros<arma::vec>(n_all);
//   arma::vec g_fun = arma::zeros<arma::vec>(n_all);
//
//   const int n_coef_m1 = coef_inf.size() - 1;
//   int npp1;
//
//   for (int i = 0; i < n; ++i) {
//
//     npp1 = 1 + i;
//
//     sig_d = omega +
//       arma::as_scalar(coef_rev.subvec(n_coef_m1 - i, n_coef_m1) * g_fun.subvec(0, i));
//
//     sigt(i) = std::pow(sig_d, 1.0 / delta);
//
//     rt(npp1) = innov(i) * sigt(i);
//     g_fun(npp1) = std::pow(std::abs(rt(npp1)) - gamma * rt(npp1), delta);
//
//   }
//
//   const arma::vec sigt_out = sigt.subvec(n - n_out, n - 1);
//   const arma::vec innov_out = innov.subvec(n - n_out, n - 1);
//
//   const arma::vec rt_out = rt.subvec(n_all - n_out, n_all - 1) + mu;
//
//   return Rcpp::List::create(_["rt"] = rt_out, _["sigt"] = sigt_out, _["etat"] = innov_out);
// }

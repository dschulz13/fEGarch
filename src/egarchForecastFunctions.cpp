#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
#include "signVec.h"

// [[Rcpp::export]]
arma::vec sigt_egarch__forecast_shortCpp(const arma::vec& x,
                               const arma::vec& et,
                               const arma::vec& sigt,
                               const double mu,
                               const double omega,
                               const arma::rowvec& phi,
                               const arma::rowvec& alpha,
                               const arma::rowvec& beta,
                               const double E_asy,
                               const double E_mag,
                               const double Elnsig2,
                               const arma::vec& powers,
                               const arma::vec& modulus,
                               const std::string& mode,
                               const int horizon) {

  const int n = x.size();
  const int p = phi.size();
  const int q = alpha.size();
  const int l = std::max(p, q);
  const int Hl = horizon + l;

  int il, il1, ilp, ilq;

  const double mod1 = modulus(0);
  const double mod2 = modulus(1);
  const double pow1 = powers.size() > 0 ? powers(0) : 1.0; // Default to 1.0 for safety
  const double pow2 = powers.size() > 1 ? powers(1) : 1.0;

  arma::vec transf1(l);
  arma::vec transf2(l);
  const arma::vec eta_i = et.subvec(n - l, n - 1);

  if (mode == "pow_pow") {
    transf1 = (signVecCpp(eta_i) * arma::pow(arma::abs(eta_i) + mod1, pow1) - E_asy) / pow1;
    transf2 = (arma::pow(arma::abs(eta_i) + mod2, pow2) - E_mag) / pow2;
  } else if (mode == "log_log") {
    transf1 = signVecCpp(eta_i) * arma::log(arma::abs(eta_i) + mod1) - E_asy;
    transf2 = arma::log(arma::abs(eta_i) + mod2) - E_mag;
  } else if (mode == "pow_log") {
    transf1 = (signVecCpp(eta_i) * arma::pow(arma::abs(eta_i) + mod1, pow1) - E_asy) / pow1;
    transf2 = arma::log(arma::abs(eta_i) + mod2) - E_mag;
  } else if (mode == "log_pow") {
    transf1 = signVecCpp(eta_i) * arma::log(arma::abs(eta_i) + mod1) - E_asy;
    transf2 = (arma::pow(arma::abs(eta_i) + mod2, pow2) - E_mag) / pow2;
  } else {
    Rcpp::stop("Invalid mode specified");
  }

  arma::vec et_asy = arma::zeros<arma::vec>(Hl);
  et_asy.subvec(0, l - 1) = transf1;

  arma::vec et_mag = arma::zeros<arma::vec>(Hl);
  et_mag.subvec(0, l - 1) = transf2;

  arma::vec lnsig2 = arma::zeros<arma::vec>(Hl);
  lnsig2.subvec(0, l - 1) = arma::log(arma::pow(sigt.subvec(n - l, n - 1), 2.0));

  const arma::rowvec alpha_r = arma::reverse(alpha);
  const arma::rowvec phi_r = arma::reverse(phi);
  const arma::rowvec beta_r = arma::reverse(beta);

  for (int i = 0; i < horizon; ++i) {
    il = i + l;
    il1 = il - 1;
    ilp = il - p;
    ilq = il - q;

    lnsig2(il) = omega +
      arma::as_scalar(alpha_r * et_asy.subvec(ilq, il1)) +
      arma::as_scalar(beta_r * et_mag.subvec(ilq, il1)) +
      arma::as_scalar(phi_r * lnsig2.subvec(ilp, il1));

  }

  return arma::exp(0.5 * lnsig2.subvec(l, Hl - 1));
}

//
// // [[Rcpp::export]]
// arma::vec sigt_egarch_longCpp(const arma::vec& x, const double mu,
//                        const arma::rowvec& coef_inf,
//                        const double kappa,
//                        const double gamma,
//                        const double E_asy,
//                        const double E_mag,
//                        const double Elnsig2,
//                        const arma::vec& powers,
//                        const arma::vec& modulus,
//                        const std::string& mode) {
//
//   const int n = x.size();
//   const int n1 = n + 1;
//   double lnsig2 = 0;
//   arma::vec series_g = arma::vec(n1).fill(0);
//
//   const arma::vec x_cent = x - mu;
//   arma::vec sigt = arma::zeros<arma::vec>(n);
//
//   const arma::rowvec coef_r = arma::reverse(coef_inf);
//   double eta_i;
//
//   int il;
//   int nm1 = n - 1;
//   const double mod1 = modulus(0);
//   const double mod2 = modulus(1);
//   const double pow1 = powers.size() > 0 ? powers(0) : 1.0; // Default to 1.0 for safety
//   const double pow2 = powers.size() > 1 ? powers(1) : 1.0;
//
//   std::function<double(double)> compute_abs_eta_pow1;
//   std::function<double(double)> compute_abs_eta_pow2;
//
//   if (mode == "pow_pow") {
//     compute_abs_eta_pow1 = [pow1, mod1, E_asy](double eta_i) {
//       return (signCpp(eta_i) * std::pow(std::abs(eta_i) + mod1, pow1) - E_asy) / pow1;
//     };
//     compute_abs_eta_pow2 = [pow2, mod2, E_mag](double eta_i) {
//       return (std::pow(std::abs(eta_i) + mod2, pow2) - E_mag) / pow2;
//     };
//   } else if (mode == "log_log") {
//     compute_abs_eta_pow1 = [mod1, E_asy](double eta_i) {
//       return signCpp(eta_i) * std::log(std::abs(eta_i) + mod1) - E_asy;
//     };
//     compute_abs_eta_pow2 = [mod2, E_mag](double eta_i) {
//       return std::log(std::abs(eta_i) + mod2) - E_mag;
//     };
//   } else if (mode == "pow_log") {
//     compute_abs_eta_pow1 = [pow1, mod1, E_asy](double eta_i) {
//       return (signCpp(eta_i) * std::pow(std::abs(eta_i) + mod1, pow1) - E_asy) / pow1;
//     };
//     compute_abs_eta_pow2 = [mod2, E_mag](double eta_i) {
//       return std::log(std::abs(eta_i) + mod2) - E_mag;
//     };
//   } else if (mode == "log_pow") {
//     compute_abs_eta_pow1 = [mod1, E_asy](double eta_i) {
//       return signCpp(eta_i) * std::log(std::abs(eta_i) + mod1) - E_asy;
//     };
//     compute_abs_eta_pow2 = [pow2, mod2, E_mag](double eta_i) {
//       return (std::pow(std::abs(eta_i) + mod2, pow2) - E_mag) / pow2;
//     };
//   } else {
//     Rcpp::stop("Invalid mode specified");
//   }
//
//   for (int i = 0; i < n; ++i) {
//
//     il = i + 1;
//
//     lnsig2 = Elnsig2 +
//       arma::as_scalar(coef_r.subvec(nm1 - i, nm1) * series_g.subvec(0, i));
//
//     sigt(i) = std::exp(lnsig2 / 2.0);
//     eta_i = x_cent(i) / sigt(i);
//
//     series_g(il) = kappa * compute_abs_eta_pow1(eta_i) + gamma * compute_abs_eta_pow2(eta_i);
//
//   }
//
//   return sigt;
//
// }


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
int signCpp(const double val) {
    return (0.0 < val) - (val < 0.0);
}

// // [[Rcpp::export]]
// double powTransfAsyCpp(const double val, const double power, const double mod,
//                          const double E) {
//     return signCpp(val) * (std::pow(std::abs(val) + mod, power) - E) / power;
// }
//
// // [[Rcpp::export]]
// double powTransfMagCpp(const double val, const double power, const double mod,
//                          const double E) {
//     return (std::pow(std::abs(val) + mod, power) - E) / power;
// }
//
// // [[Rcpp::export]]
// double logTransfAsyCpp(const double val, const double mod,
//                          const double E) {
//     return signCpp(val) * (std::log(std::abs(val) + mod) - E);
// }
//
// // [[Rcpp::export]]
// double logTransfMagCpp(const double val, const double mod,
//                          const double E) {
//     return std::log(std::abs(val) + mod) - E;
// }


// [[Rcpp::export]]
arma::vec sigt_egarch_short_pp(const arma::vec& x, const double mu,
                       const double omega,
                       const arma::rowvec& phi,
                       const arma::rowvec& alpha,
                       const arma::rowvec& beta,
                       const double E1,
                       const double E2,
                       const double Elnsig2,
                       const arma::vec& powers,
                       const arma::vec& modulus) {

  const int n = x.size();
  const int p = phi.size();
  const int q = alpha.size();
  const int l = std::max(p, q);
  const int nl = n + l;
  arma::vec lnsig2 = arma::vec(nl).fill(Elnsig2);
  double eta_i;
  arma::vec series1 = arma::vec(nl).fill(0);
  arma::vec series2 = arma::vec(nl).fill(0);

  const arma::vec x_cent = x - mu;
  arma::vec sigt = arma::zeros<arma::vec>(n);


  const arma::rowvec alpha_r = arma::reverse(alpha);
  const arma::rowvec phi_r = arma::reverse(phi);
  const arma::rowvec beta_r = arma::reverse(beta);

  int il;
  int il1;
  int ilp;
  int ilq;

  const double pow1 = powers(0);
  const double pow2 = powers(1);
  const double mod1 = modulus(0);
  const double mod2 = modulus(1);
  double abs_eta_pow1;
  double abs_eta_pow2;

  for (int i = 0; i < n; ++i) {

    il = i + l;
    il1 = il - 1;
    ilp = il - p;
    ilq = il - q;
    lnsig2(il) = omega +
      arma::as_scalar(alpha_r * series1.subvec(ilq, il1)) +
      arma::as_scalar(beta_r * series2.subvec(ilq, il1)) +
      arma::as_scalar(phi_r * lnsig2.subvec(ilp, il1));

    sigt(i) = std::exp(lnsig2(il) / 2.0);
    eta_i = x_cent(i) / sigt(i);
    abs_eta_pow1 = signCpp(eta_i) * std::pow(std::abs(eta_i) + mod1, pow1);
    abs_eta_pow2 = std::pow(std::abs(eta_i) + mod2, pow2);
    series1(il) = (abs_eta_pow1 - E1) / pow1;
    series2(il) = (abs_eta_pow2 - E2) / pow2;

  }

  return sigt;

}

// [[Rcpp::export]]
arma::vec sigt_egarch_short_loglog(const arma::vec& x, const double mu,
                       const double omega,
                       const arma::rowvec& phi,
                       const arma::rowvec& alpha,
                       const arma::rowvec& beta,
                       const double E1,
                       const double E2,
                       const double Elnsig2,
                       const arma::vec& modulus) {

  const int n = x.size();
  const int p = phi.size();
  const int q = alpha.size();
  const int l = std::max(p, q);
  const int nl = n + l;
  arma::vec lnsig2 = arma::vec(nl).fill(Elnsig2);
  double eta_i;
  arma::vec series1 = arma::vec(nl).fill(0);
  arma::vec series2 = arma::vec(nl).fill(0);

  const arma::vec x_cent = x - mu;
  arma::vec sigt = arma::zeros<arma::vec>(n);


  const arma::rowvec alpha_r = arma::reverse(alpha);
  const arma::rowvec phi_r = arma::reverse(phi);
  const arma::rowvec beta_r = arma::reverse(beta);

  int il;
  int il1;
  int ilp;
  int ilq;

  const double mod1 = modulus(0);
  const double mod2 = modulus(1);
  double abs_eta_pow1;
  double abs_eta_pow2;

  for (int i = 0; i < n; ++i) {

    il = i + l;
    il1 = il - 1;
    ilp = il - p;
    ilq = il - q;
    lnsig2(il) = omega +
      arma::as_scalar(alpha_r * series1.subvec(ilq, il1)) +
      arma::as_scalar(beta_r * series2.subvec(ilq, il1)) +
      arma::as_scalar(phi_r * lnsig2.subvec(ilp, il1));

    sigt(i) = std::exp(lnsig2(il) / 2.0);
    eta_i = x_cent(i) / sigt(i);
    abs_eta_pow1 = signCpp(eta_i) * std::log(std::abs(eta_i) + mod1);
    abs_eta_pow2 = std::log(std::abs(eta_i) + mod2);
    series1(il) = abs_eta_pow1 - E1;
    series2(il) = abs_eta_pow2 - E2;

  }

  return sigt;

}

// [[Rcpp::export]]
arma::vec sigt_egarch_short_plog(const arma::vec& x, const double mu,
                       const double omega,
                       const arma::rowvec& phi,
                       const arma::rowvec& alpha,
                       const arma::rowvec& beta,
                       const double E1,
                       const double E2,
                       const double Elnsig2,
                       const double pow1,
                       const arma::vec& modulus) {

  const int n = x.size();
  const int p = phi.size();
  const int q = alpha.size();
  const int l = std::max(p, q);
  const int nl = n + l;
  arma::vec lnsig2 = arma::vec(nl).fill(Elnsig2);
  double eta_i;
  arma::vec series1 = arma::vec(nl).fill(0);
  arma::vec series2 = arma::vec(nl).fill(0);

  const arma::vec x_cent = x - mu;
  arma::vec sigt = arma::zeros<arma::vec>(n);


  const arma::rowvec alpha_r = arma::reverse(alpha);
  const arma::rowvec phi_r = arma::reverse(phi);
  const arma::rowvec beta_r = arma::reverse(beta);

  int il;
  int il1;
  int ilp;
  int ilq;

  const double mod1 = modulus(0);
  const double mod2 = modulus(1);
  double abs_eta_pow1;
  double abs_eta_pow2;

  for (int i = 0; i < n; ++i) {

    il = i + l;
    il1 = il - 1;
    ilp = il - p;
    ilq = il - q;
    lnsig2(il) = omega +
      arma::as_scalar(alpha_r * series1.subvec(ilq, il1)) +
      arma::as_scalar(beta_r * series2.subvec(ilq, il1)) +
      arma::as_scalar(phi_r * lnsig2.subvec(ilp, il1));

    sigt(i) = std::exp(lnsig2(il) / 2.0);
    eta_i = x_cent(i) / sigt(i);
    abs_eta_pow1 = signCpp(eta_i) * std::pow(std::abs(eta_i) + mod1, pow1);
    abs_eta_pow2 = std::log(std::abs(eta_i) + mod2);
    series1(il) = (abs_eta_pow1 - E1) / pow1;
    series2(il) = abs_eta_pow2 - E2;

  }

  return sigt;

}

// [[Rcpp::export]]
arma::vec sigt_egarch_short_logp(const arma::vec& x, const double mu,
                       const double omega,
                       const arma::rowvec& phi,
                       const arma::rowvec& alpha,
                       const arma::rowvec& beta,
                       const double E1,
                       const double E2,
                       const double Elnsig2,
                       const double pow2,
                       const arma::vec& modulus) {

  const int n = x.size();
  const int p = phi.size();
  const int q = alpha.size();
  const int l = std::max(p, q);
  const int nl = n + l;
  arma::vec lnsig2 = arma::vec(nl).fill(Elnsig2);
  double eta_i;
  arma::vec series1 = arma::vec(nl).fill(0);
  arma::vec series2 = arma::vec(nl).fill(0);

  const arma::vec x_cent = x - mu;
  arma::vec sigt = arma::zeros<arma::vec>(n);


  const arma::rowvec alpha_r = arma::reverse(alpha);
  const arma::rowvec phi_r = arma::reverse(phi);
  const arma::rowvec beta_r = arma::reverse(beta);

  int il;
  int il1;
  int ilp;
  int ilq;

  const double mod1 = modulus(0);
  const double mod2 = modulus(1);
  double abs_eta_pow1;
  double abs_eta_pow2;

  for (int i = 0; i < n; ++i) {

    il = i + l;
    il1 = il - 1;
    ilp = il - p;
    ilq = il - q;
    lnsig2(il) = omega +
      arma::as_scalar(alpha_r * series1.subvec(ilq, il1)) +
      arma::as_scalar(beta_r * series2.subvec(ilq, il1)) +
      arma::as_scalar(phi_r * lnsig2.subvec(ilp, il1));

    sigt(i) = std::exp(lnsig2(il) / 2.0);
    eta_i = x_cent(i) / sigt(i);
    abs_eta_pow1 = signCpp(eta_i) * std::log(std::abs(eta_i) + mod1);
    abs_eta_pow2 = std::pow(std::abs(eta_i) + mod2, pow2);
    series1(il) = abs_eta_pow1 - E1;
    series2(il) = (abs_eta_pow2 - E2) / pow2;

  }

  return sigt;

}


////////////////////////////////////////


// [[Rcpp::export]]
arma::vec sigt_egarch_long_pp(const arma::vec& x, const double mu,
                       const arma::rowvec& alpha,
                       const arma::rowvec& beta,
                       const double E1,
                       const double E2,
                       const double Elnsig2,
                       const arma::vec& powers,
                       const arma::vec& modulus) {

  const int n = x.size();
  const int n1 = n + 1;
  arma::vec lnsig2 = arma::vec(n1).fill(0);
  arma::vec series_alpha = arma::vec(n1).fill(0);
  arma::vec series_beta = arma::vec(n1).fill(0);

  const arma::vec x_cent = x - mu;
  arma::vec sigt = arma::zeros<arma::vec>(n);

  const arma::rowvec alpha_r = arma::reverse(alpha);
  const arma::rowvec beta_r = arma::reverse(beta);
  double eta_i;

  int il;
  int nm1 = n - 1;
  const double pow1 = powers(0);
  const double pow2 = powers(1);
  const double mod1 = modulus(0);
  const double mod2 = modulus(1);
  double abs_eta_pow1;
  double abs_eta_pow2;

  for (int i = 0; i < n; ++i) {

    il = i + 1;

    lnsig2(il) = Elnsig2 +
      arma::as_scalar(alpha_r.subvec(nm1 - i, nm1) * series_alpha.subvec(0, i)) +
      arma::as_scalar(beta_r.subvec(nm1 - i, nm1) * series_beta.subvec(0, i));

    sigt(i) = std::exp(lnsig2(il) / 2.0);
    eta_i = x_cent(i) / sigt(i);
    abs_eta_pow1 = signCpp(eta_i) * std::pow((std::abs(eta_i) + mod1), pow1);
    abs_eta_pow2 = std::pow((std::abs(eta_i) + mod2), pow2);
    series_alpha(il) = (abs_eta_pow1 - E1) / pow1;
    series_beta(il) = (abs_eta_pow2 - E2) / pow2;

  }

  return sigt;

}

// [[Rcpp::export]]
arma::vec sigt_egarch_long_loglog(const arma::vec& x, const double mu,
                       const arma::rowvec& alpha,
                       const arma::rowvec& beta,
                       const double E1,
                       const double E2,
                       const double Elnsig2,
                       const arma::vec& modulus) {

  const int n = x.size();
  const int n1 = n + 1;
  arma::vec lnsig2 = arma::vec(n1).fill(0);
  arma::vec series_alpha = arma::vec(n1).fill(0);
  arma::vec series_beta = arma::vec(n1).fill(0);

  const arma::vec x_cent = x - mu;
  arma::vec sigt = arma::zeros<arma::vec>(n);

  const arma::rowvec alpha_r = arma::reverse(alpha);
  const arma::rowvec beta_r = arma::reverse(beta);
  double eta_i;

  int il;
  int nm1 = n - 1;
  const double mod1 = modulus(0);
  const double mod2 = modulus(1);
  double abs_eta_pow1;
  double abs_eta_pow2;

  for (int i = 0; i < n; ++i) {

    il = i + 1;

    lnsig2(il) = Elnsig2 +
      arma::as_scalar(alpha_r.subvec(nm1 - i, nm1) * series_alpha.subvec(0, i)) +
      arma::as_scalar(beta_r.subvec(nm1 - i, nm1) * series_beta.subvec(0, i));

    sigt(i) = std::exp(lnsig2(il) / 2.0);
    eta_i = x_cent(i) / sigt(i);
    abs_eta_pow1 = signCpp(eta_i) * std::log(std::abs(eta_i) + mod1);
    abs_eta_pow2 = std::log(std::abs(eta_i) + mod2);
    series_alpha(il) = abs_eta_pow1 - E1;
    series_beta(il) = abs_eta_pow2 - E2;

  }

  return sigt;

}

// [[Rcpp::export]]
arma::vec sigt_egarch_long_plog(const arma::vec& x, const double mu,
                       const arma::rowvec& alpha,
                       const arma::rowvec& beta,
                       const double E1,
                       const double E2,
                       const double Elnsig2,
                       const double pow1,
                       const arma::vec& modulus) {

  const int n = x.size();
  const int n1 = n + 1;
  arma::vec lnsig2 = arma::vec(n1).fill(0);
  arma::vec series_alpha = arma::vec(n1).fill(0);
  arma::vec series_beta = arma::vec(n1).fill(0);

  const arma::vec x_cent = x - mu;
  arma::vec sigt = arma::zeros<arma::vec>(n);

  const arma::rowvec alpha_r = arma::reverse(alpha);
  const arma::rowvec beta_r = arma::reverse(beta);
  double eta_i;

  int il;
  int nm1 = n - 1;
  const double mod1 = modulus(0);
  const double mod2 = modulus(1);
  double abs_eta_pow1;
  double abs_eta_pow2;

  for (int i = 0; i < n; ++i) {

    il = i + 1;

    lnsig2(il) = Elnsig2 +
      arma::as_scalar(alpha_r.subvec(nm1 - i, nm1) * series_alpha.subvec(0, i)) +
      arma::as_scalar(beta_r.subvec(nm1 - i, nm1) * series_beta.subvec(0, i));

    sigt(i) = std::exp(lnsig2(il) / 2.0);
    eta_i = x_cent(i) / sigt(i);
    abs_eta_pow1 = signCpp(eta_i) * std::pow(std::abs(eta_i) + mod1, pow1);
    abs_eta_pow2 = std::log(std::abs(eta_i) + mod2);
    series_alpha(il) = (abs_eta_pow1 - E1) / pow1;
    series_beta(il) = abs_eta_pow2 - E2;

  }

  return sigt;

}

// [[Rcpp::export]]
arma::vec sigt_egarch_long_logp(const arma::vec& x, const double mu,
                       const arma::rowvec& alpha,
                       const arma::rowvec& beta,
                       const double E1,
                       const double E2,
                       const double Elnsig2,
                       const double pow2,
                       const arma::vec& modulus) {

  const int n = x.size();
  const int n1 = n + 1;
  arma::vec lnsig2 = arma::vec(n1).fill(0);
  arma::vec series_alpha = arma::vec(n1).fill(0);
  arma::vec series_beta = arma::vec(n1).fill(0);

  const arma::vec x_cent = x - mu;
  arma::vec sigt = arma::zeros<arma::vec>(n);

  const arma::rowvec alpha_r = arma::reverse(alpha);
  const arma::rowvec beta_r = arma::reverse(beta);
  double eta_i;

  int il;
  int nm1 = n - 1;
  const double mod1 = modulus(0);
  const double mod2 = modulus(1);
  double abs_eta_pow1;
  double abs_eta_pow2;

  for (int i = 0; i < n; ++i) {

    il = i + 1;

    lnsig2(il) = Elnsig2 +
      arma::as_scalar(alpha_r.subvec(nm1 - i, nm1) * series_alpha.subvec(0, i)) +
      arma::as_scalar(beta_r.subvec(nm1 - i, nm1) * series_beta.subvec(0, i));

    sigt(i) = std::exp(lnsig2(il) / 2.0);
    eta_i = x_cent(i) / sigt(i);
    abs_eta_pow1 = signCpp(eta_i) * std::log(std::abs(eta_i) + mod1);
    abs_eta_pow2 = std::pow(std::abs(eta_i) + mod2, pow2);
    series_alpha(il) = abs_eta_pow1 - E1;
    series_beta(il) = (abs_eta_pow2 - E2) / pow2;

  }

  return sigt;

}

////////////////////////////////////////

// [[Rcpp::export]]
arma::vec sigt_filoggarchCpp(const arma::vec& x, const double mu,
                       const arma::rowvec& coef_inf,
                       const double E1,
                       const double Elnsig2) {

  const int n = x.size();
  const int n1 = n + 1;
  arma::vec lnsig2 = arma::vec(n1).fill(0);
  arma::vec epst = arma::vec(n1).fill(0);

  const arma::vec x_cent = x - mu;
  arma::vec sigt = arma::zeros<arma::vec>(n);

  const arma::rowvec coef_inf_r = arma::reverse(coef_inf);
  double eta_i;

  int il;
  int nm1 = n - 1;

  for (int i = 0; i < n; ++i) {

    il = i + 1;

    lnsig2(il) = Elnsig2 +
      arma::as_scalar(coef_inf_r.subvec(nm1 - i, nm1) * epst.subvec(0, i));

    sigt(i) = std::exp(lnsig2(il) / 2.0);
    eta_i = x_cent(i) / sigt(i);
    epst(il) = std::log(std::pow(eta_i, 2.0)) - E1;

  }

  return sigt;

}

// [[Rcpp::export]]
arma::vec sigt_mafiloggarchCpp(const arma::vec& x, const double mu,
                       const arma::rowvec& alpha,
                       const arma::rowvec& beta,
                       const double E1,
                       const double E2,
                       const double Elnsig2) {

  const int n = x.size();
  const int n1 = n + 1;
  arma::vec lnsig2 = arma::vec(n1).fill(0);
  arma::vec series_alpha = arma::vec(n1).fill(0);
  arma::vec series_beta = arma::vec(n1).fill(0);

  const arma::vec x_cent = x - mu;
  arma::vec sigt = arma::zeros<arma::vec>(n);

  const arma::rowvec alpha_r = arma::reverse(alpha);
  const arma::rowvec beta_r = arma::reverse(beta);
  double eta_i;

  int il;
  int nm1 = n - 1;
  double zeta_i;

  for (int i = 0; i < n; ++i) {

    il = i + 1;

    lnsig2(il) = Elnsig2 +
      arma::as_scalar(alpha_r.subvec(nm1 - i, nm1) * series_alpha.subvec(0, i)) +
      arma::as_scalar(beta_r.subvec(nm1 - i, nm1) * series_beta.subvec(0, i));

    sigt(i) = std::exp(lnsig2(il) / 2.0);
    eta_i = x_cent(i) / sigt(i);
    zeta_i = signCpp(eta_i) * std::log(std::abs(eta_i) + 1) - E1;
    series_alpha(il) = zeta_i;
    series_beta(il) = std::abs(zeta_i) - E2;

  }

  return sigt;

}


// [[Rcpp::export]]
arma::vec sigt_megarchCpp(const arma::vec& x, const double mu,
                       const double omega,
                       const arma::rowvec& phi,
                       const arma::rowvec& alpha,
                       const arma::rowvec& beta,
                       const double E1,
                       const double E2,
                       const double Elnsig2) {

  const int n = x.size();
  const int p = phi.size();
  const int q = alpha.size();
  const int l = std::max(p, q);
  const int nl = n + l;
  arma::vec lnsig2 = arma::vec(nl).fill(Elnsig2);
  double eta_i;
  arma::vec series1 = arma::vec(nl).fill(0);
  arma::vec series2 = arma::vec(nl).fill(0);

  const arma::vec x_cent = x - mu;
  arma::vec sigt = arma::zeros<arma::vec>(n);

  const arma::rowvec alpha_r = arma::reverse(alpha);
  const arma::rowvec phi_r = arma::reverse(phi);
  const arma::rowvec beta_r = arma::reverse(beta);

  int il;
  int il1;
  int ilp;
  int ilq;

  for (int i = 0; i < n; ++i) {

    il = i + l;
    il1 = il - 1;
    ilp = il - p;
    ilq = il - q;
    lnsig2(il) = omega +
      arma::as_scalar(alpha_r * series1.subvec(ilq, il1)) +
      arma::as_scalar(beta_r * series2.subvec(ilq, il1)) +
      arma::as_scalar(phi_r * lnsig2.subvec(ilp, il1));

    sigt(i) = std::exp(lnsig2(il) / 2.0);
    eta_i = x_cent(i) / sigt(i);
    series1(il) = signCpp(eta_i) * std::log(std::abs(eta_i) + 1.0) - E1;
    series2(il) = std::abs(eta_i) - E2;

  }

  return sigt;

}

// [[Rcpp::export]]
arma::vec sigt_mloggarchCpp(const arma::vec& x, const double mu,
                       const double omega,
                       const arma::rowvec& phi,
                       const arma::rowvec& alpha,
                       const arma::rowvec& beta,
                       const double E1,
                       const double E2,
                       const double Elnsig2) {

  const int n = x.size();
  const int p = phi.size();
  const int q = alpha.size();
  const int l = std::max(p, q);
  const int nl = n + l;
  arma::vec lnsig2 = arma::vec(nl).fill(Elnsig2);
  double eta_i;
  double zeta_i;
  arma::vec series1 = arma::vec(nl).fill(0);
  arma::vec series2 = arma::vec(nl).fill(0);

  const arma::vec x_cent = x - mu;
  arma::vec sigt = arma::zeros<arma::vec>(n);


  const arma::rowvec alpha_r = arma::reverse(alpha);
  const arma::rowvec phi_r = arma::reverse(phi);
  const arma::rowvec beta_r = arma::reverse(beta);

  int il;
  int il1;
  int ilp;
  int ilq;

  for (int i = 0; i < n; ++i) {

    il = i + l;
    il1 = il - 1;
    ilp = il - p;
    ilq = il - q;
    lnsig2(il) = omega +
      arma::as_scalar(alpha_r * series1.subvec(ilq, il1)) +
      arma::as_scalar(beta_r * series2.subvec(ilq, il1)) +
      arma::as_scalar(phi_r * lnsig2.subvec(ilp, il1));

    sigt(i) = std::exp(lnsig2(il) / 2.0);
    eta_i = x_cent(i) / sigt(i);
    zeta_i = signCpp(eta_i) * std::log(std::abs(eta_i) + 1) - E1;
    series1(il) = zeta_i;
    series2(il) = std::abs(zeta_i) - E2;

  }

  return sigt;

}

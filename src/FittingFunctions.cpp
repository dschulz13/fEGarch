#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
int signCpp(const double val) {
    return (0.0 < val) - (val < 0.0);
}

// [[Rcpp::export]]
arma::vec sigt_fiegarchCpp(const arma::vec& x, const double mu,
                       const arma::rowvec& alpha,
                       const arma::rowvec& beta,
                       const double E1,
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

  for (int i = 0; i < n; ++i) {

    il = i + 1;

    lnsig2(il) = Elnsig2 +
      arma::as_scalar(alpha_r.subvec(nm1 - i, nm1) * series_alpha.subvec(0, i)) +
      arma::as_scalar(beta_r.subvec(nm1 - i, nm1) * series_beta.subvec(0, i));

    sigt(i) = std::exp(lnsig2(il) / 2.0);
    eta_i = x_cent(i) / sigt(i);
    series_alpha(il) = eta_i;
    series_beta(il) = std::abs(eta_i) - E1;

  }

  return sigt;

}


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

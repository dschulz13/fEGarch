#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Author: Dominik Schulz
// Version: January 16, 2025

// [[Rcpp::export]]
double binomial_coef(const int n, const int k) {
  if (k == 0 || k == n) {
    return 1;
  }
  double out = double(n) / k;
  for (int K = 1; K <= (k - 1); ++K) {
    out *= double(n - K) / (k - K);
  }
  return out;
}

// [[Rcpp::export]]
arma::vec pdf_sm(const arma::vec& x, const double mu, const arma::vec& sigt, const int P) {
  const double a = std::sqrt(2.0 * (P + 1.0));
  const double K = std::pow(2.0, -2.0 * P) * binomial_coef(2 * P, P);

  arma::vec cj = arma::ones<arma::vec>(P + 1);
  for (int j = 1; j < P; ++j) {
    cj(j + 1) = 2.0 * (P - j) / ((j + 1.0) * (2 * P - j)) * cj(j);
  }

  const double c1 = a * K / 2.0;

  const arma::vec one_div_sigt = 1.0 / sigt;

  const arma::vec x_abs = arma::abs((x - mu) % one_div_sigt);

  const arma::vec c2_vec = arma::exp(-a * x_abs);

  const int xn = x.size();

  arma::mat x_mat(xn, P + 1);

  for (int i = 0; i <= P; ++i) {
    x_mat.col(i) = cj(i) * std::pow(a, i) * arma::pow(x_abs, i);
  }

  const arma::vec c3_vec = arma::sum(x_mat, 1); // sum of each row;

  return c1 * (c2_vec % c3_vec);

}

// [[Rcpp::export]]
arma::vec pdf_sm_v1(const arma::vec& x, const int P) {
  const double a = std::sqrt(2.0 * (P + 1.0));
  const double K = std::pow(2.0, -2.0 * P) * binomial_coef(2 * P, P);

  arma::vec cj = arma::ones<arma::vec>(P + 1);
  for (int j = 1; j < P; ++j) {
    cj(j + 1) = 2.0 * (P - j) / ((j + 1.0) * (2 * P - j)) * cj(j);
  }

  const double c1 = a * K / 2.0;

  const arma::vec x_abs = arma::abs(x);

  const arma::vec c2_vec = arma::exp(-a * x_abs);

  const int xn = x.size();

  arma::mat x_mat(xn, P + 1);

  for (int i = 0; i <= P; ++i) {
    x_mat.col(i) = cj(i) * std::pow(a, i) * arma::pow(x_abs, i);
  }

  const arma::vec c3_vec = arma::sum(x_mat, 1); // sum of each row;

  return c1 * (c2_vec % c3_vec);

}

// [[Rcpp::export]]
arma::vec pdf_skew_ssm(const arma::vec& x, const int P, const double skew) {
  const double inv_skew = 1.0 / skew;
  const arma::uvec Hnz = (-x >= 0);
  const arma::uvec Hz = 1 - Hnz;
  const double f1 = 2.0 / (skew + inv_skew);
  const arma::vec f2 = Hnz % pdf_sm_v1(skew * x, P);
  const arma::vec f3 = Hz % pdf_sm_v1(inv_skew * x, P);
  return f1 * (f2 + f3);
}

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Author: Dominik Schulz
// Version: January 16, 2025

// [[Rcpp::export]]
int factorial(const int k) {
  if (k == 0 || k == 1) {
    return 1;
  }
  int out = 1;
  for (int i = 2; i <= k; ++i) {
    out *= i;
  }
  return out;
}

// [[Rcpp::export]]
int binomial_coef(const int n, const int k) {
  int out = factorial(n) / (factorial(k) * factorial(n - k));
}

// [[Rcpp::export]]
arma::vec pdf_sm(const arma::vec& x, const double mu, const arma::vec& sigt, const double P, const double b) {
  const double a = (P + 1.0) / b;
  const int K = std::pow(2.0, -2.0 * P) * binomial_coef(2 * P, P);

  arma::vec cj = arma::ones<arma::vec>(P + 1);
  for (int j = 1; j < P; ++j) {
    cj(j + 1) = 2.0 * (P - j) / ((j + 1.0) * (2 * P - j));
  }

  const double c1 = a * K / 2.0;

  const arma::vec x_abs = arma::abs(x);

  const arma::vec c2 = arma::exp(-a * x_abs);

  const int xn = x.size();

  const arma::mat x_mat = arma::mat(xn, P + 1);

  for (int i = 0; i <= P; ++i) {
    x_mat.col(i) = cj(i) * std::pow(a, i) * arma::pow(x_abs, i);
  }

  const arma::vec c3_vec = arma::sum(x_mat, 1); // sum of each row;

  return c1 * c2 * c3_vec;

}

// [[Rcpp::export]]
arma::vec pdf_norm_v1(const arma::vec& rt) {
  const arma::vec eta_t2 = arma::pow(rt, 2.0);
  const double f1 = std::sqrt(2 * M_PI);
  const arma::vec f2 = -0.5 * eta_t2;
  return arma::exp(f2) / f1;
}

// [[Rcpp::export]]
arma::vec pdf_std(const arma::vec& rt, const double mu, const arma::vec& sigt,
                  const double df) {
  const arma::vec one_div_sigt = 1.0 / sigt;
  const arma::vec eta_t2 = arma::pow((rt - mu) % one_div_sigt, 2.0);
  const double dfp1div2 = (df + 1.0) / 2.0;
  const double dfm2 = df - 2.0;
  const double f1 = std::tgamma(dfp1div2);
  const arma::vec f2 = one_div_sigt / (std::tgamma(df / 2.0) * std::sqrt(M_PI * dfm2));
  const arma::vec f3 = arma::pow(1.0 + eta_t2 / (dfm2), -dfp1div2);
  return f1 * (f2 % f3);
}


// [[Rcpp::export]]
arma::vec pdf_std_v1(const arma::vec& rt, const double df) {
  const arma::vec eta_t2 = arma::pow(rt, 2.0);
  const double dfp1div2 = (df + 1.0) / 2.0;
  const double dfm2 = df - 2.0;
  const double f1 = std::tgamma(dfp1div2);
  const double f2 = std::tgamma(df / 2.0) * std::sqrt(M_PI * dfm2);
  const arma::vec f3 = arma::pow(1.0 + eta_t2 / (dfm2), -dfp1div2);
  return f1 * (f3 / f2);
}

// [[Rcpp::export]]
arma::vec pdf_ged(const arma::vec& rt, const double mu, const arma::vec& sigt,
                  const double shape) {
  const arma::vec one_div_sigt = 1.0 / sigt;
  const arma::vec abs_eta_t = arma::abs((rt - mu) % one_div_sigt);
  const double g3b = std::tgamma(3.0 / shape);
  const double g1b = std::tgamma(1.0 / shape);
  const arma::vec f1 = one_div_sigt * std::sqrt(g3b / std::pow(g1b, 3.0));
  const double f2 = shape / 2.0;
  const double f3s = std::pow(g3b / g1b, f2);
  const arma::vec f3 = -arma::pow(abs_eta_t, shape) * f3s;
  return f2 * (f1 % arma::exp(f3));
}


// [[Rcpp::export]]
arma::vec pdf_ged_v1(const arma::vec& rt, const double shape) {
  const arma::vec abs_eta_t = arma::abs(rt);
  const double g3b = std::tgamma(3.0 / shape);
  const double g1b = std::tgamma(1.0 / shape);
  const double f1 = std::sqrt(g3b / std::pow(g1b, 3.0));
  const double f2 = shape / 2.0;
  const double f3s = std::pow(g3b / g1b, f2);
  const arma::vec f3 = -arma::pow(abs_eta_t, shape) * f3s;
  return f2 * (f1 * arma::exp(f3));
}

// [[Rcpp::export]]
arma::vec pdf_skew_snorm(const arma::vec& x, const double skew) {
  const double inv_skew = 1.0 / skew;
  const arma::uvec Hnz = (-x >= 0);
  const arma::uvec Hz = 1 - Hnz;
  const double f1 = 2.0 / (skew + inv_skew);
  const arma::vec f2 = Hnz % pdf_norm_v1(skew * x);
  const arma::vec f3 = Hz % pdf_norm_v1(inv_skew * x);
  return f1 * (f2 + f3);
}

// [[Rcpp::export]]
arma::vec pdf_skew_sstd(const arma::vec& x, const double df, const double skew) {
  const double inv_skew = 1.0 / skew;
  const arma::uvec Hnz = (-x >= 0);
  const arma::uvec Hz = 1 - Hnz;
  const double f1 = 2.0 / (skew + inv_skew);
  const arma::vec f2 = Hnz % pdf_std_v1(skew * x, df);
  const arma::vec f3 = Hz % pdf_std_v1(inv_skew * x, df);
  return f1 * (f2 + f3);
}

// [[Rcpp::export]]
arma::vec pdf_skew_sged(const arma::vec& x, const double shape, const double skew) {
  const double inv_skew = 1.0 / skew;
  const arma::uvec Hnz = (-x >= 0);
  const arma::uvec Hz = 1 - Hnz;
  const double f1 = 2.0 / (skew + inv_skew);
  const arma::vec f2 = Hnz % pdf_ged_v1(skew * x, shape);
  const arma::vec f3 = Hz % pdf_ged_v1(inv_skew * x, shape);
  return f1 * (f2 + f3);
}

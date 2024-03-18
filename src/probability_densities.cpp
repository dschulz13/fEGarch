#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector pdf_norm_cpp(NumericVector x) {

  return 1.0 / (std::sqrt(2.0 * M_PI)) * exp(-0.5 * Rcpp::pow(x, 2.0));

}

// [[Rcpp::export]]
NumericVector pdf_std_cpp(NumericVector x, const double df) {
  const double sd = std::sqrt(df / (df - 2));

  return sd * std::tgamma((df + 1.0) / 2.0) / (std::sqrt(df * M_PI) * std::tgamma(df / 2.0)) * Rcpp::pow(1 + Rcpp::pow(x * sd, 2.0) / df, -(df + 1.0) / 2.0);

}

// [[Rcpp::export]]
NumericVector pdf_ged_cpp(NumericVector x, const double shape) {
  const double c1 = std::tgamma(1 / shape);
  const double scale = std::sqrt(c1 / std::tgamma(3 / shape)) * std::pow(2.0, -1.0 / shape);
  return shape / (std::pow(2.0, 1.0 + 1.0 / shape) * scale * c1) * exp(-0.5 * Rcpp::pow(abs(x) / scale, shape));
}

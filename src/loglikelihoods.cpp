#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
int sign_cpp(const double x) {
  int out;
  if (x >= 0) {
    out = 1;
  } else {
    out = -1;
  }
  return out;
}

// [[Rcpp::export]]
arma::vec obtain_volas(arma::rowvec& parameters, const arma::vec& obs, const int p, const int q, const int include_mean, const double E1, const double E2) {

  double mean;
  arma::rowvec kappa;
  arma::rowvec gamma;
  arma::rowvec beta;

  const int n = obs.size();

  int j = 0;

  if (include_mean == 1) {
    mean = parameters(j);
    j = j + 1;
  } else {
    mean = 0;
  }

  const double omega = parameters(j);
  j = j + 1;

  if (p > 0) {
    kappa = parameters.subvec(j, j + p - 1);
    j = j + p;
    gamma = parameters.subvec(j, j + p - 1);
    j = j + p;
  } else {
    kappa = arma::zeros<arma::rowvec>(1);
    gamma = arma::zeros<arma::rowvec>(1);
  }

  if (q > 0) {
    beta = parameters.subvec(j, j + q - 1);
    j = j + q;
  } else {
    beta = arma::zeros<arma::rowvec>(1);
  }

  const double power1 = parameters(j);
  const double power2 = parameters(j + 1);

  const double Elnsig2 = omega / (1.0 - arma::accu(beta)); // Theoretical expectation of ln(sig_t^2)

  const int ps = kappa.size();
  const int qs = beta.size();

  const int l = std::max(ps, qs);

  arma::vec lnsig2(l + n, arma::fill::value(Elnsig2));
  arma::vec innovs_c_1 = arma::zeros<arma::vec>(l + n);
  arma::vec innovs_c_2 = arma::zeros<arma::vec>(l + n);

  const arma::rowvec kappa_rev = arma::reverse(kappa);
  const arma::rowvec gamma_rev = arma::reverse(gamma);
  const arma::rowvec beta_rev = arma::reverse(beta);
  arma::vec sigma(n);
  double innov_i;
  double term1;
  double term2;
  double term3;

  for (int i = l; i < n + l; ++i) {
    term1 = arma::as_scalar(kappa_rev * innovs_c_1.subvec(i - p, i - 1));
    term2 = arma::as_scalar(gamma_rev * innovs_c_2.subvec(i - p, i - 1));
    term3 = arma::as_scalar(beta_rev * lnsig2.subvec(i - q, i - 1));
    lnsig2(i) = omega + term1 / power1 +
      term2 / power2 +
      term3;
    sigma(i - l) = exp(lnsig2(i) / 2);
    innov_i = (obs(i - l) - mean) / sigma(i);
    innovs_c_1(i) = sign_cpp(innov_i) * std::pow(std::abs(innov_i), power1) - E1;
    innovs_c_2(i) = std::pow(std::abs(innov_i), power2) - E2;
  }

  return sigma;

}

// [[Rcpp::export]]
double llhNORM(const arma::vec& parameters, const arma::vec& obs, const int p, const int q, const int include_mean, const double E1, const double E2) {

  double mean = 0;
  if (include_mean == 1) {
    mean = parameters(0);
  }

  const int n = obs.size();

  const arma::vec volas = obtain_volas(parameters, obs, include_mean, p, q, E1, E2);

  const arma::vec inv_volas = arma::pow(volas, -1.0);

  const double c1 = (n / -2.0) * log(2 * M_PI);
  const double c2 = -1.0 * arma::accu(arma::log(volas));
  const double c3 = -1.0 * arma::accu((obs - mean) % inv_volas);

  return c1 + c2 + c3;

}

// [[Rcpp::export]]
double llhSTD(const arma::vec& parameters, const arma::vec& obs, const int p, const int q, const int include_mean, const double E1, const double E2) {

  double mean;
  double df;
  if (include_mean == 1) {
    mean = parameters(0);
    df = parameters(4 + 2 * p + q);
  } else {
    mean = 0;
    df = parameters(3 + 2 * p + q);
  }

  const int n = obs.size();

  const arma::vec volas = obtain_volas(parameters, obs, include_mean, p, q, E1, E2);
  const arma::vec inv_volas = arma::pow(volas, -1.0);

  const double c1 = n * log(std::tgamma((df + 1.0) / 2.0) / (std::sqrt(M_PI * (df - 2.0)) * std::tgamma(df / 2.0)));
  const double c2 = -1.0 * arma::accu(arma::log(volas));
  const double c3 = -1.0 * ((df + 1.0) / 2.0) * arma::accu(arma::log(1 + arma::pow((obs - mean) % inv_volas, 2.0) / (df - 2.0)));

  return c1 + c2 + c3;

}

// [[Rcpp::export]]
double llhGED(const arma::vec& parameters, const arma::vec& obs, const int p, const int q, const int include_mean, const double E1, const double E2) {

  double mean;
  double shape;
  if (include_mean == 1) {
    mean = parameters(0);
    shape = parameters(4 + 2 * p + q);
  } else {
    mean = 0;
    shape = parameters(3 + 2 * p + q);
  }

  const int n = obs.size();

  const double g1 = tgamma(3.0 / shape);
  const double g2 = tgamma(1.0 / shape);

  const arma::vec volas = obtain_volas(parameters, obs, include_mean, p, q, E1, E2);
  const arma::vec inv_volas = arma::pow(volas, -1.0);

  const double c1 = n * (log(shape / 2.0) + 0.5 * log(g1) - 1.5 * log(g2));
  const double c2 = -1.0 * arma::accu(arma::log(volas));
  const double c3 = -1.0 * std::pow(g1 / g2, shape / 2.0) * arma::accu(arma::pow(arma::abs((obs - mean) % inv_volas), shape));

  return c1 + c2 + c3;

}


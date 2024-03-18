#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec obtain_volas(arma::vec& parameters, const arma::vec& obs, const int p, const int q, const int include_mean, const double E1, const double E2) {

  double mean;
  arma::vec kappa;
  arma::vec gamma;
  arma::vec beta;

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
    kappa = arma::zeros<arma::vec>(1);
    gamma = arma::zeros<arma::vec>(1);
  }

  if (q > 0) {
    beta = parameters.subvec(j, j + q - 1);
    j = j + q;
  } else {
    beta = arma::zeros<arma::vec>(1);
  }

  const double power1 = parameters(j);
  const double power2 = parameters(j + 1);

  const double Elnsig2 = omega / (1.0 - arma::accu(beta)); // Theoretical expectation of ln(sig_t^2)

  const int ps = kappa.size();
  const int qs = beta.size();

  const int l = std::max(ps, qs);

  arma::vec lnsig2(l + n, arma::fill::value(Elnsig2));
  arma::vec innovs_c = arma::zeros<arma::vec>(l + n);

  for (int i = l; i < n + l; ++i) {
    lnsig2(i) = omega + kappa //continue here;
  }

}

// [[Rcpp::export]]
double llhNORM(const arma::vec& parameters, const arma::vec& obs, const int p, const int q, const int include_mean, const double E1, const double E2) {

  double mean = 0;
  if (include_mean == 1) {
    mean = parameters(0);
  }

  const int n = obs.size();

  const arma::vec volas = obtain_volas(parameters, obs, include_mean, p, q, E1, E2, E_c);

  const arma::vec inv_volas = arma::pow(volas, -1.0);

  const double c1 = (n / -2.0) * log(2 * M_PI);
  const double c2 = -1.0 * arma::accu(arma::log(volas));
  const double c3 = -1.0 * arma::accu((obs - mean) % inv_volas);

  return c1 + c2 + c3;

}

// [[Rcpp::export]]
double llhSTD(const arma::vec& parameters, const arma::vec& obs, const int p, const int q, const int include_mean, const double E1, const double E2) {

  double mean = 0;
  if (include_mean == 1) {
    mean = parameters(0);
  }

  const int n = obs.size();

  const arma::vec volas = obtain_volas(parameters, obs, include_mean, p, q, E1, E2, E_c);
  const arma::vec inv_volas = arma::pow(volas, -1.0);

  const double c1 = n * log(std::tgamma((df + 1.0) / 2.0) / (std::sqrt(M_PI * (df - 2.0)) * std::tgamma(df / 2.0)));
  const double c2 = -1.0 * arma::accu(arma::log(volas));
  const double c3 = -1.0 * ((df + 1.0) / 2.0) * arma::accu(arma::log(1 + arma::pow((obs - mean) % inv_volas, 2.0) / (df - 2.0)));

  return c1 + c2 + c3;

}

// [[Rcpp::export]]
double llhGED(const arma::vec& parameters, const arma::vec& obs, const int p, const int q, const int include_mean, const double E1, const double E2) {

  double mean = 0;
  if (include_mean == 1) {
    mean = parameters(0);
  }

  const int n = obs.size();

  const double g1 = tgamma(3.0 / shape);
  const double g2 = tgamma(1.0 / shape);

  const arma::vec volas = obtain_volas(parameters, obs, include_mean, p, q, E1, E2, E_c);
  const arma::vec inv_volas = arma::pow(volas, -1.0);

  const double c1 = n * (log(shape / 2.0) + 0.5 * log(g1) - 1.5 * log(g2));
  const double c2 = -1.0 * arma::accu(arma::log(volas));
  const double c3 = -1.0 * std::pow(g1 / g2, shape / 2.0) * arma::accu(arma::pow(arma::abs((obs - mean) % inv_volas), shape));

  return c1 + c2 + c3;

}


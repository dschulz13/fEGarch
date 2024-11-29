#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
int signCpp(const double val) {
    return (0.0 < val) - (val < 0.0);
}

// [[Rcpp::export]]
arma::vec sigt_egarch_shortCpp(const arma::vec& x, const double mu,
                               const double omega,
                               const arma::rowvec& phi,
                               const arma::rowvec& alpha,
                               const arma::rowvec& beta,
                               const double E_asy,
                               const double E_mag,
                               const double Elnsig2,
                               const arma::vec& powers,
                               const arma::vec& modulus,
                               const std::string& mode) {

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

  int il, il1, ilp, ilq;

  const double mod1 = modulus(0);
  const double mod2 = modulus(1);
  const double pow1 = powers.size() > 0 ? powers(0) : 1.0; // Default to 1.0 for safety
  const double pow2 = powers.size() > 1 ? powers(1) : 1.0;

  std::function<double(double)> compute_abs_eta_pow1;
  std::function<double(double)> compute_abs_eta_pow2;

  if (mode == "pow_pow") {
    compute_abs_eta_pow1 = [pow1, mod1, E_asy](double eta_i) {
      return (signCpp(eta_i) * std::pow(std::abs(eta_i) + mod1, pow1) - E_asy) / pow1;
    };
    compute_abs_eta_pow2 = [pow2, mod2, E_mag](double eta_i) {
      return (std::pow(std::abs(eta_i) + mod2, pow2) - E_mag) / pow2;
    };
  } else if (mode == "log_log") {
    compute_abs_eta_pow1 = [mod1, E_asy](double eta_i) {
      return signCpp(eta_i) * std::log(std::abs(eta_i) + mod1) - E_asy;
    };
    compute_abs_eta_pow2 = [mod2, E_mag](double eta_i) {
      return std::log(std::abs(eta_i) + mod2) - E_mag;
    };
  } else if (mode == "pow_log") {
    compute_abs_eta_pow1 = [pow1, mod1, E_asy](double eta_i) {
      return (signCpp(eta_i) * std::pow(std::abs(eta_i) + mod1, pow1) - E_asy) / pow1;
    };
    compute_abs_eta_pow2 = [mod2, E_mag](double eta_i) {
      return std::log(std::abs(eta_i) + mod2) - E_mag;
    };
  } else if (mode == "log_pow") {
    compute_abs_eta_pow1 = [mod1, E_asy](double eta_i) {
      return signCpp(eta_i) * std::log(std::abs(eta_i) + mod1) - E_asy;
    };
    compute_abs_eta_pow2 = [pow2, mod2, E_mag](double eta_i) {
      return (std::pow(std::abs(eta_i) + mod2, pow2) - E_mag) / pow2;
    };
  } else {
    Rcpp::stop("Invalid mode specified");
  }

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

    series1(il) = compute_abs_eta_pow1(eta_i);
    series2(il) = compute_abs_eta_pow2(eta_i);
  }

  return sigt;
}


// [[Rcpp::export]]
arma::vec sigt_egarch_longCpp(const arma::vec& x, const double mu,
                       const arma::rowvec& coef_inf,
                       const double kappa,
                       const double gamma,
                       const double E_asy,
                       const double E_mag,
                       const double Elnsig2,
                       const arma::vec& powers,
                       const arma::vec& modulus,
                       const std::string& mode) {

  const int n = x.size();
  const int n1 = n + 1;
  arma::vec lnsig2 = arma::vec(n1).fill(0);
  arma::vec series_g = arma::vec(n1).fill(0);

  const arma::vec x_cent = x - mu;
  arma::vec sigt = arma::zeros<arma::vec>(n);

  const arma::rowvec coef_r = arma::reverse(coef_inf);
  double eta_i;

  int il;
  int nm1 = n - 1;
  const double mod1 = modulus(0);
  const double mod2 = modulus(1);
  const double pow1 = powers.size() > 0 ? powers(0) : 1.0; // Default to 1.0 for safety
  const double pow2 = powers.size() > 1 ? powers(1) : 1.0;

  std::function<double(double)> compute_abs_eta_pow1;
  std::function<double(double)> compute_abs_eta_pow2;

  if (mode == "pow_pow") {
    compute_abs_eta_pow1 = [pow1, mod1, E_asy](double eta_i) {
      return (signCpp(eta_i) * std::pow(std::abs(eta_i) + mod1, pow1) - E_asy) / pow1;
    };
    compute_abs_eta_pow2 = [pow2, mod2, E_mag](double eta_i) {
      return (std::pow(std::abs(eta_i) + mod2, pow2) - E_mag) / pow2;
    };
  } else if (mode == "log_log") {
    compute_abs_eta_pow1 = [mod1, E_asy](double eta_i) {
      return signCpp(eta_i) * std::log(std::abs(eta_i) + mod1) - E_asy;
    };
    compute_abs_eta_pow2 = [mod2, E_mag](double eta_i) {
      return std::log(std::abs(eta_i) + mod2) - E_mag;
    };
  } else if (mode == "pow_log") {
    compute_abs_eta_pow1 = [pow1, mod1, E_asy](double eta_i) {
      return (signCpp(eta_i) * std::pow(std::abs(eta_i) + mod1, pow1) - E_asy) / pow1;
    };
    compute_abs_eta_pow2 = [mod2, E_mag](double eta_i) {
      return std::log(std::abs(eta_i) + mod2) - E_mag;
    };
  } else if (mode == "log_pow") {
    compute_abs_eta_pow1 = [mod1, E_asy](double eta_i) {
      return signCpp(eta_i) * std::log(std::abs(eta_i) + mod1) - E_asy;
    };
    compute_abs_eta_pow2 = [pow2, mod2, E_mag](double eta_i) {
      return (std::pow(std::abs(eta_i) + mod2, pow2) - E_mag) / pow2;
    };
  } else {
    Rcpp::stop("Invalid mode specified");
  }

  for (int i = 0; i < n; ++i) {

    il = i + 1;

    lnsig2(il) = Elnsig2 +
      arma::as_scalar(coef_r.subvec(nm1 - i, nm1) * series_g.subvec(0, i));

    sigt(i) = std::exp(lnsig2(il) / 2.0);
    eta_i = x_cent(i) / sigt(i);

    series_g(il) = kappa * compute_abs_eta_pow1(eta_i) + gamma * compute_abs_eta_pow2(eta_i);

  }

  return sigt;

}

// [[Rcpp::export]]
arma::vec sigt_loggarch_long(const arma::vec& x, const double mu,
                       const arma::rowvec& coef_inf,
                       const double Elneta2,
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
    epst(il) = std::log(std::pow(eta_i, 2.0)) - Elneta2;

  }

  return sigt;

}

// [[Rcpp::export]]
arma::vec sigt_loggarch_short(const arma::vec& x, const double mu,
                       const double omega,
                       const arma::rowvec& phi,
                       const arma::rowvec& psi,
                       const double Elneta2,
                       const double Elnsig2) {

  const int n = x.size();
  const int p = phi.size();
  const int q = psi.size();
  const int l = std::max(p, q);
  const int nl = n + l;
  arma::vec lnsig2 = arma::vec(nl).fill(Elnsig2);
  double eta_i;
  arma::vec series_innov = arma::vec(nl).fill(0);

  const arma::vec x_cent = x - mu;
  arma::vec sigt = arma::zeros<arma::vec>(n);

  const arma::rowvec psi_r = arma::reverse(psi);
  const arma::rowvec phi_r = arma::reverse(phi);

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
      arma::as_scalar(psi_r * series_innov.subvec(ilq, il1)) +
      arma::as_scalar(phi_r * lnsig2.subvec(ilp, il1));

    sigt(i) = std::exp(lnsig2(il) / 2.0);
    eta_i = x_cent(i) / sigt(i);
    series_innov(il) = std::log(std::pow(eta_i, 2.0)) - Elneta2;

  }

  return sigt;

}

// [[Rcpp::export]]
arma::vec sigt_delta_fiaparch_Cpp(const arma::vec& x, const double mu,
                       const arma::rowvec& coef_inf,
                       const double gamma,
                       const double delta,
                       const double omega_sig) {

  const int n = x.size();
  const int l = n + 1000;
  arma::vec sig_delta = arma::vec(n).fill(0);

  const arma::vec x_cent = x - mu;
  const arma::vec transf = arma::pow(arma::abs(x_cent) - gamma * x_cent, delta);
  //Initialization based on sample
  //const double E_g_delta = arma::mean(transf);
  const double E_g_delta = 0.02;
  arma::vec g_fun = arma::vec(l).fill(E_g_delta);
  g_fun.subvec(1000, l - 1) = transf;

  const arma::rowvec coef_r = arma::reverse(coef_inf);

  for (int i = 0; i < n; ++i) {

    sig_delta(i) = omega_sig +
      arma::as_scalar(coef_r.subvec(l - 1000 - i, l - 1) * g_fun.subvec(0, 999 + i));

  }

  return sig_delta;

}

// [[Rcpp::export]]
arma::vec sigt_lmgarch_Cpp(const arma::vec& x, const double mu,
                       const arma::rowvec& coef_inf,
                       const double E_e2) {

  const int n = x.size();
  const int l = n + 1000;
  arma::vec sig = arma::vec(n).fill(0);

  const arma::vec x_cent = x - mu;

  //Initialization based on sample
  arma::vec e2 = arma::vec(l).fill(E_e2);
  e2.subvec(1000, l - 1) = arma::pow(x_cent, 2.0);

  const arma::rowvec coef_r = arma::reverse(coef_inf);

  for (int i = 0; i < n; ++i) {

    sig(i) = std::sqrt(arma::as_scalar(coef_r.subvec(l - 1000 - i, l - 1) * e2.subvec(0, 999 + i)));

  }

  return sig;

}


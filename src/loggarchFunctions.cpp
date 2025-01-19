#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


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
Rcpp::List loggarch_shortmemo_type_sim(const arma::vec& innov, const double omega_sig,
                          const arma::rowvec& ar_coef, const arma::rowvec& ma_coef,
                          const int n_out,
                          const double mu,
                          const double Econst) {

  const int n = innov.size();

  const int n_ar = ar_coef.size();
  const int n_ma = ma_coef.size();
  const int l = std::max(n_ar, n_ma);

  const int n_pre = l;
  const int n_all = n + n_pre;

  arma::vec lnsig2(n_all);
  lnsig2.fill(omega_sig);

  const arma::rowvec ar_rev = arma::reverse(ar_coef);
  const arma::rowvec ma_rev = arma::reverse(ma_coef);

  arma::vec xi = arma::zeros<arma::vec>(n_all);
  xi.subvec(l, n_all - 1) = arma::log(arma::pow(innov, 2.0)) - Econst;

  const double omega = omega_sig * (1 - arma::accu(ar_coef));
  const int lar = l - n_ar;
  const int lma = l - n_ma;
  const int lm1 = l - 1;
  int li;

  for (int i = 0; i < n; ++i) {

    li = lm1 + i;

    lnsig2(l + i) = omega +
      arma::as_scalar(ar_rev * lnsig2.subvec(lar + i, li)) +
      arma::as_scalar(ma_rev * xi.subvec(lma + i, li));

  }

  const arma::vec sigt = arma::exp(lnsig2.subvec(n_all - n_out, n_all - 1) / 2.0);
  const arma::vec innov_out = innov.subvec(n - n_out, n - 1);

  const arma::vec rt = mu + sigt % innov_out;

  return Rcpp::List::create(_["rt"] = rt, _["sigt"] = sigt, _["etat"] = innov_out);
}

// [[Rcpp::export]]
Rcpp::List loggarch_longmemo_type_sim(const arma::vec& innov, const double omega_sig,
                          const arma::rowvec& coef_inf,
                          const int n_out,
                          const double mu,
                          const double Econst) {

  const int n = innov.size();

  const int n_pre = 1;
  const int n_all = n + n_pre;

  arma::vec lnsig2(n);

  const arma::rowvec coef_rev = arma::reverse(coef_inf);

  arma::vec xi = arma::zeros<arma::vec>(n_all);
  xi.subvec(1, n_all - 1) = arma::log(arma::pow(innov, 2.0)) - Econst;

  const int n_coef_m1 = coef_inf.size() - 1;

  for (int i = 0; i < n; ++i) {

    lnsig2(i) = omega_sig +
      arma::as_scalar(coef_rev.subvec(n_coef_m1 - i, n_coef_m1) * xi.subvec(0, i));

  }

  const arma::vec sigt = arma::exp(lnsig2.subvec(n - n_out, n - 1) / 2.0);
  const arma::vec innov_out = innov.subvec(n - n_out, n - 1);

  const arma::vec rt = mu + sigt % innov_out;

  return Rcpp::List::create(_["rt"] = rt, _["sigt"] = sigt, _["etat"] = innov_out);
}

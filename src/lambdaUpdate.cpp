#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' @keywords internal
// [[Rcpp::export]]
arma::mat lambdaUpdate(const Rcpp::List& b_, const Rcpp::List& imat_,
                       const Rcpp::List& zt_, const Rcpp::List& pb_,
                       const Rcpp::List& v_, const arma::mat& gam,
                       const arma::vec& gam_vec, const int& q, const arma::vec& nev,
                       const Rcpp::List& h_) {

  // Updates of lambda0 (E-step and M-step)

  arma::vec haz = arma::zeros<arma::vec>(nev.n_elem);

  // loop over subjects
  for (int i=0; i<b_.size(); i++) {

    // extract matrices from lists for subject i
    arma::mat b = Rcpp::as<arma::mat>(b_[i]);
    arma::mat I = Rcpp::as<arma::mat>(imat_[i]);
    arma::mat zt = Rcpp::as<arma::mat>(zt_[i]);
    arma::vec pb = Rcpp::as<arma::vec>(pb_[i]);
    arma::vec v = Rcpp::as<arma::vec>(v_[i]);
    Rcpp::DataFrame h = Rcpp::as<Rcpp::DataFrame>(h_[i]);

    // subjects who are censored before the first failure time
    // do not contribute towards \lambda estimation
    int tj_ind = h["tj.ind"];
    if (tj_ind == 0) continue;

    arma::mat expW_new = exp((b * gam) * trans(I * zt));
    arma::mat EexpVstar = mean(expW_new.each_col() % pb, 0);
    if (q > 0) {
      EexpVstar *= arma::as_scalar(exp(v.t() * gam_vec.subvec(0, q-1)));
    }
    haz.subvec(0, EexpVstar.n_cols-1) += EexpVstar.t();

  }

  return(nev / haz);

}

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' Calculation of exp{W(tj, b)}
//'
//' @keywords internal
// [[Rcpp::export]]
List expWArma(Rcpp::List imat_, Rcpp::List zt_,
              Rcpp::List b_, arma::mat gam) {
  List expw(zt_.size());
  for (int i=0; i<zt_.size(); i++) {
    arma::mat I  = Rcpp::as<arma::mat>(imat_[i]);
    arma::mat zt = Rcpp::as<arma::mat>(zt_[i]);
    arma::mat b  = Rcpp::as<arma::mat>(b_[i]);
    expw[i] = exp(b * gam * zt * I.t());
  }
  return(expw);
}


//' Calculation of E[exp{W(tj, b)}]
//'
//' @keywords internal
// [[Rcpp::export]]
List EexpWArma(Rcpp::List w_, Rcpp::List pb_) {
  List out(w_.size());
  for (int i=0; i<w_.size(); i++) {
    arma::mat w = Rcpp::as<arma::mat>(w_[i]);
    arma::vec pb = Rcpp::as<arma::vec>(pb_[i]);
    out[i] = mean(w.each_col() % pb, 0);
  }
  return(out);
}

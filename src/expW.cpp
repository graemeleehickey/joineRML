#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//---------------------------------------------------------
// exp{W(tj, b)}
//---------------------------------------------------------

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

// expW <- expWArma(IW.fail, Zit.fail, bi.y, gamma.scale)

//---------------------------------------------------------
// E[exp{W(tj, b)}]
//---------------------------------------------------------

// [[Rcpp::export]]
List EexpWArma(Rcpp::List w_, Rcpp::List f_, 
               Rcpp::List d_) {
  List out(w_.size());
  for (int i=0; i<w_.size(); i++) {
    arma::mat w = Rcpp::as<arma::mat>(w_[i]);
    arma::vec f = Rcpp::as<arma::vec>(f_[i]);
    double d    = Rcpp::as<double>(d_[i]);
    out[i] = mean(w.each_col() % f, 0) / d;
  }
  return(out);
}

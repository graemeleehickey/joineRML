#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' Calculation of exp{W(tj, b)}
//'
//' @keywords internal
// [[Rcpp::export]]
List expWArma(Rcpp::List iz_, Rcpp::List b_, arma::mat gam, Rcpp::List h_) {

  List expw(b_.size());

  for (int i=0; i<b_.size(); i++) {

    Rcpp::DataFrame h = Rcpp::as<Rcpp::DataFrame>(h_[i]);
    arma::mat iz = Rcpp::as<arma::mat>(iz_[i]);
    arma::mat b = Rcpp::as<arma::mat>(b_[i]);
    int tj_ind = h["tj.ind"];

    expw[i] = exp(b * gam * iz);

    if (tj_ind == 0) {
      // subjects who are censored before first failure
      // time do not contribute anything
      expw[i] = arma::zeros(arma::size(as<arma::mat>(expw[i])));
    }

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

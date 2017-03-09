#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' @keywords internal
// [[Rcpp::export]]
List expWArma(const Rcpp::List& iz_, const Rcpp::List& b_, const arma::mat& gam,
              const Rcpp::List& h_) {

  // Calculation of exp{W(tj, b)}

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

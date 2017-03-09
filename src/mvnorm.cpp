#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' @keywords internal
// [[Rcpp::export]]
arma::mat mvrnormArma(const int& n, const arma::vec& mu, const arma::mat& sigma) {

  // Antithetic simulation of MVN random variables

  int ncols = sigma.n_cols;
  // +/- 1 vector (n each) for antithetic calculation
  arma::vec pmvec = arma::ones<arma::vec>(2*n);
  pmvec.tail(n) *= -1;
  arma::mat Z = repmat(arma::randn(n, ncols), 2, 1);

  return(arma::repmat(mu, 1, 2*n).t() +
         (Z.each_col() % pmvec) * arma::trimatu(arma::chol(sigma)));

}


//' @keywords internal
// [[Rcpp::export]]
List bSim(const int& n, const Rcpp::List& Mean_, const Rcpp::List& Sigma_) {

  // Get a list of MVN samples for each subject

  List b(Mean_.size());
  for(int i=0; i<Mean_.size(); i++) {
    arma::vec m = Rcpp::as<arma::vec>(Mean_[i]);
    arma::mat s = Rcpp::as<arma::mat>(Sigma_[i]);
    b[i] = mvrnormArma(n, m, s);
  }

  return(b);

}

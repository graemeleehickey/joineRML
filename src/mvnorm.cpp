#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' Antithetic simulation of MVN random variables
//'
//' @keywords internal
// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu,
                      arma::mat sigma) {
  int ncols = sigma.n_cols;
  // +/- 1 vector (n each) for antithetic calculation
  arma::vec pmvec = arma::ones<arma::vec>(2*n);
  pmvec.tail(n) *= -1;
  arma::mat Z = repmat(arma::randn(n, ncols), 2, 1);
  return(arma::repmat(mu, 1, 2*n).t() +
         (Z.each_col() % pmvec) * arma::trimatu(arma::chol(sigma)));
}


//' Get a list of MVN samples for each subject
//'
//' @keywords internal
// [[Rcpp::export]]
List bSim(int n, List Mean_, List Sigma_) {
  List b(Mean_.size());
  for(int i=0; i<Mean_.size(); i++) {
    arma::vec m = Rcpp::as<arma::vec>(Mean_[i]);
    arma::mat s = Rcpp::as<arma::mat>(Sigma_[i]);
    b[i] = mvrnormArma(n, m, s);
  }
  return(b);
}

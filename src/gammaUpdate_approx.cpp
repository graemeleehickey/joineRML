#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' Newton-Raphson updates of gamma (E-step and M-step) using an empirical
//' infomration matrix approximation
//'
//' @keywords internal
// [[Rcpp::export]]
List gammaUpdate_approx(Rcpp::List b_, Rcpp::List z_, Rcpp::List w_,
                        Rcpp::List pb_, arma::vec haz,
                        Rcpp::List v_, Rcpp::List h_, int K, int q, int nev) {

  // declare score and E[delta x v*]
  arma::mat Si = arma::zeros<arma::mat>(q+K, w_.size());
  arma::vec S = arma::zeros<arma::vec>(q+K);
  arma::mat Evstari = arma::zeros<arma::mat>(q+K, w_.size());
  arma::mat I = arma::zeros<arma::mat>(q+K, q+K);

  // loop over subjects
  for (int i=0; i<w_.size(); i++) {

    Rcpp::checkUserInterrupt();

    // extract matrices from lists for subject i
    arma::mat b = Rcpp::as<arma::mat>(b_[i]);
    arma::mat z = Rcpp::as<arma::mat>(z_[i]);
    arma::mat w = Rcpp::as<arma::mat>(w_[i]);
    arma::vec pb = Rcpp::as<arma::vec>(pb_[i]);
    arma::vec v = Rcpp::as<arma::vec>(v_[i]);
    Rcpp::DataFrame h = Rcpp::as<Rcpp::DataFrame>(h_[i]);

    // subjects who are censored before the first failure time
    // do not contribute towards \gamma estimation
    int tj_ind = h["tj.ind"];
    if (tj_ind == 0) continue;

    int nj = w.n_cols;  // number of failures upto time T_i
    int delta = h["delta"]; // delta_i

    arma::mat bzt = b * z; // b x t(z)
    arma::mat bztev = bzt % repmat(w, 1, K); // b x t(Z) . exp(v*gamma)
    arma::mat Eexpvj = (mean(w.each_col() % pb, 0)) % trans(haz.subvec(0, nj-1));
    arma::mat Eexpv = sum(Eexpvj, 1); // lambda0 x E[exp(v*gamma)]
    arma::mat hexpand = trans(repmat(haz.subvec(0, nj-1), K, 1)); // K reps of lambda0(tj)
    arma::mat outj = (mean(bztev.each_col() % pb, 0)) % hexpand;

    arma::mat Eb = mean(b.each_col() % pb, 0);

    // loop of K longitudinal outcomes
    for(int k=0; k<K; k++) {
      // E[delta x v*(T_i)]
      Evstari(q+k, i) = delta * arma::dot(z.col(nj*(k+1)-1), Eb);
      // score elements for K Zb's
      Si(q+k, i) = Evstari(q+k, i) - arma::as_scalar(sum(outj.cols(nj*k, nj*(k+1)-1), 1));

    } // end loop over outcomes k

    // elements for q V_i's (only when q > 0)
    if (q > 0) {
      // E[delta x v]
      Evstari.submat(0, i, q-1, i) = delta * v;
      // score elements
      Si.submat(0, i, q-1, i) = Evstari.submat(0, i, q-1, i) - v * Eexpv;
    }

    S += Si.col(i);
    I += arma::kron(Si.col(i), Si.col(i).t()); // see below

  } // end loop over subjects i


  // an approximate I-matrix using observed empirical information
  I = I - arma::kron(S, S.t()) / w_.size();

  return List::create(
    _["gDelta"]  = 0.5 * solve(I, S),
    _["scorei"] = Si
  );

}

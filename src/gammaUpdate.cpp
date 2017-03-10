#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' @keywords internal
// [[Rcpp::export]]
List gammaUpdate(const Rcpp::List& b_, const Rcpp::List& z_, const Rcpp::List& w_,
                 const Rcpp::List& pb_, const arma::vec& haz, const Rcpp::List& v_,
                 const Rcpp::List& h_, const int& K, const int& q, const int& nev,
                 const arma::vec& jcount) {

  // Newton-Raphson updates of gamma (E-step and M-step) using an exact observed
  // information calculation

  // declare score, E[delta x v*], and information matrix
  arma::mat Si = arma::zeros<arma::mat>(q+K, w_.size());
  arma::vec S = arma::zeros<arma::vec>(q+K);
  arma::mat Evstari = arma::zeros<arma::mat>(q+K, w_.size());
  arma::mat I = arma::zeros<arma::mat>(q+K, q+K);
  arma::mat Gammaj = arma::zeros<arma::mat>(q+K, nev);
  arma::mat Gamma = arma::zeros<arma::mat>(q+K, q+K);

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

    int nj = w.n_cols;  // number of failure times upto time T_i
    int delta = h["delta"]; // delta_i

    arma::mat Ii_int = arma::zeros<arma::mat>(q+K, q+K); // information matrix (uninitialized) for subject i
    arma::mat bzt = b * z; // b x t(z)
    arma::mat bztev = bzt % repmat(w, 1, K); // b x t(Z) . exp(v*gamma)
    arma::mat Eexpvj = (mean(w.each_col() % pb, 0)) % trans(haz.subvec(0, nj-1));
    arma::mat Eexpv = sum(Eexpvj, 1); // lambda0 x E[exp(v*gamma)]
    arma::mat hexpand = trans(repmat(haz.subvec(0, nj-1), K, 1)); // K reps of lambda0(tj)
    arma::mat outj = (mean(bztev.each_col() % pb, 0)) % hexpand;
    arma::mat bzt2ev = bzt % bztev; // [b x t(z)]^2 . exp(v*gamma)
    arma::mat Ii_int_Kdiag = (mean(bzt2ev.each_col() % pb, 0)) % hexpand;

    arma::mat Eb = mean(b.each_col() % pb, 0);

    // loop of K longitudinal outcomes
    for (int k=0; k<K; k++) {
      // E[delta x v*(T_i)]
      Evstari(q+k, i) = delta * arma::dot(z.col(nj*(k+1)-1), Eb);

      // score elements for K Zb's
      Si(q+k, i) = arma::as_scalar(sum(outj.cols(nj*k, nj*(k+1)-1), 1));

      // cross-prod (diagonal) elements for K Zb's only
      Ii_int(q+k, q+k) = arma::as_scalar(sum(Ii_int_Kdiag.cols(nj*k, nj*(k+1)-1), 1));

      // cross-prod (off-diagonal) elements for K Zb's only
      for (int k2=k+1; k2<K; k2++) {
        arma::mat bztcross = bztev.cols(nj*k, nj*(k+1)-1) % bzt.cols(nj*k2, nj*(k2+1)-1);
        Ii_int(q+k, k2+q) = arma::as_scalar(sum((mean(bztcross.each_col() % pb, 0)) % trans(haz.subvec(0, nj-1)), 1));
        Ii_int(k2+q, q+k) = Ii_int(q+k, k2+q);
      }

      // cross-prod elements for q V_i's and K Zb's
      if (q > 0) {
        for (int j=0; j<q; j++) {
           Ii_int(j, q+k) = Si(q+k, i) * v(j);
           Ii_int(q+k, j) = Ii_int(j, q+k);
        }
      }
    } // end loop over outcomes k

    // Gamma_j vectors (for cross-product later)
    Gammaj.submat(q, 0, q+K-1, nj-1) += trans(reshape(trans(outj), nj, K));

    // elements for q V_i's (only when q > 0)
    if (q > 0) {
      // E[delta x v]
      Evstari.submat(0, i, q-1, i) = delta * v;
      // score elements
      Si.submat(0, i, q-1, i) = v * Eexpv;
      // cross-prod elements
      Ii_int.submat(0, 0, q-1, q-1) = (v * trans(v)) * arma::as_scalar(Eexpv);
      // Gamma_j elements
      Gammaj.submat(0, 0, q-1, nj-1) += v * Eexpvj;
    }

    S += (Evstari.col(i) - Si.col(i)); // NB: actual score is sum(Evstari - Si)
    I += Ii_int;

  } // end loop over subjects i

  // lambda0 x Gamma_j sum term (minus from information matrix)
  for (int t=0; t<nev; t++) {
    //Rcpp::Rcout << jcount(t) << std::endl;
    Gamma += Gammaj.col(t) * trans(Gammaj.col(t)) / jcount(t);
  }

  return List::create(
    Named("gDelta")  = solve(I - Gamma, S),
    Named("scorei")  = Evstari - Si
  );

}


//' @keywords internal
// [[Rcpp::export]]
arma::mat hazHat(const Rcpp::List& w_, const Rcpp::List& pb_, const arma::vec& nev) {

  // lambda0(t) for profile score function of beta

  arma::vec haz = arma::zeros<arma::vec>(nev.n_elem);

  // loop over subjects
  for (int i=0; i<w_.size(); i++) {

    // extract matrices from lists for subject i
    arma::mat w = Rcpp::as<arma::mat>(w_[i]);
    arma::vec pb = Rcpp::as<arma::vec>(pb_[i]);

    haz.subvec(0, w.n_cols-1) += arma::trans(mean(w.each_col() % pb, 0));

  }

  return(nev / haz);

}




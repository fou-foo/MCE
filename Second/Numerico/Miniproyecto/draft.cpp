#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat cube_sum1(arma::cube c) {
  arma::mat ss(c.n_cols, c.n_slices, arma::fill::zeros);
  for(int i = 0; i < c.n_rows; i++) {
    ss += c.tube(arma::span(i), arma::span::all);
  }
  return ss;
}

// [[Rcpp::export]]
arma::mat cube_sum2(arma::cube c) {
  arma::mat ss(c.n_cols, c.n_slices, arma::fill::zeros);
  for(int i = 0; i < c.n_rows; i++) {
    for(int j = 0; j < c.n_cols; j++) {
      for(int k = 0; k < c.n_slices; k++) {
        ss(j,k) += c(i,j,k);
      }}}
  return ss;
}


// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

#include <RcppArmadillo.h>
// #include <RcppEigen.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
SEXP MatMult3(arma::mat A, arma::mat B, arma::mat C){
  arma::mat D = A * B * C;
  return Rcpp::wrap(D);
}

// [[Rcpp::export]]
SEXP MatMult2(arma::mat A, arma::mat B){
  arma::mat D = A * B;
  return Rcpp::wrap(D);
}

// [[Rcpp::export]]
SEXP NumMult2(arma::vec A, arma::vec B){
  arma::vec D = A % B;
  return Rcpp::wrap(D);
}

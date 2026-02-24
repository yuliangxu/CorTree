#ifndef POLYAGAMMA_H
#define POLYAGAMMA_H

#include <RcppArmadillo.h>

using namespace Rcpp;

// Function declaration
arma::vec arma_pgdraw(const arma::vec& b, const arma::vec& c);
NumericVector rcpp_pgdraw(NumericVector b, NumericVector c);

#endif // POLYAGAMMA_H
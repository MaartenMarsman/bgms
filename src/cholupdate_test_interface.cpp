// -----------------------------------------------------------------------------
// cholupdate_test_interface.cpp
//
// R-facing test entry for the rank-1 Cholesky downdate. Exposes the success
// flag so tests can check that a non-positive-definite downdate is reported
// to the caller (the model classes rebuild the factor when it is not).
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>

#include "math/cholupdate.h"

// -----------------------------------------------------------------------------
// test_cholesky_downdate:
//   Apply the rank-1 downdate R'R - uu' and return the resulting factor
//   together with the positive-definiteness flag.
// -----------------------------------------------------------------------------

// [[Rcpp::export(name = "test_cholesky_downdate")]]
Rcpp::List test_cholesky_downdate(arma::mat R, arma::vec u, double eps = 1e-12) {
    bool ok = cholesky_downdate(R, u, eps);
    return Rcpp::List::create(Rcpp::Named("R") = R, Rcpp::Named("ok") = ok);
}

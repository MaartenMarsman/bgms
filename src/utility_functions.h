#include <Rcpp.h>
using namespace Rcpp;

// ----------------------------------------------------------------------------|
// A c++ version of table
// ----------------------------------------------------------------------------|
Rcpp::IntegerVector table_cpp(Rcpp::IntegerVector x);

// ----------------------------------------------------------------------------|
// Remove row i and column i from a matrix
// ----------------------------------------------------------------------------|
Rcpp::NumericMatrix remove_row_col_matrix(Rcpp::NumericMatrix X, int i);

// ----------------------------------------------------------------------------|
// Add a row and column to a matrix (and fill with beta variables)
// ----------------------------------------------------------------------------|
Rcpp::NumericMatrix add_row_col_block_prob_matrix(Rcpp::NumericMatrix X,
                                                  double beta_alpha,
                                                  double beta_beta);
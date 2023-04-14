#include <Rcpp.h>
using namespace Rcpp;


// ----------------------------------------------------------------------------|
// A c++ version of table
// ----------------------------------------------------------------------------|
IntegerVector table_cpp(IntegerVector x) {
  int n = x.size();
  int m = max(x);
  IntegerVector counts(m + 1);

  for (int i = 0; i < n; i++) {
    counts[x[i]]++;
  }

  return counts;
}

// ----------------------------------------------------------------------------|
// Remove row i and column i from a matrix
// ----------------------------------------------------------------------------|
NumericMatrix remove_row_col_matrix(NumericMatrix X, int i) {
  int n = X.nrow();

  // Handle special case of matrix with only two rows and columns
  if (n == 2) {
    if(i == 0) {
      X = X(Range(1, 1), Range(1, 1));
    } else {
      X = X(Range(0, 0), Range(0, 0));
    }
    return X;
  }

  // Remove row i
  for (int j = i; j < n - 1; j++) {
    X(j, _) = X(j+1, _);  // Shift all rows below i up by one
  }
  X = X(Range(0, n - 2), _);  // Remove last row

  // Remove column i
  for (int j = i; j < n - 1; j++) {
    X(_, j) = X(_, j + 1);  // Shift all columns right of i to the left by one
  }
  X = X(_, Range(0, n - 2));  // Remove last column

  return X;
}

// ----------------------------------------------------------------------------|
// Add a row and column to a matrix (and fill with beta variables)
// ----------------------------------------------------------------------------|
NumericMatrix add_row_col_block_prob_matrix(NumericMatrix X,
                                            double beta_alpha,
                                            double beta_beta) {
  int dim = X.nrow();
  NumericMatrix Y(dim + 1, dim + 1);

  for(int r = 0; r < dim; r++) {
    for(int c = 0; c < dim; c++) {
      Y(r, c) = X(r, c);
    }
  }

  for(int i = 0; i < dim; i++) {
    Y(dim, i) = R::rbeta(beta_alpha, beta_beta);
    Y(i, dim) = Y(dim, i);
  }
  Y(dim, dim) = R::rbeta(beta_alpha, beta_beta);

  return Y;
}
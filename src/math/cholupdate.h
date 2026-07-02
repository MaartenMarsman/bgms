#pragma once

/**
 * @file cholupdate.h
 * @brief Rank-1 Cholesky update and downdate.
 *
 * Givens-rotation (update) and hyperbolic-rotation (downdate) algorithms
 * for modifying an upper-triangular Cholesky factor R in place so that
 *
 *   update:   R'R + uu' = R1'R1
 *   downdate: R'R - uu' = R1'R1
 *
 * The core routine `chol_up()` is ported from Simon Wood's mgcv package:
 *   https://github.com/cran/mgcv/blob/master/src/mat.c  (function chol_up)
 * See Golub & Van Loan (2013, 4th ed., Section 6.5.4) for the hyperbolic
 * rotation downdate algorithm.
 *
 * The two C++ wrappers below accept Armadillo matrices and forward to the
 * column-major C implementation. The sub-diagonal entries of the first two
 * columns of R are used as scratch storage during the computation but are
 * zeroed on return.
 */

#include <RcppArmadillo.h>

/**
 * Rank-1 Cholesky update: R'R + uu' = R1'R1.
 *
 * Modifies R in place using Givens rotations.
 *
 * @param R    Upper-triangular Cholesky factor (p x p, modified in place)
 * @param u    Update vector (length p, modified in place)
 * @param eps  Tolerance for near-singularity (default 1e-12)
 */
void cholesky_update(  arma::mat& R, arma::vec& u, double eps = 1e-12);

/**
 * Rank-1 Cholesky downdate: R'R - uu' = R1'R1.
 *
 * Modifies R in place using hyperbolic rotations.
 *
 * @param R    Upper-triangular Cholesky factor (p x p, modified in place)
 * @param u    Downdate vector (length p, modified in place)
 * @param eps  Tolerance for near-singularity (default 1e-12)
 * @return     true on success. false if the downdated matrix would not be
 *             positive definite; R is then left partially updated and the
 *             caller must rebuild the factor from the source matrix.
 */
bool cholesky_downdate(arma::mat& R, arma::vec& u, double eps = 1e-12);

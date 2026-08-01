#pragma once

/**
 * @file cholesky_helpers.h
 * @brief Shared algebraic helpers for Cholesky-based precision updates.
 *
 * Pure functions with no model-specific state.  Used by both GGMModel and
 * MixedMRFModel for proposal constant extraction and log-determinant
 * computation.
 */

#include <RcppArmadillo.h>
#include <array>
#include <cmath>
#include "math/explog_macros.h"

namespace cholesky_helpers {

/**
 * Log-determinant of a positive-definite matrix from its upper-triangular
 * Cholesky factor R (where Ω = R'R).
 *
 * @param R  Upper-triangular Cholesky factor.
 * @return   log|Ω| = 2 Σ log(R_ii).
 */
inline double get_log_det(const arma::mat& R) {
    return 2.0 * arma::accu(ARMA_MY_LOG(arma::vec(R.diag())));
}

/**
 * Schur complement element: A(ii,jj) − A(ii,i) A(jj,i) / A(i,i).
 *
 * Used to compute entries of the inverse of a submatrix from the full
 * covariance matrix.
 *
 * @param A   Symmetric positive-definite matrix.
 * @param i   Conditioning index.
 * @param ii  Row index of the desired element.
 * @param jj  Column index of the desired element.
 * @return    Schur complement entry.
 */
inline double compute_inv_submatrix_i(const arma::mat& A, size_t i,
                                      size_t ii, size_t jj) {
    return A(ii, jj) - A(ii, i) * A(jj, i) / A(i, i);
}

/**
 * Rank-2 matrix-determinant-lemma log-ratio log|K'| − log|K| for a precision
 * change confined to entries (i,j), (j,i), (j,j). cc11/cc12/cc22 are the 2×2
 * update-Gram-matrix entries I + Vᵀ Σ U.
 *
 * Shared by GGMModel::log_det_ratio_edge and MixedMRFModel::log_det_ratio_yy_edge,
 * which differ only in where the current precision entries come from (GGM stores
 * K directly; Mixed stores −½K), so they pass the resolved precision scalars in.
 *
 * @param covariance      Σ = K⁻¹.
 * @param i, j            Off-diagonal/diagonal indices (i < j).
 * @param precision_curr_ij, precision_curr_jj  Current K entries.
 * @param precision_prop_ij, precision_prop_jj  Proposed K entries.
 */
inline double log_det_ratio_edge_kernel(
        const arma::mat& covariance, size_t i, size_t j,
        double precision_curr_ij, double precision_prop_ij,
        double precision_curr_jj, double precision_prop_jj) {
    double Ui2 = precision_curr_ij - precision_prop_ij;
    double Uj2 = (precision_curr_jj - precision_prop_jj) / 2;

    double cc11 = covariance(j, j);
    double cc12 = 1.0 - (covariance(i, j) * Ui2 + covariance(j, j) * Uj2);
    double cc22 = Ui2 * Ui2 * covariance(i, i)
                + 2.0 * Ui2 * Uj2 * covariance(i, j)
                + Uj2 * Uj2 * covariance(j, j);

    return MY_LOG(std::abs(cc11 * cc22 - cc12 * cc12));
}

/**
 * Rank-1 specialisation of log_det_ratio_edge_kernel (Ui2 = 0): a diagonal-only
 * precision change at (j,j). Shared by GGMModel::log_det_ratio_diag and
 * MixedMRFModel::log_det_ratio_yy_diag.
 */
inline double log_det_ratio_diag_kernel(
        const arma::mat& covariance, size_t j,
        double precision_curr_jj, double precision_prop_jj) {
    double Uj2 = (precision_curr_jj - precision_prop_jj) / 2;

    double cc11 = covariance(j, j);
    double cc12 = 1.0 - covariance(j, j) * Uj2;
    double cc22 = Uj2 * Uj2 * covariance(j, j);

    return MY_LOG(std::abs(cc11 * cc22 - cc12 * cc12));
}

/**
 * Extract the six reparameterization constants for a rank-2 off-diagonal MH
 * proposal in precision space (Roverato move). Shared by GGMModel::get_constants
 * and MixedMRFModel::get_precision_constants, which differ only in where the
 * current precision entries come from (K vs −½K), passed in as scalars.
 *
 * @param logdet        log|K| (cached by the caller).
 * @param covariance    Σ = K⁻¹.
 * @param i, j          Off-diagonal/diagonal indices (i < j).
 * @param precision_ij, precision_jj  Current K entries K(i,j), K(j,j).
 * @return  {Phi_q1q, Phi_q1q1, c2, Phi_q1q1, c4, c5}.
 */
inline std::array<double, 6> precision_proposal_constants(
        double logdet, const arma::mat& covariance,
        size_t i, size_t j, double precision_ij, double precision_jj) {

    double log_adj_ii = logdet + MY_LOG(std::abs(covariance(i, i)));
    double log_adj_ij = logdet + MY_LOG(std::abs(covariance(i, j)));
    double log_adj_jj = logdet + MY_LOG(std::abs(covariance(j, j)));

    double inv_sub_jj = compute_inv_submatrix_i(covariance, i, j, j);
    double log_abs_inv_sub_jj = log_adj_ii + MY_LOG(std::abs(inv_sub_jj));

    double Phi_q1q  = (2 * std::signbit(covariance(i, j)) - 1) * MY_EXP(
        (log_adj_ij - (log_adj_jj + log_abs_inv_sub_jj) / 2)
    );
    double Phi_q1q1 = MY_EXP((log_adj_jj - log_abs_inv_sub_jj) / 2);

    std::array<double, 6> c{};
    c[0] = Phi_q1q;
    c[1] = Phi_q1q1;
    c[2] = precision_ij - Phi_q1q * Phi_q1q1;
    c[3] = Phi_q1q1;
    c[4] = precision_jj - Phi_q1q * Phi_q1q;
    c[5] = c[4] + c[2] * c[2] / (c[3] * c[3]);
    return c;
}

/** Overload computing log|K| from the Cholesky factor (callers without a cached log-determinant). */
inline std::array<double, 6> precision_proposal_constants(
        const arma::mat& cholesky, const arma::mat& covariance,
        size_t i, size_t j, double precision_ij, double precision_jj) {
    return precision_proposal_constants(get_log_det(cholesky), covariance,
                                        i, j, precision_ij, precision_jj);
}

} // namespace cholesky_helpers

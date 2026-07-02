// -----------------------------------------------------------------------------
// sbm_correction_test_interface.cpp
//
// R-facing test entries for the stochastic-block normalizing-constant
// corrections: the self-consistent per-edge slopes, the mini-TI label-move
// correction, and the corrected new-cluster collapsed marginal. With a
// constant slope curve all three have closed forms, so tests can pin the
// ported numerics without running a chain.
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>

#include "priors/edge_prior_correction.h"
#include "priors/sbm_edge_prior.h"

// -----------------------------------------------------------------------------
// test_sbm_compute_ce:
//   Self-consistent per-edge slope matrix for given labels (1-based) and
//   block probabilities.
// -----------------------------------------------------------------------------

// [[Rcpp::export(name = "test_sbm_compute_ce")]]
arma::mat test_sbm_compute_ce(
    const arma::ivec& cluster_assign,
    const arma::mat& block_probs,
    const arma::vec& fprime_density,
    const arma::vec& fprime,
    const arma::vec& quad_theta,
    const arma::vec& quad_f
) {
    SBMCorrection correction(fprime_density, fprime, quad_theta, quad_f);
    arma::uvec z = arma::conv_to<arma::uvec>::from(cluster_assign - 1);
    return compute_ce_sbm(z, block_probs, z.n_elem, correction);
}

// -----------------------------------------------------------------------------
// test_sbm_miniti_node:
//   Mini-TI normalizing-constant difference for moving one node (1-based
//   indices) between clusters.
// -----------------------------------------------------------------------------

// [[Rcpp::export(name = "test_sbm_miniti_node")]]
double test_sbm_miniti_node(
    int node,
    const arma::ivec& cluster_assign,
    const arma::mat& block_probs,
    int cur,
    int cand,
    const arma::vec& fprime_density,
    const arma::vec& fprime,
    const arma::vec& quad_theta,
    const arma::vec& quad_f
) {
    SBMCorrection correction(fprime_density, fprime, quad_theta, quad_f);
    arma::uvec z = arma::conv_to<arma::uvec>::from(cluster_assign - 1);
    arma::vec deg_base = degrees_ld_sbm(z, block_probs, z.n_elem, correction);
    return miniti_node_sbm(
        static_cast<arma::uword>(node - 1), z, block_probs, z.n_elem,
        static_cast<arma::uword>(cur - 1), static_cast<arma::uword>(cand - 1),
        deg_base, correction
    );
}

// -----------------------------------------------------------------------------
// test_sbm_miniti_removal:
//   Mini-TI edge-removal contribution for one node (1-based indices).
// -----------------------------------------------------------------------------

// [[Rcpp::export(name = "test_sbm_miniti_removal")]]
double test_sbm_miniti_removal(
    int node,
    const arma::ivec& cluster_assign,
    const arma::mat& block_probs,
    int cur,
    const arma::vec& fprime_density,
    const arma::vec& fprime,
    const arma::vec& quad_theta,
    const arma::vec& quad_f
) {
    SBMCorrection correction(fprime_density, fprime, quad_theta, quad_f);
    arma::uvec z = arma::conv_to<arma::uvec>::from(cluster_assign - 1);
    arma::vec deg_base = degrees_ld_sbm(z, block_probs, z.n_elem, correction);
    return miniti_removal_sbm(
        static_cast<arma::uword>(node - 1), z, block_probs, z.n_elem,
        static_cast<arma::uword>(cur - 1), deg_base, correction
    );
}

// -----------------------------------------------------------------------------
// test_sbm_corrected_log_marginal:
//   Corrected collapsed log-marginal for assigning a node (1-based) to a new
//   cluster.
// -----------------------------------------------------------------------------

// [[Rcpp::export(name = "test_sbm_corrected_log_marginal")]]
double test_sbm_corrected_log_marginal(
    int node,
    const arma::ivec& cluster_assign,
    const arma::imat& indicator,
    double alpha_between,
    double beta_between,
    const arma::vec& fprime_density,
    const arma::vec& fprime,
    const arma::vec& quad_theta,
    const arma::vec& quad_f
) {
    SBMCorrection correction(fprime_density, fprime, quad_theta, quad_f);
    arma::uvec z = arma::conv_to<arma::uvec>::from(cluster_assign - 1);
    arma::umat ind = arma::conv_to<arma::umat>::from(indicator);
    return corrected_log_marginal_mfm_sbm(
        z, ind, static_cast<arma::uword>(node - 1), z.n_elem,
        alpha_between, beta_between, correction
    );
}

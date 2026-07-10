#pragma once

/**
 * @file bgmCompare_state.h
 * @brief Persistent per-group sweep state for the bgmCompare sampler.
 */

#include <RcppArmadillo.h>
#include <vector>

/**
 * Persistent per-group state shared by the bgmCompare Gibbs sweeps.
 *
 * Holds the double-precision observation matrices, the group-specific
 * effective pairwise weights, and the residual (rest-score) matrices
 * `residual[g] = obs_double[g] * pairwise_group[g]`. The state is owned by
 * run_gibbs_sampler_bgmCompare() and threaded by reference through the
 * sweeps, which keep it consistent with the model parameters:
 *  - accepted pairwise-effect moves and accepted difference-indicator flips
 *    update `pairwise_group` and two residual columns per group;
 *  - accepted main-effect moves leave the state untouched;
 *  - imputation refills the observation matrices and triggers a full weight
 *    rebuild once per pass;
 *  - NUTS updates and graph re-initialisation trigger a full weight rebuild
 *    via rebuild_sweep_state_weights().
 *
 * The state also caches the per-variable, per-group log-normalizer sums of
 * the pseudolikelihood, `log_normalizer(v, g) = sum_i [bound_i + log
 * denom_i]` at the current parameters and rest scores. A variable's row is
 * valid while its main effects and its residual columns are unchanged:
 * accepted Metropolis moves store the proposed-state normalizers, accepted
 * indicator flips invalidate the affected variables, and weight rebuilds
 * invalidate every variable.
 */
struct CompareSweepState {
  arma::mat obs_double_all;              ///< Observations as double (n x V)
  std::vector<arma::mat> obs_double;     ///< Per-group observations as double (n_g x V)
  std::vector<arma::mat> pairwise_group; ///< Per-group effective pairwise weights (V x V)
  std::vector<arma::mat> residual;       ///< Per-group rest scores obs_double[g] * pairwise_group[g] (n_g x V)
  arma::mat log_normalizer;              ///< Cached pseudolikelihood log-normalizer sums (V x G)
  arma::uvec normalizer_valid;           ///< 1 if a variable's log_normalizer row is current (V)
};

/**
 * Fill the observation members of a CompareSweepState.
 *
 * Converts the integer observation matrix to double, both as a whole and
 * per group. Does not touch the weight or residual members.
 *
 * @param[out] state          Sweep state to fill
 * @param observations        Integer observation matrix (n x V)
 * @param group_indices       Group start/end row indices per group (G x 2)
 * @param num_groups          Number of groups (G)
 */
void initialize_sweep_state_observations(
    CompareSweepState& state,
    const arma::imat& observations,
    const arma::imat& group_indices,
    const int num_groups
);

/**
 * Rebuild the effective pairwise weights and residual matrices.
 *
 * Recomputes `pairwise_group[g]` from the pairwise-effect matrix and the
 * inclusion indicators, and `residual[g]` as one matrix product per group.
 * Invalidates the normalizer cache for every variable. Called after
 * wholesale parameter changes (NUTS updates, graph initialisation) and once
 * at sampler start.
 *
 * @param[in,out] state           Sweep state with observation members filled
 * @param pairwise_effects        Pairwise-effect matrix (rows = pairs, cols = groups)
 * @param pairwise_effect_indices Row index per variable pair (V x V)
 * @param inclusion_indicator     Difference inclusion indicators (V x V)
 * @param projection              Group contrast matrix (G x (G-1))
 * @param num_groups              Number of groups (G)
 */
void rebuild_sweep_state_weights(
    CompareSweepState& state,
    const arma::mat& pairwise_effects,
    const arma::imat& pairwise_effect_indices,
    const arma::imat& inclusion_indicator,
    const arma::mat& projection,
    const int num_groups
);

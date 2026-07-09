#pragma once

/**
 * @file bgmCompare_logp_and_grad.h
 * @brief Log-pseudoposterior and gradient functions for bgmCompare.
 *
 * Free functions implementing the log-pseudoposterior, its gradient, and
 * component-wise evaluations used by Metropolis and NUTS samplers in
 * the multi-group comparison model.
 *
 * Parameter columns are indexed by h: h = 0 is the overall (shared/baseline)
 * effect; h > 0 are group-difference effects projected through the contrast
 * matrix.
 */

#include <RcppArmadillo.h>
#include "priors/parameter_prior.h"


/**
 * Compute the total length of the flat parameter vector.
 *
 * Accounts for which group-difference entries are active based on the
 * current inclusion indicators.
 *
 * @param num_variables            Number of variables (V)
 * @param main_effect_indices      Start/end row indices per variable (V x 2)
 * @param pairwise_effect_indices  Row index per variable pair (V x V)
 * @param inclusion_indicator      Edge inclusion indicators (V x V)
 * @param num_categories           Number of categories per variable
 * @param is_ordinal_variable      1 = ordinal, 0 = Blume-Capel
 * @param num_groups               Number of groups (G)
 * @return Total number of active parameters
 */
arma::uword total_length(
    const int num_variables,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::imat& inclusion_indicator,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable,
    const int num_groups
);

/**
 * Compute the observed-data contribution to the gradient.
 *
 * Projects sufficient statistics into the active-parameter layout using
 * the group contrast matrix and the index maps from build_index_maps().
 *
 * @param main_effect_indices          Start/end row indices per variable (V x 2)
 * @param pairwise_effect_indices      Row index per variable pair (V x V)
 * @param projection                   Group contrast matrix (G x (G-1))
 * @param observations                 Integer observation matrix (n x V)
 * @param group_indices                Group start/end indices per group (G x 2)
 * @param num_categories               Number of categories per variable
 * @param inclusion_indicator          Edge inclusion indicators (V x V)
 * @param counts_per_category_group    Category counts per group
 * @param blume_capel_stats_group      Blume-Capel sufficient statistics per group
 * @param pairwise_stats_group         Pairwise sufficient statistics per group
 * @param num_groups                   Number of groups (G)
 * @param is_ordinal_variable          1 = ordinal, 0 = Blume-Capel
 * @param baseline_category            Reference categories for Blume-Capel variables
 * @param main_index                   Main-effect index map from build_index_maps()
 * @param pair_index                   Pairwise index map from build_index_maps()
 * @return Gradient vector (observed-data component only)
 */
arma::vec gradient_observed_active(
  const arma::imat& main_effect_indices,
  const arma::imat& pairwise_effect_indices,
  const arma::mat& projection,
  const arma::imat& observations,
  const arma::imat& group_indices,
  const arma::ivec& num_categories,
  const arma::imat& inclusion_indicator,
  const std::vector<arma::imat>& counts_per_category_group,
  const std::vector<arma::imat>& blume_capel_stats_group,
  const std::vector<arma::mat>&  pairwise_stats_group,
  const int num_groups,
  const arma::uvec& is_ordinal_variable,
  const arma::ivec& baseline_category,
  const arma::imat main_index,
  const arma::imat pair_index
);

/**
 * Compute the full gradient of the log-pseudoposterior.
 *
 * Combines observed sufficient statistics (grad_obs), expected sufficient
 * statistics (computed on-the-fly via softmax probabilities), and prior
 * gradient terms (logistic-Beta for baselines, Cauchy for differences).
 *
 * @param main_effects             Current main-effect matrix
 * @param pairwise_effects         Current pairwise-effect matrix
 * @param main_effect_indices      Start/end row indices per variable (V x 2)
 * @param pairwise_effect_indices  Row index per variable pair (V x V)
 * @param projection               Group contrast matrix (G x (G-1))
 * @param observations_double      Observations as double (n x V)
 * @param group_indices            Group start/end indices per group (G x 2)
 * @param num_categories           Number of categories per variable
 * @param counts_per_category_group Category counts per group
 * @param blume_capel_stats_group  Blume-Capel sufficient statistics per group
 * @param pairwise_stats_group     Pairwise sufficient statistics per group
 * @param num_groups               Number of groups (G)
 * @param inclusion_indicator      Edge inclusion indicators (V x V)
 * @param is_ordinal_variable      1 = ordinal, 0 = Blume-Capel
 * @param baseline_category        Reference categories for Blume-Capel variables
 * @param main_index               Main-effect index map from build_index_maps()
 * @param pair_index               Pairwise index map from build_index_maps()
 * @param grad_obs                 Pre-computed observed-data gradient
 * @param interaction_prior        Prior on baseline pairwise effects
 * @param difference_prior         Prior on group-difference parameters
 * @param threshold_prior          Prior on baseline main effects
 * @return Full gradient vector
 */
arma::vec gradient(
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::mat& projection,
    const arma::mat& observations_double,
    const arma::imat& group_indices,
    const arma::ivec& num_categories,
    const std::vector<arma::imat>& counts_per_category_group,
    const std::vector<arma::imat>& blume_capel_stats_group,
    const std::vector<arma::mat>&  pairwise_stats_group,
    const int num_groups,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const arma::imat& main_index,
    const arma::imat& pair_index,
    const arma::vec& grad_obs,
    const BaseParameterPrior& interaction_prior,
    const BaseParameterPrior& difference_prior,
    const BaseParameterPrior& threshold_prior
);

/**
 * Compute the log-pseudoposterior and its gradient in a single pass.
 *
 * Shares intermediate computations (group-specific effects, residual
 * matrices, probability vectors) to avoid redundant work during NUTS.
 *
 * @return Pair of (log-pseudoposterior value, gradient vector)
 * @see gradient() for parameter descriptions
 */
std::pair<double, arma::vec> logp_and_gradient(
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::mat& projection,
    const arma::mat& observations_double,
    const arma::imat& group_indices,
    const arma::ivec& num_categories,
    const std::vector<arma::imat>& counts_per_category_group,
    const std::vector<arma::imat>& blume_capel_stats_group,
    const std::vector<arma::mat>&  pairwise_stats_group,
    const int num_groups,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const arma::imat& main_index,
    const arma::imat& pair_index,
    const arma::vec& grad_obs,
    const BaseParameterPrior& interaction_prior,
    const BaseParameterPrior& difference_prior,
    const BaseParameterPrior& threshold_prior
);

/**
 * Log-pseudoposterior contribution of a single main-effect parameter.
 *
 * Used by element-wise Metropolis updates. Evaluates the pseudolikelihood
 * and prior for one (variable, category/par, column h) entry. Rest scores
 * are read from the maintained residual matrices; the pairwise effects do
 * not change during main-effect updates.
 *
 * @param residual_groups  Per-group rest-score matrices (n_g x V)
 * @param variable   Variable index
 * @param category   Category index (ordinal variables only)
 * @param par        Parameter index: 0 = linear, 1 = quadratic (Blume-Capel only)
 * @param h          Column index: 0 = overall baseline, >0 = group difference
 * @param normalizers_in   Optional cached per-group log-normalizer sums for
 *                         the variable (length G); skips their computation
 * @param normalizers_out  Optional output for the computed per-group
 *                         log-normalizer sums (length G)
 * @see gradient() for remaining parameter descriptions
 */
double log_pseudoposterior_main_component(
    const arma::mat& main_effects,
    const arma::imat& main_effect_indices,
    const arma::mat& projection,
    const std::vector<arma::mat>& residual_groups,
    const arma::ivec& num_categories,
    const std::vector<arma::imat>& counts_per_category_group,
    const std::vector<arma::imat>& blume_capel_stats_group,
    const int num_groups,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    int variable,
    int category,
    int par,
    int h,
    const BaseParameterPrior& difference_prior,
    const BaseParameterPrior& threshold_prior,
    const arma::vec* normalizers_in = nullptr,
    arma::vec* normalizers_out = nullptr
);

/**
 * Log-pseudoposterior contribution of a single pairwise-effect parameter.
 *
 * Uses pre-computed residual matrices adjusted by delta to avoid full
 * recomputation. Used by element-wise Metropolis updates.
 *
 * @param obs_double_groups  Per-group observation matrices converted to double
 * @param residual_matrices  Pre-computed residual matrices per group
 * @param variable1          First variable index
 * @param variable2          Second variable index
 * @param h                  Column index: 0 = overall baseline, >0 = group difference
 * @param delta              Proposed change to pairwise_effects(idx, h)
 * @param normalizers_in     Optional cached per-group log-normalizer sums for
 *                           the endpoint variables (G x 2); skips their computation
 * @param normalizers_out    Optional output for the computed per-group
 *                           log-normalizer sums (G x 2)
 * @see gradient() for remaining parameter descriptions
 */
double log_pseudoposterior_pair_component(
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::mat& projection,
    const std::vector<arma::mat>& obs_double_groups,
    const arma::ivec& num_categories,
    const std::vector<arma::mat>& pairwise_stats_group,
    const std::vector<arma::mat>& residual_matrices,
    const int num_groups,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    int variable1,
    int variable2,
    int h,
    double delta,
    const BaseParameterPrior& interaction_prior,
    const BaseParameterPrior& difference_prior,
    const arma::mat* normalizers_in = nullptr,
    arma::mat* normalizers_out = nullptr
);


/**
 * Log-pseudolikelihood ratio for toggling a variable's main-effect differences.
 *
 * Compares proposed vs. current main-effect parameters across all groups,
 * combining sufficient-statistic differences with normalizing-constant ratios.
 * Both states share the rest scores held in the maintained residual matrices.
 * Used by the Metropolis-Hastings indicator update for main effects.
 *
 * @param current_main_effects   Current main-effect matrix
 * @param proposed_main_effects  Proposed main-effect matrix
 * @param residual_groups        Per-group rest-score matrices (n_g x V)
 * @param variable               Variable whose main effect is being toggled
 * @see gradient() for remaining parameter descriptions
 */
double log_pseudolikelihood_ratio_main(
    const arma::mat& current_main_effects,
    const arma::mat& proposed_main_effects,
    const arma::imat& main_effect_indices,
    const arma::mat&  projection,
    const std::vector<arma::mat>& residual_groups,
    const arma::ivec& num_categories,
    const std::vector<arma::imat>& counts_per_category_group,
    const std::vector<arma::imat>& blume_capel_stats_group,
    const int num_groups,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const int variable
);

/**
 * Log-pseudolikelihood ratio for toggling a pairwise interaction's differences.
 *
 * Compares proposed vs. current pairwise-effect parameters for a single edge,
 * summing the data contribution and normalizing-constant ratios for both
 * endpoint variables. Current-state rest scores come from the maintained
 * residual matrices; proposed-state rest scores adjust them by the change in
 * the group-specific effective weight. Used by the Metropolis-Hastings
 * indicator update.
 *
 * @param current_pairwise_effects   Current pairwise-effect matrix
 * @param proposed_pairwise_effects  Proposed pairwise-effect matrix
 * @param obs_double_groups          Per-group observation matrices as double
 * @param residual_groups            Per-group rest-score matrices (n_g x V)
 * @param var1                       First variable index
 * @param var2                       Second variable index
 * @see gradient() for remaining parameter descriptions
 */
double log_pseudolikelihood_ratio_pairwise(
    const arma::mat& main_effects,
    const arma::mat& current_pairwise_effects,
    const arma::mat& proposed_pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::mat& projection,
    const std::vector<arma::mat>& obs_double_groups,
    const std::vector<arma::mat>& residual_groups,
    const arma::ivec& num_categories,
    const std::vector<arma::mat>& pairwise_stats_group,
    const int num_groups,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const int var1,
    const int var2
);
#pragma once

/**
 * @file bgmCompare_sampler.h
 * @brief MCMC chain runner for the bgmCompare multi-group comparison model.
 */

#include <RcppArmadillo.h>
#include "utils/common_helpers.h"
#include "models/bgmCompare/bgmCompare_output.h"
#include "priors/parameter_prior.h"
#include "priors/edge_prior.h"
#include <string>

struct SafeRNG;
class ProgressManager;

/**
 * Run a complete MCMC chain for a bgmCompare model.
 *
 * Initializes parameters and adaptation controllers, loops over warmup and
 * sampling iterations (with optional missing-data imputation, difference
 * selection, and prior updates), and collects posterior samples and
 * diagnostics into a bgmCompareOutput struct.
 *
 * @param chain_id                   Identifier for this chain
 * @param observations               Integer observation matrix (n x V)
 * @param num_groups                 Number of groups (G)
 * @param counts_per_category        Category counts per group (modified during imputation)
 * @param blume_capel_stats          Blume-Capel statistics per group (modified during imputation)
 * @param pairwise_stats             Pairwise sufficient statistics per group (modified during imputation)
 * @param num_categories             Number of categories per variable
 * @param difference_selection_alpha Alpha for difference-selection Beta-Bernoulli prior
 * @param difference_selection_beta  Beta for difference-selection Beta-Bernoulli prior
 * @param difference_prior_type      Prior family for inclusion probabilities (e.g. "Beta-Bernoulli")
 * @param iter                       Number of post-warmup iterations
 * @param warmup                     Number of warmup iterations
 * @param na_impute                  Enable missing-data imputation
 * @param missing_data_indices       (person, variable) pairs of missing entries
 * @param is_ordinal_variable        1 = ordinal, 0 = Blume-Capel
 * @param baseline_category          Reference categories for Blume-Capel variables
 * @param difference_selection       Enable group-difference indicator updates
 * @param main_difference_selection  Enable main-effect difference selection
 * @param main_effect_indices        Start/end row indices per variable (V x 2)
 * @param pairwise_effect_indices    Row index per variable pair (V x V)
 * @param target_accept              Target acceptance rate for NUTS dual averaging
 * @param nuts_max_depth             Maximum NUTS tree depth
 * @param learn_mass_matrix          Enable mass-matrix adaptation during warmup
 * @param projection                 Group contrast matrix (G x (G-1))
 * @param group_membership           Group label per observation
 * @param group_indices              Group start/end indices per group (G x 2)
 * @param interaction_index_matrix   Edge-pair index matrix
 * @param inclusion_probability      Prior inclusion probabilities (V x V)
 * @param rng                        Random number generator
 * @param update_method              Sampler type (adaptive_metropolis, nuts)
 * @param pm                         Progress manager for user feedback
 * @param interaction_prior          Prior for the baseline pairwise interaction parameters
 * @param difference_prior           Prior for the group-difference parameters;
 *                                   its scale is mutated in place when a random
 *                                   difference scale is enabled, so it must be a
 *                                   per-chain object
 * @param threshold_prior            Prior for the threshold (main effect) parameters
 * @param difference_scale_prior     Mean-1 hyperprior on the difference-scale
 *                                   multiplier u = s / s0; nullptr keeps the
 *                                   difference scale fixed
 * @param difference_edge_prior      Indicator prior for difference selection;
 *                                   for Stochastic-Block, governs off-diagonal
 *                                   (pairwise) inclusions and exposes block
 *                                   allocations via has_allocations()
 * @return bgmCompareOutput containing posterior samples and diagnostics
 */

bgmCompareOutput run_gibbs_sampler_bgmCompare(
    int chain_id,
    arma::imat observations,
    const int num_groups,
    std::vector<arma::imat>& counts_per_category,
    std::vector<arma::imat>& blume_capel_stats,
    std::vector<arma::mat>& pairwise_stats,
    const arma::ivec& num_categories,
    const double difference_selection_alpha,
    const double difference_selection_beta,
    const std::string& difference_prior_type,
    const int iter,
    const int warmup,
    const bool na_impute,
    const arma::imat& missing_data_indices,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const bool difference_selection,
    const bool main_difference_selection,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const double target_accept,
    const int nuts_max_depth,
    const bool learn_mass_matrix,
    const arma::mat& projection,
    const arma::ivec& group_membership,
    const arma::imat& group_indices,
    const arma::imat& interaction_index_matrix,
    arma::mat inclusion_probability,
    SafeRNG& rng,
    const UpdateMethod update_method,
    ProgressManager& pm,
    const BaseParameterPrior& interaction_prior,
    BaseParameterPrior& difference_prior,
    const BaseParameterPrior& threshold_prior,
    BaseEdgePrior& difference_edge_prior,
    const BaseParameterPrior* difference_scale_prior
);
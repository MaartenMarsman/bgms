#include <RcppArmadillo.h>

/**
 * Function: log_pseudolikelihood_ratio_interaction
 *
 * Computes the change in log pseudo-likelihood when a pairwise interaction parameter
 * between two variables is updated from a current value to a proposed value.
 *
 * This function evaluates:
 *  1. The direct contribution of the interaction to the joint linear predictor.
 *  2. The change in pseudo-likelihood for each affected variable (via softmax terms),
 *     accounting for the influence of the interaction on their respective rest scores.
 *
 * Inputs:
 *  - pairwise_effects: Current matrix of pairwise interaction parameters [V × V].
 *  - main_effects: Matrix of main effect (threshold) parameters [V × max_categories].
 *  - observations: Integer matrix of observed category scores [N × V].
 *  - num_categories: Vector of category counts per variable [V].
 *  - num_persons: Number of individuals (rows in the data).
 *  - variable1, variable2: Indices of the two interacting variables (0-based).
 *  - proposed_state: Proposed new interaction weight.
 *  - current_state: Current interaction weight.
 *  - residual_matrix: Matrix of residual linear predictors [N × V].
 *  - is_ordinal_variable: Logical vector indicating whether each variable is ordinal (1) or BC (0).
 *  - reference_category: Vector of reference categories per variable (used for BC variables).
 *
 * Returns:
 *  - The log pseudo-likelihood ratio:
 *      log p(y | β_proposed) - log p(y | β_current)
 */
double log_pseudolikelihood_ratio_interaction (
    const arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const int num_persons,
    const int variable1,
    const int variable2,
    const double proposed_state,
    const double current_state,
    const arma::mat& residual_matrix,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category
);


/**
 * Function: gradient_log_pseudoposterior_interaction_single
 *
 * Computes the gradient of the log pseudoposterior with respect to a single
 * interaction parameter (i, j), assuming symmetric pairwise_effects matrix.
 *
 * Inputs:
 *  - i, j: Variable indices for the interaction pair (i < j).
 *  - pairwise_effects: Matrix of interaction parameters.
 *  - main_effects: Matrix of main effect parameters.
 *  - observations: Matrix of categorical data [N × V].
 *  - num_categories: Vector with number of categories per variable.
 *  - is_ordinal_variable: Logical vector (1 = ordinal, 0 = Blume-Capel).
 *  - reference_category: Vector of reference categories for BC variables.
 *  - interaction_scale: Scale of Cauchy prior.
 *
 * Returns:
 *  - A single double value: the gradient with respect to β_ij.
 */
double gradient_log_pseudoposterior_interaction_single (
    int var1,
    int var2,
    const arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    const double interaction_scale
);
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

#include "gibbs_functions.h"
#include <Rcpp.h>

using namespace Rcpp;


// ----------------------------------------------------------------------------|
// List of to do's
// 1. Add between model moves for main effects.
// 2. Add g prior.
// 3. Add "enum" to handle threshold scenarios
// ----------------------------------------------------------------------------|

/**
 * Function: update_with_robbins_monro
 * Purpose: Performs Robbins-Monro updates for proposal standard deviations.
 * Inputs:
 *  - current_sd: Current standard deviation of the proposal.
 *  - observed_log_acceptance_probability: Log acceptance probability from the Metropolis-Hastings step.
 *  - target_acceptance_rate: Target acceptance rate.
 *  - t: Iteration number.
 *  - rm_adaptation_rate: Robbins-Monro learning rate.
 *  - rm_lower_bound: Minimum allowable standard deviation.
 *  - rm_upper_bound: Maximum allowable standard deviation.
 * Returns:
 *  - Updated proposal standard deviation, clamped within bounds.
 */
double update_with_robbins_monro (
    const double current_sd,
    const double observed_log_acceptance_probability,
    const double target_acceptance_rate,
    const double rm_lower_bound,
    const double rm_upper_bound,
    const double rm_weight
) {

  // Normalize the acceptance probability
  double observed_acceptance_probability = 1.0;
  if(observed_log_acceptance_probability < 0.0) {
    observed_acceptance_probability = std::exp(observed_log_acceptance_probability);
  }

  // Update the proposal standard deviation
  double update = current_sd +
    (observed_acceptance_probability - target_acceptance_rate) * rm_weight;

  // Handle NaN cases by resetting to default value
  if (std::isnan(update)) {
    update = 1.0; // Default proposal standard deviation
  }

  // Clamp the updated standard deviation within bounds
  return std::clamp(update, rm_lower_bound, rm_upper_bound);
}


/**
 * Function: compute_group_thresholds
 * Purpose: Computes thresholds for a specific variable in a given group.
 * Inputs:
 *  - variable: Index of the variable for which thresholds are computed.
 *  - is_ordinal_variable: Logical vector indicating if variables are ordinal.
 *  - group: Index of the group for which thresholds are computed.
 *  - num_groups: Total number of groups in the analysis.
 *  - num_categories: Matrix of category counts per variable and group.
 *  - main_effects: Matrix of main effects across variables and groups.
 *  - main_effect_indices: Indices for main effect parameters.
 *  - projection: Projection matrix for group-specific scaling.
 *  - independent_thresholds: Whether thresholds are modeled independently.
 * Outputs:
 *  - A `arma::vec` containing the computed thresholds for the variable in the group.
 */
arma::vec compute_group_thresholds(
    const int variable,
    const arma::uvec& is_ordinal_variable,
    const int group,
    const int num_groups,
    const arma::imat& num_categories,
    const arma::mat& main_effects,
    const arma::imat& main_effect_indices,
    const arma::mat& projection,
    const bool independent_thresholds
) {
  // Determine the length of the threshold vector based on whether the variable is ordinal
  int vector_length = (is_ordinal_variable[variable]) ? num_categories(variable, group) : 2;
  arma::vec GroupThresholds(vector_length);  // Output vector for thresholds

  // Base index for accessing main effects for this variable
  int base_category_index = main_effect_indices(variable, 0);

  if (!independent_thresholds) {
    // Dependent thresholds: model group differences
    int n_cats = num_categories(variable, group);
    if (is_ordinal_variable[variable]) {
      // Regular binary or ordinal variable
      arma::vec category_effects = main_effects.rows(base_category_index, base_category_index + n_cats - 1).col(0);
      arma::mat projection_contributions = projection.row(group) * main_effects.rows(base_category_index, base_category_index + n_cats - 1).cols(1, num_groups - 1).t();
      GroupThresholds = category_effects + projection_contributions.t();
    } else {
      // Blume-Capel ordinal variable

      // Extract the rows corresponding to the two threshold parameters
      arma::mat threshold_effects = main_effects.submat(base_category_index, 0, base_category_index + 1, 0);

      // Add group-specific contributions for both rows at once
      threshold_effects += projection.row(group) * main_effects.submat(base_category_index, 1, base_category_index + 1, num_groups - 1).t();

      // Assign the computed group-specific thresholds
      GroupThresholds[0] = threshold_effects(0, 0);
      GroupThresholds[1] = threshold_effects(1, 0);
    }
  } else {
    // Independent thresholds: compute separately for each group
    if (is_ordinal_variable[variable]) {
      // Regular binary or ordinal variable
      GroupThresholds = main_effects.submat(base_category_index, group, base_category_index + vector_length - 1, group);
    } else {
      // Blume-Capel ordinal variable
      GroupThresholds[0] = main_effects(base_category_index, group);
      GroupThresholds[1] = main_effects(base_category_index + 1, group);
    }
  }

  return GroupThresholds;
}


/**
 * Function: impute_missing_data_for_anova_model
 * Purpose: Imputes missing data for independent samples designs by generating new observations
 *          based on the model parameters and pseudo-likelihood.
 *
 * Inputs:
 *  - main_effects: Numeric matrix of main effects across variables and groups.
 *  - pairwise_effects: Numeric matrix of pairwise interaction effects between variables across groups.
 *  - main_effect_indices: Integer matrix mapping variable indices to main effect parameters.
 *  - pairwise_effect_indices: Integer matrix mapping variable pairs to pairwise interaction parameters.
 *  - projection: Numeric matrix representing group-specific scaling for effects.
 *  - observations: Integer matrix of observed data (individuals x variables), with missing data encoded.
 *  - num_groups: Number of groups in the analysis.
 *  - group_membership: Integer vector mapping individuals to their respective groups.
 *  - n_cat_obs: List of matrices, one per group, recording category frequencies per variable.
 *  - sufficient_blume_capel: List of matrices storing sufficient statistics for Blume-Capel variables.
 *  - num_categories: Integer matrix of category counts for each variable and group.
 *  - residual_matrix: Numeric matrix of residual effects for pseudo-likelihood calculations.
 *  - missing_data_indices: Integer matrix of indices indicating missing observations (row x column pairs).
 *  - is_ordinal_variable: Logical vector indicating whether variables are ordinal.
 *  - baseline_category: Integer vector of reference categories for Blume-Capel variables.
 *  - independent_thresholds: Boolean indicating whether thresholds are modeled independently.
 *
 * Outputs:
 *  - A List containing:
 *    - `observations`: Updated observation matrix with imputed values.
 *    - `n_cat_obs`: Updated list of category counts per group.
 *    - `sufficient_blume_capel`: Updated sufficient statistics for Blume-Capel variables.
 *    - `residual_matrix`: Updated residual effects matrix.
 */
List impute_missing_data_for_anova_model(
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::mat& projection,
    arma::imat& observations,
    const int num_groups,
    const arma::ivec& group_membership,
    List& n_cat_obs,
    List& sufficient_blume_capel,
    const arma::imat& num_categories,
    arma::mat& residual_matrix,
    const arma::imat& missing_data_indices,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    bool independent_thresholds
) {
  const int num_variables = observations.n_cols;
  const int num_missings = missing_data_indices.n_rows;
  const int max_num_categories = arma::max(arma::vectorise(num_categories));

  arma::vec category_response_probabilities(max_num_categories + 1);
  double exponent, rest_score, cumsum, u, GroupInteraction;
  int score, person, variable, new_observation, old_observation, gr, int_index;

  //Impute missing data
  for(int missing = 0; missing < num_missings; missing++) {
    // Identify the observation to impute
    person = missing_data_indices(missing, 0);
    variable = missing_data_indices(missing, 1);
    gr = group_membership[person];

    // Compute thresholds for the variable in the given group
    arma::vec GroupThresholds = compute_group_thresholds(variable,
                                                         is_ordinal_variable,
                                                         gr,
                                                         num_groups,
                                                         num_categories,
                                                         main_effects,
                                                         main_effect_indices,
                                                         projection,
                                                         independent_thresholds);

    // Generate a new observation based on the model
    rest_score = residual_matrix(person, variable);
    if(is_ordinal_variable[variable] == true) {
      // For regular binary or ordinal variables
      cumsum = 1.0;
      category_response_probabilities[0] = 1.0;
      for(int category = 1; category <= num_categories(variable, gr); category++) {
        exponent = GroupThresholds(category - 1);
        exponent += category * rest_score;
        cumsum += std::exp(exponent);
        category_response_probabilities[category] = cumsum;
      }
    } else {
      // For Blume-Capel variables
      cumsum = 0.0;
      for(int category = 0; category <= num_categories(variable, gr); category++) {
        exponent = GroupThresholds[0] * category;
        exponent += GroupThresholds[1] *
          (category - baseline_category[variable]) *
          (category - baseline_category[variable]);
        exponent += category * rest_score;
        cumsum += std::exp(exponent);
        category_response_probabilities[category] = cumsum;
      }
    }

    // Sample a new value based on computed probabilities
    u = cumsum * R::unif_rand();
    score = 0;
    while (u > category_response_probabilities[score]) {
      score++;
    }
    new_observation = score;
    old_observation = observations(person, variable);

    if(old_observation != new_observation) {
      // Update raw observations
      observations(person, variable) = new_observation;

      // Update category counts or sufficient statistics
      if(is_ordinal_variable[variable] == true) {
        arma::imat n_cat_obs_gr = n_cat_obs[gr];
        n_cat_obs_gr(old_observation, variable)--;
        n_cat_obs_gr(new_observation, variable)++;
        n_cat_obs[gr] = n_cat_obs_gr;
      } else {
        arma::imat sufficient_blume_capel_gr = sufficient_blume_capel[gr];
        sufficient_blume_capel_gr(0, variable) -= old_observation;
        sufficient_blume_capel_gr(0, variable) += new_observation;
        sufficient_blume_capel_gr(1, variable) -=
          (old_observation - baseline_category[variable]) *
          (old_observation - baseline_category[variable]);
        sufficient_blume_capel_gr(1, variable) +=
          (new_observation - baseline_category[variable]) *
          (new_observation - baseline_category[variable]);
        sufficient_blume_capel[gr] = sufficient_blume_capel_gr;
      }

      // Update rest scores
      for(int vertex = 0; vertex < num_variables; vertex++) {
        int_index = pairwise_effect_indices(vertex, variable);
        GroupInteraction = pairwise_effects(int_index, 0);
        for(int h = 0; h < num_groups - 1; h++) {
          GroupInteraction += projection(gr, h) * pairwise_effects(int_index, h + 1);
        }
        residual_matrix(person, vertex) -= old_observation * GroupInteraction;
        residual_matrix(person, vertex) += new_observation * GroupInteraction;
      }
    }
  }

  return List::create(Named("observations") = observations,
                      Named("n_cat_obs") = n_cat_obs,
                      Named("sufficient_blume_capel") = sufficient_blume_capel,
                      Named("residual_matrix") = residual_matrix);
}


/**
 * Function: log_pseudolikelihood_ratio_interaction
 * Purpose: Computes the log pseudo-likelihood ratio for proposed vs. current interaction states
 *          for a multi-group independent samples design.
 *
 * Inputs:
 *  - main_effects: Numeric matrix of main effects for variables across groups.
 *  - main_effect_indices: Integer matrix mapping variables to main effect parameters.
 *  - projection: Numeric matrix of group-specific scaling factors for interactions.
 *  - observations: Integer matrix of observed data (individuals x variables).
 *  - num_groups: Number of groups in the analysis.
 *  - group_indices: Integer matrix containing start and end indices for each group in `observations`.
 *  - num_categories: Integer matrix indicating the number of categories for each variable and group.
 *  - independent_thresholds: Boolean indicating whether thresholds are modeled independently.
 *  - num_persons: Total number of individuals in the dataset.
 *  - variable1: Index of the first variable involved in the interaction.
 *  - variable2: Index of the second variable involved in the interaction.
 *  - proposed_state: Proposed value for the interaction parameter.
 *  - current_state: Current value of the interaction parameter.
 *  - residual_matrix: Numeric matrix of residual effects for pseudo-likelihood calculations.
 *  - is_ordinal_variable: Logical vector indicating whether each variable is ordinal.
 *  - baseline_category: Integer vector indicating the reference categories for Blume-Capel variables.
 *
 * Outputs:
 *  - Returns the log pseudo-likelihood ratio for the proposed vs. current interaction states.
 */
double log_pseudolikelihood_ratio_interaction(
    const arma::mat& main_effects,
    const arma::imat& main_effect_indices,
    const arma::mat& projection,
    const arma::imat& observations,
    const int num_groups,
    const arma::imat& group_indices,
    const arma::imat& num_categories,
    const bool independent_thresholds,
    const int num_persons,
    const int variable1,
    const int variable2,
    const double proposed_state,
    const double current_state,
    const arma::mat& residual_matrix,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category
) {
  // Compute the difference in interaction states
  double delta_state = 2.0 * (proposed_state - current_state);
  double pseudolikelihood_ratio = 0.0;
  int n_cats_v1, n_cats_v2, obs_score1, obs_score2, var, obs, n_cats, score;
  double obs_current, obs_proposed, rest_score, bound, denominator_prop;
  double denominator_curr, exponent;

  // Loop over groups to compute contributions to the pseudo-likelihood ratio
  for(int gr = 0; gr < num_groups; gr++) {
    // Retrieve thresholds for the two variables in the current group
    arma::vec GroupThresholds_v1 = compute_group_thresholds (
      variable1, is_ordinal_variable, gr, num_groups, num_categories, main_effects,
      main_effect_indices, projection, independent_thresholds);

    arma::vec GroupThresholds_v2 = compute_group_thresholds (
      variable2, is_ordinal_variable, gr, num_groups, num_categories, main_effects,
      main_effect_indices, projection, independent_thresholds);

    n_cats_v1 = num_categories(variable1, gr);
    n_cats_v2 = num_categories(variable2, gr);

    // Loop over individuals in the current group
    for(int person = group_indices(gr, 0); person <= group_indices(gr, 1); person++) {
      obs_score1 = observations(person, variable1);
      obs_score2 = observations(person, variable2);

      // Contribution from the interaction term
      pseudolikelihood_ratio += obs_score1 * obs_score2 * delta_state;

      // Compute contributions for each variable (variable1 and variable2)
      for (int variable = 1; variable <= 2; variable++) {
        // Assign variable-specific data
        var = (variable == 1) ? variable1 : variable2;
        obs = (variable == 1) ? obs_score2 : obs_score1;
        obs_current = obs * current_state;
        obs_proposed = obs * proposed_state;

        n_cats = (variable == 1) ? n_cats_v1 : n_cats_v2;
        arma::vec& GroupThresholds = (variable == 1) ? GroupThresholds_v1 : GroupThresholds_v2;

        // Compute the rest score and bounds
        rest_score = residual_matrix(person, var) - obs_current;
        bound = (rest_score > 0) ? n_cats * rest_score : 0.0;
        denominator_prop = 0.0;
        denominator_curr = 0.0;

        // Compute pseudo-likelihood terms
        if(is_ordinal_variable[var]) {
          // Regular binary or ordinal MRF variable
          denominator_prop += std::exp(-bound);
          denominator_curr += std::exp(-bound);
          for(int cat = 0; cat < n_cats; cat++) {
            score = cat + 1;
            exponent = GroupThresholds[cat] + score * rest_score - bound;
            denominator_prop += std::exp(exponent + score * obs_proposed);
            denominator_curr += std::exp(exponent + score * obs_current);
          }
        } else {
          // Blume-Capel ordinal MRF variable
          for(int cat = 0; cat <= n_cats; cat++) {
            exponent = GroupThresholds[0] * cat;
            exponent += GroupThresholds[1] *
              (cat - baseline_category[var]) *
              (cat - baseline_category[var]);
            exponent += cat * rest_score - bound;
            denominator_prop += std::exp(exponent + cat * obs_proposed);
            denominator_curr += std::exp(exponent + cat * obs_current);
          }
        }

        // Update the pseudo-likelihood ratio
        pseudolikelihood_ratio -= std::log(denominator_prop);
        pseudolikelihood_ratio += std::log(denominator_curr);
      }
    }
  }
  return pseudolikelihood_ratio;
}

/**
 * Function: metropolis_interaction
 * Purpose: Uses the Metropolis-Hastings algorithm to sample from the full conditional distribution
 *          of a pairwise interaction parameter, updating its value and the associated proposal standard deviation.
 *
 * Inputs:
 *  - main_effects: Numeric matrix of main effects across variables and groups.
 *  - pairwise_effects: Numeric matrix of pairwise effects across variable pairs and groups.
 *  - main_effect_indices: Integer matrix mapping variables to main effect parameters.
 *  - pairwise_effect_indices: Integer matrix mapping variable pairs to pairwise interaction parameters.
 *  - projection: Numeric matrix representing group-specific scaling factors for interactions.
 *  - observations: Integer matrix of observed data (individuals x variables).
 *  - num_groups: Number of groups in the analysis.
 *  - group_indices: Integer matrix indicating the start and end indices for each group in `observations`.
 *  - num_categories: Integer matrix indicating the number of categories for each variable and group.
 *  - independent_thresholds: Boolean indicating whether thresholds are modeled independently.
 *  - num_persons: Total number of individuals in the dataset.
 *  - residual_matrix: Numeric matrix of residual effects for pseudo-likelihood calculations.
 *  - is_ordinal_variable: Logical vector indicating whether each variable is ordinal.
 *  - baseline_category: Integer vector indicating the reference categories for Blume-Capel variables.
 *  - proposal_sd_pairwise_effects: Numeric matrix of proposal standard deviations for pairwise interaction parameters.
 *  - interaction_scale: Scale parameter for the prior distribution of interaction terms.
 *  - num_variables: Total number of variables in the analysis.
 *  - rm_adaptation_rate: Robbins-Monro learning rate for adapting proposal standard deviations.
 *  - target_acceptance_rate: Target log acceptance rate for the Metropolis-Hastings updates.
 *  - t: Current iteration number in the sampling procedure.
 *  - rm_lower_bound: Lower bound for the proposal standard deviations.
 *  - rm_upper_bound: Upper bound for the proposal standard deviations.
 *
 * Outputs:
 *  - Updates `pairwise_effects` with new sampled values.
 *  - Updates `residual_matrix` with adjusted residual effects.
 *  - Updates `proposal_sd_pairwise_effects` with adjusted proposal standard deviations using Robbins-Monro updates.
 */
void metropolis_interaction(
    const arma::mat& main_effects,
    arma::mat& pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::mat& projection,
    const arma::imat& observations,
    const int num_groups,
    const arma::imat& group_indices,
    const arma::imat& num_categories,
    const bool independent_thresholds,
    const int num_persons,
    arma::mat& residual_matrix,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    arma::mat& proposal_sd_pairwise_effects,
    double interaction_scale,
    const int num_variables,
    double rm_adaptation_rate,
    double target_acceptance_rate,
    const int t,
    double rm_lower_bound,
    double rm_upper_bound
) {
  double log_acceptance_probability, current_state, proposed_state, U, state_difference;
  double exp_neg_log_t_rm_adaptation_rate = std::exp(-std::log(t) * rm_adaptation_rate);// Precompute Robbins-Monro decay term
  int int_index, obs1, obs2;

  // Iterate over all pairs of variables for interaction updates
  for(int variable1 = 0; variable1 <  num_variables - 1; variable1++) {
    for(int variable2 = variable1 + 1; variable2 <  num_variables; variable2++) {
      int_index = pairwise_effect_indices(variable1, variable2);

      // Retrieve the current state of the interaction parameter
      current_state = pairwise_effects(int_index, 0);

      // Propose a new state from a normal density centered at the current state
      proposed_state = R::rnorm(current_state, proposal_sd_pairwise_effects(int_index, 0));

      // Compute the log pseudo-likelihood ratio for proposed vs. current states
      log_acceptance_probability = log_pseudolikelihood_ratio_interaction(
        main_effects, main_effect_indices, projection, observations, num_groups,
        group_indices, num_categories, independent_thresholds, num_persons,
        variable1, variable2, proposed_state, current_state, residual_matrix,
        is_ordinal_variable, baseline_category);

      // Add prior probabilities for the interaction parameter
      log_acceptance_probability += R::dcauchy(proposed_state, 0.0, interaction_scale, true);
      log_acceptance_probability -= R::dcauchy(current_state, 0.0, interaction_scale, true);

      // Metropolis-Hastings acceptance step
      U = R::unif_rand();
      if(std::log(U) < log_acceptance_probability) {
        // Update the interaction parameter
        state_difference = proposed_state - current_state;
        pairwise_effects(int_index, 0) = proposed_state;

        // Update the residual matrix to reflect the new interaction parameter
        for (int person = 0; person < num_persons; person++) {
          obs1 = observations(person, variable1);
          obs2 = observations(person, variable2);
          residual_matrix(person, variable1) += obs2 * state_difference;
          residual_matrix(person, variable2) += obs1 * state_difference;
        }
      }

      // Robbins-Monro update to adapt the proposal standard deviation
      proposal_sd_pairwise_effects(int_index, 0) = update_with_robbins_monro(
        proposal_sd_pairwise_effects(int_index, 0), log_acceptance_probability, target_acceptance_rate, rm_lower_bound,
        rm_upper_bound, exp_neg_log_t_rm_adaptation_rate);
    }
  }
}


/**
 * Function: log_pseudolikelihood_ratio_pairwise_difference
 * Purpose: Computes the log pseudo-likelihood ratio for proposed vs. current
 *          values of a pairwise interaction parameter difference in a multi-group
 *          independent samples design.
 *
 * Inputs:
 *  - main_effects: Numeric matrix of main effects for variables across groups.
 *  - main_effect_indices: Integer matrix mapping variables to main effect parameters.
 *  - projection: Numeric matrix representing group-specific scaling factors for interactions.
 *  - observations: Integer matrix of observed data (individuals x variables).
 *  - num_groups: Number of groups in the analysis.
 *  - group_indices: Integer matrix indicating the start and end indices for each group in `observations`.
 *  - num_categories: Integer matrix indicating the number of categories for each variable and group.
 *  - independent_thresholds: Boolean indicating whether thresholds are modeled independently.
 *  - num_persons: Total number of individuals in the dataset.
 *  - variable1: Index of the first variable in the interaction pair.
 *  - variable2: Index of the second variable in the interaction pair.
 *  - h: Index of the group difference being modeled.
 *  - proposed_state: Proposed value for the interaction parameter difference.
 *  - current_state: Current value of the interaction parameter difference.
 *  - residual_matrix: Numeric matrix of residual effects for pseudo-likelihood calculations.
 *  - is_ordinal_variable: Logical vector indicating whether each variable is ordinal.
 *  - baseline_category: Integer vector indicating the reference categories for Blume-Capel variables.
 *
 * Outputs:
 *  - Returns the log pseudo-likelihood ratio for the proposed vs. current values
 *    of the pairwise interaction parameter difference.
 */
double log_pseudolikelihood_ratio_pairwise_difference(
    const arma::mat& main_effects,
    const arma::imat& main_effect_indices,
    const arma::mat& projection,
    const arma::imat& observations,
    const int num_groups,
    const arma::imat& group_indices,
    const arma::imat& num_categories,
    const bool independent_thresholds,
    const int num_persons,
    const int variable1,
    const int variable2,
    const int h,
    const double proposed_state,
    const double current_state,
    const arma::mat& residual_matrix,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category
) {
  double pseudolikelihood_ratio = 0.0;
  double delta_state = 2.0 * (proposed_state - current_state); // Change in parameter value


  double obs_proposed_p1, obs_current_p1, obs_proposed_p2, obs_current_p2;
  double obs_proposed_p, obs_current_p, denominator_prop, denominator_curr;
  double rest_score, bound, P, delta_state_group, exponent;

  int n_cats_v1, n_cats_v2, n_cats, obs_score1, obs_score2, score, var;

  // Loop over groups
  for(int gr = 0; gr < num_groups; gr++) {
    // Compute thresholds for both variables in the current group
    arma::vec GroupThresholds_v1 = compute_group_thresholds(
      variable1, is_ordinal_variable, gr, num_groups, num_categories, main_effects,
      main_effect_indices, projection, independent_thresholds);

    arma::vec GroupThresholds_v2 = compute_group_thresholds(
      variable2, is_ordinal_variable, gr, num_groups, num_categories, main_effects,
      main_effect_indices, projection, independent_thresholds);

    // Scale the delta_state by the projection factor for this group
    P = projection(gr, h);
    delta_state_group = delta_state * P;

    // Cache the number of categories for both variables
    n_cats_v1 = num_categories(variable1, gr);
    n_cats_v2 = num_categories(variable2, gr);

    // Iterate over all individuals in the current group
    for(int person = group_indices(gr, 0); person <= group_indices(gr, 1); person++) {
      // Cache observation scores and scaled terms
      obs_score1 = observations(person, variable1);
      obs_score2 = observations(person, variable2);
      obs_proposed_p1 = obs_score2 * proposed_state * P;
      obs_current_p1 = obs_score2 * current_state * P;
      obs_proposed_p2 = obs_score1 * proposed_state * P;
      obs_current_p2 = obs_score1 * current_state * P;

      // Contribution from the interaction term
      pseudolikelihood_ratio +=  obs_score1 * obs_score2 * delta_state_group;

      // Process each variable in the interaction pair
      for (int variable = 1; variable <= 2; variable++) {
        var = (variable == 1) ? variable1 : variable2;
        n_cats = (variable == 1) ? n_cats_v1 : n_cats_v2;
        arma::vec& GroupThresholds = (variable == 1) ? GroupThresholds_v1 : GroupThresholds_v2;
        obs_proposed_p = (variable == 1) ? obs_proposed_p1 : obs_proposed_p2;
        obs_current_p = (variable == 1) ? obs_current_p1 : obs_current_p2;

        // Compute the rest score and bound
        rest_score = residual_matrix(person, var) - obs_current_p;
        bound = (rest_score > 0) ? n_cats * rest_score : 0.0;

        // Compute denominators based on whether the variable is ordinal
        if (is_ordinal_variable[var]) {
          // Binary or ordinal MRF variable
          denominator_prop = std::exp(-bound);
          denominator_curr = std::exp(-bound);
          for (int cat = 0; cat < n_cats; cat++) {
            score = cat + 1;
            exponent = GroupThresholds[cat] + score * rest_score - bound;
            denominator_prop += std::exp(exponent + score * obs_proposed_p);
            denominator_curr += std::exp(exponent + score * obs_current_p);
          }
        } else {
          // Blume-Capel ordinal MRF variable
          denominator_prop = 0.0;
          denominator_curr = 0.0;
          for (int cat = 0; cat <= n_cats; cat++) {
            exponent = GroupThresholds[0] * cat +
              GroupThresholds[1] * (cat - baseline_category[var]) *
              (cat - baseline_category[var]) +
              cat * rest_score - bound;
            denominator_prop += std::exp(exponent + cat * obs_proposed_p);
            denominator_curr += std::exp(exponent + cat * obs_current_p);
          }
        }

        // Update the pseudo-likelihood ratio
        pseudolikelihood_ratio -= std::log(denominator_prop);
        pseudolikelihood_ratio += std::log(denominator_curr);
      }
    }
  }
  return pseudolikelihood_ratio;
}


/**
 * Function: metropolis_pairwise_difference
 * Purpose: Uses the Metropolis-Hastings algorithm to sample from the full conditional
 *          distribution of pairwise interaction parameter differences for multiple groups,
 *          updating their values and the associated proposal standard deviations.
 *
 * Inputs:
 *  - main_effects: Numeric matrix of main effects across variables and groups.
 *  - pairwise_effects: Numeric matrix of pairwise interaction parameter differences across groups.
 *  - main_effect_indices: Integer matrix mapping variables to main effect parameters.
 *  - pairwise_effect_indices: Integer matrix mapping variable pairs to pairwise interaction parameters.
 *  - projection: Numeric matrix representing group-specific scaling factors for interactions.
 *  - observations: Integer matrix of observed data (individuals x variables).
 *  - num_groups: Number of groups in the analysis.
 *  - group_indices: Integer matrix indicating the start and end indices for each group in `observations`.
 *  - num_categories: Integer matrix indicating the number of categories for each variable and group.
 *  - independent_thresholds: Boolean indicating whether thresholds are modeled independently.
 *  - inclusion_indicator: Integer matrix indicating whether pairwise differences are included in the model.
 *  - num_persons: Total number of individuals in the dataset.
 *  - residual_matrix: Numeric matrix of residual effects for pseudo-likelihood calculations.
 *  - is_ordinal_variable: Logical vector indicating whether each variable is ordinal.
 *  - baseline_category: Integer vector indicating the reference categories for Blume-Capel variables.
 *  - proposal_sd_pairwise_effects: Numeric matrix of proposal standard deviations for pairwise interaction parameters.
 *  - pairwise_difference_scale: Scale parameter for the prior distribution of pairwise differences.
 *  - num_variables: Total number of variables in the analysis.
 *  - rm_adaptation_rate: Robbins-Monro learning rate for adapting proposal standard deviations.
 *  - target_acceptance_rate: Target log acceptance rate for the Metropolis-Hastings updates.
 *  - t: Current iteration number in the sampling procedure.
 *  - rm_lower_bound: Lower bound for the proposal standard deviations.
 *  - rm_upper_bound: Upper bound for the proposal standard deviations.
 *
 * Outputs:
 *  - Updates `pairwise_effects` with new sampled values.
 *  - Updates `residual_matrix` with adjusted residual effects.
 *  - Updates `proposal_sd_pairwise_effects` with adjusted proposal standard deviations using Robbins-Monro updates.
 */
void metropolis_pairwise_difference(
    const arma::mat& main_effects,
    arma::mat& pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::mat& projection,
    const arma::imat& observations,
    const int num_groups,
    const arma::imat& group_indices,
    const arma::imat& num_categories,
    const bool independent_thresholds,
    const arma::imat& inclusion_indicator,
    const int num_persons,
    arma::mat& residual_matrix,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    arma::mat& proposal_sd_pairwise_effects,
    const double pairwise_difference_scale,
    const int num_variables,
    const double rm_adaptation_rate,
    const double target_acceptance_rate,
    const int t,
    const double rm_lower_bound,
    const double rm_upper_bound
) {
  double exp_neg_log_t_rm_adaptation_rate = std::exp(-std::log(t) * rm_adaptation_rate); // Precompute Robbins-Monro decay term
  int int_index, obs1, obs2;
  double current_state, proposed_state,log_acceptance_probability, U;
  double state_difference;

  // Iterate over all variable pairs
  for(int variable1 = 0; variable1 <  num_variables - 1; variable1++) {
    for(int variable2 = variable1 + 1; variable2 <  num_variables; variable2++) {
      if (inclusion_indicator(variable1, variable2) == 1) {
        int_index = pairwise_effect_indices(variable1, variable2);

        // Iterate over all groups (excluding the reference group)
        for(int h = 1; h < num_groups; h++) {
          // Retrieve the current state and propose a new state
          current_state = pairwise_effects(int_index, h);
          proposed_state = R::rnorm(current_state, proposal_sd_pairwise_effects(int_index, h));

          // Compute the log pseudo-likelihood ratio for the proposed vs. current state
          log_acceptance_probability = log_pseudolikelihood_ratio_pairwise_difference(
            main_effects, main_effect_indices, projection, observations, num_groups,
            group_indices, num_categories, independent_thresholds, num_persons,
            variable1, variable2, h - 1, proposed_state, current_state,
            residual_matrix, is_ordinal_variable, baseline_category);

          // Add prior probabilities for the pairwise difference parameter
          log_acceptance_probability += R::dcauchy(proposed_state, 0.0, pairwise_difference_scale, true);
          log_acceptance_probability -= R::dcauchy(current_state, 0.0, pairwise_difference_scale, true);

          // Metropolis-Hastings acceptance step
          U = R::unif_rand();
          if(std::log(U) < log_acceptance_probability) {
            pairwise_effects(int_index, h) = proposed_state;

            // Update the rest matrix to reflect the new pairwise difference
            for(int gr = 0; gr < num_groups; gr++) {

              state_difference = (proposed_state - current_state) * projection(gr, h - 1);
              for(int person = group_indices(gr, 0); person <= group_indices(gr, 1); person++) {
                obs1 = observations(person, variable1);
                obs2 = observations(person, variable2);
                residual_matrix(person, variable1) += obs2 * state_difference;
                residual_matrix(person, variable2) += obs1 * state_difference;
              }
            }
          }

          // Robbins-Monro update for the proposal standard deviation
          proposal_sd_pairwise_effects(int_index, h) = update_with_robbins_monro(
            proposal_sd_pairwise_effects(int_index, h), log_acceptance_probability, target_acceptance_rate, rm_lower_bound,
            rm_upper_bound, exp_neg_log_t_rm_adaptation_rate);
        }
      }
    }
  }
}


/**
 * Function: log_pseudolikelihood_ratio_pairwise_differences
 * Purpose:
 *   Computes the log pseudo-likelihood ratio for proposed vs. current pairwise
 *   difference states in a Bayesian ANOVA.
 *
 * Inputs:
 *   - main_effects: arma::mat of main effects for all variables and groups.
 *   - main_effect_indices: arma::imat mapping variables to category indices.
 *   - projection: arma::mat representing group-specific scaling.
 *   - observations: arma::imat of observed data (individuals by variables).
 *   - num_groups: Total number of groups in the analysis.
 *   - group_indices: arma::imat specifying group-wise start and end indices for individuals.
 *   - num_categories: arma::imat containing category counts for each variable and group.
 *   - independent_thresholds: Boolean flag for whether thresholds are modeled independently.
 *   - num_persons: Total number of individuals in the analysis.
 *   - variable1: Index of the first variable in the pair.
 *   - variable2: Index of the second variable in the pair.
 *   - proposed_states: arma::vec of proposed pairwise difference states for all groups.
 *   - current_states: arma::vec of current pairwise difference states for all groups.
 *   - residual_matrix: arma::mat of residuals used for pseudo-likelihood calculations.
 *   - is_ordinal_variable: arma::uvec indicating whether variables are ordinal.
 *   - baseline_category: arma::ivec specifying reference categories for each variable.
 *
 * Outputs:
 *   - A double representing the log pseudo-likelihood ratio for the proposed
 *     vs. current states.
 *
 * Details:
 *   - Iterates over all groups and individuals to compute contributions to the
 *     pseudo-likelihood ratio.
 *   - Handles ordinal and Blume-Capel variables separately, adjusting for their
 *     specific modeling requirements.
 */
double log_pseudolikelihood_ratio_pairwise_differences(
    const arma::mat& main_effects,
    const arma::imat& main_effect_indices,
    const arma::mat& projection,
    const arma::imat& observations,
    const int num_groups,
    const arma::imat& group_indices,
    const arma::imat& num_categories,
    const bool independent_thresholds,
    const int num_persons,
    const int variable1,
    const int variable2,
    const arma::vec& proposed_states,
    const arma::vec& current_states,
    const arma::mat& residual_matrix,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category
) {

  double pseudolikelihood_ratio = 0.0;

  // Loop over groups
  for (int gr = 0; gr < num_groups; gr++) {
    // Compute thresholds for both variables
    arma::vec GroupThresholds_v1 = compute_group_thresholds(
      variable1, is_ordinal_variable, gr, num_groups, num_categories, main_effects,
      main_effect_indices, projection, independent_thresholds);

    arma::vec GroupThresholds_v2 = compute_group_thresholds(
      variable2, is_ordinal_variable, gr, num_groups, num_categories, main_effects,
      main_effect_indices, projection, independent_thresholds);

    // Compute current and proposed group interactions
    double current_group_interaction = 0.0;
    double proposed_group_interaction = 0.0;
    for (int h = 0; h < num_groups - 1; h++) {
      current_group_interaction += current_states[h] * projection(gr, h);
      proposed_group_interaction += proposed_states[h] * projection(gr, h);
    }
    double delta_state_group = 2.0 * (proposed_group_interaction - current_group_interaction);

    // Cache the number of categories for both variables
    int n_cats_v1 = num_categories(variable1, gr);
    int n_cats_v2 = num_categories(variable2, gr);

    // Iterate over all individuals in the current group
    for (int person = group_indices(gr, 0); person <= group_indices(gr, 1); person++) {
      int obs_score1 = observations(person, variable1);
      int obs_score2 = observations(person, variable2);
      double obs_proposed_p1 = obs_score2 * proposed_group_interaction;
      double obs_current_p1 = obs_score2 * current_group_interaction;
      double obs_proposed_p2 = obs_score1 * proposed_group_interaction;
      double obs_current_p2 = obs_score1 * current_group_interaction;

      // Contribution from the interaction term
      pseudolikelihood_ratio += obs_score1 * obs_score2 * delta_state_group;

      // Process each variable in the interaction pair
      for (int variable = 1; variable <= 2; variable++) {
        int var = (variable == 1) ? variable1 : variable2;
        int n_cats = (variable == 1) ? n_cats_v1 : n_cats_v2;
        arma::vec& GroupThresholds = (variable == 1) ? GroupThresholds_v1 : GroupThresholds_v2;
        double obs_proposed_p = (variable == 1) ? obs_proposed_p1 : obs_proposed_p2;
        double obs_current_p = (variable == 1) ? obs_current_p1 : obs_current_p2;

        double rest_score = residual_matrix(person, var) - obs_current_p;
        double bound = (rest_score > 0) ? n_cats * rest_score : 0.0;

        double denominator_curr = 0.0, denominator_prop = 0.0;

        // Compute denominators
        if (is_ordinal_variable[var]) {
          denominator_prop += std::exp(-bound);
          denominator_curr += std::exp(-bound);

          for (int cat = 0; cat < n_cats; cat++) {
            int score = cat + 1;
            double exponent = GroupThresholds[cat] + rest_score * score - bound;
            denominator_curr += std::exp(exponent + score * obs_current_p);
            denominator_prop += std::exp(exponent + score * obs_proposed_p);
          }
        } else {
          for (int cat = 0; cat <= n_cats; cat++) {
            double exponent = GroupThresholds[0] * cat +
              GroupThresholds[1] * (cat - baseline_category[var]) *
              (cat - baseline_category[var]) +
              rest_score * cat - bound;
            denominator_curr += std::exp(exponent + cat * obs_current_p);
            denominator_prop += std::exp(exponent + cat * obs_proposed_p);
          }
        }

        // Update log pseudo-likelihood ratio
        pseudolikelihood_ratio += std::log(denominator_curr) - std::log(denominator_prop);
      }
    }
  }

  return pseudolikelihood_ratio;
}


/**
 * Function: metropolis_pairwise_difference_between_model
 * Purpose:
 *   Implements a between-model Metropolis-Hastings algorithm to update
 *   the inclusion inclusion_indicator and pairwise differences for variable pairs.
 *
 * Inputs:
 *   - inclusion_probability_difference: arma::mat of inclusion probabilities
 *                                       for pairwise differences.
 *   - index: arma::imat mapping pairwise differences to variable indices.
 *   - main_effects: arma::mat of main effects for all variables and groups.
 *   - pairwise_effects: arma::mat of pairwise effects for all variable pairs and groups.
 *   - main_effect_indices: arma::imat mapping variables to category indices.
 *   - pairwise_effect_indices: arma::imat mapping variable pairs to pairwise effect indices.
 *   - projection: arma::mat for group-specific scaling.
 *   - observations: arma::imat of observed data (individuals by variables).
 *   - num_groups: Total number of groups in the analysis.
 *   - group_indices: arma::imat specifying group-wise start and end indices for individuals.
 *   - num_categories: arma::imat containing category counts for each variable and group.
 *   - independent_thresholds: Boolean flag for whether thresholds are modeled independently.
 *   - inclusion_indicator: arma::imat indicating active pairwise differences.
 *   - num_persons: Total number of individuals in the analysis.
 *   - residual_matrix: arma::mat of residuals used for pseudo-likelihood calculations.
 *   - is_ordinal_variable: arma::uvec indicating whether variables are ordinal.
 *   - baseline_category: arma::ivec specifying reference categories for each variable.
 *   - proposal_sd_pairwise_effects: arma::mat of proposal standard deviations for pairwise differences.
 *   - pairwise_difference_scale: Double representing the scale of the prior distribution
 *                                for pairwise differences.
 *   - num_variables: Total number of variables in the analysis.
 *   - rm_adaptation_rate: Robbins-Monro learning rate parameter.
 *   - target_acceptance_rate: Target acceptance probability for Metropolis-Hastings updates.
 *   - t: Current iteration number.
 *   - rm_lower_bound: Lower bound for proposal standard deviations.
 *   - rm_upper_bound: Upper bound for proposal standard deviations.
 *   - num_pairwise: Total number of pairwise differences.
 *
 * Outputs:
 *   - Updates `inclusion_indicator`, `pairwise_effects`, and `residual_matrix` to reflect
 *     accepted proposals for pairwise differences.
 */
void metropolis_pairwise_difference_between_model(
    const arma::mat& inclusion_probability_difference,
    const arma::imat& index,
    const arma::mat& main_effects,
    arma::mat& pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::mat& projection,
    const arma::imat& observations,
    const int num_groups,
    const arma::imat& group_indices,
    const arma::imat& num_categories,
    const bool independent_thresholds,
    arma::imat& inclusion_indicator,
    const int num_persons,
    arma::mat& residual_matrix,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    arma::mat& proposal_sd_pairwise_effects,
    const double pairwise_difference_scale,
    const int num_variables,
    const double rm_adaptation_rate,
    const double target_acceptance_rate,
    const int t,
    const double rm_lower_bound,
    const double rm_upper_bound,
    const int num_pairwise
) {
  // Vectors to store current and proposed states for pairwise differences
  arma::vec proposed_states(num_groups - 1);
  arma::vec current_states(num_groups - 1);

  // Loop over all pairwise differences
  for (int cntr = 0; cntr < num_pairwise; cntr++) {
    int variable1 = index(cntr, 1);
    int variable2 = index(cntr, 2);
    int int_index = pairwise_effect_indices(variable1, variable2);

    double log_acceptance_probability = 0.0;

    // Loop over groups to process current and proposed states
    for (int h = 1; h < num_groups; h++) {
      double current_state = pairwise_effects(int_index, h);
      current_states[h - 1] = current_state;

      double proposed_state = 0.0;

      // Update log probabilities based on the inclusion inclusion_indicator
      if (inclusion_indicator(variable1, variable2) == 1) {
        // Difference is included
        log_acceptance_probability -= R::dcauchy(current_state, 0.0, pairwise_difference_scale, true);
        log_acceptance_probability += R::dnorm(current_state, proposed_state, proposal_sd_pairwise_effects(int_index, h), true);
      } else {
        // Propose a new state
        proposed_state = R::rnorm(current_state, proposal_sd_pairwise_effects(int_index, h));
        log_acceptance_probability += R::dcauchy(proposed_state, 0.0, pairwise_difference_scale, true);
        log_acceptance_probability -= R::dnorm(proposed_state, current_state, proposal_sd_pairwise_effects(int_index, h), true);
      }

      proposed_states[h - 1] = proposed_state;
    }

    // Update log probability for the inclusion inclusion_indicator
    if (inclusion_indicator(variable1, variable2) == 1) {
      log_acceptance_probability -= std::log(inclusion_probability_difference(variable1, variable2));
      log_acceptance_probability += std::log(1 - inclusion_probability_difference(variable1, variable2));
    } else {
      log_acceptance_probability += std::log(inclusion_probability_difference(variable1, variable2));
      log_acceptance_probability -= std::log(1 - inclusion_probability_difference(variable1, variable2));
    }

    // Compute log pseudo-likelihood ratio
    log_acceptance_probability += log_pseudolikelihood_ratio_pairwise_differences(
      main_effects, main_effect_indices, projection, observations,
      num_groups, group_indices, num_categories, independent_thresholds, num_persons, variable1,
      variable2, proposed_states, current_states, residual_matrix, is_ordinal_variable, baseline_category);

    // Metropolis-Hastings acceptance step
    double U = R::unif_rand();
    if (std::log(U) < log_acceptance_probability) {
      // Update inclusion inclusion_indicator
      inclusion_indicator(variable1, variable2) = 1 - inclusion_indicator(variable1, variable2);
      inclusion_indicator(variable2, variable1) = inclusion_indicator(variable1, variable2);

      // Update pairwise effects and rest matrix
      for (int h = 1; h < num_groups; h++) {
        pairwise_effects(int_index, h) = proposed_states[h - 1];
      }

      // Update residuals in the rest matrix
      for (int gr = 0; gr < num_groups; ++gr) {
        double state_difference = 0.0;
        for (int h = 0; h < num_groups - 1; h++) {
          state_difference += (proposed_states[h] - current_states[h]) * projection(gr, h);
        }

        for (int person = group_indices(gr, 0); person <= group_indices(gr, 1); person++) {
          int obs1 = observations(person, variable1);
          int obs2 = observations(person, variable2);
          residual_matrix(person, variable1) += obs2 * state_difference;
          residual_matrix(person, variable2) += obs1 * state_difference;
        }
      }
    }
  }
}


/**
 * Function: metropolis_threshold_regular
 * Purpose: Uses the Metropolis-Hastings algorithm to sample from the full conditional
 *          distribution of overall category threshold parameters for regular ordinal variables.
 *
 * Inputs:
 *  - main_effects: Numeric matrix containing the main effect parameters for all variables and groups.
 *  - main_effect_indices: Integer matrix mapping variables to their main effect parameter indices.
 *  - projection: Numeric matrix specifying group-specific scaling for thresholds.
 *  - observations: Integer matrix containing the observed data (individuals x variables).
 *  - num_groups: Total number of groups in the analysis.
 *  - group_indices: Integer matrix specifying start and end indices for each group in `observations`.
 *  - num_categories: Integer matrix specifying the number of categories for each variable and group.
 *  - num_persons: Total number of individuals in the dataset.
 *  - residual_matrix: Numeric matrix storing residual effects for pseudo-likelihood calculations.
 *  - n_cat_obs: List of integer matrices, where each matrix tracks the number of observations
 *               for each category of a variable in a specific group.
 *  - prior_threshold_alpha: Hyperparameter for the prior distribution of threshold parameters (shape parameter).
 *  - prior_threshold_beta: Hyperparameter for the prior distribution of threshold parameters (rate parameter).
 *  - variable: Index of the variable whose thresholds are being updated.
 *
 * Outputs:
 *  - Updates the `main_effects` matrix with new sampled threshold parameters for the specified variable.
 */
void metropolis_threshold_regular(
    arma::mat& main_effects,
    const arma::imat& main_effect_indices,
    const arma::mat& projection,
    const arma::imat& observations,
    const int num_groups,
    const arma::imat& group_indices,
    const arma::imat& num_categories,
    const int num_persons,
    const arma::mat& residual_matrix,
    const List& n_cat_obs,
    const double prior_threshold_alpha,
    const double prior_threshold_beta,
    const int variable
) {
  arma::vec q(num_persons); // Intermediate storage for pseudo-likelihood calculations
  arma::vec r(num_persons); // Intermediate storage for pseudo-likelihood calculations

  // Cache the number of categories and main effect index for the variable
  int n_cats = num_categories(variable, 0); // Number of categories for the variable
  int cat_index = main_effect_indices(variable, 0); // Index in `main_effects` for this variable

  arma::vec GroupThresholds(n_cats); // Vector to store group-specific thresholds

  // Iterate over each category threshold for the variable
  for(int category = 0; category < n_cats; category++) {
    double current_state = main_effects(cat_index + category, 0); // Current state of the threshold
    double exp_current = std::exp(current_state); // Exponentiated current threshold
    double c = (prior_threshold_alpha + prior_threshold_beta) / (1 + exp_current); // Initial value for c

    // Compute group-specific thresholds and contributions to `q` and `r`
    for (int gr = 0; gr < num_groups; gr++) {
      // Update thresholds for the current group
      for (int cat = 0; cat < n_cats; cat++) {
        double threshold = main_effects(cat_index + cat, 0);
        for (int h = 1; h < num_groups; h++) {
          threshold += projection(gr, h - 1) * main_effects(cat_index + cat, h);
        }
        GroupThresholds[cat] = threshold;
      }

      // Subtract the current category's base threshold
      GroupThresholds[category] -= main_effects(cat_index + category, 0);

      // Compute `q` and `r` for each person in the group
      for (int person = group_indices(gr, 0); person <= group_indices(gr, 1); person++) {
        double rest_score = residual_matrix(person, variable);
        double q_person = 1.0; // Initialize q for the person
        for (int cat = 0; cat < n_cats; ++cat) {
          if (cat != category) {
            double exponent = GroupThresholds[cat] + (cat + 1) * rest_score;
            q_person += std::exp(exponent);
          }
        }
        double exponent_r = GroupThresholds[category] + (category + 1) * rest_score;
        double r_person = std::exp(exponent_r);
        q[person] = q_person;
        r[person] = r_person;
        c += r_person / (q_person + r_person * exp_current);
      }
    }

    // Update `c` to include the prior's contribution
    double tmp = num_persons + prior_threshold_alpha + prior_threshold_beta - exp_current * c;
    c /= tmp;

    // Generate a proposed state using a generalized beta-prime proposal
    double a = prior_threshold_alpha;
    double b = num_persons + prior_threshold_beta;
    for(int gr = 0; gr < num_groups; gr++) {
      arma::imat n_cat_obs_gr = n_cat_obs[gr];
      a += n_cat_obs_gr(category + 1, variable); // Update a with observations in the category
      b -= n_cat_obs_gr(category + 1, variable); // Update b with remaining observations
    }

    // Sample from the beta distribution
    double tmp_beta = R::rbeta(a, b);
    double proposed_state = std::log(tmp_beta / (1  - tmp_beta) / c);
    double exp_proposed = std::exp(proposed_state);

    // Compute the log acceptance probability
    double log_acceptance_probability = 0.0;

    // Compute pseudo-likelihood ratio
    for (int gr = 0; gr < num_groups; ++gr) {
      for (int person = group_indices(gr, 0); person <= group_indices(gr, 1); ++person) {
        log_acceptance_probability += std::log(q[person] + r[person] * exp_current);
        log_acceptance_probability -= std::log(q[person] + r[person] * exp_proposed);
      }
    }

    // Add prior density ratio
    log_acceptance_probability -= (prior_threshold_alpha + prior_threshold_beta) * std::log(1 + exp_proposed);
    log_acceptance_probability += (prior_threshold_alpha + prior_threshold_beta) * std::log(1 + exp_current);

    // Add proposal density ratio
    log_acceptance_probability -= (a + b) * std::log(1 + c * exp_current);
    log_acceptance_probability += (a + b) * std::log(1 + c * exp_proposed);

    // Perform Metropolis-Hastings acceptance step
    double U = std::log(R::unif_rand());
    if(U < log_acceptance_probability) {
      main_effects(cat_index + category, 0) = proposed_state; // Accept the proposal
    }
  }
}


/**
 * Function: log_pseudolikelihood_ratio_main_difference
 * Purpose: Computes the log pseudo-likelihood ratio for a proposed change in
 *          a category threshold difference for an independent samples design.
 *
 * Inputs:
 *  - main_effects: Numeric matrix containing the main effect parameters for all variables and groups.
 *  - main_effect_indices: Integer matrix mapping variables to their main effect parameter indices.
 *  - projection: Numeric matrix specifying group-specific scaling for thresholds.
 *  - observations: Integer matrix containing the observed data (individuals x variables).
 *  - num_groups: Total number of groups in the analysis.
 *  - group_indices: Integer matrix specifying start and end indices for each group in `observations`.
 *  - num_categories: Integer matrix specifying the number of categories for each variable and group.
 *  - num_persons: Total number of individuals in the dataset.
 *  - residual_matrix: Numeric matrix storing residual effects for pseudo-likelihood calculations.
 *  - n_cat_obs: List of integer matrices tracking the number of observations for each category
 *               of a variable in a specific group.
 *  - variable: Index of the variable whose category threshold is being updated.
 *  - category: Index of the specific category threshold being updated for the variable.
 *  - h: Index for the projection matrix, corresponding to the group difference.
 *  - proposed_state: Proposed value for the category threshold difference.
 *  - current_state: Current value of the category threshold difference.
 *
 * Outputs:
 *  - Returns the log pseudo-likelihood ratio for the proposed versus current state.
 */
double log_pseudolikelihood_ratio_main_difference(
    const arma::mat& main_effects,
    const arma::imat& main_effect_indices,
    const arma::mat& projection,
    const arma::imat& observations,
    const int num_groups,
    const arma::imat& group_indices,
    const arma::imat& num_categories,
    const int num_persons,
    const arma::mat& residual_matrix,
    const List& n_cat_obs,
    const int variable,
    const int category,
    const int h,
    double proposed_state,
    double current_state
) {
  double pseudolikelihood_ratio = 0.0; // Initialize the log pseudo-likelihood ratio
  double delta_state = proposed_state - current_state; // Difference between proposed and current states
  int cat_index = main_effect_indices(variable, 0); // Index for the main effects of the variable


  // Loop over all groups
  for(int gr = 0; gr < num_groups; gr++) {
    int n_cats = num_categories(variable, gr); // Number of categories for the variable in this group
    arma::vec current_thresholds(n_cats); // Store current thresholds for the group
    arma::vec proposed_thresholds(n_cats); // Store proposed thresholds for the group
    double P = projection(gr, h); // Group-specific projection scaling factor

    // Compute current and proposed thresholds for all categories
    for(int cat = 0; cat < n_cats; cat++) {
      int full_cat_index = cat_index + cat;
      double threshold = main_effects(full_cat_index, 0);
      for (int hh = 1; hh < num_groups; ++hh) {
        threshold += projection(gr, hh - 1) * main_effects(full_cat_index, hh);
      }
      current_thresholds[cat] = threshold;
      proposed_thresholds[cat] = threshold;
    }

    // Adjust the threshold for the specific category
    proposed_thresholds[category] -= P * current_state;
    proposed_thresholds[category] += P * proposed_state;

    // Add the contribution from delta_state based on observations
    arma::imat n_cat_obs_gr = n_cat_obs[gr];
    pseudolikelihood_ratio += delta_state * P * n_cat_obs_gr(category + 1, variable);

    // Loop over all persons in the group
    for(int person = group_indices(gr, 0); person <= group_indices(gr, 1); person++) {
      double rest_score = residual_matrix(person, variable); // Compute residual score
      double bound = (rest_score > 0) ? n_cats * rest_score : 0.0;

      // Compute the denominators for proposed and current thresholds
      double denominator_proposed = std::exp(-bound);
      double denominator_current = std::exp(-bound);
      for(int cat = 0; cat < n_cats; cat++) {
        double exponent = (cat + 1) * rest_score - bound;
        denominator_proposed += std::exp(exponent + proposed_thresholds[cat]);
        denominator_current += std::exp(exponent + current_thresholds[cat]);
      }

      // Update the pseudo-likelihood ratio with log-likelihood differences
      pseudolikelihood_ratio -= std::log(denominator_proposed);
      pseudolikelihood_ratio += std::log(denominator_current);
    }
  }

  return pseudolikelihood_ratio;
}


/**
 * Function: metropolis_main_difference_regular
 * Purpose: Implements the Metropolis-Hastings (MH) algorithm to sample from
 *          the full conditional of the category threshold difference parameter
 *          for regular ordinal variables in an independent samples design.
 *
 * Inputs:
 *  - main_effects: Numeric matrix containing the main effect parameters for all variables and groups.
 *  - main_effect_indices: Integer matrix mapping variables to their main effect parameter indices.
 *  - projection: Numeric matrix specifying group-specific scaling for thresholds.
 *  - observations: Integer matrix containing the observed data (individuals x variables).
 *  - num_groups: Total number of groups in the analysis.
 *  - group_indices: Integer matrix specifying start and end indices for each group in `observations`.
 *  - num_categories: Integer matrix specifying the number of categories for each variable and group.
 *  - num_persons: Total number of individuals in the dataset.
 *  - residual_matrix: Numeric matrix storing residual effects for pseudo-likelihood calculations.
 *  - n_cat_obs: List of integer matrices tracking the number of observations for each category
 *               of a variable in a specific group.
 *  - variable: Index of the variable being updated.
 *  - inclusion_indicator: Integer matrix indicating active variables for the analysis.
 *  - proposal_sd_main_effects: Numeric matrix specifying proposal standard deviations for category thresholds.
 *  - main_difference_scale: Scale parameter for the Cauchy prior on threshold differences.
 *  - rm_adaptation_rate: Robbins-Monro learning rate for adaptive proposal standard deviations.
 *  - target_acceptance_rate: Target log acceptance rate for Metropolis-Hastings updates.
 *  - t: Current iteration number of the sampler.
 *  - rm_lower_bound: Minimum allowable standard deviation for proposals.
 *  - rm_upper_bound: Maximum allowable standard deviation for proposals.
 *
 * Outputs:
 *  - Updates `main_effects` with sampled values for threshold differences.
 *  - Updates `proposal_sd_main_effects` with adjusted proposal standard deviations using Robbins-Monro updates.
 */
void metropolis_main_difference_regular(
    arma::mat& main_effects,
    const arma::imat& main_effect_indices,
    const arma::mat& projection,
    const arma::imat& observations,
    const int num_groups,
    const arma::imat& group_indices,
    const arma::imat& num_categories,
    const int num_persons,
    const arma::mat& residual_matrix,
    const List& n_cat_obs,
    const int variable,
    const arma::imat& inclusion_indicator,
    arma::mat& proposal_sd_main_effects,
    const double main_difference_scale,
    const double rm_adaptation_rate,
    const double target_acceptance_rate,
    const int t,
    const double rm_lower_bound,
    const double rm_upper_bound
) {
  // Precompute Robbins-Monro term
  double exp_neg_log_t_rm_adaptation_rate = std::exp(-std::log(t) * rm_adaptation_rate);

  // Check if the variable is active
  if (inclusion_indicator(variable, variable) != 1) {
    return; // Skip if variable is inactive
  }

  // Look up and cache the base category index for this variable
  int base_cat_index = main_effect_indices(variable, 0);

  // Loop over all categories for this variable
  int n_cats = num_categories(variable, 0);
  for (int category = 0; category < n_cats; category++) {
    int cat_index = base_cat_index + category;

    // Loop over groups (starting from h = 1)
    for (int h = 1; h < num_groups; h++) {
      double current_state = main_effects(cat_index, h); // Current threshold difference
      double proposed_state = R::rnorm(current_state, proposal_sd_main_effects(cat_index, h)); // Propose a new state

      // Compute log pseudo-likelihood ratio for proposed vs current state
      double log_acceptance_probability = log_pseudolikelihood_ratio_main_difference(
        main_effects, main_effect_indices, projection, observations, num_groups,
        group_indices, num_categories, num_persons, residual_matrix, n_cat_obs,
        variable, category, h - 1, proposed_state, current_state);

      // Add contributions from the Cauchy prior
      log_acceptance_probability += R::dcauchy(proposed_state, 0.0, main_difference_scale, true);
      log_acceptance_probability -= R::dcauchy(current_state, 0.0, main_difference_scale, true);

      // Metropolis-Hastings acceptance step
      double U = R::unif_rand();
      if (std::log(U) < log_acceptance_probability) {
        main_effects(cat_index, h) = proposed_state; // Accept the proposed state
      }

      // Robbins-Monro update to the proposal standard deviation
      proposal_sd_main_effects(cat_index, h) = update_with_robbins_monro(
        proposal_sd_main_effects(cat_index, h), log_acceptance_probability, target_acceptance_rate,
        rm_lower_bound, rm_upper_bound, exp_neg_log_t_rm_adaptation_rate);
    }
  }
}


/**
 * Function: log_pseudolikelihood_ratio_thresholds_blumecapel
 * Purpose: Computes the log pseudo-likelihood ratio for proposed versus current
 *          category threshold parameters in the Blume-Capel model for an independent samples design.
 *
 * Inputs:
 *  - linear_current: Current value of the linear coefficient for the threshold.
 *  - quadratic_current: Current value of the quadratic coefficient for the threshold.
 *  - linear_proposed: Proposed value of the linear coefficient for the threshold.
 *  - quadratic_proposed: Proposed value of the quadratic coefficient for the threshold.
 *  - variable: Index of the variable for which thresholds are being updated.
 *  - baseline_category: Integer vector specifying the reference category for each variable.
 *  - main_effects: Numeric matrix of main effect parameters for all variables and groups.
 *  - main_effect_indices: Integer matrix mapping variables to their main effect parameter indices.
 *  - projection: Numeric matrix specifying group-specific scaling for thresholds.
 *  - sufficient_blume_capel: List of integer matrices containing sufficient statistics
 *                            for the Blume-Capel model, including linear and quadratic terms.
 *  - num_persons: Total number of individuals in the dataset.
 *  - num_groups: Total number of groups in the analysis.
 *  - group_indices: Integer matrix specifying start and end indices for each group in `observations`.
 *  - residual_matrix: Numeric matrix storing residual effects for pseudo-likelihood calculations.
 *  - num_categories: Integer matrix specifying the number of categories for each variable and group.
 *
 * Outputs:
 *  - Returns the computed log pseudo-likelihood ratio for the proposed versus current parameters.
 */
double log_pseudolikelihood_ratio_thresholds_blumecapel(
    double linear_current,
    double quadratic_current,
    double linear_proposed,
    double quadratic_proposed,
    const int variable,
    const arma::ivec& baseline_category,
    const arma::mat& main_effects,
    const arma::imat& main_effect_indices,
    const arma::mat& projection,
    const List& sufficient_blume_capel,
    const int num_persons,
    const int num_groups,
    const arma::imat& group_indices,
    const arma::mat& residual_matrix,
    const arma::imat& num_categories
) {
  // Variables for bounds, scores, and likelihood components
  double lbound, bound;
  double rest_score, numerator, denominator, exponent;
  double pseudolikelihood_ratio = 0.0;
  int linear_score, quadratic_score;
  int cat_index = main_effect_indices(variable, 0);

  // Loop over all groups
  for(int gr = 0; gr < num_groups; gr++) {
    // Precompute constant terms for the numerator and denominator
    arma::vec constant_numerator (num_categories(variable, gr) + 1);
    arma::vec constant_denominator (num_categories(variable, gr) + 1);

    for(int category = 0; category <= num_categories(variable, gr); category++) {
      linear_score = category;
      quadratic_score = (category - baseline_category[variable]) *
        (category - baseline_category[variable]);

      // Linear and quadratic contributions for current and proposed states
      constant_numerator[category] = linear_current * linear_score;
      constant_numerator[category] += quadratic_current * quadratic_score;
      constant_denominator[category] = linear_proposed * linear_score ;
      constant_denominator[category] += quadratic_proposed * quadratic_score;

      // Add group-specific contributions
      for(int h = 1; h < num_groups; h++) {
        double P = projection(gr, h - 1);
        double m1 =  main_effects(cat_index, h);
        double m2 =  main_effects(cat_index + 1, h);
        constant_numerator[category] += P * m1 * linear_score;
        constant_numerator[category] += P * m2 * quadratic_score;
        constant_denominator[category] += P * m1 * linear_score;
        constant_denominator[category] += P * m2 * quadratic_score;
      }
    }

    // Precompute bounds for numerical stability
    double tmp_num = max(constant_numerator);
    double tmp_den = max(constant_denominator);
    if(tmp_num > 0) {
      if(tmp_num > tmp_den) {
        lbound = tmp_num;
      } else {
        lbound = tmp_den;
      }
    } else {
      lbound = 0.0;
    }

    arma::imat sufficient_blume_capel_gr = sufficient_blume_capel[gr];

    // Add contributions from sufficient statistics
    pseudolikelihood_ratio += sufficient_blume_capel_gr(0, variable) * linear_proposed;
    pseudolikelihood_ratio += sufficient_blume_capel_gr(1, variable) * quadratic_proposed;
    pseudolikelihood_ratio -= sufficient_blume_capel_gr(0, variable) * linear_current;
    pseudolikelihood_ratio -= sufficient_blume_capel_gr(1, variable) * quadratic_current;

    // Loop over individuals in the group
    for(int person = group_indices(gr, 0); person <= group_indices(gr, 1); person++) {
      rest_score = residual_matrix(person, variable);
      if(rest_score > 0) {
        bound = num_categories(variable, gr) * rest_score + lbound;
      } else {
        bound = lbound;
      }

      // Compute the likelihood contributions
      numerator = 0.0;
      denominator = 0.0;
      for(int category = 0; category <= num_categories[variable]; category ++) {
        exponent = category * rest_score - bound;
        numerator += std::exp(constant_numerator[category] + exponent);
        denominator += std::exp(constant_denominator[category] + exponent);
      }
      pseudolikelihood_ratio += std::log(numerator);
      pseudolikelihood_ratio -= std::log(denominator);
    }
  }

  return pseudolikelihood_ratio;
}


/**
 * Function: metropolis_threshold_blumecapel
 * Purpose: Samples from the full conditional distribution of the overall category
 *          threshold parameters for Blume-Capel ordinal variables using the
 *          Metropolis-Hastings (MH) algorithm with Robbins-Monro adaptive updates.
 *
 * Inputs:
 *  - main_effects: Numeric matrix of main effect parameters for all variables and groups.
 *  - main_effect_indices: Integer matrix mapping variables to their main effect parameter indices.
 *  - projection: Numeric matrix specifying group-specific scaling for thresholds.
 *  - num_categories: Integer matrix specifying the number of categories for each variable and group.
 *  - sufficient_blume_capel: List of sufficient statistics for the Blume-Capel model,
 *                            including linear and quadratic terms.
 *  - num_persons: Total number of individuals in the dataset.
 *  - num_groups: Total number of groups in the analysis.
 *  - group_indices: Integer matrix specifying start and end indices for each group in `observations`.
 *  - variable: Index of the variable for which thresholds are being updated.
 *  - baseline_category: Integer vector specifying the reference category for each variable.
 *  - prior_threshold_alpha: Hyperparameter for the Blume-Capel threshold prior distribution.
 *  - prior_threshold_beta: Hyperparameter for the Blume-Capel threshold prior distribution.
 *  - residual_matrix: Numeric matrix storing residual effects for pseudo-likelihood calculations.
 *  - proposal_sd_main_effects: Numeric matrix storing the proposal standard deviations for the parameters.
 *  - rm_adaptation_rate: Robbins-Monro learning rate.
 *  - target_acceptance_rate: Target log acceptance probability for the MH algorithm.
 *  - t: Current iteration number of the sampler.
 *  - rm_lower_bound: Lower bound for the proposal standard deviation.
 *  - rm_upper_bound: Upper bound for the proposal standard deviation.
 *
 * Outputs:
 *  - Updates `main_effects` for the linear and quadratic threshold parameters.
 *  - Updates `proposal_sd_main_effects` with Robbins-Monro adjustments for the proposal standard deviations.
 */
void metropolis_threshold_blumecapel(
    arma::mat& main_effects,
    const arma::imat& main_effect_indices,
    const arma::mat& projection,
    const arma::imat& num_categories,
    const List& sufficient_blume_capel,
    const int num_persons,
    const int num_groups,
    const arma::imat& group_indices,
    const int variable,
    const arma::ivec& baseline_category,
    const double prior_threshold_alpha,
    const double prior_threshold_beta,
    const arma::mat& residual_matrix,
    arma::mat& proposal_sd_main_effects,
    const double rm_adaptation_rate,
    const double target_acceptance_rate,
    const int t,
    const double rm_lower_bound,
    const double rm_upper_bound
) {
  // Precompute Robbins-Monro adjustment term
  double exp_neg_log_t_rm_adaptation_rate = std::exp(-std::log(t) * rm_adaptation_rate);

  // Adaptive Metropolis procedure for the linear Blume-Capel model
  int cat_index = main_effect_indices(variable, 0); // Base index for the variable
  double current_state = main_effects(cat_index, 0);
  double proposed_state = R::rnorm(current_state, proposal_sd_main_effects(cat_index, 0));

  // Compute log acceptance probability
  double linear_current = current_state;
  double quadratic_current = main_effects(cat_index + 1, 0);
  double linear_proposed = proposed_state;
  double quadratic_proposed = main_effects(cat_index + 1, 0);

  double log_acceptance_probability = log_pseudolikelihood_ratio_thresholds_blumecapel(
    linear_current, quadratic_current, linear_proposed, quadratic_proposed,
    variable, baseline_category, main_effects, main_effect_indices, projection,
    sufficient_blume_capel, num_persons, num_groups, group_indices, residual_matrix,
    num_categories);

  // Add prior ratio to log acceptance probability
  log_acceptance_probability += prior_threshold_alpha * (proposed_state - current_state);
  log_acceptance_probability += (prior_threshold_alpha + prior_threshold_beta) * std::log(1 + std::exp(current_state));
  log_acceptance_probability -= (prior_threshold_alpha + prior_threshold_beta) * std::log(1 + std::exp(proposed_state));

  // Perform Metropolis-Hastings acceptance step
  double U = R::unif_rand();
  if(std::log(U) < log_acceptance_probability) {
    main_effects(cat_index, 0) = proposed_state;
  }

  // Robbins-Monro update for proposal standard deviation
  proposal_sd_main_effects(cat_index, 0) = update_with_robbins_monro(
    proposal_sd_main_effects(cat_index, 0), log_acceptance_probability, target_acceptance_rate, rm_lower_bound,
    rm_upper_bound, exp_neg_log_t_rm_adaptation_rate);

  // Adaptive Metropolis procedure for the quadratic Blume-Capel model
  current_state = main_effects(cat_index + 1, 0);
  proposed_state = R::rnorm(current_state, proposal_sd_main_effects(cat_index + 1, 0));

  // Compute log acceptance probability
  linear_current = main_effects(cat_index, 0);
  quadratic_current = current_state;
  linear_proposed =  main_effects(cat_index, 0);
  quadratic_proposed =  proposed_state;

  log_acceptance_probability = log_pseudolikelihood_ratio_thresholds_blumecapel(
    linear_current, quadratic_current, linear_proposed, quadratic_proposed,
    variable, baseline_category, main_effects, main_effect_indices, projection,
    sufficient_blume_capel, num_persons, num_groups, group_indices, residual_matrix,
    num_categories);

  // Add prior ratio to log acceptance probability
  log_acceptance_probability += prior_threshold_alpha * (proposed_state - current_state);
  log_acceptance_probability += (prior_threshold_alpha + prior_threshold_beta) * std::log(1 + std::exp(current_state));
  log_acceptance_probability -= (prior_threshold_alpha + prior_threshold_beta) * std::log(1 + std::exp(proposed_state));

  // Perform Metropolis-Hastings acceptance step
  U = R::unif_rand();
  if(std::log(U) < log_acceptance_probability) {
    main_effects(cat_index + 1, 0) = proposed_state;
  }

  // Robbins-Monro update for proposal standard deviation
  proposal_sd_main_effects(cat_index + 1, 0) = update_with_robbins_monro(
    proposal_sd_main_effects(cat_index + 1, 0), log_acceptance_probability, target_acceptance_rate, rm_lower_bound,
    rm_upper_bound, exp_neg_log_t_rm_adaptation_rate);
}


/**
 * Function: log_pseudolikelihood_ratio_main_difference_blumecapel
 * Purpose: Computes the log pseudo-likelihood ratio (proposed against current)
 *          for the differences in two category threshold parameters in the
 *          Blume-Capel model.
 *
 * Inputs:
 *  - linear_current: Current value of the linear parameter for the Blume-Capel model.
 *  - quadratic_current: Current value of the quadratic parameter for the Blume-Capel model.
 *  - linear_proposed: Proposed value of the linear parameter for the Blume-Capel model.
 *  - quadratic_proposed: Proposed value of the quadratic parameter for the Blume-Capel model.
 *  - variable: Index of the variable for which thresholds are being updated.
 *  - h: Index of the group being updated in the multi-group design.
 *  - baseline_category: Integer vector specifying the reference category for each variable.
 *  - main_effects: Numeric matrix of main effect parameters for all variables and groups.
 *  - main_effect_indices: Integer matrix mapping variables to their main effect parameter indices.
 *  - projection: Numeric matrix specifying group-specific scaling for thresholds.
 *  - sufficient_blume_capel: List of sufficient statistics for the Blume-Capel model,
 *                            including linear and quadratic terms.
 *  - num_persons: Total number of individuals in the dataset.
 *  - num_groups: Total number of groups in the analysis.
 *  - group_indices: Integer matrix specifying start and end indices for each group in `observations`.
 *  - residual_matrix: Numeric matrix storing residual effects for pseudo-likelihood calculations.
 *  - num_categories: Integer matrix specifying the number of categories for each variable and group.
 *
 * Outputs:
 *  - Returns the computed log pseudo-likelihood ratio.
 */
double log_pseudolikelihood_ratio_main_difference_blumecapel(
    const double linear_current,
    const double quadratic_current,
    const double linear_proposed,
    const double quadratic_proposed,
    const int variable,
    const int h,
    const arma::ivec& baseline_category,
    const arma::mat& main_effects,
    const arma::imat& main_effect_indices,
    const arma::mat& projection,
    const List& sufficient_blume_capel,
    const int num_persons,
    const int num_groups,
    const arma::imat& group_indices,
    const arma::mat& residual_matrix,
    const arma::imat& num_categories
) {
  double lbound, bound; // Variables for numerical bounds to ensure stability
  double rest_score, numerator, denominator, exponent; // Variables for pseudo-likelihood computations
  double pseudolikelihood_ratio = 0.0; // Accumulated log pseudo-likelihood ratio
  int linear_score, quadratic_score; // Scores for linear and quadratic terms
  int cat_index = main_effect_indices(variable, 0); // Base index for variable categories

  // Loop over all groups
  for(int gr = 0; gr < num_groups; gr++) {
    // Precomputed constants for numerator and denominator for each category
    arma::vec constant_numerator (num_categories(variable, gr) + 1);
    arma::vec constant_denominator (num_categories(variable, gr) + 1);
    double P = projection(gr, h); // Scaling factor for group-specific adjustments

    // Pre-compute terms for all categories in the current group
    for(int category = 0; category <= num_categories(variable, gr); category++) {
      // Compute linear and quadratic scores for the current category
      linear_score = category;
      quadratic_score = (category - baseline_category[variable]) *
        (category - baseline_category[variable]);

      // Initialize numerator and denominator contributions with main effects
      constant_numerator[category] = main_effects(cat_index, 0) * linear_score;
      constant_numerator[category] +=  main_effects(cat_index + 1, 0) * quadratic_score;
      constant_denominator[category] =  main_effects(cat_index, 0) * linear_score ;
      constant_denominator[category] += main_effects(cat_index + 1, 0) * quadratic_score;

      // Add group-specific contributions from projections
      for(int hh = 1; hh < num_groups; hh++) {
        double P = projection(gr, hh - 1);
        double m1 =  main_effects(cat_index, hh);
        double m2 =  main_effects(cat_index + 1, hh);
        constant_numerator[category] += P * m1 * linear_score;
        constant_numerator[category] += P * m2 * quadratic_score;
        constant_denominator[category] += P * m1 * linear_score;
        constant_denominator[category] += P * m2 * quadratic_score;
      }

      // Adjust denominator with proposed changes for this group
      constant_denominator[category] -= P * main_effects(cat_index, h + 1) * linear_score;
      constant_denominator[category] -= P * main_effects(cat_index + 1, h + 1) * quadratic_score;
      constant_denominator[category] += P * linear_proposed * linear_score;
      constant_denominator[category] += P * quadratic_proposed * quadratic_score;
    }

    // Compute numerical bounds for stability
    double tmp_num = max(constant_numerator);
    double tmp_den = max(constant_denominator);
    if(tmp_num > 0) {
      if(tmp_num > tmp_den) {
        lbound = tmp_num;
      } else {
        lbound = tmp_den;
      }
    } else {
      lbound = 0.0;
    }

    // Add contributions from sufficient statistics
    arma::imat sufficient_blume_capel_gr = sufficient_blume_capel[gr];
    double sp0 = sufficient_blume_capel_gr(0, variable) * P;
    double sp1 = sufficient_blume_capel_gr(1, variable) * P;
    pseudolikelihood_ratio += sp0 * linear_proposed;
    pseudolikelihood_ratio += sp1 * quadratic_proposed;
    pseudolikelihood_ratio -= sp0 * linear_current;
    pseudolikelihood_ratio -= sp1 * quadratic_current;

    // Process each person in the group to compute pseudo-likelihood terms
    for(int person = group_indices(gr, 0); person <= group_indices(gr, 1); person++) {
      rest_score = residual_matrix(person, variable);
      if(rest_score > 0) {
        bound = num_categories(variable, gr) * rest_score + lbound;
      } else {
        bound = lbound;
      }

      // Initialize numerator and denominator
      numerator = 0.0;
      denominator = 0.0;

      // Compute category-specific contributions
      for(int category = 0; category <= num_categories[variable]; category ++) {
        exponent = category * rest_score - bound;
        numerator += std::exp(constant_numerator[category] + exponent);
        denominator += std::exp(constant_denominator[category] + exponent);
      }

      // Update the pseudo-likelihood ratio with log probabilities
      pseudolikelihood_ratio += std::log(numerator);
      pseudolikelihood_ratio -= std::log(denominator);
    }
  }

  return pseudolikelihood_ratio; // Return the computed pseudo-likelihood ratio
}


/**
 * Function: metropolis_main_difference_blumecapel
 * Purpose:
 *   Perform Metropolis-Hastings sampling for the full-conditional distribution
 *   of the group-specific category threshold difference parameters for a Blume-Capel
 *   ordinal variable in an ANOVA model.
 *
 * Inputs:
 *   - main_effects: A matrix of main effects for all variables and groups.
 *   - main_effect_indices: An integer matrix mapping variables to their category indices.
 *   - projection: A projection matrix mapping group-specific effects.
 *   - num_categories: An integer matrix containing the number of categories for
 *                    each variable across groups.
 *   - sufficient_blume_capel: A list of sufficient statistics for Blume-Capel
 *                             variables per group.
 *   - num_persons: Total number of individuals in the dataset.
 *   - num_groups: Total number of groups in the analysis.
 *   - group_indices: A matrix containing the start and end indices for individuals
 *                  in each group.
 *   - variable: Index of the variable being updated.
 *   - baseline_category: A vector indicating the reference category for each variable.
 *   - prior_threshold_alpha: Shape parameter of the prior distribution for thresholds.
 *   - prior_threshold_beta: Rate parameter of the prior distribution for thresholds.
 *   - residual_matrix: A matrix of residual scores for the pseudo-likelihood calculation.
 *   - inclusion_indicator: A binary matrix indicating whether a variable is active for sampling.
 *   - proposal_sd_main_effects: A matrix of proposal standard deviations for the main effects.
 *   - rm_adaptation_rate: Robbins-Monro learning rate parameter for adaptive Metropolis-Hastings.
 *   - target_acceptance_rate: Target acceptance probability for Robbins-Monro.
 *   - t: Current iteration number in the sampling procedure.
 *   - rm_lower_bound: Lower bound for the proposal standard deviation.
 *   - rm_upper_bound: Upper bound for the proposal standard deviation.
 *
 * Outputs:
 *   - Updates `main_effects` to reflect accepted proposals for the linear and
 *     quadratic threshold differences.
 *   - Updates `proposal_sd_main_effects` for the linear and quadratic parameters using
 *     Robbins-Monro updates.
 */
void metropolis_main_difference_blumecapel(
    arma::mat& main_effects,
    const arma::imat& main_effect_indices,
    const arma::mat& projection,
    const arma::imat& num_categories,
    const List& sufficient_blume_capel,
    const int num_persons,
    const int num_groups,
    const arma::imat& group_indices,
    const int variable,
    const arma::ivec& baseline_category,
    const double prior_threshold_alpha,
    const double prior_threshold_beta,
    arma::mat& residual_matrix,
    const arma::imat& inclusion_indicator,
    arma::mat& proposal_sd_main_effects,
    const double rm_adaptation_rate,
    const double target_acceptance_rate,
    const int t,
    const double rm_lower_bound,
    const double rm_upper_bound
) {
  double log_acceptance_probability, U; // Log probability and random uniform value for MH step
  double current_state, proposed_state; // Current and proposed parameter values
  double exp_neg_log_t_rm_adaptation_rate = std::exp(-std::log(t) * rm_adaptation_rate); // Precomputed Robbins-Monro term

  // Check if the variable is active for sampling
  if(inclusion_indicator(variable, variable) == 1) {

    // Loop over the group-specific difference effects
    for(int h = 1; h < num_groups; h++) {

      // Adaptive Metropolis procedure for the linear Blume-Capel parameter
      int cat_index = main_effect_indices(variable, 0); // Base index for the variable's categories
      current_state = main_effects(cat_index, h); // Current value for linear parameter
      proposed_state = R::rnorm(current_state, proposal_sd_main_effects(cat_index, h)); // Propose new value

      // Compute log pseudo-likelihood ratio for the proposed and current states
      double linear_current = current_state;
      double quadratic_current = main_effects(cat_index + 1, h);
      double linear_proposed = proposed_state;
      double quadratic_proposed = main_effects(cat_index + 1, h);

      log_acceptance_probability = log_pseudolikelihood_ratio_main_difference_blumecapel(
        linear_current, quadratic_current, linear_proposed, quadratic_proposed,
        variable, h - 1, baseline_category, main_effects, main_effect_indices,
        projection, sufficient_blume_capel, num_persons, num_groups, group_indices,
        residual_matrix, num_categories);

      // Add prior contributions to the log probability
      log_acceptance_probability += prior_threshold_alpha * (proposed_state - current_state);
      log_acceptance_probability += (prior_threshold_alpha + prior_threshold_beta) * std::log(1 + std::exp(current_state));
      log_acceptance_probability -= (prior_threshold_alpha + prior_threshold_beta) * std::log(1 + std::exp(proposed_state));

      // Metropolis-Hastings acceptance step
      U = R::unif_rand();
      if(std::log(U) < log_acceptance_probability) {
        main_effects(cat_index, h) = proposed_state; // Accept proposed value
      }

      // Robbins-Monro update for the proposal standard deviation
      proposal_sd_main_effects(cat_index, h) = update_with_robbins_monro(
        proposal_sd_main_effects(cat_index, h), log_acceptance_probability, target_acceptance_rate, rm_lower_bound,
        rm_upper_bound, exp_neg_log_t_rm_adaptation_rate);


      // Adaptive Metropolis procedure for the quadratic Blume-Capel parameter
      current_state = main_effects(cat_index + 1, h); // Current value for quadratic parameter
      proposed_state = R::rnorm(current_state, proposal_sd_main_effects(cat_index + 1, h)); // Propose new value

      // Compute log pseudo-likelihood ratio for the proposed and current states
      linear_current = main_effects(cat_index, h);
      quadratic_current = current_state;
      linear_proposed =  main_effects(cat_index, h);
      quadratic_proposed =  proposed_state;

      log_acceptance_probability = log_pseudolikelihood_ratio_main_difference_blumecapel(
        linear_current, quadratic_current, linear_proposed, quadratic_proposed,
        variable, h - 1, baseline_category, main_effects, main_effect_indices,
        projection, sufficient_blume_capel, num_persons, num_groups, group_indices,
        residual_matrix, num_categories);

      // Add prior contributions to the log probability
      log_acceptance_probability += prior_threshold_alpha * (proposed_state - current_state);
      log_acceptance_probability += (prior_threshold_alpha + prior_threshold_beta) * std::log(1 + std::exp(current_state));
      log_acceptance_probability -= (prior_threshold_alpha + prior_threshold_beta) * std::log(1 + std::exp(proposed_state));

      // Metropolis-Hastings acceptance step
      U = R::unif_rand();
      if(std::log(U) < log_acceptance_probability) {
        main_effects(cat_index + 1, h) = proposed_state;
      }

      // Robbins-Monro update for the proposal standard deviation
      proposal_sd_main_effects(cat_index + 1, h) = update_with_robbins_monro(
        proposal_sd_main_effects(cat_index + 1, h), log_acceptance_probability, target_acceptance_rate, rm_lower_bound,
        rm_upper_bound, exp_neg_log_t_rm_adaptation_rate);
    }
  }
}


/**
 * Function: metropolis_thresholds_regular_free
 * Purpose:
 *   Performs Metropolis-Hastings sampling for the full-conditional distribution
 *   of threshold parameters for a regular binary or ordinal variable in a
 *   free-threshold independent samples design.
 *
 * Inputs:
 *   - main_effects: A matrix of main effects for all variables and groups.
 *   - main_effect_indices: An integer matrix mapping variables to their category indices.
 *   - observations: An integer matrix containing observed data for individuals
 *                   by variables.
 *   - num_groups: Total number of groups in the analysis.
 *   - group_indices: A matrix containing start and end indices for individuals
 *                  in each group.
 *   - num_categories: A matrix containing the number of categories for each
 *                    variable and group.
 *   - residual_matrix: A matrix of residual scores for the pseudo-likelihood
 *                  calculations.
 *   - n_cat_obs: A list of matrices containing category-specific counts for
 *                each group.
 *   - prior_threshold_alpha: Shape parameter for the prior distribution.
 *   - prior_threshold_beta: Rate parameter for the prior distribution.
 *   - variable: Index of the variable being updated.
 *   - group: Index of the group being updated.
 *
 * Outputs:
 *   - Updates the `main_effects` matrix to reflect accepted proposals for
 *     threshold parameters.
 */
void metropolis_thresholds_regular_free(
    arma::mat& main_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& observations,
    const int num_groups,
    const arma::imat& group_indices,
    const arma::imat& num_categories,
    arma::mat& residual_matrix,
    const List& n_cat_obs,
    const double prior_threshold_alpha,
    const double prior_threshold_beta,
    const int variable,
    const int group
) {

  // Number of persons in the group
  int num_persons = group_indices(group, 1) - group_indices(group, 0) + 1;

  // Cache group-specific category observations
  arma::imat n_cat_obs_gr = n_cat_obs[group];

  // Base category index for the variable
  int base_cat_index = main_effect_indices(variable, 0);
  int n_cats = num_categories(variable, group);

  // Temporary storage for pseudo-likelihood elements
  arma::vec q(num_persons);
  arma::vec r(num_persons);

  // Loop over categories
  for(int category = 0; category < n_cats; category++) {
    double current_state = main_effects(base_cat_index + category, group);
    double exp_current = std::exp(current_state);

    // Initialize scaling factor `c` for the generalized beta-prime proposal
    double c = (prior_threshold_alpha + prior_threshold_beta) / (1 + exp_current);

    // Compute `q` and `r` for each person in the group
    for(int person = 0; person < num_persons; person++) {
      double rest_score = residual_matrix(group_indices(group, 0) + person, variable);

      // Calculate pseudo-likelihood component `q`
      double q_person = 1.0;
      for(int cat = 0; cat < n_cats; cat++) {
        if(cat != category) {
          q_person += std::exp(main_effects(base_cat_index + cat, group) + (cat + 1) * rest_score);
        }
      }

      // Calculate pseudo-likelihood component `r`
      double r_person = std::exp((category + 1) * rest_score);

      // Store results for each person
      q[person] = q_person;
      r[person] = r_person;

      // Update scaling factor `c`
      c += r_person / (q_person + r_person * exp_current);
    }

    // Finalize `c` by accounting for prior contributions
    c /= (num_persons + prior_threshold_alpha + prior_threshold_beta - exp_current * c);

    // Propose a new state using the generalized beta-prime distribution
    double a = n_cat_obs_gr(category + 1, variable) + prior_threshold_alpha;
    double b = num_persons + prior_threshold_beta - n_cat_obs_gr(category + 1, variable);
    double tmp = R::rbeta(a, b);
    double proposed_state = std::log(tmp / ((1 - tmp) * c));
    double exp_proposed = std::exp(proposed_state);


    // Compute the log acceptance probability
    double log_acceptance_probability = 0.0;

    // Add pseudo-likelihood ratio contributions
    for (int person = 0; person < num_persons; ++person) {
      log_acceptance_probability += std::log(q[person] + r[person] * exp_current);
      log_acceptance_probability -= std::log(q[person] + r[person] * exp_proposed);
    }

    // Add prior ratio contributions
    log_acceptance_probability -= (prior_threshold_alpha + prior_threshold_beta) * std::log(1 + exp_proposed);
    log_acceptance_probability += (prior_threshold_alpha + prior_threshold_beta) * std::log(1 + exp_current);

    // Add proposal ratio contributions
    log_acceptance_probability -= (a + b) * std::log(1 + c * exp_current);
    log_acceptance_probability += (a + b) * std::log(1 + c * exp_proposed);

    // Perform the Metropolis-Hastings acceptance step
    double U = std::log(R::unif_rand());
    if(U < log_acceptance_probability) {
      // Update the main effects matrix if the proposal is accepted
      main_effects(base_cat_index + category, group) = proposed_state;
    }
  }
}

/**
 * Function: metropolis_thresholds_blumecapel_free
 * Purpose:
 *   Performs Metropolis-Hastings sampling for the full-conditional distribution
 *   of the threshold parameters (linear and quadratic) for a Blume-Capel ordinal
 *   variable in a free-threshold independent samples design.
 *
 * Inputs:
 *   - main_effects: A matrix of main effects for all variables and groups.
 *   - main_effect_indices: An integer matrix mapping variables to their category indices.
 *   - observations: An integer matrix containing observed data for individuals
 *                   by variables.
 *   - num_groups: Total number of groups in the analysis.
 *   - group_indices: A matrix containing start and end indices for individuals
 *                  in each group.
 *   - baseline_category: A vector specifying the reference category for each variable.
 *   - num_categories: A matrix containing the number of categories for each
 *                    variable and group.
 *   - sufficient_blume_capel: A list of matrices containing sufficient statistics
 *                             for each group.
 *   - num_persons: Total number of individuals in the dataset.
 *   - residual_matrix: A matrix of residual scores for pseudo-likelihood calculations.
 *   - n_cat_obs: A list of matrices containing category-specific counts for
 *                each group and variable.
 *   - prior_threshold_alpha: Shape parameter for the prior distribution.
 *   - prior_threshold_beta: Rate parameter for the prior distribution.
 *   - variable: Index of the variable being updated.
 *   - group: Index of the group being updated.
 *   - proposal_sd_main_effects: A matrix of proposal standard deviations for the
 *                       Metropolis-Hastings updates.
 *   - rm_adaptation_rate: Robbins-Monro learning rate parameter.
 *   - target_acceptance_rate: Target acceptance rate for
 *                                         Metropolis-Hastings updates.
 *   - t: Current iteration number.
 *   - rm_lower_bound: Lower bound for Robbins-Monro updates to proposal standard deviations.
 *   - rm_upper_bound: Upper bound for Robbins-Monro updates to proposal standard deviations.
 *
 * Outputs:
 *   - Updates the `main_effects` matrix to reflect accepted proposals for
 *     the linear and quadratic threshold parameters.
 *   - Adjusts the `proposal_sd_main_effects` matrix using Robbins-Monro adaptive updates.
 */
void metropolis_thresholds_blumecapel_free(
    arma::mat& main_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& observations,
    const int num_groups,
    const arma::imat& group_indices,
    const arma::ivec& baseline_category,
    const arma::imat& num_categories,
    const List& sufficient_blume_capel,
    const int num_persons,
    arma::mat& residual_matrix,
    const List& n_cat_obs,
    const double prior_threshold_alpha,
    const double prior_threshold_beta,
    const int variable,
    const int group,
    arma::mat& proposal_sd_main_effects,
    const double rm_adaptation_rate,
    const double target_acceptance_rate,
    const int t,
    const double rm_lower_bound,
    const double rm_upper_bound
) {

  // Robbins-Monro term precomputed for efficiency
  double exp_neg_log_t_rm_adaptation_rate = std::exp(-std::log(t) * rm_adaptation_rate);

  // Retrieve sufficient statistics for the current group
  arma::mat sufficient_blume_capel_group = sufficient_blume_capel[group];

  // Preallocate vectors for constants used in likelihood calculations
  arma::vec constant_numerator(num_categories(variable, group) + 1);
  arma::vec constant_denominator(num_categories(variable, group) + 1);

  // Adaptive Metropolis procedure for the linear Blume-Capel parameter
  int cat_index = main_effect_indices(variable, 0);
  double current_state = main_effects(cat_index, group);
  double proposed_state = R::rnorm(current_state, proposal_sd_main_effects(cat_index, group));

  // Difference between proposed and current state
  double difference = proposed_state - current_state;

  // Precompute constants for likelihood ratio calculation
  for(int category = 0; category <= num_categories(variable, group); category ++) {
    double exponent = main_effects(cat_index + 1, group) *
      (category - baseline_category[variable]) *
      (category - baseline_category[variable]);
    constant_numerator[category] = current_state * category + exponent;
    constant_denominator[category] = proposed_state * category + exponent;
  }

  // Precompute bounds for numerical stability
  double tmp_n = max(constant_numerator);
  double tmp_d = max(constant_denominator);
  double lbound = 0.0;
  if(tmp_n > 0) {
    if(tmp_n > tmp_d) {
      lbound = tmp_n;
    } else {
      lbound = tmp_d;
    }
  }

  // Initialize log acceptance probability with prior contribution
  double log_acceptance_probability = prior_threshold_alpha * difference;
  log_acceptance_probability += sufficient_blume_capel_group(0, variable) * difference;

  // Loop over individuals in the group to compute likelihood ratio
  for(int person = group_indices(group, 0); person <= group_indices(group, 1); person++) {
    double rest_score = residual_matrix(person, variable);
    double bound = lbound;
    if(rest_score > 0) {
      bound += num_categories(variable, group) * rest_score;
    }

    // Compute likelihood numerator and denominator
    double numerator = std::exp(constant_numerator[0] - bound);
    double denominator = std::exp(constant_denominator[0] - bound);
    for(int score = 1; score <= num_categories(variable, group); score++) {
      double exponent = score * rest_score - bound;
      numerator += std::exp(constant_numerator[score] + exponent);
      denominator += std::exp(constant_denominator[score] + exponent);
    }
    log_acceptance_probability += std::log(numerator);
    log_acceptance_probability -= std::log(denominator);
  }

  // Add prior ratio to log acceptance probability
  log_acceptance_probability += (prior_threshold_alpha + prior_threshold_beta) * std::log(1 + std::exp(current_state));
  log_acceptance_probability -= (prior_threshold_alpha + prior_threshold_beta) * std::log(1 + std::exp(proposed_state));

  // Metropolis acceptance step
  double U = R::unif_rand();
  if(std::log(U) < log_acceptance_probability) {
    main_effects(cat_index, group) = proposed_state;
  }

  // Robbins-Monro adaptive update for proposal standard deviation
  proposal_sd_main_effects(cat_index, group) = update_with_robbins_monro(
    proposal_sd_main_effects(cat_index, group), log_acceptance_probability, target_acceptance_rate, rm_lower_bound,
    rm_upper_bound, exp_neg_log_t_rm_adaptation_rate);


  // Adaptive Metropolis procedure for the quadratic Blume-Capel parameter
  current_state = main_effects(cat_index + 1, group);
  proposed_state = R::rnorm(current_state, proposal_sd_main_effects(cat_index + 1, group));
  difference = proposed_state - current_state;

  // Recompute constants for quadratic term
  for(int category = 0; category <= num_categories(variable, group); category ++) {
    double exponent = main_effects(cat_index, group) * category;
    int score = (category - baseline_category[variable]) *
      (category - baseline_category[variable]);

    constant_numerator[category] = current_state * score + exponent;
    constant_denominator[category] = proposed_state * score + exponent;
  }

  tmp_n = max(constant_numerator);
  tmp_d = max(constant_denominator);
  if(tmp_n > 0) {
    if(tmp_n > tmp_d) {
      lbound = tmp_n;
    } else {
      lbound = tmp_d;
    }
  } else {
    lbound = 0.0;
  }

  log_acceptance_probability = prior_threshold_alpha * difference;
  log_acceptance_probability += sufficient_blume_capel_group(1, variable) * difference;

  // Loop over individuals to compute likelihood ratio for quadratic term
  for(int person = group_indices(group, 0); person <= group_indices(group, 1); person++) {
    double rest_score = residual_matrix(person, variable);
    double bound = lbound;
    if(rest_score > 0) {
      bound += num_categories[variable] * rest_score;
    }

    double numerator = std::exp(constant_numerator[0] - bound);
    double denominator = std::exp(constant_denominator[0] - bound);

    for(int score = 1; score <= num_categories(variable, group); score ++) {
      double exponent = score * rest_score - bound;
      numerator += std::exp(constant_numerator[score] + exponent);
      denominator += std::exp(constant_denominator[score] + exponent);
    }

    log_acceptance_probability += std::log(numerator);
    log_acceptance_probability -= std::log(denominator);
  }
  log_acceptance_probability += (prior_threshold_alpha + prior_threshold_beta) * std::log(1 + std::exp(current_state));
  log_acceptance_probability -= (prior_threshold_alpha + prior_threshold_beta) * std::log(1 + std::exp(proposed_state));

  U = R::unif_rand();
  if(std::log(U) < log_acceptance_probability) {
    main_effects(cat_index + 1, group) = proposed_state;
  }

  // Robbins-Monro update for the proposal standard deviation
  proposal_sd_main_effects(cat_index + 1, group) = update_with_robbins_monro(
    proposal_sd_main_effects(cat_index + 1, group), log_acceptance_probability, target_acceptance_rate, rm_lower_bound,
    rm_upper_bound, exp_neg_log_t_rm_adaptation_rate);

}


/**
 * Computes the log pseudo-likelihood ratio for a regular ordinal variable
 * when comparing two models (proposed vs current) for group-level main
 * difference parameters.
 *
 * @param current_main_effects Matrix of current main effects for all categories and groups.
 * @param proposed_main_effects Matrix of proposed main effects for all categories and groups.
 * @param projection Matrix of group-specific projection weights.
 * @param num_groups Number of groups in the data.
 * @param group_indices Matrix of group start and end indices for individuals.
 * @param num_categories Matrix of the number of categories per variable and group.
 * @param num_persons Total number of individuals.
 * @param residual_matrix Matrix of residual scores for individuals and variables.
 * @param n_cat_obs List of category observation counts for each group.
 * @param variable Index of the variable being evaluated.
 * @return The log pseudo-likelihood ratio comparing the proposed and current models.
 */
double log_pseudolikelihood_ratio_main_difference_regular_between_model(
    const arma::mat& current_main_effects,
    const arma::mat& proposed_main_effects,
    const arma::mat& projection,
    const int num_groups,
    const arma::imat& group_indices,
    const arma::imat& num_categories,
    const arma::mat& residual_matrix,
    const List& n_cat_obs,
    const int variable
) {
  double pseudolikelihood_ratio = 0.0; // Initialize the log pseudo-likelihood ratio
  int num_cats = num_categories(variable, 0); // Number of categories for the variable
  arma::vec current_thresholds(num_cats); // Store current thresholds for the group
  arma::vec proposed_thresholds(num_cats); // Store proposed thresholds for the group

  // Loop over all groups to compute the contribution to the log pseudo-likelihood ratio
  for(int gr = 0; gr < num_groups; gr++) {
    // Retrieve group-specific category observation counts
    arma::imat n_cat_obs_gr = n_cat_obs[gr];

    // Compute category thresholds for each category within the group
    for(int cat = 0; cat < num_cats; cat++) {
      // Compute current and proposed category thresholds
      double current_threshold = current_main_effects(cat, 0);
      double proposed_threshold = proposed_main_effects(cat, 0);
      for (int h = 1; h < num_groups; ++h) {
        current_threshold += projection(gr, h - 1) * current_main_effects(cat, h);
        proposed_threshold += projection(gr, h - 1) * proposed_main_effects(cat, h);
      }
      current_thresholds[cat] = current_threshold;
      proposed_thresholds[cat] = proposed_threshold;

      // Update the pseudo-likelihood ratio with contributions from the observations
      double delta_state = proposed_threshold - current_threshold;
      pseudolikelihood_ratio += delta_state * n_cat_obs_gr(cat + 1, variable);
    }

    // Loop over individuals within the group
    for(int person = group_indices(gr, 0); person <= group_indices(gr, 1); person++) {
      // Compute the rest score (residual) for the current variable
      double rest_score = residual_matrix(person, variable);

      // Compute numerical bounds for stability
      double bound = (rest_score > 0) ? num_cats * rest_score : 0.0;

      // Initialize denominator terms for the proposed and current pseudolikelihoods
      double denominator_proposed = std::exp(-bound);
      double denominator_current = std::exp(-bound);

      // Add contributions from each category to the denominators
      for(int cat = 0; cat < num_cats; cat++) {
        double exponent = (cat + 1) * rest_score - bound;
        denominator_proposed += std::exp(exponent + proposed_thresholds[cat]);
        denominator_current += std::exp(exponent + current_thresholds[cat]);
      }

      // Update the pseudo-likelihood ratio with log-likelihood differences
      pseudolikelihood_ratio -= std::log(denominator_proposed);
      pseudolikelihood_ratio += std::log(denominator_current);
    }
  }

  return pseudolikelihood_ratio; // Return the computed log pseudo-likelihood ratio
}


/**
 * Perform Metropolis-Hastings sampling for group-level main difference parameters
 * for regular ordinal variables between models (current vs. proposed).
 *
 * @param inclusion_indicator Matrix indicating inclusion for variables.
 * @param inclusion_probability_difference Matrix of inclusion probabilities for variable differences.
 * @param main_effects Matrix of current main effects for all variables and groups.
 * @param main_effect_indices Matrix of indices for main effect parameters.
 * @param observations Matrix of observations for all individuals and variables.
 * @param num_groups Number of groups in the model.
 * @param group_indices Matrix of start and end indices for individuals in each group.
 * @param num_categories Matrix indicating the number of categories for each variable and group.
 * @param num_persons Total number of individuals in the dataset.
 * @param residual_matrix Matrix of residual scores for individuals and variables.
 * @param n_cat_obs List of category observation counts for each group.
 * @param prior_threshold_alpha Alpha parameter for prior threshold.
 * @param prior_threshold_beta Beta parameter for prior threshold.
 * @param variable Index of the variable being updated.
 * @param group Index of the group being updated.
 * @param proposal_sd_main_effects Matrix of proposal standard deviations for main effects.
 * @param main_difference_scale Scale parameter for the Cauchy prior on main differences.
 * @param projection Matrix of projection weights for groups.
 */
void metropolis_main_difference_regular_between_model(
    arma::imat& inclusion_indicator,
    const arma::mat& inclusion_probability_difference,
    arma::mat& main_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& observations,
    const int num_groups,
    const arma::imat& group_indices,
    const arma::imat& num_categories,
    const arma::mat& residual_matrix,
    const List& n_cat_obs,
    const int variable,
    const arma::mat& proposal_sd_main_effects,
    const double main_difference_scale,
    const arma::mat& projection
) {
  // Number of rows for this variable's main effects
  int num_row = main_effect_indices(variable, 1) -  main_effect_indices(variable, 0) + 1;
  arma::mat proposed_main_effects(num_row, num_groups);
  arma::mat current_main_effects(num_row, num_groups);

  // Retrieve the current and proposed indicators for inclusion
  int current_indicator = inclusion_indicator(variable, variable);
  int proposed_indicator = 1 - current_indicator;

  int main_effect_index = main_effect_indices(variable, 0);
  int num_cats = num_categories(variable, 0);
  double log_prob = 0.0;

  // Populate the current and proposed main effects matrices
  for(int cat = 0; cat < num_cats; cat++) {
    proposed_main_effects(cat, 0) = main_effects(main_effect_index + cat, 0);
    current_main_effects(cat, 0) = main_effects(main_effect_index + cat, 0);

    if(current_indicator == 1) {
      for(int h = 1; h < num_groups; h++) {
        double proposal_sd = proposal_sd_main_effects(main_effect_index + cat, h);
        double current_state = main_effects(main_effect_index + cat, h);
        double proposed_state = 0.0;
        current_main_effects(cat, h) = current_state;
        proposed_main_effects(cat, h) = proposed_state;

        // Update log acceptance probability with prior and proposal contributions
        log_prob -= R::dcauchy(current_state, 0.0, main_difference_scale, true);
        log_prob += R::dnorm(current_state, proposed_state, proposal_sd, true);
      }
    } else {
      for(int h = 1; h < num_groups; h++) {
        double proposal_sd = proposal_sd_main_effects(main_effect_index + cat, h);
        double current_state = 0.0;
        double proposed_state = R::rnorm(current_state, proposal_sd);
        current_main_effects(cat, h) = current_state;
        proposed_main_effects(cat, h) = proposed_state;

        // Update log acceptance probability with prior and proposal contributions
        log_prob += R::dcauchy(proposed_state, 0.0, main_difference_scale, true);
        log_prob -= R::dnorm(proposed_state, current_state, proposal_sd, true);
      }
    }
  }

  // Add pseudo-likelihood ratio contribution
  log_prob += log_pseudolikelihood_ratio_main_difference_regular_between_model(
    current_main_effects, proposed_main_effects, projection, num_groups,
    group_indices, num_categories, residual_matrix, n_cat_obs, variable
  );

  // Add prior inclusion odds contributions
  if(current_indicator == 1) {
    log_prob -= std::log(inclusion_probability_difference(variable, variable));
    log_prob += std::log(1 - inclusion_probability_difference(variable, variable));
  } else {
    log_prob += std::log(inclusion_probability_difference(variable, variable));
    log_prob -= std::log(1 - inclusion_probability_difference(variable, variable));
  }

  // Perform Metropolis-Hastings step
  double U = R::unif_rand();
  if(std::log(U) < log_prob) {
    inclusion_indicator(variable, variable) = proposed_indicator;
    for(int cat = 0; cat < num_cats; cat++) {
      for(int h = 1; h < num_groups; h++) {
        main_effects(main_effect_index + cat, h) = proposed_main_effects(cat, h);
      }
    }
  }
}


double log_pseudolikelihood_ratio_main_difference_blume_capel_between_model(
    const arma::mat& current_main_effects,
    const arma::mat& proposed_main_effects,
    const arma::mat& projection,
    const arma::ivec& baseline_category,
    const List& sufficient_blume_capel,
    const int num_groups,
    const arma::imat& group_indices,
    const arma::imat& num_categories,
    const int num_persons,
    const arma::mat& residual_matrix,
    const int variable
) {
  double pseudolikelihood_ratio = 0.0; // Initialize the log pseudo-likelihood ratio
  int num_cats = num_categories(variable, 0); // Number of categories for the variable
  arma::vec current_parameters(2); // Store current parameters for the group
  arma::vec proposed_parameters(2); // Store proposed parameters for the group

  // Loop over all groups
  for(int gr = 0; gr < num_groups; gr++) {
    arma::imat sufficient_statistics = sufficient_blume_capel[gr]; // Group-specific data

    // Compute current and proposed values for the linear Blume-Capel parameter
    double current_linear = current_main_effects(0, 0);
    double proposed_linear = proposed_main_effects(0, 0);
    for (int h = 1; h < num_groups; ++h) {
      current_linear += projection(gr, h - 1) * current_main_effects(0, h);
      proposed_linear += projection(gr, h - 1) * proposed_main_effects(0, h);
    }
    current_parameters[0] = current_linear;
    proposed_parameters[0] = proposed_linear;
    // Add the contribution from delta_state based on observations
    double delta_state = proposed_linear - current_linear;
    pseudolikelihood_ratio += delta_state * sufficient_statistics(0, variable);

    // Compute current and proposed values for the quadratic Blume-Capel parameter
    double current_quadratic = current_main_effects(0, 0);
    double proposed_quadratic = proposed_main_effects(0, 0);
    for (int h = 1; h < num_groups; ++h) {
      current_quadratic += projection(gr, h - 1) * current_main_effects(1, h);
      proposed_quadratic += projection(gr, h - 1) * proposed_main_effects(1, h);
    }
    current_parameters[1] = current_quadratic;
    proposed_parameters[1] = proposed_quadratic;
    // Add the contribution from delta_state based on observations
    delta_state = proposed_quadratic - current_quadratic;
    pseudolikelihood_ratio += delta_state * sufficient_statistics(1, variable);

    // Precomputed constants for numerator and denominator for each category
    arma::vec current_denominator_constant (num_cats + 1);
    arma::vec proposed_denominator_constant (num_cats + 1);

    // Pre-compute terms for all categories in the current group
    for(int cat = 0; cat <= num_cats; cat++) {
      // Compute linear and quadratic scores for the current cat
      int linear_score = cat;
      int quadratic_score = (cat - baseline_category[variable]) *
                            (cat - baseline_category[variable]);

      // Initialize numerator and denominator contributions with main effects
      current_denominator_constant[cat] = current_parameters[0] * linear_score;
      current_denominator_constant[cat] += current_parameters[1] * quadratic_score;
      proposed_denominator_constant[cat] = proposed_parameters[0] * linear_score ;
      proposed_denominator_constant[cat] += proposed_parameters[1] * quadratic_score;
    }

    // Compute numerical bounds for stability
    double tmp_den1 = max(current_denominator_constant);
    double tmp_max = max(proposed_denominator_constant);
    if(tmp_den1 > tmp_max) {
      tmp_max = tmp_den1;
    }
    double lbound = 0.0;
    if(tmp_max > 0.0) {
      lbound = tmp_max;
    }

    // Loop over all persons in the group
    for(int person = group_indices(gr, 0); person <= group_indices(gr, 1); person++) {
      double rest_score = residual_matrix(person, variable); // Compute residual score
      double bound = (rest_score > 0) ? lbound + num_cats * rest_score : lbound;

      // Compute the denominators for proposed and current thresholds
      double denominator_proposed = 0.0;
      double denominator_current = 0.0;

      // Compute category-specific contributions
      for(int cat = 0; cat <= num_categories[variable]; cat++) {
        double exponent = cat * rest_score - bound;
        denominator_proposed += std::exp(proposed_denominator_constant[cat] + exponent);
        denominator_current += std::exp(current_denominator_constant[cat] + exponent);
      }

      // Update the pseudo-likelihood ratio with log-likelihood differences
      pseudolikelihood_ratio -= std::log(denominator_proposed);
      pseudolikelihood_ratio += std::log(denominator_current);
    }
  }

  return pseudolikelihood_ratio;
}

void metropolis_main_difference_blume_capel_between_model(
    arma::imat& inclusion_indicator,
    const arma::mat& inclusion_probability_difference,
    arma::mat& main_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& observations,
    const int num_groups,
    const arma::imat& group_indices,
    const arma::ivec& baseline_category,
    const arma::imat& num_categories,
    const List& sufficient_blume_capel,
    const int num_persons,
    const arma::mat& residual_matrix,
    const List& n_cat_obs,
    const double prior_threshold_alpha,
    const double prior_threshold_beta,
    const int variable,
    const arma::mat& proposal_sd_main_effects,
    const double rm_adaptation_rate,
    const double target_acceptance_rate,
    const int t,
    const double rm_lower_bound,
    const double rm_upper_bound,
    const arma::uvec& is_ordinal_variable,
    const double main_difference_scale,
    const arma::mat projection
) {
  arma::mat proposed_main_effects(2, num_groups);
  arma::mat current_main_effects(2, num_groups);
  int current_inclusion_indicator = inclusion_indicator(variable, variable);
  int proposed_inclusion_indicator = 1 - current_inclusion_indicator;
  int main_effect_index = main_effect_indices(variable, 0);
  double log_prob = 0.0;

  proposed_main_effects(0, 0) = main_effects(main_effect_index, 0);
  current_main_effects(0, 0) = main_effects(main_effect_index, 0);
  proposed_main_effects(1, 0) = main_effects(main_effect_index + 1, 0);
  current_main_effects(1, 0) = main_effects(main_effect_index + 1, 0);

  if(inclusion_indicator(variable, variable) == 1) {
    for(int h = 1; h < num_groups; h++) {
      // First, the linear parameter
      double proposal_sd = proposal_sd_main_effects(main_effect_index, h);
      double current_state = main_effects(main_effect_index, h);
      double proposed_state = 0.0;
      current_main_effects(0, h) = current_state;
      proposed_main_effects(0, h) = proposed_state;
      //Prior and proposal ratio
      log_prob -= R::dcauchy(proposed_state, 0.0, main_difference_scale, true);
      log_prob += R::dnorm(proposed_state, current_state, proposal_sd, true);

      // Second, the quadratic parameter
      proposal_sd = proposal_sd_main_effects(main_effect_index + 1, h);
      current_state = 0.0;
      proposed_state = R::rnorm(current_state, proposal_sd);
      current_main_effects(1, h) = current_state;
      proposed_main_effects(1, h) = proposed_state;
      //Prior and proposal ratio
      log_prob += R::dcauchy(proposed_state, 0.0, main_difference_scale, true);
      log_prob -= R::dnorm(proposed_state, current_state, proposal_sd, true);
    }
  } else {
    for(int h = 1; h < num_groups; h++) {
      // First, the linear parameter
      double proposal_sd = proposal_sd_main_effects(main_effect_index, h);
      double current_state = 0.0;
      double proposed_state = R::rnorm(current_state, proposal_sd);
      current_main_effects(0, h) = current_state;
      proposed_main_effects(0, h) = proposed_state;
      //Prior and proposal ratio
      log_prob += R::dcauchy(proposed_state, 0.0, main_difference_scale, true);
      log_prob -= R::dnorm(proposed_state, current_state, proposal_sd, true);

      // Second, the quadratic parameter
      proposal_sd = proposal_sd_main_effects(main_effect_index + 1, h);
      current_state = 0.0;
      proposed_state = R::rnorm(current_state, proposal_sd);
      current_main_effects(1, h) = current_state;
      proposed_main_effects(1, h) = proposed_state;
      //Prior and proposal ratio
      log_prob += R::dcauchy(proposed_state, 0.0, main_difference_scale, true);
      log_prob -= R::dnorm(proposed_state, current_state, proposal_sd, true);
    }
  }

  // Pseudolikelihood ratio
  log_prob += log_pseudolikelihood_ratio_main_difference_blume_capel_between_model(
    current_main_effects, proposed_main_effects, projection,
    baseline_category, sufficient_blume_capel, num_groups, group_indices,
    num_categories, num_persons, residual_matrix, variable);


  // Prior odds
  if(current_inclusion_indicator == 1) {
    log_prob -= std::log(inclusion_probability_difference(variable, variable));
    log_prob += std::log(1-inclusion_probability_difference(variable, variable));
  } else {
    log_prob += std::log(inclusion_probability_difference(variable, variable));
    log_prob -= std::log(1-inclusion_probability_difference(variable, variable));
  }

  double U = R::unif_rand();
  if(std::log(U) < log_prob) {
    inclusion_indicator(variable, variable) = proposed_inclusion_indicator;
    for(int h = 1; h < num_groups; h++) {
      main_effects(main_effect_index, h) = proposed_main_effects(0, h);
      main_effects(main_effect_index + 1, h) = proposed_main_effects(1, h);
    }
  }
}

/**
 * Function: gibbs_step_gm
 * Purpose:
 *   Executes a single Gibbs sampling step to update grarm_adaptation_ratecal model parameters
 *   for a Bayesian parameter comparison across multiple groups.
 *
 * Inputs:
 *   - main_effects: arma::mat containing the main effects across variables and groups.
 *   - main_effect_indices: arma::imat mapping variables to their corresponding category indices.
 *   - pairwise_effects: arma::mat for pairwise effects between variables.
 *   - pairwise_effect_indices: arma::imat mapping variable pairs to their indices in `pairwise_effects`.
 *   - projection: arma::mat representing the group-specific scaling projection.
 *   - num_categories: arma::imat containing the number of categories for each variable and group.
 *   - observations: arma::imat with observed data for individuals by variables.
 *   - num_persons: Total number of individuals in the dataset.
 *   - num_groups: Total number of groups.
 *   - group_indices: arma::imat specifying start and end indices for individuals in each group.
 *   - n_cat_obs: List of category-specific counts for each group and variable.
 *   - sufficient_blume_capel: List of sufficient statistics for Blume-Capel variables in each group.
 *   - residual_matrix: arma::mat of residual scores used in pseudo-likelihood calculations.
 *   - independent_thresholds: Boolean indicating if thresholds are modeled independently for groups.
 *   - is_ordinal_variable: arma::uvec indicating if variables are ordinal.
 *   - baseline_category: arma::ivec specifying reference categories for each variable.
 *   - inclusion_indicator: arma::imat indicating active pairwise and main differences.
 *   - inclusion_probability_difference: arma::mat for inclusion probabilities of pairwise differences.
 *   - proposal_sd_main_effects: arma::mat of proposal standard deviations for main effects.
 *   - proposal_sd_pairwise_effects: arma::mat of proposal standard deviations for pairwise effects.
 *   - interaction_scale: Scale parameter for pairwise interaction priors.
 *   - main_difference_scale: Scale parameter for main difference priors.
 *   - pairwise_difference_scale: Scale parameter for pairwise difference priors.
 *   - prior_threshold_alpha: Shape parameter for the threshold prior distribution.
 *   - prior_threshold_beta: Rate parameter for the threshold prior distribution.
 *   - rm_adaptation_rate: Robbins-Monro learning rate parameter.
 *   - target_acceptance_rate: Target acceptance probability for Metropolis-Hastings updates.
 *   - t: Current iteration number.
 *   - rm_lower_bound: Lower bound for Robbins-Monro updates.
 *   - rm_upper_bound: Upper bound for Robbins-Monro updates.
 *   - difference_selection: Boolean indicating whether to perform difference selection.
 *   - num_pairwise: Total number of pairwise effects.
 *   - num_variables: Total number of variables in the analysis.
 *
 * Outputs:
 *   - Updated model parameters, including:
 *       - `inclusion_indicator`: Updated inclusion indicators for pairwise differences.
 *       - `main_effects`: Updated main effects matrix.
 *       - `pairwise_effects`: Updated pairwise effects matrix.
 *       - `residual_matrix`: Updated residual scores.
 *       - `proposal_sd_main_effects`: Updated proposal standard deviations for main effects.
 *       - `proposal_sd_pairwise_effects`: Updated proposal standard deviations for pairwise effects.
 *
 * Steps:
 *   1. Update pairwise interaction parameters using Metropolis-Hastings.
 *   2. Update the selection of pairwise differences using Metropolis-Hastings (between model).
 *   3. Update pairwise differences using Metropolis-Hastings (within model).
 *   4. Update thresholds (main effects) for each variable:
 *       - If `independent_thresholds`, update thresholds separately for each group.
 *       - Otherwise, model group differences in thresholds.
 */
List gibbs_step_gm(
    arma::mat& main_effects,
    const arma::imat& main_effect_indices,
    arma::mat& pairwise_effects,
    const arma::imat& pairwise_effect_indices,
    const arma::mat& projection,
    const arma::imat& num_categories,
    const arma::imat& observations,
    const int num_persons,
    const int num_groups,
    const arma::imat& group_indices,
    const List& n_cat_obs,
    const List& sufficient_blume_capel,
    arma::mat& residual_matrix,
    const bool independent_thresholds,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    arma::imat& inclusion_indicator,
    const arma::mat& inclusion_probability_difference,
    arma::mat& proposal_sd_main_effects,
    arma::mat& proposal_sd_pairwise_effects,
    const double interaction_scale,
    const double main_difference_scale,
    const double pairwise_difference_scale,
    const double prior_threshold_alpha,
    const double prior_threshold_beta,
    const double rm_adaptation_rate,
    const double target_acceptance_rate,
    const int t,
    const double rm_lower_bound,
    const double rm_upper_bound,
    const bool difference_selection,
    const int num_pairwise,
    const int num_variables,
    const arma::imat& index
) {

  arma::mat proposal_sd_blumecapel_group(num_variables, 2);

  // Step 1: Update pairwise interaction parameters
  metropolis_interaction(
    main_effects, pairwise_effects, main_effect_indices,
    pairwise_effect_indices, projection, observations, num_groups, group_indices,
    num_categories, independent_thresholds, num_persons, residual_matrix,
    is_ordinal_variable, baseline_category, proposal_sd_pairwise_effects,
    interaction_scale, num_variables, rm_adaptation_rate, target_acceptance_rate, t, rm_lower_bound, rm_upper_bound);


  //  Step 2: Update the selection of pairwise differences
  if(difference_selection) {
    metropolis_pairwise_difference_between_model(
      inclusion_probability_difference, index, main_effects, pairwise_effects,
      main_effect_indices, pairwise_effect_indices, projection, observations, num_groups,
      group_indices, num_categories, independent_thresholds, inclusion_indicator, num_persons,
      residual_matrix, is_ordinal_variable, baseline_category, proposal_sd_pairwise_effects,
      pairwise_difference_scale, num_variables, rm_adaptation_rate,
      target_acceptance_rate, t, rm_lower_bound, rm_upper_bound, num_pairwise);
  }

  // Step 3: Update pairwise differences
  metropolis_pairwise_difference(
    main_effects, pairwise_effects, main_effect_indices, pairwise_effect_indices, projection,
    observations, num_groups, group_indices, num_categories, independent_thresholds,
    inclusion_indicator, num_persons, residual_matrix, is_ordinal_variable, baseline_category,
    proposal_sd_pairwise_effects, pairwise_difference_scale, num_variables, rm_adaptation_rate, target_acceptance_rate, t, rm_lower_bound,
    rm_upper_bound);

  // Step 4: Update thresholds
  for (int variable = 0; variable < num_variables; variable++) {
    if (independent_thresholds) {
      // Independent thresholds: Update separately for each group
      for (int group = 0; group < num_groups; ++group) {
        if (is_ordinal_variable[variable]) {
          metropolis_thresholds_regular_free(
            main_effects, main_effect_indices, observations, num_groups, group_indices,
            num_categories, residual_matrix, n_cat_obs, prior_threshold_alpha,
            prior_threshold_beta, variable, group);
        } else {
          metropolis_thresholds_blumecapel_free(
            main_effects, main_effect_indices, observations, num_groups, group_indices,
            baseline_category, num_categories, sufficient_blume_capel, num_persons,
            residual_matrix, n_cat_obs, prior_threshold_alpha, prior_threshold_beta, variable,
            group, proposal_sd_main_effects, rm_adaptation_rate, target_acceptance_rate, t,
            rm_lower_bound, rm_upper_bound);
        }
      }
    } else {
      // Thresholds modeled with group differences
      if (is_ordinal_variable[variable]) {
        metropolis_threshold_regular(
          main_effects, main_effect_indices, projection, observations, num_groups,
          group_indices, num_categories, num_persons, residual_matrix, n_cat_obs,
          prior_threshold_alpha, prior_threshold_beta, variable);

        if(difference_selection) {
          metropolis_main_difference_regular_between_model(
            inclusion_indicator, inclusion_probability_difference, main_effects,
            main_effect_indices, observations, num_groups, group_indices,
            num_categories, residual_matrix, n_cat_obs, variable,
            proposal_sd_main_effects, main_difference_scale, projection
          );
        }

        metropolis_main_difference_regular(
          main_effects, main_effect_indices, projection, observations, num_groups,
          group_indices, num_categories, num_persons, residual_matrix, n_cat_obs,
          variable, inclusion_indicator, proposal_sd_main_effects, main_difference_scale, rm_adaptation_rate,
          target_acceptance_rate, t, rm_lower_bound, rm_upper_bound);
      } else {
        metropolis_threshold_blumecapel(
          main_effects, main_effect_indices, projection, num_categories, sufficient_blume_capel,
          num_persons, num_groups, group_indices, variable, baseline_category,
          prior_threshold_alpha, prior_threshold_beta, residual_matrix, proposal_sd_main_effects, rm_adaptation_rate,
          target_acceptance_rate, t, rm_lower_bound, rm_upper_bound);

        if(difference_selection) {
          //...
        }

        metropolis_main_difference_blumecapel(
          main_effects, main_effect_indices, projection, num_categories, sufficient_blume_capel,
          num_persons, num_groups, group_indices, variable, baseline_category,
          prior_threshold_alpha, prior_threshold_beta, residual_matrix, inclusion_indicator,
          proposal_sd_main_effects, rm_adaptation_rate, target_acceptance_rate, t,
          rm_lower_bound, rm_upper_bound);
      }
    }
  }

  // Return updated parameters
  return List::create(Named("inclusion_indicator") = inclusion_indicator,
                      Named("main_effects") = main_effects,
                      Named("pairwise_effects") = pairwise_effects,
                      Named("residual_matrix") = residual_matrix,
                      Named("proposal_sd_main_effects") = proposal_sd_main_effects,
                      Named("proposal_sd_pairwise_effects") = proposal_sd_pairwise_effects);
}

/**
 * Function: compare_anova_gibbs_sampler
 * Purpose:
 *   Executes the Gibbs sampling process for Bayesian parameter comparisons across multiple groups.
 *   Updates main effects, pairwise effects, thresholds, and inclusion probabilities iteratively.
 *
 * Inputs:
 *   - observations: arma::imat of observed data (individuals by variables).
 *   - main_effect_indices: arma::imat mapping variables to their category indices.
 *   - pairwise_effect_indices: arma::imat mapping variable pairs to pairwise effect indices.
 *   - projection: arma::mat representing group-specific scaling.
 *   - num_categories: arma::imat containing category counts for each variable and group.
 *   - num_groups: Total number of groups in the analysis.
 *   - group_indices: arma::imat specifying group-wise start and end indices for individuals.
 *   - interaction_scale: Scale parameter for pairwise interaction priors.
 *   - pairwise_difference_scale: Scale parameter for pairwise difference priors.
 *   - main_difference_scale: Scale parameter for main difference priors.
 *   - pairwise_difference_prior: Type of prior for pairwise differences (e.g., "Beta-Bernoulli").
 *   - main_difference_prior: Type of prior for main differences (e.g., "Beta-Bernoulli").
 *   - inclusion_probability_difference: arma::mat for inclusion probabilities of differences.
 *   - pairwise_beta_bernoulli_alpha: Alpha parameter for Beta-Bernoulli prior on pairwise differences.
 *   - pairwise_beta_bernoulli_beta: Beta parameter for Beta-Bernoulli prior on pairwise differences.
 *   - main_beta_bernoulli_alpha: Alpha parameter for Beta-Bernoulli prior on main differences.
 *   - main_beta_bernoulli_beta: Beta parameter for Beta-Bernoulli prior on main differences.
 *   - Index: arma::imat for randomization during Gibbs sampling iterations.
 *   - iter: Total number of Gibbs sampling iterations.
 *   - burnin: Number of burn-in iterations.
 *   - n_cat_obs: List of category-specific observation counts for each group.
 *   - sufficient_blume_capel: List of sufficient statistics for Blume-Capel variables by group.
 *   - prior_threshold_alpha: Shape parameter for the threshold prior distribution.
 *   - prior_threshold_beta: Rate parameter for the threshold prior distribution.
 *   - na_impute: Boolean flag indicating whether to impute missing data.
 *   - missing_data_indices: arma::imat of indices for missing data entries.
 *   - is_ordinal_variable: arma::uvec indicating whether variables are ordinal.
 *   - baseline_category: arma::ivec specifying reference categories for each variable.
 *   - independent_thresholds: Boolean flag for whether thresholds are modeled independently across groups.
 *   - save: Boolean flag for saving results during Gibbs sampling (default: false).
 *   - display_progress: Boolean flag for displaying a progress bar during sampling (default: false).
 *   - difference_selection: Boolean flag for whether to perform difference selection (default: true).
 *
 * Outputs:
 *   - A List containing:
 *       - `main`: Matrix of updated main effects after sampling.
 *       - `pairwise`: Matrix of updated pairwise effects after sampling.
 *       - `inclusion_indicator`: Matrix of updated inclusion probabilities after sampling.
 *
 * Steps:
 *   1. Initialize matrices for main effects, pairwise effects, inclusion probabilities, and residuals.
 *   2. Iterate through burn-in and sampling iterations:
 *       a. Impute missing data if `na_impute` is true.
 *       b. Perform parameter updates via Gibbs sampling using helper functions for:
 *           - Pairwise interaction parameters.
 *           - Pairwise differences.
 *           - Thresholds (either independent or grouped).
 *       c. Update inclusion probabilities if difference selection is enabled.
 *       d. Save running averages of parameters after burn-in.
 *   3. Return the updated parameters.
 */
// [[Rcpp::export]]
List compare_anova_gibbs_sampler(
    arma::imat& observations,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::mat& projection,
    const arma::imat& num_categories,
    const int num_groups,
    const arma::imat& group_indices,
    const double interaction_scale,
    const double pairwise_difference_scale,
    const double main_difference_scale,
    const String& pairwise_difference_prior,
    const String& main_difference_prior,
    arma::mat& inclusion_probability_difference,
    const double pairwise_beta_bernoulli_alpha,
    const double pairwise_beta_bernoulli_beta,
    const double main_beta_bernoulli_alpha,
    const double main_beta_bernoulli_beta,
    const arma::imat& Index,
    const int iter,
    const int burnin,
    List& n_cat_obs,
    List& sufficient_blume_capel,
    const double prior_threshold_alpha,
    const double prior_threshold_beta,
    const bool na_impute,
    const arma::imat& missing_data_indices,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const bool independent_thresholds,
    const bool save_main = false,
    const bool save_pairwise = false,
    const bool save_indicator = false,
    const bool display_progress = false,
    bool difference_selection = true
) {

  // Step 1: Initialize parameters and helper matrices
  int num_variables = observations.n_cols;
  int num_persons = observations.n_rows;
  int num_main = 1 + main_effect_indices(num_variables - 1, 1);
  int num_pairwise = Index.n_rows;

  // Initialize person-group inclusion_indicator
  arma::ivec group_membership (num_persons);
  for(int group = 0; group < num_groups; group++) {
    for(int person = group_indices(group, 0); person < group_indices(group, 1) + 1; person++) {
      group_membership[person] = group;
    }
  }

  // Initialize model parameters
  arma::mat main_effects(num_main, num_groups);
  arma::mat pairwise_effects(num_pairwise, num_groups);
  arma::imat inclusion_indicator(num_variables, num_variables);
  std::fill(inclusion_indicator.begin(), inclusion_indicator.end(), 1);

  // Adaptive Metropolis proposal standard deviations
  arma::mat proposal_sd_main_effects(num_main, num_groups);
  arma::mat proposal_sd_pairwise_effects(num_pairwise, num_groups);
  std::fill(proposal_sd_main_effects.begin(), proposal_sd_main_effects.end(), 1.0);
  std::fill(proposal_sd_pairwise_effects.begin(), proposal_sd_pairwise_effects.end(), 1.0);

  // Robbins-Monro parameters
  double rm_adaptation_rate = 0.75, target_acceptance_rate = 0.234;
  int rm_scale = group_indices(0, 1) - group_indices(0, 0);
  for(int gr = 1; gr < num_groups; gr++) {
    int tmp = group_indices(gr, 1) - group_indices(gr, 0);
    if(tmp > rm_scale)
      rm_scale = tmp;
  }

  double rm_lower_bound = 1.0 / static_cast<double>(rm_scale);
  double rm_upper_bound = 2.0;

  //Randomized index for the pairwise updates ----------------------------------
  arma::ivec v = seq(0, num_pairwise - 1);
  arma::ivec order(num_pairwise);
  arma::imat index(num_pairwise, 3);

  // Rest matrix for pseudo-likelihoods
  arma::mat residual_matrix(num_persons, num_variables);

  // Output matrices
  arma::mat posterior_mean_main(num_main, num_groups);
  arma::mat posterior_mean_pairwise(num_pairwise, num_groups);
  arma::mat posterior_mean_indicator(num_variables, num_variables);
  std::fill(posterior_mean_indicator.begin(), posterior_mean_indicator.end(), 1);

  // Allocate matrices conditionally to save memory
  arma::mat* main_effect_samples = nullptr;
  arma::mat* pairwise_effect_samples = nullptr;
  arma::mat* inclusion_indicator_samples = nullptr;

  // Conditionally allocate based on save flags
  if (save_main) {
    main_effect_samples = new arma::mat(iter, num_main * num_groups);
  }
  if (save_pairwise) {
    pairwise_effect_samples = new arma::mat(iter, num_pairwise * num_groups);
  }
  if (save_indicator) {
    inclusion_indicator_samples = new arma::mat(iter, num_pairwise + num_variables);
  }

  // For difference selection we use a burnin in two phases
  bool enable_difference_selection = difference_selection; // Store the original difference selection input
  int total_burnin = burnin * (enable_difference_selection ? 2 : 1); // Compute the total burn-in duration
  difference_selection = false; // Flag to enable difference selection after the initial burn-in phase

  // Progress bar
  Progress p(iter + total_burnin, display_progress);

  // Step 2: Gibbs sampling loop
  for (int iteration = 0; iteration < iter + total_burnin; ++iteration) {
    if (Progress::check_abort()) {
      return List::create(Named("main") = posterior_mean_main,
                          Named("pairwise") = posterior_mean_pairwise,
                          Named("inclusion_indicator") = posterior_mean_indicator);
    }
    p.increment();

    // Enable difference selection at the midpoint of burn-in
    if (enable_difference_selection && iteration == burnin) {
      difference_selection = true;
    }

    // Create a random ordering of pairwise effects for updating
    arma::uvec order = arma::randperm(v.n_elem);
    //sample(v, num_pairwise, false, R_NilValue);

    for(int cntr = 0; cntr < num_pairwise; cntr++) {
      index(cntr, 0) = Index(order[cntr], 0);
      index(cntr, 1) = Index(order[cntr], 1);
      index(cntr, 2) = Index(order[cntr], 2);
    }

    // Handle missing data if required
    if (na_impute) {
      List impute_out = impute_missing_data_for_anova_model (
        main_effects, pairwise_effects, main_effect_indices, pairwise_effect_indices,
        projection, observations, num_groups, group_membership,
        n_cat_obs, sufficient_blume_capel, num_categories, residual_matrix,
        missing_data_indices, is_ordinal_variable, baseline_category, independent_thresholds);

      arma::imat observations_tmp = impute_out["observations"];
      List n_cat_obs_tmp = impute_out["n_cat_obs"];
      List sufficient_blume_capel_tmp = impute_out["sufficient_blume_capel"];
      arma::mat residual_matrix_tmp = impute_out["residual_matrix"];

      // Reassign to original variables
      observations = observations_tmp;
      n_cat_obs = n_cat_obs_tmp;
      sufficient_blume_capel = sufficient_blume_capel_tmp;
      residual_matrix = residual_matrix_tmp;
    }

    // Perform parameter updates
    List step_out = gibbs_step_gm (
      main_effects, main_effect_indices, pairwise_effects, pairwise_effect_indices, projection,
      num_categories, observations, num_persons, num_groups, group_indices, n_cat_obs,
      sufficient_blume_capel, residual_matrix, independent_thresholds, is_ordinal_variable,
      baseline_category, inclusion_indicator, inclusion_probability_difference,
      proposal_sd_main_effects, proposal_sd_pairwise_effects, interaction_scale,
      main_difference_scale, pairwise_difference_scale, prior_threshold_alpha,
      prior_threshold_beta, rm_adaptation_rate, target_acceptance_rate, iteration, rm_lower_bound, rm_upper_bound,
      difference_selection, num_pairwise, num_variables, index);

    arma::imat indicator_tmp = step_out["inclusion_indicator"];
    arma::mat main_effects_tmp = step_out["main_effects"];
    arma::mat pairwise_effects_tmp = step_out["pairwise_effects"];
    arma::mat residual_matrix_tmp = step_out["residual_matrix"];
    arma::mat proposal_sd_pairwise_tmp = step_out["proposal_sd_pairwise_effects"];
    arma::mat proposal_sd_main_tmp = step_out["proposal_sd_main_effects"];

    // Update matrices from Gibbs step output
    inclusion_indicator = indicator_tmp;
    main_effects = main_effects_tmp;
    pairwise_effects = pairwise_effects_tmp;
    residual_matrix = residual_matrix_tmp;
    proposal_sd_pairwise_effects = proposal_sd_pairwise_tmp;
    proposal_sd_main_effects = proposal_sd_main_tmp;

    // Update inclusion probabilities
    if (difference_selection) {
      int sumG = 0;

      if (pairwise_difference_prior == "Beta-Bernoulli") {
        // Update pairwise inclusion probabilities
        for (int i = 0; i < num_variables - 1; ++i) {
          for (int j = i + 1; j < num_variables; ++j) {
            sumG += inclusion_indicator(i, j);
          }
        }
        double prob = R::rbeta(pairwise_beta_bernoulli_alpha + sumG,
                               pairwise_beta_bernoulli_beta + num_pairwise - sumG);
        std::fill(inclusion_probability_difference.begin(), inclusion_probability_difference.end(), prob);
      }
    }

    // Save sampled states after burn-in if the respective saving option is enabled
    if (iteration >= total_burnin) {
      int iter_adj = iteration - total_burnin + 1;
      for (int col = 0; col < num_groups; ++col) {
        for (int row = 0; row < num_main; ++row) {
          posterior_mean_main(row, col) = (posterior_mean_main(row, col) * (iter_adj - 1) + main_effects(row, col)) /
            static_cast<double>(iter_adj);
        }
      }
      for (int col = 0; col < num_groups; ++col) {
        for (int row = 0; row < num_pairwise; ++row) {
          posterior_mean_pairwise(row, col) = (posterior_mean_pairwise(row, col) * (iter_adj - 1) + pairwise_effects(row, col)) /
            static_cast<double>(iter_adj);
        }
      }
      if(difference_selection) {
        for (int i = 0; i < num_variables - 1; ++i) {
          for (int j = i; j < num_variables; ++j) {
            posterior_mean_indicator(i, j) = (posterior_mean_indicator(i, j) * (iter_adj - 1) + inclusion_indicator(i, j)) /
              static_cast<double>(iter_adj);
            posterior_mean_indicator(j, i) = posterior_mean_indicator(i, j);
          }
        }
        int i = num_variables - 1;
        posterior_mean_indicator(i, i) = (posterior_mean_indicator(i, i) * (iter_adj - 1) + inclusion_indicator(i, i)) /
          static_cast<double>(iter_adj);
      }
      if(save_main) {
        int cntr = 0;
        for (int col = 0; col < num_groups; ++col) {
          for (int row = 0; row < num_main; ++row) {
            (*main_effect_samples)(iter_adj - 1, cntr) = main_effects(row, col);
            cntr++;
          }
        }
      }
      if(save_pairwise) {
        int cntr = 0;
        for (int col = 0; col < num_groups; ++col) {
          for (int row = 0; row < num_pairwise; ++row) {
            (*pairwise_effect_samples)(iter_adj - 1, cntr) = pairwise_effects(row, col);
            cntr++;
          }
        }
      }
      if(difference_selection) {
        if(save_indicator) {
          int cntr = 0;
          for (int i = 0; i < num_variables - 1; ++i) {
            for (int j = i; j < num_variables; ++j) {
              (*inclusion_indicator_samples)(iter_adj - 1, cntr) = inclusion_indicator(i, j);
              cntr++;
            }
          }
        }
      }
    }
  }

  // Compile the output based on saving options
  List output = List::create(Named("posterior_mean_main") = posterior_mean_main,
                             Named("posterior_mean_pairwise") = posterior_mean_pairwise,
                             Named("posterior_mean_indicator") = posterior_mean_indicator);

  // Cleanup at the end
  if (save_main) {
    output["main_effect_samples"] = *main_effect_samples;
    delete main_effect_samples;  // Free allocated memory
  }
  if (save_pairwise) {
    output["pairwise_effect_samples"] = *pairwise_effect_samples;
    delete pairwise_effect_samples;
  }
  if (save_indicator) {
    output["inclusion_indicator_samples"] = *inclusion_indicator_samples;
    delete inclusion_indicator_samples;
  }

  return output;
}
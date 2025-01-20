// [[Rcpp::depends(RcppProgress)]]
#include <Rcpp.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include "gibbs_functions.h"
using namespace Rcpp;


// ----------------------------------------------------------------------------|
// List of to do's
// 1. Add between model moves.
// 2. Use "index" to vary between-model updates to pairwise effects.
// 3. Add g prior.
// 4. Output matrices (resizing based on ``save'' and ``independent_thresholds'')
// 5. Add "enum" to handle threshold scenarios
// ----------------------------------------------------------------------------|

/**
 * Function: robbins_monro_update
 * Purpose: Performs Robbins-Monro updates for proposal standard deviations.
 * Inputs:
 *  - current_sd: Current standard deviation of the proposal.
 *  - observed_log_acceptance_probability: Log acceptance probability from the Metropolis-Hastings step.
 *  - target_acceptance_probability: Target acceptance rate.
 *  - t: Iteration number.
 *  - phi: Robbins-Monro learning rate.
 *  - epsilon_lo: Minimum allowable standard deviation.
 *  - epsilon_hi: Maximum allowable standard deviation.
 * Returns:
 *  - Updated proposal standard deviation, clamped within bounds.
 */
double robbins_monro_update(double current_sd,
                            double observed_log_acceptance_probability,
                            double target_acceptance_probability,
                            double epsilon_lo,
                            double epsilon_hi,
                            double rm_weight) {

  // Normalize the acceptance probability
  double observed_acceptance_probability = 1.0;
  if(observed_log_acceptance_probability < 0) {
    observed_acceptance_probability = std::exp(observed_log_acceptance_probability);
  }

  // Update the proposal standard deviation
  double update = current_sd +
    (observed_acceptance_probability - target_acceptance_probability) * rm_weight;

  // Handle NaN cases by resetting to default value
  if (std::isnan(update)) {
    update = 1.0; // Default proposal standard deviation
  }

  // Clamp the updated standard deviation within bounds
  return std::clamp(update, epsilon_lo, epsilon_hi);
}


/**
 * Function: group_thresholds_for_variable
 * Purpose: Computes thresholds for a specific variable in a given group.
 * Inputs:
 *  - variable: Index of the variable for which thresholds are computed.
 *  - ordinal_variable: Logical vector indicating if variables are ordinal.
 *  - group: Index of the group for which thresholds are computed.
 *  - no_groups: Total number of groups in the analysis.
 *  - no_categories: Matrix of category counts per variable and group.
 *  - main_effects: Matrix of main effects across variables and groups.
 *  - main_index: Indices for main effect parameters.
 *  - projection: Projection matrix for group-specific scaling.
 *  - independent_thresholds: Whether thresholds are modeled independently.
 * Outputs:
 *  - A `NumericVector` containing the computed thresholds for the variable in the group.
 */
NumericVector group_thresholds_for_variable (int variable,
                                             LogicalVector ordinal_variable,
                                             int group,
                                             int no_groups,
                                             IntegerMatrix no_categories,
                                             NumericMatrix main_effects,
                                             IntegerMatrix main_index,
                                             NumericMatrix projection,
                                             bool independent_thresholds) {
  // Determine the length of the threshold vector based on whether the variable is ordinal
  int vector_length = (ordinal_variable[variable]) ? no_categories(variable, group) : 2;
  NumericVector GroupThresholds(vector_length);

  // Base index for accessing main effects for this variable
  int base_category_index = main_index(variable, 0);

  if (!independent_thresholds) {
    // Compute thresholds for dependent case (i.e., model group differences)
    int n_cats = no_categories(variable, group);
    if(ordinal_variable[variable]) {
      // Regular binary or ordinal variable
      for (int category = 0; category < n_cats; category++) {
        int category_index = base_category_index + category;
        double threshold = main_effects(category_index, 0);

        for (int h = 0; h < no_groups - 1; h++) {
          threshold += projection(group, h) * main_effects(category_index, h + 1);
        }

        GroupThresholds[category] = threshold;
      }
    } else {
      // Blume-Capel ordinal variable
      double threshold0 = main_effects(base_category_index, 0);
      double threshold1 = main_effects(base_category_index + 1, 0);

      for (int h = 0; h < no_groups - 1; h++) {
        threshold0 += projection(group, h) * main_effects(base_category_index, h + 1);
        threshold1 += projection(group, h) * main_effects(base_category_index + 1, h + 1);
      }

      GroupThresholds[0] = threshold0;
      GroupThresholds[1] = threshold1;
    }
  } else {
    // Compute thresholds for groups independently
    if (ordinal_variable[variable]) {
      // Regular binary or ordinal variable
      int n_cats = no_categories(variable, group);
      for (int category = 0; category < n_cats; category++) {
        int category_index = base_category_index + category;
        GroupThresholds[category] = main_effects(category_index, group);
      }
    } else {
      // Blume-Capel ordinal variable
      GroupThresholds[0] = main_effects(base_category_index, group);
      GroupThresholds[1] = main_effects(base_category_index + 1, group);
    }
  }

  return GroupThresholds;
}


/**
 * Function: compare_anova_impute_missing_data
 * Purpose: Imputes missing data for independent samples designs by generating new observations
 *          based on the model parameters and pseudo-likelihood.
 *
 * Inputs:
 *  - main_effects: Numeric matrix of main effects across variables and groups.
 *  - pairwise_effects: Numeric matrix of pairwise interaction effects between variables across groups.
 *  - main_index: Integer matrix mapping variable indices to main effect parameters.
 *  - pairwise_index: Integer matrix mapping variable pairs to pairwise interaction parameters.
 *  - projection: Numeric matrix representing group-specific scaling for effects.
 *  - observations: Integer matrix of observed data (individuals x variables), with missing data encoded.
 *  - no_groups: Number of groups in the analysis.
 *  - person_group_indicator: Integer vector mapping individuals to their respective groups.
 *  - n_cat_obs: List of matrices, one per group, recording category frequencies per variable.
 *  - sufficient_blume_capel: List of matrices storing sufficient statistics for Blume-Capel variables.
 *  - no_categories: Integer matrix of category counts for each variable and group.
 *  - rest_matrix: Numeric matrix of residual effects for pseudo-likelihood calculations.
 *  - missing_index: Integer matrix of indices indicating missing observations (row x column pairs).
 *  - ordinal_variable: Logical vector indicating whether variables are ordinal.
 *  - reference_category: Integer vector of reference categories for Blume-Capel variables.
 *  - independent_thresholds: Boolean indicating whether thresholds are modeled independently.
 *
 * Outputs:
 *  - A List containing:
 *    - `observations`: Updated observation matrix with imputed values.
 *    - `n_cat_obs`: Updated list of category counts per group.
 *    - `sufficient_blume_capel`: Updated sufficient statistics for Blume-Capel variables.
 *    - `rest_matrix`: Updated residual effects matrix.
 */
List compare_anova_impute_missing_data(NumericMatrix main_effects,
                                       NumericMatrix pairwise_effects,
                                       IntegerMatrix main_index,
                                       IntegerMatrix pairwise_index,
                                       NumericMatrix projection,
                                       IntegerMatrix observations,
                                       int no_groups,
                                       IntegerVector person_group_indicator,
                                       List n_cat_obs,
                                       List sufficient_blume_capel,
                                       IntegerMatrix no_categories,
                                       NumericMatrix rest_matrix,
                                       IntegerMatrix missing_index,
                                       LogicalVector ordinal_variable,
                                       IntegerVector reference_category,
                                       bool independent_thresholds) {
  int no_variables = observations.ncol();
  int no_missings = missing_index.nrow();
  int max_no_categories = max(no_categories);

  NumericVector probabilities(max_no_categories + 1);
  double exponent, rest_score, cumsum, u;
  int score, person, variable, new_observation, old_observation;

  //Impute missing data
  for(int missing = 0; missing < no_missings; missing++) {
    // Identify the observation to impute
    person = missing_index(missing, 0);
    variable = missing_index(missing, 1);
    int gr = person_group_indicator[person];

    // Compute thresholds for the variable in the given group
    NumericVector GroupThresholds = group_thresholds_for_variable(variable,
                                                                  ordinal_variable,
                                                                  gr,
                                                                  no_groups,
                                                                  no_categories,
                                                                  main_effects,
                                                                  main_index,
                                                                  projection,
                                                                  independent_thresholds);

    // Generate a new observation based on the model
    rest_score = rest_matrix(person, variable);
    if(ordinal_variable[variable] == true) {
      // For regular binary or ordinal variables
      cumsum = 1.0;
      probabilities[0] = 1.0;
      for(int category = 1; category <= no_categories(variable, gr); category++) {
        exponent = GroupThresholds(category - 1);
        exponent += category * rest_score;
        cumsum += std::exp(exponent);
        probabilities[category] = cumsum;
      }
    } else {
      // For Blume-Capel variables
      cumsum = 0.0;
      for(int category = 0; category <= no_categories(variable, gr); category++) {
        exponent = GroupThresholds[0] * category;
        exponent += GroupThresholds[1] *
          (category - reference_category[variable]) *
          (category - reference_category[variable]);
        exponent += category * rest_score;
        cumsum += std::exp(exponent);
        probabilities[category] = cumsum;
      }
    }

    // Sample a new value based on computed probabilities
    u = cumsum * R::unif_rand();
    score = 0;
    while (u > probabilities[score]) {
      score++;
    }
    new_observation = score;
    old_observation = observations(person, variable);

    if(old_observation != new_observation) {
      // Update raw observations
      observations(person, variable) = new_observation;

      // Update category counts or sufficient statistics
      if(ordinal_variable[variable] == true) {
        IntegerMatrix n_cat_obs_gr = n_cat_obs[gr];
        n_cat_obs_gr(old_observation, variable)--;
        n_cat_obs_gr(new_observation, variable)++;
        n_cat_obs[gr] = n_cat_obs_gr;
      } else {
        IntegerMatrix sufficient_blume_capel_gr = sufficient_blume_capel[gr];
        sufficient_blume_capel_gr(0, variable) -= old_observation;
        sufficient_blume_capel_gr(0, variable) += new_observation;
        sufficient_blume_capel_gr(1, variable) -=
          (old_observation - reference_category[variable]) *
          (old_observation - reference_category[variable]);
        sufficient_blume_capel_gr(1, variable) +=
          (new_observation - reference_category[variable]) *
          (new_observation - reference_category[variable]);
        sufficient_blume_capel[gr] = sufficient_blume_capel_gr;
      }

      // Update rest scores
      for(int vertex = 0; vertex < no_variables; vertex++) {
        int int_index = pairwise_index(vertex, variable);
        double GroupInteraction = pairwise_effects(int_index, 0);
        for(int h = 0; h < no_groups - 1; h++) {
          GroupInteraction += projection(gr, h) * pairwise_effects(int_index, h + 1);
        }
        rest_matrix(person, vertex) -= old_observation * GroupInteraction;
        rest_matrix(person, vertex) += new_observation * GroupInteraction;
      }
    }
  }

  return List::create(Named("observations") = observations,
                      Named("n_cat_obs") = n_cat_obs,
                      Named("sufficient_blume_capel") = sufficient_blume_capel,
                      Named("rest_matrix") = rest_matrix);
}


/**
 * Function: compare_anova_log_pseudolikelihood_ratio_interaction
 * Purpose: Computes the log pseudo-likelihood ratio for proposed vs. current interaction states
 *          for a multi-group independent samples design.
 *
 * Inputs:
 *  - main_effects: Numeric matrix of main effects for variables across groups.
 *  - main_index: Integer matrix mapping variables to main effect parameters.
 *  - projection: Numeric matrix of group-specific scaling factors for interactions.
 *  - observations: Integer matrix of observed data (individuals x variables).
 *  - no_groups: Number of groups in the analysis.
 *  - group_index: Integer matrix containing start and end indices for each group in `observations`.
 *  - no_categories: Integer matrix indicating the number of categories for each variable and group.
 *  - independent_thresholds: Boolean indicating whether thresholds are modeled independently.
 *  - no_persons: Total number of individuals in the dataset.
 *  - variable1: Index of the first variable involved in the interaction.
 *  - variable2: Index of the second variable involved in the interaction.
 *  - proposed_state: Proposed value for the interaction parameter.
 *  - current_state: Current value of the interaction parameter.
 *  - rest_matrix: Numeric matrix of residual effects for pseudo-likelihood calculations.
 *  - ordinal_variable: Logical vector indicating whether each variable is ordinal.
 *  - reference_category: Integer vector indicating the reference categories for Blume-Capel variables.
 *
 * Outputs:
 *  - Returns the log pseudo-likelihood ratio for the proposed vs. current interaction states.
 */
double compare_anova_log_pseudolikelihood_ratio_interaction(NumericMatrix main_effects,
                                                            IntegerMatrix main_index,
                                                            NumericMatrix projection,
                                                            IntegerMatrix observations,
                                                            int no_groups,
                                                            IntegerMatrix group_index,
                                                            IntegerMatrix no_categories,
                                                            bool independent_thresholds,
                                                            int no_persons,
                                                            int variable1,
                                                            int variable2,
                                                            double proposed_state,
                                                            double current_state,
                                                            NumericMatrix rest_matrix,
                                                            LogicalVector ordinal_variable,
                                                            IntegerVector reference_category) {
  // Compute the difference in interaction states
  double delta_state = 2.0 * (proposed_state - current_state);
  double pseudolikelihood_ratio = 0.0;

  // Loop over groups to compute contributions to the pseudo-likelihood ratio
  for(int gr = 0; gr < no_groups; gr++) {
    // Retrieve thresholds for the two variables in the current group
    NumericVector GroupThresholds_v1 = group_thresholds_for_variable (
      variable1, ordinal_variable, gr, no_groups, no_categories, main_effects,
      main_index, projection, independent_thresholds);

    NumericVector GroupThresholds_v2 = group_thresholds_for_variable (
      variable2, ordinal_variable, gr, no_groups, no_categories, main_effects,
      main_index, projection, independent_thresholds);

    int n_cats_v1 = no_categories(variable1, gr);
    int n_cats_v2 = no_categories(variable2, gr);

    // Loop over individuals in the current group
    for(int person = group_index(gr, 0); person <= group_index(gr, 1); person++) {
      int obs_score1 = observations(person, variable1);
      int obs_score2 = observations(person, variable2);

      // Contribution from the interaction term
      pseudolikelihood_ratio += obs_score1 * obs_score2 * delta_state;

      // Compute contributions for each variable (variable1 and variable2)
      for (int variable = 1; variable <= 2; variable++) {
        // Assign variable-specific data
        int var = (variable == 1) ? variable1 : variable2;
        int obs = (variable == 1) ? obs_score2 : obs_score1;
        double obs_current = obs * current_state;
        double obs_proposed = obs * proposed_state;

        int n_cats = (variable == 1) ? n_cats_v1 : n_cats_v2;
        NumericVector& GroupThresholds = (variable == 1) ? GroupThresholds_v1 : GroupThresholds_v2;

        // Compute the rest score and bounds
        double rest_score = rest_matrix(person, var) - obs_current;
        double bound = (rest_score > 0) ? n_cats * rest_score : 0.0;
        double denominator_prop = 0.0;
        double denominator_curr = 0.0;

        // Compute pseudo-likelihood terms
        if(ordinal_variable[var]) {
          // Regular binary or ordinal MRF variable
          denominator_prop += std::exp(-bound);
          denominator_curr += std::exp(-bound);
          for(int cat = 0; cat < n_cats; cat++) {
            int score = cat + 1;
            double exponent = GroupThresholds[cat] + score * rest_score - bound;
            denominator_prop += std::exp(exponent + score * obs_proposed);
            denominator_curr += std::exp(exponent + score * obs_current);
          }
        } else {
          // Blume-Capel ordinal MRF variable
          for(int cat = 0; cat <= n_cats; cat++) {
            double exponent = GroupThresholds[0] * cat;
            exponent += GroupThresholds[1] *
              (cat - reference_category[var]) *
              (cat - reference_category[var]);
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
 * Function: compare_anova_metropolis_interaction
 * Purpose: Uses the Metropolis-Hastings algorithm to sample from the full conditional distribution
 *          of a pairwise interaction parameter, updating its value and the associated proposal standard deviation.
 *
 * Inputs:
 *  - main_effects: Numeric matrix of main effects across variables and groups.
 *  - pairwise_effects: Numeric matrix of pairwise effects across variable pairs and groups.
 *  - main_index: Integer matrix mapping variables to main effect parameters.
 *  - pairwise_index: Integer matrix mapping variable pairs to pairwise interaction parameters.
 *  - projection: Numeric matrix representing group-specific scaling factors for interactions.
 *  - observations: Integer matrix of observed data (individuals x variables).
 *  - no_groups: Number of groups in the analysis.
 *  - group_index: Integer matrix indicating the start and end indices for each group in `observations`.
 *  - no_categories: Integer matrix indicating the number of categories for each variable and group.
 *  - independent_thresholds: Boolean indicating whether thresholds are modeled independently.
 *  - no_persons: Total number of individuals in the dataset.
 *  - rest_matrix: Numeric matrix of residual effects for pseudo-likelihood calculations.
 *  - ordinal_variable: Logical vector indicating whether each variable is ordinal.
 *  - reference_category: Integer vector indicating the reference categories for Blume-Capel variables.
 *  - proposal_sd_pairwise: Numeric matrix of proposal standard deviations for pairwise interaction parameters.
 *  - interaction_scale: Scale parameter for the prior distribution of interaction terms.
 *  - no_variables: Total number of variables in the analysis.
 *  - phi: Robbins-Monro learning rate for adapting proposal standard deviations.
 *  - target_acceptance_probability: Target log acceptance rate for the Metropolis-Hastings updates.
 *  - t: Current iteration number in the sampling procedure.
 *  - epsilon_lo: Lower bound for the proposal standard deviations.
 *  - epsilon_hi: Upper bound for the proposal standard deviations.
 *
 * Outputs:
 *  - Updates `pairwise_effects` with new sampled values.
 *  - Updates `rest_matrix` with adjusted residual effects.
 *  - Updates `proposal_sd_pairwise` with adjusted proposal standard deviations using Robbins-Monro updates.
 */
void compare_anova_metropolis_interaction(NumericMatrix main_effects,
                                          NumericMatrix pairwise_effects,
                                          IntegerMatrix main_index,
                                          IntegerMatrix pairwise_index,
                                          NumericMatrix projection,
                                          IntegerMatrix observations,
                                          int no_groups,
                                          IntegerMatrix group_index,
                                          IntegerMatrix no_categories,
                                          bool independent_thresholds,
                                          int no_persons,
                                          NumericMatrix rest_matrix,
                                          LogicalVector ordinal_variable,
                                          IntegerVector reference_category,
                                          NumericMatrix proposal_sd_pairwise,
                                          double interaction_scale,
                                          int no_variables,
                                          double phi,
                                          double target_acceptance_probability,
                                          int t,
                                          double epsilon_lo,
                                          double epsilon_hi) {
  double log_prob;
  double exp_neg_log_t_phi = std::exp(-std::log(t) * phi);// Precompute Robbins-Monro decay term

  // Iterate over all pairs of variables for interaction updates
  for(int variable1 = 0; variable1 <  no_variables - 1; variable1++) {
    for(int variable2 = variable1 + 1; variable2 <  no_variables; variable2++) {
      int int_index = pairwise_index(variable1, variable2);

      // Retrieve the current state of the interaction parameter
      double current_state = pairwise_effects(int_index, 0);

      // Propose a new state from a normal density centered at the current state
      double proposed_state = R::rnorm(current_state, proposal_sd_pairwise(int_index, 0));

      // Compute the log pseudo-likelihood ratio for proposed vs. current states
      log_prob = compare_anova_log_pseudolikelihood_ratio_interaction(
        main_effects, main_index, projection, observations, no_groups,
        group_index, no_categories, independent_thresholds, no_persons,
        variable1, variable2, proposed_state, current_state, rest_matrix,
        ordinal_variable, reference_category);

      // Add prior probabilities for the interaction parameter
      log_prob += R::dcauchy(proposed_state, 0.0, interaction_scale, true);
      log_prob -= R::dcauchy(current_state, 0.0, interaction_scale, true);

      // Metropolis-Hastings acceptance step
      double U = R::unif_rand();
      if(std::log(U) < log_prob) {
        // Update the interaction parameter
        double state_difference = proposed_state - current_state;
        pairwise_effects(int_index, 0) = proposed_state;

        // Update the rest matrix to reflect the new interaction parameter
        for(int person = 0; person < no_persons; person++) {
          double obs1 = observations(person, variable1);
          double obs2 = observations(person, variable2);
          rest_matrix(person, variable1) += obs2 * state_difference;
          rest_matrix(person, variable2) += obs1 * state_difference;
        }
      }

      // Robbins-Monro update to adapt the proposal standard deviation
      proposal_sd_pairwise(int_index, 0) = robbins_monro_update(
        proposal_sd_pairwise(int_index, 0), log_prob, target_acceptance_probability, epsilon_lo,
        epsilon_hi, exp_neg_log_t_phi);
    }
  }
}


/**
 * Function: compare_anova_log_pseudolikelihood_ratio_pairwise_difference
 * Purpose: Computes the log pseudo-likelihood ratio for proposed vs. current
 *          values of a pairwise interaction parameter difference in a multi-group
 *          independent samples design.
 *
 * Inputs:
 *  - main_effects: Numeric matrix of main effects for variables across groups.
 *  - main_index: Integer matrix mapping variables to main effect parameters.
 *  - projection: Numeric matrix representing group-specific scaling factors for interactions.
 *  - observations: Integer matrix of observed data (individuals x variables).
 *  - no_groups: Number of groups in the analysis.
 *  - group_index: Integer matrix indicating the start and end indices for each group in `observations`.
 *  - no_categories: Integer matrix indicating the number of categories for each variable and group.
 *  - independent_thresholds: Boolean indicating whether thresholds are modeled independently.
 *  - no_persons: Total number of individuals in the dataset.
 *  - variable1: Index of the first variable in the interaction pair.
 *  - variable2: Index of the second variable in the interaction pair.
 *  - h: Index of the group difference being modeled.
 *  - proposed_state: Proposed value for the interaction parameter difference.
 *  - current_state: Current value of the interaction parameter difference.
 *  - rest_matrix: Numeric matrix of residual effects for pseudo-likelihood calculations.
 *  - ordinal_variable: Logical vector indicating whether each variable is ordinal.
 *  - reference_category: Integer vector indicating the reference categories for Blume-Capel variables.
 *
 * Outputs:
 *  - Returns the log pseudo-likelihood ratio for the proposed vs. current values
 *    of the pairwise interaction parameter difference.
 */
double compare_anova_log_pseudolikelihood_ratio_pairwise_difference(NumericMatrix main_effects,
                                                                    IntegerMatrix main_index,
                                                                    NumericMatrix projection,
                                                                    IntegerMatrix observations,
                                                                    int no_groups,
                                                                    IntegerMatrix group_index,
                                                                    IntegerMatrix no_categories,
                                                                    bool independent_thresholds,
                                                                    int no_persons,
                                                                    int variable1,
                                                                    int variable2,
                                                                    int h,
                                                                    double proposed_state,
                                                                    double current_state,
                                                                    NumericMatrix rest_matrix,
                                                                    LogicalVector ordinal_variable,
                                                                    IntegerVector reference_category) {
  double pseudolikelihood_ratio = 0.0;
  double delta_state = 2.0 * (proposed_state - current_state); // Change in parameter value

  // Iterate over all groups
  double denominator_prop, denominator_curr;

  // Loop over groups
  for(int gr = 0; gr < no_groups; gr++) {
    // Compute thresholds for both variables in the current group
    NumericVector GroupThresholds_v1 = group_thresholds_for_variable(
      variable1, ordinal_variable, gr, no_groups, no_categories, main_effects,
      main_index, projection, independent_thresholds);

    NumericVector GroupThresholds_v2 = group_thresholds_for_variable(
      variable2, ordinal_variable, gr, no_groups, no_categories, main_effects,
      main_index, projection, independent_thresholds);

    // Scale the delta_state by the projection factor for this group
    double P = projection(gr, h);
    double delta_state_group = delta_state * P;

    // Cache the number of categories for both variables
    int n_cats_v1 = no_categories(variable1, gr);
    int n_cats_v2 = no_categories(variable2, gr);

    // Iterate over all individuals in the current group
    for(int person = group_index(gr, 0); person <= group_index(gr, 1); person++) {
      // Cache observation scores and scaled terms
      int obs_score1 = observations(person, variable1);
      int obs_score2 = observations(person, variable2);
      double obs_proposed_p1 = obs_score2 * proposed_state * P;
      double obs_current_p1 = obs_score2 * current_state * P;
      double obs_proposed_p2 = obs_score1 * proposed_state * P;
      double obs_current_p2 = obs_score1 * current_state * P;

      // Contribution from the interaction term
      pseudolikelihood_ratio +=  obs_score1 * obs_score2 * delta_state_group;

      // Process each variable in the interaction pair
      for (int variable = 1; variable <= 2; variable++) {
        int var = (variable == 1) ? variable1 : variable2;
        int n_cats = (variable == 1) ? n_cats_v1 : n_cats_v2;
        NumericVector& GroupThresholds = (variable == 1) ? GroupThresholds_v1 : GroupThresholds_v2;
        double obs_proposed_p = (variable == 1) ? obs_proposed_p1 : obs_proposed_p2;
        double obs_current_p = (variable == 1) ? obs_current_p1 : obs_current_p2;

        // Compute the rest score and bound
        double rest_score = rest_matrix(person, var) - obs_current_p;
        double bound = (rest_score > 0) ? n_cats * rest_score : 0.0;

        // Compute denominators based on whether the variable is ordinal
        if (ordinal_variable[var]) {
          // Binary or ordinal MRF variable
          denominator_prop = std::exp(-bound);
          denominator_curr = std::exp(-bound);
          for (int cat = 0; cat < n_cats; cat++) {
            int score = cat + 1;
            double exponent = GroupThresholds[cat] + score * rest_score - bound;
            denominator_prop += std::exp(exponent + score * obs_proposed_p);
            denominator_curr += std::exp(exponent + score * obs_current_p);
          }
        } else {
          // Blume-Capel ordinal MRF variable
          denominator_prop = 0.0;
          denominator_curr = 0.0;
          for (int cat = 0; cat <= n_cats; cat++) {
            double exponent = GroupThresholds[0] * cat +
              GroupThresholds[1] * (cat - reference_category[var]) *
              (cat - reference_category[var]) +
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
 * Function: compare_anova_metropolis_pairwise_difference
 * Purpose: Uses the Metropolis-Hastings algorithm to sample from the full conditional
 *          distribution of pairwise interaction parameter differences for multiple groups,
 *          updating their values and the associated proposal standard deviations.
 *
 * Inputs:
 *  - main_effects: Numeric matrix of main effects across variables and groups.
 *  - pairwise_effects: Numeric matrix of pairwise interaction parameter differences across groups.
 *  - main_index: Integer matrix mapping variables to main effect parameters.
 *  - pairwise_index: Integer matrix mapping variable pairs to pairwise interaction parameters.
 *  - projection: Numeric matrix representing group-specific scaling factors for interactions.
 *  - observations: Integer matrix of observed data (individuals x variables).
 *  - no_groups: Number of groups in the analysis.
 *  - group_index: Integer matrix indicating the start and end indices for each group in `observations`.
 *  - no_categories: Integer matrix indicating the number of categories for each variable and group.
 *  - independent_thresholds: Boolean indicating whether thresholds are modeled independently.
 *  - indicator: Integer matrix indicating whether pairwise differences are included in the model.
 *  - no_persons: Total number of individuals in the dataset.
 *  - rest_matrix: Numeric matrix of residual effects for pseudo-likelihood calculations.
 *  - ordinal_variable: Logical vector indicating whether each variable is ordinal.
 *  - reference_category: Integer vector indicating the reference categories for Blume-Capel variables.
 *  - proposal_sd_pairwise: Numeric matrix of proposal standard deviations for pairwise interaction parameters.
 *  - pairwise_difference_scale: Scale parameter for the prior distribution of pairwise differences.
 *  - no_variables: Total number of variables in the analysis.
 *  - phi: Robbins-Monro learning rate for adapting proposal standard deviations.
 *  - target_acceptance_probability: Target log acceptance rate for the Metropolis-Hastings updates.
 *  - t: Current iteration number in the sampling procedure.
 *  - epsilon_lo: Lower bound for the proposal standard deviations.
 *  - epsilon_hi: Upper bound for the proposal standard deviations.
 *
 * Outputs:
 *  - Updates `pairwise_effects` with new sampled values.
 *  - Updates `rest_matrix` with adjusted residual effects.
 *  - Updates `proposal_sd_pairwise` with adjusted proposal standard deviations using Robbins-Monro updates.
 */
void compare_anova_metropolis_pairwise_difference(NumericMatrix main_effects,
                                                  NumericMatrix pairwise_effects,
                                                  IntegerMatrix main_index,
                                                  IntegerMatrix pairwise_index,
                                                  NumericMatrix projection,
                                                  IntegerMatrix observations,
                                                  int no_groups,
                                                  IntegerMatrix group_index,
                                                  IntegerMatrix no_categories,
                                                  bool independent_thresholds,
                                                  IntegerMatrix indicator,
                                                  int no_persons,
                                                  NumericMatrix rest_matrix,
                                                  LogicalVector ordinal_variable,
                                                  IntegerVector reference_category,
                                                  NumericMatrix proposal_sd_pairwise,
                                                  double pairwise_difference_scale,
                                                  int no_variables,
                                                  double phi,
                                                  double target_acceptance_probability,
                                                  int t,
                                                  double epsilon_lo,
                                                  double epsilon_hi) {
  double exp_neg_log_t_phi = std::exp(-std::log(t) * phi); // Precompute Robbins-Monro decay term

  // Iterate over all variable pairs
  for(int variable1 = 0; variable1 <  no_variables - 1; variable1++) {
    for(int variable2 = variable1 + 1; variable2 <  no_variables; variable2++) {
      if (indicator(variable1, variable2) == 1) {
        int int_index = pairwise_index(variable1, variable2);

        // Iterate over all groups (excluding the reference group)
        for(int h = 1; h < no_groups; h++) {
          // Retrieve the current state and propose a new state
          double current_state = pairwise_effects(int_index, h);
          double proposed_state = R::rnorm(current_state, proposal_sd_pairwise(int_index, h));

          // Compute the log pseudo-likelihood ratio for the proposed vs. current state
          double log_prob = compare_anova_log_pseudolikelihood_ratio_pairwise_difference(
            main_effects, main_index, projection, observations, no_groups,
            group_index, no_categories, independent_thresholds, no_persons,
            variable1, variable2, h - 1, proposed_state, current_state,
            rest_matrix, ordinal_variable, reference_category);

          // Add prior probabilities for the pairwise difference parameter
          log_prob += R::dcauchy(proposed_state, 0.0, pairwise_difference_scale, true);
          log_prob -= R::dcauchy(current_state, 0.0, pairwise_difference_scale, true);

          // Metropolis-Hastings acceptance step
          double U = R::unif_rand();
          if(std::log(U) < log_prob) {
            pairwise_effects(int_index, h) = proposed_state;

            // Update the rest matrix to reflect the new pairwise difference
            for(int gr = 0; gr < no_groups; gr++) {
              double state_difference = (proposed_state - current_state) * projection(gr, h - 1);
              for(int person = group_index(gr, 0); person <= group_index(gr, 1); person++) {
                int obs1 = observations(person, variable1);
                int obs2 = observations(person, variable2);
                rest_matrix(person, variable1) += obs2 * state_difference;
                rest_matrix(person, variable2) += obs1 * state_difference;
              }
            }
          }

          // Robbins-Monro update for the proposal standard deviation
          proposal_sd_pairwise(int_index, h) = robbins_monro_update(
            proposal_sd_pairwise(int_index, h), log_prob, target_acceptance_probability, epsilon_lo,
            epsilon_hi, exp_neg_log_t_phi);
        }
      }
    }
  }
}


/**
 * Function: compare_anova_log_pseudolikelihood_ratio_pairwise_differences
 * Purpose:
 *   Computes the log pseudo-likelihood ratio for proposed vs. current pairwise
 *   difference states in a Bayesian ANOVA.
 *
 * Inputs:
 *   - main_effects: NumericMatrix of main effects for all variables and groups.
 *   - main_index: IntegerMatrix mapping variables to category indices.
 *   - pairwise_effects: NumericMatrix of pairwise effects for all variable pairs and groups.
 *   - pairwise_index: IntegerMatrix mapping variable pairs to pairwise effect indices.
 *   - projection: NumericMatrix representing group-specific scaling.
 *   - observations: IntegerMatrix of observed data (individuals by variables).
 *   - no_groups: Total number of groups in the analysis.
 *   - group_index: IntegerMatrix specifying group-wise start and end indices for individuals.
 *   - no_categories: IntegerMatrix containing category counts for each variable and group.
 *   - independent_thresholds: Boolean flag for whether thresholds are modeled independently.
 *   - no_persons: Total number of individuals in the analysis.
 *   - variable1: Index of the first variable in the pair.
 *   - variable2: Index of the second variable in the pair.
 *   - proposed_states: NumericVector of proposed pairwise difference states for all groups.
 *   - current_states: NumericVector of current pairwise difference states for all groups.
 *   - rest_matrix: NumericMatrix of residuals used for pseudo-likelihood calculations.
 *   - ordinal_variable: LogicalVector indicating whether variables are ordinal.
 *   - reference_category: IntegerVector specifying reference categories for each variable.
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
double compare_anova_log_pseudolikelihood_ratio_pairwise_differences(NumericMatrix main_effects,
                                                                     IntegerMatrix main_index,
                                                                     NumericMatrix pairwise_effects,
                                                                     IntegerMatrix pairwise_index,
                                                                     NumericMatrix projection,
                                                                     IntegerMatrix observations,
                                                                     int no_groups,
                                                                     IntegerMatrix group_index,
                                                                     IntegerMatrix no_categories,
                                                                     bool independent_thresholds,
                                                                     int no_persons,
                                                                     int variable1,
                                                                     int variable2,
                                                                     NumericVector proposed_states,
                                                                     NumericVector current_states,
                                                                     NumericMatrix rest_matrix,
                                                                     LogicalVector ordinal_variable,
                                                                     IntegerVector reference_category) {

  double pseudolikelihood_ratio = 0.0;

  // Loop over groups
  for (int gr = 0; gr < no_groups; ++gr) {
    // Compute thresholds for both variables
    NumericVector GroupThresholds_v1 = group_thresholds_for_variable(
      variable1, ordinal_variable, gr, no_groups, no_categories, main_effects,
      main_index, projection, independent_thresholds);

    NumericVector GroupThresholds_v2 = group_thresholds_for_variable(
      variable2, ordinal_variable, gr, no_groups, no_categories, main_effects,
      main_index, projection, independent_thresholds);

    // Compute current and proposed group interactions
    double current_group_interaction = 0.0;
    double proposed_group_interaction = 0.0;
    for (int h = 0; h < no_groups - 1; h++) {
      current_group_interaction += current_states[h] * projection(gr, h);
      proposed_group_interaction += proposed_states[h] * projection(gr, h);
    }
    double delta_state_group = 2.0 * (proposed_group_interaction - current_group_interaction);

    // Cache the number of categories for both variables
    int n_cats_v1 = no_categories(variable1, gr);
    int n_cats_v2 = no_categories(variable2, gr);

    // Iterate over all individuals in the current group
    for (int person = group_index(gr, 0); person <= group_index(gr, 1); ++person) {
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
        NumericVector& GroupThresholds = (variable == 1) ? GroupThresholds_v1 : GroupThresholds_v2;
        double obs_proposed_p = (variable == 1) ? obs_proposed_p1 : obs_proposed_p2;
        double obs_current_p = (variable == 1) ? obs_current_p1 : obs_current_p2;


        double rest_score = rest_matrix(person, var) - obs_current_p;
        double bound = (rest_score > 0) ? n_cats * rest_score : 0.0;

        double denominator_curr = 0.0, denominator_prop = 0.0;

        // Compute denominators
        if (ordinal_variable[var]) {
          for (int cat = 0; cat < n_cats; cat++) {
            int score = cat + 1;
            double exponent = GroupThresholds[cat] + rest_score * score - bound;
            denominator_curr += std::exp(exponent + score * obs_current_p);
            denominator_prop += std::exp(exponent + score * obs_proposed_p);
          }
        } else {
          for (int cat = 0; cat <= n_cats; cat++) {
            double exponent = GroupThresholds[0] * cat +
                              GroupThresholds[1] * (cat - reference_category[var]) *
                              (cat - reference_category[var]) +
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
 * Function: compare_anova_metropolis_pairwise_difference_between_model
 * Purpose:
 *   Implements a between-model Metropolis-Hastings algorithm to update
 *   the inclusion indicator and pairwise differences for variable pairs.
 *
 * Inputs:
 *   - inclusion_probability_difference: NumericMatrix of inclusion probabilities
 *                                       for pairwise differences.
 *   - index: IntegerMatrix mapping pairwise differences to variable indices.
 *   - main_effects: NumericMatrix of main effects for all variables and groups.
 *   - pairwise_effects: NumericMatrix of pairwise effects for all variable pairs and groups.
 *   - main_index: IntegerMatrix mapping variables to category indices.
 *   - pairwise_index: IntegerMatrix mapping variable pairs to pairwise effect indices.
 *   - projection: NumericMatrix for group-specific scaling.
 *   - observations: IntegerMatrix of observed data (individuals by variables).
 *   - no_groups: Total number of groups in the analysis.
 *   - group_index: IntegerMatrix specifying group-wise start and end indices for individuals.
 *   - no_categories: IntegerMatrix containing category counts for each variable and group.
 *   - independent_thresholds: Boolean flag for whether thresholds are modeled independently.
 *   - indicator: IntegerMatrix indicating active pairwise differences.
 *   - no_persons: Total number of individuals in the analysis.
 *   - rest_matrix: NumericMatrix of residuals used for pseudo-likelihood calculations.
 *   - ordinal_variable: LogicalVector indicating whether variables are ordinal.
 *   - reference_category: IntegerVector specifying reference categories for each variable.
 *   - proposal_sd_pairwise: NumericMatrix of proposal standard deviations for pairwise differences.
 *   - pairwise_difference_scale: Double representing the scale of the prior distribution
 *                                for pairwise differences.
 *   - no_variables: Total number of variables in the analysis.
 *   - phi: Robbins-Monro learning rate parameter.
 *   - target_acceptance_probability: Target acceptance probability for Metropolis-Hastings updates.
 *   - t: Current iteration number.
 *   - epsilon_lo: Lower bound for proposal standard deviations.
 *   - epsilon_hi: Upper bound for proposal standard deviations.
 *   - no_pairwise: Total number of pairwise differences.
 *
 * Outputs:
 *   - Updates `indicator`, `pairwise_effects`, and `rest_matrix` to reflect
 *     accepted proposals for pairwise differences.
 */
void compare_anova_metropolis_pairwise_difference_between_model(NumericMatrix inclusion_probability_difference,
                                                                IntegerMatrix index,
                                                                NumericMatrix main_effects,
                                                                NumericMatrix pairwise_effects,
                                                                IntegerMatrix main_index,
                                                                IntegerMatrix pairwise_index,
                                                                NumericMatrix projection,
                                                                IntegerMatrix observations,
                                                                int no_groups,
                                                                IntegerMatrix group_index,
                                                                IntegerMatrix no_categories,
                                                                bool independent_thresholds,
                                                                IntegerMatrix indicator,
                                                                int no_persons,
                                                                NumericMatrix rest_matrix,
                                                                LogicalVector ordinal_variable,
                                                                IntegerVector reference_category,
                                                                NumericMatrix proposal_sd_pairwise,
                                                                double pairwise_difference_scale,
                                                                int no_variables,
                                                                double phi,
                                                                double target_acceptance_probability,
                                                                int t,
                                                                double epsilon_lo,
                                                                double epsilon_hi,
                                                                int no_pairwise) {
// Vectors to store current and proposed states for pairwise differences
  NumericVector proposed_states(no_groups - 1);
  NumericVector current_states(no_groups - 1);

  // Loop over all pairwise differences
  for (int cntr = 0; cntr < no_pairwise; ++cntr) {
    int variable1 = index(cntr, 1);
    int variable2 = index(cntr, 2);
    int int_index = pairwise_index(variable1, variable2);

    double log_prob = 0.0;

    // Loop over groups to process current and proposed states
    for (int h = 1; h < no_groups; h++) {
      double current_state = pairwise_effects(int_index, h);
      current_states[h - 1] = current_state;

      double proposed_state = 0.0;

      // Update log probabilities based on the inclusion indicator
      if (indicator(variable1, variable2) == 1) {
        // Difference is included
        log_prob -= R::dcauchy(current_state, 0.0, pairwise_difference_scale, true);
        log_prob += R::dnorm(current_state, proposed_state, proposal_sd_pairwise(int_index, h), true);
      } else {
        // Propose a new state
        proposed_state = R::rnorm(current_state, proposal_sd_pairwise(int_index, h));
        log_prob += R::dcauchy(proposed_state, 0.0, pairwise_difference_scale, true);
        log_prob -= R::dnorm(proposed_state, current_state, proposal_sd_pairwise(int_index, h), true);
      }

      proposed_states[h - 1] = proposed_state;
    }

    // Update log probability for the inclusion indicator
    if (indicator(variable1, variable2) == 1) {
      log_prob -= std::log(inclusion_probability_difference(variable1, variable2));
      log_prob += std::log(1 - inclusion_probability_difference(variable1, variable2));
    } else {
      log_prob += std::log(inclusion_probability_difference(variable1, variable2));
      log_prob -= std::log(1 - inclusion_probability_difference(variable1, variable2));
    }

    // Compute log pseudo-likelihood ratio
    log_prob += compare_anova_log_pseudolikelihood_ratio_pairwise_differences(
      main_effects, main_index, pairwise_effects, pairwise_index, projection, observations,
      no_groups, group_index, no_categories, independent_thresholds, no_persons, variable1,
      variable2, proposed_states, current_states, rest_matrix, ordinal_variable, reference_category);

    // Metropolis-Hastings acceptance step
    double U = R::unif_rand();
    if (std::log(U) < log_prob) {
      // Update inclusion indicator
      indicator(variable1, variable2) = 1 - indicator(variable1, variable2);
      indicator(variable2, variable1) = indicator(variable1, variable2);

      // Update pairwise effects and rest matrix
      for (int h = 1; h < no_groups; ++h) {
        pairwise_effects(int_index, h) = proposed_states[h - 1];
      }

      // Update residuals in the rest matrix
      for (int gr = 0; gr < no_groups; ++gr) {
        double state_difference = 0.0;
        for (int h = 0; h < no_groups - 1; ++h) {
          state_difference += (proposed_states[h] - current_states[h]) * projection(gr, h);
        }
        for (int person = group_index(gr, 0); person <= group_index(gr, 1); ++person) {
          int obs1 = observations(person, variable1);
          int obs2 = observations(person, variable2);
          rest_matrix(person, variable1) += obs2 * state_difference;
          rest_matrix(person, variable2) += obs1 * state_difference;
        }
      }
    }
  }
}


/**
 * Function: compare_anova_metropolis_threshold_regular
 * Purpose: Uses the Metropolis-Hastings algorithm to sample from the full conditional
 *          distribution of overall category threshold parameters for regular ordinal variables.
 *
 * Inputs:
 *  - main_effects: Numeric matrix containing the main effect parameters for all variables and groups.
 *  - main_index: Integer matrix mapping variables to their main effect parameter indices.
 *  - projection: Numeric matrix specifying group-specific scaling for thresholds.
 *  - observations: Integer matrix containing the observed data (individuals x variables).
 *  - no_groups: Total number of groups in the analysis.
 *  - group_index: Integer matrix specifying start and end indices for each group in `observations`.
 *  - no_categories: Integer matrix specifying the number of categories for each variable and group.
 *  - no_persons: Total number of individuals in the dataset.
 *  - rest_matrix: Numeric matrix storing residual effects for pseudo-likelihood calculations.
 *  - n_cat_obs: List of integer matrices, where each matrix tracks the number of observations
 *               for each category of a variable in a specific group.
 *  - threshold_alpha: Hyperparameter for the prior distribution of threshold parameters (shape parameter).
 *  - threshold_beta: Hyperparameter for the prior distribution of threshold parameters (rate parameter).
 *  - variable: Index of the variable whose thresholds are being updated.
 *
 * Outputs:
 *  - Updates the `main_effects` matrix with new sampled threshold parameters for the specified variable.
 */
void compare_anova_metropolis_threshold_regular(NumericMatrix main_effects,
                                                IntegerMatrix main_index,
                                                NumericMatrix projection,
                                                IntegerMatrix observations,
                                                int no_groups,
                                                IntegerMatrix group_index,
                                                IntegerMatrix no_categories,
                                                int no_persons,
                                                NumericMatrix rest_matrix,
                                                List n_cat_obs,
                                                double threshold_alpha,
                                                double threshold_beta,
                                                int variable) {
  NumericVector q(no_persons); // Intermediate storage for pseudo-likelihood calculations
  NumericVector r(no_persons); // Intermediate storage for pseudo-likelihood calculations

  // Cache the number of categories and main effect index for the variable
  int n_cats = no_categories(variable, 0); // Number of categories for the variable
  int cat_index = main_index(variable, 0); // Index in `main_effects` for this variable

  NumericVector GroupThresholds(n_cats); // Vector to store group-specific thresholds

  // Iterate over each category threshold for the variable
  for(int category = 0; category < n_cats; category++) {
    double current_state = main_effects(cat_index + category, 0); // Current state of the threshold
    double exp_current = std::exp(current_state); // Exponentiated current threshold
    double c = (threshold_alpha + threshold_beta) / (1 + exp_current); // Initial value for c

    // Compute group-specific thresholds and contributions to `q` and `r`
    for (int gr = 0; gr < no_groups; gr++) {
      // Update thresholds for the current group
      for (int cat = 0; cat < n_cats; cat++) {
        double threshold = main_effects(cat_index + cat, 0);
        for (int h = 1; h < no_groups; h++) {
          threshold += projection(gr, h - 1) * main_effects(cat_index + cat, h);
        }
        GroupThresholds[cat] = threshold;
      }

      // Subtract the current category's base threshold
      GroupThresholds[category] -= main_effects(cat_index + category, 0);

      // Compute `q` and `r` for each person in the group
      for (int person = group_index(gr, 0); person <= group_index(gr, 1); person++) {
        double rest_score = rest_matrix(person, variable);
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
    double tmp = no_persons + threshold_alpha + threshold_beta - exp_current * c;
    c /= tmp;

    // Generate a proposed state using a generalized beta-prime proposal
    double a = threshold_alpha;
    double b = no_persons + threshold_beta;
    for(int gr = 0; gr < no_groups; gr++) {
      IntegerMatrix n_cat_obs_gr = n_cat_obs[gr];
      a += n_cat_obs_gr(category + 1, variable); // Update a with observations in the category
      b -= n_cat_obs_gr(category + 1, variable); // Update b with remaining observations
    }

    // Sample from the beta distribution
    double tmp_beta = R::rbeta(a, b);
    double proposed_state = std::log(tmp_beta / (1  - tmp_beta) / c);
    double exp_proposed = std::exp(proposed_state);

    // Compute the log acceptance probability
    double log_prob = 0.0;

    // Compute pseudo-likelihood ratio
    for (int gr = 0; gr < no_groups; ++gr) {
      for (int person = group_index(gr, 0); person <= group_index(gr, 1); ++person) {
        log_prob += std::log(q[person] + r[person] * exp_current);
        log_prob -= std::log(q[person] + r[person] * exp_proposed);
      }
    }

    // Add prior density ratio
    log_prob -= (threshold_alpha + threshold_beta) * std::log(1 + exp_proposed);
    log_prob += (threshold_alpha + threshold_beta) * std::log(1 + exp_current);

    // Add proposal density ratio
    log_prob -= (a + b) * std::log(1 + c * exp_current);
    log_prob += (a + b) * std::log(1 + c * exp_proposed);

    // Perform Metropolis-Hastings acceptance step
    double U = std::log(R::unif_rand());
    if(U < log_prob) {
      main_effects(cat_index + category, 0) = proposed_state; // Accept the proposal
    }
  }
}


/**
 * Function: compare_anova_log_pseudolikelihood_ratio_main_difference
 * Purpose: Computes the log pseudo-likelihood ratio for a proposed change in
 *          a category threshold difference for an independent samples design.
 *
 * Inputs:
 *  - main_effects: Numeric matrix containing the main effect parameters for all variables and groups.
 *  - main_index: Integer matrix mapping variables to their main effect parameter indices.
 *  - projection: Numeric matrix specifying group-specific scaling for thresholds.
 *  - observations: Integer matrix containing the observed data (individuals x variables).
 *  - no_groups: Total number of groups in the analysis.
 *  - group_index: Integer matrix specifying start and end indices for each group in `observations`.
 *  - no_categories: Integer matrix specifying the number of categories for each variable and group.
 *  - no_persons: Total number of individuals in the dataset.
 *  - rest_matrix: Numeric matrix storing residual effects for pseudo-likelihood calculations.
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
double compare_anova_log_pseudolikelihood_ratio_main_difference(NumericMatrix main_effects,
                                                                IntegerMatrix main_index,
                                                                NumericMatrix projection,
                                                                IntegerMatrix observations,
                                                                int no_groups,
                                                                IntegerMatrix group_index,
                                                                IntegerMatrix no_categories,
                                                                int no_persons,
                                                                NumericMatrix rest_matrix,
                                                                List n_cat_obs,
                                                                int variable,
                                                                int category,
                                                                int h,
                                                                double proposed_state,
                                                                double current_state) {
  double pseudolikelihood_ratio = 0.0; // Initialize the log pseudo-likelihood ratio
  double delta_state = proposed_state - current_state; // Difference between proposed and current states
  int cat_index = main_index(variable, 0); // Index for the main effects of the variable


  // Loop over all groups
  for(int gr = 0; gr < no_groups; gr++) {
    int n_cats = no_categories(variable, gr); // Number of categories for the variable in this group
    NumericVector current_thresholds(n_cats); // Store current thresholds for the group
    NumericVector proposed_thresholds(n_cats); // Store proposed thresholds for the group
    double P = projection(gr, h); // Group-specific projection scaling factor

    // Compute current and proposed thresholds for all categories
    for(int cat = 0; cat < n_cats; cat++) {
      int full_cat_index = cat_index + cat;
      double threshold = main_effects(full_cat_index, 0);
      for (int hh = 1; hh < no_groups; ++hh) {
        threshold += projection(gr, hh - 1) * main_effects(full_cat_index, hh);
      }
      current_thresholds[cat] = threshold;
      proposed_thresholds[cat] = threshold;
    }

    // Adjust the threshold for the specific category
    proposed_thresholds[category] -= P * current_state;
    proposed_thresholds[category] += P * proposed_state;

    // Add the contribution from delta_state based on observations
    IntegerMatrix n_cat_obs_gr = n_cat_obs[gr];
    pseudolikelihood_ratio += delta_state * P * n_cat_obs_gr(category + 1, variable);

    // Loop over all persons in the group
    for(int person = group_index(gr, 0); person <= group_index(gr, 1); person++) {
      double rest_score = rest_matrix(person, variable); // Compute residual score
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
 * Function: compare_anova_metropolis_main_difference_regular
 * Purpose: Implements the Metropolis-Hastings (MH) algorithm to sample from
 *          the full conditional of the category threshold difference parameter
 *          for regular ordinal variables in an independent samples design.
 *
 * Inputs:
 *  - main_effects: Numeric matrix containing the main effect parameters for all variables and groups.
 *  - main_index: Integer matrix mapping variables to their main effect parameter indices.
 *  - projection: Numeric matrix specifying group-specific scaling for thresholds.
 *  - observations: Integer matrix containing the observed data (individuals x variables).
 *  - no_groups: Total number of groups in the analysis.
 *  - group_index: Integer matrix specifying start and end indices for each group in `observations`.
 *  - no_categories: Integer matrix specifying the number of categories for each variable and group.
 *  - no_persons: Total number of individuals in the dataset.
 *  - rest_matrix: Numeric matrix storing residual effects for pseudo-likelihood calculations.
 *  - n_cat_obs: List of integer matrices tracking the number of observations for each category
 *               of a variable in a specific group.
 *  - variable: Index of the variable being updated.
 *  - indicator: Integer matrix indicating active variables for the analysis.
 *  - proposal_sd_main: Numeric matrix specifying proposal standard deviations for category thresholds.
 *  - main_difference_scale: Scale parameter for the Cauchy prior on threshold differences.
 *  - phi: Robbins-Monro learning rate for adaptive proposal standard deviations.
 *  - target_acceptance_probability: Target log acceptance rate for Metropolis-Hastings updates.
 *  - t: Current iteration number of the sampler.
 *  - epsilon_lo: Minimum allowable standard deviation for proposals.
 *  - epsilon_hi: Maximum allowable standard deviation for proposals.
 *
 * Outputs:
 *  - Updates `main_effects` with sampled values for threshold differences.
 *  - Updates `proposal_sd_main` with adjusted proposal standard deviations using Robbins-Monro updates.
 */
void compare_anova_metropolis_main_difference_regular(NumericMatrix main_effects,
                                                      IntegerMatrix main_index,
                                                      NumericMatrix projection,
                                                      IntegerMatrix observations,
                                                      int no_groups,
                                                      IntegerMatrix group_index,
                                                      IntegerMatrix no_categories,
                                                      int no_persons,
                                                      NumericMatrix rest_matrix,
                                                      List n_cat_obs,
                                                      int variable,
                                                      IntegerMatrix indicator,
                                                      NumericMatrix proposal_sd_main,
                                                      double main_difference_scale,
                                                      double phi,
                                                      double target_acceptance_probability,
                                                      int t,
                                                      double epsilon_lo,
                                                      double epsilon_hi) {
  // Precompute Robbins-Monro term
  double exp_neg_log_t_phi = std::exp(-std::log(t) * phi);

  // Check if the variable is active
  if (indicator(variable, variable) != 1) {
    return; // Skip if variable is inactive
  }

  // Look up and cache the base category index for this variable
  int base_cat_index = main_index(variable, 0);

  // Loop over all categories for this variable
  int n_cats = no_categories(variable, 0);
  for (int category = 0; category < n_cats; category++) {
    int cat_index = base_cat_index + category;

    // Loop over groups (starting from h = 1)
    for (int h = 1; h < no_groups; h++) {
      double current_state = main_effects(cat_index, h); // Current threshold difference
      double proposed_state = R::rnorm(current_state, proposal_sd_main(cat_index, h)); // Propose a new state

      // Compute log pseudo-likelihood ratio for proposed vs current state
      double log_prob = compare_anova_log_pseudolikelihood_ratio_main_difference(
        main_effects, main_index, projection, observations, no_groups,
        group_index, no_categories, no_persons, rest_matrix, n_cat_obs,
        variable, category, h - 1, proposed_state, current_state);

      // Add contributions from the Cauchy prior
      log_prob += R::dcauchy(proposed_state, 0.0, main_difference_scale, true);
      log_prob -= R::dcauchy(current_state, 0.0, main_difference_scale, true);

      // Metropolis-Hastings acceptance step
      double U = R::unif_rand();
      if (std::log(U) < log_prob) {
        main_effects(cat_index, h) = proposed_state; // Accept the proposed state
      }

      // Robbins-Monro update to the proposal standard deviation
      proposal_sd_main(cat_index, h) = robbins_monro_update(
        proposal_sd_main(cat_index, h), log_prob, target_acceptance_probability,
        epsilon_lo, epsilon_hi, exp_neg_log_t_phi);
    }
  }
}


/**
 * Function: compare_anova_log_pseudolikelihood_ratio_thresholds_blumecapel
 * Purpose: Computes the log pseudo-likelihood ratio for proposed versus current
 *          category threshold parameters in the Blume-Capel model for an independent samples design.
 *
 * Inputs:
 *  - linear_current: Current value of the linear coefficient for the threshold.
 *  - quadratic_current: Current value of the quadratic coefficient for the threshold.
 *  - linear_proposed: Proposed value of the linear coefficient for the threshold.
 *  - quadratic_proposed: Proposed value of the quadratic coefficient for the threshold.
 *  - variable: Index of the variable for which thresholds are being updated.
 *  - reference_category: Integer vector specifying the reference category for each variable.
 *  - main_effects: Numeric matrix of main effect parameters for all variables and groups.
 *  - main_index: Integer matrix mapping variables to their main effect parameter indices.
 *  - projection: Numeric matrix specifying group-specific scaling for thresholds.
 *  - sufficient_blume_capel: List of integer matrices containing sufficient statistics
 *                            for the Blume-Capel model, including linear and quadratic terms.
 *  - no_persons: Total number of individuals in the dataset.
 *  - no_groups: Total number of groups in the analysis.
 *  - group_index: Integer matrix specifying start and end indices for each group in `observations`.
 *  - rest_matrix: Numeric matrix storing residual effects for pseudo-likelihood calculations.
 *  - no_categories: Integer matrix specifying the number of categories for each variable and group.
 *
 * Outputs:
 *  - Returns the computed log pseudo-likelihood ratio for the proposed versus current parameters.
 */
double compare_anova_log_pseudolikelihood_ratio_thresholds_blumecapel(double linear_current,
                                                                      double quadratic_current,
                                                                      double linear_proposed,
                                                                      double quadratic_proposed,
                                                                      int variable,
                                                                      IntegerVector reference_category,
                                                                      NumericMatrix main_effects,
                                                                      IntegerMatrix main_index,
                                                                      NumericMatrix projection,
                                                                      List sufficient_blume_capel,
                                                                      int no_persons,
                                                                      int no_groups,
                                                                      IntegerMatrix group_index,
                                                                      NumericMatrix rest_matrix,
                                                                      IntegerMatrix no_categories) {
  // Variables for bounds, scores, and likelihood components
  double lbound, bound;
  double rest_score, numerator, denominator, exponent;
  double pseudolikelihood_ratio = 0.0;
  int linear_score, quadratic_score;
  int cat_index = main_index(variable, 0);

  // Loop over all groups
  for(int gr = 0; gr < no_groups; gr++) {
    // Precompute constant terms for the numerator and denominator
    NumericVector constant_numerator (no_categories(variable, gr) + 1);
    NumericVector constant_denominator (no_categories(variable, gr) + 1);

    for(int category = 0; category <= no_categories(variable, gr); category++) {
      linear_score = category;
      quadratic_score = (category - reference_category[variable]) *
        (category - reference_category[variable]);

      // Linear and quadratic contributions for current and proposed states
      constant_numerator[category] = linear_current * linear_score;
      constant_numerator[category] += quadratic_current * quadratic_score;
      constant_denominator[category] = linear_proposed * linear_score ;
      constant_denominator[category] += quadratic_proposed * quadratic_score;

      // Add group-specific contributions
      for(int h = 1; h < no_groups; h++) {
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

    IntegerMatrix sufficient_blume_capel_gr = sufficient_blume_capel[gr];

    // Add contributions from sufficient statistics
    pseudolikelihood_ratio += sufficient_blume_capel_gr(0, variable) * linear_proposed;
    pseudolikelihood_ratio += sufficient_blume_capel_gr(1, variable) * quadratic_proposed;
    pseudolikelihood_ratio -= sufficient_blume_capel_gr(0, variable) * linear_current;
    pseudolikelihood_ratio -= sufficient_blume_capel_gr(1, variable) * quadratic_current;

    // Loop over individuals in the group
    for(int person = group_index(gr, 0); person <= group_index(gr, 1); person++) {
      rest_score = rest_matrix(person, variable);
      if(rest_score > 0) {
        bound = no_categories(variable, gr) * rest_score + lbound;
      } else {
        bound = lbound;
      }

      // Compute the likelihood contributions
      numerator = 0.0;
      denominator = 0.0;
      for(int category = 0; category <= no_categories[variable]; category ++) {
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
 * Function: compare_anova_metropolis_threshold_blumecapel
 * Purpose: Samples from the full conditional distribution of the overall category
 *          threshold parameters for Blume-Capel ordinal variables using the
 *          Metropolis-Hastings (MH) algorithm with Robbins-Monro adaptive updates.
 *
 * Inputs:
 *  - main_effects: Numeric matrix of main effect parameters for all variables and groups.
 *  - main_index: Integer matrix mapping variables to their main effect parameter indices.
 *  - projection: Numeric matrix specifying group-specific scaling for thresholds.
 *  - no_categories: Integer matrix specifying the number of categories for each variable and group.
 *  - sufficient_blume_capel: List of sufficient statistics for the Blume-Capel model,
 *                            including linear and quadratic terms.
 *  - no_persons: Total number of individuals in the dataset.
 *  - no_groups: Total number of groups in the analysis.
 *  - group_index: Integer matrix specifying start and end indices for each group in `observations`.
 *  - variable: Index of the variable for which thresholds are being updated.
 *  - reference_category: Integer vector specifying the reference category for each variable.
 *  - threshold_alpha: Hyperparameter for the Blume-Capel threshold prior distribution.
 *  - threshold_beta: Hyperparameter for the Blume-Capel threshold prior distribution.
 *  - rest_matrix: Numeric matrix storing residual effects for pseudo-likelihood calculations.
 *  - proposal_sd_main: Numeric matrix storing the proposal standard deviations for the parameters.
 *  - phi: Robbins-Monro learning rate.
 *  - target_acceptance_probability: Target log acceptance probability for the MH algorithm.
 *  - t: Current iteration number of the sampler.
 *  - epsilon_lo: Lower bound for the proposal standard deviation.
 *  - epsilon_hi: Upper bound for the proposal standard deviation.
 *
 * Outputs:
 *  - Updates `main_effects` for the linear and quadratic threshold parameters.
 *  - Updates `proposal_sd_main` with Robbins-Monro adjustments for the proposal standard deviations.
 */
void compare_anova_metropolis_threshold_blumecapel(NumericMatrix main_effects,
                                                   IntegerMatrix main_index,
                                                   NumericMatrix projection,
                                                   IntegerMatrix no_categories,
                                                   List sufficient_blume_capel,
                                                   int no_persons,
                                                   int no_groups,
                                                   IntegerMatrix group_index,
                                                   int variable,
                                                   IntegerVector reference_category,
                                                   double threshold_alpha,
                                                   double threshold_beta,
                                                   NumericMatrix rest_matrix,
                                                   NumericMatrix proposal_sd_main,
                                                   double phi,
                                                   double target_acceptance_probability,
                                                   int t,
                                                   double epsilon_lo,
                                                   double epsilon_hi) {
  // Precompute Robbins-Monro adjustment term
  double exp_neg_log_t_phi = std::exp(-std::log(t) * phi);

  // Adaptive Metropolis procedure for the linear Blume-Capel model
  int cat_index = main_index(variable, 0); // Base index for the variable
  double current_state = main_effects(cat_index, 0);
  double proposed_state = R::rnorm(current_state, proposal_sd_main(cat_index, 0));

  // Compute log acceptance probability
  double linear_current = current_state;
  double quadratic_current = main_effects(cat_index + 1, 0);
  double linear_proposed = proposed_state;
  double quadratic_proposed = main_effects(cat_index + 1, 0);

  double log_prob = compare_anova_log_pseudolikelihood_ratio_thresholds_blumecapel(
    linear_current, quadratic_current, linear_proposed, quadratic_proposed,
    variable, reference_category, main_effects, main_index, projection,
    sufficient_blume_capel, no_persons, no_groups, group_index, rest_matrix,
    no_categories);

  // Add prior ratio to log acceptance probability
  log_prob += threshold_alpha * (proposed_state - current_state);
  log_prob += (threshold_alpha + threshold_beta) * std::log(1 + std::exp(current_state));
  log_prob -= (threshold_alpha + threshold_beta) * std::log(1 + std::exp(proposed_state));

  // Perform Metropolis-Hastings acceptance step
  double U = R::unif_rand();
  if(std::log(U) < log_prob) {
    main_effects(cat_index, 0) = proposed_state;
  }

  // Robbins-Monro update for proposal standard deviation
  proposal_sd_main(cat_index, 0) = robbins_monro_update(
    proposal_sd_main(cat_index, 0), log_prob, target_acceptance_probability, epsilon_lo,
    epsilon_hi, exp_neg_log_t_phi);

  // Adaptive Metropolis procedure for the quadratic Blume-Capel model
  current_state = main_effects(cat_index + 1, 0);
  proposed_state = R::rnorm(current_state, proposal_sd_main(cat_index + 1, 0));

  // Compute log acceptance probability
  linear_current = main_effects(cat_index, 0);
  quadratic_current = current_state;
  linear_proposed =  main_effects(cat_index, 0);
  quadratic_proposed =  proposed_state;

  log_prob = compare_anova_log_pseudolikelihood_ratio_thresholds_blumecapel(
    linear_current, quadratic_current, linear_proposed, quadratic_proposed,
    variable, reference_category, main_effects, main_index, projection,
    sufficient_blume_capel, no_persons, no_groups, group_index, rest_matrix,
    no_categories);

  // Add prior ratio to log acceptance probability
  log_prob += threshold_alpha * (proposed_state - current_state);
  log_prob += (threshold_alpha + threshold_beta) * std::log(1 + std::exp(current_state));
  log_prob -= (threshold_alpha + threshold_beta) * std::log(1 + std::exp(proposed_state));

  // Perform Metropolis-Hastings acceptance step
  U = R::unif_rand();
  if(std::log(U) < log_prob) {
    main_effects(cat_index + 1, 0) = proposed_state;
  }

  // Robbins-Monro update for proposal standard deviation
  proposal_sd_main(cat_index + 1, 0) = robbins_monro_update(
    proposal_sd_main(cat_index + 1, 0), log_prob, target_acceptance_probability, epsilon_lo,
    epsilon_hi, exp_neg_log_t_phi);
}


/**
 * Function: compare_anova_log_pseudolikelihood_ratio_main_difference_blumecapel
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
 *  - reference_category: Integer vector specifying the reference category for each variable.
 *  - main_effects: Numeric matrix of main effect parameters for all variables and groups.
 *  - main_index: Integer matrix mapping variables to their main effect parameter indices.
 *  - projection: Numeric matrix specifying group-specific scaling for thresholds.
 *  - sufficient_blume_capel: List of sufficient statistics for the Blume-Capel model,
 *                            including linear and quadratic terms.
 *  - no_persons: Total number of individuals in the dataset.
 *  - no_groups: Total number of groups in the analysis.
 *  - group_index: Integer matrix specifying start and end indices for each group in `observations`.
 *  - rest_matrix: Numeric matrix storing residual effects for pseudo-likelihood calculations.
 *  - no_categories: Integer matrix specifying the number of categories for each variable and group.
 *
 * Outputs:
 *  - Returns the computed log pseudo-likelihood ratio.
 */
double compare_anova_log_pseudolikelihood_ratio_main_difference_blumecapel(double linear_current,
                                                                           double quadratic_current,
                                                                           double linear_proposed,
                                                                           double quadratic_proposed,
                                                                           int variable,
                                                                           int h,
                                                                           IntegerVector reference_category,
                                                                           NumericMatrix main_effects,
                                                                           IntegerMatrix main_index,
                                                                           NumericMatrix projection,
                                                                           List sufficient_blume_capel,
                                                                           int no_persons,
                                                                           int no_groups,
                                                                           IntegerMatrix group_index,
                                                                           NumericMatrix rest_matrix,
                                                                           IntegerMatrix no_categories) {
  double lbound, bound; // Variables for numerical bounds to ensure stability
  double rest_score, numerator, denominator, exponent; // Variables for pseudo-likelihood computations
  double pseudolikelihood_ratio = 0.0; // Accumulated log pseudo-likelihood ratio
  int linear_score, quadratic_score; // Scores for linear and quadratic terms
  int cat_index = main_index(variable, 0); // Base index for variable categories

  // Loop over all groups
  for(int gr = 0; gr < no_groups; gr++) {
    // Precomputed constants for numerator and denominator for each category
    NumericVector constant_numerator (no_categories(variable, gr) + 1);
    NumericVector constant_denominator (no_categories(variable, gr) + 1);
    double P = projection(gr, h); // Scaling factor for group-specific adjustments

    // Pre-compute terms for all categories in the current group
    for(int category = 0; category <= no_categories(variable, gr); category++) {
      // Compute linear and quadratic scores for the current category
      linear_score = category;
      quadratic_score = (category - reference_category[variable]) *
        (category - reference_category[variable]);

      // Initialize numerator and denominator contributions with main effects
      constant_numerator[category] = main_effects(cat_index, 0) * linear_score;
      constant_numerator[category] +=  main_effects(cat_index + 1, 0) * quadratic_score;
      constant_denominator[category] =  main_effects(cat_index, 0) * linear_score ;
      constant_denominator[category] += main_effects(cat_index + 1, 0) * quadratic_score;

      // Add group-specific contributions from projections
      for(int hh = 1; hh < no_groups; hh++) {
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
    IntegerMatrix sufficient_blume_capel_gr = sufficient_blume_capel[gr];
    double sp0 = sufficient_blume_capel_gr(0, variable) * P;
    double sp1 = sufficient_blume_capel_gr(1, variable) * P;
    pseudolikelihood_ratio += sp0 * linear_proposed;
    pseudolikelihood_ratio += sp1 * quadratic_proposed;
    pseudolikelihood_ratio -= sp0 * linear_current;
    pseudolikelihood_ratio -= sp1 * quadratic_current;

    // Process each person in the group to compute pseudo-likelihood terms
    for(int person = group_index(gr, 0); person <= group_index(gr, 1); person++) {
      rest_score = rest_matrix(person, variable);
      if(rest_score > 0) {
        bound = no_categories(variable, gr) * rest_score + lbound;
      } else {
        bound = lbound;
      }

      // Initialize numerator and denominator
      numerator = 0.0;
      denominator = 0.0;

      // Compute category-specific contributions
      for(int category = 0; category <= no_categories[variable]; category ++) {
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
 * Function: compare_anova_metropolis_main_difference_blumecapel
 * Purpose:
 *   Perform Metropolis-Hastings sampling for the full-conditional distribution
 *   of the group-specific category threshold difference parameters for a Blume-Capel
 *   ordinal variable in an ANOVA model.
 *
 * Inputs:
 *   - main_effects: A matrix of main effects for all variables and groups.
 *   - main_index: An integer matrix mapping variables to their category indices.
 *   - projection: A projection matrix mapping group-specific effects.
 *   - no_categories: An integer matrix containing the number of categories for
 *                    each variable across groups.
 *   - sufficient_blume_capel: A list of sufficient statistics for Blume-Capel
 *                             variables per group.
 *   - no_persons: Total number of individuals in the dataset.
 *   - no_groups: Total number of groups in the analysis.
 *   - group_index: A matrix containing the start and end indices for individuals
 *                  in each group.
 *   - variable: Index of the variable being updated.
 *   - reference_category: A vector indicating the reference category for each variable.
 *   - threshold_alpha: Shape parameter of the prior distribution for thresholds.
 *   - threshold_beta: Rate parameter of the prior distribution for thresholds.
 *   - rest_matrix: A matrix of residual scores for the pseudo-likelihood calculation.
 *   - indicator: A binary matrix indicating whether a variable is active for sampling.
 *   - proposal_sd_main: A matrix of proposal standard deviations for the main effects.
 *   - phi: Robbins-Monro learning rate parameter for adaptive Metropolis-Hastings.
 *   - target_acceptance_probability: Target acceptance probability for Robbins-Monro.
 *   - t: Current iteration number in the sampling procedure.
 *   - epsilon_lo: Lower bound for the proposal standard deviation.
 *   - epsilon_hi: Upper bound for the proposal standard deviation.
 *
 * Outputs:
 *   - Updates `main_effects` to reflect accepted proposals for the linear and
 *     quadratic threshold differences.
 *   - Updates `proposal_sd_main` for the linear and quadratic parameters using
 *     Robbins-Monro updates.
 */
void compare_anova_metropolis_main_difference_blumecapel(NumericMatrix main_effects,
                                                         IntegerMatrix main_index,
                                                         NumericMatrix projection,
                                                         IntegerMatrix no_categories,
                                                         List sufficient_blume_capel,
                                                         int no_persons,
                                                         int no_groups,
                                                         IntegerMatrix group_index,
                                                         int variable,
                                                         IntegerVector reference_category,
                                                         double threshold_alpha,
                                                         double threshold_beta,
                                                         NumericMatrix rest_matrix,
                                                         IntegerMatrix indicator,
                                                         NumericMatrix proposal_sd_main,
                                                         double phi,
                                                         double target_acceptance_probability,
                                                         int t,
                                                         double epsilon_lo,
                                                         double epsilon_hi) {
  double log_prob, U; // Log probability and random uniform value for MH step
  double current_state, proposed_state; // Current and proposed parameter values
  double exp_neg_log_t_phi = std::exp(-std::log(t) * phi); // Precomputed Robbins-Monro term

  // Check if the variable is active for sampling
  if(indicator(variable, variable) == 1) {

    // Loop over the group-specific difference effects
    for(int h = 1; h < no_groups; h++) {

      // Adaptive Metropolis procedure for the linear Blume-Capel parameter
      int cat_index = main_index(variable, 0); // Base index for the variable's categories
      current_state = main_effects(cat_index, h); // Current value for linear parameter
      proposed_state = R::rnorm(current_state, proposal_sd_main(cat_index, h)); // Propose new value

      // Compute log pseudo-likelihood ratio for the proposed and current states
      double linear_current = current_state;
      double quadratic_current = main_effects(cat_index + 1, h);
      double linear_proposed = proposed_state;
      double quadratic_proposed = main_effects(cat_index + 1, h);

      log_prob = compare_anova_log_pseudolikelihood_ratio_main_difference_blumecapel(
        linear_current, quadratic_current, linear_proposed, quadratic_proposed,
        variable, h - 1, reference_category, main_effects, main_index,
        projection, sufficient_blume_capel, no_persons, no_groups, group_index,
        rest_matrix, no_categories);

      // Add prior contributions to the log probability
      log_prob += threshold_alpha * (proposed_state - current_state);
      log_prob += (threshold_alpha + threshold_beta) * std::log(1 + std::exp(current_state));
      log_prob -= (threshold_alpha + threshold_beta) * std::log(1 + std::exp(proposed_state));

      // Metropolis-Hastings acceptance step
      U = R::unif_rand();
      if(std::log(U) < log_prob) {
        main_effects(cat_index, h) = proposed_state; // Accept proposed value
      }

      // Robbins-Monro update for the proposal standard deviation
      proposal_sd_main(cat_index, h) = robbins_monro_update(
        proposal_sd_main(cat_index, h), log_prob, target_acceptance_probability, epsilon_lo,
        epsilon_hi, exp_neg_log_t_phi);


      // Adaptive Metropolis procedure for the quadratic Blume-Capel parameter
      current_state = main_effects(cat_index + 1, h); // Current value for quadratic parameter
      proposed_state = R::rnorm(current_state, proposal_sd_main(cat_index + 1, h)); // Propose new value

      // Compute log pseudo-likelihood ratio for the proposed and current states
      linear_current = main_effects(cat_index, h);
      quadratic_current = current_state;
      linear_proposed =  main_effects(cat_index, h);
      quadratic_proposed =  proposed_state;

      log_prob = compare_anova_log_pseudolikelihood_ratio_main_difference_blumecapel(
        linear_current, quadratic_current, linear_proposed, quadratic_proposed,
        variable, h - 1, reference_category, main_effects, main_index,
        projection, sufficient_blume_capel, no_persons, no_groups, group_index,
        rest_matrix, no_categories);

      // Add prior contributions to the log probability
      log_prob += threshold_alpha * (proposed_state - current_state);
      log_prob += (threshold_alpha + threshold_beta) * std::log(1 + std::exp(current_state));
      log_prob -= (threshold_alpha + threshold_beta) * std::log(1 + std::exp(proposed_state));

      // Metropolis-Hastings acceptance step
      U = R::unif_rand();
      if(std::log(U) < log_prob) {
        main_effects(cat_index + 1, h) = proposed_state;
      }

      // Robbins-Monro update for the proposal standard deviation
      proposal_sd_main(cat_index + 1, h) = robbins_monro_update(
        proposal_sd_main(cat_index + 1, h), log_prob, target_acceptance_probability, epsilon_lo,
        epsilon_hi, exp_neg_log_t_phi);
    }
  }
}


/**
 * Function: compare_anova_metropolis_thresholds_regular_free
 * Purpose:
 *   Performs Metropolis-Hastings sampling for the full-conditional distribution
 *   of threshold parameters for a regular binary or ordinal variable in a
 *   free-threshold independent samples design.
 *
 * Inputs:
 *   - main_effects: A matrix of main effects for all variables and groups.
 *   - main_index: An integer matrix mapping variables to their category indices.
 *   - observations: An integer matrix containing observed data for individuals
 *                   by variables.
 *   - no_groups: Total number of groups in the analysis.
 *   - group_index: A matrix containing start and end indices for individuals
 *                  in each group.
 *   - no_categories: A matrix containing the number of categories for each
 *                    variable and group.
 *   - rest_matrix: A matrix of residual scores for the pseudo-likelihood
 *                  calculations.
 *   - n_cat_obs: A list of matrices containing category-specific counts for
 *                each group.
 *   - threshold_alpha: Shape parameter for the prior distribution.
 *   - threshold_beta: Rate parameter for the prior distribution.
 *   - variable: Index of the variable being updated.
 *   - group: Index of the group being updated.
 *
 * Outputs:
 *   - Updates the `main_effects` matrix to reflect accepted proposals for
 *     threshold parameters.
 */
void compare_anova_metropolis_thresholds_regular_free(NumericMatrix main_effects,
                                                      IntegerMatrix main_index,
                                                      IntegerMatrix observations,
                                                      int no_groups,
                                                      IntegerMatrix group_index,
                                                      IntegerMatrix no_categories,
                                                      NumericMatrix rest_matrix,
                                                      List n_cat_obs,
                                                      double threshold_alpha,
                                                      double threshold_beta,
                                                      int variable,
                                                      int group) {

  // Number of persons in the group
  int no_persons = group_index(group, 1) - group_index(group, 0) + 1;

  // Cache group-specific category observations
  IntegerMatrix n_cat_obs_gr = n_cat_obs[group];

  // Base category index for the variable
  int base_cat_index = main_index(variable, 0);
  int n_cats = no_categories(variable, group);

  // Temporary storage for pseudo-likelihood elements
  NumericVector q(no_persons);
  NumericVector r(no_persons);

  // Loop over categories
  for(int category = 0; category < n_cats; category++) {
    double current_state = main_effects(base_cat_index + category, group);
    double exp_current = std::exp(current_state);

    // Initialize scaling factor `c` for the generalized beta-prime proposal
    double c = (threshold_alpha + threshold_beta) / (1 + exp_current);

    // Compute `q` and `r` for each person in the group
    for(int person = 0; person < no_persons; person++) {
      double rest_score = rest_matrix(group_index(group, 0) + person, variable);

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
    c /= (no_persons + threshold_alpha + threshold_beta - exp_current * c);

    // Propose a new state using the generalized beta-prime distribution
    double a = n_cat_obs_gr(category + 1, variable) + threshold_alpha;
    double b = no_persons + threshold_beta - n_cat_obs_gr(category + 1, variable);
    double tmp = R::rbeta(a, b);
    double proposed_state = std::log(tmp / ((1 - tmp) * c));
    double exp_proposed = std::exp(proposed_state);


    // Compute the log acceptance probability
    double log_prob = 0.0;

    // Add pseudo-likelihood ratio contributions
    for (int person = 0; person < no_persons; ++person) {
      log_prob += std::log(q[person] + r[person] * exp_current);
      log_prob -= std::log(q[person] + r[person] * exp_proposed);
    }

    // Add prior ratio contributions
    log_prob -= (threshold_alpha + threshold_beta) * std::log(1 + exp_proposed);
    log_prob += (threshold_alpha + threshold_beta) * std::log(1 + exp_current);

    // Add proposal ratio contributions
    log_prob -= (a + b) * std::log(1 + c * exp_current);
    log_prob += (a + b) * std::log(1 + c * exp_proposed);

    // Perform the Metropolis-Hastings acceptance step
    double U = std::log(R::unif_rand());
    if(U < log_prob) {
      // Update the main effects matrix if the proposal is accepted
      main_effects(base_cat_index + category, group) = proposed_state;
    }
  }
}

/**
 * Function: compare_anova_metropolis_thresholds_blumecapel_free
 * Purpose:
 *   Performs Metropolis-Hastings sampling for the full-conditional distribution
 *   of the threshold parameters (linear and quadratic) for a Blume-Capel ordinal
 *   variable in a free-threshold independent samples design.
 *
 * Inputs:
 *   - main_effects: A matrix of main effects for all variables and groups.
 *   - main_index: An integer matrix mapping variables to their category indices.
 *   - observations: An integer matrix containing observed data for individuals
 *                   by variables.
 *   - no_groups: Total number of groups in the analysis.
 *   - group_index: A matrix containing start and end indices for individuals
 *                  in each group.
 *   - reference_category: A vector specifying the reference category for each variable.
 *   - no_categories: A matrix containing the number of categories for each
 *                    variable and group.
 *   - sufficient_blume_capel: A list of matrices containing sufficient statistics
 *                             for each group.
 *   - no_persons: Total number of individuals in the dataset.
 *   - rest_matrix: A matrix of residual scores for pseudo-likelihood calculations.
 *   - n_cat_obs: A list of matrices containing category-specific counts for
 *                each group and variable.
 *   - threshold_alpha: Shape parameter for the prior distribution.
 *   - threshold_beta: Rate parameter for the prior distribution.
 *   - variable: Index of the variable being updated.
 *   - group: Index of the group being updated.
 *   - proposal_sd_main: A matrix of proposal standard deviations for the
 *                       Metropolis-Hastings updates.
 *   - phi: Robbins-Monro learning rate parameter.
 *   - target_acceptance_probability: Target acceptance rate for
 *                                         Metropolis-Hastings updates.
 *   - t: Current iteration number.
 *   - epsilon_lo: Lower bound for Robbins-Monro updates to proposal standard deviations.
 *   - epsilon_hi: Upper bound for Robbins-Monro updates to proposal standard deviations.
 *
 * Outputs:
 *   - Updates the `main_effects` matrix to reflect accepted proposals for
 *     the linear and quadratic threshold parameters.
 *   - Adjusts the `proposal_sd_main` matrix using Robbins-Monro adaptive updates.
 */
void compare_anova_metropolis_thresholds_blumecapel_free(NumericMatrix main_effects,
                                                         IntegerMatrix main_index,
                                                         IntegerMatrix observations,
                                                         int no_groups,
                                                         IntegerMatrix group_index,
                                                         IntegerVector reference_category,
                                                         IntegerMatrix no_categories,
                                                         List sufficient_blume_capel,
                                                         int no_persons,
                                                         NumericMatrix rest_matrix,
                                                         List n_cat_obs,
                                                         double threshold_alpha,
                                                         double threshold_beta,
                                                         int variable,
                                                         int group,
                                                         NumericMatrix proposal_sd_main,
                                                         double phi,
                                                         double target_acceptance_probability,
                                                         int t,
                                                         double epsilon_lo,
                                                         double epsilon_hi) {

  // Robbins-Monro term precomputed for efficiency
  double exp_neg_log_t_phi = std::exp(-std::log(t) * phi);

  // Retrieve sufficient statistics for the current group
  NumericMatrix sufficient_blume_capel_group = sufficient_blume_capel[group];

  // Preallocate vectors for constants used in likelihood calculations
  NumericVector constant_numerator(no_categories(variable, group) + 1);
  NumericVector constant_denominator(no_categories(variable, group) + 1);

  // Adaptive Metropolis procedure for the linear Blume-Capel parameter
  int cat_index = main_index(variable, 0);
  double current_state = main_effects(cat_index, group);
  double proposed_state = R::rnorm(current_state, proposal_sd_main(cat_index, group));

  // Difference between proposed and current state
  double difference = proposed_state - current_state;

  // Precompute constants for likelihood ratio calculation
  for(int category = 0; category <= no_categories(variable, group); category ++) {
    double exponent = main_effects(cat_index + 1, group) *
                      (category - reference_category[variable]) *
                      (category - reference_category[variable]);
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
  double log_prob = threshold_alpha * difference;
  log_prob += sufficient_blume_capel_group(0, variable) * difference;

  // Loop over individuals in the group to compute likelihood ratio
  for(int person = group_index(group, 0); person <= group_index(group, 1); person++) {
    double rest_score = rest_matrix(person, variable);
    double bound = lbound;
    if(rest_score > 0) {
      bound += no_categories(variable, group) * rest_score;
    }

    // Compute likelihood numerator and denominator
    double numerator = std::exp(constant_numerator[0] - bound);
    double denominator = std::exp(constant_denominator[0] - bound);
    for(int score = 1; score <= no_categories(variable, group); score++) {
      double exponent = score * rest_score - bound;
      numerator += std::exp(constant_numerator[score] + exponent);
      denominator += std::exp(constant_denominator[score] + exponent);
    }
    log_prob += std::log(numerator);
    log_prob -= std::log(denominator);
  }

  // Add prior ratio to log acceptance probability
  log_prob += (threshold_alpha + threshold_beta) * std::log(1 + std::exp(current_state));
  log_prob -= (threshold_alpha + threshold_beta) * std::log(1 + std::exp(proposed_state));

  // Metropolis acceptance step
  double U = R::unif_rand();
  if(std::log(U) < log_prob) {
    main_effects(cat_index, group) = proposed_state;
  }

  // Robbins-Monro adaptive update for proposal standard deviation
  proposal_sd_main(cat_index, group) = robbins_monro_update(
    proposal_sd_main(cat_index, group), log_prob, target_acceptance_probability, epsilon_lo,
    epsilon_hi, exp_neg_log_t_phi);


  // Adaptive Metropolis procedure for the quadratic Blume-Capel parameter
  current_state = main_effects(cat_index + 1, group);
  proposed_state = R::rnorm(current_state, proposal_sd_main(cat_index + 1, group));
  difference = proposed_state - current_state;

  // Recompute constants for quadratic term
  for(int category = 0; category <= no_categories(variable, group); category ++) {
    double exponent = main_effects(cat_index, group) * category;
    int score = (category - reference_category[variable]) *
                (category - reference_category[variable]);

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

  log_prob = threshold_alpha * difference;
  log_prob += sufficient_blume_capel_group(1, variable) * difference;

  // Loop over individuals to compute likelihood ratio for quadratic term
  for(int person = group_index(group, 0); person <= group_index(group, 1); person++) {
    double rest_score = rest_matrix(person, variable);
    double bound = lbound;
    if(rest_score > 0) {
      bound += no_categories[variable] * rest_score;
    }

    double numerator = std::exp(constant_numerator[0] - bound);
    double denominator = std::exp(constant_denominator[0] - bound);

    for(int score = 1; score <= no_categories(variable, group); score ++) {
      double exponent = score * rest_score - bound;
      numerator += std::exp(constant_numerator[score] + exponent);
      denominator += std::exp(constant_denominator[score] + exponent);
    }

    log_prob += std::log(numerator);
    log_prob -= std::log(denominator);
  }
  log_prob += (threshold_alpha + threshold_beta) * std::log(1 + std::exp(current_state));
  log_prob -= (threshold_alpha + threshold_beta) * std::log(1 + std::exp(proposed_state));

  U = R::unif_rand();
  if(std::log(U) < log_prob) {
    main_effects(cat_index + 1, group) = proposed_state;
  }

  // Robbins-Monro update for the proposal standard deviation
  proposal_sd_main(cat_index + 1, group) = robbins_monro_update(
    proposal_sd_main(cat_index + 1, group), log_prob, target_acceptance_probability, epsilon_lo,
    epsilon_hi, exp_neg_log_t_phi);

}


/**
 * Function: compare_anova_gibbs_step_gm
 * Purpose:
 *   Executes a single Gibbs sampling step to update graphical model parameters
 *   for a Bayesian parameter comparison across multiple groups.
 *
 * Inputs:
 *   - main_effects: NumericMatrix containing the main effects across variables and groups.
 *   - main_index: IntegerMatrix mapping variables to their corresponding category indices.
 *   - pairwise_effects: NumericMatrix for pairwise effects between variables.
 *   - pairwise_index: IntegerMatrix mapping variable pairs to their indices in `pairwise_effects`.
 *   - projection: NumericMatrix representing the group-specific scaling projection.
 *   - no_categories: IntegerMatrix containing the number of categories for each variable and group.
 *   - observations: IntegerMatrix with observed data for individuals by variables.
 *   - no_persons: Total number of individuals in the dataset.
 *   - no_groups: Total number of groups.
 *   - group_index: IntegerMatrix specifying start and end indices for individuals in each group.
 *   - n_cat_obs: List of category-specific counts for each group and variable.
 *   - sufficient_blume_capel: List of sufficient statistics for Blume-Capel variables in each group.
 *   - rest_matrix: NumericMatrix of residual scores used in pseudo-likelihood calculations.
 *   - independent_thresholds: Boolean indicating if thresholds are modeled independently for groups.
 *   - ordinal_variable: LogicalVector indicating if variables are ordinal.
 *   - reference_category: IntegerVector specifying reference categories for each variable.
 *   - indicator: IntegerMatrix indicating active pairwise and main differences.
 *   - inclusion_probability_difference: NumericMatrix for inclusion probabilities of pairwise differences.
 *   - proposal_sd_main: NumericMatrix of proposal standard deviations for main effects.
 *   - proposal_sd_pairwise: NumericMatrix of proposal standard deviations for pairwise effects.
 *   - interaction_scale: Scale parameter for pairwise interaction priors.
 *   - main_difference_scale: Scale parameter for main difference priors.
 *   - pairwise_difference_scale: Scale parameter for pairwise difference priors.
 *   - threshold_alpha: Shape parameter for the threshold prior distribution.
 *   - threshold_beta: Rate parameter for the threshold prior distribution.
 *   - phi: Robbins-Monro learning rate parameter.
 *   - target_acceptance_probability: Target acceptance probability for Metropolis-Hastings updates.
 *   - t: Current iteration number.
 *   - epsilon_lo: Lower bound for Robbins-Monro updates.
 *   - epsilon_hi: Upper bound for Robbins-Monro updates.
 *   - difference_selection: Boolean indicating whether to perform difference selection.
 *   - no_pairwise: Total number of pairwise effects.
 *   - no_variables: Total number of variables in the analysis.
 *
 * Outputs:
 *   - Updated model parameters, including:
 *       - `indicator`: Updated inclusion indicators for pairwise differences.
 *       - `main_effects`: Updated main effects matrix.
 *       - `pairwise_effects`: Updated pairwise effects matrix.
 *       - `rest_matrix`: Updated residual scores.
 *       - `proposal_sd_main`: Updated proposal standard deviations for main effects.
 *       - `proposal_sd_pairwise`: Updated proposal standard deviations for pairwise effects.
 *
 * Steps:
 *   1. Update pairwise interaction parameters using Metropolis-Hastings.
 *   2. Update the selection of pairwise differences using Metropolis-Hastings (between model).
 *   3. Update pairwise differences using Metropolis-Hastings (within model).
 *   4. Update thresholds (main effects) for each variable:
 *       - If `independent_thresholds`, update thresholds separately for each group.
 *       - Otherwise, model group differences in thresholds.
 */
List compare_anova_gibbs_step_gm(NumericMatrix main_effects,
                                 IntegerMatrix main_index,
                                 NumericMatrix pairwise_effects,
                                 IntegerMatrix pairwise_index,
                                 NumericMatrix projection,
                                 IntegerMatrix no_categories,
                                 IntegerMatrix observations,
                                 int no_persons,
                                 int no_groups,
                                 IntegerMatrix group_index,
                                 List n_cat_obs,
                                 List sufficient_blume_capel,
                                 NumericMatrix rest_matrix,
                                 bool independent_thresholds,
                                 LogicalVector ordinal_variable,
                                 IntegerVector reference_category,
                                 IntegerMatrix indicator,
                                 NumericMatrix inclusion_probability_difference,
                                 NumericMatrix proposal_sd_main,
                                 NumericMatrix proposal_sd_pairwise,
                                 double interaction_scale,
                                 double main_difference_scale,
                                 double pairwise_difference_scale,
                                 double threshold_alpha,
                                 double threshold_beta,
                                 double phi,
                                 double target_acceptance_probability,
                                 int t,
                                 double epsilon_lo,
                                 double epsilon_hi,
                                 bool difference_selection,
                                 int no_pairwise,
                                 int no_variables,
                                 IntegerMatrix index) {

  NumericMatrix proposal_sd_blumecapel_group(no_variables, 2);

  // Step 1: Update pairwise interaction parameters
  compare_anova_metropolis_interaction(
    main_effects, pairwise_effects, main_index,
    pairwise_index, projection, observations, no_groups, group_index,
    no_categories, independent_thresholds, no_persons, rest_matrix,
    ordinal_variable, reference_category, proposal_sd_pairwise,
    interaction_scale, no_variables, phi, target_acceptance_probability, t, epsilon_lo, epsilon_hi);


  //  Step 2: Update the selection of pairwise differences
  if(difference_selection == true) {
    compare_anova_metropolis_pairwise_difference_between_model(
      inclusion_probability_difference, index, main_effects, pairwise_effects,
      main_index, pairwise_index, projection, observations, no_groups,
      group_index, no_categories, independent_thresholds, indicator, no_persons,
      rest_matrix, ordinal_variable, reference_category, proposal_sd_pairwise,
      pairwise_difference_scale, no_variables, phi,
      target_acceptance_probability, t, epsilon_lo, epsilon_hi, no_pairwise);
  }

  // Step 3: Update pairwise differences
  compare_anova_metropolis_pairwise_difference(
    main_effects, pairwise_effects, main_index, pairwise_index, projection,
    observations, no_groups, group_index, no_categories, independent_thresholds,
    indicator, no_persons, rest_matrix, ordinal_variable, reference_category,
    proposal_sd_pairwise, pairwise_difference_scale, no_variables, phi, target_acceptance_probability, t, epsilon_lo,
    epsilon_hi);

  // Step 4: Update thresholds
  for (int variable = 0; variable < no_variables; variable++) {
    if (independent_thresholds) {
      // Independent thresholds: Update separately for each group
      for (int group = 0; group < no_groups; ++group) {
        if (ordinal_variable[variable]) {
          compare_anova_metropolis_thresholds_regular_free(
            main_effects, main_index, observations, no_groups, group_index,
            no_categories, rest_matrix, n_cat_obs, threshold_alpha,
            threshold_beta, variable, group);
        } else {
          compare_anova_metropolis_thresholds_blumecapel_free(
            main_effects, main_index, observations, no_groups, group_index,
            reference_category, no_categories, sufficient_blume_capel, no_persons,
            rest_matrix, n_cat_obs, threshold_alpha, threshold_beta, variable,
            group, proposal_sd_main, phi, target_acceptance_probability, t,
            epsilon_lo, epsilon_hi);
        }
      }
    } else {
      // Thresholds modeled with group differences
      if (ordinal_variable[variable]) {
        compare_anova_metropolis_threshold_regular(
          main_effects, main_index, projection, observations, no_groups,
          group_index, no_categories, no_persons, rest_matrix, n_cat_obs,
          threshold_alpha, threshold_beta, variable);

        compare_anova_metropolis_main_difference_regular(
          main_effects, main_index, projection, observations, no_groups,
          group_index, no_categories, no_persons, rest_matrix, n_cat_obs,
          variable, indicator, proposal_sd_main, main_difference_scale, phi,
          target_acceptance_probability, t, epsilon_lo, epsilon_hi);
      } else {
        compare_anova_metropolis_threshold_blumecapel(
          main_effects, main_index, projection, no_categories, sufficient_blume_capel,
          no_persons, no_groups, group_index, variable, reference_category,
          threshold_alpha, threshold_beta, rest_matrix, proposal_sd_main, phi,
          target_acceptance_probability, t, epsilon_lo, epsilon_hi);

        compare_anova_metropolis_main_difference_blumecapel(
          main_effects, main_index, projection, no_categories, sufficient_blume_capel,
          no_persons, no_groups, group_index, variable, reference_category,
          threshold_alpha, threshold_beta, rest_matrix, indicator,
          proposal_sd_main, phi, target_acceptance_probability, t,
          epsilon_lo, epsilon_hi);
      }
    }
  }

  // Return updated parameters
  return List::create(Named("indicator") = indicator,
                      Named("main_effects") = main_effects,
                      Named("pairwise_effects") = pairwise_effects,
                      Named("rest_matrix") = rest_matrix,
                      Named("proposal_sd_main") = proposal_sd_main,
                      Named("proposal_sd_pairwise") = proposal_sd_pairwise);
}

/**
 * Function: compare_anova_gibbs_sampler
 * Purpose:
 *   Executes the Gibbs sampling process for Bayesian parameter comparisons across multiple groups.
 *   Updates main effects, pairwise effects, thresholds, and inclusion probabilities iteratively.
 *
 * Inputs:
 *   - observations: IntegerMatrix of observed data (individuals by variables).
 *   - main_index: IntegerMatrix mapping variables to their category indices.
 *   - pairwise_index: IntegerMatrix mapping variable pairs to pairwise effect indices.
 *   - projection: NumericMatrix representing group-specific scaling.
 *   - no_categories: IntegerMatrix containing category counts for each variable and group.
 *   - no_groups: Total number of groups in the analysis.
 *   - group_index: IntegerMatrix specifying group-wise start and end indices for individuals.
 *   - interaction_scale: Scale parameter for pairwise interaction priors.
 *   - pairwise_difference_scale: Scale parameter for pairwise difference priors.
 *   - main_difference_scale: Scale parameter for main difference priors.
 *   - pairwise_difference_prior: Type of prior for pairwise differences (e.g., "Beta-Bernoulli").
 *   - main_difference_prior: Type of prior for main differences (e.g., "Beta-Bernoulli").
 *   - inclusion_probability_difference: NumericMatrix for inclusion probabilities of differences.
 *   - pairwise_beta_bernoulli_alpha: Alpha parameter for Beta-Bernoulli prior on pairwise differences.
 *   - pairwise_beta_bernoulli_beta: Beta parameter for Beta-Bernoulli prior on pairwise differences.
 *   - main_beta_bernoulli_alpha: Alpha parameter for Beta-Bernoulli prior on main differences.
 *   - main_beta_bernoulli_beta: Beta parameter for Beta-Bernoulli prior on main differences.
 *   - Index: IntegerMatrix for randomization during Gibbs sampling iterations.
 *   - iter: Total number of Gibbs sampling iterations.
 *   - burnin: Number of burn-in iterations.
 *   - n_cat_obs: List of category-specific observation counts for each group.
 *   - sufficient_blume_capel: List of sufficient statistics for Blume-Capel variables by group.
 *   - threshold_alpha: Shape parameter for the threshold prior distribution.
 *   - threshold_beta: Rate parameter for the threshold prior distribution.
 *   - na_impute: Boolean flag indicating whether to impute missing data.
 *   - missing_index: IntegerMatrix of indices for missing data entries.
 *   - ordinal_variable: LogicalVector indicating whether variables are ordinal.
 *   - reference_category: IntegerVector specifying reference categories for each variable.
 *   - independent_thresholds: Boolean flag for whether thresholds are modeled independently across groups.
 *   - save: Boolean flag for saving results during Gibbs sampling (default: false).
 *   - display_progress: Boolean flag for displaying a progress bar during sampling (default: false).
 *   - difference_selection: Boolean flag for whether to perform difference selection (default: true).
 *
 * Outputs:
 *   - A List containing:
 *       - `main`: Matrix of updated main effects after sampling.
 *       - `pairwise`: Matrix of updated pairwise effects after sampling.
 *       - `indicator`: Matrix of updated inclusion probabilities after sampling.
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
List compare_anova_gibbs_sampler(IntegerMatrix observations,
                                 IntegerMatrix main_index,
                                 IntegerMatrix pairwise_index,
                                 NumericMatrix projection,
                                 IntegerMatrix no_categories,
                                 int no_groups,
                                 IntegerMatrix group_index,
                                 double interaction_scale,
                                 double pairwise_difference_scale,
                                 double main_difference_scale,
                                 String pairwise_difference_prior,
                                 String main_difference_prior,
                                 NumericMatrix inclusion_probability_difference,
                                 double pairwise_beta_bernoulli_alpha,
                                 double pairwise_beta_bernoulli_beta,
                                 double main_beta_bernoulli_alpha,
                                 double main_beta_bernoulli_beta,
                                 IntegerMatrix Index,
                                 int iter,
                                 int burnin,
                                 List n_cat_obs,
                                 List sufficient_blume_capel,
                                 double threshold_alpha,
                                 double threshold_beta,
                                 bool na_impute,
                                 IntegerMatrix missing_index,
                                 LogicalVector ordinal_variable,
                                 IntegerVector reference_category,
                                 bool independent_thresholds,
                                 bool save = false,
                                 bool display_progress = false,
                                 bool difference_selection = true) {

  // Step 1: Initialize parameters and helper matrices
  int no_variables = observations.ncol();
  int no_persons = observations.nrow();
  int no_main = 1 + main_index(no_variables - 1, 1);
  int no_pairwise = no_variables * (no_variables - 1) / 2;

  // Initialize person-group indicator
  IntegerVector person_group_indicator (no_persons);
  for(int group = 0; group < no_groups; group++) {
    for(int person = group_index(group, 0); person < group_index(group, 1) + 1; person++) {
      person_group_indicator[person] = group;
    }
  }

  // Initialize model parameters
  NumericMatrix main_effects(no_main, no_groups);
  NumericMatrix pairwise_effects(no_pairwise, no_groups);
  IntegerMatrix indicator(no_variables, no_variables);
  std::fill(indicator.begin(), indicator.end(), 1);

  // Adaptive Metropolis proposal standard deviations
  NumericMatrix proposal_sd_main(no_main, no_groups);
  NumericMatrix proposal_sd_pairwise(no_pairwise, no_groups);
  std::fill(proposal_sd_main.begin(), proposal_sd_main.end(), 1.0);
  std::fill(proposal_sd_pairwise.begin(), proposal_sd_pairwise.end(), 1.0);

  // Robbins-Monro parameters
  double phi = 0.75, target_acceptance_probability = 0.234;
  int rm_scale = group_index(0, 1) - group_index(0, 0);
  for(int gr = 1; gr < no_groups; gr++) {
    int tmp = group_index(gr, 1) - group_index(gr, 0);
    if(tmp > rm_scale)
      rm_scale = tmp;
  }

  double epsilon_lo = 1.0 / static_cast<double>(rm_scale);
  double epsilon_hi = 2.0;

  //Randomized index for the pairwise updates ----------------------------------
  IntegerVector v = seq(0, no_pairwise - 1);
  IntegerVector order(no_pairwise);
  IntegerMatrix index(no_pairwise, 3);

  // Rest matrix for pseudo-likelihoods
  NumericMatrix rest_matrix(no_persons, no_variables);

  // Output matrices
  NumericMatrix out_main(no_main, no_groups);
  NumericMatrix out_pairwise(no_pairwise, no_groups);
  NumericMatrix out_indicator(no_variables, no_variables);

  // Progress bar
  Progress p(iter + burnin, display_progress);


  // Store the original difference selection input
  bool enable_difference_selection = difference_selection;

  // Compute the total burn-in duration
  int total_burnin = burnin * (enable_difference_selection ? 2 : 1);

  // Flag to enable difference selection after the initial burn-in phase
  difference_selection = false;

  // Step 2: Gibbs sampling loop
  for (int iteration = 0; iteration < iter + total_burnin; ++iteration) {
    if (Progress::check_abort()) {
      return List::create(Named("main") = out_main,
                          Named("pairwise") = out_pairwise,
                          Named("indicator") = out_indicator);

    }
    p.increment();

    // Enable difference selection at the midpoint of burn-in
    if (enable_difference_selection && iteration == burnin) {
      difference_selection = true;
    }

    // Create a random ordering of pairwise effects for updating
    order = sample(v,
                   no_pairwise,
                   false,
                   R_NilValue);

    for(int cntr = 0; cntr < no_pairwise; cntr++) {
      index(cntr, 0) = Index(order[cntr], 0);
      index(cntr, 1) = Index(order[cntr], 1);
      index(cntr, 2) = Index(order[cntr], 2);
    }

    // Handle missing data if required
    if (na_impute) {
      List impute_out = compare_anova_impute_missing_data(
        main_effects, pairwise_effects, main_index, pairwise_index,
        projection, observations, no_groups, person_group_indicator,
        n_cat_obs, sufficient_blume_capel, no_categories, rest_matrix,
        missing_index, ordinal_variable, reference_category, independent_thresholds);

      IntegerMatrix observations_tmp = impute_out["observations"];
      List n_cat_obs_tmp = impute_out["n_cat_obs"];
      List sufficient_blume_capel_tmp = impute_out["sufficient_blume_capel"];
      NumericMatrix rest_matrix_tmp = impute_out["rest_matrix"];

      // Reassign to original variables
      observations = observations_tmp;
      n_cat_obs = n_cat_obs_tmp;
      sufficient_blume_capel = sufficient_blume_capel_tmp;
      rest_matrix = rest_matrix_tmp;
    }

    // Perform parameter updates
    List step_out = compare_anova_gibbs_step_gm(
      main_effects, main_index, pairwise_effects, pairwise_index, projection,
      no_categories, observations, no_persons, no_groups, group_index, n_cat_obs,
      sufficient_blume_capel, rest_matrix, independent_thresholds, ordinal_variable,
      reference_category, indicator, inclusion_probability_difference,
      proposal_sd_main, proposal_sd_pairwise, interaction_scale,
      main_difference_scale, pairwise_difference_scale, threshold_alpha,
      threshold_beta, phi, target_acceptance_probability, iteration, epsilon_lo, epsilon_hi,
      difference_selection, no_pairwise, no_variables, index);

    IntegerMatrix indicator_tmp = step_out["indicator"];
    NumericMatrix main_effects_tmp = step_out["main_effects"];
    NumericMatrix pairwise_effects_tmp = step_out["pairwise_effects"];
    NumericMatrix rest_matrix_tmp = step_out["rest_matrix"];
    NumericMatrix proposal_sd_pairwise_tmp = step_out["proposal_sd_pairwise"];
    NumericMatrix proposal_sd_main_tmp = step_out["proposal_sd_main"];

    // Update matrices from Gibbs step output
    indicator = indicator_tmp;
    main_effects = main_effects_tmp;
    pairwise_effects = pairwise_effects_tmp;
    rest_matrix = rest_matrix_tmp;
    proposal_sd_pairwise = proposal_sd_pairwise_tmp;
    proposal_sd_main = proposal_sd_main_tmp;


    // Update inclusion probabilities
    if (difference_selection) {
      int sumG = 0;

      if (pairwise_difference_prior == "Beta-Bernoulli") {
        // Update pairwise inclusion probabilities
        for (int i = 0; i < no_variables - 1; ++i) {
          for (int j = i + 1; j < no_variables; ++j) {
            sumG += indicator(i, j);
          }
        }
        double prob = R::rbeta(pairwise_beta_bernoulli_alpha + sumG,
                               pairwise_beta_bernoulli_beta + no_pairwise - sumG);
        std::fill(inclusion_probability_difference.begin(),
                  inclusion_probability_difference.end(), prob);
      }

      if (main_difference_prior == "Beta-Bernoulli") {
        // Update main inclusion probabilities
        sumG = 0;
        for (int i = 0; i < no_variables; ++i) {
          sumG += indicator(i, i);
        }
        double prob = R::rbeta(main_beta_bernoulli_alpha + sumG,
                               main_beta_bernoulli_beta + no_variables - sumG);
        for (int i = 0; i < no_variables; ++i) {
          inclusion_probability_difference(i, i) = prob;
        }
      }
    }

    // Save results after burn-in
    if (iteration >= total_burnin) {
      int iter_adj = iteration - total_burnin + 1;
      for (int col = 0; col < no_groups; ++col) {
        for (int row = 0; row < no_main; ++row) {
          out_main(row, col) = (out_main(row, col) * (iter_adj - 1) + main_effects(row, col)) /
            static_cast<double>(iter_adj);
        }
      }
      for (int col = 0; col < no_groups; ++col) {
        for (int row = 0; row < no_pairwise; ++row) {
          out_pairwise(row, col) = (out_pairwise(row, col) * (iter_adj - 1) + pairwise_effects(row, col)) /
            static_cast<double>(iter_adj);
        }
      }
      for (int i = 0; i < no_variables - 1; ++i) {
        for (int j = i + 1; j < no_variables; ++j) {
          out_indicator(i, j) = (out_indicator(i, j) * (iter_adj - 1) + indicator(i, j)) /
            static_cast<double>(iter_adj);
          out_indicator(j, i) = out_indicator(i, j);
        }
      }
    }
  }

  return List::create(Named("main") = out_main,
                      Named("pairwise") = out_pairwise,
                      Named("indicator") = out_indicator);
}
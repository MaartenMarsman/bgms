// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppProgress)]]
#include <Rcpp.h>
#include "gibbs_functions_edge_prior.h"
#include <progress.hpp>
#include <progress_bar.hpp>
using namespace Rcpp;

/**
 * Function: impute_missing_data
 * Purpose:
 *   Imputes missing observations for a graphical model with ordinal and Blume-Capel variables
 *   using their full conditional pseudo-posterior distributions.
 *
 * Inputs:
 *   - interactions: Pairwise interaction matrix between variables.
 *   - thresholds: Threshold parameters for ordinal or Blume-Capel variables.
 *   - observations: Matrix of observed data (individuals x variables), updated in-place.
 *   - num_obs_categories: Matrix of counts for each category per variable (for ordinal variables).
 *   - sufficient_blume_capel: Sufficient statistics (linear and squared sums) for BC variables.
 *   - num_categories: Number of categories per variable.
 *   - residual_matrix: Matrix of linear predictors excluding the current variable, updated in-place.
 *   - missing_index: Matrix of (row, variable) indices for missing observations.
 *   - is_ordinal_variable: Indicator vector specifying whether a variable is ordinal or Blume-Capel.
 *   - reference_category: Reference category used in Blume-Capel variables.
 *
 * Outputs:
 *   - A List with updated:
 *     - `observations`: Imputed observation matrix.
 *     - `num_obs_categories`: Updated category counts (ordinal only).
 *     - `sufficient_blume_capel`: Updated sufficient statistics for BC variables.
 *     - `residual_matrix`: Updated pseudo-likelihood contributions.
 */
List impute_missing_data(
    const arma::mat& interactions,
    const arma::mat& thresholds,
    arma::imat& observations,
    arma::imat& num_obs_categories,
    arma::imat& sufficient_blume_capel,
    const arma::ivec& num_categories,
    arma::mat& residual_matrix,
    const arma::imat& missing_index,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category
) {

  const int num_variables = observations.n_cols;
  const int num_missings = missing_index.n_rows;
  const int max_num_categories = num_categories.max();

  arma::vec probabilities(max_num_categories + 1);

  // Loop over all missing values
  for(int miss = 0; miss < num_missings; miss++) {
    // Identify the observation to impute
    const int person = missing_index(miss, 0);
    const int variable = missing_index(miss, 1);

    const double rest_score = residual_matrix(person, variable);
    const int num_cats = num_categories[variable];
    const bool is_ordinal = is_ordinal_variable[variable];

    double cumsum = 0.0;

    // Generate a new observation based on the model
    if(is_ordinal) {
      // For regular binary or ordinal variables
      cumsum += 1.0;
      probabilities[0] = 1.0;
      for(int category = 0; category < num_cats; category++) {
        int score = category + 1;
        double exponent = thresholds(variable, category) + score * rest_score;
        cumsum += std::exp(exponent);
        probabilities[score] = cumsum;
      }
    } else {
      // For Blume-Capel variables
      const int ref = reference_category[variable];
      double exponent = thresholds(variable, 1) * ref * ref;
      cumsum = std::exp(exponent);
      cumsum = std::exp(exponent);
      probabilities[0] = cumsum;
      for(int category = 0; category < num_cats; category++) {
        int score = category + 1;
        int centered = score - ref;
        exponent = thresholds(variable, 0) * score +
          thresholds(variable, 1) * centered * centered +
          score * rest_score;
        cumsum += std::exp(exponent);
        probabilities[score] = cumsum;
      }
    }

    // Sample a new value based on computed probabilities
    double u = cumsum * R::unif_rand();
    int sampled_score = 0;
    while (u > probabilities[sampled_score]) sampled_score++;

    int new_obs = sampled_score;
    int old_obs = observations(person, variable);
    if(old_obs != new_obs) {
      // Update raw observations
      observations(person, variable) = new_obs;

      // Update category counts or sufficient statistics
      if(is_ordinal) {
        num_obs_categories(old_obs, variable)--;
        num_obs_categories(new_obs, variable)++;
      } else {
        const int ref = reference_category[variable];
        sufficient_blume_capel(0, variable) += (new_obs - old_obs);
        sufficient_blume_capel(1, variable) +=
          (new_obs - ref) * (new_obs - ref) -
          (old_obs - ref) * (old_obs - ref);
      }

      // Update rest scores
      for(int v = 0; v < num_variables; v++) {
        double delta = (new_obs - old_obs) * interactions(v, variable);
        residual_matrix(person, v) += delta;
      }
    }
  }

  return List::create(
    Named("observations") = observations,
    Named("num_obs_categories") = num_obs_categories,
    Named("sufficient_blume_capel") = sufficient_blume_capel,
    Named("residual_matrix") = residual_matrix);
}

// ----------------------------------------------------------------------------|
// MH algorithm to sample from the full-conditional of the threshold parameters
//   for a regular binary or ordinal variable
// ----------------------------------------------------------------------------|
void metropolis_thresholds_regular(
    arma::mat& thresholds,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& num_obs_categories,
    const int no_persons,
    const int variable,
    const double threshold_alpha,
    const double threshold_beta,
    const arma::mat& residual_matrix
) {

  arma::vec g(no_persons);
  arma::vec q(no_persons);

  double log_prob, rest_score;
  double a, b, c;
  double tmp;
  double current_state, proposed_state;
  double U;
  double exp_current, exp_proposed;

  for(int category = 0; category < num_categories[variable]; category++) {
    current_state = thresholds(variable, category);
    exp_current = std::exp(current_state);
    c = (threshold_alpha + threshold_beta) / (1 + exp_current);
    for(int person = 0; person < no_persons; person++) {
      g[person] = 1.0;
      q[person] = 1.0;
      rest_score = residual_matrix(person, variable);
      for(int cat = 0; cat < num_categories[variable]; cat++) {
        if(cat != category) {
          g[person] += std::exp(thresholds(variable, cat) +
            (cat + 1) * rest_score);
        }
      }
      q[person] = std::exp((category + 1) * rest_score);
      c +=  q[person] / (g[person] + q[person] * exp_current);
    }
    c = c / ((no_persons + threshold_alpha + threshold_beta) -
      exp_current * c);

    //Proposal is generalized beta-prime.
    a = num_obs_categories(category + 1, variable) + threshold_alpha;
    b = no_persons + threshold_beta - num_obs_categories(category + 1, variable);
    tmp = R::rbeta(a, b);
    proposed_state = std::log(tmp / (1  - tmp) / c);
    exp_proposed = exp(proposed_state);

    //Compute log_acceptance probability for Metropolis.
    //First, we use g and q above to compute the ratio of pseudolikelihoods
    log_prob = 0;
    for(int person = 0; person < no_persons; person++) {
      log_prob += std::log(g[person] + q[person] * exp_current);
      log_prob -= std::log(g[person] + q[person] * exp_proposed);
    }
    //Second, we add the ratio of prior probabilities
    log_prob -= (threshold_alpha + threshold_beta) *
      std::log(1 + exp_proposed);
    log_prob += (threshold_alpha + threshold_beta) *
      std::log(1 + exp_current);
    //Third, we add the ratio of proposals
    log_prob -= (a + b) * std::log(1 + c * exp_current);
    log_prob += (a + b) * std::log(1 + c * exp_proposed);

    U = std::log(R::unif_rand());
    if(U < log_prob) {
      thresholds(variable, category) = proposed_state;
    }
  }
}

// ----------------------------------------------------------------------------|
// Adaptive Metropolis algorithm to sample from the full-conditional of the
//   threshold parameters for a Blume-Capel ordinal variable
// ----------------------------------------------------------------------------|
void metropolis_thresholds_blumecapel(
    arma::mat& thresholds,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& sufficient_blume_capel,
    const int no_persons,
    const int variable,
    const arma::ivec& reference_category,
    const double threshold_alpha,
    const double threshold_beta,
    const arma::mat& residual_matrix,
    arma::mat& proposal_sd_blumecapel,
    const double phi,
    const double target_ar,
    const int t,
    const double epsilon_lo,
    const double epsilon_hi
) {

  double log_prob, U;
  double current_state, proposed_state, difference;
  double numerator, denominator;
  double lbound, bound, exponent, rest_score;
  arma::vec constant_numerator (num_categories[variable] + 1);
  arma::vec constant_denominator (num_categories[variable] + 1);

  //----------------------------------------------------------------------------
  //Adaptive Metropolis for the linear Blume-Capel parameter
  //----------------------------------------------------------------------------
  current_state = thresholds(variable, 0);
  proposed_state = R::rnorm(current_state, proposal_sd_blumecapel(variable, 0));

  //Precompute terms for the log acceptance probability ------------------------
  difference = proposed_state - current_state;

  for(int category = 0; category < num_categories[variable] + 1; category ++) {
    exponent = thresholds(variable, 1) *
      (category - reference_category[variable]) *
      (category - reference_category[variable]);
    constant_numerator[category] = current_state * category + exponent;
    constant_denominator[category] = proposed_state * category + exponent;
  }
  double tmp_n = max(constant_numerator);
  double tmp_d = max(constant_denominator);
  if(tmp_n > 0) {
    if(tmp_n > tmp_d) {
      lbound = tmp_n;
    } else {
      lbound = tmp_d;
    }
  } else {
    lbound = 0.0;
  }

  //Compute the log acceptance probability -------------------------------------
  log_prob = threshold_alpha * difference;
  log_prob += sufficient_blume_capel(0, variable) * difference;

  for(int person = 0; person < no_persons; person++) {
    rest_score = residual_matrix(person, variable);
    if(rest_score > 0) {
      bound = num_categories[variable] * rest_score + lbound;
    } else {
      bound = lbound;
    }
    numerator = std::exp(constant_numerator[0] - bound);
    denominator = std::exp(constant_denominator[0] - bound);
    for(int category = 0; category < num_categories[variable]; category ++) {
      exponent = (category + 1) * rest_score - bound;
      numerator += std::exp(constant_numerator[category + 1] + exponent);
      denominator += std::exp(constant_denominator[category + 1] + exponent);
    }
    log_prob += std::log(numerator);
    log_prob -= std::log(denominator);
  }

  log_prob += (threshold_alpha + threshold_beta) *
    std::log(1 + std::exp(current_state));
  log_prob -= (threshold_alpha + threshold_beta) *
    std::log(1 + std::exp(proposed_state));

  U = R::unif_rand();
  if(std::log(U) < log_prob) {
    thresholds(variable, 0) = proposed_state;
  }

  //Robbins-Monro update of the proposal variance ------------------------------
  if(log_prob > 0) {
    log_prob = 1;
  } else {
    log_prob = std::exp(log_prob);
  }

  double update_proposal_sd = proposal_sd_blumecapel(variable, 0) +
    (log_prob - target_ar) * std::exp(-log(t) * phi);

  if(std::isnan(update_proposal_sd) == true) {
    update_proposal_sd = 1.0;
  }

  update_proposal_sd = std::clamp(update_proposal_sd, epsilon_lo, epsilon_hi);
  proposal_sd_blumecapel(variable, 0) = update_proposal_sd;

  //----------------------------------------------------------------------------
  //Adaptive Metropolis for the quadratic Blume-Capel parameter
  //----------------------------------------------------------------------------
  current_state = thresholds(variable, 1);
  proposed_state = R::rnorm(current_state, proposal_sd_blumecapel(variable, 1));

  //Precompute terms for the log acceptance probability ------------------------
  difference = proposed_state - current_state;

  for(int category = 0; category < num_categories[variable] + 1; category ++) {
    exponent = thresholds(variable, 0) * category;
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

  //Compute the log acceptance probability -------------------------------------
  log_prob = threshold_alpha * difference;
  log_prob += sufficient_blume_capel(1, variable) * difference;

  for(int person = 0; person < no_persons; person++) {
    rest_score = residual_matrix(person, variable);
    if(rest_score > 0) {
      bound = num_categories[variable] * rest_score + lbound;
    } else {
      bound = lbound;
    }

    numerator = std::exp(constant_numerator[0] - bound);
    denominator = std::exp(constant_denominator[0] - bound);

    for(int category = 0; category < num_categories[variable]; category ++) {
      exponent = (category + 1) * rest_score - bound;
      numerator += std::exp(constant_numerator[category + 1] + exponent);
      denominator += std::exp(constant_denominator[category + 1] + exponent);
    }

    log_prob += std::log(numerator);
    log_prob -= std::log(denominator);
  }
  log_prob += (threshold_alpha + threshold_beta) *
    std::log(1 + std::exp(current_state));
  log_prob -= (threshold_alpha + threshold_beta) *
    std::log(1 + std::exp(proposed_state));

  U = R::unif_rand();
  if(std::log(U) < log_prob) {
    thresholds(variable, 1) = proposed_state;
  }

  //Robbins-Monro update of the proposal variance ------------------------------
  if(log_prob > 0) {
    log_prob = 1;
  } else {
    log_prob = std::exp(log_prob);
  }

  update_proposal_sd = proposal_sd_blumecapel(variable, 1) +
    (log_prob - target_ar) * std::exp(-log(t) * phi);

  if(std::isnan(update_proposal_sd) == true) {
    update_proposal_sd = 1.0;
  }

  update_proposal_sd = std::clamp(update_proposal_sd, epsilon_lo, epsilon_hi);
  proposal_sd_blumecapel(variable, 1) = update_proposal_sd;

}

// ----------------------------------------------------------------------------|
// The log pseudolikelihood ratio [proposed against current] for an interaction
// ----------------------------------------------------------------------------|
double log_pseudolikelihood_ratio_interaction(
    const arma::mat& interactions,
    const arma::mat& thresholds,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const int no_persons,
    const int variable1,
    const int variable2,
    const double proposed_state,
    const double current_state,
    const arma::mat& residual_matrix,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category
) {
  double rest_score, bound;
  double pseudolikelihood_ratio = 0.0;
  double denominator_prop, denominator_curr, exponent;
  int score, obs_score1, obs_score2;

  double delta_state = proposed_state - current_state;

  for(int person = 0; person < no_persons; person++) {
    obs_score1 = observations(person, variable1);
    obs_score2 = observations(person, variable2);

    pseudolikelihood_ratio += 2 * obs_score1 * obs_score2 * delta_state;

    //variable 1 log pseudolikelihood ratio
    rest_score = residual_matrix(person, variable1) -
      obs_score2 * interactions(variable2, variable1);

    if(rest_score > 0) {
      bound = num_categories[variable1] * rest_score;
    } else {
      bound = 0.0;
    }

    if(is_ordinal_variable[variable1] == true) {
      //Regular binary or ordinal MRF variable ---------------------------------
      denominator_prop = std::exp(-bound);
      denominator_curr = std::exp(-bound);
      for(int category = 0; category < num_categories[variable1]; category++) {
        score = category + 1;
        exponent = thresholds(variable1, category) +
          score * rest_score -
          bound;
        denominator_prop +=
          std::exp(exponent + score * obs_score2 * proposed_state);
        denominator_curr +=
          std::exp(exponent + score * obs_score2 * current_state);
      }
    } else {
      //Blume-Capel ordinal MRF variable ---------------------------------------
      denominator_prop = 0.0;
      denominator_curr = 0.0;
      for(int category = 0; category < num_categories[variable1] + 1; category++) {
        exponent = thresholds(variable1, 0) * category;
        exponent += thresholds(variable1, 1) *
          (category - reference_category[variable1]) *
          (category - reference_category[variable1]);
        exponent+= category * rest_score - bound;
        denominator_prop +=
          std::exp(exponent + category * obs_score2 * proposed_state);
        denominator_curr +=
          std::exp(exponent + category * obs_score2 * current_state);
      }
    }
    pseudolikelihood_ratio -= std::log(denominator_prop);
    pseudolikelihood_ratio += std::log(denominator_curr);

    //variable 2 log pseudolikelihood ratio
    rest_score = residual_matrix(person, variable2) -
      obs_score1 * interactions(variable1, variable2);

    if(rest_score > 0) {
      bound = num_categories[variable2] * rest_score;
    } else {
      bound = 0.0;
    }

    if(is_ordinal_variable[variable2] == true) {
      //Regular binary or ordinal MRF variable ---------------------------------
      denominator_prop = std::exp(-bound);
      denominator_curr = std::exp(-bound);
      for(int category = 0; category < num_categories[variable2]; category++) {
        score = category + 1;
        exponent = thresholds(variable2, category) +
          score * rest_score -
          bound;
        denominator_prop +=
          std::exp(exponent + score * obs_score1 * proposed_state);
        denominator_curr +=
          std::exp(exponent + score * obs_score1 * current_state);
      }
    } else {
      //Blume-Capel ordinal MRF variable ---------------------------------------
      denominator_prop = 0.0;
      denominator_curr = 0.0;
      for(int category = 0; category < num_categories[variable2] + 1; category++) {
        exponent = thresholds(variable2, 0) * category;
        exponent += thresholds(variable2, 1) *
          (category - reference_category[variable2]) *
          (category - reference_category[variable2]);
        exponent+=  category * rest_score - bound;
        denominator_prop +=
          std::exp(exponent + category * obs_score1 * proposed_state);
        denominator_curr +=
          std::exp(exponent + category * obs_score1 * current_state);
      }
    }
    pseudolikelihood_ratio -= std::log(denominator_prop);
    pseudolikelihood_ratio += std::log(denominator_curr);
  }
  return pseudolikelihood_ratio;
}

// ----------------------------------------------------------------------------|
// MH algorithm to sample from the full-conditional of the active interaction
//  parameters for Bayesian edge selection
// ----------------------------------------------------------------------------|
void metropolis_interactions(
    arma::mat& interactions,
    const arma::mat& thresholds,
    const arma::imat& indicator,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    arma::mat& proposal_sd,
    const double interaction_scale,
    const int no_persons,
    const int no_variables,
    arma::mat& residual_matrix,
    const double phi,
    const double target_ar,
    const int t,
    const double epsilon_lo,
    const double epsilon_hi,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category
) {
  double proposed_state;
  double current_state;
  double log_prob;
  double U;

  for(int variable1 = 0; variable1 <  no_variables - 1; variable1++) {
    for(int variable2 = variable1 + 1; variable2 <  no_variables; variable2++) {
      if(indicator(variable1, variable2) == 1) {
        current_state = interactions(variable1, variable2);
        proposed_state = R::rnorm(current_state, proposal_sd(variable1, variable2));

        log_prob = log_pseudolikelihood_ratio_interaction(interactions,
                                                          thresholds,
                                                          observations,
                                                          num_categories,
                                                          no_persons,
                                                          variable1,
                                                          variable2,
                                                          proposed_state,
                                                          current_state,
                                                          residual_matrix,
                                                          is_ordinal_variable,
                                                          reference_category);
        log_prob += R::dcauchy(proposed_state, 0.0, interaction_scale, true);
        log_prob -= R::dcauchy(current_state, 0.0, interaction_scale, true);

        U = R::unif_rand();
        if(std::log(U) < log_prob) {
          double state_diff = proposed_state - current_state;
          interactions(variable1, variable2) = proposed_state;
          interactions(variable2, variable1) = proposed_state;

          //Update the matrix of rest scores
          for(int person = 0; person < no_persons; person++) {
            residual_matrix(person, variable1) += observations(person, variable2) *
              state_diff;
            residual_matrix(person, variable2) += observations(person, variable1) *
              state_diff;
          }
        }

        if(log_prob > 0) {
          log_prob = 1;
        } else {
          log_prob = std::exp(log_prob);
        }

        double update_proposal_sd = proposal_sd(variable1, variable2) +
          (log_prob - target_ar) * std::exp(-log(t) * phi);

        if(std::isnan(update_proposal_sd) == true) {
          update_proposal_sd = 1.0;
        }

        update_proposal_sd = std::clamp(update_proposal_sd, epsilon_lo, epsilon_hi);
        proposal_sd(variable1, variable2) = update_proposal_sd;

      }
    }
  }
}

// ----------------------------------------------------------------------------|
// MH algorithm to sample from the full-conditional of an edge + interaction
//  pair for Bayesian edge selection
// ----------------------------------------------------------------------------|
void metropolis_edge_interaction_pair(
    arma::mat& interactions,
    const arma::mat& thresholds,
    arma::imat& indicator,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::mat& proposal_sd,
    const double interaction_scale,
    const arma::imat& index,
    const int no_interactions,
    const int no_persons,
    arma::mat& residual_matrix,
    const arma::mat& theta,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category
) {
  double proposed_state;
  double current_state;
  double log_prob;
  double U;

  int variable1;
  int variable2;

  for(int cntr = 0; cntr < no_interactions; cntr ++) {
    variable1 = index(cntr, 1);
    variable2 = index(cntr, 2);

    current_state = interactions(variable1, variable2);

    if(indicator(variable1, variable2) == 0) {
      proposed_state = R::rnorm(current_state, proposal_sd(variable1, variable2));
    } else {
      proposed_state = 0.0;
    }

    log_prob = log_pseudolikelihood_ratio_interaction(interactions,
                                                      thresholds,
                                                      observations,
                                                      num_categories,
                                                      no_persons,
                                                      variable1,
                                                      variable2,
                                                      proposed_state,
                                                      current_state,
                                                      residual_matrix,
                                                      is_ordinal_variable,
                                                      reference_category);

    if(indicator(variable1, variable2) == 0) {
      log_prob += R::dcauchy(proposed_state, 0.0, interaction_scale, true);
      log_prob -= R::dnorm(proposed_state,
                           current_state,
                           proposal_sd(variable1, variable2),
                           true);

      log_prob += log(theta(variable1, variable2) / (1 - theta(variable1, variable2)));
    } else {
      log_prob -= R::dcauchy(current_state, 0.0, interaction_scale, true);
      log_prob += R::dnorm(current_state,
                           proposed_state,
                           proposal_sd(variable1, variable2),
                           true);

      log_prob -= log(theta(variable1, variable2) / (1 - theta(variable1, variable2)));
    }

    U = R::unif_rand();
    if(std::log(U) < log_prob) {
      indicator(variable1, variable2) = 1 - indicator(variable1, variable2);
      indicator(variable2, variable1) = 1 - indicator(variable2, variable1);

      interactions(variable1, variable2) = proposed_state;
      interactions(variable2, variable1) = proposed_state;

      double state_diff = proposed_state - current_state;

      //Update the matrix of rest scores ---------------------------------------
      for(int person = 0; person < no_persons; person++) {
        residual_matrix(person, variable1) += observations(person, variable2) *
          state_diff;
        residual_matrix(person, variable2) += observations(person, variable1) *
          state_diff;
      }
    }
  }
}

// ----------------------------------------------------------------------------|
// A Gibbs step for graphical model parameters for Bayesian edge selection
// ----------------------------------------------------------------------------|
List gibbs_step_gm(
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const double interaction_scale,
    arma::mat& proposal_sd,
    arma::mat& proposal_sd_blumecapel,
    const arma::imat& index,
    const arma::imat& num_obs_categories,
    const arma::imat& sufficient_blume_capel,
    const double threshold_alpha,
    const double threshold_beta,
    const int no_persons,
    const int no_variables,
    const int no_interactions,
    const int no_thresholds,
    const int max_num_categories,
    arma::imat& indicator,
    arma::mat& interactions,
    arma::mat& thresholds,
    arma::mat& residual_matrix,
    const arma::mat& theta,
    const  double phi,
    const double target_ar,
    const int t,
    const double epsilon_lo,
    const double epsilon_hi,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    const bool edge_selection
) {

  if(edge_selection == true) {
    //Between model move (update edge indicators and interaction parameters)
    metropolis_edge_interaction_pair(interactions,
                                     thresholds,
                                     indicator,
                                     observations,
                                     num_categories,
                                     proposal_sd,
                                     interaction_scale,
                                     index,
                                     no_interactions,
                                     no_persons,
                                     residual_matrix,
                                     theta,
                                     is_ordinal_variable,
                                     reference_category);
  }

  //Within model move (update interaction parameters)
  metropolis_interactions(interactions,
                          thresholds,
                          indicator,
                          observations,
                          num_categories,
                          proposal_sd,
                          interaction_scale,
                          no_persons,
                          no_variables,
                          residual_matrix,
                          phi,
                          target_ar,
                          t,
                          epsilon_lo,
                          epsilon_hi,
                          is_ordinal_variable,
                          reference_category);

  //Update threshold parameters
  for(int variable = 0; variable < no_variables; variable++) {
    if(is_ordinal_variable[variable] == true) {
      metropolis_thresholds_regular(thresholds,
                                    observations,
                                    num_categories,
                                    num_obs_categories,
                                    no_persons,
                                    variable,
                                    threshold_alpha,
                                    threshold_beta,
                                    residual_matrix);
    } else {
      metropolis_thresholds_blumecapel(thresholds,
                                       observations,
                                       num_categories,
                                       sufficient_blume_capel,
                                       no_persons,
                                       variable,
                                       reference_category,
                                       threshold_alpha,
                                       threshold_beta,
                                       residual_matrix,
                                       proposal_sd_blumecapel,
                                       phi,
                                       target_ar,
                                       t,
                                       epsilon_lo,
                                       epsilon_hi);
    }
  }

  return List::create(Named("indicator") = indicator,
                      Named("interactions") = interactions,
                      Named("thresholds") = thresholds,
                      Named("residual_matrix") = residual_matrix,
                      Named("proposal_sd") = proposal_sd);
}

// ----------------------------------------------------------------------------|
// The Gibbs sampler for Bayesian edge selection
// ----------------------------------------------------------------------------|
// [[Rcpp::export]]
List gibbs_sampler(
    arma::imat& observations,
    arma::imat& indicator,
    arma::mat& interactions,
    arma::mat& thresholds,
    const arma::ivec& num_categories,
    const double interaction_scale,
    arma::mat& proposal_sd,
    arma::mat& proposal_sd_blumecapel,
    const String& edge_prior,
    arma::mat& theta,
    const double beta_bernoulli_alpha,
    const double beta_bernoulli_beta,
    const double dirichlet_alpha,
    const double lambda,
    const arma::imat& Index,
    const int iter,
    const int burnin,
    arma::imat& num_obs_categories,
    arma::imat& sufficient_blume_capel,
    const double threshold_alpha,
    const double threshold_beta,
    const bool na_impute,
    const arma::imat& missing_index,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    const bool save = false,
    const bool display_progress = false,
    bool edge_selection = true
) {
  int cntr;
  int no_variables = observations.n_cols;
  int no_persons = observations.n_rows;
  int no_interactions = Index.n_rows;
  int no_thresholds = sum(num_categories);
  int max_num_categories = max(num_categories);
  arma::ivec K_values;  // To store sampled K values


  arma::ivec v = seq(0, no_interactions - 1);
  arma::ivec order(no_interactions);
  arma::imat index(no_interactions, 3);

  //Parameters of adaptive proposals -------------------------------------------
  double phi = 0.75;
  double target_ar = 0.234;
  double epsilon_lo = 1.0 / static_cast<double>(no_persons);
  double epsilon_hi = 2.0;

  //The resizing based on ``save'' could probably be prettier ------------------
  int nrow = no_variables;
  int ncol_edges = no_variables;
  int ncol_thresholds = max_num_categories;

  if(save == true) {
    nrow = iter;
    ncol_edges= no_interactions;
    ncol_thresholds = no_thresholds;
  }

  arma::mat out_pairwise_effects(nrow, ncol_edges);
  arma::mat out_main_effects(nrow, ncol_thresholds);

  if(edge_selection == false) {
    for(int variable1 = 0; variable1 < no_variables - 1; variable1++) {
      for(int variable2 = variable1 + 1; variable2 < no_variables; variable2++) {
        indicator(variable1, variable2) = 1;
        indicator(variable2, variable1) = 1;
      }
    }
    nrow = 1;
    ncol_edges = 1;
  }
  arma::mat out_indicator(nrow, ncol_edges);

  arma::mat residual_matrix(no_persons, no_variables);
  for(int variable1 = 0; variable1 < no_variables; variable1++) {
    for(int person = 0; person < no_persons; person++) {
      for(int variable2 = 0; variable2 < no_variables; variable2++) {
        residual_matrix(person, variable1) +=
          observations(person, variable2) * interactions(variable2, variable1);
      }
    }
  }

  //Variable declaration edge prior
  arma::ivec cluster_allocations(no_variables);
  arma::mat cluster_prob(1, 1);
  arma::vec log_Vn(1);

  // store the allocation indices for each iteration
  arma::mat out_allocations(iter, no_variables);

  // if(edge_prior == "Stochastic-Block") { // Initial Configuration of the cluster allocations
  //   cluster_allocations[0] = 0;
  //   cluster_allocations[1] = 1;
  //   for(int i = 2; i < no_variables; i++) {
  //     double U = R::unif_rand();
  //     if(U > 0.5){
  //       cluster_allocations[i] = 1;
  //     } else {
  //       cluster_allocations[i] = 0;
  //     }
  //   }
  //
  //   cluster_prob = block_probs_mfm_sbm(cluster_allocations,
  //                                      indicator,
  //                                      no_variables,
  //                                      beta_bernoulli_alpha,
  //                                      beta_bernoulli_beta);
  //
  //   for(int i = 0; i < no_variables - 1; i++) {
  //     for(int j = i + 1; j < no_variables; j++) {
  //       theta(i, j) = cluster_prob(cluster_allocations[i], cluster_allocations[j]);
  //       theta(j, i) = cluster_prob(cluster_allocations[i], cluster_allocations[j]);
  //     }
  //   }
  //
  //   log_Vn = compute_Vn_mfm_sbm(no_variables,
  //                               dirichlet_alpha,
  //                               no_variables + 10,
  //                               lambda);
  // }

  //The Gibbs sampler ----------------------------------------------------------
  //First, we do burn-in iterations---------------------------------------------

  //When edge_selection = true we do 2 * burnin iterations. The first burnin
  // iterations without selection to ensure good starting values, and proposal
  // calibration. The second burnin iterations with selection.

  int first_burnin = burnin;
  int second_burnin = 0;
  if(edge_selection == true)
    second_burnin = burnin;
  bool input_edge_selection = edge_selection;
  edge_selection = false;

  //Progress bar ---------------------------------------------------------------
  Progress p(iter + first_burnin + second_burnin, display_progress);

  for(int iteration = 0; iteration < first_burnin + second_burnin; iteration++) {
    if(iteration >= first_burnin) {
      edge_selection = input_edge_selection;
    }
    if (Progress::check_abort()) {
      return List::create(Named("indicator") = out_indicator,
                          Named("interactions") = out_pairwise_effects,
                          Named("thresholds") = out_main_effects);
    }
    p.increment();

    //Update interactions and model (between model move) -----------------------
    // Create a random ordering of pairwise effects for updating
    arma::uvec order = arma::randperm(v.n_elem);

    for(int cntr = 0; cntr < no_interactions; cntr++) {
      index(cntr, 0) = Index(order[cntr], 0);
      index(cntr, 1) = Index(order[cntr], 1);
      index(cntr, 2) = Index(order[cntr], 2);
    }

    if(na_impute == true) {
      List out = impute_missing_data(interactions,
                                     thresholds,
                                     observations,
                                     num_obs_categories,
                                     sufficient_blume_capel,
                                     num_categories,
                                     residual_matrix,
                                     missing_index,
                                     is_ordinal_variable,
                                     reference_category);

      arma::imat observations = out["observations"];
      arma::imat num_obs_categories = out["num_obs_categories"];
      arma::imat sufficient_blume_capel = out["sufficient_blume_capel"];
      arma::mat residual_matrix = out["residual_matrix"];
    }

    List out = gibbs_step_gm(observations,
                             num_categories,
                             interaction_scale,
                             proposal_sd,
                             proposal_sd_blumecapel,
                             index,
                             num_obs_categories,
                             sufficient_blume_capel,
                             threshold_alpha,
                             threshold_beta,
                             no_persons,
                             no_variables,
                             no_interactions,
                             no_thresholds,
                             max_num_categories,
                             indicator,
                             interactions,
                             thresholds,
                             residual_matrix,
                             theta,
                             phi,
                             target_ar,
                             iteration + 1,
                             epsilon_lo,
                             epsilon_hi,
                             is_ordinal_variable,
                             reference_category,
                             edge_selection);

    arma::imat indicator = out["indicator"];
    arma::mat interactions = out["interactions"];
    arma::mat thresholds = out["thresholds"];
    arma::mat residual_matrix = out["residual_matrix"];
    arma::mat proposal_sd = out["proposal_sd"];

    if(edge_selection == true) {
      if(edge_prior == "Beta-Bernoulli") {
        int sumG = 0;
        for(int i = 0; i < no_variables - 1; i++) {
          for(int j = i + 1; j < no_variables; j++) {
            sumG += indicator(i, j);
          }
        }
        double probability = R::rbeta(beta_bernoulli_alpha + sumG,
                                      beta_bernoulli_beta + no_interactions - sumG);

        for(int i = 0; i < no_variables - 1; i++) {
          for(int j = i + 1; j < no_variables; j++) {
            theta(i, j) = probability;
            theta(j, i) = probability;
          }
        }
      }
      // if(edge_prior == "Stochastic-Block") {
      //   cluster_allocations = block_allocations_mfm_sbm(cluster_allocations,
      //                                                   no_variables,
      //                                                   log_Vn,
      //                                                   cluster_prob,
      //                                                   indicator,
      //                                                   dirichlet_alpha,
      //                                                   beta_bernoulli_alpha,
      //                                                   beta_bernoulli_beta);
      //
      //
      //   cluster_prob = block_probs_mfm_sbm(cluster_allocations,
      //                                      indicator,
      //                                      no_variables,
      //                                      beta_bernoulli_alpha,
      //                                      beta_bernoulli_beta);
      //
      //   for(int i = 0; i < no_variables - 1; i++) {
      //     for(int j = i + 1; j < no_variables; j++) {
      //       theta(i, j) = cluster_prob(cluster_allocations[i], cluster_allocations[j]);
      //       theta(j, i) = cluster_prob(cluster_allocations[i], cluster_allocations[j]);
      //     }
      //   }
      // }
    }
  }
  //To ensure that edge_selection is reinstated to the input value -------------
  edge_selection = input_edge_selection;

  //The post burn-in iterations ------------------------------------------------
  for(int iteration = 0; iteration < iter; iteration++) {
    if (Progress::check_abort()) {
      if(edge_selection == true && edge_prior == "Stochastic-Block") {
        return List::create(Named("indicator") = out_indicator,
                            Named("interactions") = out_pairwise_effects,
                            Named("thresholds") = out_main_effects,
                            Named("allocations") = out_allocations);
      } if(edge_selection == true && edge_prior != "Stochastic-Block") {
        return List::create(Named("indicator") = out_indicator,
                            Named("interactions") = out_pairwise_effects,
                            Named("thresholds") = out_main_effects);
      }
      else {
        return List::create(Named("interactions") = out_pairwise_effects,
                            Named("thresholds") = out_main_effects);
      }
    }
    p.increment();

    //Update interactions and model (between model move) -----------------------
    // Create a random ordering of pairwise effects for updating
    arma::uvec order = arma::randperm(v.n_elem);

    for(int cntr = 0; cntr < no_interactions; cntr++) {
      index(cntr, 0) = Index(order[cntr], 0);
      index(cntr, 1) = Index(order[cntr], 1);
      index(cntr, 2) = Index(order[cntr], 2);
    }

    if(na_impute == true) {
      List out = impute_missing_data(interactions,
                                     thresholds,
                                     observations,
                                     num_obs_categories,
                                     sufficient_blume_capel,
                                     num_categories,
                                     residual_matrix,
                                     missing_index,
                                     is_ordinal_variable,
                                     reference_category);

      arma::imat observations = out["observations"];
      arma::imat num_obs_categories = out["num_obs_categories"];
      arma::imat sufficient_blume_capel = out["sufficient_blume_capel"];
      arma::mat residual_matrix = out["residual_matrix"];
    }

    List out = gibbs_step_gm(observations,
                             num_categories,
                             interaction_scale,
                             proposal_sd,
                             proposal_sd_blumecapel,
                             index,
                             num_obs_categories,
                             sufficient_blume_capel,
                             threshold_alpha,
                             threshold_beta,
                             no_persons,
                             no_variables,
                             no_interactions,
                             no_thresholds,
                             max_num_categories,
                             indicator,
                             interactions,
                             thresholds,
                             residual_matrix,
                             theta,
                             phi,
                             target_ar,
                             iteration + 1,
                             epsilon_lo,
                             epsilon_hi,
                             is_ordinal_variable,
                             reference_category,
                             edge_selection);

    arma::imat indicator = out["indicator"];
    arma::mat interactions = out["interactions"];
    arma::mat thresholds = out["thresholds"];
    arma::mat residual_matrix = out["residual_matrix"];
    arma::mat proposal_sd = out["proposal_sd"];

    if(edge_selection == true) {
      if(edge_prior == "Beta-Bernoulli") {
        int sumG = 0;
        for(int i = 0; i < no_variables - 1; i++) {
          for(int j = i + 1; j < no_variables; j++) {
            sumG += indicator(i, j);
          }
        }
        double probability = R::rbeta(beta_bernoulli_alpha + sumG,
                                      beta_bernoulli_beta + no_interactions - sumG);

        for(int i = 0; i < no_variables - 1; i++) {
          for(int j = i + 1; j < no_variables; j++) {
            theta(i, j) = probability;
            theta(j, i) = probability;
          }
        }
      }
      // if(edge_prior == "Stochastic-Block") {
      //   cluster_allocations = block_allocations_mfm_sbm(cluster_allocations,
      //                                                   no_variables,
      //                                                   log_Vn,
      //                                                   cluster_prob,
      //                                                   indicator,
      //                                                   dirichlet_alpha,
      //                                                   beta_bernoulli_alpha,
      //                                                   beta_bernoulli_beta);
      //
      //
      //   cluster_prob = block_probs_mfm_sbm(cluster_allocations,
      //                                      indicator,
      //                                      no_variables,
      //                                      beta_bernoulli_alpha,
      //                                      beta_bernoulli_beta);
      //
      //
      //   for(int i = 0; i < no_variables - 1; i++) {
      //     for(int j = i + 1; j < no_variables; j++) {
      //       theta(i, j) = cluster_prob(cluster_allocations[i], cluster_allocations[j]);
      //       theta(j, i) = cluster_prob(cluster_allocations[i], cluster_allocations[j]);
      //     }
      //   }
      // }
    }


    //Output -------------------------------------------------------------------
    if(save == true) {
      //Save raw samples -------------------------------------------------------
      cntr = 0;
      for(int variable1 = 0; variable1 < no_variables - 1; variable1++) {
        for(int variable2 = variable1 + 1; variable2 < no_variables; variable2++) {
          if(edge_selection == true) {
            out_indicator(iteration, cntr) = indicator(variable1, variable2);
          }
          out_pairwise_effects(iteration, cntr) = interactions(variable1, variable2);
          cntr++;
        }
      }
      cntr = 0;
      for(int variable = 0; variable < no_variables; variable++) {
        if(is_ordinal_variable[variable] == true) {
          for(int category = 0; category < num_categories[variable]; category++) {
            out_main_effects(iteration, cntr) = thresholds(variable, category);
            cntr++;
          }
        } else {
          out_main_effects(iteration, cntr) = thresholds(variable, 0);
          cntr++;
          out_main_effects(iteration, cntr) = thresholds(variable, 1);
          cntr++;
        }
      }
    } else {
      //Compute running averages -----------------------------------------------
      for(int variable1 = 0; variable1 < no_variables - 1; variable1++) {
        for(int variable2 = variable1 + 1; variable2 < no_variables; variable2++) {
          if(edge_selection == true) {
            out_indicator(variable1, variable2) *= iteration;
            out_indicator(variable1, variable2) += indicator(variable1, variable2);
            out_indicator(variable1, variable2) /= iteration + 1;
            out_indicator(variable2, variable1) = out_indicator(variable1, variable2);
          }

          out_pairwise_effects(variable1, variable2) *= iteration;
          out_pairwise_effects(variable1, variable2) += interactions(variable1, variable2);
          out_pairwise_effects(variable1, variable2) /= iteration + 1;
          out_pairwise_effects(variable2, variable1) = out_pairwise_effects(variable1, variable2);
        }

        if(is_ordinal_variable[variable1] == true) {
          for(int category = 0; category < num_categories[variable1]; category++) {
            out_main_effects(variable1, category) *= iteration;
            out_main_effects(variable1, category) += thresholds(variable1, category);
            out_main_effects(variable1, category) /= iteration + 1;
          }
        } else {
          out_main_effects(variable1, 0) *= iteration;
          out_main_effects(variable1, 0) += thresholds(variable1, 0);
          out_main_effects(variable1, 0) /= iteration + 1;
          out_main_effects(variable1, 1) *= iteration;
          out_main_effects(variable1, 1) += thresholds(variable1, 1);
          out_main_effects(variable1, 1) /= iteration + 1;
        }
      }
      if(is_ordinal_variable[no_variables - 1] == true) {
        for(int category = 0; category < num_categories[no_variables - 1]; category++) {
          out_main_effects(no_variables - 1, category) *= iteration;
          out_main_effects(no_variables - 1, category) += thresholds(no_variables - 1, category);
          out_main_effects(no_variables - 1, category) /= iteration + 1;
        }
      } else {
        out_main_effects(no_variables - 1, 0) *= iteration;
        out_main_effects(no_variables - 1, 0) += thresholds(no_variables - 1, 0);
        out_main_effects(no_variables - 1, 0) /= iteration + 1;
        out_main_effects(no_variables - 1, 1) *= iteration;
        out_main_effects(no_variables - 1, 1) += thresholds(no_variables - 1, 1);
        out_main_effects(no_variables - 1, 1) /= iteration + 1;
      }
    }

    for(int i = 0; i < no_variables; i++) {
      out_allocations(iteration, i) = cluster_allocations[i] + 1;
    }
  }


  if(edge_selection == true) {
    if(edge_prior == "Stochastic-Block"){
      return List::create(Named("indicator") = out_indicator,
                          Named("pairwise_effects") = out_pairwise_effects,
                          Named("main_effects") = out_main_effects,
                          Named("allocations") = out_allocations);
    } else {
      return List::create(Named("indicator") = out_indicator,
                          Named("pairwise_effects") = out_pairwise_effects,
                          Named("main_effects") = out_main_effects);
    }
  } else {
    return List::create(Named("pairwise_effects") = out_pairwise_effects,
                        Named("main_effects") = out_main_effects);
  }
}
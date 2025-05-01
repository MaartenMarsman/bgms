// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppProgress)]]
#include <Rcpp.h>
#include "gibbs_functions_edge_prior.h"
#include <progress.hpp>
#include <progress_bar.hpp>
using namespace Rcpp;



/**
 * Adapts the log step size using dual averaging during MCMC burn-in.
 *
 * Implements the dual averaging algorithm to adaptively update the MALA step size
 * toward a target acceptance rate. Only used during burn-in.
 *
 * Inputs:
 *  - acceptance_probability: Observed acceptance rate at the current iteration.
 *  - iteration: Current iteration number (starting from 1).
 *  - state: Vector of length 3, modified in-place:
 *      * state[0] = current log step size
 *      * state[1] = running average of log step size
 *      * state[2] = running average of acceptance error
 *
 * Modifies:
 *  - `state` vector in-place to update the log step size parameters.
 */
inline void update_step_size_with_dual_averaging(
    const double acceptance_probability,
    const int iteration,
    arma::vec& state
) {
  constexpr double initial_step_size = 0.01;
  const double target_log_step_size = std::log(10.0 * initial_step_size);
  constexpr int stabilization_offset = 10;

  double& log_step_size = state[0];
  double& log_step_size_avg = state[1];
  double& acceptance_error_avg = state[2];

  const double adjusted_iter = iteration + stabilization_offset;
  const double error = 0.574 - acceptance_probability;

  acceptance_error_avg = (1.0 - 1.0 / adjusted_iter) * acceptance_error_avg +
    (1.0 / adjusted_iter) * error;

  log_step_size = target_log_step_size - std::sqrt(static_cast<double>(iteration)) / 0.05 * acceptance_error_avg;

  const double weight = std::pow(static_cast<double>(iteration), -0.75);
  log_step_size_avg = weight * log_step_size + (1.0 - weight) * log_step_size_avg;
}



/**
 * Function: update_step_size_with_robbins_monro
 *
 * Performs Robbins-Monro adaptation of the step size on the log scale.
 * Applies exponential decay to gradually reduce adaptation influence over time.
 *
 * Inputs:
 *  - acceptance_probability: Observed acceptance rate at the current iteration.
 *  - iteration: Current iteration number (starting from 1).
 *  - step_size_mala: Step size to update (in-place).
 *
 * Notes:
 *  - Uses a fixed target acceptance of 0.574 (optimal for MALA).
 *  - Decay rate is 0.75 for stable adaptation.
 */
inline void update_step_size_with_robbins_monro(
    const double acceptance_probability,
    const int iteration,
    double& step_size_mala
) {
  constexpr double target_acceptance = 0.574;
  constexpr double decay_rate = 0.75;

  const double error = acceptance_probability - target_acceptance;
  const double decay = std::pow(static_cast<double>(iteration), -decay_rate);

  double log_step_size = std::log(step_size_mala);
  log_step_size += error * decay;
  step_size_mala = std::exp(log_step_size);
}



/**
 * Function: update_proposal_sd_with_robbins_monro
 * Purpose: Performs Robbins-Monro updates for proposal standard deviations.
 *
 * Inputs:
 *  - current_sd: Current standard deviation of the proposal.
 *  - observed_log_acceptance_probability: Log acceptance probability from the Metropolis-Hastings step.
 *  - rm_weight: Robbins-Monro adaptation weight (e.g. iteration^{-0.75}).
 *
 * Returns:
 *  - Updated proposal standard deviation, clamped within bounds.
 */
inline double update_proposal_sd_with_robbins_monro (
    const double current_sd,
    const double observed_log_acceptance_probability,
    const double rm_weight
) {
  constexpr double target_acceptance = 0.234;
  constexpr double rm_lower_bound = 0.001;
  constexpr double rm_upper_bound = 2.0;

  // Normalize the acceptance probability
  double observed_acceptance_probability = 1.0;
  if (observed_log_acceptance_probability < 0.0) {
    observed_acceptance_probability = std::exp(observed_log_acceptance_probability);
  }

  // Robbins-Monro update step
  double updated_sd = current_sd +
    (observed_acceptance_probability - target_acceptance) * rm_weight;

  // Handle NaNs robustly
  if (std::isnan(updated_sd)) {
    updated_sd = 1.0;
  }

  return std::clamp(updated_sd, rm_lower_bound, rm_upper_bound);
}



/**
 * Function: count_num_main_effects
 *
 * Computes the total number of main effect (threshold) parameters across variables.
 *
 * Ordinal variables contribute one parameter per category.
 * Blume-Capel variables contribute exactly two parameters.
 *
 * Inputs:
 *  - num_categories: Vector of category counts per variable.
 *  - is_ordinal_variable: Logical vector (0 or 1) indicating which variables are ordinal.
 *
 * Returns:
 *  - Total number of main effect parameters to estimate.
 */
inline int count_num_main_effects(
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable
) {
  int n_params = 0;
  for (int i = 0; i < num_categories.n_elem; ++i) {
    n_params += is_ordinal_variable[i] ? num_categories[i] : 2;
  }
  return n_params;
}



/**
 * Function: vectorize_thresholds
 *
 * Converts a matrix of main effect parameters into a flat vector for optimization.
 * Respects the structure of ordinal vs. Blume-Capel variables.
 *
 * Inputs:
 *  - main_effects: Matrix of main effect parameters.
 *  - num_categories: Category count per variable.
 *  - is_ordinal_variable: Logical vector indicating ordinal (1) or Blume-Capel (0).
 *
 * Returns:
 *  - A flat vector of all parameters, in row-major order.
 *    Ordinal: one parameter per category;
 *    Blume-Capel: two parameters (linear and quadratic).
 */
arma::vec vectorize_thresholds(
    const arma::mat& main_effects,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable
) {
  const int num_parameters = count_num_main_effects(num_categories, is_ordinal_variable);
  arma::vec vector(num_parameters);
  int offset = 0;

  for (int variable = 0; variable < main_effects.n_rows; ++variable) {
    const int num_pars = is_ordinal_variable(variable) ? num_categories(variable) : 2;
    vector.subvec(offset, offset + num_pars - 1) =
      main_effects.row(variable).cols(0, num_pars - 1).t();
    offset += num_pars;
  }

  return vector;
}



/**
 * Function: reconstruct_threshold_matrix
 *
 * Reconstructs a threshold matrix from a flat vector of parameters,
 * reversing the operation of `vectorize_thresholds()`.
 *
 * Inputs:
 *  - vector: Flattened vector of threshold parameters.
 *  - num_categories: Vector of category counts per variable.
 *  - is_ordinal_variable: Logical vector (1 for ordinal, 0 for Blume-Capel).
 *
 * Returns:
 *  - A matrix of main effect parameters, with each row corresponding to a variable.
 *    Entries beyond the actual number of parameters per variable are zero-filled.
 */
arma::mat reconstruct_threshold_matrix(
    const arma::vec& vector,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable
) {
  const int num_variables = num_categories.n_elem;
  const int max_categories = num_categories.max();

  arma::mat matrix(num_variables, max_categories, arma::fill::zeros);

  int offset = 0;
  for (int variable = 0; variable < num_variables; ++variable) {
    const int num_pars = is_ordinal_variable[variable] ? num_categories[variable] : 2;
    matrix.row(variable).cols(0, num_pars - 1) =
      vector.subvec(offset, offset + num_pars - 1).t();
    offset += num_pars;
  }

  return matrix;
}



/**
 * Function: impute_missing_values_for_graphical_model
 *
 * Imputes missing values in the observation matrix using the current model parameters.
 * For each missing entry, a category is sampled from the posterior predictive distribution,
 * and if the value changes, all dependent structures are updated accordingly.
 *
 * Inputs:
 *  - pairwise_effects: matrix of interaction weights.
 *  - main_effects: matrix of threshold parameters (main effects).
 *  - missing_index: matrix of (person, variable) indices for missing entries.
 *  - num_categories: vector of category counts per variable.
 *  - is_ordinal_variable: logical vector (0/1) indicating ordinal vs. Blume-Capel.
 *  - reference_category: vector of reference categories per variable.
 *
 * Updates (in-place):
 *  - observations: matrix of imputed categorical scores.
 *  - num_obs_categories: counts of observed categories per variable.
 *  - sufficient_blume_capel: sufficient stats for Blume-Capel thresholds.
 *  - residual_matrix: updated linear predictors after imputation.
 */
void impute_missing_values_for_graphical_model (
    const arma::mat& pairwise_effects,
    const arma::mat& main_effects,
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

  arma::vec category_probabilities(max_num_categories + 1);

  for (int miss = 0; miss < num_missings; miss++) {
    const int person = missing_index(miss, 0);
    const int variable = missing_index(miss, 1);

    const double rest_score = residual_matrix(person, variable);
    const int num_cats = num_categories(variable);
    const bool is_ordinal = is_ordinal_variable(variable);

    double cumsum = 0.0;

    if (is_ordinal) {
      // Compute cumulative unnormalized probabilities for ordinal variable
      cumsum = 1.0;
      category_probabilities[0] = cumsum;
      for (int cat = 0; cat < num_cats; cat++) {
        const int score = cat + 1;
        const double exponent = main_effects(variable, cat) + score * rest_score;
        cumsum += std::exp(exponent);
        category_probabilities[score] = cumsum;
      }
    } else {
      // Compute probabilities for Blume-Capel variable
      const int ref = reference_category(variable);

      cumsum = std::exp(main_effects(variable, 1) * ref * ref);
      category_probabilities[0] = cumsum;

      for (int cat = 0; cat < num_cats; cat++) {
        const int score = cat + 1;
        const int centered = score - ref;
        const double exponent =
          main_effects(variable, 0) * score +
          main_effects(variable, 1) * centered * centered +
          score * rest_score;
        cumsum += std::exp(exponent);
        category_probabilities[score] = cumsum;
      }
    }

    // Sample from categorical distribution via inverse transform
    const double u = R::unif_rand() * cumsum;
    int sampled_score = 0;
    while (u > category_probabilities[sampled_score]) {
      sampled_score++;
    }

    const int new_value = sampled_score;
    const int old_value = observations(person, variable);

    if (new_value != old_value) {
      // Update observation matrix
      observations(person, variable) = new_value;

      if (is_ordinal) {
        num_obs_categories(old_value, variable)--;
        num_obs_categories(new_value, variable)++;
      } else {
        const int ref = reference_category(variable);
        const int delta = new_value - old_value;
        const int delta_sq =
          (new_value - ref) * (new_value - ref) -
          (old_value - ref) * (old_value - ref);

        sufficient_blume_capel(0, variable) += delta;
        sufficient_blume_capel(1, variable) += delta_sq;
      }

      // Update residuals across all variables
      for (int var = 0; var < num_variables; var++) {
        const double delta_score = (new_value - old_value) * pairwise_effects(var, variable);
        residual_matrix(person, var) += delta_score;
      }
    }
  }
}


/**
 * Function: log_pseudoposterior_thresholds
 *
 * Computes the log pseudo-posterior for all main effect parameters.
 *
 * Inputs:
 *  - main_effects: Matrix of threshold parameters.
 *  - residual_matrix: Residual scores for each observation and variable.
 *  - num_categories: Number of categories per variable.
 *  - num_obs_categories: Observed category count matrix.
 *  - sufficient_blume_capel: Sufficient statistics for Blume-Capel variables.
 *  - reference_category: Reference category per variable (for Blume-Capel).
 *  - is_ordinal_variable: Logical vector (1 = ordinal, 0 = Blume-Capel).
 *  - threshold_alpha, threshold_beta: Prior hyperparameters.
 *
 * Returns:
 *  - Scalar log pseudo-posterior.
 */
double log_pseudoposterior_thresholds(
    const arma::mat& main_effects,
    const arma::mat& residual_matrix,
    const arma::ivec& num_categories,
    const arma::imat& num_obs_categories,
    const arma::imat& sufficient_blume_capel,
    const arma::ivec& reference_category,
    const arma::uvec& is_ordinal_variable,
    const double threshold_alpha = 1.0,
    const double threshold_beta = 1.0
) {
  const int num_variables = residual_matrix.n_cols;
  const int num_persons = residual_matrix.n_rows;

  double log_posterior = 0.0;

  for (int variable = 0; variable < num_variables; ++variable) {
    const int num_cats = num_categories(variable);

    if (is_ordinal_variable(variable)) {
      // Regular ordinal variable
      for (int cat = 0; cat < num_cats; ++cat) {
        const double theta = main_effects(variable, cat);
        log_posterior += theta * (num_obs_categories(cat + 1, variable) + threshold_alpha);
        log_posterior -= std::log1p(std::exp(theta)) * (threshold_alpha + threshold_beta);
      }

      for (int person = 0; person < num_persons; ++person) {
        const double rest_score = residual_matrix(person, variable);
        const double bound = num_cats * rest_score;

        double denom = std::exp(-bound);
        for (int cat = 0; cat < num_cats; ++cat) {
          const double exponent = main_effects(variable, cat) + (cat + 1) * rest_score - bound;
          denom += std::exp(exponent);
        }

        log_posterior -= bound + std::log(denom);
      }

    } else {
      // Blume-Capel variable
      const double theta_lin = main_effects(variable, 0);
      const double theta_quad = main_effects(variable, 1);

      log_posterior += theta_lin * (sufficient_blume_capel(0, variable) + threshold_alpha);
      log_posterior -= std::log1p(std::exp(theta_lin)) * (threshold_alpha + threshold_beta);

      log_posterior += theta_quad * (sufficient_blume_capel(1, variable) + threshold_alpha);
      log_posterior -= std::log1p(std::exp(theta_quad)) * (threshold_alpha + threshold_beta);

      const int ref = reference_category(variable);

      for (int person = 0; person < num_persons; ++person) {
        const double rest_score = residual_matrix(person, variable);
        const double bound = num_cats * rest_score;

        double denom = 0.0;
        for (int cat = 0; cat <= num_cats; ++cat) {
          const int centered = cat - ref;
          const double exponent =
            theta_lin * cat +
            theta_quad * centered * centered +
            cat * rest_score - bound;
          denom += std::exp(exponent);
        }

        log_posterior -= bound + std::log(denom);
      }
    }
  }

  return log_posterior;
}



/**
 * Function: gradient_log_pseudoposterior_thresholds
 *
 * Computes the gradient of the log pseudo-posterior for threshold parameters.
 * Includes both data likelihood and logistic-Beta prior components.
 *
 * Inputs:
 *  - main_effects: Matrix of threshold parameters.
 *  - residual_matrix: Residual scores for each observation and variable.
 *  - num_categories: Number of categories per variable.
 *  - num_obs_categories: Observed category count matrix.
 *  - sufficient_blume_capel: Sufficient statistics for Blume-Capel variables.
 *  - reference_category: Reference category per variable (for Blume-Capel).
 *  - is_ordinal_variable: Logical vector (1 = ordinal, 0 = Blume-Capel).
 *  - threshold_alpha, threshold_beta: Prior hyperparameters.
 *
 * Returns:
 *  - A flat gradient vector in the same order as vectorized thresholds.
 */
arma::vec gradient_log_pseudoposterior_thresholds(
    const arma::mat& main_effects,
    const arma::mat& residual_matrix,
    const arma::ivec& num_categories,
    const arma::imat& num_obs_categories,
    const arma::imat& sufficient_blume_capel,
    const arma::ivec& reference_category,
    const arma::uvec& is_ordinal_variable,
    const double threshold_alpha = 1.0,
    const double threshold_beta = 1.0
) {
  const int num_variables = residual_matrix.n_cols;
  const int num_persons = residual_matrix.n_rows;
  const int num_parameters = count_num_main_effects(num_categories, is_ordinal_variable);

  arma::vec gradient(num_parameters, arma::fill::zeros);
  int offset = 0;

  for (int variable = 0; variable < num_variables; ++variable) {
    const int num_cats = num_categories(variable);

    if (is_ordinal_variable(variable)) {
      const double max_theta = main_effects.row(variable).max();

      for (int cat = 0; cat < num_cats; ++cat) {
        gradient(offset + cat) = num_obs_categories(cat + 1, variable);
      }

      for (int person = 0; person < num_persons; ++person) {
        const double rest_score = residual_matrix(person, variable);
        const double bound = max_theta + num_cats * rest_score;

        double denom = std::exp(-bound);
        arma::vec numerators(num_cats, arma::fill::zeros);

        for (int cat = 0; cat < num_cats; ++cat) {
          const double exponent = main_effects(variable, cat) + (cat + 1) * rest_score - bound;
          numerators(cat) = std::exp(exponent);
          denom += numerators(cat);
        }

        for (int cat = 0; cat < num_cats; ++cat) {
          gradient(offset + cat) -= numerators(cat) / denom;
        }
      }

      for (int cat = 0; cat < num_cats; ++cat) {
        const double theta = main_effects(variable, cat);
        const double p = 1.0 / (1.0 + std::exp(-theta));
        gradient(offset + cat) += threshold_alpha - (threshold_alpha + threshold_beta) * p;
      }

      offset += num_cats;

    } else {
      // Blume-Capel variable
      const int ref = reference_category(variable);
      const double theta_lin = main_effects(variable, 0);
      const double theta_quad = main_effects(variable, 1);

      gradient(offset)     = sufficient_blume_capel(0, variable);
      gradient(offset + 1) = sufficient_blume_capel(1, variable);

      for (int person = 0; person < num_persons; ++person) {
        const double rest_score = residual_matrix(person, variable);
        const double bound = num_cats * rest_score;

        double denom = std::exp(theta_quad * ref * ref - bound);
        double sum_lin = 0.0;
        double sum_quad = ref * ref * denom;

        for (int cat = 0; cat < num_cats; ++cat) {
          const int score = cat + 1;
          const int centered = score - ref;

          const double exponent =
            theta_lin * score +
            theta_quad * centered * centered +
            score * rest_score - bound;

          const double weight = std::exp(exponent);
          sum_lin  += weight * score;
          sum_quad += weight * centered * centered;
          denom    += weight;
        }

        gradient(offset)     -= sum_lin / denom;
        gradient(offset + 1) -= sum_quad / denom;
      }

      for (int i = 0; i < 2; ++i) {
        const double theta = main_effects(variable, i);
        const double p = 1.0 / (1.0 + std::exp(-theta));
        gradient(offset + i) += threshold_alpha - (threshold_alpha + threshold_beta) * p;
      }

      offset += 2;
    }
  }

  return gradient;
}



/**
 * Function: update_thresholds_with_adaptive_mala
 *
 * Performs a MALA update of threshold parameters with adaptive step size tuning.
 * Applies dual averaging during burn-in and Robbins-Monro afterward.
 *
 * Inputs:
 *  - main_effects: Matrix of threshold parameters (updated in-place).
 *  - step_size_mala: MALA step size (updated in-place).
 *  - residual_matrix: Residual scores for each observation and variable.
 *  - num_categories: Number of categories per variable.
 *  - num_obs_categories: Observed category count matrix.
 *  - sufficient_blume_capel: Sufficient statistics for Blume-Capel variables.
 *  - reference_category: Reference category per variable (for Blume-Capel).
 *  - is_ordinal_variable: Logical vector (1 = ordinal, 0 = Blume-Capel).
 *  - iteration: Current iteration number.
 *  - burnin: Total number of burn-in iterations.
 *  - dual_averaging_state: Dual averaging state vector [log_eps, log_eps_avg, H_bar] (updated in-place).
 *  - threshold_alpha, threshold_beta: Prior hyperparameters.
 *
 * Modifies:
 *  - main_effects
 *  - step_size_mala
 *  - dual_averaging_state
 */
void update_thresholds_with_adaptive_mala(
    arma::mat& main_effects,
    double& step_size_mala,
    const arma::mat& residual_matrix,
    const arma::ivec& num_categories,
    const arma::imat& num_obs_categories,
    const arma::imat& sufficient_blume_capel,
    const arma::ivec& reference_category,
    const arma::uvec& is_ordinal_variable,
    const int iteration,
    const int burnin,
    arma::vec& dual_averaging_state,
    const double threshold_alpha,
    const double threshold_beta
) {
  arma::vec flat_theta = vectorize_thresholds(main_effects, num_categories, is_ordinal_variable);

  arma::vec grad = gradient_log_pseudoposterior_thresholds(
    main_effects, residual_matrix, num_categories, num_obs_categories,
    sufficient_blume_capel, reference_category, is_ordinal_variable,
    threshold_alpha, threshold_beta
  );

  const double log_post_current = log_pseudoposterior_thresholds(
    main_effects, residual_matrix, num_categories, num_obs_categories,
    sufficient_blume_capel, reference_category, is_ordinal_variable,
    threshold_alpha, threshold_beta
  );

  const double sqrt_step = std::sqrt(step_size_mala);
  arma::vec proposal = flat_theta + 0.5 * step_size_mala * grad + sqrt_step * arma::randn(flat_theta.n_elem);
  arma::mat proposed_thresholds = reconstruct_threshold_matrix(proposal, num_categories, is_ordinal_variable);

  const double log_post_proposal = log_pseudoposterior_thresholds(
    proposed_thresholds, residual_matrix, num_categories, num_obs_categories,
    sufficient_blume_capel, reference_category, is_ordinal_variable,
    threshold_alpha, threshold_beta
  );

  arma::vec grad_proposal = gradient_log_pseudoposterior_thresholds(
    proposed_thresholds, residual_matrix, num_categories, num_obs_categories,
    sufficient_blume_capel, reference_category, is_ordinal_variable,
    threshold_alpha, threshold_beta
  );

  const arma::vec forward_mean = flat_theta + 0.5 * step_size_mala * grad;
  const arma::vec backward_mean = proposal + 0.5 * step_size_mala * grad_proposal;

  const double log_forward = -0.5 / step_size_mala * arma::accu(arma::square(proposal - forward_mean));
  const double log_backward = -0.5 / step_size_mala * arma::accu(arma::square(flat_theta - backward_mean));
  const double log_acceptance = log_post_proposal + log_backward - log_post_current - log_forward;

  // Accept proposal
  if (std::log(R::unif_rand()) < log_acceptance) {
    main_effects = proposed_thresholds;
    flat_theta = proposal;
  }

  const double accept_prob = std::min(1.0, std::exp(log_acceptance));

  if (iteration <= burnin) {
    update_step_size_with_dual_averaging(accept_prob, iteration, dual_averaging_state);
    step_size_mala = std::exp(dual_averaging_state[0]);
  } else {
    update_step_size_with_robbins_monro(accept_prob, iteration - burnin, step_size_mala);
  }
}



/**
 * Function: update_regular_thresholds_with_metropolis
 *
 * Performs a Metropolis-Hastings update for each threshold of an ordinal variable.
 * Uses a generalized beta-prime proposal and logistic-Beta prior.
 *
 * Inputs:
 *  - main_effects: Matrix of thresholds (updated in-place).
 *  - observations: Matrix of category scores.
 *  - num_categories: Vector of number of categories per variable.
 *  - num_obs_categories: Count matrix of observed scores.
 *  - no_persons: Number of individuals.
 *  - variable: Index of the variable being updated.
 *  - threshold_alpha, threshold_beta: Prior hyperparameters.
 *  - residual_matrix: Residual scores for each observation and variable.
 *
 * Modifies:
 *  - main_effects (only for the specified variable)
 */
void update_regular_thresholds_with_metropolis(
    arma::mat& main_effects,
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
  const int num_cats = num_categories(variable);

  for (int category = 0; category < num_cats; ++category) {
    double current = main_effects(variable, category);
    double exp_current = std::exp(current);
    double c = (threshold_alpha + threshold_beta) / (1.0 + exp_current);

    for (int person = 0; person < no_persons; ++person) {
      double rest_score = residual_matrix(person, variable);
      double denom = 1.0;
      double numer = std::exp((category + 1) * rest_score);

      for (int cat = 0; cat < num_cats; ++cat) {
        if (cat != category) {
          denom += std::exp(main_effects(variable, cat) + (cat + 1) * rest_score);
        }
      }

      g(person) = denom;
      q(person) = numer;
      c += numer / (denom + numer * exp_current);
    }

    c /= (no_persons + threshold_alpha + threshold_beta - exp_current * c);

    // Sample from generalized beta-prime proposal
    double a = num_obs_categories(category + 1, variable) + threshold_alpha;
    double b = no_persons + threshold_beta - num_obs_categories(category + 1, variable);
    double tmp = R::rbeta(a, b);
    double proposed = std::log(tmp / (1.0 - tmp) / c);
    double exp_proposed = std::exp(proposed);

    // Compute MH acceptance probability
    double log_acceptance_probability = 0.0;
    for (int person = 0; person < no_persons; ++person) {
      log_acceptance_probability += std::log(g(person) + q(person) * exp_current);
      log_acceptance_probability -= std::log(g(person) + q(person) * exp_proposed);
    }

    log_acceptance_probability -= (threshold_alpha + threshold_beta) * std::log1p(exp_proposed);
    log_acceptance_probability += (threshold_alpha + threshold_beta) * std::log1p(exp_current);
    log_acceptance_probability -= (a + b) * std::log1p(c * exp_current);
    log_acceptance_probability += (a + b) * std::log1p(c * exp_proposed);

    if (std::log(R::unif_rand()) < log_acceptance_probability) {
      main_effects(variable, category) = proposed;
    }
  }
}



/**
 * Function: update_blumecapel_thresholds_with_adaptive_metropolis
 *
 * Performs an adaptive Metropolis update of the Blume-Capel threshold parameters
 * (linear and quadratic) for a single variable, with Robbins-Monro tuning.
 *
 * Inputs:
 *  - main_effects: Matrix of threshold parameters (updated in-place).
 *  - observations: Matrix of categorical scores.
 *  - num_categories: Number of categories per variable.
 *  - sufficient_blume_capel: Sufficient statistics for Blume-Capel variables.
 *  - no_persons: Number of observations.
 *  - variable: Index of the variable being updated.
 *  - reference_category: Reference category per variable.
 *  - threshold_alpha, threshold_beta: Prior hyperparameters.
 *  - residual_matrix: Residual scores.
 *  - proposal_sd_blumecapel: Matrix of proposal SDs for each variable (updated in-place).
 *  - exp_neg_log_t_rm_adaptation_rate: Robbins-Monro adaptation weight.
 *
 * Modifies:
 *  - main_effects (for the given variable)
 *  - proposal_sd_blumecapel
 */
void update_blumecapel_thresholds_with_adaptive_metropolis(
    arma::mat& main_effects,
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
    const double exp_neg_log_t_rm_adaptation_rate
) {
  const int num_cats = num_categories(variable);
  const int ref = reference_category(variable);

  auto update_parameter = [&](int param_index) {
    const double current = main_effects(variable, param_index);
    const double proposed = R::rnorm(current, proposal_sd_blumecapel(variable, param_index));
    const double diff = proposed - current;

    arma::vec numer_current(num_cats + 1);
    arma::vec numer_proposed(num_cats + 1);

    for (int cat = 0; cat <= num_cats; ++cat) {
      const int centered = cat - ref;

      if (param_index == 0) {
        double quad_term = main_effects(variable, 1) * (centered * centered);
        numer_current(cat) = current * cat + quad_term;
        numer_proposed(cat) = proposed * cat + quad_term;
      } else {
        double lin_term = main_effects(variable, 0) * cat;
        numer_current(cat) = current * (centered * centered) + lin_term;
        numer_proposed(cat) = proposed * (centered * centered) + lin_term;
      }
    }

    const double max_curr = numer_current.max();
    const double max_prop = numer_proposed.max();
    const double lbound = (max_curr > 0.0 || max_prop > 0.0) ? std::max(max_curr, max_prop) : 0.0;

    double log_accept = threshold_alpha * diff +
      static_cast<double>(sufficient_blume_capel(param_index, variable)) * diff;

    for (int person = 0; person < no_persons; ++person) {
      const double rest_score = residual_matrix(person, variable);
      const double bound = (rest_score > 0.0) ? num_cats * rest_score + lbound : lbound;

      double denom_curr = std::exp(numer_current(0) - bound);
      double denom_prop = std::exp(numer_proposed(0) - bound);

      for (int cat = 0; cat < num_cats; ++cat) {
        const double score_term = (cat + 1) * rest_score - bound;
        denom_curr += std::exp(numer_current(cat + 1) + score_term);
        denom_prop += std::exp(numer_proposed(cat + 1) + score_term);
      }

      log_accept += std::log(denom_curr) - std::log(denom_prop);
    }

    // Prior ratio (logistic-Beta)
    log_accept += (threshold_alpha + threshold_beta) *
      (std::log1p(std::exp(current)) - std::log1p(std::exp(proposed)));

    // Metropolis step
    if (std::log(R::unif_rand()) < log_accept) {
      main_effects(variable, param_index) = proposed;
    }

    // Robbins-Monro adaptation
    proposal_sd_blumecapel(variable, param_index) =
      update_proposal_sd_with_robbins_monro(
        proposal_sd_blumecapel(variable, param_index),
        log_accept,
        exp_neg_log_t_rm_adaptation_rate
      );
  };

  update_parameter(0); // linear
  update_parameter(1); // quadratic
}



/**
 * Function: log_pseudolikelihood_ratio_interaction
 *
 * Computes the log pseudo-likelihood ratio between a proposed and current interaction value.
 *
 * Inputs:
 *  - pairwise_effects: Current matrix of pairwise interaction parameters.
 *  - main_effects: Matrix of threshold parameters.
 *  - observations: Matrix of category scores.
 *  - num_categories: Number of categories per variable.
 *  - no_persons: Number of observations.
 *  - variable1, variable2: Indices of the interacting variables.
 *  - proposed_state: Proposed new value for the interaction.
 *  - current_state: Current value of the interaction.
 *  - residual_matrix: Residual scores matrix.
 *  - is_ordinal_variable: Logical vector indicating which variables are ordinal.
 *  - reference_category: Reference category per variable (Blume-Capel).
 *
 * Returns:
 *  - log p(y | β_proposed) - log p(y | β_current)
 */
double log_pseudolikelihood_ratio_interaction(
    const arma::mat& pairwise_effects,
    const arma::mat& main_effects,
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
  double log_ratio = 0.0;
  const double delta = proposed_state - current_state;

  for (int person = 0; person < no_persons; ++person) {
    const int score1 = observations(person, variable1);
    const int score2 = observations(person, variable2);

    log_ratio += 2.0 * score1 * score2 * delta;

    // ---- Variable 1 contribution ----
    {
      double rest_score = residual_matrix(person, variable1) - score2 * current_state;
      double bound = rest_score > 0.0 ? num_categories(variable1) * rest_score : 0.0;

      double denom_curr = 0.0;
      double denom_prop = 0.0;

      if (is_ordinal_variable(variable1)) {
        denom_curr += std::exp(-bound);
        denom_prop += std::exp(-bound);

        for (int cat = 0; cat < num_categories(variable1); ++cat) {
          double exponent = main_effects(variable1, cat) + (cat + 1) * rest_score;
          denom_curr += std::exp(exponent + (cat + 1) * score2 * current_state - bound);
          denom_prop += std::exp(exponent + (cat + 1) * score2 * proposed_state - bound);
        }
      } else {
        const int ref = reference_category(variable1);

        for (int cat = 0; cat <= num_categories(variable1); ++cat) {
          const int centered = cat - ref;
          const double exponent =
            main_effects(variable1, 0) * cat +
            main_effects(variable1, 1) * centered * centered +
            cat * rest_score - bound;

          denom_curr += std::exp(exponent + cat * score2 * current_state);
          denom_prop += std::exp(exponent + cat * score2 * proposed_state);
        }
      }

      log_ratio += std::log(denom_curr) - std::log(denom_prop);
    }

    // ---- Variable 2 contribution ----
    {
      double rest_score = residual_matrix(person, variable2) - score1 * current_state;
      double bound = rest_score > 0.0 ? num_categories(variable2) * rest_score : 0.0;

      double denom_curr = 0.0;
      double denom_prop = 0.0;

      if (is_ordinal_variable(variable2)) {
        denom_curr += std::exp(-bound);
        denom_prop += std::exp(-bound);

        for (int cat = 0; cat < num_categories(variable2); ++cat) {
          double exponent = main_effects(variable2, cat) + (cat + 1) * rest_score;
          denom_curr += std::exp(exponent + (cat + 1) * score1 * current_state - bound);
          denom_prop += std::exp(exponent + (cat + 1) * score1 * proposed_state - bound);
        }
      } else {
        const int ref = reference_category(variable2);

        for (int cat = 0; cat <= num_categories(variable2); ++cat) {
          const int centered = cat - ref;
          const double exponent =
            main_effects(variable2, 0) * cat +
            main_effects(variable2, 1) * centered * centered +
            cat * rest_score - bound;

          denom_curr += std::exp(exponent + cat * score1 * current_state);
          denom_prop += std::exp(exponent + cat * score1 * proposed_state);
        }
      }

      log_ratio += std::log(denom_curr) - std::log(denom_prop);
    }
  }

  return log_ratio;
}



/**
 * Function: update_interactions_with_adaptive_metropolis
 *
 * Performs adaptive Metropolis-Hastings updates for all active pairwise interactions.
 *
 * For each pair (i, j) with inclusion_indicator(i, j) == 1, proposes a new interaction strength.
 * If accepted, updates the interaction matrix and residual matrix.
 *
 * Inputs:
 *  - pairwise_effects: Current matrix of interaction parameters (updated in-place).
 *  - main_effects: Matrix of main effect (threshold) parameters.
 *  - inclusion_indicator: Binary matrix indicating active interactions.
 *  - observations: Matrix of category scores.
 *  - num_categories: Number of categories per variable.
 *  - proposal_sd_pairwise_effects: Matrix of proposal standard deviations (updated in-place).
 *  - interaction_scale: Scale parameter for the Cauchy prior.
 *  - num_persons: Number of observations.
 *  - num_variables: Number of variables.
 *  - residual_matrix: Matrix of residual scores (updated in-place).
 *  - exp_neg_log_t_rm_adaptation_rate: Robbins-Monro adaptation weight.
 *  - is_ordinal_variable: Logical vector (1 = ordinal, 0 = Blume-Capel).
 *  - reference_category: Reference category per variable (Blume-Capel).
 *
 * Modifies:
 *  - pairwise_effects
 *  - residual_matrix
 *  - proposal_sd_pairwise_effects
 */
void update_interactions_with_adaptive_metropolis(
    arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    const arma::imat& inclusion_indicator,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    arma::mat& proposal_sd_pairwise_effects,
    const double interaction_scale,
    const int num_persons,
    const int num_variables,
    arma::mat& residual_matrix,
    const double exp_neg_log_t_rm_adaptation_rate,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category
) {
  for (int variable1 = 0; variable1 < num_variables - 1; ++variable1) {
    for (int variable2 = variable1 + 1; variable2 < num_variables; ++variable2) {
      if (inclusion_indicator(variable1, variable2) == 1) {
        const double current_state = pairwise_effects(variable1, variable2);
        const double proposed_state = R::rnorm(current_state, proposal_sd_pairwise_effects(variable1, variable2));

        double log_acceptance = log_pseudolikelihood_ratio_interaction(
          pairwise_effects, main_effects, observations, num_categories, num_persons,
          variable1, variable2, proposed_state, current_state,
          residual_matrix, is_ordinal_variable, reference_category
        );

        // Prior ratio (Cauchy)
        log_acceptance += R::dcauchy(proposed_state, 0.0, interaction_scale, true);
        log_acceptance -= R::dcauchy(current_state, 0.0, interaction_scale, true);

        if (std::log(R::unif_rand()) < log_acceptance) {
          const double delta = proposed_state - current_state;
          pairwise_effects(variable1, variable2) = proposed_state;
          pairwise_effects(variable2, variable1) = proposed_state;

          for (int person = 0; person < num_persons; ++person) {
            residual_matrix(person, variable1) += observations(person, variable2) * delta;
            residual_matrix(person, variable2) += observations(person, variable1) * delta;
          }
        }

        // Robbins-Monro proposal SD update
        proposal_sd_pairwise_effects(variable1, variable2) =
          update_proposal_sd_with_robbins_monro(
            proposal_sd_pairwise_effects(variable1, variable2),
            log_acceptance,
            exp_neg_log_t_rm_adaptation_rate
          );
      }
    }
  }
}



/**
 * Function: update_indicator_interaction_pair_with_metropolis
 *
 * Metropolis-Hastings update of pairwise inclusion indicators for a predefined set of edges.
 *
 * Inputs:
 *  - pairwise_effects: Matrix of interaction weights (updated in-place).
 *  - main_effects: Matrix of main effect (threshold) parameters.
 *  - indicator: Matrix of edge inclusion flags (updated in-place).
 *  - observations: Matrix of category scores.
 *  - num_categories: Number of categories per variable.
 *  - proposal_sd: Matrix of proposal standard deviations for pairwise effects.
 *  - interaction_scale: Scale parameter for the Cauchy prior.
 *  - index: List of interaction pairs to update.
 *  - no_interactions: Number of interaction pairs.
 *  - no_persons: Number of observations.
 *  - residual_matrix: Residual scores matrix (updated in-place).
 *  - theta: Matrix of prior inclusion probabilities.
 *  - is_ordinal_variable: Logical vector indicating variable type.
 *  - reference_category: Reference category per variable (Blume-Capel).
 *
 * Modifies:
 *  - indicator
 *  - pairwise_effects
 *  - residual_matrix
 */
void update_indicator_interaction_pair_with_metropolis(
    arma::mat& pairwise_effects,
    const arma::mat& main_effects,
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
  for (int cntr = 0; cntr < no_interactions; ++cntr) {
    const int variable1 = index(cntr, 1);
    const int variable2 = index(cntr, 2);

    const double current_state = pairwise_effects(variable1, variable2);
    double proposed_state;

    bool proposing_addition = (indicator(variable1, variable2) == 0);
    proposed_state = proposing_addition
    ? R::rnorm(current_state, proposal_sd(variable1, variable2))
      : 0.0;

    double log_accept = log_pseudolikelihood_ratio_interaction(
      pairwise_effects, main_effects, observations, num_categories, no_persons,
      variable1, variable2, proposed_state, current_state,
      residual_matrix, is_ordinal_variable, reference_category
    );

    if (proposing_addition) {
      log_accept += R::dcauchy(proposed_state, 0.0, interaction_scale, true);
      log_accept -= R::dnorm(proposed_state, current_state, proposal_sd(variable1, variable2), true);
      log_accept += std::log(theta(variable1, variable2)) - std::log(1.0 - theta(variable1, variable2));
    } else {
      log_accept -= R::dcauchy(current_state, 0.0, interaction_scale, true);
      log_accept += R::dnorm(current_state, proposed_state, proposal_sd(variable1, variable2), true);
      log_accept -= std::log(theta(variable1, variable2)) - std::log(1.0 - theta(variable1, variable2));
    }

    if (std::log(R::unif_rand()) < log_accept) {
      const int updated_indicator = 1 - indicator(variable1, variable2);
      indicator(variable1, variable2) = updated_indicator;
      indicator(variable2, variable1) = updated_indicator;

      pairwise_effects(variable1, variable2) = proposed_state;
      pairwise_effects(variable2, variable1) = proposed_state;

      const double delta = proposed_state - current_state;

      for (int person = 0; person < no_persons; ++person) {
        residual_matrix(person, variable1) += observations(person, variable2) * delta;
        residual_matrix(person, variable2) += observations(person, variable1) * delta;
      }
    }
  }
}



/**
 * Function: gibbs_update_step_for_graphical_model_parameters
 *
 * Performs a single Gibbs update step for all graphical model parameters.
 *
 * This step updates:
 *  - Main effects (threshold parameters), using either MALA or Metropolis-Hastings
 *  - Pairwise interaction parameters using adaptive Metropolis-Hastings
 *  - Edge inclusion indicators (if edge selection is enabled)
 *
 * Adaptive tuning is applied via Robbins-Monro and dual averaging where applicable.
 *
 * Inputs:
 *  - All model parameters, update indices, residuals, and proposal SDs.
 *  - Current iteration and burn-in info for adaptive tuning.
 *  - Logical flags for edge selection and MALA usage.
 *
 * Modifies:
 *  - main_effects
 *  - pairwise_effects
 *  - inclusion_indicator (if edge selection is active)
 *  - residual_matrix
 *  - proposal_sd_main_effects
 *  - proposal_sd_pairwise_effects
 *  - step_size_mala and dual_averaging_state (if MALA is used)
 */
void gibbs_update_step_for_graphical_model_parameters(
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const double interaction_scale,
    arma::mat& proposal_sd_pairwise_effects,
    arma::mat& proposal_sd_main_effects,
    const arma::imat& update_index,
    const arma::imat& num_obs_categories,
    const arma::imat& sufficient_blume_capel,
    const double threshold_alpha,
    const double threshold_beta,
    const int num_persons,
    const int num_variables,
    const int num_pairwise,
    const int num_main,
    const int max_categories,
    arma::imat& inclusion_indicator,
    arma::mat& pairwise_effects,
    arma::mat& main_effects,
    arma::mat& residual_matrix,
    const arma::mat& theta,
    const double rm_decay_rate,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    const bool edge_selection,
    double& step_size_mala,
    const int iteration,
    arma::vec& dual_averaging_state,
    const int total_burnin,
    const bool use_mala_for_main_effects
) {
  const double exp_neg_log_t_rm_adaptation_rate =
    std::exp(-std::log(static_cast<double>(iteration)) * rm_decay_rate);

  // 1. Update inclusion indicators (edge selection)
  if (edge_selection) {
    update_indicator_interaction_pair_with_metropolis(
      pairwise_effects, main_effects, inclusion_indicator, observations,
      num_categories, proposal_sd_pairwise_effects, interaction_scale,
      update_index, num_pairwise, num_persons, residual_matrix,
      theta, is_ordinal_variable, reference_category
    );
  }

  // 2. Update interaction parameters
  update_interactions_with_adaptive_metropolis(
    pairwise_effects, main_effects, inclusion_indicator, observations,
    num_categories, proposal_sd_pairwise_effects, interaction_scale,
    num_persons, num_variables, residual_matrix,
    exp_neg_log_t_rm_adaptation_rate, is_ordinal_variable, reference_category
  );

  // 3. Update main effect (threshold) parameters
  if (use_mala_for_main_effects) {
    update_thresholds_with_adaptive_mala(
      main_effects, step_size_mala, residual_matrix, num_categories,
      num_obs_categories, sufficient_blume_capel, reference_category,
      is_ordinal_variable, iteration, total_burnin, dual_averaging_state,
      threshold_alpha, threshold_beta
    );
  } else {
    for (int variable = 0; variable < num_variables; ++variable) {
      if (is_ordinal_variable(variable)) {
        update_regular_thresholds_with_metropolis(
          main_effects, observations, num_categories, num_obs_categories,
          num_persons, variable, threshold_alpha, threshold_beta,
          residual_matrix
        );
      } else {
        update_blumecapel_thresholds_with_adaptive_metropolis(
          main_effects, observations, num_categories, sufficient_blume_capel,
          num_persons, variable, reference_category, threshold_alpha,
          threshold_beta, residual_matrix, proposal_sd_main_effects,
          exp_neg_log_t_rm_adaptation_rate
        );
      }
    }
  }
}



/**
 * Function: run_gibbs_sampler_for_bgm
 *
 * Runs the full Gibbs sampling algorithm for a graphical model with ordinal
 * and/or Blume-Capel variables, including optional edge selection and imputation.
 *
 * Each iteration updates:
 *  - Main effects (thresholds) using MALA or Metropolis-Hastings
 *  - Pairwise interactions using adaptive Metropolis-Hastings
 *  - Edge inclusion indicators (if edge selection is enabled)
 *  - Missing values (if imputation is enabled)
 *
 * Proposal distributions are tuned adaptively using Robbins-Monro and dual averaging.
 *
 * Inputs:
 *  - observations: Matrix of categorical scores (with possible missing entries).
 *  - num_categories: Number of categories per variable.
 *  - interaction_scale: Scale of the Cauchy prior for pairwise effects.
 *  - edge_prior: Prior type: "Bernoulli", "Beta-Bernoulli", or "Stochastic-Block".
 *  - theta: Matrix of inclusion probabilities under the edge prior.
 *  - beta_bernoulli_alpha, beta_bernoulli_beta: Beta prior parameters (if applicable).
 *  - dirichlet_alpha, lambda: Parameters for the Stochastic Block Model prior.
 *  - Index: Matrix of interaction edge indices.
 *  - iter, burnin: Number of MCMC iterations and burn-in iterations.
 *  - num_obs_categories: Observed score counts by category.
 *  - sufficient_blume_capel: Sufficient statistics for Blume-Capel variables.
 *  - threshold_alpha, threshold_beta: Hyperparameters for the main effect priors.
 *  - na_impute: Whether to impute missing values.
 *  - missing_index: Matrix of missing (row, col) positions.
 *  - is_ordinal_variable: Logical indicator for ordinal (1) vs. Blume-Capel (0).
 *  - reference_category: Reference category per variable (for Blume-Capel).
 *  - save_main, save_pairwise, save_indicator: Whether to save MCMC samples.
 *  - display_progress: Show progress bar during sampling.
 *  - use_mala_for_main_effects: If TRUE, uses MALA; otherwise Metropolis.
 *
 * Returns:
 *  - A List with:
 *      * "main": Posterior mean of main effects
 *      * "pairwise": Posterior mean of pairwise effects
 *      * "inclusion_indicator": Posterior mean of edge indicators
 *  - If enabled:
 *      * "main_samples", "pairwise_samples", "inclusion_indicator_samples"
 *  - If SBM prior is used:
 *      * "allocations": Cluster membership per iteration
 */
// [[Rcpp::export]]
List run_gibbs_sampler_for_bgm(
    arma::imat& observations,
    const arma::ivec& num_categories,
    const double interaction_scale,
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
    const bool save_main = false,
    const bool save_pairwise = false,
    const bool save_indicator = false,
    const bool display_progress = false,
    bool edge_selection = true,
    bool use_mala_for_main_effects = false
) {
  const int num_variables = observations.n_cols;
  const int num_persons = observations.n_rows;
  const int max_num_categories = num_categories.max();
  const int num_pairwise = Index.n_rows;

  arma::mat main_effects(num_variables, max_num_categories, arma::fill::zeros);
  arma::mat pairwise_effects(num_variables, num_variables, arma::fill::zeros);
  arma::imat inclusion_indicator(num_variables, num_variables, arma::fill::ones);

  arma::mat posterior_mean_main(num_variables, max_num_categories, arma::fill::zeros);
  arma::mat posterior_mean_pairwise(num_variables, num_variables, arma::fill::zeros);
  arma::mat posterior_mean_indicator(num_variables, num_variables, arma::fill::zeros);
  arma::mat residual_matrix(num_persons, num_variables, arma::fill::zeros);

  const int num_main = count_num_main_effects(num_categories, is_ordinal_variable);
  arma::mat* main_effect_samples = nullptr;
  arma::mat* pairwise_effect_samples = nullptr;
  arma::imat* indicator_samples = nullptr;

  if (save_main) main_effect_samples = new arma::mat(iter, num_main);
  if (save_pairwise) pairwise_effect_samples = new arma::mat(iter, num_pairwise);
  if (save_indicator) indicator_samples = new arma::imat(iter, num_pairwise);

  arma::mat proposal_sd_blumecapel(num_main, 2, arma::fill::ones);
  arma::mat proposal_sd_pairwise_effects(num_variables, num_variables, arma::fill::ones);

  double step_size_mala = 0.01;
  arma::vec dual_averaging_state(3, arma::fill::zeros);
  const double decay_rate_robins_monro = 0.75;

  arma::uvec v = arma::regspace<arma::uvec>(0, num_pairwise - 1);
  arma::uvec order(num_pairwise);
  arma::imat index(num_pairwise, 3);

  arma::uvec K_values;
  arma::uvec cluster_allocations(num_variables);
  arma::mat cluster_prob(1, 1);
  arma::vec log_Vn(1);
  arma::imat out_allocations(iter, num_variables);

  if (edge_prior == "Stochastic-Block") {
    cluster_allocations[0] = 0;
    cluster_allocations[1] = 1;
    for (int i = 2; i < num_variables; i++) {
      cluster_allocations[i] = (R::unif_rand() > 0.5) ? 1 : 0;
    }

    cluster_prob = block_probs_mfm_sbm(
      cluster_allocations, arma::conv_to<arma::umat>::from(inclusion_indicator),
      num_variables, beta_bernoulli_alpha, beta_bernoulli_beta
    );

    for (int i = 0; i < num_variables - 1; i++) {
      for (int j = i + 1; j < num_variables; j++) {
        theta(i, j) = cluster_prob(cluster_allocations[i], cluster_allocations[j]);
        theta(j, i) = theta(i, j);
      }
    }

    log_Vn = compute_Vn_mfm_sbm(num_variables, dirichlet_alpha, num_variables + 10, lambda);
  }

  bool enable_edge_selection = edge_selection;
  const int total_burnin = burnin * (enable_edge_selection ? 2 : 1);
  edge_selection = false;
  const int total_iter = total_burnin + iter;
  Progress p(total_iter, display_progress);

  for (int iteration = 0; iteration < total_iter; iteration++) {
    if (Progress::check_abort()) {
      return List::create(
        Named("main") = posterior_mean_main,
        Named("pairwise") = posterior_mean_pairwise,
        Named("inclusion_indicator") = posterior_mean_indicator
      );
    }
    p.increment();

    if (enable_edge_selection && iteration == burnin) edge_selection = true;

    order = arma::randperm(num_pairwise);
    for (int i = 0; i < num_pairwise; i++) {
      index.row(i) = Index.row(order(i));
    }

    if (na_impute) {
      impute_missing_values_for_graphical_model(
        pairwise_effects, main_effects, observations, num_obs_categories,
        sufficient_blume_capel, num_categories, residual_matrix,
        missing_index, is_ordinal_variable, reference_category
      );
    }

    gibbs_update_step_for_graphical_model_parameters(
      observations, num_categories, interaction_scale,
      proposal_sd_pairwise_effects, proposal_sd_blumecapel, index,
      num_obs_categories, sufficient_blume_capel, threshold_alpha,
      threshold_beta, num_persons, num_variables, num_pairwise, num_main,
      max_num_categories, inclusion_indicator, pairwise_effects, main_effects,
      residual_matrix, theta, decay_rate_robins_monro, is_ordinal_variable,
      reference_category, edge_selection, step_size_mala, iteration,
      dual_averaging_state, total_burnin, use_mala_for_main_effects
    );

    if (edge_selection) {
      if (edge_prior == "Beta-Bernoulli") {
        int num_edges_included = 0;
        for (int i = 0; i < num_variables - 1; i++)
          for (int j = i + 1; j < num_variables; j++)
            num_edges_included += inclusion_indicator(i, j);

        double prob = R::rbeta(
          beta_bernoulli_alpha + num_edges_included,
          beta_bernoulli_beta + num_pairwise - num_edges_included
        );

        for (int i = 0; i < num_variables - 1; i++)
          for (int j = i + 1; j < num_variables; j++)
            theta(i, j) = theta(j, i) = prob;

      } else if (edge_prior == "Stochastic-Block") {
        cluster_allocations = block_allocations_mfm_sbm(
          cluster_allocations, num_variables, log_Vn, cluster_prob,
          arma::conv_to<arma::umat>::from(inclusion_indicator), dirichlet_alpha,
          beta_bernoulli_alpha, beta_bernoulli_beta
        );

        cluster_prob = block_probs_mfm_sbm(
          cluster_allocations,
          arma::conv_to<arma::umat>::from(inclusion_indicator), num_variables,
          beta_bernoulli_alpha, beta_bernoulli_beta
        );

        for (int i = 0; i < num_variables - 1; i++) {
          for (int j = i + 1; j < num_variables; j++) {
            theta(i, j) = theta(j, i) = cluster_prob(cluster_allocations[i], cluster_allocations[j]);
          }
        }
      }
    }

    if (iteration >= total_burnin) {
      int iter_adj = iteration - total_burnin + 1;

      posterior_mean_main = (posterior_mean_main * (iter_adj - 1) + main_effects) / iter_adj;
      posterior_mean_pairwise = (posterior_mean_pairwise * (iter_adj - 1) + pairwise_effects) / iter_adj;

      if (edge_selection) {
        posterior_mean_indicator = (posterior_mean_indicator * (iter_adj - 1) +
          arma::conv_to<arma::mat>::from(inclusion_indicator)) / iter_adj;
      }

      int sample_index = iteration - total_burnin;

      if (save_main) {
        arma::vec vectorized_main = vectorize_thresholds(main_effects, num_categories, is_ordinal_variable);
        main_effect_samples->row(sample_index) = vectorized_main.t();
      }

      if (save_pairwise) {
        arma::vec vectorized_pairwise(num_pairwise);
        for (int i = 0; i < num_pairwise; ++i) {
          vectorized_pairwise(i) = pairwise_effects(Index(i, 1), Index(i, 2));
        }
        pairwise_effect_samples->row(sample_index) = vectorized_pairwise.t();
      }

      if (save_indicator) {
        arma::ivec vectorized_indicator(num_pairwise);
        for (int i = 0; i < num_pairwise; ++i) {
          vectorized_indicator(i) = inclusion_indicator(Index(i, 1), Index(i, 2));
        }
        indicator_samples->row(sample_index) = vectorized_indicator.t();
      }

      if (edge_prior == "Stochastic-Block") {
        for (int j = 0; j < num_variables; ++j)
          out_allocations(sample_index, j) = cluster_allocations[j] + 1;
      }
    }
  }

  List out = List::create(
    Named("main") = posterior_mean_main,
    Named("pairwise") = posterior_mean_pairwise,
    Named("inclusion_indicator") = posterior_mean_indicator
  );

  if (save_main) {
    out["main_samples"] = *main_effect_samples;
    delete main_effect_samples;
  }
  if (save_pairwise) {
    out["pairwise_samples"] = *pairwise_effect_samples;
    delete pairwise_effect_samples;
  }
  if (save_indicator) {
    out["inclusion_indicator_samples"] = *indicator_samples;
    delete indicator_samples;
  }
  if (edge_prior == "Stochastic-Block") {
    out["allocations"] = out_allocations;
  }

  return out;
}
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
inline void update_step_size_with_dual_averaging (
    const double initial_step_size,
    const double acceptance_probability,
    const int iteration,
    arma::vec& state
) {
  const double target_log_step_size = std::log (10.0 * initial_step_size);
  constexpr int stabilization_offset = 10;

  double& log_step_size = state[0];
  double& log_step_size_avg = state[1];
  double& acceptance_error_avg = state[2];

  const double adjusted_iter = iteration + stabilization_offset;
  const double error = 0.574 - acceptance_probability;

  acceptance_error_avg = (1.0 - 1.0 / adjusted_iter) * acceptance_error_avg +
    (1.0 / adjusted_iter) * error;

  log_step_size = target_log_step_size - std::sqrt (static_cast<double> (iteration)) / 0.05 * acceptance_error_avg;

  const double weight = std::pow (static_cast<double> (iteration), -0.75);
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
inline void update_step_size_with_robbins_monro (
    const double acceptance_probability,
    const int iteration,
    double& step_size_mala
) {
  constexpr double target_acceptance = 0.574;
  constexpr double decay_rate = 0.75;

  const double error = acceptance_probability - target_acceptance;
  const double decay = std::pow (static_cast<double> (iteration), -decay_rate);

  double log_step_size = std::log (step_size_mala);
  log_step_size += error * decay;
  step_size_mala = std::exp (log_step_size);
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
    observed_acceptance_probability = std::exp (observed_log_acceptance_probability);
  }

  // Robbins-Monro update step
  double updated_sd = current_sd +
    (observed_acceptance_probability - target_acceptance) * rm_weight;

  // Handle NaNs robustly
  if (std::isnan (updated_sd)) {
    updated_sd = 1.0;
  }

  return std::clamp (updated_sd, rm_lower_bound, rm_upper_bound);
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
inline int count_num_main_effects (
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable
) {
  int n_params = 0;
  for (int i = 0; i < num_categories.n_elem; i++) {
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
arma::vec vectorize_thresholds (
    const arma::mat& main_effects,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable
) {
  const int num_parameters = count_num_main_effects (num_categories, is_ordinal_variable);
  arma::vec vector (num_parameters);
  int offset = 0;

  for (int variable = 0; variable < main_effects.n_rows; variable++) {
    const int num_pars = is_ordinal_variable (variable) ? num_categories(variable) : 2;
    vector.subvec (offset, offset + num_pars - 1) =
      main_effects.row (variable).cols (0, num_pars - 1).t ();
    offset += num_pars;
  }

  return vector;
}



/**
 * Function: unvectorize_thresholds
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
arma::mat unvectorize_thresholds (
    const arma::vec& vector,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable
) {
  const int num_variables = num_categories.n_elem;
  const int max_categories = num_categories.max ();

  arma::mat matrix (num_variables, max_categories, arma::fill::zeros);

  int offset = 0;
  for (int variable = 0; variable < num_variables; variable++) {
    const int num_pars = is_ordinal_variable[variable] ? num_categories[variable] : 2;
    matrix.row (variable).cols (0, num_pars - 1) =
      vector.subvec (offset, offset + num_pars - 1).t ();
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
  const int max_num_categories = num_categories.max ();

  arma::vec category_probabilities (max_num_categories + 1);

  for (int miss = 0; miss < num_missings; miss++) {
    const int person = missing_index (miss, 0);
    const int variable = missing_index (miss, 1);

    const double rest_score = residual_matrix (person, variable);
    const int num_cats = num_categories (variable);
    const bool is_ordinal = is_ordinal_variable (variable);

    double cumsum = 0.0;

    if (is_ordinal) {
      // Compute cumulative unnormalized probabilities for ordinal variable
      cumsum = 1.0;
      category_probabilities[0] = cumsum;
      for (int cat = 0; cat < num_cats; cat++) {
        const int score = cat + 1;
        const double exponent = main_effects (variable, cat) + score * rest_score;
        cumsum += std::exp (exponent);
        category_probabilities[score] = cumsum;
      }
    } else {
      // Compute probabilities for Blume-Capel variable
      const int ref = reference_category (variable);

      cumsum = std::exp (main_effects (variable, 1) * ref * ref);
      category_probabilities[0] = cumsum;

      for (int cat = 0; cat < num_cats; cat++) {
        const int score = cat + 1;
        const int centered = score - ref;
        const double exponent =
          main_effects (variable, 0) * score +
          main_effects (variable, 1) * centered * centered +
          score * rest_score;
        cumsum += std::exp (exponent);
        category_probabilities[score] = cumsum;
      }
    }

    // Sample from categorical distribution via inverse transform
    const double u = R::unif_rand () * cumsum;
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
double log_pseudoposterior_thresholds (
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

  auto log_beta_prior = [&](double theta) {
    return theta * threshold_alpha - std::log1p (std::exp (theta)) * (threshold_alpha + threshold_beta);
  };

  for (int variable = 0; variable < num_variables; variable++) {
    const int num_cats = num_categories(variable);

    if (is_ordinal_variable(variable)) {
      // Prior contribution + sufficient statistic
      for (int cat = 0; cat < num_cats; cat++) {
        const double theta = main_effects(variable, cat);
        log_posterior += theta * num_obs_categories(cat + 1, variable);
        log_posterior += log_beta_prior (theta);
      }

      // Vectorized likelihood contribution
      // For each person, we compute the unnormalized log-likelihood denominator:
      //   denom = exp (-bound) + sum_c exp (theta_c + (c+1) * rest_score - bound)
      // Where:
      //   - rest_score is the summed interaction score excluding the variable itself
      //   - bound = num_cats * rest_score (for numerical stability)
      //   - theta_c is the threshold parameter for category c (0-based)
      arma::vec rest_score = residual_matrix.col (variable);                     // rest scores for all persons
      arma::vec bound = num_cats * rest_score;                                  // numerical bound vector
      arma::vec denom = arma::exp (-bound);                                      // initialize with base term
      arma::vec theta = main_effects.row (variable).cols (0, num_cats - 1).t ();   // threshold parameters

      for (int cat = 0; cat < num_cats; cat++) {
        arma::vec exponent = theta(cat) + (cat + 1) * rest_score - bound;       // exponent per person
        denom += arma::exp (exponent);                                           // accumulate exp terms
      }

      // We then compute the total log-likelihood contribution as:
      //   log_posterior -= bound + log (denom), summed over all persons
      log_posterior -= arma::accu (bound + arma::log (denom));                    // total contribution
    } else {
      const double theta_lin = main_effects(variable, 0);
      const double theta_quad = main_effects(variable, 1);
      const int ref = reference_category(variable);

      // Prior contribution + sufficient statistic
      log_posterior += theta_lin * sufficient_blume_capel(0, variable);
      log_posterior += log_beta_prior(theta_lin);
      log_posterior += theta_quad * sufficient_blume_capel(1, variable);
      log_posterior += log_beta_prior(theta_quad);;

      // Vectorized likelihood contribution
      // For each person, we compute the unnormalized log-likelihood denominator:
      //   denom = sum_c exp (θ_lin * c + θ_quad * (c - ref)^2 + c * rest_score - bound)
      // Where:
      //   - θ_lin, θ_quad are linear and quadratic thresholds
      //   - ref is the reference category (used for centering)
      //   - bound = num_cats * rest_score (stabilizes exponentials)
      arma::vec rest_score = residual_matrix.col(variable);                     // rest scores for all persons
      arma::vec bound = num_cats * rest_score;                                  // numerical bound vector
      arma::vec denom(num_persons, arma::fill::zeros);                          // initialize denominator
      for (int cat = 0; cat <= num_cats; cat++) {
        int centered = cat - ref;                                               // centered category
        double quad_term = theta_quad * centered * centered;                    // precompute quadratic term
        double lin_term = theta_lin * cat;                                      // precompute linear term

        arma::vec exponent = lin_term + quad_term + cat * rest_score - bound;
        denom += arma::exp (exponent);                                           // accumulate over categories
      }

      // The final log-likelihood contribution is then:
      //   log_posterior -= bound + log (denom), summed over all persons
      log_posterior -= arma::accu (bound + arma::log (denom));                    // total contribution
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
arma::vec gradient_log_pseudoposterior_thresholds (
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

  for (int variable = 0; variable < num_variables; variable++) {
    const int num_cats = num_categories(variable);

    if (is_ordinal_variable(variable)) {
      for (int cat = 0; cat < num_cats; cat++) {
        gradient(offset + cat) = num_obs_categories(cat + 1, variable);
      }

      // Vectorized computation of expected category counts
      //
      // For each person, we compute softmax-like probabilities over categories:
      //   probs[p, c] = exp (θ_c + (c+1) * rest_score_p - bound_p)
      // where:
      //   - θ_c is the threshold for category c
      //   - bound_p = max(θ) + num_cats * rest_score_p (numerical stabilization)
      // We normalize across categories and subtract expected counts from the gradient.

      arma::vec rest_score = residual_matrix.col(variable);                     // rest scores per person
      arma::vec theta = main_effects.row(variable).cols(0, num_cats - 1).t();   // thresholds for current variable
      arma::vec bound = theta.max() + num_cats * rest_score;                    // vector of bounds per person

      arma::mat exponents(num_persons, num_cats);                               // log unnormalized probabilities
      for (int cat = 0; cat < num_cats; cat++) {
        exponents.col(cat) = theta(cat) + (cat + 1) * rest_score - bound;
      }

      arma::mat probs = arma::exp (exponents);                                   // unnormalized probabilities
      arma::vec denom = arma::sum(probs, 1) + arma::exp (-bound);                // normalization constants per person

      // Accumulate gradient contributions by subtracting expected counts
      for (int cat = 0; cat < num_cats; cat++) {
        arma::vec normalized = probs.col(cat) / denom;                          // normalized prob for category
        gradient(offset + cat) -= arma::accu (normalized);                       // accumulate gradient
      }

      // Compute prior contribution to gradient
      for (int cat = 0; cat < num_cats; cat++) {
        const double p = 1.0 / (1.0 + std::exp (-theta(cat)));
        gradient(offset + cat) += threshold_alpha - (threshold_alpha + threshold_beta) * p;
      }

      offset += num_cats;

    } else {
      // Blume-Capel variable
      const int ref = reference_category(variable);
      const double theta_lin = main_effects(variable, 0);
      const double theta_quad = main_effects(variable, 1);

      gradient(offset) = sufficient_blume_capel(0, variable);
      gradient(offset + 1) = sufficient_blume_capel(1, variable);

      // Vectorized computation of expected statistics
      //
      // For each person, compute the expected sufficient statistics:
      //   - sum_lin  = E[score]
      //   - sum_quad = E[(score - ref)^2]
      // where probabilities are softmax-like terms derived from θ_lin and θ_quad.
      //
      // This replaces the nested loop with vectorized accumulation over categories.
      arma::vec rest_score = residual_matrix.col(variable);                     // Residuals per person
      arma::vec bound = num_cats * rest_score;                                  // Stabilization bound
      arma::vec denom = arma::exp (theta_quad * ref * ref - bound);              // Initial term at score = 0

      arma::vec sum_lin(num_persons, arma::fill::zeros);                        // E[score]
      arma::vec sum_quad = ref * ref * denom;                                   // E[(score - ref)^2], starts at score = 0

      for (int cat = 0; cat < num_cats; cat++) {
        int score = cat + 1;
        int centered = score - ref;

        double lin_term = theta_lin * score;
        double quad_term = theta_quad * centered * centered;

        arma::vec exponent = lin_term + quad_term + score * rest_score - bound;
        arma::vec weight = arma::exp (exponent);                                 // Unnormalized probabilities

        sum_lin += weight * score;                                              // Accumulate score-weighted terms
        sum_quad += weight * centered * centered;                               // Accumulate centered^2-weighted terms
        denom += weight;                                                        // Update normalization constant
      }

      // Finalize the gradient updates
      gradient(offset) -= arma::accu (sum_lin / denom);                          // Gradient for θ_lin
      gradient(offset + 1) -= arma::accu (sum_quad / denom);                     // Gradient for θ_quad


      // Compute prior contribution to gradient
      for (int i = 0; i < 2; i++) {
        const double theta = main_effects(variable, i);
        const double p = 1.0 / (1.0 + std::exp (-theta));
        gradient(offset + i) += threshold_alpha - (threshold_alpha + threshold_beta) * p;
      }

      offset += 2;
    }
  }

  return gradient;
}



/**
 * Function: find_reasonable_initial_step_size_mala_thresholds
 *
 * Heuristically selects an initial step size for MALA updates of the threshold parameters
 * by adapting the acceptance rate toward a target value (~0.574). The procedure
 * increases or decreases the log step size until the Metropolis-Hastings acceptance
 * crosses the threshold, using a method adapted from Algorithm 4 in:
 *
 *   Hoffman, M.D. & Gelman, A. (2014).
 *   The No-U-Turn Sampler: Adaptively Setting Path Lengths in Hamiltonian Monte Carlo.
 *   Journal of Machine Learning Research, 15, 1593–1623.
 *
 * Inputs:
 *  - main_effects: Initial threshold matrix.
 *  - residual_matrix: Current residual predictor matrix.
 *  - num_categories, num_obs_categories: Category structure for each variable.
 *  - sufficient_blume_capel, reference_category: Data model components.
 *  - is_ordinal_variable: Logical vector indicating ordinal predictors.
 *  - threshold_alpha, threshold_beta: Prior hyperparameters for thresholds.
 *
 * Returns:
 *  - A positive scalar step size tuned for MALA with target acceptance ~0.574.
 *
 * Notes:
 *  - Uses log-space doubling/halving; returns 0.01 if the search fails.
 *  - Assumes Euclidean (identity) preconditioning.
 */
double find_reasonable_initial_step_size_mala_thresholds (
    arma::mat main_effects,
    const arma::mat& residual_matrix,
    const arma::ivec& num_categories,
    const arma::imat& num_obs_categories,
    const arma::imat& sufficient_blume_capel,
    const arma::ivec& reference_category,
    const arma::uvec& is_ordinal_variable,
    const double threshold_alpha,
    const double threshold_beta
) {
  constexpr double initial_step_size = 0.1;
  constexpr double target_acceptance = 0.574;
  constexpr int max_attempts = 20;
  constexpr double max_log_step = 10.0;

  double log_step_size = std::log (initial_step_size);

  // Current state and log posterior
  const arma::vec current_state = vectorize_thresholds(main_effects, num_categories, is_ordinal_variable);
  const arma::vec current_grad = gradient_log_pseudoposterior_thresholds(
    main_effects, residual_matrix, num_categories, num_obs_categories,
    sufficient_blume_capel, reference_category, is_ordinal_variable,
    threshold_alpha, threshold_beta
  );
  const double current_log_post = log_pseudoposterior_thresholds(
    main_effects, residual_matrix, num_categories, num_obs_categories,
    sufficient_blume_capel, reference_category, is_ordinal_variable,
    threshold_alpha, threshold_beta
  );

  int direction = 0;
  double accept_prob = 0.0;

  // Propose initial MALA step: θ' = θ + ½ε ∇log p(θ) + √ε * N(0, I)
  {
    const double step_size = std::exp (log_step_size);
    const double sqrt_step = std::sqrt(step_size);
    const arma::vec proposed_state = current_state + 0.5 * step_size * current_grad +
      sqrt_step * arma::randn(current_state.n_elem);
    const arma::mat proposed_main_effects = unvectorize_thresholds(
      proposed_state, num_categories, is_ordinal_variable
    );

    const double proposed_log_post = log_pseudoposterior_thresholds(
      proposed_main_effects, residual_matrix, num_categories, num_obs_categories,
      sufficient_blume_capel, reference_category, is_ordinal_variable,
      threshold_alpha, threshold_beta
    );
    const arma::vec proposed_grad = gradient_log_pseudoposterior_thresholds(
      proposed_main_effects, residual_matrix, num_categories, num_obs_categories,
      sufficient_blume_capel, reference_category, is_ordinal_variable,
      threshold_alpha, threshold_beta
    );

    const arma::vec forward = current_state + 0.5 * step_size * current_grad;
    const arma::vec backward = proposed_state + 0.5 * step_size * proposed_grad;

    const double log_fwd = -0.5 / step_size * arma::accu (arma::square(proposed_state - forward));
    const double log_bwd = -0.5 / step_size * arma::accu (arma::square(current_state - backward));
    const double log_accept = proposed_log_post + log_bwd - current_log_post - log_fwd;

    accept_prob = std::min(1.0, std::exp (log_accept));
    if (std::abs(accept_prob - target_acceptance) < 0.1) {
      return std::exp (log_step_size);
    }
    direction = (accept_prob > target_acceptance) ? 1 : -1;
  }

  // Log-scale search for step size that brackets the target acceptance rate
  for (int attempt = 0; attempt < max_attempts; attempt++) {
    log_step_size += direction;
    const double step_size = std::exp (log_step_size);
    const double sqrt_step = std::sqrt(step_size);

    const arma::vec proposed_state = current_state + 0.5 * step_size * current_grad +
      sqrt_step * arma::randn(current_state.n_elem);
    const arma::mat proposed_main_effects = unvectorize_thresholds(
      proposed_state, num_categories, is_ordinal_variable
    );

    const double proposed_log_post = log_pseudoposterior_thresholds(
      proposed_main_effects, residual_matrix, num_categories, num_obs_categories,
      sufficient_blume_capel, reference_category, is_ordinal_variable,
      threshold_alpha, threshold_beta
    );
    const arma::vec proposed_grad = gradient_log_pseudoposterior_thresholds(
      proposed_main_effects, residual_matrix, num_categories, num_obs_categories,
      sufficient_blume_capel, reference_category, is_ordinal_variable,
      threshold_alpha, threshold_beta
    );

    const arma::vec forward = current_state + 0.5 * step_size * current_grad;
    const arma::vec backward = proposed_state + 0.5 * step_size * proposed_grad;

    const double log_fwd = -0.5 / step_size * arma::accu (arma::square(proposed_state - forward));
    const double log_bwd = -0.5 / step_size * arma::accu (arma::square(current_state - backward));
    const double log_accept = proposed_log_post + log_bwd - current_log_post - log_fwd;

    const double new_accept_prob = std::min(1.0, std::exp (log_accept));

    // Exit if acceptance flips across the target
    if ((direction == 1 && new_accept_prob < target_acceptance) ||
        (direction == -1 && new_accept_prob > target_acceptance)) {
      break;
    }

    accept_prob = new_accept_prob;

    if (std::abs(log_step_size) > max_log_step) {
      Rcpp::Rcout << "Warning: Step size search failed. Falling back to default (0.01)." << std::endl;
      return 0.01;
    }
  }

  return std::exp (log_step_size);
}



/**
 * Function: initialize_fisher_preconditioner
 *
 * Initializes the square root of the inverse Fisher information matrix,
 * based on the outer product of the gradient. This serves as the starting
 * preconditioner for Fisher-MALA after burn-in.
 *
 * Based on:
 *   Titsias, M.K. (2023).
 *   Gradient-Based MCMC Using Preconditioned Langevin Dynamics.
 *   JMLR, 24(216):1–40.
 *
 * Inputs:
 *  - grad: Gradient vector at the initial state.
 *
 * Returns:
 *  - Square root of the inverse Fisher approximation (d × d matrix).
 *
 * Notes:
 *  - A damping parameter is used to regularize the preconditioner.
 */
inline arma::mat initialize_fisher_preconditioner(
    const arma::vec& grad
) {
  constexpr double damping_par = 10.0;

  const int dim = grad.n_elem;
  const double inner = arma::dot (grad, grad);
  const arma::mat outer = grad * grad.t();

  // Shrinkage ratio: balances curvature and damping
  const double shrinkage = 1.0 / (1.0 + std::sqrt(damping_par / (damping_par + inner)));

  // Fisher inverse root approximation (Titsias 2023, Eq. 12 form)
  arma::mat sqrt_inv_fisher =
    (1.0 / std::sqrt(damping_par)) *
    (arma::eye(dim, dim) - shrinkage * outer / (damping_par + inner));

  return sqrt_inv_fisher;
}



/**
 * Function: update_fisher_preconditioner
 *
 * Performs a rank-1 update of the inverse Fisher matrix approximation,
 * as described in:
 *
 *   Titsias, M.K. (2023).
 *   Gradient-Based MCMC Using Preconditioned Langevin Dynamics.
 *   JMLR, 24(216):1–40.
 *
 * Inputs:
 *  - sqrt_inv_fisher: Current square root inverse Fisher matrix (updated in-place).
 *  - score_diff: Rao-Blackwellized score difference vector (∝ proposed_grad - current_grad).
 *
 * Notes:
 *  - The method maintains an efficient low-rank structure (O(d²)).
 *  - Intended to be used after burn-in.
 */
inline void update_fisher_preconditioner(
    arma::mat& sqrt_inv_fisher,
    const arma::vec& score_diff
) {
  // Transform score difference into scaled preconditioned space
  const arma::vec phi = sqrt_inv_fisher.t() * score_diff;
  const double inner = arma::dot (phi, phi);
  const arma::mat outer = phi * phi.t();

  // Compute the Titsias update scaling factor
  const double r = 1.0 / (1.0 + std::sqrt(1.0 / (1.0 + inner)));

  // Apply rank-1 update to the inverse Fisher root
  sqrt_inv_fisher -= r * sqrt_inv_fisher * outer / (1.0 + inner);
}



/**
 * Function: update_thresholds_with_fisher_mala
 *
 * Performs a Fisher-preconditioned MALA update of the threshold parameters.
 * Uses dual averaging during burn-in and Robbins-Monro adaptation afterward.
 * Proposal is based on inverse Fisher matrix approximation following:
 *
 *   Titsias, M.K. (2023).
 *   Gradient-Based MCMC Using Preconditioned Langevin Dynamics.
 *   Journal of Machine Learning Research, 24(216):1–40.
 *
 * Modifies:
 *  - main_effects, step_size, dual_averaging_state, sqrt_inv_fisher
 */
void update_thresholds_with_fisher_mala (
    arma::mat& main_effects,
    double& step_size,
    const arma::mat& residual_matrix,
    const arma::ivec& num_categories,
    const arma::imat& num_obs_categories,
    const arma::imat& sufficient_blume_capel,
    const arma::ivec& reference_category,
    const arma::uvec& is_ordinal_variable,
    const int iteration,
    const int total_burnin,
    arma::vec& dual_averaging_state,
    arma::mat& sqrt_inv_fisher,
    const double threshold_alpha,
    const double threshold_beta,
    const double initial_step_size
) {
  // --- Compute current parameter vector and its gradient ---
  const arma::vec current_state = vectorize_thresholds (
    main_effects, num_categories, is_ordinal_variable
  );

  const arma::vec current_grad = gradient_log_pseudoposterior_thresholds (
    main_effects, residual_matrix, num_categories, num_obs_categories,
    sufficient_blume_capel, reference_category, is_ordinal_variable,
    threshold_alpha, threshold_beta
  );

  const double log_post = log_pseudoposterior_thresholds (
    main_effects, residual_matrix, num_categories, num_obs_categories,
    sufficient_blume_capel, reference_category, is_ordinal_variable,
    threshold_alpha, threshold_beta
  );

  const int dim = current_state.n_elem;

  // --- If burn-in just ended, initialize Fisher matrix using gradient
  if (iteration == total_burnin) {
    sqrt_inv_fisher = initialize_fisher_preconditioner (current_grad);
  }

  // --- Set inverse Fisher matrix: identity during warm-up, adapted afterward
  arma::mat inv_fisher;
  if (iteration >= total_burnin) {
    inv_fisher = sqrt_inv_fisher * sqrt_inv_fisher.t();
  } else {
    inv_fisher.eye(dim, dim);
  }

  // --- Construct Fisher-preconditioned MALA proposal ---
  const double trace_inv_fisher = arma::trace(inv_fisher);
  const double scaled_step_size = step_size / (trace_inv_fisher / dim);
  const double sqrt_step = std::sqrt(scaled_step_size);

  // Drift (mean shift) and stochastic noise
  const arma::vec proposal_drift = 0.5 * scaled_step_size * inv_fisher * current_grad;
  const arma::vec noise = sqrt_inv_fisher * arma::randn(dim);

  const arma::vec proposed_state = current_state + proposal_drift + sqrt_step * noise;

  // --- Map vector back to main_effects matrix
  const arma::mat proposed_main_effects = unvectorize_thresholds (
    proposed_state, num_categories, is_ordinal_variable
  );
  // --- Evaluate proposed log posterior and gradient
  const arma::vec proposed_grad = gradient_log_pseudoposterior_thresholds (
    proposed_main_effects, residual_matrix, num_categories, num_obs_categories,
    sufficient_blume_capel, reference_category, is_ordinal_variable,
    threshold_alpha, threshold_beta
  );

  const double log_post_prop = log_pseudoposterior_thresholds (
    proposed_main_effects, residual_matrix, num_categories, num_obs_categories,
    sufficient_blume_capel, reference_category, is_ordinal_variable,
    threshold_alpha, threshold_beta
  );

  // --- Compute MALA acceptance correction (Titsias 2023, Prop. 1) ---
  const arma::vec forward_proposal_residual =
    proposed_state - current_state -
    0.25 * scaled_step_size * inv_fisher * current_grad;
  const arma::vec reverse_proposal_residual =
    current_state - proposed_state -
    0.25 * scaled_step_size * inv_fisher * proposed_grad;

  const double log_forward = 0.5 * arma::dot (forward_proposal_residual, current_grad);
  const double log_backward = 0.5 * arma::dot (reverse_proposal_residual, proposed_grad);

  const double log_accept = log_post_prop - log_post + (log_backward - log_forward);
  const double accept_prob = std::min(1.0, std::exp (log_accept));

  // --- Accept or reject proposed move ---
  if (std::log (R::unif_rand()) < log_accept) {
    main_effects = proposed_main_effects;
  }

  // --- Update step size and Fisher matrix ---
  if (iteration < total_burnin) {
    // During warm-up: dual averaging adaptation
    update_step_size_with_dual_averaging (
      initial_step_size, accept_prob, iteration + 1, dual_averaging_state
    );
    step_size = std::exp (dual_averaging_state[1]);
  } else {
    // After warm-up: Robbins-Monro + Fisher preconditioner update
    update_step_size_with_robbins_monro (
      accept_prob, iteration - total_burnin + 1, step_size
    );

    const arma::vec score_diff =
      std::sqrt(accept_prob) * (proposed_grad - current_grad);

    update_fisher_preconditioner(sqrt_inv_fisher, score_diff);
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
void update_regular_thresholds_with_metropolis (
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

  for (int category = 0; category < num_cats; category++) {
    double current = main_effects(variable, category);
    double exp_current = std::exp (current);
    double c = (threshold_alpha + threshold_beta) / (1.0 + exp_current);

    for (int person = 0; person < no_persons; person++) {
      double rest_score = residual_matrix(person, variable);
      double denom = 1.0;
      double numer = std::exp ((category + 1) * rest_score);

      for (int cat = 0; cat < num_cats; cat++) {
        if (cat != category) {
          denom += std::exp (main_effects(variable, cat) + (cat + 1) * rest_score);
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
    double proposed = std::log (tmp / (1.0 - tmp) / c);
    double exp_proposed = std::exp (proposed);

    // Compute MH acceptance probability
    double log_acceptance_probability = 0.0;
    for (int person = 0; person < no_persons; person++) {
      log_acceptance_probability += std::log (g(person) + q(person) * exp_current);
      log_acceptance_probability -= std::log (g(person) + q(person) * exp_proposed);
    }

    log_acceptance_probability -= (threshold_alpha + threshold_beta) * std::log1p(exp_proposed);
    log_acceptance_probability += (threshold_alpha + threshold_beta) * std::log1p(exp_current);
    log_acceptance_probability -= (a + b) * std::log1p(c * exp_current);
    log_acceptance_probability += (a + b) * std::log1p(c * exp_proposed);

    if (std::log (R::unif_rand()) < log_acceptance_probability) {
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
void update_blumecapel_thresholds_with_adaptive_metropolis (
    arma::mat& main_effects,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& sufficient_blume_capel,
    const int num_persons,
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

  // --- Define helper for prior contribution
  auto log_beta_prior_diff = [&](double curr, double prop) {
    return (threshold_alpha + threshold_beta) *
      (std::log1p(std::exp (curr)) - std::log1p(std::exp (prop)));
  };

  // --- Update each threshold parameter: 0 = linear, 1 = quadratic
  for (int param = 0; param < 2; param++) {
    double& proposal_sd = proposal_sd_blumecapel(variable, param);
    double current = main_effects(variable, param);
    double proposed = R::rnorm(current, proposal_sd);
    double diff = proposed - current;

    arma::vec numer_current(num_cats + 1);
    arma::vec numer_proposed(num_cats + 1);

    // --- Step 1: Construct numerators for softmax (for all categories)
    for (int cat = 0; cat <= num_cats; cat++) {
      int centered = cat - ref;
      if (param == 0) {
        // Linear update
        double quad = main_effects(variable, 1) * centered * centered;
        numer_current(cat) = current * cat + quad;
        numer_proposed(cat) = proposed * cat + quad;
      } else {
        // Quadratic update
        double lin = main_effects(variable, 0) * cat;
        numer_current(cat) = current * centered * centered + lin;
        numer_proposed(cat) = proposed * centered * centered + lin;
      }
    }

    // --- Step 2: Compute lbound for numerical stability
    double max_curr = numer_current.max();
    double max_prop = numer_proposed.max();
    double lbound = (max_curr > 0.0 || max_prop > 0.0) ? std::max(max_curr, max_prop) : 0.0;

    // --- Step 3: Likelihood ratio
    // Accumulate log acceptance probability based on change in pseudo-likelihood

    // Contribution from sufficient statistics and prior
    double log_accept = diff * (threshold_alpha + sufficient_blume_capel(param, variable));

    // Vectorized likelihood ratio for all persons:
    //
    // For each person, compute:
    //   log p(y_i | proposed) - log p(y_i | current)
    //   using softmax-style normalization for categorical probabilities.
    //
    // The bound stabilizes exponentials across categories and persons.
    arma::vec rest_score = residual_matrix.col(variable);                       // Person-wise residuals
    arma::vec bound = arma::max(rest_score, arma::zeros<arma::vec>(num_persons)) * num_cats + lbound;

    arma::vec denom_curr = arma::exp (numer_current(0) - bound);                 // Score = 0 contribution
    arma::vec denom_prop = arma::exp (numer_proposed(0) - bound);

    for (int cat = 0; cat < num_cats; cat++) {
      arma::vec score_term = (cat + 1) * rest_score - bound;

      // Compute exponentials for each category and add to denominator
      denom_curr += arma::exp (numer_current(cat + 1) + score_term);
      denom_prop += arma::exp (numer_proposed(cat + 1) + score_term);
    }

    // Accumulate the person-wise log ratio contributions
    log_accept += arma::accu (arma::log (denom_curr) - arma::log (denom_prop));

    // --- Step 4: Add prior ratio
    log_accept += log_beta_prior_diff(current, proposed);

    // --- Step 5: Metropolis accept/reject
    if (std::log (R::unif_rand()) < log_accept) {
      main_effects(variable, param) = proposed;
    }

    // --- Step 6: Robbins-Monro proposal adaptation
    proposal_sd = update_proposal_sd_with_robbins_monro (
      proposal_sd, log_accept, exp_neg_log_t_rm_adaptation_rate
    );
  }
}



/**
 * Function: gradient_log_pseudoposterior_interactions
 *
 * Computes the gradient of the log pseudoposterior with respect to the
 * active interaction parameters. This is used in MALA updates of the
 * pairwise interaction matrix.
 *
 * Inputs:
 *  - pairwise_effects: Symmetric matrix of interaction parameters [V × V].
 *  - main_effects: Matrix of main effect (threshold) parameters [V × max_categories].
 *  - observations: Matrix of ordinal and Blume-Capel scores (row = individual, col = variable).
 *  - num_categories: Vector of number of categories per variable.
 *  - inclusion_indicator: Symmetric binary matrix indicating which interactions are active.
 *  - is_ordinal_variable: Logical vector indicating ordinal (1) or Blume-Capel (0) variables.
 *  - reference_category: Vector of reference categories for BC variables.
 *  - interaction_scale: Cauchy prior scale on interaction weights.
 *
 * Returns:
 *  - Gradient vector for all pairwise interactions (in upper-triangle order).
 */
arma::vec gradient_log_pseudoposterior_interactions (
    const arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    const double interaction_scale
) {
  const int num_variables = observations.n_cols;
  const int num_observations = observations.n_rows;
  const int num_interactions = (num_variables * (num_variables - 1)) / 2;

  arma::umat index_matrix(num_variables, num_variables);
  int counter = 0;
  for(int var1 = 0; var1 < num_variables-1; var1++) {
    for(int var2 = var1 + 1; var2 < num_variables; var2++) {
      index_matrix(var1, var2) = counter;
      counter++;
    }
  }

  arma::vec gradient (num_interactions, arma::fill::zeros);

  for (int var = 0; var < num_variables; var++) {
    int num_cats = num_categories (var);
    arma::mat score_weights (num_observations, num_cats, arma::fill::zeros);    //First column would be zero

    arma::vec rest_scores = observations * pairwise_effects.col (var);
    arma::vec denominator = arma::zeros (num_observations);
    arma::vec bounds = arma::max (rest_scores, arma::zeros<arma::vec> (num_observations)) * num_cats;

    if (is_ordinal_variable (var)) {
      denominator += arma::exp (-bounds);
      for (int category = 0; category < num_cats; category++) {
        arma::vec exponent = main_effects (var, category) + (category + 1) * rest_scores - bounds;
        arma::vec weight = arma::exp (exponent);
        denominator += weight;
        score_weights.col(category) = (category + 1) * weight;
      }
    } else {
      const int ref_cat = reference_category (var);
      // Zero category
      double quad_term = main_effects (var, 1) * ref_cat * ref_cat;
      arma::vec exponent = quad_term - bounds;
      denominator = arma::exp (exponent);
      for (int category = 1; category <= num_cats; category++) {
        int centered_cat = category - ref_cat;
        double lin_term = main_effects (var, 0) * category;
        double quad_term = main_effects (var, 1) * centered_cat * centered_cat;
        arma::vec exponent = lin_term + quad_term + category * rest_scores - bounds;
        arma::vec weight = arma::exp (exponent);
        denominator += weight;
        score_weights.col(category - 1) = category * weight;
      }
    }
    score_weights.each_col() /= denominator;

    for(int var2 = 0; var2 < num_variables; var2++) {
      if (inclusion_indicator (var, var2) == 0)
        continue;

      arma::vec expected_value(num_observations, arma::fill::zeros);
      if(var == var2)
        continue;

      arma::ivec xv = observations.col(var2);

      for (int category = 0; category < num_cats; category++) {
        expected_value += score_weights.col(category) % xv;
      }

      int location = (var < var2) ? index_matrix(var, var2) : index_matrix(var2, var);

      gradient(location) -= arma::accu(expected_value);
    }
  }

  for(int var1 = 0; var1 < num_variables; var1++) {
    for (int var2 = var1 + 1; var2 < num_variables; var2++) {
      if (inclusion_indicator (var1, var2) == 0)
        continue;

      int location = index_matrix(var1, var2);

      gradient (location) += 2.0 * arma::dot (observations.col(var1), observations.col(var2));

      // ---- Gradient contribution from Cauchy prior
      const double effect = pairwise_effects (var1, var2);
      gradient (location) -= 2.0 * effect / (effect * effect + interaction_scale * interaction_scale);
    }
  }

  return gradient;
}


/**
 * Function: log_pseudoposterior_interactions
 *
 * Computes the log of the pseudoposterior distribution over interaction parameters.
 * Includes:
 *  - The pairwise quadratic sufficient statistic: trace(X * B * Xᵗ)
 *  - The (unnormalized) nodewise log-likelihood contributions
 *  - The log Cauchy prior for included interaction terms
 *
 * Inputs:
 *  - pairwise_effects: Symmetric matrix of interaction parameters.
 *  - main_effects: Matrix of main effect (threshold) parameters.
 *  - observations: Matrix of ordinal and BC scores.
 *  - num_categories: Vector of number of categories for each variable.
 *  - inclusion_indicator: Binary matrix indicating active interactions.
 *  - is_ordinal_variable: Logical vector: 1 = ordinal, 0 = Blume-Capel.
 *  - reference_category: Reference category index per variable (for BC).
 *  - interaction_scale: Cauchy prior scale for interaction terms.
 *
 * Returns:
 *  - Scalar log pseudoposterior (double)
 */
double log_pseudoposterior_interactions (
    const arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    const double interaction_scale
) {
  const int num_variables = observations.n_cols;
  const int num_observations = observations.n_rows;

  // Convert to double matrix for trace calculation
  arma::mat real_observations = arma::conv_to<arma::mat>::from (observations);

  // Leading term: trace(X * B * X^T)
  double log_pseudo_likelihood = arma::trace (real_observations * pairwise_effects * real_observations.t ());

  for (int var = 0; var < num_variables; var++) {
    int num_categories_var = num_categories (var);

    // Compute rest score: contribution from other variables
    arma::vec rest_scores = observations * pairwise_effects.col (var);
    arma::vec bounds = arma::max (rest_scores, arma::zeros<arma::vec> (num_observations)) * num_categories_var;
    arma::vec denominator = arma::zeros (num_observations);

    if (is_ordinal_variable (var)) {
      // Ordinal variable: denominator includes exp (-bounds) + exp over categories
      denominator += arma::exp ( -bounds );
      for (int category = 0; category < num_categories_var; category++) {
        arma::vec exponent = main_effects (var, category) + (category + 1) * rest_scores - bounds;
        denominator += arma::exp (exponent);
      }
    } else {
      // Binary/categorical variable: quadratic + linear term
      const int ref_cat = reference_category (var);
      for (int category = 0; category <= num_categories_var; category++) {
        int centered_cat = category - ref_cat;
        double lin_term = main_effects (var, 0) * category;
        double quad_term = main_effects (var, 1) * centered_cat * centered_cat;
        arma::vec exponent = lin_term + quad_term + category * rest_scores - bounds;
        denominator += arma::exp (exponent);
      }
    }

    // Subtract log partition function and bounds adjustment
    log_pseudo_likelihood -= arma::accu (arma::log (denominator));
    log_pseudo_likelihood -= arma::accu (bounds);
  }

  // Add Cauchy prior terms for included pairwise effects
  for (int var1 = 0; var1 < num_variables - 1; var1++) {
    for (int var2 = var1 + 1; var2 < num_variables; var2++) {
      if (inclusion_indicator (var1, var2) == 1) {
        log_pseudo_likelihood += R::dcauchy (pairwise_effects (var1, var2), 0.0, interaction_scale, true);
      }
    }
  }

  return log_pseudo_likelihood;
}



/**
 * Function: find_reasonable_initial_step_size_mala_interactions
 *
 * Heuristically selects an initial step size for MALA updates of the interaction parameters
 * by adapting the acceptance rate toward a target value (~0.574).
 *
 * Inputs:
 *  - pairwise_effects: Matrix of pairwise interaction parameters.
 *  - main_effects: Matrix of main effect (threshold) parameters.
 *  - observations: Matrix of categorical scores.
 *  - num_categories: Number of categories per variable.
 *  - inclusion_indicator: Symmetric binary matrix of active interactions.
 *  - is_ordinal_variable: Logical vector (1 = ordinal, 0 = Blume-Capel).
 *  - reference_category: Vector of reference categories (for BC variables).
 *  - interaction_scale: Prior scale for Cauchy prior on interactions.
 *
 * Returns:
 *  - A positive step size value that gives roughly 0.574 acceptance rate.
 */
double find_reasonable_initial_step_size_mala_interactions (
    arma::mat pairwise_effects,
    const arma::mat& main_effects,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    const double interaction_scale
) {
  const double target_acceptance = 0.574;
  const double initial_step_size = 0.1;
  const int max_attempts = 20;

  const int num_variables = pairwise_effects.n_rows;
  const int num_interactions = (num_variables * (num_variables - 1)) / 2;

  double log_step_size = std::log (initial_step_size);

  // --- Step 1: Extract current interaction state as vector
  arma::vec current_state (num_interactions, arma::fill::zeros);
  int interaction_index = -1;
  for (int var1 = 0; var1 < num_variables - 1; var1++) {
    for (int var2 = var1 + 1; var2 < num_variables; var2++) {
      interaction_index++;
      if (inclusion_indicator (var1, var2) == 1) {
        current_state (interaction_index) = pairwise_effects (var1, var2);
      }
    }
  }

  // --- Step 2: Evaluate gradient and posterior at current state
  arma::vec gradient = gradient_log_pseudoposterior_interactions (
    pairwise_effects, main_effects, observations, num_categories,
    inclusion_indicator, is_ordinal_variable, reference_category,
    interaction_scale
  );

  double log_post_current = log_pseudoposterior_interactions (
    pairwise_effects, main_effects, observations, num_categories,
    inclusion_indicator, is_ordinal_variable, reference_category,
    interaction_scale
  );

  int direction = 0;
  double accept_prob = 0.0;

  // --- Step 3: Exponential step-size search loop
  for (int attempt = 0; attempt < max_attempts; attempt++) {
    double step_size = std::exp (log_step_size);

    // Generate Langevin proposal
    arma::vec noise = arma::randn<arma::vec> (num_interactions);
    arma::vec proposal = current_state +
      0.5 * step_size * gradient + std::sqrt (step_size) * noise;

    // Rebuild symmetric proposal matrix
    arma::mat proposal_matrix = pairwise_effects;
    interaction_index = -1;
    for (int var1 = 0; var1 < num_variables - 1; var1++) {
      for (int var2 = var1 + 1; var2 < num_variables; var2++) {
        interaction_index++;
        if (inclusion_indicator (var1, var2) == 1) {
          proposal_matrix (var1, var2) = proposal (interaction_index);
          proposal_matrix (var2, var1) = proposal (interaction_index);
        }
      }
    }

    // Evaluate posterior and gradient at proposed state
    double log_post_proposal = log_pseudoposterior_interactions (
      proposal_matrix, main_effects, observations, num_categories,
      inclusion_indicator, is_ordinal_variable, reference_category,
      interaction_scale
    );

    arma::vec gradient_prop = gradient_log_pseudoposterior_interactions (
      proposal_matrix, main_effects, observations, num_categories,
      inclusion_indicator, is_ordinal_variable, reference_category,
      interaction_scale
    );

    // Compute log proposal densities
    arma::vec forward = current_state + 0.5 * step_size * gradient;
    arma::vec backward = proposal + 0.5 * step_size * gradient_prop;

    double log_q_forward = -0.5 / step_size * arma::accu (arma::square (proposal - forward));
    double log_q_reverse = -0.5 / step_size * arma::accu (arma::square (current_state - backward));

    double log_accept = log_post_proposal + log_q_reverse - log_post_current - log_q_forward;
    accept_prob = std::min (1.0, std::exp (log_accept));

    // --- Step 4: Decide direction based on first attempt
    if (attempt == 0) {
      direction = (accept_prob > target_acceptance) ? 1 : -1;
    } else {
      if ((direction == 1 && accept_prob < target_acceptance) ||
          (direction == -1 && accept_prob > target_acceptance)) {
        break;
      }
    }

    // Adjust log step size in chosen direction
    log_step_size += direction;
  }

  return std::exp (log_step_size);
}




/**
 * Function: update_interactions_with_mala
 *
 * Performs a blockwise MALA update of the interaction parameters using a fixed
 * or adaptive step size. Uses full-length vectors over all pairwise interactions,
 * with inactive elements set to zero.
 *
 * Inputs:
 *  - pairwise_effects: Symmetric matrix of interaction parameters (updated in-place).
 *  - residual_matrix: Residual matrix (updated in-place if proposal accepted).
 *  - main_effects: Matrix of main effect (threshold) parameters.
 *  - observations: Matrix of observed scores.
 *  - num_categories: Vector of number of categories per variable.
 *  - inclusion_indicator: Binary matrix indicating active interactions.
 *  - is_ordinal_variable: Logical vector: 1 = ordinal, 0 = Blume-Capel.
 *  - reference_category: Vector of reference categories for BC variables.
 *  - interaction_scale: Scale parameter for Cauchy prior on interactions.
 *  - step_size_interactions: Current step size (updated if adaptive).
 *  - initial_step_size_interactions: Initial step size used during dual averaging.
 *  - iteration: Current MCMC iteration (0-based).
 *  - total_burnin: Total number of burn-in iterations.
 *  - dual_averaging_state: Vector (length 3) tracking dual averaging state.
 *
 * Modifies:
 *  - pairwise_effects
 *  - residual_matrix
 *  - step_size_interactions (if in burn-in or adaptation phase)
 *  - dual_averaging_state (during burn-in)
 */
void update_interactions_with_mala (
    arma::mat& pairwise_effects,
    arma::mat& residual_matrix,
    const arma::mat& main_effects,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    const double interaction_scale,
    double& step_size_interactions,
    const double initial_step_size_interactions,
    const int iteration,
    const int total_burnin,
    arma::vec& dual_averaging_state
) {
  const int num_variables = pairwise_effects.n_rows;
  const int num_interactions = (num_variables * (num_variables - 1)) / 2;

  // --- Flatten current interaction matrix to vector
  arma::vec current_state (num_interactions, arma::fill::zeros);
  int interaction_index = -1;
  for (int var1 = 0; var1 < num_variables - 1; var1++) {
    for (int var2 = var1 + 1; var2 < num_variables; var2++) {
      interaction_index++;
      if (inclusion_indicator (var1, var2) == 1) {
        current_state (interaction_index) = pairwise_effects (var1, var2);
      }
    }
  }

  // --- Compute gradient and posterior at current state
  arma::vec grad_current = gradient_log_pseudoposterior_interactions (
    pairwise_effects, main_effects, observations, num_categories,
    inclusion_indicator, is_ordinal_variable, reference_category, interaction_scale
  );

  double log_post_current = log_pseudoposterior_interactions (
    pairwise_effects, main_effects, observations, num_categories,
    inclusion_indicator, is_ordinal_variable, reference_category, interaction_scale
  );

  // --- Generate Langevin proposal
  arma::vec noise = arma::randn<arma::vec> (num_interactions);
  arma::vec proposal = current_state +
    0.5 * step_size_interactions * grad_current +
    std::sqrt (step_size_interactions) * noise;

  // --- Build symmetric proposal matrix
  arma::mat proposal_matrix = pairwise_effects;
  interaction_index = -1;
  for (int var1 = 0; var1 < num_variables - 1; var1++) {
    for (int var2 = var1 + 1; var2 < num_variables; var2++) {
      interaction_index++;
      if (inclusion_indicator (var1, var2) == 1) {
        proposal_matrix (var1, var2) = proposal (interaction_index);
        proposal_matrix (var2, var1) = proposal (interaction_index);
      }
    }
  }

  // --- Evaluate posterior and gradient at proposed state
  double log_post_proposal = log_pseudoposterior_interactions (
    proposal_matrix, main_effects, observations, num_categories,
    inclusion_indicator, is_ordinal_variable, reference_category, interaction_scale
  );

  arma::vec grad_proposed = gradient_log_pseudoposterior_interactions (
    proposal_matrix, main_effects, observations, num_categories,
    inclusion_indicator, is_ordinal_variable, reference_category, interaction_scale
  );

  // --- Compute MH log acceptance ratio
  arma::vec forward_mean = current_state + 0.5 * step_size_interactions * grad_current;
  arma::vec reverse_mean = proposal + 0.5 * step_size_interactions * grad_proposed;

  double log_q_forward = -0.5 / step_size_interactions * arma::accu (arma::square (proposal - forward_mean));
  double log_q_reverse = -0.5 / step_size_interactions * arma::accu (arma::square (current_state - reverse_mean));

  double log_accept = log_post_proposal + log_q_reverse - log_post_current - log_q_forward;

  double logu = std::log (R::runif (0.0, 1.0));
  bool accepted = (logu < log_accept);
  double acceptance_prob = std::min (1.0, std::exp (log_accept));

  // --- Adapt step size
  if (iteration < total_burnin) {
    update_step_size_with_dual_averaging (
        initial_step_size_interactions,
        acceptance_prob,
        iteration + 1,
        dual_averaging_state
    );
    step_size_interactions = std::exp (dual_averaging_state (1));
  } else {
    update_step_size_with_robbins_monro (
        acceptance_prob,
        iteration - total_burnin + 1,
        step_size_interactions
    );
  }

  // --- Accept and update interaction matrix and residuals
  if (accepted) {
    interaction_index = -1;
    for (int var1 = 0; var1 < num_variables - 1; var1++) {
      for (int var2 = var1 + 1; var2 < num_variables; var2++) {
        interaction_index++;
        if (inclusion_indicator (var1, var2) == 1) {
          double delta = proposal (interaction_index) - pairwise_effects (var1, var2);
          pairwise_effects (var1, var2) = proposal (interaction_index);
          pairwise_effects (var2, var1) = proposal (interaction_index);
          residual_matrix.col (var1) += arma::conv_to<arma::vec>::from (observations.col (var2)) * delta;
          residual_matrix.col (var2) += arma::conv_to<arma::vec>::from (observations.col (var1)) * delta;
        }
      }
    }
  }
}



/**
 * Function: compute_log_likelihood_ratio_for_variable
 *
 * Computes the log pseudo-likelihood ratio contribution for a single variable,
 * comparing a proposed vs. current interaction value. This is used to evaluate
 * Metropolis-Hastings updates to a pairwise interaction parameter.
 *
 * The function is vectorized over persons and supports both ordinal and
 * Blume-Capel variables.
 *
 * Inputs:
 *  - variable: Index of the variable whose likelihood contribution is evaluated.
 *  - interacting_score: Vector of category scores for the interacting variable (one per person).
 *  - proposed_state: Proposed interaction value.
 *  - current_state: Current interaction value.
 *  - main_effects: Matrix of threshold parameters [variables × categories].
 *  - num_categories: Vector with number of categories per variable.
 *  - residual_matrix: Current matrix of residual predictors (one column per variable).
 *  - observations: Data matrix of categorical scores (only used for row count).
 *  - is_ordinal_variable: Logical vector (1 if ordinal, 0 if Blume-Capel).
 *  - reference_category: Reference category per variable (for BC variables).
 *
 * Returns:
 *  - The total log pseudo-likelihood ratio for the given variable, summed over all persons.
 */
double compute_log_likelihood_ratio_for_variable (
    int variable,
    const arma::ivec& interacting_score,
    double proposed_state,
    double current_state,
    const arma::mat& main_effects,
    const arma::ivec& num_categories,
    const arma::mat& residual_matrix,
    const arma::imat& observations,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category
) {
  // Convert interaction score vector to double precision
  arma::vec interaction = arma::conv_to<arma::vec>::from (interacting_score);

  const int num_persons = residual_matrix.n_rows;
  const int num_categories_var = num_categories (variable);

  // Compute adjusted linear predictors without the current interaction
  arma::vec rest_scores = residual_matrix.col (variable) - interaction * current_state;

  // Stability bound for softmax (scaled by number of categories)
  arma::vec bounds = arma::max (rest_scores, arma::zeros<arma::vec> (num_persons)) * num_categories_var;

  arma::vec denom_current = arma::zeros (num_persons);
  arma::vec denom_proposed = arma::zeros (num_persons);

  if (is_ordinal_variable (variable)) {
    // Ordinal model: initialize with underflow-protection constant
    denom_current += arma::exp ( -bounds );
    denom_proposed += arma::exp ( -bounds );

    for (int category = 0; category < num_categories_var; category++) {
      arma::vec exponent = main_effects (variable, category) + (category + 1) * rest_scores;
      denom_current += arma::exp (exponent + (category + 1) * interaction * current_state - bounds);
      denom_proposed += arma::exp (exponent + (category + 1) * interaction * proposed_state - bounds);
    }

  } else {
    // Binary or categorical variable: linear + quadratic score
    const int ref_cat = reference_category (variable);

    for (int category = 0; category <= num_categories_var; category++) {
      int centered = category - ref_cat;
      double lin_term = main_effects (variable, 0) * category;
      double quad_term = main_effects (variable, 1) * centered * centered;
      arma::vec exponent = lin_term + quad_term + category * rest_scores - bounds;

      denom_current += arma::exp (exponent + category * interaction * current_state);
      denom_proposed += arma::exp (exponent + category * interaction * proposed_state);
    }
  }

  // Accumulated log-likelihood difference across persons
  return arma::accu (arma::log (denom_current) - arma::log (denom_proposed));
}



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
 *  - no_persons: Number of individuals (rows in the data).
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

  // Extract score vectors for both variables across all persons
  arma::ivec score1 = observations.col(variable1);
  arma::ivec score2 = observations.col(variable2);

  // (1) Direct interaction contribution to the linear predictor:
  //     Δβ × ∑(score1_i × score2_i) for all persons i
  log_ratio += 2.0 * arma::dot (score1, score2) * delta;

  // (2) Change in pseudo-likelihood for variable1 due to the update in its interaction with variable2
  log_ratio += compute_log_likelihood_ratio_for_variable (
    variable1, score2, proposed_state, current_state, main_effects,
    num_categories, residual_matrix, observations, is_ordinal_variable,
    reference_category
  );

  // (3) Symmetric change for variable2 due to its interaction with variable1
  log_ratio += compute_log_likelihood_ratio_for_variable (
    variable2, score1, proposed_state, current_state, main_effects,
    num_categories, residual_matrix, observations, is_ordinal_variable,
    reference_category
  );

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
void update_interactions_with_adaptive_metropolis (
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
  for (int variable1 = 0; variable1 < num_variables - 1; variable1++) {
    for (int variable2 = variable1 + 1; variable2 < num_variables; variable2++) {
      if (inclusion_indicator(variable1, variable2) == 1) {
        // Sample proposal from Gaussian centered at current state
        const double current_state = pairwise_effects(variable1, variable2);
        const double proposed_state = R::rnorm(current_state, proposal_sd_pairwise_effects(variable1, variable2));

        // Compute log acceptance ratio: data + prior
        double log_acceptance = log_pseudolikelihood_ratio_interaction(
          pairwise_effects, main_effects, observations, num_categories, num_persons,
          variable1, variable2, proposed_state, current_state,
          residual_matrix, is_ordinal_variable, reference_category
        );

        // Add symmetric Cauchy prior ratio
        log_acceptance += R::dcauchy(proposed_state, 0.0, interaction_scale, true);
        log_acceptance -= R::dcauchy(current_state, 0.0, interaction_scale, true);

        // Accept proposal with MH step
        if (std::log (R::unif_rand()) < log_acceptance) {
          const double delta = proposed_state - current_state;

          // Update effect matrix symmetrically
          pairwise_effects(variable1, variable2) = proposed_state;
          pairwise_effects(variable2, variable1) = proposed_state;

          // Vectorized update of residual matrix for both variables
          residual_matrix.col(variable1) += arma::conv_to<arma::vec>::from(observations.col(variable2)) * delta;
          residual_matrix.col(variable2) += arma::conv_to<arma::vec>::from(observations.col(variable1)) * delta;
        }

        // Robbins-Monro adaptation of proposal SD
        proposal_sd_pairwise_effects(variable1, variable2) =
          update_proposal_sd_with_robbins_monro (
            proposal_sd_pairwise_effects(variable1, variable2), log_acceptance,
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
void update_indicator_interaction_pair_with_metropolis (
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
  for (int cntr = 0; cntr < no_interactions; cntr++) {
    const int variable1 = index(cntr, 1);
    const int variable2 = index(cntr, 2);

    const double current_state = pairwise_effects(variable1, variable2);

    // Propose a new state: either add a new edge or remove an existing one
    const bool proposing_addition = (indicator(variable1, variable2) == 0);
    const double proposed_state = proposing_addition ? R::rnorm(current_state, proposal_sd(variable1, variable2)) : 0.0;

    // Compute log pseudo-likelihood ratio
    double log_accept = log_pseudolikelihood_ratio_interaction (
      pairwise_effects, main_effects, observations, num_categories, no_persons,
      variable1, variable2, proposed_state, current_state, residual_matrix,
      is_ordinal_variable, reference_category
    );

    // Add prior ratio and proposal correction
    const double theta_ij = theta(variable1, variable2);
    const double sd = proposal_sd(variable1, variable2);

    if (proposing_addition) {
      log_accept += R::dcauchy(proposed_state, 0.0, interaction_scale, true);
      log_accept -= R::dnorm(proposed_state, current_state, sd, true);
      log_accept += std::log (theta_ij) - std::log (1.0 - theta_ij);
    } else {
      log_accept -= R::dcauchy(current_state, 0.0, interaction_scale, true);
      log_accept += R::dnorm(current_state, proposed_state, sd, true);
      log_accept -= std::log (theta_ij) - std::log (1.0 - theta_ij);
    }

    // Metropolis-Hastings accept step
    if (std::log (R::unif_rand()) < log_accept) {
      const int updated_indicator = 1 - indicator(variable1, variable2);
      indicator(variable1, variable2) = updated_indicator;
      indicator(variable2, variable1) = updated_indicator;

      pairwise_effects(variable1, variable2) = proposed_state;
      pairwise_effects(variable2, variable1) = proposed_state;

      const double delta = proposed_state - current_state;

      // Vectorized residual update
      residual_matrix.col(variable1) += arma::conv_to<arma::vec>::from(observations.col(variable2)) * delta;
      residual_matrix.col(variable2) += arma::conv_to<arma::vec>::from(observations.col(variable1)) * delta;
    }
  }
}



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
) {
  const int num_persons = observations.n_rows;

  // Extract observed score vectors for each variable
  arma::vec x_var1 = arma::conv_to<arma::vec>::from (observations.col (var1));
  arma::vec x_var2 = arma::conv_to<arma::vec>::from (observations.col (var2));

  // First-order gradient from data
  double gradient = 2.0 * arma::dot (x_var1, x_var2);

  // --- Contribution from var1
  int num_categories_var1 = num_categories (var1);
  arma::vec rest_scores_var1 = observations * pairwise_effects.col (var1);  // β_{var1,var1} = 0
  arma::vec numerator_var1 (num_persons, arma::fill::zeros);
  arma::vec denominator_var1 (num_persons, arma::fill::zeros);
  arma::vec bounds_var1 = arma::max (rest_scores_var1, arma::zeros<arma::vec> (num_persons)) * num_categories_var1;

  if (is_ordinal_variable (var1)) {
    denominator_var1 += arma::exp ( -bounds_var1 );
    for (int category = 0; category < num_categories_var1; category++) {
      arma::vec exponent = main_effects (var1, category) + (category + 1) * rest_scores_var1 - bounds_var1;
      arma::vec weight = arma::exp (exponent);
      denominator_var1 += weight;
      numerator_var1 += (category + 1) * x_var2 % weight;
    }
  } else {
    const int ref_cat = reference_category (var1);
    for (int category = 0; category <= num_categories_var1; category++) {
      int centered = category - ref_cat;
      double lin_term = main_effects (var1, 0) * category;
      double quad_term = main_effects (var1, 1) * centered * centered;
      arma::vec exponent = lin_term + quad_term + category * rest_scores_var1 - bounds_var1;
      arma::vec weight = arma::exp (exponent);
      denominator_var1 += weight;
      numerator_var1 += category * x_var2 % weight;
    }
  }

  gradient -= arma::accu (numerator_var1 / denominator_var1);

  // --- Contribution from var2
  int num_categories_var2 = num_categories (var2);
  arma::vec rest_scores_var2 = observations * pairwise_effects.col (var2);
  arma::vec numerator_var2 (num_persons, arma::fill::zeros);
  arma::vec denominator_var2 (num_persons, arma::fill::zeros);
  arma::vec bounds_var2 = arma::max (rest_scores_var2, arma::zeros<arma::vec> (num_persons)) * num_categories_var2;

  if (is_ordinal_variable (var2)) {
    denominator_var2 += arma::exp ( -bounds_var2 );
    for (int category = 0; category < num_categories_var2; category++) {
      arma::vec exponent = main_effects (var2, category) + (category + 1) * rest_scores_var2 - bounds_var2;
      arma::vec weight = arma::exp (exponent);
      denominator_var2 += weight;
      numerator_var2 += (category + 1) * x_var1 % weight;
    }
  } else {
    const int ref_cat = reference_category (var2);
    for (int category = 0; category <= num_categories_var2; category++) {
      int centered = category - ref_cat;
      double lin_term = main_effects (var2, 0) * category;
      double quad_term = main_effects (var2, 1) * centered * centered;
      arma::vec exponent = lin_term + quad_term + category * rest_scores_var2 - bounds_var2;
      arma::vec weight = arma::exp (exponent);
      denominator_var2 += weight;
      numerator_var2 += category * x_var1 % weight;
    }
  }

  gradient -= arma::accu (numerator_var2 / denominator_var2);

  // --- Cauchy prior derivative
  double beta = pairwise_effects (var1, var2);
  gradient -= 2.0 * beta / (beta * beta + interaction_scale * interaction_scale);

  return gradient;
}



/**
 * Function: update_indicator_interaction_pair_with_mala
 *
 * Proposes and accepts/rejects inclusion of pairwise interaction terms in a sparse graphical model
 * using a Metropolis-adjusted Langevin algorithm (MALA). Operates one pair at a time across a list
 * of possible interactions.
 *
 * Inputs:
 *  - pairwise_effects: Matrix of pairwise interaction effects (updated in-place).
 *  - main_effects: Matrix of main effect (threshold) parameters.
 *  - indicator: Binary matrix indicating inclusion of pairwise interaction terms (updated in-place).
 *  - observations: Matrix of observed scores.
 *  - num_categories: Vector with number of categories per variable.
 *  - step_size_mala_pairwise: Step size used for MALA proposals.
 *  - interaction_scale: Scale parameter of the Cauchy prior on interaction effects.
 *  - index: Matrix listing all candidate pairs for interaction updates.
 *  - no_interactions: Number of candidate interaction pairs to consider.
 *  - num_persons: Number of observations (individuals).
 *  - residual_matrix: Matrix of residuals (updated if proposal accepted).
 *  - theta: Matrix of posterior inclusion probabilities for pairwise interactions.
 *  - is_ordinal_variable: Indicator vector specifying which variables are ordinal.
 *  - reference_category: Vector of reference categories for binary/categorical variables.
 *
 * Modifies:
 *  - pairwise_effects
 *  - indicator
 *  - residual_matrix (only if proposal accepted)
 */
void update_indicator_interaction_pair_with_mala (
    arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    arma::imat& indicator,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const double step_size_mala_pairwise,
    const double interaction_scale,
    const arma::imat& index,
    const int no_interactions,
    const int num_persons,
    arma::mat& residual_matrix,
    const arma::mat& theta,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category
) {
  for (int cntr = 0; cntr < no_interactions; cntr++) {
    const int variable1 = index(cntr, 1);
    const int variable2 = index(cntr, 2);

    // Determine if we are proposing to add (if currently absent)
    const bool proposing_addition = (indicator(variable1, variable2) == 0);
    const double current_state = pairwise_effects(variable1, variable2);
    double proposed_state = 0.0;
    double log_accept = 0.0;
    const double theta_ij = theta(variable1, variable2);

    if (proposing_addition) {
      // Compute gradient of log-pseudo-posterior for interaction term
      double grad = gradient_log_pseudoposterior_interaction_single (
        variable1, variable2, pairwise_effects, main_effects, observations,
        num_categories, is_ordinal_variable, reference_category,
        interaction_scale
      );

      // MALA proposal: Langevin step forward
      double sd = std::sqrt(step_size_mala_pairwise);
      double noise = R::rnorm(0.0, sd);
      double forward_mean = current_state + 0.5 * step_size_mala_pairwise * grad;
      proposed_state = forward_mean + noise;
      log_accept -= R::dnorm(proposed_state, forward_mean, sd, true);

      // Cauchy prior on interaction effect
      log_accept += R::dcauchy(proposed_state, 0.0, interaction_scale, true);

      // Prior inclusion probability
      log_accept += std::log (theta_ij) - std::log (1.0 - theta_ij);
    } else {

      // MALA proposal: Langevin step backward
      arma::mat proposed_matrix = pairwise_effects;
      proposed_matrix(variable1, variable2) = proposed_state;
      proposed_matrix(variable2, variable1) = proposed_state;
      double grad = gradient_log_pseudoposterior_interaction_single (
        variable1, variable2, proposed_matrix, main_effects, observations,
        num_categories, is_ordinal_variable, reference_category,
        interaction_scale
      );
      double sd = std::sqrt(step_size_mala_pairwise);
      double backward_mean = proposed_state + 0.5 * step_size_mala_pairwise * grad;
      log_accept += R::dnorm(current_state, backward_mean, sd, true);

      // Cauchy prior on interaction effect
      log_accept -= R::dcauchy(current_state, 0.0, interaction_scale, true);
      // Prior inclusion probability
      log_accept -= std::log (theta_ij) - std::log (1.0 - theta_ij);
    }

    log_accept += log_pseudolikelihood_ratio_interaction (
      pairwise_effects, main_effects, observations, num_categories, num_persons,
      variable1, variable2, proposed_state, current_state, residual_matrix,
      is_ordinal_variable, reference_category
    );

    // Metropolis-Hastings accept step
    if (std::log (R::unif_rand()) < log_accept) {
      const int new_value = 1 - indicator(variable1, variable2);
      indicator(variable1, variable2) = new_value;
      indicator(variable2, variable1) = new_value;

      const double delta = proposed_state - current_state;
      pairwise_effects(variable1, variable2) = proposed_state;
      pairwise_effects(variable2, variable1) = proposed_state;

      residual_matrix.col(variable1) += arma::conv_to<arma::vec>::from(observations.col(variable2)) * delta;
      residual_matrix.col(variable2) += arma::conv_to<arma::vec>::from(observations.col(variable1)) * delta;
    }
  }
}



/**
 * Function: gibbs_update_step_for_graphical_model_parameters
 *
 * Performs a single iteration of the Gibbs sampler for a graphical model.
 *
 * This update step includes:
 *  - Edge inclusion updates (if selection is enabled)
 *  - Pairwise interaction parameter updates via Metropolis-Hastings
 *  - Main effect (threshold) updates using MALA or Metropolis
 *
 * Proposal standard deviations are adapted using Robbins-Monro or dual averaging.
 *
 * Inputs:
 *  - observations: Matrix of observed categorical scores.
 *  - num_categories: Vector of category counts per variable.
 *  - interaction_scale: Scale for the Cauchy prior on interactions.
 *  - proposal_sd_pairwise_effects: Matrix of proposal SDs for pairwise parameters (updated in-place).
 *  - proposal_sd_main_effects: Matrix of proposal SDs for main effect parameters (updated in-place).
 *  - update_index: Edge index list to use for edge selection.
 *  - num_obs_categories: Matrix of observed category counts.
 *  - sufficient_blume_capel: Sufficient statistics for BC parameters.
 *  - threshold_alpha, threshold_beta: Hyperparameters for threshold priors.
 *  - num_persons: Number of observations.
 *  - num_variables: Number of variables.
 *  - num_pairwise: Number of candidate interaction pairs.
 *  - num_main: Number of main effect parameters.
 *  - inclusion_indicator: Symmetric matrix of active edges (updated in-place).
 *  - pairwise_effects: Matrix of interaction parameters (updated in-place).
 *  - main_effects: Matrix of main effect parameters (updated in-place).
 *  - residual_matrix: Matrix of residuals (updated in-place).
 *  - theta: Matrix of prior inclusion probabilities.
 *  - rm_decay_rate: Robbins-Monro decay rate.
 *  - is_ordinal_variable: Logical vector (1 = ordinal, 0 = Blume-Capel).
 *  - reference_category: Reference categories for BC variables.
 *  - edge_selection: Whether to update inclusion indicators this iteration.
 *  - step_size_mala: MALA step size (updated if MALA is active).
 *  - iteration: Current iteration number.
 *  - dual_averaging_state: Vector of dual averaging state [log_eps, log_eps_avg, H_bar] (updated).
 *  - total_burnin: Total number of burn-in iterations.
 *  - use_mala: Whether to use MALA for updating.
 *
 * Updates (in-place):
 *  - inclusion_indicator
 *  - pairwise_effects
 *  - main_effects
 *  - residual_matrix
 *  - step_size_mala
 *  - proposal_sd_main_effects
 *  - proposal_sd_pairwise_effects
 */
void gibbs_update_step_for_graphical_model_parameters (
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const double interaction_scale,
    arma::mat& proposal_sd_pairwise_effects,
    arma::mat& proposal_sd_blumecapel,
    const arma::imat& index,
    const arma::imat& num_obs_categories,
    const arma::imat& sufficient_blume_capel,
    const double threshold_alpha,
    const double threshold_beta,
    const int num_persons,
    const int num_variables,
    const int num_pairwise,
    const int num_main,
    arma::imat& inclusion_indicator,
    arma::mat& pairwise_effects,
    arma::mat& main_effects,
    arma::mat& residual_matrix,
    const arma::mat& theta,
    const double rm_decay_rate,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    const bool edge_selection,
    double& step_size_mala_main,
    const int iteration,
    arma::vec& dual_averaging_main,
    const int total_burnin,
    const bool use_mala,
    const double initial_step_size_mala_main,
    arma::mat& sqrt_inv_fisher,
    double& step_size_mala_pairwise,
    arma::vec& dual_averaging_pairwise,
    const double initial_step_size_mala_pairwise
) {
  // Robbins-Monro learning rate (e.g., iteration^-0.75)
  const double exp_neg_log_t_rm_adaptation_rate =
    std::exp (-std::log (static_cast<double>(iteration)) * rm_decay_rate);

  // Step 1: Edge selection via MH indicator updates (if enabled)
  if (edge_selection) {
    if(use_mala) {
      update_indicator_interaction_pair_with_mala (
        pairwise_effects, main_effects, inclusion_indicator, observations,
        num_categories, step_size_mala_pairwise, interaction_scale,
        index, num_pairwise, num_persons, residual_matrix,
        theta, is_ordinal_variable, reference_category
      );
    } else {
      update_indicator_interaction_pair_with_metropolis (
        pairwise_effects, main_effects, inclusion_indicator, observations,
        num_categories, proposal_sd_pairwise_effects, interaction_scale,
        index, num_pairwise, num_persons, residual_matrix,
        theta, is_ordinal_variable, reference_category
      );
    }
  }

  // Step 2: Update interaction weights for active edges
  if (use_mala) {
    update_interactions_with_mala (
        pairwise_effects, residual_matrix, main_effects, observations,
        num_categories, inclusion_indicator, is_ordinal_variable,
        reference_category, interaction_scale,
        step_size_mala_pairwise, initial_step_size_mala_pairwise,
        iteration, total_burnin, dual_averaging_pairwise
    );
  } else {
    update_interactions_with_adaptive_metropolis (
        pairwise_effects, main_effects, inclusion_indicator, observations,
        num_categories, proposal_sd_pairwise_effects, interaction_scale,
        num_persons, num_variables, residual_matrix,
        exp_neg_log_t_rm_adaptation_rate, is_ordinal_variable,
        reference_category
    );
  }

  // Step 3: Update main effect (threshold) parameters
  if (use_mala) {
    update_thresholds_with_fisher_mala (
      main_effects, step_size_mala_main, residual_matrix, num_categories,
      num_obs_categories, sufficient_blume_capel, reference_category,
      is_ordinal_variable, iteration, total_burnin, dual_averaging_main,
      sqrt_inv_fisher, threshold_alpha, threshold_beta,
      initial_step_size_mala_main
    );
  } else {
    for (int variable = 0; variable < num_variables; variable++) {
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
          threshold_beta, residual_matrix, proposal_sd_blumecapel,
          exp_neg_log_t_rm_adaptation_rate
        );
      }
    }
  }
}



/**
 * Function: run_gibbs_sampler_for_bgm
 *
 * Runs the full Gibbs sampler for a graphical model with ordinal and/or Blume-Capel variables.
 * Optionally performs edge selection using a Beta-Bernoulli or Stochastic Block prior.
 *
 * During each iteration, the algorithm:
 *  - Updates missing values (if imputation enabled)
 *  - Updates main effects (thresholds) via MALA or Metropolis-Hastings
 *  - Updates interaction parameters via Metropolis-Hastings
 *  - Updates inclusion indicators if edge selection is active
 *  - Adapts proposal variances using Robbins-Monro
 *  - Saves MCMC samples or updates running averages
 *
 * Inputs:
 *  - observations: Imputed categorical data.
 *  - num_categories: Number of categories per variable.
 *  - interaction_scale: Scale for Cauchy prior on pairwise interactions.
 *  - edge_prior: Type of edge prior ("Beta-Bernoulli" or "Stochastic-Block").
 *  - theta: Matrix of edge inclusion probabilities (updated if SBM/Beta-Bernoulli).
 *  - beta_bernoulli_alpha, beta_bernoulli_beta: Beta prior parameters for edge inclusion.
 *  - dirichlet_alpha, lambda: SBM prior parameters.
 *  - Index: Matrix of pairwise edge indices (E × 3).
 *  - iter: Number of post-burn-in iterations.
 *  - burnin: Number of burn-in iterations.
 *  - num_obs_categories: Category counts for each variable.
 *  - sufficient_blume_capel: Sufficient statistics for Blume-Capel variables.
 *  - threshold_alpha, threshold_beta: Prior parameters for logistic-Beta.
 *  - na_impute: Whether to impute missing values.
 *  - missing_index: List of missing value locations.
 *  - is_ordinal_variable: Logical vector (1 = ordinal, 0 = Blume-Capel).
 *  - reference_category: Reference category for centering (for Blume-Capel).
 *  - save_main, save_pairwise, save_indicator: Whether to save MCMC samples.
 *  - display_progress: Show progress bar during sampling.
 *  - edge_selection: Whether to enable edge selection during burn-in.
 *  - use_mala: Use MALA updates if true.
 *
 * Returns:
 *  - List containing:
 *    - "main": Posterior mean of main effects
 *    - "pairwise": Posterior mean of interaction effects
 *    - "inclusion_indicator": Posterior mean of inclusion matrix
 *    - (optional) "main_samples", "pairwise_samples", "inclusion_indicator_samples"
 *    - (optional) "allocations": Cluster allocations (if SBM used)
 */
// [[Rcpp::export]]
List run_gibbs_sampler_for_bgm (
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
    bool use_mala = false
) {
  // --- Setup: dimensions and storage structures
  const int num_variables = observations.n_cols;
  const int num_persons = observations.n_rows;
  const int max_num_categories = num_categories.max();
  const int num_pairwise = Index.n_rows;

  // Initialize model parameter matrices
  arma::mat main_effects(num_variables, max_num_categories, arma::fill::zeros);
  arma::mat pairwise_effects(num_variables, num_variables, arma::fill::zeros);
  arma::imat inclusion_indicator(num_variables, num_variables, arma::fill::ones);

  // Posterior mean accumulators
  arma::mat posterior_mean_main(num_variables, max_num_categories, arma::fill::zeros);
  arma::mat posterior_mean_pairwise(num_variables, num_variables, arma::fill::zeros);
  arma::mat posterior_mean_indicator(num_variables, num_variables, arma::fill::zeros);

  // Residuals used in pseudo-likelihood computation
  arma::mat residual_matrix(num_persons, num_variables, arma::fill::zeros);

  // Allocate optional storage for MCMC samples
  const int num_main = count_num_main_effects(num_categories, is_ordinal_variable);
  arma::mat* main_effect_samples = nullptr;
  arma::mat* pairwise_effect_samples = nullptr;
  arma::imat* indicator_samples = nullptr;

  if (save_main) main_effect_samples = new arma::mat(iter, num_main);
  if (save_pairwise) pairwise_effect_samples = new arma::mat(iter, num_pairwise);
  if (save_indicator) indicator_samples = new arma::imat(iter, num_pairwise);

  // Initialize proposal SDs and MALA tracking
  arma::mat proposal_sd_blumecapel(num_main, 2, arma::fill::ones);
  arma::mat proposal_sd_pairwise_effects(num_variables, num_variables, arma::fill::ones);

  double step_size_mala_main = 0.01;
  double step_size_mala_pairwise = 0.01;
  arma::vec dual_averaging_main(3, arma::fill::zeros);
  arma::vec dual_averaging_pairwise(3, arma::fill::zeros);
  double initial_step_size_mala_main = 0.01;
  double initial_step_size_mala_pairwise = 0.01;
  const double rm_decay_rate = 0.75;

  // Edge update shuffling setup
  arma::uvec v = arma::regspace<arma::uvec>(0, num_pairwise - 1);
  arma::uvec order(num_pairwise);
  arma::imat index(num_pairwise, 3);

  // SBM-specific structures
  arma::uvec K_values;
  arma::uvec cluster_allocations(num_variables);
  arma::mat cluster_prob(1, 1);
  arma::vec log_Vn(1);
  arma::imat out_allocations(iter, num_variables);

  // --- Initialize SBM prior if applicable
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

  // --- Optional MALA warmup stage (step size tuning only)
  arma::mat sqrt_inv_fisher(num_main, num_main, arma::fill::eye);  // Identity at start
  if (use_mala) {
    if (na_impute) {
      impute_missing_values_for_graphical_model(
        pairwise_effects, main_effects, observations, num_obs_categories,
        sufficient_blume_capel, num_categories, residual_matrix,
        missing_index, is_ordinal_variable, reference_category
      );
    }

    // Warmup phase for MALA: find initial step size (per Hoffman & Gelman, 2014)
    // Uses one MALA step to tune the log step size before dual averaging.
    initial_step_size_mala_main = find_reasonable_initial_step_size_mala_thresholds (
      main_effects, residual_matrix, num_categories, num_obs_categories,
      sufficient_blume_capel, reference_category, is_ordinal_variable,
      threshold_alpha, threshold_beta
    );
    initial_step_size_mala_pairwise = find_reasonable_initial_step_size_mala_interactions(
      pairwise_effects, main_effects, observations, num_categories,
      inclusion_indicator, is_ordinal_variable, reference_category, interaction_scale
    );

    step_size_mala_main = initial_step_size_mala_main;
    step_size_mala_pairwise = initial_step_size_mala_pairwise;

    dual_averaging_main[0] = std::log (step_size_mala_main);
    dual_averaging_pairwise[0] = std::log (step_size_mala_pairwise);
  }

  // --- Set up total number of iterations (burn-in + sampling)
  bool enable_edge_selection = edge_selection;
  int total_burnin = burnin * (enable_edge_selection ? 2 : 1);
  edge_selection = false;
  const int total_iter = total_burnin + iter;
  Progress p(total_iter, display_progress);


  // --- Main Gibbs sampling loop
  for (int iteration = 0; iteration < total_iter; iteration++) {
    if (Progress::check_abort()) {
      return List::create(
        Named("main") = posterior_mean_main,
        Named("pairwise") = posterior_mean_pairwise,
        Named("inclusion_indicator") = posterior_mean_indicator
      );
    }
    p.increment();

    // Re-enable edge selection halfway through burn-in
    if (enable_edge_selection && iteration == burnin) edge_selection = true;

    // Shuffle update order of edge indices
    order = arma::randperm(num_pairwise);
    for (int i = 0; i < num_pairwise; i++) {
      index.row(i) = Index.row(order(i));
    }

    // Optional imputation
    if (na_impute) {
      impute_missing_values_for_graphical_model(
        pairwise_effects, main_effects, observations, num_obs_categories,
        sufficient_blume_capel, num_categories, residual_matrix,
        missing_index, is_ordinal_variable, reference_category
      );
    }

    // Main Gibbs update step for parameters
    gibbs_update_step_for_graphical_model_parameters(
      observations, num_categories, interaction_scale,
      proposal_sd_pairwise_effects, proposal_sd_blumecapel, index,
      num_obs_categories, sufficient_blume_capel, threshold_alpha, threshold_beta,
      num_persons, num_variables, num_pairwise, num_main,
      inclusion_indicator, pairwise_effects, main_effects, residual_matrix, theta,
      rm_decay_rate, is_ordinal_variable, reference_category, edge_selection,
      step_size_mala_main, iteration, dual_averaging_main, total_burnin,
      use_mala, initial_step_size_mala_main, sqrt_inv_fisher,
      step_size_mala_pairwise, dual_averaging_pairwise,
      initial_step_size_mala_pairwise
    );

    // --- Update edge probabilities under the prior (if edge selection is active)
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

    // --- Save samples and update posterior means
    if (iteration >= total_burnin) {
      int iter_adj = iteration - total_burnin + 1;

      // Running posterior means
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
        for (int i = 0; i < num_pairwise; i++) {
          vectorized_pairwise(i) = pairwise_effects(Index(i, 1), Index(i, 2));
        }
        pairwise_effect_samples->row(sample_index) = vectorized_pairwise.t();
      }

      if (save_indicator) {
        arma::ivec vectorized_indicator(num_pairwise);
        for (int i = 0; i < num_pairwise; i++) {
          vectorized_indicator(i) = inclusion_indicator(Index(i, 1), Index(i, 2));
        }
        indicator_samples->row(sample_index) = vectorized_indicator.t();
      }

      if (edge_prior == "Stochastic-Block") {
        for (int j = 0; j < num_variables; j++)
          out_allocations(sample_index, j) = cluster_allocations[j] + 1;
      }
    }
  }

  // --- Final output
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
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppProgress)]]
#include <Rcpp.h>
#include "gibbs_functions_edge_prior.h"
#include <progress.hpp>
#include <progress_bar.hpp>
using namespace Rcpp;



/**
 * Function: update_step_size_with_dual_averaging
 * Purpose:
 *   Adapts the log step size in MCMC (e.g., MALA or HMC) using the dual averaging algorithm.
 *   This method adjusts the step size toward a target acceptance probability using an exponentially
 *   weighted average of the acceptance error. The adaptation is based on the algorithm described in:
 *
 *     Hoffman, M. D., & Gelman, A. (2014).
 *     The No-U-Turn Sampler: Adaptively Setting Path Lengths in Hamiltonian Monte Carlo.
 *     JMLR, 15(1), 1593–1623.
 *
 * Inputs:
 *   - acceptance_probability: The acceptance rate at the current iteration.
 *   - iteration: The current iteration number (starting from 1).
 *   - state: Vector of length 3, representing the dual averaging state. Updated in-place.
 *       - state[0] = log_step_size (current)
 *       - state[1] = log_step_size_avg (running average)
 *       - state[2] = acceptance_error_avg (running avg of acceptance error)
 *   - target_log_step_size: Desired log step size target (typically log(10 * initial_step_size)).
 *   - stabilization_offset: Constant to reduce adaptation sensitivity early on (e.g., 10).
 *
 * Outputs:
 *   - Updates the `state` vector in-place.
 *
 * Usage:
 *   - Used inside the adaptive Metropolis-adjusted Langevin algorithm (MALA) for updating threshold parameters.
 *   - Called during the burn-in phase inside `update_thresholds_with_adaptive_mala()` to adaptively tune the step size.
 */
inline void update_step_size_with_dual_averaging (
    const double acceptance_probability,
    const arma::uword iteration,
    arma::vec& state,
    const double target_log_step_size,
    const arma::uword stabilization_offset
) {
  double& log_step_size = state[0];
  double& log_step_size_avg = state[1];
  double& acceptance_error_avg = state[2];

  const double target_acceptance = 0.574;
  const double gamma = 0.05;
  const double kappa = 0.75;

  const double iter_double = static_cast<double>(iteration);
  const double adjusted_iter = iter_double + static_cast<double>(stabilization_offset);
  const double error = target_acceptance - acceptance_probability;

  // Update running average of acceptance error
  acceptance_error_avg = (1.0 - 1.0 / adjusted_iter) * acceptance_error_avg +
    1.0 / adjusted_iter * error;

  // Update log step size
  log_step_size = target_log_step_size - std::sqrt(iter_double) / gamma * acceptance_error_avg;

  // Update running average of log step size
  const double weight = std::pow(iter_double, -kappa);
  log_step_size_avg = weight * log_step_size + (1.0 - weight) * log_step_size_avg;
}



/**
 * Robbins-Monro update of the step size for MALA proposals using exponential decay.
 *
 * This function adapts the step size on the log scale to target a specific
 * acceptance probability. It uses a Robbins-Monro rule with polynomially
 * decaying learning rate:
 *
 *     log(step_size_mala) ← log(step_size_mala) + (accept_prob - target) * iteration^{-decay_rate_robins_monro}
 *
 * This ensures diminishing adaptation over time and avoids instability
 * in later iterations. The target acceptance rate (0.574) is optimal
 * for MALA under typical assumptions.
 *
 * Parameters:
 *  - acceptance_probability: Observed acceptance rate at the current iteration.
 *  - iteration: Current iteration count (starts at 1).
 *  - step_size_mala: The current step size (updated in-place).
 *
 * Notes:
 *  - decay_rate_robins_monro is fixed at 0.75 (typical value for stable decay).
 *  - This version assumes the step size is always positive.
 */
inline void update_step_size_with_robbins_monro(
    const double acceptance_probability,
    const arma::uword iteration,
    double& step_size_mala
) {
  const double target_acceptance = 0.574;
  const double decay_rate_robins_monro = 0.75;
  const double error = acceptance_probability - target_acceptance;
  const double decay = std::pow(static_cast<double>(iteration), -decay_rate_robins_monro);

  double log_step_size = std::log(step_size_mala);
  log_step_size += error * decay;
  step_size_mala = std::exp(log_step_size);
}



/**
 * Counts the total number of threshold parameters across all variables.
 *
 * Ordinal variables contribute one threshold per category.
 * Blume-Capel variables contribute two parameters (linear and quadratic).
 *
 * Parameters:
 *  - num_categories: Vector of category counts per variable.
 *  - is_ordinal_variable: Boolean vector indicating which variables are ordinal.
 *
 * Returns:
 *  - Total number of threshold parameters to be estimated.
 */
inline arma::uword count_num_main_effects(
    const arma::uvec& num_categories,
    const arma::uvec& is_ordinal_variable
) {
  arma::uword n_params = 0;
  for (arma::uword variable = 0; variable < num_categories.n_elem; variable++) {
    n_params += is_ordinal_variable[variable] ? num_categories[variable] : 2;
  }
  return n_params;
}



/**
 * Function: vectorize_thresholds
 * Purpose:
 *   Converts a matrix of threshold parameters into a flat vector for optimization,
 *   respecting the structure of ordinal vs. Blume-Capel variables.
 *
 * Inputs:
 *   - main_effects: Matrix of threshold parameters [variables × categories].
 *   - num_categories: Vector of the number of categories for each variable.
 *   - is_ordinal_variable: Logical vector indicating if a variable is ordinal (1) or Blume-Capel (0).
 *
 * Outputs:
 *   - Returns a flat vector of all threshold parameters, in order:
 *       - Ordinal variables: one value per category (length = num_categories[var])
 *       - Blume-Capel variables: exactly two values (linear and quadratic)
 *
 * Usage:
 *   - Used when preparing thresholds for MALA proposals or storing samples.
 */
arma::vec vectorize_thresholds(
    const arma::mat& main_effects,
    const arma::uvec& num_categories,
    const arma::uvec& is_ordinal_variable
) {
  const arma::uword num_parameters = count_num_main_effects(num_categories, is_ordinal_variable);
  arma::vec vector(num_parameters);
  arma::uword offset = 0;

  for (arma::uword variable = 0; variable < main_effects.n_rows; variable++) {
    const arma::uword num_pars = is_ordinal_variable(variable) ? num_categories(variable) : 2;
    vector.subvec(offset, offset + num_pars - 1) = main_effects.row(variable).cols(0, num_pars - 1).t();
    offset += num_pars;
  }
  return vector;
}



/**
 * Function: reconstruct_threshold_matrix
 * Purpose:
 *   Reconstructs a threshold matrix from a flat vector of parameters,
 *   reversing the operation of `vectorize_thresholds`.
 *
 * Inputs:
 *   - flat_vector: Vector containing all threshold parameters in a flattened form.
 *   - num_categories: Vector of the number of categories per variable.
 *   - is_ordinal_variable: Logical vector indicating if a variable is ordinal (1) or Blume-Capel (0).
 *
 * Outputs:
 *   - Returns a matrix of thresholds [variables × max(num_categories)], with unused entries (if any) uninitialized.
 *     Each row contains:
 *       - Ordinal: num_categories[var] entries
 *       - Blume-Capel: 2 entries (linear and quadratic thresholds)
 *
 * Usage:
 *   - Used to convert parameter vectors (e.g., from MALA or sampling storage) back to structured threshold matrices.
 */
arma::mat reconstruct_threshold_matrix(
    const arma::vec& vector,
    const arma::uvec& num_categories,
    const arma::uvec& is_ordinal_variable
) {
  const arma::uword num_variables = num_categories.n_elem;
  const arma::uword max_categories = num_categories.max();
  arma::mat matrix(num_variables, max_categories, arma::fill::zeros); // May contain unused entries

  arma::uword offset = 0;
  for (arma::uword variable = 0; variable < num_variables; variable++) {
    arma::uword num_pars = is_ordinal_variable[variable] ? num_categories[variable] : 2;
    matrix.row(variable).cols(0, num_pars - 1) = vector.subvec(offset, offset + num_pars - 1).t();
    offset += num_pars;
  }

  return matrix;
}



/**
 * Function: update_proposal_sd_with_robbins_monro
 * Purpose: Performs Robbins-Monro updates for proposal standard deviations.
 *
 * Inputs:
 *  - current_sd: Current standard deviation of the proposal.
 *  - observed_log_acceptance_probability: Log acceptance probability from the Metropolis-Hastings step.
 *  - target_acceptance_rate: Target acceptance rate (e.g. 0.234 or 0.574).
 *  - rm_lower_bound: Minimum allowable standard deviation.
 *  - rm_upper_bound: Maximum allowable standard deviation.
 *  - rm_weight: Robbins-Monro adaptation weight (e.g. iteration^{-0.75}).
 *
 * Returns:
 *  - Updated proposal standard deviation, clamped within bounds.
 */
inline double update_proposal_sd_with_robbins_monro (
    const double current_sd,
    const double observed_log_acceptance_probability,
    const double target_acceptance_rate,
    const double rm_lower_bound,
    const double rm_upper_bound,
    const double rm_weight
) {
  // Normalize the acceptance probability
  double observed_acceptance_probability = 1.0;
  if (observed_log_acceptance_probability < 0.0) {
    observed_acceptance_probability = std::exp(observed_log_acceptance_probability);
  }

  // Robbins-Monro update step
  double updated_sd = current_sd +
    (observed_acceptance_probability - target_acceptance_rate) * rm_weight;

  // Handle NaNs robustly
  if (std::isnan(updated_sd)) {
    updated_sd = 1.0;  // Safe default
  }

  // Clamp the value within specified bounds
  return std::clamp(updated_sd, rm_lower_bound, rm_upper_bound);
}



/**
 * Function: impute_missing_values_for_graphical_model
 *
 * Imputes missing values in the data matrix using current model parameters and updates
 * sufficient statistics and residual matrix accordingly.
 *
 * Inputs:
 *  - pairwise_effects: current matrix of interaction weights.
 *  - main_effects: current matrix of threshold parameters.
 *  - missing_index: matrix of indices (i, j) for missing values.
 *  - num_categories: vector of category counts per variable.
 *  - is_ordinal_variable: indicator for ordinal vs. Blume-Capel.
 *  - reference_category: vector of reference categories.
 *
 * Updates (in-place):
 *  - observations
 *  - num_obs_categories
 *  - sufficient_blume_capel
 *  - residual_matrix
 */
void impute_missing_values_for_graphical_model (
    const arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    arma::umat& observations,
    arma::umat& num_obs_categories,
    arma::umat& sufficient_blume_capel,
    const arma::uvec& num_categories,
    arma::mat& residual_matrix,
    const arma::umat& missing_index,
    const arma::uvec& is_ordinal_variable,
    const arma::uvec& reference_category
) {
  const arma::uword num_variables = observations.n_cols;
  const arma::uword num_missings = missing_index.n_rows;
  const arma::uword max_num_categories = num_categories.max();

  arma::vec category_probabilities(max_num_categories + 1);

  for (arma::uword miss = 0; miss < num_missings; miss++) {
    const arma::uword person = missing_index(miss, 0);
    const arma::uword variable = missing_index(miss, 1);

    const double rest_score = residual_matrix(person, variable);
    const arma::uword num_cats = num_categories(variable);
    const bool is_ordinal = is_ordinal_variable(variable);

    double cumsum = 0.0;

    if (is_ordinal) {
      // Compute cumulative unnormalized probabilities for ordinal variable
      cumsum = 1.0;
      category_probabilities[0] = cumsum;
      for (arma::uword cat = 0; cat < num_cats; cat++) {
        const arma::uword score = cat + 1;
        const double exponent = main_effects(variable, cat) + score * rest_score;
        cumsum += std::exp(exponent);
        category_probabilities[score] = cumsum;
      }
    } else {
      // Compute probabilities for Blume-Capel variable
      const arma::uword ref = reference_category(variable);

      cumsum = std::exp(main_effects(variable, 1) * ref * ref);
      category_probabilities[0] = cumsum;

      for (arma::uword cat = 0; cat < num_cats; cat++) {
        const arma::uword score = cat + 1;
        const arma::sword centered = static_cast<arma::sword>(score) - static_cast<arma::sword>(ref);
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
    arma::uword sampled_score = 0;
    while (u > category_probabilities[sampled_score]) {
      sampled_score++;
    }

    const arma::uword new_value = sampled_score;
    const arma::uword old_value = observations(person, variable);

    if (new_value != old_value) {
      // Update observation matrix
      observations(person, variable) = new_value;

      if (is_ordinal) {
        num_obs_categories(old_value, variable)--;
        num_obs_categories(new_value, variable)++;
      } else {
        const arma::sword ref = reference_category(variable);
        const arma::sword delta = static_cast<arma::sword>(new_value) - static_cast<arma::sword>(old_value);
        const arma::sword delta_sq =
          (new_value - ref) * (new_value - ref) -
          (old_value - ref) * (old_value - ref);

        sufficient_blume_capel(0, variable) += delta;
        sufficient_blume_capel(1, variable) += delta_sq;
      }

      // Update residuals across all variables
      for (arma::uword v = 0; v < num_variables; v++) {
        const double delta_score = (static_cast<double>(new_value) - old_value) * pairwise_effects(v, variable);
        residual_matrix(person, v) += delta_score;
      }
    }
  }
}


/**
 * Function: log_pseudoposterior_thresholds
 * Purpose:
 *   Computes the log pseudo-posterior for threshold parameters in a graphical model
 *   combining contributions from the data likelihood and logistic-Beta priors.
 *   Supports both regular ordinal and Blume-Capel variables.
 *
 * Inputs:
 *   - main_effects: Matrix of threshold parameters [variables × categories].
 *   - residual_matrix: Matrix of residual predictors (individuals × variables).
 *   - num_categories: Vector of the number of score levels for each variable.
 *   - num_obs_categories: Matrix of observed category counts [category × variable].
 *   - sufficient_blume_capel: Matrix of sufficient statistics for Blume-Capel variables [2 × variable].
 *   - reference_category: Vector of reference (centering) categories for each variable.
 *   - is_ordinal_variable: Logical vector (1 if ordinal, 0 if Blume-Capel).
 *   - threshold_alpha: Hyperparameter alpha for logistic-Beta prior (default = 1.0).
 *   - threshold_beta: Hyperparameter beta for logistic-Beta prior (default = 1.0).
 *
 * Outputs:
 *   - Returns the total log pseudo-posterior (data likelihood + priors) as a scalar `double`.
 *
 * Usage:
 *   - Called inside `update_thresholds_with_adaptive_mala()` and Metropolis updates for thresholds.
 */
double log_pseudoposterior_thresholds (
    const arma::mat& main_effects,
    const arma::mat& residual_matrix,
    const arma::uvec& num_categories,
    const arma::umat& num_obs_categories,
    const arma::umat& sufficient_blume_capel,
    const arma::uvec& reference_category,
    const arma::uvec& is_ordinal_variable,
    const double threshold_alpha = 1.0,
    const double threshold_beta = 1.0
) {
  const arma::uword num_variables = residual_matrix.n_cols;
  const arma::uword num_persons = residual_matrix.n_rows;

  double log_posterior = 0.0;

  for (arma::uword variable = 0; variable < num_variables; variable++) {
    const arma::uword num_cats = num_categories(variable);

    if (is_ordinal_variable(variable)) {
      // Regular ordinal variable
      for (arma::uword cat = 0; cat < num_cats; cat++) {
        const double theta = main_effects(variable, cat);
        log_posterior += theta * (num_obs_categories(cat + 1, variable) + threshold_alpha);
        log_posterior -= std::log1p(std::exp(theta)) * (threshold_alpha + threshold_beta);
      }

      for (arma::uword person = 0; person < num_persons; person++) {
        const double rest_score = residual_matrix(person, variable);
        const double bound = num_cats * rest_score;

        double denominator = std::exp(-bound);
        for (arma::uword cat = 0; cat < num_cats; cat++) {
          const double exponent = main_effects(variable, cat) + (cat + 1) * rest_score - bound;
          denominator += std::exp(exponent);
        }

        log_posterior -= bound + std::log(denominator);
      }

    } else {
      // Blume-Capel variable
      const double theta_lin = main_effects(variable, 0);
      const double theta_quad = main_effects(variable, 1);

      log_posterior += theta_lin * (sufficient_blume_capel(0, variable) + threshold_alpha);
      log_posterior -= std::log1p(std::exp(theta_lin)) * (threshold_alpha + threshold_beta);

      log_posterior += theta_quad * (sufficient_blume_capel(1, variable) + threshold_alpha);
      log_posterior -= std::log1p(std::exp(theta_quad)) * (threshold_alpha + threshold_beta);

      const int ref = static_cast<int>(reference_category(variable));

      for (arma::uword person = 0; person < num_persons; person++) {
        const double rest_score = residual_matrix(person, variable);
        const double bound = num_cats * rest_score;

        double denominator = 0.0;
        for (int cat = 0; cat <= num_cats; cat++) {
          const int centered = cat - ref;
          const double exponent =
            theta_lin * cat +
            theta_quad * centered * centered +
            cat * rest_score - bound;
          denominator += std::exp(exponent);
        }

        log_posterior -= bound + std::log(denominator);
      }
    }
  }

  return log_posterior;
}


/**
 * Computes the gradient of the log pseudo-posterior with respect to threshold parameters.
 *
 * Supports both ordinal and Blume-Capel variables. The gradient includes:
 *   - Data pseudo-likelihood contribution (via expected counts)
 *   - Logistic-Beta prior contributions on the transformed thresholds
 *
 * Parameters:
 *  - main_effects: [variables × categories] matrix of current threshold parameters.
 *  - residual_matrix: [persons × variables] matrix of residual scores.
 *  - num_categories: Vector with number of score levels per variable.
 *  - num_obs_categories: [category × variable] matrix of observed category counts.
 *  - sufficient_blume_capel: [2 × variable] matrix of sufficient stats.
 *  - reference_category: Reference (centering) category per variable.
 *  - is_ordinal_variable: 1 = ordinal, 0 = Blume-Capel.
 *  - threshold_alpha, threshold_beta: Hyperparameters for the logistic-Beta prior.
 *
 * Returns:
 *  - A vector containing gradients, flattened to match threshold vector layout.
 */
arma::vec gradient_log_pseudoposterior_thresholds(
    const arma::mat& main_effects,
    const arma::mat& residual_matrix,
    const arma::uvec& num_categories,
    const arma::umat& num_obs_categories,
    const arma::umat& sufficient_blume_capel,
    const arma::uvec& reference_category,
    const arma::uvec& is_ordinal_variable,
    const double threshold_alpha = 1.0,
    const double threshold_beta = 1.0
) {
  const arma::uword num_variables = residual_matrix.n_cols;
  const arma::uword num_persons = residual_matrix.n_rows;
  const arma::uword num_parameters = count_num_main_effects(num_categories, is_ordinal_variable);

  arma::vec gradient(num_parameters, arma::fill::zeros);
  arma::uword offset = 0;

  for (arma::uword variable = 0; variable < num_variables; variable++) {
    const arma::uword num_cats = num_categories(variable);

    if (is_ordinal_variable(variable)) {
      // Ordinal variable
      const double max_threshold = main_effects.row(variable).max();

      for (arma::uword cat = 0; cat < num_cats; cat++) {
        gradient(offset + cat) = num_obs_categories(cat + 1, variable);
      }

      for (arma::uword person = 0; person < num_persons; person++) {
        const double rest_score = residual_matrix(person, variable);
        const double bound = max_threshold + num_cats * rest_score;

        double denom = std::exp(-bound);
        arma::vec numerators(num_cats, arma::fill::zeros);

        for (arma::uword cat = 0; cat < num_cats; cat++) {
          const double exponent = main_effects(variable, cat) + (cat + 1) * rest_score - bound;
          numerators(cat) = std::exp(exponent);
          denom += numerators(cat);
        }

        for (arma::uword cat = 0; cat < num_cats; cat++) {
          gradient(offset + cat) -= numerators(cat) / denom;
        }
      }

      for (arma::uword cat = 0; cat < num_cats; cat++) {
        const double theta = main_effects(variable, cat);
        const double p = 1.0 / (1.0 + std::exp(-theta));
        gradient(offset + cat) += threshold_alpha - (threshold_alpha + threshold_beta) * p;
      }

      offset += num_cats;

    } else {
      // Blume-Capel variable
      const int ref = static_cast<int>(reference_category(variable));
      const double theta_lin = main_effects(variable, 0);
      const double theta_quad = main_effects(variable, 1);

      gradient(offset)     = sufficient_blume_capel(0, variable);
      gradient(offset + 1) = sufficient_blume_capel(1, variable);

      for (arma::uword person = 0; person < num_persons; person++) {
        const double rest_score = residual_matrix(person, variable);
        const double bound = num_cats * rest_score;

        double denom = std::exp(theta_quad * ref * ref - bound);
        double sum_lin = 0.0;
        double sum_quad = ref * ref * denom;

        for (arma::uword cat = 0; cat < num_cats; cat++) {
          const arma::uword score = cat + 1;
          const int centered = static_cast<int>(score) - ref;

          const double exponent =
            theta_lin * score +
            theta_quad * centered * centered +
            score * rest_score - bound;

          const double weight = std::exp(exponent);
          sum_lin += weight * score;
          sum_quad += weight * centered * centered;

          denom += weight;
        }

        gradient(offset)     -= sum_lin / denom;
        gradient(offset + 1) -= sum_quad / denom;
      }

      for (arma::uword i = 0; i < 2; i++) {
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
 * Performs a MALA (Metropolis-adjusted Langevin Algorithm) update of threshold parameters
 * with adaptive step size tuning. Supports both ordinal and Blume-Capel variables.
 *
 * During burn-in, step size is adapted using dual averaging. After burn-in,
 * Robbins-Monro is used.
 *
 * Parameters:
 *  - main_effects: [variables × max_categories] matrix of current thresholds (updated in-place).
 *  - step_size_mala: MALA step size (updated in-place).
 *  - residual_matrix: [persons × variables] matrix of residual scores.
 *  - num_categories: Number of categories per variable.
 *  - num_obs_categories: Count of observed categories [category × variable].
 *  - sufficient_blume_capel: Sufficient stats for BC variables [2 × variable].
 *  - reference_category: Reference category for centering.
 *  - is_ordinal_variable: Indicator for ordinal (1) or Blume-Capel (0) per variable.
 *  - iteration: Current iteration.
 *  - burnin: Number of burn-in iterations.
 *  - dual_averaging_state: [log_eps, log_eps_bar, H_bar] vector (updated in-place).
 *  - target_log_step_size: Target log step size (typically log(10 * ε₀)).
 *  - threshold_alpha, threshold_beta: Hyperparameters for logistic-Beta prior.
 */
void update_thresholds_with_adaptive_mala (
    arma::mat& main_effects,
    double& step_size_mala,
    const arma::mat& residual_matrix,
    const arma::uvec& num_categories,
    const arma::umat& num_obs_categories,
    const arma::umat& sufficient_blume_capel,
    const arma::uvec& reference_category,
    const arma::uvec& is_ordinal_variable,
    const arma::uword iteration,
    const arma::uword burnin,
    arma::vec& dual_averaging_state,
    const double target_log_step_size,
    const double threshold_alpha,
    const double threshold_beta
) {
  // Flatten parameters
  arma::vec flat_theta = vectorize_thresholds(main_effects, num_categories, is_ordinal_variable);

  // Evaluate current gradient and log posterior
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

  // Propose new theta
  const double sqrt_step = std::sqrt(step_size_mala);
  arma::vec proposal = flat_theta + 0.5 * step_size_mala * grad + sqrt_step * arma::randn(flat_theta.n_elem);
  arma::mat proposed_thresholds = reconstruct_threshold_matrix(proposal, num_categories, is_ordinal_variable);

  // Evaluate proposal log posterior and gradient
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

  // Compute forward and backward proposal densities
  const arma::vec forward_mean = flat_theta + 0.5 * step_size_mala * grad;
  const arma::vec backward_mean = proposal + 0.5 * step_size_mala * grad_proposal;

  const double log_forward = -0.5 / step_size_mala * arma::accu(arma::square(proposal - forward_mean));
  const double log_backward = -0.5 / step_size_mala * arma::accu(arma::square(flat_theta - backward_mean));

  // Metropolis acceptance probability
  const double log_acceptance = log_post_proposal + log_backward - log_post_current - log_forward;

  // Accept/reject
  if (std::log(R::unif_rand()) < log_acceptance) {
    main_effects = proposed_thresholds;
    flat_theta = proposal;
  }

  // Compute acceptance probability (bounded at 1)
  const double accept_prob = std::min(1.0, std::exp(log_acceptance));

  // Adapt step size
  if (iteration <= burnin) {
    update_step_size_with_dual_averaging(
      accept_prob, iteration, dual_averaging_state,
      target_log_step_size, 10
    );
    step_size_mala = std::exp(dual_averaging_state[1]);
  } else {
    update_step_size_with_robbins_monro(accept_prob, iteration - burnin, step_size_mala);
  }
}



/**
 * Function: update_regular_thresholds_with_metropolis
 * Purpose:
 *   Performs a Metropolis-Hastings update for threshold parameters of ordinal variables.
 *   Each threshold is updated one at a time using a generalized beta-prime proposal,
 *   with acceptance determined by the pseudo-likelihood and logistic-Beta prior.
 *
 * Inputs:
 *   - main_effects: Matrix of threshold parameters [variables × categories]; updated in-place.
 *   - observations: Matrix of observed scores [persons × variables].
 *   - num_categories: Vector of number of categories per variable.
 *   - num_obs_categories: Matrix of counts per (category, variable).
 *   - no_persons: Number of individuals in the data.
 *   - variable: Index of the variable being updated.
 *   - threshold_alpha: Alpha parameter of the logistic-Beta prior.
 *   - threshold_beta: Beta parameter of the logistic-Beta prior.
 *   - residual_matrix: Matrix of linear predictors excluding current variable.
 *
 * Outputs:
 *   - Updates `main_effects(variable, category)` in-place for each category of the specified variable.
 *
 * Usage:
 *   - Called during Gibbs updates when using non-gradient MH threshold proposals.
 */
void update_regular_thresholds_with_metropolis (
    arma::mat& main_effects,
    const arma::umat& observations,
    const arma::uvec& num_categories,
    const arma::umat& num_obs_categories,
    const arma::uword no_persons,
    const arma::uword variable,
    const double threshold_alpha,
    const double threshold_beta,
    const arma::mat& residual_matrix
) {
  arma::vec g(no_persons);
  arma::vec q(no_persons);

  for (arma::uword category = 0; category < num_categories[variable]; category++) {
    double current = main_effects(variable, category);
    double exp_current = std::exp(current);
    double c = (threshold_alpha + threshold_beta) / (1 + exp_current);

    // Compute proposal scaling constant `c`
    for (arma::uword person = 0; person < no_persons; person++) {
      double rest_score = residual_matrix(person, variable);
      double denom = 1.0;  // base of sum (g)
      double numer = std::exp((category + 1) * rest_score);  // q

      for (arma::uword cat = 0; cat < num_categories[variable]; cat++) {
        if (cat != category) {
          denom += std::exp(main_effects(variable, cat) + (cat + 1) * rest_score);
        }
      }

      g[person] = denom;
      q[person] = numer;
      c += q[person] / (g[person] + q[person] * exp_current);
    }

    c /= ((no_persons + threshold_alpha + threshold_beta) -
      exp_current * c);

    // Sample from generalized beta-prime proposal
    double a = num_obs_categories(category + 1, variable) + threshold_alpha;
    double b = no_persons + threshold_beta - num_obs_categories(category + 1, variable);
    double tmp = R::rbeta(a, b);
    double proposed = std::log(tmp / (1.0 - tmp) / c);
    double exp_proposed = std::exp(proposed);

    // Compute MH acceptance probability
    double log_acceptance_probability = 0.0;
    for (arma::uword person = 0; person < no_persons; person++) {
      log_acceptance_probability += std::log(g[person] + q[person] * exp_current);
      log_acceptance_probability -= std::log(g[person] + q[person] * exp_proposed);
    }

    // Add prior ratio (logistic-Beta)
    log_acceptance_probability -= (threshold_alpha + threshold_beta) * std::log(1 + exp_proposed);
    log_acceptance_probability += (threshold_alpha + threshold_beta) * std::log(1 + exp_current);

    // Add proposal ratio (generalized beta-prime)
    log_acceptance_probability -= (a + b) * std::log(1 + c * exp_current);
    log_acceptance_probability += (a + b) * std::log(1 + c * exp_proposed);

    // Metropolis step
    if (std::log(R::unif_rand()) < log_acceptance_probability) {
      main_effects(variable, category) = proposed;
    }
  }
}



/**
 * Function: update_blumecapel_thresholds_with_adaptive_metropolis
 *
 * Purpose:
 * Performs an adaptive Metropolis update of the Blume-Capel threshold parameters
 * (linear and quadratic components) for a single variable, using a log-pseudolikelihood
 * and a logistic-Beta prior. Robbins-Monro adaptation is used to tune proposal variances.
 *
 * Inputs:
 *  - main_effects: [variables × 2] matrix of threshold parameters (linear and quadratic).
 *  - observations: [persons × variables] matrix of ordinal responses.
 *  - num_categories: Vector of number of categories per variable.
 *  - sufficient_blume_capel: [2 × variables] matrix of sufficient statistics.
 *  - no_persons: Number of observations (rows in observations).
 *  - variable: Index of the variable being updated (0-based).
 *  - reference_category: Vector of reference categories per variable.
 *  - threshold_alpha, threshold_beta: Shape parameters of the Beta-prime prior.
 *  - residual_matrix: [persons × variables] matrix of linear predictors.
 *  - proposal_sd_blumecapel: [variables × 2] matrix of proposal standard deviations.
 *  - exp_neg_log_t_rm_adaptation_rate: Robbins-Monro weight: exp(-log(iter) * rate).
 *  - target_acceptance_rate_mh: Target acceptance rate for Metropolis (e.g., 0.234).
 *  - rm_lower_bound, rm_upper_bound: Bounds for clamping proposal SD.
 *
 * Updates:
 *  - main_effects: Updated in-place for the specified variable.
 *  - proposal_sd_blumecapel: Updated in-place using Robbins-Monro rule.
 */
void update_blumecapel_thresholds_with_adaptive_metropolis (
    arma::mat& main_effects,
    const arma::umat& observations,
    const arma::uvec& num_categories,
    const arma::umat& sufficient_blume_capel,
    const arma::uword no_persons,
    const arma::uword variable,
    const arma::uvec& reference_category,
    const double threshold_alpha,
    const double threshold_beta,
    const arma::mat& residual_matrix,
    arma::mat& proposal_sd_blumecapel,
    const double exp_neg_log_t_rm_adaptation_rate,
    const double target_acceptance_rate_mh,
    const double rm_lower_bound,
    const double rm_upper_bound
) {
  const arma::uword num_cats = num_categories(variable);
  const int ref = static_cast<int>(reference_category(variable));

  // Helper lambda to update either the linear (0) or quadratic (1) threshold
  auto update_parameter = [&](arma::uword param_index) {
    const double current = main_effects(variable, param_index);
    const double proposed = R::rnorm(current, proposal_sd_blumecapel(variable, param_index));
    const double diff = proposed - current;

    arma::vec numer_current(num_cats + 1), numer_proposed(num_cats + 1);

    for (arma::uword cat = 0; cat <= num_cats; cat++) {
      const int centered = static_cast<int>(cat) - ref;
      if (param_index == 0) { // linear update
        double quad_term = main_effects(variable, 1) * (centered * centered);
        numer_current(cat) = current * cat + quad_term;
        numer_proposed(cat) = proposed * cat + quad_term;
      } else { // quadratic update
        double lin_term = main_effects(variable, 0) * cat;
        numer_current(cat) = current * (centered * centered) + lin_term;
        numer_proposed(cat) = proposed * (centered * centered) + lin_term;
      }
    }

    const double max_curr = numer_current.max();
    const double max_prop = numer_proposed.max();
    const double lbound = (max_curr > 0 || max_prop > 0) ? std::max(max_curr, max_prop) : 0.0;

    double log_acceptance_probability = threshold_alpha * diff +
      static_cast<double>(sufficient_blume_capel(param_index, variable)) * diff;

    for (arma::uword person = 0; person < no_persons; person++) {
      const double rest_score = residual_matrix(person, variable);
      const double bound = (rest_score > 0) ? num_cats * rest_score + lbound : lbound;

      double denom_curr = std::exp(numer_current(0) - bound);
      double denom_prop = std::exp(numer_proposed(0) - bound);

      for (arma::uword cat = 0; cat < num_cats; cat++) {
        const double score_term = (cat + 1) * rest_score - bound;
        denom_curr += std::exp(numer_current(cat + 1) + score_term);
        denom_prop += std::exp(numer_proposed(cat + 1) + score_term);
      }

      log_acceptance_probability += std::log(denom_curr) - std::log(denom_prop);
    }

    // Prior ratio
    log_acceptance_probability += (threshold_alpha + threshold_beta) * (
      std::log1p(std::exp(current)) - std::log1p(std::exp(proposed))
    );

    // Metropolis accept/reject
    const double logu = std::log(R::unif_rand());
    if (logu < log_acceptance_probability) {
      main_effects(variable, param_index) = proposed;
    }

    // Robbins-Monro adaptation
    proposal_sd_blumecapel(variable, param_index) =
      update_proposal_sd_with_robbins_monro(
        proposal_sd_blumecapel(variable, param_index),
        log_acceptance_probability,
        target_acceptance_rate_mh,
        rm_lower_bound,
        rm_upper_bound,
        exp_neg_log_t_rm_adaptation_rate
      );
  };

  // Update both parameters
  update_parameter(0); // linear
  update_parameter(1); // quadratic
}



/**
 * Computes the log pseudo-likelihood ratio between a proposed and current value
 * of an interaction parameter between two variables.
 *
 * This is used in Metropolis-Hastings updates for interaction parameters.
 *
 * Parameters:
 *  - pairwise_effects: [V × V] matrix of interaction parameters.
 *  - main_effects: [V × max_categories] matrix of threshold parameters.
 *  - observations: [N × V] matrix of observed scores (integer encoded).
 *  - num_categories: Number of categories per variable.
 *  - no_persons: Number of observations (rows in observations matrix).
 *  - variable1, variable2: Indices of the two interacting variables.
 *  - proposed_state: Proposed new interaction value.
 *  - current_state: Current interaction value.
 *  - residual_matrix: [N × V] matrix of residual scores (linear predictors).
 *  - is_ordinal_variable: Indicator vector (1 = ordinal, 0 = Blume-Capel).
 *  - reference_category: Centering category per variable (for Blume-Capel).
 *
 * Returns:
 *  - The log pseudo-likelihood ratio: log p(y | β_proposed) - log p(y | β_current)
 */
double log_pseudolikelihood_ratio_interaction (
    const arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    const arma::umat& observations,
    const arma::uvec& num_categories,
    const arma::uword no_persons,
    const arma::uword variable1,
    const arma::uword variable2,
    const double proposed_state,
    const double current_state,
    const arma::mat& residual_matrix,
    const arma::uvec& is_ordinal_variable,
    const arma::uvec& reference_category
) {
  double log_ratio = 0.0;
  const double delta = proposed_state - current_state;

  for (arma::uword person = 0; person < no_persons; person++) {
    const arma::uword score1 = observations(person, variable1);
    const arma::uword score2 = observations(person, variable2);

    // Linear interaction contribution to likelihood
    log_ratio += 2.0 * score1 * score2 * delta;

    // Variable 1: log-likelihood change
    {
      double rest_score = residual_matrix(person, variable1) - score2 * current_state;
      double bound = rest_score > 0 ? num_categories(variable1) * rest_score : 0.0;


      double denom_current = 0.0;
      double denom_proposed = 0.0;

      if (is_ordinal_variable(variable1)) {
        denom_current += std::exp(-bound);
        denom_proposed += std::exp(-bound);
        for (arma::uword cat = 0; cat < num_categories(variable1); cat++) {
          const double exponent = main_effects(variable1, cat) + (cat + 1) * rest_score;
          denom_current += std::exp(exponent + (cat + 1) * score2 * current_state - bound);
          denom_proposed += std::exp(exponent + (cat + 1) * score2 * proposed_state - bound);
        }
      } else {
        const arma::uword ref = reference_category(variable1);
        for (arma::uword cat = 0; cat <= num_categories(variable1); cat++) {
          const int centered = static_cast<int>(cat) - static_cast<int>(ref);
          const double exponent =
            main_effects(variable1, 0) * cat +
            main_effects(variable1, 1) * centered * centered +
            cat * rest_score - bound;
          denom_current += std::exp(exponent + cat * score2 * current_state);
          denom_proposed += std::exp(exponent + cat * score2 * proposed_state);
        }
      }

      log_ratio += std::log(denom_current) - std::log(denom_proposed);
    }

    // Variable 2: log-likelihood change (symmetric)
    {
      double rest_score = residual_matrix(person, variable2) - score1 * current_state;
      double bound = rest_score > 0 ? num_categories(variable2) * rest_score : 0.0;

      double denom_current = 0.0;
      double denom_proposed = 0.0;

      if (is_ordinal_variable(variable2)) {
        denom_current += std::exp(-bound);
        denom_proposed += std::exp(-bound);
        for (arma::uword cat = 0; cat < num_categories(variable2); cat++) {
          const double exponent = main_effects(variable2, cat) + (cat + 1) * rest_score;
          denom_current += std::exp(exponent + (cat + 1) * score1 * current_state - bound);
          denom_proposed += std::exp(exponent + (cat + 1) * score1 * proposed_state - bound);
        }
      } else {
        const arma::uword ref = reference_category(variable2);
        for (arma::uword cat = 0; cat <= num_categories(variable2); cat++) {
          const int centered = static_cast<int>(cat) - static_cast<int>(ref);
          const double exponent =
            main_effects(variable2, 0) * cat +
            main_effects(variable2, 1) * centered * centered +
            cat * rest_score - bound;
          denom_current += std::exp(exponent + cat * score1 * current_state);
          denom_proposed += std::exp(exponent + cat * score1 * proposed_state);
        }
      }

      log_ratio += std::log(denom_current) - std::log(denom_proposed);
    }
  }

  return log_ratio;
}



/**
 * Function: update_interactions_with_adaptive_metropolis
 *
 * Performs adaptive Metropolis-Hastings updates for all active pairwise interactions.
 *
 * For each pair (i, j) where `inclusion_indicator(i, j) == 1`, proposes a new value for the
 * interaction strength from a Gaussian proposal distribution. The proposal is accepted or
 * rejected based on the log-pseudo-likelihood ratio and a symmetric Cauchy prior.
 *
 * Proposal standard deviations are updated using Robbins-Monro adaptation.
 *
 * Inputs:
 *  - pairwise_effects: current matrix of pairwise interaction parameters.
 *  - main_effects: matrix of main effect parameters (used in the likelihood).
 *  - inclusion_indicator: binary matrix indicating active pairwise effects.
 *  - observations: matrix of category scores.
 *  - num_categories: number of categories for each variable.
 *  - proposal_sd_pairwise_effects: matrix of proposal SDs (adapted in-place).
 *  - interaction_scale: scale parameter for the Cauchy prior.
 *  - num_persons: number of observations.
 *  - num_variables: number of variables.
 *  - residual_matrix: matrix of linear predictors (updated in-place if proposal accepted).
 *  - exp_neg_log_t_rm_adaptation_rate: Robbins-Monro weight: exp(-log(t) * rate).
 *  - target_acceptance_rate_mh: target MH acceptance rate.
 *  - rm_lower_bound, rm_upper_bound: bounds for Robbins-Monro proposal SDs.
 *  - is_ordinal_variable: indicator for ordinal vs. Blume-Capel.
 *  - reference_category: reference category per variable.
 *
 * Updates (in-place):
 *  - pairwise_effects
 *  - residual_matrix
 *  - proposal_sd_pairwise_effects
 */
void update_interactions_with_adaptive_metropolis (
    arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    const arma::umat& inclusion_indicator,
    const arma::umat& observations,
    const arma::uvec& num_categories,
    arma::mat& proposal_sd_pairwise_effects,
    const double interaction_scale,
    const arma::uword num_persons,
    const arma::uword num_variables,
    arma::mat& residual_matrix,
    const double exp_neg_log_t_rm_adaptation_rate,
    const double target_acceptance_rate_mh,
    const double rm_lower_bound,
    const double rm_upper_bound,
    const arma::uvec& is_ordinal_variable,
    const arma::uvec& reference_category
) {
  for (arma::uword variable1 = 0; variable1 < num_variables - 1; variable1++) {
    for (arma::uword variable2 = variable1 + 1; variable2 < num_variables; variable2++) {
      if (inclusion_indicator(variable1, variable2) == 1) {
        const double current_state = pairwise_effects(variable1, variable2);
        const double proposed_state = R::rnorm(current_state, proposal_sd_pairwise_effects(variable1, variable2));

        double log_acceptance_probability = log_pseudolikelihood_ratio_interaction(
          pairwise_effects, main_effects, observations, num_categories, num_persons,
          variable1, variable2, proposed_state, current_state,
          residual_matrix, is_ordinal_variable, reference_category
        );

        // Add log-prior difference (Cauchy prior)
        log_acceptance_probability += R::dcauchy(proposed_state, 0.0, interaction_scale, true);
        log_acceptance_probability -= R::dcauchy(current_state, 0.0, interaction_scale, true);

        // Accept/reject
        const double U = R::unif_rand();
        if (std::log(U) < log_acceptance_probability) {
          double state_diff = proposed_state - current_state;
          pairwise_effects(variable1, variable2) = proposed_state;
          pairwise_effects(variable2, variable1) = proposed_state;

          // Update residual matrix
          for (arma::uword person = 0; person < num_persons; person++) {
            residual_matrix(person, variable1) += observations(person, variable2) * state_diff;
            residual_matrix(person, variable2) += observations(person, variable1) * state_diff;
          }
        }

        // Robbins-Monro update
        proposal_sd_pairwise_effects(variable1, variable2) =
          update_proposal_sd_with_robbins_monro (proposal_sd_pairwise_effects(variable1, variable2),
          log_acceptance_probability, target_acceptance_rate_mh, rm_lower_bound, rm_upper_bound,
          exp_neg_log_t_rm_adaptation_rate
        );
      }
    }
  }
}



/**
 * Function: update_indicator_interaction_pair_with_metropolis
 *
 * Performs a Metropolis-Hastings update of the edge inclusion indicators
 * for a predefined set of interaction pairs. For each pair (i, j),
 * this function proposes toggling the inclusion indicator (add or remove the edge),
 * and updates the interaction parameter and residual matrix accordingly.
 *
 * The proposal is evaluated using the pseudo-likelihood ratio and prior inclusion probabilities.
 *
 * Inputs:
 *  - interactions: current matrix of interaction weights.
 *  - thresholds: current matrix of main effect parameters.
 *  - indicator: matrix of edge inclusion indicators (updated in-place).
 *  - observations: matrix of observed scores per person and variable.
 *  - num_categories: vector of category counts per variable.
 *  - proposal_sd: matrix of proposal SDs for pairwise interactions.
 *  - interaction_scale: scale parameter for the Cauchy prior.
 *  - index: matrix of edge indices (E × 3): [index, i, j].
 *  - no_interactions: number of edge pairs (E).
 *  - no_persons: number of individuals in the data.
 *  - residual_matrix: matrix of linear predictors (updated on acceptance).
 *  - theta: matrix of inclusion probabilities under the edge prior.
 *  - is_ordinal_variable: logical vector indicating ordinal vs. Blume-Capel variables.
 *  - reference_category: vector of reference categories.
 *
 * Updates (in-place):
 *  - indicator
 *  - interactions
 *  - residual_matrix
 */
void update_indicator_interaction_pair_with_metropolis (
    arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    arma::umat& indicator,
    const arma::umat& observations,
    const arma::uvec& num_categories,
    const arma::mat& proposal_sd,
    const double interaction_scale,
    const arma::umat& index,
    const arma::uword no_interactions,
    const arma::uword no_persons,
    arma::mat& residual_matrix,
    const arma::mat& theta,
    const arma::uvec& is_ordinal_variable,
    const arma::uvec& reference_category
) {
  for (arma::uword cntr = 0; cntr < no_interactions; cntr++) {
    const arma::uword variable1 = index(cntr, 1);
    const arma::uword variable2 = index(cntr, 2);

    const double current_state = pairwise_effects(variable1, variable2);
    double proposed_state;

    bool proposing_addition = (indicator(variable1, variable2) == 0);
    if (proposing_addition) {
      proposed_state = R::rnorm(current_state, proposal_sd(variable1, variable2));
    } else {
      proposed_state = 0.0;
    }

    double log_acceptance_probability = log_pseudolikelihood_ratio_interaction(
      pairwise_effects, main_effects, observations, num_categories, no_persons,
      variable1, variable2, proposed_state, current_state,
      residual_matrix, is_ordinal_variable, reference_category
    );

    if (proposing_addition) {
      // Prior and proposal corrections for proposing an addition
      log_acceptance_probability += R::dcauchy(proposed_state, 0.0, interaction_scale, true);
      log_acceptance_probability -= R::dnorm(proposed_state, current_state, proposal_sd(variable1, variable2), true);
      log_acceptance_probability += std::log(theta(variable1, variable2)) - std::log(1.0 - theta(variable1, variable2));
    } else {
      // Prior and proposal corrections for proposing a removal
      log_acceptance_probability -= R::dcauchy(current_state, 0.0, interaction_scale, true);
      log_acceptance_probability += R::dnorm(current_state, proposed_state, proposal_sd(variable1, variable2), true);
      log_acceptance_probability -= std::log(theta(variable1, variable2)) - std::log(1.0 - theta(variable1, variable2));
    }

    if (std::log(R::unif_rand()) < log_acceptance_probability) {
      indicator(variable1, variable2) = 1 - indicator(variable1, variable2);
      indicator(variable2, variable1) = indicator(variable1, variable2);

      pairwise_effects(variable1, variable2) = proposed_state;
      pairwise_effects(variable2, variable1) = proposed_state;

      const double state_diff = proposed_state - current_state;

      for (arma::uword person = 0; person < no_persons; person++) {
        residual_matrix(person, variable1) += observations(person, variable2) * state_diff;
        residual_matrix(person, variable2) += observations(person, variable1) * state_diff;
      }
    }
  }
}



/**
 * Function: gibbs_update_step_for_graphical_model_parameters
 *
 * Performs a single Gibbs update step for all model parameters in the graphical model.
 *
 * This includes:
 * - Adaptive Metropolis-Hastings updates for pairwise effects.
 * - MALA or Metropolis updates for main effects (ordinal and Blume-Capel).
 * - Edge selection via indicator sampling (if enabled).
 *
 * Step sizes are tuned using Robbins-Monro or dual averaging. All updates are done in-place.
 *
 * Inputs:
 *  - observations: matrix of category scores.
 *  - num_categories: number of categories for each variable.
 *  - interaction_scale: scale parameter for the Cauchy prior.
 *  - proposal_sd_pairwise_effects: matrix of proposal SDs for pairwise effects.
 *  - proposal_sd_main_effects: matrix of proposal SDs for main effects.
 *  - update_index: index matrix of interaction pairs to update.
 *  - num_obs_categories: matrix of observed category counts.
 *  - sufficient_blume_capel: sufficient stats for Blume-Capel updates.
 *  - threshold_alpha, threshold_beta: prior hyperparameters for main effects.
 *  - num_persons: number of observations.
 *  - num_variables: number of variables.
 *  - num_pairwise: number of interaction pairs.
 *  - num_main: number of main effect parameters.
 *  - max_categories: maximum number of response categories.
 *  - theta: matrix of edge probabilities under the edge prior.
 *  - rm_decay_rate: Robbins-Monro decay rate.
 *  - target_acceptance_rate_mh: target MH acceptance rate (e.g., 0.234).
 *  - rm_lower_bound, rm_upper_bound: bounds for Robbins-Monro proposals.
 *  - is_ordinal_variable: indicator for ordinal vs. Blume-Capel variables.
 *  - reference_category: vector of reference categories.
 *  - edge_selection: whether edge selection is active this iteration.
 *  - step_size_mala: current MALA step size (updated via dual averaging).
 *  - iteration: current iteration number.
 *  - dual_averaging_state: state vector for dual averaging [log_eps, log_eps_avg, H_bar].
 *  - target_log_step_size: target log step size for MALA.
 *  - total_burnin: total number of burn-in iterations.
 *  - use_mala_for_main_effects: whether MALA is enabled for main effect updates.
 *
 * Updates (in-place):
 *  - inclusion_indicator
 *  - pairwise_effects
 *  - main_effects
 *  - residual_matrix
 *  - proposal_sd_main_effects
 *  - proposal_sd_pairwise_effects
 *  - step_size_mala (if MALA is active)
 */
void gibbs_update_step_for_graphical_model_parameters (
    const arma::umat& observations,
    const arma::uvec& num_categories,
    const double interaction_scale,
    arma::mat& proposal_sd_pairwise_effects,
    arma::mat& proposal_sd_main_effects,
    const arma::umat& update_index,
    const arma::umat& num_obs_categories,
    const arma::umat& sufficient_blume_capel,
    const double threshold_alpha,
    const double threshold_beta,
    const arma::uword num_persons,
    const arma::uword num_variables,
    const arma::uword num_pairwise,
    const arma::uword num_main,
    const arma::uword max_categories,
    arma::umat& inclusion_indicator,
    arma::mat& pairwise_effects,
    arma::mat& main_effects,
    arma::mat& residual_matrix,
    const arma::mat& theta,
    const double rm_decay_rate,
    const double target_acceptance_rate_mh,
    const double rm_lower_bound,
    const double rm_upper_bound,
    const arma::uvec& is_ordinal_variable,
    const arma::uvec& reference_category,
    const bool edge_selection,
    double& step_size_mala,
    const arma::uword iteration,
    arma::vec& dual_averaging_state,
    const double target_log_step_size,
    const arma::uword total_burnin,
    const bool use_mala_for_main_effects
) {

  const double exp_neg_log_t_rm_adaptation_rate =
    std::exp(-std::log(static_cast<double>(iteration)) * rm_decay_rate);

  //Between model move (update edge indicators and interaction parameters)
  if (edge_selection) {
    update_indicator_interaction_pair_with_metropolis (
        pairwise_effects, main_effects, inclusion_indicator, observations,
        num_categories, proposal_sd_pairwise_effects, interaction_scale,
        update_index, num_pairwise, num_persons, residual_matrix,
        theta, is_ordinal_variable, reference_category
      );
  }

  //Within model move (update interaction parameters)
  update_interactions_with_adaptive_metropolis (
    pairwise_effects, main_effects, inclusion_indicator, observations,
    num_categories, proposal_sd_pairwise_effects, interaction_scale,
    num_persons, num_variables, residual_matrix,
    exp_neg_log_t_rm_adaptation_rate, target_acceptance_rate_mh, rm_lower_bound,
    rm_upper_bound, is_ordinal_variable, reference_category
  );

  //Update threshold parameters
  if (use_mala_for_main_effects) {
    update_thresholds_with_adaptive_mala(
      main_effects, step_size_mala, residual_matrix, num_categories,
      num_obs_categories, sufficient_blume_capel, reference_category,
      is_ordinal_variable, iteration, total_burnin, dual_averaging_state,
      target_log_step_size, threshold_alpha, threshold_beta
    );
  } else {
    for (arma::uword variable = 0; variable < num_variables; ++variable) {
      if (is_ordinal_variable(variable)) {
        update_regular_thresholds_with_metropolis (
          main_effects, observations, num_categories, num_obs_categories,
          num_persons, variable, threshold_alpha, threshold_beta,
          residual_matrix
        );
      } else {
        update_blumecapel_thresholds_with_adaptive_metropolis(
          main_effects, observations, num_categories, sufficient_blume_capel,
          num_persons, variable, reference_category, threshold_alpha,
          threshold_beta, residual_matrix, proposal_sd_main_effects,
          exp_neg_log_t_rm_adaptation_rate, target_acceptance_rate_mh,
          rm_lower_bound, rm_upper_bound
        );
      }
    }
  }
}



/**
 * Function: run_gibbs_sampler_for_bgm
 *
 * Runs the Gibbs sampling algorithm for a Markov random field model with
 * ordinal and/or Blume-Capel variables, including optional edge selection.
 *
 * The sampler alternates between updating:
 * - Main effects (threshold parameters) using MALA or Metropolis
 * - Pairwise interactions using adaptive Metropolis-Hastings
 * - Inclusion indicators (edges) if edge selection is enabled
 * - Imputation of missing values if applicable
 *
 * Adaptive tuning of proposal distributions is performed using Robbins-Monro or dual averaging.
 *
 * Inputs:
 *  - observations: data matrix with categorical scores.
 *  - num_categories: number of categories per variable.
 *  - interaction_scale: scale of Cauchy prior for pairwise interactions.
 *  - proposal_sd: initial proposal SDs for pairwise effects.
 *  - proposal_sd_blumecapel: initial proposal SDs for Blume-Capel main effects.
 *  - edge_prior: prior type ("Bernoulli", "Beta-Bernoulli", "Stochastic-Block").
 *  - theta: matrix of inclusion probabilities under the edge prior.
 *  - beta_bernoulli_alpha, beta_bernoulli_beta: shape parameters for Beta-Bernoulli prior.
 *  - dirichlet_alpha, lambda: parameters for the SBM prior (if used).
 *  - Index: list of interaction indices (row id, i, j).
 *  - iter, burnin: number of sampling and burn-in iterations.
 *  - num_obs_categories: counts of observed categories.
 *  - sufficient_blume_capel: sufficient stats for Blume-Capel updates.
 *  - threshold_alpha, threshold_beta: prior hyperparameters for main effects.
 *  - na_impute: whether to impute missing values.
 *  - missing_index: matrix of (i, j) missing locations.
 *  - is_ordinal_variable: indicator for ordinal vs. Blume-Capel variables.
 *  - reference_category: reference categories for Blume-Capel variables.
 *  - save_main, save_pairwise, save_indicator: flags to store MCMC samples.
 *  - display_progress: whether to show progress bar.
 *  - use_mala_for_main_effects: use MALA for threshold updates.
 *
 * Returns:
 *  - A List containing posterior mean estimates:
 *      * "main": average thresholds
 *      * "pairwise": average interaction weights
 *      * "inclusion_indicator": average edge inclusion matrix
 *  - If sampling is enabled:
 *      * "main_samples", "pairwise_samples", "inclusion_indicator_samples"
 *  - If edge_prior is SBM:
 *      * "allocations": sampled cluster assignments over iterations
 */
// [[Rcpp::export]]
List run_gibbs_sampler_for_bgm (
    arma::umat& observations,
    const arma::uvec& num_categories,
    const double interaction_scale,
    const String& edge_prior,
    arma::mat& theta,
    const double beta_bernoulli_alpha,
    const double beta_bernoulli_beta,
    const double dirichlet_alpha,
    const double lambda,
    const arma::umat& Index,
    const arma::uword iter,
    const arma::uword burnin,
    arma::umat& num_obs_categories,
    arma::umat& sufficient_blume_capel,
    const double threshold_alpha,
    const double threshold_beta,
    const bool na_impute,
    const arma::umat& missing_index,
    const arma::uvec& is_ordinal_variable,
    const arma::uvec& reference_category,
    const bool save_main = false,
    const bool save_pairwise = false,
    const bool save_indicator = false,
    const bool display_progress = false,
    bool edge_selection = true,
    bool use_mala_for_main_effects = false
) {

  // Step 1: Initialize parameters, proposal settings, and storage
  const arma::uword num_variables = observations.n_cols;
  const arma::uword num_persons = observations.n_rows;
  const arma::uword max_num_categories = num_categories.max();
  const arma::uword num_pairwise = Index.n_rows;

  // Model parameters
  arma::mat main_effects(num_variables, max_num_categories, arma::fill::zeros);
  arma::mat pairwise_effects(num_variables, num_variables, arma::fill::zeros);
  arma::umat inclusion_indicator(num_variables, num_variables, arma::fill::ones);

  // Posterior mean accumulation
  arma::mat posterior_mean_main(num_variables, max_num_categories, arma::fill::zeros);
  arma::mat posterior_mean_pairwise(num_variables, num_variables, arma::fill::zeros);
  arma::mat posterior_mean_indicator(num_variables, num_variables, arma::fill::zeros);

  // Residual matrix
  arma::mat residual_matrix(num_persons, num_variables, arma::fill::zeros);

  // Optional sample storage
  const arma::uword num_main = count_num_main_effects(num_categories, is_ordinal_variable);
  arma::mat* main_effect_samples = nullptr;
  arma::mat* pairwise_effect_samples = nullptr;
  arma::umat* indicator_samples = nullptr;

  if (save_main) {
    main_effect_samples = new arma::mat(iter, num_main);
  }
  if (save_pairwise) {
    pairwise_effect_samples = new arma::mat(iter, num_pairwise);
  }
  if (save_indicator) {
    indicator_samples = new arma::umat(iter, num_pairwise);  // Edge indicators
  }

  // Adaptive Metropolis proposal standard deviations
  arma::mat proposal_sd_blumecapel(num_main, 2, arma::fill::ones);
  arma::mat proposal_sd_pairwise_effects(num_variables, num_variables, arma::fill::ones);

  // Adaptive step size and proposal tuning
  const double initial_step_size = 0.01;
  double step_size_mala = initial_step_size;
  const double target_log_step_size = std::log(10.0 * initial_step_size);
  arma::vec dual_averaging_state(3, arma::fill::zeros);  // [log_eps, log_eps_avg, H_bar]

  const double target_acceptance_rate_mh = 0.234;

  const double decay_rate_robins_monro = 0.75;
  const double rm_lower_bound = 1.0 / static_cast<double>(num_persons);
  const double rm_upper_bound = 2.0;

  // Edge update indexing
  arma::uvec v = arma::regspace<arma::uvec>(0, num_pairwise - 1);
  arma::uvec order(num_pairwise);
  arma::umat index(num_pairwise, 3);

  // SBM model setup (only if used)
  arma::uvec K_values;
  arma::uvec cluster_allocations(num_variables);
  arma::mat cluster_prob(1, 1);
  arma::vec log_Vn(1);
  arma::mat out_allocations(iter, num_variables);

  if (edge_prior == "Stochastic-Block") {
    // Initialize random cluster assignments
    cluster_allocations[0] = 0;
    cluster_allocations[1] = 1;
    for (arma::uword i = 2; i < num_variables; i++) {
      cluster_allocations[i] = (R::unif_rand() > 0.5) ? 1 : 0;
    }

    // Compute initial cluster probabilities
    cluster_prob = block_probs_mfm_sbm(
      cluster_allocations, inclusion_indicator,
      num_variables, beta_bernoulli_alpha, beta_bernoulli_beta
    );

    for (arma::uword i = 0; i < num_variables - 1; i++) {
      for (arma::uword j = i + 1; j < num_variables; j++) {
        theta(i, j) = cluster_prob(cluster_allocations[i], cluster_allocations[j]);
        theta(j, i) = theta(i, j);
      }
    }

    log_Vn = compute_Vn_mfm_sbm(
      num_variables, dirichlet_alpha, num_variables + 10, lambda
    );
  }

  // Progress bar
  bool enable_edge_selection = edge_selection;
  const arma::uword total_burnin = burnin * (enable_edge_selection ? 2 : 1);
  edge_selection = false;  // Disabled during first burn-in phase

  const arma::uword total_iter = total_burnin + iter;
  Progress p(total_iter, display_progress);


  // Step 2: Gibbs sampling loop
  for (arma::uword iteration = 0; iteration < total_iter; iteration++) {
    if (Progress::check_abort()) {
      return List::create(
        Named("main") = posterior_mean_main,
        Named("pairwise") = posterior_mean_pairwise,
        Named("inclusion_indicator") = posterior_mean_indicator
      );
    }
    p.increment();

    // Re-enable edge selection after first burn-in stage
    if (enable_edge_selection && iteration == burnin) {
      edge_selection = true;
    }

    // Shuffle edge update order
    order = arma::randperm(num_pairwise);
    for (arma::uword i = 0; i < num_pairwise; i++) {
      index.row(i) = Index.row(order(i));
    }

    // Missing value imputation (if enabled)
    if (na_impute) {
      impute_missing_values_for_graphical_model (
        pairwise_effects, main_effects, observations, num_obs_categories,
        sufficient_blume_capel, num_categories, residual_matrix,
        missing_index, is_ordinal_variable, reference_category
      );
    }

    // Gibbs update step for model parameters
    gibbs_update_step_for_graphical_model_parameters (
      observations, num_categories, interaction_scale,
      proposal_sd_pairwise_effects, proposal_sd_blumecapel, index,
      num_obs_categories, sufficient_blume_capel, threshold_alpha,
      threshold_beta, num_persons, num_variables, num_pairwise, num_main,
      max_num_categories, inclusion_indicator, pairwise_effects, main_effects,
      residual_matrix, theta, decay_rate_robins_monro, target_acceptance_rate_mh,
      rm_lower_bound, rm_upper_bound, is_ordinal_variable, reference_category,
      edge_selection, step_size_mala, iteration, dual_averaging_state,
      target_log_step_size, total_burnin, use_mala_for_main_effects
    );


    // Update edge inclusion probabilities if edge selection is enabled
    if (edge_selection) {
      if (edge_prior == "Beta-Bernoulli") {
        arma::uword num_edges_included = 0;
        for (arma::uword i = 0; i < num_variables - 1; i++) {
          for (arma::uword j = i + 1; j < num_variables; j++) {
            num_edges_included += inclusion_indicator(i, j);
          }
        }

        double prob = R::rbeta(
          beta_bernoulli_alpha + num_edges_included,
          beta_bernoulli_beta + num_pairwise - num_edges_included
        );

        for (arma::uword i = 0; i < num_variables - 1; i++) {
          for (arma::uword j = i + 1; j < num_variables; j++) {
            theta(i, j) = theta(j, i) = prob;
          }
        }

      } else if (edge_prior == "Stochastic-Block") {
        cluster_allocations = block_allocations_mfm_sbm(
          cluster_allocations, num_variables, log_Vn, cluster_prob,
          inclusion_indicator, dirichlet_alpha,
          beta_bernoulli_alpha, beta_bernoulli_beta
        );

        cluster_prob = block_probs_mfm_sbm(
          cluster_allocations, inclusion_indicator, num_variables,
          beta_bernoulli_alpha, beta_bernoulli_beta
        );

        for (arma::uword i = 0; i < num_variables - 1; i++) {
          for (arma::uword j = i + 1; j < num_variables; j++) {
            theta(i, j) = theta(j, i) = cluster_prob(cluster_allocations[i], cluster_allocations[j]);
          }
        }
      }
    }

    // Save samples or update running posterior means after burn-in
    if (iteration >= total_burnin) {
      arma::uword iter_adj = iteration - total_burnin + 1;

      // Update running posterior mean: main effects
      posterior_mean_main = (posterior_mean_main * (iter_adj - 1) + main_effects) / static_cast<double>(iter_adj);

      // Update running posterior mean: pairwise effects
      posterior_mean_pairwise = (posterior_mean_pairwise * (iter_adj - 1) + pairwise_effects) / static_cast<double>(iter_adj);

      if(edge_selection) {
        // Update running posterior mean: inclusion indicator
        posterior_mean_indicator = (posterior_mean_indicator * (iter_adj - 1) +
          arma::conv_to<arma::mat>::from(inclusion_indicator)) / static_cast<double>(iter_adj);
      }

      arma::uword sample_index = iteration - total_burnin;

      if (save_main) {
        arma::vec vectorized_main = vectorize_thresholds(main_effects, num_categories, is_ordinal_variable);
        (*main_effect_samples).row(sample_index) = vectorized_main.t();
      }

      if (save_pairwise) {
        arma::vec vectorized_pairwise(num_pairwise);
        for (arma::uword i = 0; i < num_pairwise; ++i) {
          vectorized_pairwise(i) = pairwise_effects(Index(i, 1), Index(i, 2));
        }
        (*pairwise_effect_samples).row(sample_index) = vectorized_pairwise.t();
      }

      if (save_indicator) {
        arma::uvec vectorized_indicator(num_pairwise);
        for (arma::uword i = 0; i < num_pairwise; ++i) {
          vectorized_indicator(i) = inclusion_indicator(Index(i, 1), Index(i, 2));
        }
        (*indicator_samples).row(sample_index) = vectorized_indicator.t();
      }

      // Save SBM cluster allocations
      if (edge_prior == "Stochastic-Block") {
        for (arma::uword j = 0; j < num_variables; ++j)
          out_allocations(sample_index, j) = cluster_allocations[j] + 1;
      }
    }
  }


  // Step 3: Return results
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
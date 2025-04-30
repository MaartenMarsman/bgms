// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppProgress)]]
#include <Rcpp.h>
#include "gibbs_functions_edge_prior.h"
#include <progress.hpp>
#include <progress_bar.hpp>
using namespace Rcpp;


/**
 * Function: dual_averaging_update
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
 *   - log_step_size_target: Desired log step size target (typically log(10 * initial_step_size)).
 *   - stabilization_offset: Constant to reduce adaptation sensitivity early on (e.g., 10).
 *
 * Outputs:
 *   - Updates the `state` vector in-place.
 *
 * Usage:
 *   - Used inside the adaptive Metropolis-adjusted Langevin algorithm (MALA) for updating threshold parameters.
 *   - Called during the burn-in phase inside `adamala_thresholds()` to adaptively tune the step size.
 */
inline void dual_averaging_update (
    const double acceptance_probability,
    const arma::uword iteration,
    arma::vec& state,
    const double log_step_size_target,
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
  log_step_size = log_step_size_target - std::sqrt(iter_double) / gamma * acceptance_error_avg;

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
 *     log(step_size) ← log(step_size) + (accept_prob - target) * t^{-phi}
 *
 * This ensures diminishing adaptation over time and avoids instability
 * in later iterations. The target acceptance rate (0.574) is optimal
 * for MALA under typical assumptions.
 *
 * Parameters:
 *  - acceptance_probability: Observed acceptance rate at the current iteration.
 *  - iteration: Current iteration count (starts at 1).
 *  - step_size: The current step size (updated in-place).
 *
 * Notes:
 *  - phi is fixed at 0.75 (typical value for stable decay).
 *  - This version assumes the step size is always positive.
 */
inline void robbins_monro_update(
    const double acceptance_probability,
    const arma::uword iteration,
    double& step_size
) {
  const double target_acceptance = 0.574;
  const double phi = 0.75;
  const double error = acceptance_probability - target_acceptance;
  const double decay = std::pow(static_cast<double>(iteration), -phi);

  double log_step_size = std::log(step_size);
  log_step_size += error * decay;

  step_size = std::exp(log_step_size);
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
inline arma::uword count_threshold_parameters(
    const arma::uvec& num_categories,
    const arma::uvec& is_ordinal_variable
) {
  arma::uword n_params = 0;
  for (arma::uword variable = 0; variable < num_categories.n_elem; variable++) {
    if (is_ordinal_variable(variable)) {
      n_params += num_categories(variable);
    } else {
      n_params += 2; // linear + quadratic for Blume-Capel
    }
  }
  return n_params;
}


/**
 * Function: flatten_thresholds
 * Purpose:
 *   Converts a matrix of threshold parameters into a flat vector for optimization,
 *   respecting the structure of ordinal vs. Blume-Capel variables.
 *
 * Inputs:
 *   - thresholds: Matrix of threshold parameters [variables × categories].
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
arma::vec flatten_thresholds(
    const arma::mat& thresholds,
    const arma::uvec& num_categories,
    const arma::uvec& is_ordinal_variable
) {
  const arma::uword num_parameters = count_threshold_parameters(num_categories, is_ordinal_variable);
  arma::vec flattened(num_parameters);

  arma::uword offset = 0;

  for (arma::uword var = 0; var < thresholds.n_rows; var++) {
    const arma::uword num_cats = is_ordinal_variable(var) ? num_categories(var) : 2;

    for (arma::uword j = 0; j < num_cats; j++) {
      flattened(offset++) = thresholds(var, j);
    }
  }

  return flattened;
}


/**
 * Function: unflatten_thresholds
 * Purpose:
 *   Reconstructs a threshold matrix from a flat vector of parameters,
 *   reversing the operation of `flatten_thresholds`.
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
arma::mat unflatten_thresholds(
    const arma::vec& flat_vector,
    const arma::uvec& num_categories,
    const arma::uvec& is_ordinal_variable
) {
  const arma::uword num_variables = num_categories.n_elem;
  const arma::uword max_categories = num_categories.max();

  arma::mat thresholds(num_variables, max_categories); // May contain unused entries

  arma::uword offset = 0;

  for (arma::uword var = 0; var < num_variables; var++) {
    const arma::uword num_cats = is_ordinal_variable(var) ? num_categories(var) : 2;

    for (arma::uword j = 0; j < num_cats; j++) {
      thresholds(var, j) = flat_vector(offset++);
    }
  }

  return thresholds;
}


/**
 * Function: impute_missing_data
 * Purpose:
 *   Imputes missing values in a graphical model with ordinal and Blume-Capel variables,
 *   using their full conditional pseudo-posterior distributions.
 *
 * Inputs:
 *   - interactions: Matrix of pairwise interaction weights between variables.
 *   - thresholds: Threshold parameters for each variable (linear and quadratic).
 *   - observations: Matrix of observed values (individuals × variables); updated in-place.
 *   - num_obs_categories: Category count matrix for ordinal variables [category × variable]; updated in-place.
 *   - sufficient_blume_capel: Sufficient statistics for BC variables [2 × variable]; updated in-place.
 *   - num_categories: Vector of the number of categories for each variable.
 *   - residual_matrix: Linear predictors excluding the current variable; updated in-place.
 *   - missing_index: Matrix of (row, col) indices of missing values.
 *   - is_ordinal_variable: Indicator for whether a variable is ordinal (1) or Blume-Capel (0).
 *   - reference_category: Reference (centering) category for each Blume-Capel variable.
 *
 * Outputs:
 *   - A List with updated:
 *       - `observations`: Updated observation matrix.
 *       - `num_obs_categories`: Updated category counts (ordinal only).
 *       - `sufficient_blume_capel`: Updated sufficient statistics (BC only).
 *       - `residual_matrix`: Updated pseudo-likelihood residuals.
 *
 * Usage:
 *   - Called within the Gibbs sampler during missing data imputation steps.
 */
List impute_missing_data (
    const arma::mat& interactions,
    const arma::mat& thresholds,
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
        const double exponent = thresholds(variable, cat) + score * rest_score;
        cumsum += std::exp(exponent);
        category_probabilities[score] = cumsum;
      }
    } else {
      // Compute probabilities for Blume-Capel variable
      const arma::uword ref = reference_category(variable);

      cumsum = std::exp(thresholds(variable, 1) * ref * ref);
      category_probabilities[0] = cumsum;

      for (arma::uword cat = 0; cat < num_cats; cat++) {
        const arma::uword score = cat + 1;
        const arma::sword centered = static_cast<arma::sword>(score) - static_cast<arma::sword>(ref);
        const double exponent =
          thresholds(variable, 0) * score +
          thresholds(variable, 1) * centered * centered +
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
        const double delta_score = (static_cast<double>(new_value) - old_value) * interactions(v, variable);
        residual_matrix(person, v) += delta_score;
      }
    }
  }

  return List::create(
    Named("observations") = observations,
    Named("num_obs_categories") = num_obs_categories,
    Named("sufficient_blume_capel") = sufficient_blume_capel,
    Named("residual_matrix") = residual_matrix
  );
}


/**
 * Function: log_pseudoposterior_thresholds
 * Purpose:
 *   Computes the log pseudo-posterior for threshold parameters in a graphical model
 *   combining contributions from the data likelihood and logistic-Beta priors.
 *   Supports both regular ordinal and Blume-Capel variables.
 *
 * Inputs:
 *   - thresholds: Matrix of threshold parameters [variables × categories].
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
 *   - Called inside `adamala_thresholds()` and Metropolis updates for thresholds.
 */
double log_pseudoposterior_thresholds (
    const arma::mat& thresholds,
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
        const double theta = thresholds(variable, cat);
        log_posterior += theta * (num_obs_categories(cat + 1, variable) + threshold_alpha);
        log_posterior -= std::log1p(std::exp(theta)) * (threshold_alpha + threshold_beta);
      }

      for (arma::uword person = 0; person < num_persons; person++) {
        const double rest_score = residual_matrix(person, variable);
        const double bound = num_cats * rest_score;

        double denominator = std::exp(-bound);
        for (arma::uword cat = 0; cat < num_cats; cat++) {
          const double exponent = thresholds(variable, cat) + (cat + 1) * rest_score - bound;
          denominator += std::exp(exponent);
        }

        log_posterior -= bound + std::log(denominator);
      }

    } else {
      // Blume-Capel variable
      const double theta_lin = thresholds(variable, 0);
      const double theta_quad = thresholds(variable, 1);

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
 *  - thresholds: [variables × categories] matrix of current threshold parameters.
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
arma::vec gradient_thresholds_pseudoposterior(
    const arma::mat& thresholds,
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
  const arma::uword num_parameters = count_threshold_parameters(num_categories, is_ordinal_variable);

  arma::vec gradient(num_parameters, arma::fill::zeros);
  arma::uword offset = 0;

  for (arma::uword variable = 0; variable < num_variables; variable++) {
    const arma::uword num_cats = num_categories(variable);

    if (is_ordinal_variable(variable)) {
      // Ordinal variable
      const double max_threshold = thresholds.row(variable).max();

      for (arma::uword cat = 0; cat < num_cats; cat++) {
        gradient(offset + cat) = num_obs_categories(cat + 1, variable);
      }

      for (arma::uword person = 0; person < num_persons; person++) {
        const double rest_score = residual_matrix(person, variable);
        const double bound = max_threshold + num_cats * rest_score;

        double denom = std::exp(-bound);
        arma::vec numerators(num_cats, arma::fill::zeros);

        for (arma::uword cat = 0; cat < num_cats; cat++) {
          const double exponent = thresholds(variable, cat) + (cat + 1) * rest_score - bound;
          numerators(cat) = std::exp(exponent);
          denom += numerators(cat);
        }

        for (arma::uword cat = 0; cat < num_cats; cat++) {
          gradient(offset + cat) -= numerators(cat) / denom;
        }
      }

      for (arma::uword cat = 0; cat < num_cats; cat++) {
        const double theta = thresholds(variable, cat);
        const double p = 1.0 / (1.0 + std::exp(-theta));
        gradient(offset + cat) += threshold_alpha - (threshold_alpha + threshold_beta) * p;
      }

      offset += num_cats;

    } else {
      // Blume-Capel variable
      const int ref = static_cast<int>(reference_category(variable));
      const double theta_lin = thresholds(variable, 0);
      const double theta_quad = thresholds(variable, 1);

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
        const double theta = thresholds(variable, i);
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
 *  - thresholds: [variables × max_categories] matrix of current thresholds (updated in-place).
 *  - step_size: MALA step size (updated in-place).
 *  - residual_matrix: [persons × variables] matrix of residual scores.
 *  - num_categories: Number of categories per variable.
 *  - num_obs_categories: Count of observed categories [category × variable].
 *  - sufficient_blume_capel: Sufficient stats for BC variables [2 × variable].
 *  - reference_category: Reference category for centering.
 *  - is_ordinal_variable: Indicator for ordinal (1) or Blume-Capel (0) per variable.
 *  - t: Current iteration.
 *  - burnin: Number of burn-in iterations.
 *  - dual_averaging_state: [log_eps, log_eps_bar, H_bar] vector (updated in-place).
 *  - log_step_size_target: Target log step size (typically log(10 * ε₀)).
 *  - threshold_alpha, threshold_beta: Hyperparameters for logistic-Beta prior.
 */
void adamala_thresholds (
    arma::mat& thresholds,
    double& step_size,
    const arma::mat& residual_matrix,
    const arma::uvec& num_categories,
    const arma::umat& num_obs_categories,
    const arma::umat& sufficient_blume_capel,
    const arma::uvec& reference_category,
    const arma::uvec& is_ordinal_variable,
    const arma::uword t,
    const arma::uword burnin,
    arma::vec& dual_averaging_state,
    const double log_step_size_target,
    const double threshold_alpha,
    const double threshold_beta
) {
  // Flatten parameters
  arma::vec flat_theta = flatten_thresholds(thresholds, num_categories, is_ordinal_variable);

  // Evaluate current gradient and log posterior
  arma::vec grad = gradient_thresholds_pseudoposterior(
    thresholds, residual_matrix, num_categories, num_obs_categories,
    sufficient_blume_capel, reference_category, is_ordinal_variable,
    threshold_alpha, threshold_beta
  );

  const double log_post_current = log_pseudoposterior_thresholds(
    thresholds, residual_matrix, num_categories, num_obs_categories,
    sufficient_blume_capel, reference_category, is_ordinal_variable,
    threshold_alpha, threshold_beta
  );

  // Propose new theta
  const double sqrt_step = std::sqrt(step_size);
  arma::vec proposal = flat_theta + 0.5 * step_size * grad + sqrt_step * arma::randn(flat_theta.n_elem);
  arma::mat proposed_thresholds = unflatten_thresholds(proposal, num_categories, is_ordinal_variable);

  // Evaluate proposal log posterior and gradient
  const double log_post_proposal = log_pseudoposterior_thresholds(
    proposed_thresholds, residual_matrix, num_categories, num_obs_categories,
    sufficient_blume_capel, reference_category, is_ordinal_variable,
    threshold_alpha, threshold_beta
  );

  arma::vec grad_proposal = gradient_thresholds_pseudoposterior(
    proposed_thresholds, residual_matrix, num_categories, num_obs_categories,
    sufficient_blume_capel, reference_category, is_ordinal_variable,
    threshold_alpha, threshold_beta
  );

  // Compute forward and backward proposal densities
  const arma::vec forward_mean = flat_theta + 0.5 * step_size * grad;
  const arma::vec backward_mean = proposal + 0.5 * step_size * grad_proposal;

  const double log_forward = -0.5 / step_size * arma::accu(arma::square(proposal - forward_mean));
  const double log_backward = -0.5 / step_size * arma::accu(arma::square(flat_theta - backward_mean));

  // Metropolis acceptance probability
  const double log_acceptance = log_post_proposal + log_backward - log_post_current - log_forward;

  // Accept/reject
  if (std::log(R::unif_rand()) < log_acceptance) {
    thresholds = proposed_thresholds;
    flat_theta = proposal;
  }

  // Compute acceptance probability (bounded at 1)
  const double accept_prob = std::min(1.0, std::exp(log_acceptance));

  // Adapt step size
  if (t <= burnin) {
    dual_averaging_update(
      accept_prob, t, dual_averaging_state,
      log_step_size_target, 10
    );
    step_size = std::exp(dual_averaging_state[1]);
  } else {
    robbins_monro_update(accept_prob, t - burnin, step_size);
  }
}


/**
 * Function: metropolis_thresholds_regular
 * Purpose:
 *   Performs a Metropolis-Hastings update for threshold parameters of ordinal variables.
 *   Each threshold is updated one at a time using a generalized beta-prime proposal,
 *   with acceptance determined by the pseudo-likelihood and logistic-Beta prior.
 *
 * Inputs:
 *   - thresholds: Matrix of threshold parameters [variables × categories]; updated in-place.
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
 *   - Updates `thresholds(variable, category)` in-place for each category of the specified variable.
 *
 * Usage:
 *   - Called during Gibbs updates when using non-gradient MH threshold proposals.
 */
void metropolis_thresholds_regular (
    arma::mat& thresholds,
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
    double current = thresholds(variable, category);
    double exp_current = std::exp(current);
    double c = (threshold_alpha + threshold_beta) / (1 + exp_current);

    // Compute proposal scaling constant `c`
    for (arma::uword person = 0; person < no_persons; person++) {
      double rest_score = residual_matrix(person, variable);
      double denom = 1.0;  // base of sum (g)
      double numer = std::exp((category + 1) * rest_score);  // q

      for (arma::uword cat = 0; cat < num_categories[variable]; cat++) {
        if (cat != category) {
          denom += std::exp(thresholds(variable, cat) + (cat + 1) * rest_score);
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
    double log_prob = 0.0;
    for (arma::uword person = 0; person < no_persons; person++) {
      log_prob += std::log(g[person] + q[person] * exp_current);
      log_prob -= std::log(g[person] + q[person] * exp_proposed);
    }

    // Add prior ratio (logistic-Beta)
    log_prob -= (threshold_alpha + threshold_beta) * std::log(1 + exp_proposed);
    log_prob += (threshold_alpha + threshold_beta) * std::log(1 + exp_current);

    // Add proposal ratio (generalized beta-prime)
    log_prob -= (a + b) * std::log(1 + c * exp_current);
    log_prob += (a + b) * std::log(1 + c * exp_proposed);

    // Metropolis step
    if (std::log(R::unif_rand()) < log_prob) {
      thresholds(variable, category) = proposed;
    }
  }
}



/**
 * Performs a Metropolis update of the Blume-Capel threshold parameters for a single variable.
 *
 * For each Blume-Capel variable, both the linear and quadratic threshold parameters
 * are updated in turn using a log pseudo-likelihood and logistic-Beta prior.
 *
 * Step sizes are adapted with Robbins-Monro based on the observed acceptance probability.
 *
 * Parameters:
 *  - thresholds: [variables × 2] matrix of threshold parameters (linear and quadratic).
 *  - observations: [persons × variables] (not used directly, passed for consistency).
 *  - num_categories: Vector of number of categories per variable.
 *  - sufficient_blume_capel: [2 × variables] matrix of sufficient statistics.
 *  - no_persons: Number of observations.
 *  - variable: Index of the variable being updated.
 *  - reference_category: Vector of reference categories per variable.
 *  - threshold_alpha, threshold_beta: Prior hyperparameters.
 *  - residual_matrix: [persons × variables] matrix of rest scores.
 *  - proposal_sd_blumecapel: [variables × 2] proposal SD matrix (updated in-place).
 *  - phi: Robbins-Monro decay rate.
 *  - target_ar: Target acceptance rate.
 *  - t: Current iteration (starting at 1).
 *  - epsilon_lo, epsilon_hi: Bounds on proposal SD.
 */
void metropolis_thresholds_blumecapel(
    arma::mat& thresholds,
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
    const double phi,
    const double target_ar,
    const arma::uword t,
    const double epsilon_lo,
    const double epsilon_hi
) {
  const arma::uword num_cats = num_categories(variable);
  const int ref = static_cast<int>(reference_category(variable));

  // Helper lambda to update either the linear (0) or quadratic (1) threshold
  auto update_parameter = [&](arma::uword param_index) {
    const double current = thresholds(variable, param_index);
    const double proposed = R::rnorm(current, proposal_sd_blumecapel(variable, param_index));
    const double diff = proposed - current;

    arma::vec numer_current(num_cats + 1), numer_proposed(num_cats + 1);

    for (arma::uword cat = 0; cat <= num_cats; cat++) {
      const int centered = static_cast<int>(cat) - ref;
      if (param_index == 0) { // linear update
        double quad_term = thresholds(variable, 1) * (centered * centered);
        numer_current(cat) = current * cat + quad_term;
        numer_proposed(cat) = proposed * cat + quad_term;
      } else { // quadratic update
        double lin_term = thresholds(variable, 0) * cat;
        numer_current(cat) = current * (centered * centered) + lin_term;
        numer_proposed(cat) = proposed * (centered * centered) + lin_term;
      }
    }

    const double max_curr = numer_current.max();
    const double max_prop = numer_proposed.max();
    const double lbound = (max_curr > 0 || max_prop > 0) ? std::max(max_curr, max_prop) : 0.0;

    double log_prob = threshold_alpha * diff +
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

      log_prob += std::log(denom_curr) - std::log(denom_prop);
    }

    // Prior ratio
    log_prob += (threshold_alpha + threshold_beta) * (
      std::log1p(std::exp(current)) - std::log1p(std::exp(proposed))
    );

    // Metropolis accept/reject
    const double logu = std::log(R::unif_rand());
    if (logu < log_prob) {
      thresholds(variable, param_index) = proposed;
    }

    // Robbins-Monro adaptation
    double acc_prob = (log_prob > 0.0) ? 1.0 : std::exp(log_prob);
    double updated_sd = proposal_sd_blumecapel(variable, param_index) +
      (acc_prob - target_ar) * std::pow(static_cast<double>(t), -phi);

    if (std::isnan(updated_sd)) updated_sd = 1.0;

    proposal_sd_blumecapel(variable, param_index) =
      std::clamp(updated_sd, epsilon_lo, epsilon_hi);
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
 *  - interactions: [V × V] matrix of interaction parameters.
 *  - thresholds: [V × max_categories] matrix of threshold parameters.
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
double log_pseudolikelihood_ratio_interaction(
    const arma::mat& interactions,
    const arma::mat& thresholds,
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
          const double exponent = thresholds(variable1, cat) + (cat + 1) * rest_score;
          denom_current += std::exp(exponent + (cat + 1) * score2 * current_state - bound);
          denom_proposed += std::exp(exponent + (cat + 1) * score2 * proposed_state - bound);
        }
      } else {
        const arma::uword ref = reference_category(variable1);
        for (arma::uword cat = 0; cat <= num_categories(variable1); cat++) {
          const int centered = static_cast<int>(cat) - static_cast<int>(ref);
          const double exponent =
            thresholds(variable1, 0) * cat +
            thresholds(variable1, 1) * centered * centered +
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
          const double exponent = thresholds(variable2, cat) + (cat + 1) * rest_score;
          denom_current += std::exp(exponent + (cat + 1) * score1 * current_state - bound);
          denom_proposed += std::exp(exponent + (cat + 1) * score1 * proposed_state - bound);
        }
      } else {
        const arma::uword ref = reference_category(variable2);
        for (arma::uword cat = 0; cat <= num_categories(variable2); cat++) {
          const int centered = static_cast<int>(cat) - static_cast<int>(ref);
          const double exponent =
            thresholds(variable2, 0) * cat +
            thresholds(variable2, 1) * centered * centered +
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
 * Performs Metropolis-Hastings updates for all active pairwise interaction parameters.
 *
 * For each pair (i, j) where `indicator(i, j) == 1`, a proposal is drawn and
 * accepted or rejected based on the pseudo-likelihood ratio and Cauchy prior.
 * The proposal SD is adapted using Robbins-Monro updates.
 *
 * Parameters:
 *  - interactions: [V × V] symmetric matrix of current interaction values (updated in-place).
 *  - thresholds: [V × max_categories] matrix of threshold parameters.
 *  - indicator: [V × V] binary matrix indicating active interactions.
 *  - observations: [N × V] matrix of observed category scores.
 *  - num_categories: Vector of category counts per variable.
 *  - proposal_sd: [V × V] matrix of proposal standard deviations (updated).
 *  - interaction_scale: Scale parameter for the Cauchy prior.
 *  - no_persons: Number of observations.
 *  - no_variables: Number of variables (V).
 *  - residual_matrix: [N × V] matrix of current residual scores (updated on acceptance).
 *  - phi: Robbins-Monro adaptation decay exponent.
 *  - target_ar: Target acceptance rate (e.g., 0.234).
 *  - t: Current iteration (starting from 1).
 *  - epsilon_lo, epsilon_hi: Bounds on adapted proposal SDs.
 *  - is_ordinal_variable: Vector indicating ordinal (1) vs. Blume-Capel (0).
 *  - reference_category: Reference category per variable for BC modeling.
 */
void metropolis_interactions(
    arma::mat& interactions,
    const arma::mat& thresholds,
    const arma::umat& indicator,
    const arma::umat& observations,
    const arma::uvec& num_categories,
    arma::mat& proposal_sd,
    const double interaction_scale,
    const arma::uword no_persons,
    const arma::uword no_variables,
    arma::mat& residual_matrix,
    const double phi,
    const double target_ar,
    const arma::uword t,
    const double epsilon_lo,
    const double epsilon_hi,
    const arma::uvec& is_ordinal_variable,
    const arma::uvec& reference_category
) {
  for (arma::uword variable1 = 0; variable1 < no_variables - 1; variable1++) {
    for (arma::uword variable2 = variable1 + 1; variable2 < no_variables; variable2++) {
      if (indicator(variable1, variable2) == 1) {
        const double current_state = interactions(variable1, variable2);
        const double proposed_state = R::rnorm(current_state, proposal_sd(variable1, variable2));

        double log_prob = log_pseudolikelihood_ratio_interaction(
          interactions, thresholds, observations, num_categories, no_persons,
          variable1, variable2, proposed_state, current_state,
          residual_matrix, is_ordinal_variable, reference_category
        );

        // Add log-prior difference (Cauchy prior)
        log_prob += R::dcauchy(proposed_state, 0.0, interaction_scale, true);
        log_prob -= R::dcauchy(current_state, 0.0, interaction_scale, true);

        // Accept/reject
        const double U = R::unif_rand();
        if (std::log(U) < log_prob) {
          double state_diff = proposed_state - current_state;
          interactions(variable1, variable2) = proposed_state;
          interactions(variable2, variable1) = proposed_state;

          // Update residual matrix
          for (arma::uword person = 0; person < no_persons; person++) {
            residual_matrix(person, variable1) += observations(person, variable2) * state_diff;
            residual_matrix(person, variable2) += observations(person, variable1) * state_diff;
          }
        }

        // Robbins-Monro update of proposal SD
        double acc_prob = (log_prob > 0.0) ? 1.0 : std::exp(log_prob);
        double updated_sd = proposal_sd(variable1, variable2) +
          (acc_prob - target_ar) * std::pow(static_cast<double>(t), -phi);

        if (std::isnan(updated_sd)) updated_sd = 1.0;

        proposal_sd(variable1, variable2) = std::clamp(updated_sd, epsilon_lo, epsilon_hi);
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
  double proposed_state;
  double current_state;
  double log_prob;
  double U;

  arma::uword variable1;
  arma::uword variable2;

  for(arma::uword cntr = 0; cntr < no_interactions; cntr ++) {
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
      for(arma::uword person = 0; person < no_persons; person++) {
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
    const arma::umat& observations,
    const arma::uvec& num_categories,
    const double interaction_scale,
    arma::mat& proposal_sd,
    arma::mat& proposal_sd_blumecapel,
    const arma::umat& index,
    const arma::umat& num_obs_categories,
    const arma::umat& sufficient_blume_capel,
    const double threshold_alpha,
    const double threshold_beta,
    const arma::uword no_persons,
    const arma::uword no_variables,
    const arma::uword no_interactions,
    const arma::uword no_thresholds,
    const arma::uword max_num_categories,
    arma::umat& indicator,
    arma::mat& interactions,
    arma::mat& thresholds,
    arma::mat& residual_matrix,
    const arma::mat& theta,
    const  double phi,
    const double target_ar,
    const double epsilon_lo,
    const double epsilon_hi,
    const arma::uvec& is_ordinal_variable,
    const arma::uvec& reference_category,
    const bool edge_selection,
    double& step_size,
    const arma::uword t,
    arma::vec& dual_averaging_state,
    const double log_step_size_target,
    const arma::uword burnin_iters,
    const bool mala = false
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
  if(mala) {
    adamala_thresholds(
      thresholds, step_size, residual_matrix, num_categories, num_obs_categories,
      sufficient_blume_capel, reference_category, is_ordinal_variable, t,
      burnin_iters, dual_averaging_state, log_step_size_target,
      threshold_alpha, threshold_beta
    );

  } else {
    for(arma::uword variable = 0; variable < no_variables; variable++) {
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
    arma::umat& observations,
    arma::umat& indicator,
    arma::mat& interactions,
    arma::mat& thresholds,
    const arma::uvec& num_categories,
    const double interaction_scale,
    arma::mat& proposal_sd,
    arma::mat& proposal_sd_blumecapel,
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
    const bool save = false,
    const bool display_progress = false,
    bool edge_selection = true,
    bool mala = false
) {
  arma::uword cntr;
  arma::uword no_variables = observations.n_cols;
  arma::uword no_persons = observations.n_rows;
  arma::uword no_interactions = Index.n_rows;
  arma::uword no_thresholds = sum(num_categories);
  arma::uword max_num_categories = max(num_categories);
  arma::uvec K_values;  // To store sampled K values

  arma::uvec v = arma::regspace<arma::uvec>(0, no_interactions - 1);
  arma::uvec order(no_interactions);
  arma::umat index(no_interactions, 3);

  //Parameters of adaptive proposals -------------------------------------------
  double initial_step_size = 0.01;
  double step_size = initial_step_size;
  double log_step_size_target = std::log(10.0 * initial_step_size);
  double phi = 0.75;
  double target_ar = 0.234;
  double epsilon_lo = 1.0 / static_cast<double>(no_persons);
  double epsilon_hi = 2.0;
  arma::vec dual_averaging_state = arma::zeros<arma::vec>(3);
  // Indexes:
  // [0] = log_step_size
  // [1] = log_step_size_avg
  // [2] = acceptance_error_avg

  //The resizing based on ``save'' could probably be prettier ------------------
  arma::uword nrow = no_variables;
  arma::uword ncol_edges = no_variables;
  arma::uword ncol_thresholds = max_num_categories;

  if(save == true) {
    nrow = iter;
    ncol_edges= no_interactions;
    ncol_thresholds = no_thresholds;
  }

  arma::mat out_pairwise_effects(nrow, ncol_edges);
  arma::mat out_main_effects(nrow, ncol_thresholds);

  if(edge_selection == false) {
    for(arma::uword variable1 = 0; variable1 < no_variables - 1; variable1++) {
      for(arma::uword variable2 = variable1 + 1; variable2 < no_variables; variable2++) {
        indicator(variable1, variable2) = 1;
        indicator(variable2, variable1) = 1;
      }
    }
    nrow = 1;
    ncol_edges = 1;
  }
  arma::mat out_indicator(nrow, ncol_edges);

  arma::mat residual_matrix(no_persons, no_variables);
  for(arma::uword variable1 = 0; variable1 < no_variables; variable1++) {
    for(arma::uword person = 0; person < no_persons; person++) {
      for(arma::uword variable2 = 0; variable2 < no_variables; variable2++) {
        residual_matrix(person, variable1) +=
          observations(person, variable2) * interactions(variable2, variable1);
      }
    }
  }

  //Variable declaration edge prior
  arma::uvec cluster_allocations(no_variables);
  arma::mat cluster_prob(1, 1);
  arma::vec log_Vn(1);

  // store the allocation indices for each iteration
  arma::mat out_allocations(iter, no_variables);

  if(edge_prior == "Stochastic-Block") { // Initial Configuration of the cluster allocations
    cluster_allocations[0] = 0;
    cluster_allocations[1] = 1;
    for(arma::uword i = 2; i < no_variables; i++) {
      double U = R::unif_rand();
      if(U > 0.5){
        cluster_allocations[i] = 1;
      } else {
        cluster_allocations[i] = 0;
      }
    }

    cluster_prob = block_probs_mfm_sbm (
      cluster_allocations, indicator, no_variables, beta_bernoulli_alpha,
      beta_bernoulli_beta);

    for(arma::uword i = 0; i < no_variables - 1; i++) {
      for(arma::uword j = i + 1; j < no_variables; j++) {
        theta(i, j) = cluster_prob(cluster_allocations[i], cluster_allocations[j]);
        theta(j, i) = cluster_prob(cluster_allocations[i], cluster_allocations[j]);
      }
    }

    log_Vn = compute_Vn_mfm_sbm (
      no_variables, dirichlet_alpha, no_variables + 10, lambda);
  }

  //The Gibbs sampler ----------------------------------------------------------
  //First, we do burn-in iterations---------------------------------------------

  //When edge_selection = true we do 2 * burnin iterations. The first burnin
  // iterations without selection to ensure good starting values, and proposal
  // calibration. The second burnin iterations with selection.

  arma::uword first_burnin = burnin;
  arma::uword second_burnin = 0;
  if(edge_selection == true)
    second_burnin = burnin;
  bool input_edge_selection = edge_selection;
  edge_selection = false;
  arma::uword burnin_iters = first_burnin + second_burnin;


  //Progress bar ---------------------------------------------------------------
  Progress p(iter + first_burnin + second_burnin, display_progress);

  for(arma::uword iteration = 0; iteration < first_burnin + second_burnin; iteration++) {
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

    for(arma::uword cntr = 0; cntr < no_interactions; cntr++) {
      index(cntr, 0) = Index(order[cntr], 0);
      index(cntr, 1) = Index(order[cntr], 1);
      index(cntr, 2) = Index(order[cntr], 2);
    }

    if(na_impute == true) {
      List out = impute_missing_data (
        interactions, thresholds, observations, num_obs_categories,
        sufficient_blume_capel, num_categories, residual_matrix,
        missing_index, is_ordinal_variable, reference_category);

      arma::umat observations = out["observations"];
      arma::umat num_obs_categories = out["num_obs_categories"];
      arma::umat sufficient_blume_capel = out["sufficient_blume_capel"];
      arma::mat residual_matrix = out["residual_matrix"];
    }

    List out = gibbs_step_gm (
      observations, num_categories, interaction_scale, proposal_sd,
      proposal_sd_blumecapel, index, num_obs_categories, sufficient_blume_capel,
      threshold_alpha, threshold_beta, no_persons,  no_variables,
      no_interactions, no_thresholds, max_num_categories, indicator,
      interactions, thresholds, residual_matrix, theta, phi, target_ar,
      epsilon_lo, epsilon_hi, is_ordinal_variable,
      reference_category, edge_selection, step_size, iteration + 1,
      dual_averaging_state, log_step_size_target, burnin_iters, mala);

    arma::umat indicator = out["indicator"];
    arma::mat interactions = out["interactions"];
    arma::mat thresholds = out["thresholds"];
    arma::mat residual_matrix = out["residual_matrix"];
    arma::mat proposal_sd = out["proposal_sd"];

    if(edge_selection == true) {
      if(edge_prior == "Beta-Bernoulli") {
        arma::uword sumG = 0;
        for(arma::uword i = 0; i < no_variables - 1; i++) {
          for(arma::uword j = i + 1; j < no_variables; j++) {
            sumG += indicator(i, j);
          }
        }
        double probability = R::rbeta(beta_bernoulli_alpha + sumG,
                                      beta_bernoulli_beta + no_interactions - sumG);

        for(arma::uword i = 0; i < no_variables - 1; i++) {
          for(arma::uword j = i + 1; j < no_variables; j++) {
            theta(i, j) = probability;
            theta(j, i) = probability;
          }
        }
      }
      if(edge_prior == "Stochastic-Block") {
        cluster_allocations = block_allocations_mfm_sbm (
          cluster_allocations, no_variables, log_Vn, cluster_prob, indicator,
          dirichlet_alpha, beta_bernoulli_alpha, beta_bernoulli_beta);

        cluster_prob = block_probs_mfm_sbm (
          cluster_allocations, indicator, no_variables, beta_bernoulli_alpha,
          beta_bernoulli_beta);

        for(arma::uword i = 0; i < no_variables - 1; i++) {
          for(arma::uword j = i + 1; j < no_variables; j++) {
            theta(i, j) = cluster_prob(cluster_allocations[i], cluster_allocations[j]);
            theta(j, i) = cluster_prob(cluster_allocations[i], cluster_allocations[j]);
          }
        }
      }
    }
  }
  //To ensure that edge_selection is reinstated to the input value -------------
  edge_selection = input_edge_selection;

  //The post burn-in iterations ------------------------------------------------
  for(arma::uword iteration = 0; iteration < iter; iteration++) {
    if (Progress::check_abort()) {
      if(edge_selection == true && edge_prior == "Stochastic-Block") {
        return List::create(Named("indicator") = out_indicator,
                            Named("pairwise_effects") = out_pairwise_effects,
                            Named("main_effects") = out_main_effects,
                            Named("allocations") = out_allocations);
      } if(edge_selection == true && edge_prior != "Stochastic-Block") {
        return List::create(Named("indicator") = out_indicator,
                            Named("pairwise_effects") = out_pairwise_effects,
                            Named("main_effects") = out_main_effects);
      }
      else {
        return List::create(Named("pairwise_effects") = out_pairwise_effects,
                            Named("main_effects") = out_main_effects);
      }
    }
    p.increment();

    //Update interactions and model (between model move) -----------------------
    // Create a random ordering of pairwise effects for updating
    arma::uvec order = arma::randperm(v.n_elem);

    for(arma::uword cntr = 0; cntr < no_interactions; cntr++) {
      index(cntr, 0) = Index(order[cntr], 0);
      index(cntr, 1) = Index(order[cntr], 1);
      index(cntr, 2) = Index(order[cntr], 2);
    }

    if(na_impute == true) {
      List out = impute_missing_data (
        interactions, thresholds, observations, num_obs_categories,
        sufficient_blume_capel, num_categories, residual_matrix,
        missing_index, is_ordinal_variable, reference_category);

      arma::umat observations = out["observations"];
      arma::umat num_obs_categories = out["num_obs_categories"];
      arma::umat sufficient_blume_capel = out["sufficient_blume_capel"];
      arma::mat residual_matrix = out["residual_matrix"];
    }

    List out = gibbs_step_gm (
      observations, num_categories, interaction_scale, proposal_sd,
      proposal_sd_blumecapel, index, num_obs_categories, sufficient_blume_capel,
      threshold_alpha, threshold_beta, no_persons,  no_variables,
      no_interactions, no_thresholds, max_num_categories, indicator,
      interactions, thresholds, residual_matrix, theta, phi, target_ar,
      epsilon_lo, epsilon_hi, is_ordinal_variable,
      reference_category, edge_selection, step_size, iteration + 1,
      dual_averaging_state, log_step_size_target, burnin_iters, mala);

    arma::umat indicator = out["indicator"];
    arma::mat interactions = out["interactions"];
    arma::mat thresholds = out["thresholds"];
    arma::mat residual_matrix = out["residual_matrix"];
    arma::mat proposal_sd = out["proposal_sd"];

    if(edge_selection == true) {
      if(edge_prior == "Beta-Bernoulli") {
        arma::uword sumG = 0;
        for(arma::uword i = 0; i < no_variables - 1; i++) {
          for(arma::uword j = i + 1; j < no_variables; j++) {
            sumG += indicator(i, j);
          }
        }
        double probability = R::rbeta(beta_bernoulli_alpha + sumG,
                                      beta_bernoulli_beta + no_interactions - sumG);

        for(arma::uword i = 0; i < no_variables - 1; i++) {
          for(arma::uword j = i + 1; j < no_variables; j++) {
            theta(i, j) = probability;
            theta(j, i) = probability;
          }
        }
      }
      if(edge_prior == "Stochastic-Block") {
        cluster_allocations = block_allocations_mfm_sbm (
          cluster_allocations, no_variables, log_Vn, cluster_prob, indicator,
          dirichlet_alpha, beta_bernoulli_alpha, beta_bernoulli_beta);

        cluster_prob = block_probs_mfm_sbm (
          cluster_allocations, indicator, no_variables, beta_bernoulli_alpha,
          beta_bernoulli_beta);

        for(arma::uword i = 0; i < no_variables - 1; i++) {
          for(arma::uword j = i + 1; j < no_variables; j++) {
            theta(i, j) = cluster_prob(cluster_allocations[i], cluster_allocations[j]);
            theta(j, i) = cluster_prob(cluster_allocations[i], cluster_allocations[j]);
          }
        }
      }
    }


    //Output -------------------------------------------------------------------
    if(save == true) {
      //Save raw samples -------------------------------------------------------
      cntr = 0;
      for(arma::uword variable1 = 0; variable1 < no_variables - 1; variable1++) {
        for(arma::uword variable2 = variable1 + 1; variable2 < no_variables; variable2++) {
          if(edge_selection == true) {
            out_indicator(iteration, cntr) = indicator(variable1, variable2);
          }
          out_pairwise_effects(iteration, cntr) = interactions(variable1, variable2);
          cntr++;
        }
      }
      cntr = 0;
      for(arma::uword variable = 0; variable < no_variables; variable++) {
        if(is_ordinal_variable[variable] == true) {
          for(arma::uword category = 0; category < num_categories[variable]; category++) {
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
      for(arma::uword variable1 = 0; variable1 < no_variables - 1; variable1++) {
        for(arma::uword variable2 = variable1 + 1; variable2 < no_variables; variable2++) {
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
          for(arma::uword category = 0; category < num_categories[variable1]; category++) {
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
        for(arma::uword category = 0; category < num_categories[no_variables - 1]; category++) {
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

    for(arma::uword i = 0; i < no_variables; i++) {
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
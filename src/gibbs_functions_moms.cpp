#include <RcppArmadillo.h>

#include <Rcpp.h>
#include "gibbs_functions.h"
using namespace Rcpp;


double hessian_log_pseudoposterior_interaction_single (
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
  double hessian = 0.0;

  // --- Contribution from var1
  int num_categories_var1 = num_categories (var1);
  arma::vec rest_scores_var1 = observations * pairwise_effects.col (var1);  // β_{var1,var1} = 0
  arma::vec numerator_var1_E (num_persons, arma::fill::zeros);
  arma::vec denominator_var1 (num_persons, arma::fill::zeros);
  arma::vec numerator_var1_E2 (num_persons, arma::fill::zeros);
  arma::vec bounds_var1 = arma::max (rest_scores_var1, arma::zeros<arma::vec> (num_persons)) * num_categories_var1;

  if (is_ordinal_variable (var1)) {
    denominator_var1 += arma::exp ( -bounds_var1 );
    for (int category = 0; category < num_categories_var1; category++) {
      arma::vec exponent = main_effects (var1, category) + (category + 1) * rest_scores_var1 - bounds_var1;
      arma::vec weight = arma::exp (exponent);
      denominator_var1 += weight;
      numerator_var1_E += (category + 1) * x_var2 % weight;
      numerator_var1_E2 += (category + 1) * (category + 1) * x_var2 % x_var2 % weight;
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
      numerator_var1_E += category * x_var2 % weight;
      numerator_var1_E2 += category * category * x_var2 % x_var2 % weight;
    }
  }
  //- E((XiXj)^2)
  hessian -= arma::accu (numerator_var1_E2 / denominator_var1);

  //+E(XiXj)^2
  arma::vec expectation = numerator_var1_E / denominator_var1;
  hessian += arma::accu(arma::square(expectation));

  // --- Contribution from var2
  int num_categories_var2 = num_categories (var2);
  arma::vec rest_scores_var2 = observations * pairwise_effects.col (var2);
  arma::vec numerator_var2_E (num_persons, arma::fill::zeros);
  arma::vec numerator_var2_E2 (num_persons, arma::fill::zeros);
  arma::vec denominator_var2 (num_persons, arma::fill::zeros);
  arma::vec bounds_var2 = arma::max (rest_scores_var2, arma::zeros<arma::vec> (num_persons)) * num_categories_var2;

  if (is_ordinal_variable (var2)) {
    denominator_var2 += arma::exp ( -bounds_var2 );
    for (int category = 0; category < num_categories_var2; category++) {
      arma::vec exponent = main_effects (var2, category) + (category + 1) * rest_scores_var2 - bounds_var2;
      arma::vec weight = arma::exp (exponent);
      denominator_var2 += weight;
      numerator_var2_E += (category + 1) * x_var1 % weight;
      numerator_var2_E2 += (category + 1) * (category + 1) * x_var1 % x_var1 % weight;
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
      numerator_var2_E += category * x_var1 % weight;
      numerator_var2_E2 += category * category * x_var1 % x_var1 % weight;
    }
  }

  //- E((XiXj)^2)
  hessian -= arma::accu (numerator_var2_E2 / denominator_var2);

  //+E(XiXj)^2
  expectation = numerator_var2_E / denominator_var2;
  hessian += arma::accu(arma::square(expectation));


  // --- Cauchy prior derivative
  double beta = pairwise_effects (var1, var2) * pairwise_effects (var1, var2);
  double s = interaction_scale * interaction_scale;
  hessian += 2.0 * (beta - s) / ((beta + s) * (beta + s));

  return hessian;
}



double log_pseudoposterior_interactions_single (
    const int var1,
    const int var2,
    const arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    const double interaction_scale
) {
  const int num_observations = observations.n_rows;

  // Convert to double matrix for trace calculation
  arma::mat real_observations = arma::conv_to<arma::mat>::from (observations);

  // Leading term: trace(X * B * X^T)
  double log_pseudo_posterior = arma::trace (real_observations * pairwise_effects * real_observations.t ());

  int var = var1;
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
  log_pseudo_posterior -= arma::accu (arma::log (denominator));
  log_pseudo_posterior -= arma::accu (bounds);

  var = var2;
  num_categories_var = num_categories (var);

  // Compute rest score: contribution from other variables
  rest_scores = observations * pairwise_effects.col (var);
  bounds = arma::max (rest_scores, arma::zeros<arma::vec> (num_observations)) * num_categories_var;
  denominator = arma::zeros (num_observations);

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
  log_pseudo_posterior -= arma::accu (arma::log (denominator));
  log_pseudo_posterior -= arma::accu (bounds);

  // Add Cauchy prior terms for included pairwise effects
  if (inclusion_indicator (var1, var2) == 1) {
    log_pseudo_posterior += R::dcauchy (pairwise_effects (var1, var2), 0.0, interaction_scale, true);
  }

  return log_pseudo_posterior;
}



double posterior_inclusion_probability (
    arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    const arma::imat& inclusion_indicator,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const int num_persons,
    const int variable1,
    const int variable2,
    const arma::mat& residual_matrix,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    const double interaction_scale,
    const arma::mat& inclusion_probability
) {
  // Step 1 Optimization of pseudoposterior w.r.t. the interaction
  double x = pairwise_effects(variable1, variable2);

  const int max_steps = 10;
  const double tolerance = std::sqrt (std::numeric_limits<double>::epsilon ());
  double hessian_at_x;
  arma::mat optimized_matrix = pairwise_effects;
  optimized_matrix(variable1, variable2) = x;
  optimized_matrix(variable2, variable1) = x;

  // Use Newton-Raphson to find the mode
  for (int t = 0; t < max_steps; t++) {
    double gradient_at_x = gradient_log_pseudoposterior_interaction_single (
      variable1, variable2, optimized_matrix, main_effects, observations,
      num_categories, is_ordinal_variable, reference_category, interaction_scale
    );

    hessian_at_x = hessian_log_pseudoposterior_interaction_single (
      variable1, variable2, optimized_matrix, main_effects, observations,
      num_categories, is_ordinal_variable, reference_category, interaction_scale
    );

    double x_new = x - gradient_at_x / hessian_at_x;

    if (std::abs(x_new - x) < tolerance) {
      x = x_new;
      break;
    }
    x = x_new;
    optimized_matrix(variable1, variable2) = x_new;
    optimized_matrix(variable2, variable1) = x_new;
  }

  if(std::abs(x - pairwise_effects(variable1, variable2)) > 5) {
    // Newton-Raphson probably diverged
    x = pairwise_effects(variable1, variable2);
  }

  // Step 2: Laplace Approximation of integral
  arma::imat one_indicator = inclusion_indicator;
  one_indicator(variable1, variable2) = 1;
  one_indicator(variable2, variable1) = 1;
  double log_post_end = log_pseudoposterior_interactions_single (
    variable1, variable2, optimized_matrix, main_effects, observations,
    num_categories, one_indicator, is_ordinal_variable,
    reference_category, interaction_scale
  );
  const double I = log_post_end + 0.5 * log(2.0 * M_PI) - 0.5 * log(-hessian_at_x);


  // Step 3: Posterior inclusion probability
  arma::mat null_matrix = pairwise_effects;
  null_matrix(variable1, variable2) = 0.0;
  null_matrix(variable2, variable1) = 0.0;
  arma::imat null_indicator = inclusion_indicator;
  null_indicator(variable1, variable2) = 0;
  null_indicator(variable2, variable1) = 0;

  const double I0 = log_pseudoposterior_interactions_single (
    variable1, variable2, null_matrix, main_effects, observations,
    num_categories, null_indicator, is_ordinal_variable,
    reference_category, interaction_scale
  );
  const double theta = inclusion_probability(variable1, variable2);
  double bound = (I > I0) ? I : I0;
  bound = (bound > 0) ? bound : 0.0;
  double eI = std::exp(I - bound);
  double eI0 = std::exp(I0 - bound);

  double posterior_inclusion_probability = eI * theta / (eI * theta + eI0 * (1 - theta));

  return posterior_inclusion_probability;
}



// [[Rcpp::export]]
void update_indicator_interaction_pair_with_metropolis_and_marginal_moms (
    arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    arma::imat& indicator,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::mat& proposal_sd,
    const double interaction_scale,
    const arma::imat& index,
    const int num_interactions,
    const int num_persons,
    arma::mat& residual_matrix,
    const arma::mat& inclusion_probability,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category
) {
  for (int cntr = 0; cntr < num_interactions; cntr++) {
    const int variable1 = index(cntr, 1);
    const int variable2 = index(cntr, 2);

    double p = posterior_inclusion_probability (
      pairwise_effects, main_effects, indicator, observations,
      num_categories, num_persons, variable1, variable2, residual_matrix,
      is_ordinal_variable,reference_category, interaction_scale,
      inclusion_probability
    );

    const double current_state = pairwise_effects(variable1, variable2);

    if ( R::unif_rand() < p) {
      indicator(variable1, variable2) = 1;
      indicator(variable2, variable1) = 1;


      const double proposed_state = R::rnorm(current_state, proposal_sd(variable1, variable2));

      // Compute log acceptance ratio: data + prior
      double log_acceptance = log_pseudolikelihood_ratio_interaction (
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

    } else {
      indicator(variable1, variable2) = 0;
      indicator(variable2, variable1) = 0;

      pairwise_effects(variable1, variable2) = 0.0;
      pairwise_effects(variable2, variable1) = 0.0;

      const double delta = 0.0 - current_state;

      // Vectorized residual update
      residual_matrix.col(variable1) += arma::conv_to<arma::vec>::from(observations.col(variable2)) * delta;
      residual_matrix.col(variable2) += arma::conv_to<arma::vec>::from(observations.col(variable1)) * delta;
    }
  }
}
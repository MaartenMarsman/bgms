// Test interface for the bgmCompare log-pseudoposterior and gradient.
//
// Mirrors ggm_test_logp_and_gradient() and mixed_test_logp_and_gradient():
// exposes the flat-parameter-space (value, gradient) pair to R so it can be
// checked against an independent reference implementation. Not used by the
// sampler; the sampler calls logp_and_gradient() directly.
//
// The flat layout is the sampler's own: overall mains, overall pairs, active
// main differences, active pair differences -- see
// vectorize_model_parameters_bgmcompare().

#include <RcppArmadillo.h>
#include "models/bgmCompare/bgmCompare_helper.h"
#include "models/bgmCompare/bgmCompare_logp_and_grad.h"
#include "priors/parameter_prior.h"
#include "utils/common_helpers.h"

// [[Rcpp::export]]
Rcpp::List bgmCompare_test_logp_and_gradient(
    const arma::vec& params,
    const arma::imat& observations,
    const arma::imat& group_indices,
    int num_groups,
    const std::vector<arma::imat>& counts_per_category,
    const std::vector<arma::imat>& blume_capel_stats,
    const std::vector<arma::mat>&  pairwise_stats,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::imat& inclusion_indicator,
    const arma::mat& projection,
    double pairwise_scale = 1.0,
    double difference_scale = 1.0,
    double main_alpha = 1.0,
    double main_beta = 1.0,
    std::string interaction_prior_type = "cauchy",
    std::string difference_prior_type = "cauchy",
    std::string threshold_prior_type = "beta-prime",
    double threshold_scale = 1.0)
{
  const int num_variables = observations.n_cols;
  const int num_main = count_num_main_effects(
    num_categories, is_ordinal_variable
  );
  const int num_pair = num_variables * (num_variables - 1) / 2;

  auto interaction_prior = create_parameter_prior(
    interaction_prior_type, pairwise_scale
  );
  auto difference_prior = create_parameter_prior(
    difference_prior_type, difference_scale
  );
  auto threshold_prior = create_parameter_prior(
    threshold_prior_type, threshold_scale, main_alpha, main_beta
  );

  arma::mat main_effects(num_main, num_groups, arma::fill::zeros);
  arma::mat pairwise_effects(num_pair, num_groups, arma::fill::zeros);

  auto index_maps = build_index_maps(
    main_effects, pairwise_effects, inclusion_indicator,
    main_effect_indices, pairwise_effect_indices,
    num_categories, is_ordinal_variable
  );
  const arma::imat& main_index = index_maps.first;
  const arma::imat& pair_index = index_maps.second;

  const arma::uword expected = total_length(
    num_variables, main_effect_indices, pairwise_effect_indices,
    inclusion_indicator, num_categories, is_ordinal_variable, num_groups
  );
  if (params.n_elem != expected) {
    Rcpp::stop(
      "params has length %d but the active parameter vector has length %d.",
      static_cast<int>(params.n_elem), static_cast<int>(expected)
    );
  }

  unvectorize_model_parameters_bgmcompare(
    params, main_effects, pairwise_effects, inclusion_indicator,
    main_effect_indices, pairwise_effect_indices, num_groups,
    num_categories, is_ordinal_variable
  );

  const arma::vec grad_obs = gradient_observed_active(
    main_effect_indices, pairwise_effect_indices, projection,
    observations, group_indices, num_categories, inclusion_indicator,
    counts_per_category, blume_capel_stats, pairwise_stats,
    num_groups, is_ordinal_variable, baseline_category,
    main_index, pair_index
  );

  const arma::mat observations_double = arma::conv_to<arma::mat>::from(
    observations
  );

  auto result = logp_and_gradient(
    main_effects, pairwise_effects, main_effect_indices,
    pairwise_effect_indices, projection, observations_double, group_indices,
    num_categories, counts_per_category, blume_capel_stats, pairwise_stats,
    num_groups, inclusion_indicator, is_ordinal_variable, baseline_category,
    main_index, pair_index, grad_obs,
    *interaction_prior, *difference_prior, *threshold_prior
  );

  return Rcpp::List::create(
    Rcpp::Named("value") = result.first,
    Rcpp::Named("gradient") = Rcpp::wrap(result.second)
  );
}

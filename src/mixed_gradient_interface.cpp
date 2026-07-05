// Test interface for the mixed MRF gradient engine.
//
// Exposes the theta-space logp_and_gradient to R for validation.

#include <RcppArmadillo.h>
#include "models/mixed/mixed_mrf_model.h"
#include "priors/parameter_prior.h"

// [[Rcpp::export]]
Rcpp::List mixed_test_logp_and_gradient(
    const arma::vec& params,
    const arma::imat& discrete_observations,
    const arma::mat& continuous_observations,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const arma::imat& edge_indicators,
    double pairwise_scale,
    double main_alpha = 1.0,
    double main_beta = 1.0,
    std::string interaction_prior_type = "cauchy",
    std::string threshold_prior_type = "beta-prime",
    double threshold_scale = 1.0,
    std::string means_prior_type = "normal",
    double means_scale = 1.0,
    std::string diagonal_prior_type = "gamma",
    double diagonal_shape = 1.0,
    double diagonal_rate = 1.0)
{
    size_t p = discrete_observations.n_cols;
    size_t q = continuous_observations.n_cols;
    size_t total = p + q;

    arma::mat inc_prob(total, total, arma::fill::value(0.5));
    bool edge_selection = false;

    MixedMRFModel model(
        discrete_observations, continuous_observations,
        num_categories, is_ordinal_variable, baseline_category,
        inc_prob, edge_indicators, edge_selection,
        create_parameter_prior(interaction_prior_type, pairwise_scale),
        create_parameter_prior(threshold_prior_type, threshold_scale, main_alpha, main_beta),
        create_parameter_prior(means_prior_type, means_scale),
        create_scale_prior(diagonal_prior_type, diagonal_shape, diagonal_rate),
        42
    );

    auto result = model.logp_and_gradient(params);

    return Rcpp::List::create(
        Rcpp::Named("value") = result.first,
        Rcpp::Named("gradient") = Rcpp::wrap(result.second)
    );
}

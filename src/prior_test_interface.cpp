// Test interface for the polymorphic parameter prior classes.
//
// Exposes logp and grad to R for unit testing.

#include <RcppArmadillo.h>
#include "priors/parameter_prior.h"
#include "models/ggm/graph_constraint_structure.h"
#include "models/ggm/ggm_gradient.h"


// [[Rcpp::export]]
Rcpp::List test_parameter_prior(
    const std::string& type,
    double x,
    double scale = 1.0,
    double alpha = 0.5,
    double beta = 0.5
) {
    auto prior = create_parameter_prior(type, scale, alpha, beta);

    return Rcpp::List::create(
        Rcpp::Named("logp") = prior->logp(x),
        Rcpp::Named("grad") = prior->grad(x)
    );
}


// [[Rcpp::export]]
Rcpp::List test_scale_prior(
    const std::string& type,
    double x,
    double shape = 1.0,
    double rate = 1.0
) {
    auto prior = create_scale_prior(type, shape, rate);

    return Rcpp::List::create(
        Rcpp::Named("logp") = prior->logp(x),
        Rcpp::Named("grad") = prior->grad(x)
    );
}


// [[Rcpp::export]]
Rcpp::List ggm_test_logp_and_gradient_prior(
    const arma::vec& theta,
    const arma::mat& suf_stat,
    int n,
    const arma::imat& edge_indicators,
    const std::string& interaction_prior_type = "cauchy",
    double interaction_scale = 1.0,
    double interaction_alpha = 0.5,
    double interaction_beta = 0.5,
    const std::string& diagonal_prior_type = "gamma",
    double diagonal_shape = 1.0,
    double diagonal_rate = 1.0
) {
    GraphConstraintStructure cs;
    cs.build(edge_indicators);

    auto ip = create_parameter_prior(
        interaction_prior_type, interaction_scale,
        interaction_alpha, interaction_beta);
    auto dp = create_scale_prior(diagonal_prior_type, diagonal_shape, diagonal_rate);

    GGMGradientEngine engine;
    engine.rebuild(cs, static_cast<size_t>(n), suf_stat, *ip, *dp);

    auto result = engine.logp_and_gradient(theta);

    return Rcpp::List::create(
        Rcpp::Named("value") = result.first,
        Rcpp::Named("gradient") = Rcpp::wrap(result.second)
    );
}


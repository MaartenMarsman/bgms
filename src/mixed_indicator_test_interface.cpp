// -----------------------------------------------------------------------------
// mixed_indicator_test_interface.cpp
//
// R-facing test entry for the mixed-MRF edge-indicator bookkeeping. Builds a
// small model, runs its edge-indicator sweeps, and returns the internal
// (p+q) x (p+q) indicator matrix. The stochastic block edge prior reads full
// columns of this matrix, so tests use this entry to check that both
// triangles stay in sync under the discrete, continuous, and cross moves.
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>

#include "models/mixed/mixed_mrf_model.h"
#include "priors/parameter_prior.h"

// -----------------------------------------------------------------------------
// test_mixed_edge_indicator_matrix:
//   Construct a MixedMRFModel with all-ones initial indicators (the
//   run_sampler default), run n_sweeps edge-indicator sweeps with edge
//   selection active, and return the internal edge-indicator matrix.
// -----------------------------------------------------------------------------

// [[Rcpp::export(name = "test_mixed_edge_indicator_matrix")]]
arma::imat test_mixed_edge_indicator_matrix(
    const arma::imat& discrete_observations,
    const arma::mat& continuous_observations,
    const arma::ivec& num_categories,
    int n_sweeps,
    int seed
) {
    const int p = discrete_observations.n_cols;
    const int q = continuous_observations.n_cols;
    const int d = p + q;

    arma::uvec is_ordinal(p, arma::fill::ones);
    arma::ivec baseline(p, arma::fill::zeros);
    arma::mat inclusion_probability(d, d);
    inclusion_probability.fill(0.5);
    arma::imat initial_indicators(d, d, arma::fill::ones);

    MixedMRFModel model(
        discrete_observations, continuous_observations,
        num_categories, is_ordinal, baseline,
        inclusion_probability, initial_indicators,
        /*edge_selection=*/true,
        create_parameter_prior("cauchy", 2.5, NA_REAL, NA_REAL),
        create_parameter_prior("beta-prime", 1.0, 0.5, 0.5),
        create_parameter_prior("normal", 1.0, NA_REAL, NA_REAL),
        create_scale_prior("gamma", 1.0, 1.0),
        seed
    );
    model.set_edge_selection_active(true);

    for (int iter = 0; iter < n_sweeps; ++iter) {
        model.prepare_iteration();
        model.update_edge_indicators();
    }

    return model.get_edge_indicators();
}

// -----------------------------------------------------------------------------
// ggm_gibbs_test_interface.cpp
//
// R-facing test entries for the GGM row-block Gibbs sampler. Kept separate from
// the sampler dispatch (chain_runner) and from ggm_gradient_interface.cpp (which
// exposes the NUTS gradient): this file only drives the conjugate within-step in
// isolation so its full-conditional math can be validated without the
// between-step, edge moves, or warmup machinery.
// -----------------------------------------------------------------------------

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "models/ggm/ggm_model.h"
#include "priors/parameter_prior.h"

// -----------------------------------------------------------------------------
// ggm_test_gibbs_sweep:
//   Construct a GGM with a Normal slab and Gamma(alpha=1) diagonal prior, then
//   run n_sweeps full row-block Gibbs sweeps via GGMModel::do_one_gibbs_step()
//   -- the same sweep the GibbsSampler drives in a real run. The graph is fixed
//   (no edge selection). Returns the p x p precision matrix K after every sweep.
// -----------------------------------------------------------------------------

// [[Rcpp::export(name = "ggm_test_gibbs_sweep")]]
Rcpp::List ggm_test_gibbs_sweep(
    const arma::mat& suf_stat,
    int n,
    const arma::imat& edge_indicators,
    double pairwise_scale,
    double gamma_shape,
    double gamma_rate,
    int n_sweeps,
    int seed)
{
    const int p = edge_indicators.n_rows;
    auto ip = create_parameter_prior("normal", pairwise_scale);
    auto dp = create_scale_prior("gamma", gamma_shape, gamma_rate);

    arma::mat inc_prob(p, p, arma::fill::value(0.5));
    GGMModel model(n, suf_stat, inc_prob, edge_indicators,
                   /*edge_selection=*/false, std::move(ip), std::move(dp));
    model.set_determinant_tilt(0.0);
    model.set_seed(seed);

    if (!model.row_block_gibbs_eligible()) {
        Rcpp::stop("ggm_test_gibbs_sweep: model is not eligible "
                   "(requires Normal slab, Gamma alpha=1, delta=0).");
    }

    arma::cube K_samples(p, p, n_sweeps);
    for (int s = 0; s < n_sweeps; ++s) {
        model.do_one_gibbs_step(s);
        K_samples.slice(s) = model.get_precision_matrix();
    }

    return Rcpp::List::create(
        Rcpp::Named("K_samples") = K_samples
    );
}

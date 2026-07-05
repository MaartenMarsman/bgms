// Test interface for the GGM gradient engine and RATTLE projection.
//
// Exposes logp_and_gradient and forward_map to R for validation.
// Also exposes sample_ggm_prior() for sampling from the GGM prior
// from the GGM prior using NUTS.

#include <RcppArmadillo.h>
#include "models/ggm/graph_constraint_structure.h"
#include "models/ggm/ggm_gradient.h"
#include "models/ggm/ggm_model.h"
#include "mcmc/algorithms/leapfrog.h"
#include "mcmc/algorithms/nuts.h"
#include "rng/rng_utils.h"
#include "priors/parameter_prior.h"
#include "priors/edge_prior.h"
#include "utils/common_helpers.h"
#include "mcmc/execution/sampler_config.h"
#include "mcmc/execution/chain_runner.h"
#include "mcmc/execution/chain_result.h"
#include "utils/progress_manager.h"

// [[Rcpp::export]]
Rcpp::List ggm_test_logp_and_gradient(
    const arma::vec& theta,
    const arma::mat& suf_stat,
    int n,
    const arma::imat& edge_indicators,
    double pairwise_scale)
{
    GraphConstraintStructure cs;
    cs.build(edge_indicators);

    CauchyPrior ip(pairwise_scale);
    GammaScalePrior dp(1.0, 1.0);
    GGMGradientEngine engine;
    engine.rebuild(cs, static_cast<size_t>(n), suf_stat, ip, dp);

    auto result = engine.logp_and_gradient(theta);

    return Rcpp::List::create(
        Rcpp::Named("value") = result.first,
        Rcpp::Named("gradient") = Rcpp::wrap(result.second)
    );
}

// [[Rcpp::export]]
Rcpp::List ggm_test_forward_map(
    const arma::vec& theta,
    const arma::imat& edge_indicators)
{
    GraphConstraintStructure cs;
    cs.build(edge_indicators);

    // Minimal engine just for forward map (no suf_stat needed)
    arma::mat dummy_S(edge_indicators.n_rows, edge_indicators.n_rows, arma::fill::zeros);
    CauchyPrior ip(1.0);
    GammaScalePrior dp(1.0, 1.0);
    GGMGradientEngine engine;
    engine.rebuild(cs, 100, dummy_S, ip, dp);

    ForwardMapResult fm = engine.forward_map(theta);

    return Rcpp::List::create(
        Rcpp::Named("Phi") = Rcpp::wrap(fm.Phi),
        Rcpp::Named("K") = Rcpp::wrap(fm.K),
        Rcpp::Named("log_det_jacobian") = fm.log_det_jacobian,
        Rcpp::Named("psi") = Rcpp::wrap(fm.psi)
    );
}



// -----------------------------------------------------------------------------
// sample_ggm_prior: Sample from the GGM prior (K_yy = -K/2 partial-association
//                   scale; output reported as K) using NUTS.
// -----------------------------------------------------------------------------
// Uses the Cholesky parameterization with NUTS. By setting n=0 and S=0,
// the likelihood vanishes and the sampler targets the prior:
//   -K_ij/2 | graph ~ Cauchy(0, scale) or Normal(0, scale)       (included edges)
//   K_ij = 0                                                      (excluded edges)
//   K_ii/2 ~ Gamma(gamma_shape, gamma_rate)                       (diagonal)
//
// Note on convention: both `interaction_prior` and `diagonal_prior` are
// applied on the partial-association scale K_yy = -K/2 (consistent with the
// mixed-MRF continuous block), NOT on the precision entries directly.
// A `Normal(0, scale)` prior therefore constrains K_yy off-diagonals with
// standard deviation `scale`, equivalent to Normal(0, 2*scale) on K_{ij}.
// A `Gamma(shape, rate)` diagonal prior is applied to -K_yy_{ii} = K_{ii}/2,
// so the implied prior on K_{ii} is Gamma(shape, rate/2).
//
// edge_indicators: p x p integer matrix with 1 = edge included, 0 = excluded.
//   Defaults to all-ones (full graph). For edge selection SBC, pass the
//   graph drawn from the edge prior so K is sampled conditional on that graph.
//
// The prior is the product of the element-wise priors plus the
// Jacobian from theta -> K induced by the Cholesky parameterization.
//
// Uses the same NUTS infrastructure (windowed warmup, dual-averaging step-size
// adaptation, Welford diagonal mass-matrix adaptation) as sample_ggm and bgm.
// -----------------------------------------------------------------------------

// [[Rcpp::export(name = "sample_ggm_prior_cpp")]]
Rcpp::List sample_ggm_prior(
    int p,
    int n_samples,
    int n_warmup = 1000,
    double pairwise_scale = 2.5,
    const std::string& interaction_prior_type = "cauchy",
    const std::string& scale_prior_type = "gamma",
    double gamma_shape = 1.0,
    double gamma_rate = 1.0,
    double step_size = 0.1,
    int max_depth = 10,
    int seed = 1,
    bool verbose = true,
    Rcpp::Nullable<Rcpp::IntegerMatrix> edge_indicators_nullable = R_NilValue,
    double delta = 0.0)
{
    // Build edge indicators (default: full graph, no constraints)
    arma::imat edge_indicators;
    if(edge_indicators_nullable.isNotNull()) {
        edge_indicators = Rcpp::as<arma::imat>(
            Rcpp::IntegerMatrix(edge_indicators_nullable.get()));
    } else {
        edge_indicators.ones(p, p);
    }

    auto ip = create_parameter_prior(interaction_prior_type, pairwise_scale);
    auto dp = create_scale_prior(scale_prior_type, gamma_shape, gamma_rate);

    // Create model with n=0, S=0 so likelihood is flat (prior-only).
    // edge_selection=false ensures the graph stays fixed throughout:
    // no update_edge_indicators() calls, required for RATTLE correctness
    // when the graph is sparse.
    arma::mat suf_stat(p, p, arma::fill::zeros);
    arma::mat inc_prob(p, p, arma::fill::value(0.5));
    GGMModel model(0, suf_stat, inc_prob, edge_indicators,
                   false, std::move(ip), std::move(dp));
    model.set_determinant_tilt(delta);

    // Configure NUTS with the standard windowed warmup (dual averaging +
    // diagonal Welford mass-matrix adaptation). This is the same path used
    // by sample_ggm / bgm, so prior draws match the sampler used for the
    // full posterior.
    SamplerConfig config;
    config.sampler_type = "nuts";
    config.no_iter = n_samples;
    config.no_warmup = n_warmup;
    config.edge_selection = false;
    config.seed = seed;
    config.target_acceptance = 0.8;
    config.max_tree_depth = max_depth;
    config.initial_step_size = step_size;
    config.na_impute = false;

    // Edge prior is required by the API but is never consulted because
    // edge_selection=false. Use a cheap Bernoulli placeholder.
    auto edge_prior_obj = create_edge_prior(
        Bernoulli, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);

    // Progress manager: progress_type 2 = full bar (when verbose), 0 = silent.
    int progress_type = verbose ? 2 : 0;
    ProgressManager pm(1, n_samples, n_warmup, 10, progress_type,
                       verbose, R_NilValue);

    // Run a single chain.
    std::vector<ChainResult> results = run_mcmc_sampler(
        model, *edge_prior_obj, config, /*no_chains=*/1, /*no_threads=*/1, pm);
    pm.finish();

    if(results.empty() || results[0].error) {
        std::string msg = results.empty() ? "empty result" : results[0].error_msg;
        Rcpp::stop("sample_ggm_prior: chain failed (%s)", msg.c_str());
    }

    // Samples are stored as the upper triangle of K (p(p+1)/2 elements per
    // iteration, column-wise iteration index): see GGMModel::extract_upper_triangle.
    // Layout: for i in 0..p-1, for j in i..p-1, entry K(i,j).
    const arma::mat& samples = results[0].samples;  // dim x no_iter

    int n_edges = p * (p - 1) / 2;
    arma::mat K_offdiag_samples(n_samples, n_edges);
    arma::mat K_diag_samples(n_samples, p);

    for(int s = 0; s < n_samples; ++s) {
        const arma::vec& upper = samples.col(s);
        size_t e = 0;
        int off_idx = 0;
        for(int i = 0; i < p; ++i) {
            for(int j = i; j < p; ++j) {
                if(i == j) {
                    K_diag_samples(s, i) = upper(e);
                } else {
                    K_offdiag_samples(s, off_idx++) = upper(e);
                }
                ++e;
            }
        }
    }

    // Build column names
    Rcpp::CharacterVector offdiag_names(n_edges);
    int idx = 0;
    for(int i = 0; i < p - 1; ++i) {
        for(int j = i + 1; j < p; ++j) {
            offdiag_names[idx++] = "K_" + std::to_string(i + 1) + "_" +
                                   std::to_string(j + 1);
        }
    }

    Rcpp::CharacterVector diag_names(p);
    for(int i = 0; i < p; ++i) {
        diag_names[i] = "K_" + std::to_string(i + 1) + "_" +
                        std::to_string(i + 1);
    }

    return Rcpp::List::create(
        Rcpp::Named("K_offdiag") = K_offdiag_samples,
        Rcpp::Named("K_diag") = K_diag_samples,
        Rcpp::Named("offdiag_names") = offdiag_names,
        Rcpp::Named("diag_names") = diag_names,
        Rcpp::Named("step_size") = config.initial_step_size,
        Rcpp::Named("edge_indicators") = Rcpp::wrap(edge_indicators)
    );
}

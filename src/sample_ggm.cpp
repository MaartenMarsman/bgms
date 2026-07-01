#include <vector>
#include <memory>
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <tbb/global_control.h>

#include "models/ggm/ggm_model.h"
#include "utils/progress_manager.h"
#include "utils/common_helpers.h"
#include "priors/edge_prior.h"
#include "priors/parameter_prior.h"
#include "mcmc/execution/chain_result.h"
#include "mcmc/execution/chain_runner.h"
#include "mcmc/execution/sampler_config.h"

// [[Rcpp::export]]
Rcpp::List sample_ggm(
    const Rcpp::List& inputFromR,
    const arma::mat& prior_inclusion_prob,
    const arma::imat& initial_edge_indicators,
    const int no_iter,
    const int no_warmup,
    const int no_chains,
    const bool edge_selection,
    const std::string& sampler_type,
    const int seed,
    const int no_threads,
    const int progress_type,
    SEXP progress_callback = R_NilValue,
    const std::string& edge_prior = "Bernoulli",
    const double beta_bernoulli_alpha = 1.0,
    const double beta_bernoulli_beta = 1.0,
    const double beta_bernoulli_alpha_between = 1.0,
    const double beta_bernoulli_beta_between = 1.0,
    const double dirichlet_alpha = 1.0,
    const double lambda = 1.0,
    const double target_acceptance = 0.8,
    const int max_tree_depth = 10,
    const bool na_impute = false,
    const Rcpp::Nullable<Rcpp::IntegerMatrix> missing_index_nullable = R_NilValue,
    const double delta = 0.0
) {

    // Create parameter priors from R input
    double pairwise_scale = Rcpp::as<double>(inputFromR["pairwise_scale"]);
    std::string ipt_str = inputFromR.containsElementNamed("interaction_prior_type")
        ? Rcpp::as<std::string>(inputFromR["interaction_prior_type"]) : "cauchy";
    double ia = inputFromR.containsElementNamed("interaction_alpha")
        ? Rcpp::as<double>(inputFromR["interaction_alpha"]) : NA_REAL;
    double ib = inputFromR.containsElementNamed("interaction_beta")
        ? Rcpp::as<double>(inputFromR["interaction_beta"]) : NA_REAL;

    auto interaction_prior = create_parameter_prior(
        ipt_str, pairwise_scale, ia, ib);

    // Scale prior on precision diagonal
    std::string spt_str = inputFromR.containsElementNamed("scale_prior_type")
        ? Rcpp::as<std::string>(inputFromR["scale_prior_type"]) : "gamma";
    double s_shape = inputFromR.containsElementNamed("scale_shape")
        ? Rcpp::as<double>(inputFromR["scale_shape"]) : 1.0;
    double s_rate = inputFromR.containsElementNamed("scale_rate")
        ? Rcpp::as<double>(inputFromR["scale_rate"]) : 1.0;

    auto diagonal_prior = create_scale_prior(spt_str, s_shape, s_rate);

    // Create model from R input
    GGMModel model = createGGMModelFromR(
        inputFromR, prior_inclusion_prob, initial_edge_indicators,
        edge_selection, std::move(interaction_prior),
        std::move(diagonal_prior), na_impute);

    // Forward target_accept to the model's MH proposal-SD tuner.
    //   - Under "adaptive-metropolis": user's target_accept goes through
    //     directly (default 0.44 = componentwise RW MH optimum).
    //   - Under "nuts": user's target_accept (default 0.80) is the
    //     HMC step-size dual-averaging target and should NOT govern the
    //     between-model MH proposal SDs, which are still 1-D componentwise
    //     RW MH. Hardcode 0.44 there to keep stage-3b RM on the right
    //     fixed point.
    //   - Under "gibbs": the within-step draw is exact (no MH tuning); the
    //     target is inert, so use the same 0.44 fixed point.
    const double mh_target =
        (sampler_type == "adaptive-metropolis") ? target_acceptance : 0.44;
    model.set_metropolis_target_accept(mh_target);

    // Determinant-tilt prior on |K|: shifts both NUTS and MH targets by
    // delta * log|K|. delta = 0 is the default (untilted). Consumed by
    // both gradient paths and all four MH ratios in GGMModel.
    model.set_determinant_tilt(delta);

    // The row-block Gibbs sampler covers a Normal or Cauchy slab on the
    // off-diagonals and a Gamma prior on the precision diagonal. Fail fast
    // with a clear message rather than let update_row_block_gibbs cast a
    // mismatched prior.
    if (sampler_type == "gibbs" && !model.row_block_gibbs_eligible()) {
        Rcpp::stop(
            "update_method = \"gibbs\" needs a Normal or Cauchy interaction "
            "(slab) prior and a Gamma scale prior on the precision diagonal. "
            "The current priors do not meet this; use another update method "
            "or adjust the priors.");
    }

    // Edge selection with the Gibbs sampler keeps the Roverato between-step,
    // which reads the slab density directly. A Cauchy slab there would be
    // inconsistent with the scale-mixture the within-step uses, so that
    // combination is not yet supported.
    if (sampler_type == "gibbs" && edge_selection && ipt_str == "cauchy") {
        Rcpp::stop(
            "update_method = \"gibbs\" with edge selection needs a Normal slab; "
            "a Cauchy slab with edge selection is not yet supported. Use a "
            "Normal slab or set edge_selection = FALSE.");
    }
    if (sampler_type == "gibbs" && edge_selection && std::abs(s_shape - 1.0) > 1e-12) {
        Rcpp::stop(
            "update_method = \"gibbs\" with edge selection needs the Gamma "
            "scale prior shape = 1; a shape != 1 with edge selection is not yet "
            "supported. Set the scale prior shape to 1 or edge_selection = FALSE.");
    }

    // Set up missing data imputation (same pattern as OMRF)
    if (na_impute && missing_index_nullable.isNotNull()) {
        arma::imat missing_index = Rcpp::as<arma::imat>(
            Rcpp::IntegerMatrix(missing_index_nullable.get()));
        model.set_missing_data(missing_index);
    }

    // Configure sampler
    SamplerConfig config;
    config.sampler_type = sampler_type;
    config.no_iter = no_iter;
    config.no_warmup = no_warmup;
    config.edge_selection = edge_selection;
    config.seed = seed;
    config.target_acceptance = target_acceptance;
    config.max_tree_depth = max_tree_depth;
    config.na_impute = na_impute;

    // Set up progress manager
    ProgressManager pm(no_chains, no_iter, no_warmup, 50, progress_type, true, progress_callback);

    // Create edge prior
    EdgePrior edge_prior_enum = edge_prior_from_string(edge_prior);
    auto edge_prior_obj = create_edge_prior(
        edge_prior_enum,
        beta_bernoulli_alpha, beta_bernoulli_beta,
        beta_bernoulli_alpha_between, beta_bernoulli_beta_between,
        dirichlet_alpha, lambda
    );

    // Run MCMC using unified infrastructure
    std::vector<ChainResult> results = run_mcmc_sampler(
        model, *edge_prior_obj, config, no_chains, no_threads, pm);

    // Convert to R list format
    Rcpp::List output = convert_results_to_list(results);

    pm.finish();

    return output;
}
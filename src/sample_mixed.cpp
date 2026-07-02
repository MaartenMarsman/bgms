// sample_mixed.cpp - R interface for Mixed MRF model sampling
//
// Uses the unified MCMC runner infrastructure to sample from models with
// both discrete (ordinal / Blume-Capel) and continuous variables.
// Supports MH and NUTS samplers, with optional edge selection.
#include <vector>
#include <memory>
#include <RcppArmadillo.h>

#include "models/mixed/mixed_mrf_model.h"
#include "priors/parameter_prior.h"
#include "utils/progress_manager.h"
#include "utils/common_helpers.h"
#include "priors/edge_prior.h"
#include "mcmc/execution/chain_result.h"
#include "mcmc/execution/chain_runner.h"
#include "mcmc/execution/sampler_config.h"

// R-exported function to sample from a Mixed MRF model.
//
// @param inputFromR              List with model specification:
//                                  discrete_observations (integer matrix n x p),
//                                  continuous_observations (numeric matrix n x q),
//                                  num_categories (integer vector, length p),
//                                  is_ordinal_variable (integer vector, length p),
//                                  baseline_category (integer vector, length p),
//                                  main_alpha, main_beta, pairwise_scale (doubles)
// @param prior_inclusion_prob    Prior inclusion probabilities ((p+q) x (p+q) matrix)
// @param initial_edge_indicators Initial edge indicators ((p+q) x (p+q) integer matrix)
// @param no_iter                 Number of post-warmup iterations
// @param no_warmup               Number of warmup iterations
// @param no_chains               Number of parallel chains
// @param edge_selection          Whether to do edge selection (spike-and-slab)
// @param seed                    Random seed
// @param no_threads              Number of threads for parallel execution
// @param progress_type           Progress bar type
// @param progress_callback       R function (SEXP) called as callback(completed, total) at regular intervals, or R_NilValue
// @param edge_prior              Edge prior type
// @param beta_bernoulli_alpha         Beta-Bernoulli alpha hyperparameter
// @param beta_bernoulli_beta          Beta-Bernoulli beta hyperparameter
// @param beta_bernoulli_alpha_between SBM between-cluster alpha
// @param beta_bernoulli_beta_between  SBM between-cluster beta
// @param dirichlet_alpha         Dirichlet alpha for SBM
// @param lambda                  Lambda for SBM
// @param sampler_type            Sampler type string ("adaptive-metropolis" or "nuts")
// @param target_acceptance       Target acceptance rate for gradient-based samplers
// @param max_tree_depth          Maximum tree depth for NUTS
// @param na_impute               Whether to impute missing data
// @param missing_index_discrete  Matrix of missing discrete indices (n_miss x 2, 0-based)
// @param missing_index_continuous Matrix of missing continuous indices (n_miss x 2, 0-based)
// @param delta                   Determinant-tilt exponent on |Kyy|
// @param edge_prior_correction   Normalizing-constant correction list for the
//                                hierarchical edge priors (see
//                                R/correction_tables.R), or R_NilValue
//
// @return List with per-chain results including samples and diagnostics
// [[Rcpp::export]]
Rcpp::List sample_mixed_mrf(
    const Rcpp::List& inputFromR,
    const arma::mat& prior_inclusion_prob,
    const arma::imat& initial_edge_indicators,
    const int no_iter,
    const int no_warmup,
    const int no_chains,
    const bool edge_selection,
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
    const std::string& sampler_type = "adaptive-metropolis",
    const double target_acceptance = 0.80,
    const int max_tree_depth = 10,
    const bool na_impute = false,
    const Rcpp::Nullable<Rcpp::IntegerMatrix> missing_index_discrete_nullable = R_NilValue,
    const Rcpp::Nullable<Rcpp::IntegerMatrix> missing_index_continuous_nullable = R_NilValue,
    const double delta = 0.0,
    const Rcpp::Nullable<Rcpp::List> edge_prior_correction = R_NilValue
) {
    // Extract model inputs from R list
    arma::imat discrete_obs = Rcpp::as<arma::imat>(inputFromR["discrete_observations"]);
    arma::mat continuous_obs = Rcpp::as<arma::mat>(inputFromR["continuous_observations"]);
    arma::ivec num_categories = Rcpp::as<arma::ivec>(inputFromR["num_categories"]);
    arma::uvec is_ordinal = Rcpp::as<arma::uvec>(inputFromR["is_ordinal_variable"]);
    arma::ivec baseline_cat = Rcpp::as<arma::ivec>(inputFromR["baseline_category"]);

    // Create parameter priors from R input
    double pairwise_scale = Rcpp::as<double>(inputFromR["pairwise_scale"]);
    std::string ipt_str = inputFromR.containsElementNamed("interaction_prior_type")
        ? Rcpp::as<std::string>(inputFromR["interaction_prior_type"]) : "cauchy";
    double ia = inputFromR.containsElementNamed("interaction_alpha")
        ? Rcpp::as<double>(inputFromR["interaction_alpha"]) : NA_REAL;
    double ib = inputFromR.containsElementNamed("interaction_beta")
        ? Rcpp::as<double>(inputFromR["interaction_beta"]) : NA_REAL;
    auto interaction_prior = create_parameter_prior(ipt_str, pairwise_scale, ia, ib);

    // Threshold prior
    std::string tpt_str = inputFromR.containsElementNamed("threshold_prior_type")
        ? Rcpp::as<std::string>(inputFromR["threshold_prior_type"]) : "beta-prime";
    double ta = inputFromR.containsElementNamed("main_alpha")
        ? Rcpp::as<double>(inputFromR["main_alpha"]) : 0.5;
    double tb = inputFromR.containsElementNamed("main_beta")
        ? Rcpp::as<double>(inputFromR["main_beta"]) : 0.5;
    double ts = inputFromR.containsElementNamed("threshold_scale")
        ? Rcpp::as<double>(inputFromR["threshold_scale"]) : 1.0;
    auto threshold_prior = create_parameter_prior(tpt_str, ts, ta, tb);

    // Means prior (continuous means)
    std::string mpt_str = inputFromR.containsElementNamed("means_prior_type")
        ? Rcpp::as<std::string>(inputFromR["means_prior_type"]) : "normal";
    double ms = inputFromR.containsElementNamed("means_scale")
        ? Rcpp::as<double>(inputFromR["means_scale"]) : 1.0;
    double ma = inputFromR.containsElementNamed("means_alpha")
        ? Rcpp::as<double>(inputFromR["means_alpha"]) : NA_REAL;
    double mb = inputFromR.containsElementNamed("means_beta")
        ? Rcpp::as<double>(inputFromR["means_beta"]) : NA_REAL;
    auto means_prior = create_parameter_prior(mpt_str, ms, ma, mb);

    // Diagonal prior (precision diagonal)
    std::string spt_str = inputFromR.containsElementNamed("scale_prior_type")
        ? Rcpp::as<std::string>(inputFromR["scale_prior_type"]) : "gamma";
    double s_shape = inputFromR.containsElementNamed("scale_shape")
        ? Rcpp::as<double>(inputFromR["scale_shape"]) : 1.0;
    double s_rate = inputFromR.containsElementNamed("scale_rate")
        ? Rcpp::as<double>(inputFromR["scale_rate"]) : 1.0;
    auto diagonal_prior = create_scale_prior(spt_str, s_shape, s_rate);

    // Create model
    MixedMRFModel model(
        discrete_obs, continuous_obs,
        num_categories, is_ordinal, baseline_cat,
        prior_inclusion_prob, initial_edge_indicators,
        edge_selection,
        std::move(interaction_prior), std::move(threshold_prior),
        std::move(means_prior), std::move(diagonal_prior),
        seed
    );

    // Forward target_accept to the model's MH proposal-SD tuner.
    //   - Under "adaptive-metropolis": user's target_accept goes through
    //     directly (default 0.44 = componentwise RW MH optimum).
    //   - Under "nuts": user's target_accept (default 0.80) is the
    //     HMC step-size dual-averaging target and should NOT govern the
    //     between-model MH proposal SDs, which are still 1-D componentwise
    //     RW MH. Hardcode 0.44 there to keep stage-3b RM on the right
    //     fixed point.
    const double mh_target = (sampler_type == "nuts") ? 0.44 : target_acceptance;
    model.set_metropolis_target_accept(mh_target);

    // Determinant-tilt prior on |Kyy|: shifts both NUTS and MH targets by
    // delta * log|Kyy|. delta = 0 is the default (untilted). Consumed by
    // both gradient paths and the continuous-block MH ratios in
    // MixedMRFModel.
    model.set_determinant_tilt_yy(delta);

    // Set up missing data imputation
    if(na_impute) {
        arma::imat missing_disc, missing_cont;
        if(missing_index_discrete_nullable.isNotNull()) {
            missing_disc = Rcpp::as<arma::imat>(
                Rcpp::IntegerMatrix(missing_index_discrete_nullable.get()));
        }
        if(missing_index_continuous_nullable.isNotNull()) {
            missing_cont = Rcpp::as<arma::imat>(
                Rcpp::IntegerMatrix(missing_index_continuous_nullable.get()));
        }
        model.set_missing_data(missing_disc, missing_cont);
    }

    // Create edge prior
    EdgePrior edge_prior_enum = edge_prior_from_string(edge_prior);
    auto edge_prior_obj = create_edge_prior(
        edge_prior_enum,
        beta_bernoulli_alpha, beta_bernoulli_beta,
        beta_bernoulli_alpha_between, beta_bernoulli_beta_between,
        dirichlet_alpha, lambda
    );

    // Attach the normalizing-constant correction for the |Kyy| tilt so the
    // hyperparameter updates target the corrected conditionals. The list's
    // is_continuous mask restricts the tilt terms to continuous-continuous
    // pairs; the table is built for the continuous block alone.
    attach_edge_prior_correction(
        edge_prior_obj.get(), edge_prior_correction, "sample_mixed_mrf");

    // Configure sampler
    SamplerConfig config;
    config.sampler_type = sampler_type;
    config.no_iter = no_iter;
    config.no_warmup = no_warmup;
    config.edge_selection = edge_selection;
    config.seed = seed;
    config.na_impute = na_impute;
    config.target_acceptance = target_acceptance;
    config.max_tree_depth = max_tree_depth;

    // Set up progress manager
    ProgressManager pm(no_chains, no_iter, no_warmup, 50, progress_type, true, progress_callback);

    // Run MCMC using unified infrastructure
    std::vector<ChainResult> results = run_mcmc_sampler(
        model, *edge_prior_obj, config, no_chains, no_threads, pm);

    // Convert to R list format
    Rcpp::List output = convert_results_to_list(results);

    pm.finish();

    return output;
}

// sample_omrf.cpp - R interface for OMRF model sampling
//
// Uses the unified MCMC runner infrastructure to sample from OMRF models.
// Supports MH and NUTS samplers with optional edge selection.
#include <vector>
#include <memory>
#include <RcppArmadillo.h>

#include "models/omrf/omrf_model.h"
#include "utils/progress_manager.h"
#include "utils/common_helpers.h"
#include "priors/edge_prior.h"
#include "priors/parameter_prior.h"
#include "mcmc/execution/chain_result.h"
#include "mcmc/execution/chain_runner.h"
#include "mcmc/execution/sampler_config.h"

// R-exported function to sample from an OMRF model
//
// @param inputFromR          List with model specification
// @param prior_inclusion_prob Prior inclusion probabilities (p x p matrix)
// @param initial_edge_indicators Initial edge indicators (p x p integer matrix)
// @param no_iter             Number of post-warmup iterations
// @param no_warmup           Number of warmup iterations
// @param no_chains           Number of parallel chains
// @param edge_selection      Whether to do edge selection (spike-and-slab)
// @param sampler_type        "adaptive-metropolis" or "nuts"
// @param seed                Random seed
// @param no_threads          Number of threads for parallel execution
// @param progress_type       Progress bar type
// @param progress_callback   R function (SEXP) called as callback(completed, total) at regular intervals, or R_NilValue
// @param edge_prior          Edge prior type: "Bernoulli", "Beta-Bernoulli", "Stochastic-Block"
// @param na_impute           Whether to impute missing data
// @param missing_index       Matrix of missing data indices (n_missing x 2, 0-based)
// @param beta_bernoulli_alpha     Beta-Bernoulli alpha hyperparameter
// @param beta_bernoulli_beta      Beta-Bernoulli beta hyperparameter
// @param beta_bernoulli_alpha_between SBM between-cluster alpha
// @param beta_bernoulli_beta_between  SBM between-cluster beta
// @param dirichlet_alpha     Dirichlet alpha for SBM
// @param lambda              Lambda for SBM
// @param target_acceptance   Target acceptance rate for NUTS (default: 0.8)
// @param max_tree_depth      Maximum tree depth for NUTS (default: 10)
//
// @return List with per-chain results including samples and diagnostics
// [[Rcpp::export]]
Rcpp::List sample_omrf(
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
    const bool na_impute = false,
    const Rcpp::Nullable<Rcpp::IntegerMatrix> missing_index_nullable = R_NilValue,
    const double beta_bernoulli_alpha = 1.0,
    const double beta_bernoulli_beta = 1.0,
    const double beta_bernoulli_alpha_between = 1.0,
    const double beta_bernoulli_beta_between = 1.0,
    const double dirichlet_alpha = 1.0,
    const double lambda = 1.0,
    const double target_acceptance = 0.8,
    const int max_tree_depth = 10,
    const bool learn_mass_matrix = true,
    const Rcpp::Nullable<Rcpp::NumericMatrix> pairwise_scaling_factors_nullable = R_NilValue
) {
    // Create parameter priors from R input
    double pairwise_scale = Rcpp::as<double>(inputFromR["pairwise_scale"]);
    std::string ipt_str = inputFromR.containsElementNamed("interaction_prior_type")
        ? Rcpp::as<std::string>(inputFromR["interaction_prior_type"]) : "cauchy";
    double ia = inputFromR.containsElementNamed("interaction_alpha")
        ? Rcpp::as<double>(inputFromR["interaction_alpha"]) : NA_REAL;
    double ib = inputFromR.containsElementNamed("interaction_beta")
        ? Rcpp::as<double>(inputFromR["interaction_beta"]) : NA_REAL;
    auto interaction_prior = create_parameter_prior(ipt_str, pairwise_scale, ia, ib);

    std::string tpt_str = inputFromR.containsElementNamed("threshold_prior_type")
        ? Rcpp::as<std::string>(inputFromR["threshold_prior_type"]) : "beta-prime";
    double ta = inputFromR.containsElementNamed("main_alpha")
        ? Rcpp::as<double>(inputFromR["main_alpha"]) : 0.5;
    double tb = inputFromR.containsElementNamed("main_beta")
        ? Rcpp::as<double>(inputFromR["main_beta"]) : 0.5;
    double ts = inputFromR.containsElementNamed("threshold_scale")
        ? Rcpp::as<double>(inputFromR["threshold_scale"]) : 1.0;
    auto threshold_prior = create_parameter_prior(tpt_str, ts, ta, tb);

    // Create model from R input
    OMRFModel model = createOMRFModelFromR(
        inputFromR, prior_inclusion_prob, initial_edge_indicators,
        std::move(interaction_prior), std::move(threshold_prior),
        edge_selection);

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

    // Set pairwise scaling factors (if provided)
    if (pairwise_scaling_factors_nullable.isNotNull()) {
        arma::mat sf = Rcpp::as<arma::mat>(
            Rcpp::NumericMatrix(pairwise_scaling_factors_nullable.get()));
        model.set_pairwise_scaling_factors(sf);
    }

    // Set up missing data imputation
    if (na_impute && missing_index_nullable.isNotNull()) {
        arma::imat missing_index = Rcpp::as<arma::imat>(
            Rcpp::IntegerMatrix(missing_index_nullable.get()));
        model.set_missing_data(missing_index);
    }

    // Create edge prior
    EdgePrior edge_prior_enum = edge_prior_from_string(edge_prior);
    auto edge_prior_obj = create_edge_prior(
        edge_prior_enum,
        beta_bernoulli_alpha, beta_bernoulli_beta,
        beta_bernoulli_alpha_between, beta_bernoulli_beta_between,
        dirichlet_alpha, lambda
    );

    // Configure sampler
    SamplerConfig config;
    config.sampler_type = sampler_type;
    config.no_iter = no_iter;
    config.no_warmup = no_warmup;
    config.edge_selection = edge_selection;
    config.seed = seed;
    config.target_acceptance = target_acceptance;
    config.max_tree_depth = max_tree_depth;
    config.learn_mass_matrix = learn_mass_matrix;
    config.na_impute = na_impute;

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

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

// R-exported function to sample from an OMRF model. Takes the model
// specification list, p x p prior inclusion probabilities and initial edge
// indicators, iteration/warmup/chain counts, the sampler type
// ("adaptive-metropolis" or "nuts") with target acceptance and max tree depth,
// the edge prior ("Bernoulli", "Beta-Bernoulli", "Stochastic-Block") with its
// Beta-Bernoulli/SBM hyperparameters, missing-data options (missing_index:
// n_missing x 2, 0-based), seed, thread count, and progress settings
// (progress_callback is called as callback(completed, total), or R_NilValue).
// Returns a list of per-chain results with samples and diagnostics.
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
    const Rcpp::Nullable<Rcpp::List> initial_parameters = R_NilValue,
    const Rcpp::Nullable<Rcpp::NumericVector> initial_step_sizes = R_NilValue,
    const Rcpp::Nullable<Rcpp::List> initial_inv_mass = R_NilValue
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

    // Optional random interaction-slab-scale hyperprior. An empty type string
    // means the scale is fixed (the default). Otherwise the multiplier u = s/s0
    // follows a mean-1 gamma/exponential (shape/rate from R), and the scale
    // becomes a sampled parameter.
    std::string ispt_str = inputFromR.containsElementNamed("interaction_scale_prior_type")
        ? Rcpp::as<std::string>(inputFromR["interaction_scale_prior_type"]) : "";
    if (!ispt_str.empty()) {
        double is_shape = Rcpp::as<double>(inputFromR["interaction_scale_shape"]);
        double is_rate = Rcpp::as<double>(inputFromR["interaction_scale_rate"]);
        model.enable_random_interaction_scale(
            create_scale_prior(ispt_str, is_shape, is_rate));
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

    // Optional per-chain warm start (final state of a previous fit).
    std::vector<arma::vec> init_params;
    if (initial_parameters.isNotNull()) {
        Rcpp::List ip(initial_parameters.get());
        init_params.reserve(ip.size());
        for (int c = 0; c < ip.size(); ++c) {
            init_params.push_back(Rcpp::as<arma::vec>(ip[c]));
        }
    }
    std::vector<double> init_step_sizes;
    if (initial_step_sizes.isNotNull()) {
        init_step_sizes = Rcpp::as<std::vector<double>>(
            Rcpp::NumericVector(initial_step_sizes.get()));
    }
    std::vector<arma::vec> init_inv_mass;
    if (initial_inv_mass.isNotNull()) {
        Rcpp::List im(initial_inv_mass.get());
        init_inv_mass.reserve(im.size());
        for (int c = 0; c < im.size(); ++c) {
            init_inv_mass.push_back(Rcpp::as<arma::vec>(im[c]));
        }
    }

    // A warm-start list is per chain: an empty list means cold, otherwise it
    // must carry exactly one entry per chain (each is indexed by chain id).
    auto require_per_chain = [&](std::size_t n, const char* what) {
        if (n != 0 && n != static_cast<std::size_t>(no_chains)) {
            Rcpp::stop("%s has %d entries but there are %d chains; a warm-start "
                       "list must be empty or one entry per chain.",
                       what, static_cast<int>(n), no_chains);
        }
    };
    require_per_chain(init_params.size(), "initial_parameters");
    require_per_chain(init_step_sizes.size(), "initial_step_sizes");
    require_per_chain(init_inv_mass.size(), "initial_inv_mass");

    // Run MCMC using unified infrastructure
    std::vector<ChainResult> results = run_mcmc_sampler(
        model, *edge_prior_obj, config, no_chains, no_threads, pm,
        init_params, init_step_sizes, init_inv_mass);

    // Convert to R list format
    Rcpp::List output = convert_results_to_list(results);

    pm.finish();

    return output;
}

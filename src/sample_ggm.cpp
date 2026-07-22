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
    const bool learn_mass_matrix = true,
    const bool na_impute = false,
    const Rcpp::Nullable<Rcpp::IntegerMatrix> missing_index_nullable = R_NilValue,
    const double delta = 0.0,
    const Rcpp::Nullable<Rcpp::List> edge_prior_correction = R_NilValue,
    const Rcpp::Nullable<Rcpp::List> zratio_spec = R_NilValue
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

    // Forward target_accept to the model's between-model MH proposal-SD tuner.
    // Only adaptive-metropolis and nuts run the componentwise RW edge move that
    // consumes it; the gibbs within-step and its full-conditional edge move are
    // exact and tune nothing, so the gibbs path must not set an MH target.
    //   - Under "adaptive-metropolis": user's target_accept goes through
    //     directly (default 0.44 = componentwise RW MH optimum).
    //   - Under "nuts": user's target_accept (default 0.80) is the
    //     HMC step-size dual-averaging target and should NOT govern the
    //     between-model MH proposal SDs, which are still 1-D componentwise
    //     RW MH. Use 0.44 there to keep stage-3b RM on the right fixed point.
    if (sampler_type != "gibbs") {
        const double mh_target =
            (sampler_type == "adaptive-metropolis") ? target_acceptance : 0.44;
        model.set_metropolis_target_accept(mh_target);
    }

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


    // Hierarchical prior specification: attach the per-edge Z-ratio engine
    // so the between-edge moves target p(K | Gamma) = rho_Gamma(K)/Z(Gamma).
    // The constants are resolved at R spec-build (zratio_constants); each
    // chain clone deep-copies the engine with its cache.
    int zratio_window = 0;
    int zratio_gauge_sweeps = 0;
    if (zratio_spec.isNotNull()) {
        Rcpp::List zs(zratio_spec.get());
        auto engine = std::make_shared<ZRatioEngine>(
            Rcpp::as<arma::vec>(zs["addc"]),
            Rcpp::as<arma::vec>(zs["tg"]),
            Rcpp::as<arma::vec>(zs["ihat"]),
            Rcpp::as<arma::vec>(zs["ghat"]),
            Rcpp::as<arma::vec>(zs["wt"]),
            Rcpp::as<double>(zs["psi0"]));
        if (zs.containsElementNamed("calibration_window")) {
            zratio_window = Rcpp::as<int>(zs["calibration_window"]);
        }
        if (zs.containsElementNamed("gauge_sweeps")) {
            zratio_gauge_sweeps = Rcpp::as<int>(zs["gauge_sweeps"]);
        }
        bool zr_cauchy = zs.containsElementNamed("slab") &&
            Rcpp::as<std::string>(zs["slab"]) == "cauchy";
        const double zr_delta = Rcpp::as<double>(zs["delta"]);
        const double zr_eta = Rcpp::as<double>(zs["eta"]);
        const double zr_alpha = zs.containsElementNamed("alpha")
            ? Rcpp::as<double>(zs["alpha"]) : 1.0;
        // Option-B surfaces (built once in R at the analysis eta) replace the
        // online OLS correction: attach them so log_zratio decomposes each
        // block and sums per-component surface moments. run_sampler zeros the
        // calibration window when the surface is present, so the OLS path stays
        // in place but dormant.
        if (zs.containsElementNamed("surface") && !Rf_isNull(zs["surface"])) {
            Rcpp::List zsurf(zs["surface"]);
            engine->set_surface(surface_family_from_list(zsurf["cn"]),
                                surface_family_from_list(zsurf["bip"]));
        }
        // The rng pointer is rebound per chain clone by GGMModel.
        if (zratio_window > 0) {
            engine->enable_calibration(zr_delta, zr_eta, nullptr, 100, 30, 1.0,
                                       6, zr_cauchy, zr_alpha, 100);
        } else {
            // No warm-up calibration (pre-packed constants): still hand the
            // engine the standardized-cell prior params so the trust gauge's
            // block-local reference can sample.
            engine->set_oracle_params(zr_delta, zr_eta, nullptr, 300, 30,
                                      zr_cauchy, zr_alpha);
        }
        model.set_zratio_engine(std::move(engine));
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
    config.learn_mass_matrix = learn_mass_matrix;
    config.na_impute = na_impute;
    config.zratio_calibration_window = zratio_window;
    config.zratio_gauge_sweeps = zratio_gauge_sweeps;

    // Set up progress manager
    ProgressManager pm(no_chains, no_iter, no_warmup + zratio_window, 50, progress_type, true, progress_callback);

    // Create edge prior
    EdgePrior edge_prior_enum = edge_prior_from_string(edge_prior);
    auto edge_prior_obj = create_edge_prior(
        edge_prior_enum,
        beta_bernoulli_alpha, beta_bernoulli_beta,
        beta_bernoulli_alpha_between, beta_bernoulli_beta_between,
        dirichlet_alpha, lambda
    );

    // Attach the normalizing-constant correction (curves built from the
    // tilted prior sampler at fit setup) so the hyperparameter updates
    // target the corrected conditionals.
    attach_edge_prior_correction(
        edge_prior_obj.get(), edge_prior_correction, "sample_ggm");

    // Run MCMC using unified infrastructure
    std::vector<ChainResult> results = run_mcmc_sampler(
        model, *edge_prior_obj, config, no_chains, no_threads, pm);

    // Convert to R list format
    Rcpp::List output = convert_results_to_list(results);

    pm.finish();

    return output;
}
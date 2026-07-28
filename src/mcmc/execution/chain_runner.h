#pragma once

#include <vector>
#include <memory>
#include <limits>
#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include "models/base_model.h"
#include "mcmc/execution/chain_result.h"
#include "priors/edge_prior.h"
#include "utils/progress_manager.h"
#include "mcmc/execution/sampler_config.h"
#include "mcmc/samplers/sampler_base.h"
#include "mcmc/execution/warmup_schedule.h"


/** Which concrete sampler a run uses. */
enum class SamplerKind { NUTS, AdaptiveMetropolis, Gibbs };

/**
 * Behavioral descriptor for a sampler type. resolve_sampler_spec is the single
 * place that decodes the config.sampler_type string; all downstream code reads
 * these fields instead of re-comparing the string.
 */
struct SamplerSpec {
    SamplerKind kind;
    bool learn_sd;    ///< NUTS dual-averaging step-size adaptation during warmup
    bool nuts_diag;   ///< reserve/store NUTS per-iteration diagnostics
    bool am_diag;     ///< reserve/store adaptive-Metropolis diagnostics
};

/**
 * Decode a sampler-type string into a SamplerSpec.
 *
 * @param sampler_type  "nuts" or "adaptive-metropolis".
 * @return Descriptor with the concrete kind and its derived behavior flags.
 */
SamplerSpec resolve_sampler_spec(const std::string& sampler_type);

/**
 * Create a sampler of the given kind.
 *
 * @param kind     Concrete sampler kind (from resolve_sampler_spec).
 * @param config   Sampler configuration (step size, tree depth, etc.)
 * @param schedule Shared warmup schedule for adaptation staging
 * @return Owning pointer to a concrete SamplerBase subclass
 */
std::unique_ptr<SamplerBase> create_sampler(SamplerKind kind, const SamplerConfig& config, WarmupSchedule& schedule);

/**
 * Run a single MCMC chain (warmup + sampling)
 *
 * @param chain_result  Output container for samples and diagnostics
 * @param model         Model to sample from (state is mutated in place)
 * @param edge_prior    Edge prior for edge-selection updates
 * @param config        Sampler configuration
 * @param chain_id      Chain identifier (0-based)
 * @param pm            Progress manager for user-facing status updates
 */
void run_mcmc_chain(
    ChainResult& chain_result,
    BaseModel& model,
    BaseEdgePrior& edge_prior,
    const SamplerConfig& config,
    int chain_id,
    ProgressManager& pm,
    double warm_step_size = std::numeric_limits<double>::quiet_NaN(),
    const arma::vec& warm_inv_mass = arma::vec()
);


/**
 * MCMCChainRunner - TBB worker for parallel chain execution
 *
 * Each chain gets its own model clone and edge prior. The worker dispatches
 * chains to threads via RcppParallel::parallelFor.
 */
struct MCMCChainRunner : public RcppParallel::Worker {
    std::vector<ChainResult>& results_;
    std::vector<std::unique_ptr<BaseModel>>& models_;
    std::vector<std::unique_ptr<BaseEdgePrior>>& edge_priors_;
    const SamplerConfig& config_;
    ProgressManager& pm_;
    const std::vector<double>& warm_step_sizes_;
    const std::vector<arma::vec>& warm_inv_masses_;

    MCMCChainRunner(
        std::vector<ChainResult>& results,
        std::vector<std::unique_ptr<BaseModel>>& models,
        std::vector<std::unique_ptr<BaseEdgePrior>>& edge_priors,
        const SamplerConfig& config,
        ProgressManager& pm,
        const std::vector<double>& warm_step_sizes,
        const std::vector<arma::vec>& warm_inv_masses
    ) :
        results_(results),
        models_(models),
        edge_priors_(edge_priors),
        config_(config),
        pm_(pm),
        warm_step_sizes_(warm_step_sizes),
        warm_inv_masses_(warm_inv_masses)
    {}

    void operator()(std::size_t begin, std::size_t end);
};


/**
 * Run multi-chain MCMC (parallel or sequential based on thread count)
 *
 * @param model       Prototype model (cloned per chain)
 * @param edge_prior  Prototype edge prior (cloned per chain)
 * @param config      Sampler configuration
 * @param no_chains   Number of chains to run
 * @param no_threads  Number of threads (1 = sequential)
 * @param pm          Progress manager for user-facing status updates
 * @param initial_parameters  Optional per-chain warm-start storage vectors
 *                            (empty = cold start). Applied to each chain's
 *                            model clone via set_storage_vectorized_parameters.
 * @param initial_step_sizes  Optional per-chain warm-start NUTS step sizes
 *                            (empty = heuristic). Passed to set_warm_step_size.
 * @param initial_inv_mass    Optional per-chain warm-start diagonal metric
 *                            (empty = cold). Passed to set_warm_inv_mass; NUTS
 *                            holds it fixed with no windowed re-adaptation.
 * @return Vector of ChainResult, one per chain
 */
std::vector<ChainResult> run_mcmc_sampler(
    BaseModel& model,
    BaseEdgePrior& edge_prior,
    const SamplerConfig& config,
    int no_chains,
    int no_threads,
    ProgressManager& pm,
    const std::vector<arma::vec>& initial_parameters = {},
    const std::vector<double>& initial_step_sizes = {},
    const std::vector<arma::vec>& initial_inv_mass = {}
);

/**
 * Convert chain results to an Rcpp::List for return to R
 *
 * @param results  Vector of completed chain results
 * @return Named list with samples, diagnostics, and metadata
 */
Rcpp::List convert_results_to_list(const std::vector<ChainResult>& results);

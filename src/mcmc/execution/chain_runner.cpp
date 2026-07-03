#include "mcmc/execution/chain_runner.h"

#include <exception>
#include <stdexcept>
#include <tbb/global_control.h>
#include "mcmc/samplers/nuts_sampler.h"
#include "mcmc/samplers/metropolis_sampler.h"
#include "mcmc/samplers/gibbs_sampler.h"


namespace {

// Store this sample's NUTS per-iteration diagnostics, if the chain is collecting
// them and the step produced a NUTSDiagnostics object. Keeps the dynamic_cast
// and the wide store_nuts_diagnostics call out of the main sampling loop.
void store_nuts_diagnostics_if_present(ChainResult& chain_result, int sample_index,
                                       const SamplerBase& sampler, const StepResult& result) {
    if (!chain_result.has_nuts_diagnostics || !sampler.has_nuts_diagnostics()) return;
    auto* diag = dynamic_cast<NUTSDiagnostics*>(result.diagnostics.get());
    if (diag) {
        chain_result.store_nuts_diagnostics(sample_index, diag->tree_depth, diag->divergent,
                                            diag->non_reversible, diag->energy, diag->accept_prob);
    }
}

}  // namespace


SamplerSpec resolve_sampler_spec(const std::string& sampler_type) {
    if (sampler_type == "nuts") {
        return SamplerSpec{SamplerKind::NUTS, /*learn_sd=*/true, /*nuts_diag=*/true, /*am_diag=*/false};
    } else if (sampler_type == "adaptive-metropolis") {
        return SamplerSpec{SamplerKind::AdaptiveMetropolis, /*learn_sd=*/false, /*nuts_diag=*/false, /*am_diag=*/true};
    } else if (sampler_type == "gibbs") {
        return SamplerSpec{SamplerKind::Gibbs, /*learn_sd=*/false, /*nuts_diag=*/false, /*am_diag=*/false};
    } else {
        // std::runtime_error rather than Rcpp::stop: this runs on worker
        // threads, where constructing an Rcpp exception is not safe.
        throw std::runtime_error("Unknown sampler_type: '" + sampler_type + "'");
    }
}

std::unique_ptr<SamplerBase> create_sampler(SamplerKind kind, const SamplerConfig& config, WarmupSchedule& schedule) {
    switch (kind) {
        case SamplerKind::NUTS:
            return std::make_unique<NUTSSampler>(config, schedule);
        case SamplerKind::AdaptiveMetropolis:
            return std::make_unique<MetropolisSampler>(config, schedule);
        case SamplerKind::Gibbs:
            return std::make_unique<GibbsSampler>(config, schedule);
    }
    throw std::runtime_error("Unhandled SamplerKind");  // unreachable: kind comes from resolve_sampler_spec
}


void run_mcmc_chain(
    ChainResult& chain_result,
    BaseModel& model,
    BaseEdgePrior& edge_prior,
    const SamplerConfig& config,
    const int chain_id,
    ProgressManager& pm
) {
    chain_result.chain_id = chain_id + 1;

    // Construct warmup schedule (shared by runner and sampler)
    const SamplerSpec spec = resolve_sampler_spec(config.sampler_type);
    WarmupSchedule schedule(config.no_warmup, config.edge_selection, spec.learn_sd,
                            /*select_during_warmup=*/spec.kind == SamplerKind::Gibbs);

    auto sampler = create_sampler(spec.kind, config, schedule);

    // Initialize sampler (step-size heuristic) before the main loop
    sampler->initialize(model);

    const int total_iter = config.no_warmup + config.no_iter;

    // ---- Main MCMC loop (warmup + sampling) ----
    for (int iter = 0; iter < total_iter; ++iter) {

        // Per-iteration preparation (e.g., shuffle edge order)
        model.prepare_iteration();

        // Optional missing-data imputation
        if (config.na_impute && model.has_missing_data()) {
            model.impute_missing();
        }

        // Edge selection
        if (schedule.selection_enabled(iter) && model.has_edge_selection()) {
            if (iter == schedule.stage3c_start) {
                model.set_edge_selection_active(true);
            }
            model.update_edge_indicators();
        }

        // Main parameter update — adaptation is internal to sampler
        StepResult result = sampler->step(model, iter);

        // Stage 3b: proposal-SD tuning
        model.tune_proposal_sd(iter, schedule);

        // Edge prior update
        if (schedule.selection_enabled(iter) && model.has_edge_selection()) {
            edge_prior.update(
                model.get_edge_indicators(),
                model.get_inclusion_probability(),
                model.get_num_variables(),
                model.get_num_pairwise(),
                model.get_rng()
            );
        }

        // Store samples (only during sampling phase)
        if (schedule.sampling(iter)) {
            int sample_index = iter - config.no_warmup;

            store_nuts_diagnostics_if_present(chain_result, sample_index, *sampler, result);

            if (chain_result.has_am_diagnostics) {
                chain_result.store_am_diagnostics(sample_index, result.accept_prob);
            }

            chain_result.store_sample(sample_index, model.get_storage_vectorized_parameters());

            if (chain_result.has_indicators) {
                chain_result.store_indicators(sample_index, model.get_vectorized_indicator_parameters());
            }

            if (chain_result.has_allocations && edge_prior.has_allocations()) {
                chain_result.store_allocations(sample_index, edge_prior.get_allocations());
            }

            if (chain_result.has_inclusion_parameter) {
                chain_result.store_inclusion_parameter(
                    sample_index, edge_prior.get_inclusion_parameter());
            }
        }

        pm.update(chain_id);
        if (pm.shouldExit()) {
            chain_result.userInterrupt = true;
            return;
        }
    }

}


void MCMCChainRunner::operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
        ChainResult& chain_result = results_[i];
        BaseModel& model = *models_[i];
        BaseEdgePrior& edge_prior = *edge_priors_[i];
        model.set_seed(config_.seed + static_cast<int>(i));

        try {
            run_mcmc_chain(chain_result, model, edge_prior, config_, static_cast<int>(i), pm_);
        } catch (std::exception& e) {
            chain_result.error = true;
            chain_result.error_msg = e.what();
        } catch (...) {
            chain_result.error = true;
            chain_result.error_msg = "Unknown error";
        }
    }
}


std::vector<ChainResult> run_mcmc_sampler(
    BaseModel& model,
    BaseEdgePrior& edge_prior,
    const SamplerConfig& config,
    const int no_chains,
    const int no_threads,
    ProgressManager& pm
) {
    const SamplerSpec spec = resolve_sampler_spec(config.sampler_type);
    const bool has_nuts_diag = spec.nuts_diag;
    const bool has_am_diag = spec.am_diag;
    const bool has_sbm_alloc = edge_prior.has_allocations() ||
        (config.edge_selection && dynamic_cast<StochasticBlockEdgePrior*>(&edge_prior) != nullptr);

    std::vector<ChainResult> results(no_chains);
    for (int c = 0; c < no_chains; ++c) {
        results[c].reserve(model.storage_dimension(), config.no_iter);

        if (config.edge_selection) {
            size_t n_edges = model.get_vectorized_indicator_parameters().n_elem;
            results[c].reserve_indicators(n_edges, config.no_iter);
        }

        if (has_sbm_alloc) {
            results[c].reserve_allocations(model.get_num_variables(), config.no_iter);
        }

        if (config.edge_selection && edge_prior.has_inclusion_parameter()) {
            results[c].reserve_inclusion_parameter(config.no_iter);
        }

        if (has_nuts_diag) {
            results[c].reserve_nuts_diagnostics(config.no_iter);
        }

        if (has_am_diag) {
            results[c].reserve_am_diagnostics(config.no_iter);
        }
    }

    if (no_threads > 1) {
        std::vector<std::unique_ptr<BaseModel>> models;
        std::vector<std::unique_ptr<BaseEdgePrior>> edge_priors;
        models.reserve(no_chains);
        edge_priors.reserve(no_chains);
        for (int c = 0; c < no_chains; ++c) {
            models.push_back(model.clone());
            models[c]->set_seed(config.seed + c);
            edge_priors.push_back(edge_prior.clone());
        }

        MCMCChainRunner runner(results, models, edge_priors, config, pm);
        tbb::global_control control(tbb::global_control::max_allowed_parallelism, no_threads);
        RcppParallel::parallelFor(0, static_cast<size_t>(no_chains), runner);

    } else {
        model.set_seed(config.seed);
        for (int c = 0; c < no_chains; ++c) {
            auto chain_model = model.clone();
            chain_model->set_seed(config.seed + c);
            auto chain_edge_prior = edge_prior.clone();
            run_mcmc_chain(results[c], *chain_model, *chain_edge_prior, config, c, pm);
        }
    }

    return results;
}


Rcpp::List convert_results_to_list(const std::vector<ChainResult>& results) {
    Rcpp::List output(results.size());

    for (size_t i = 0; i < results.size(); ++i) {
        const ChainResult& chain = results[i];
        Rcpp::List chain_list;

        chain_list["chain_id"] = chain.chain_id;

        if (chain.error) {
            chain_list["error"] = true;
            chain_list["error_msg"] = chain.error_msg;
        } else {
            chain_list["error"] = false;
            chain_list["samples"] = chain.samples;
            chain_list["userInterrupt"] = chain.userInterrupt;

            if (chain.has_indicators) {
                chain_list["indicator_samples"] = chain.indicator_samples;
            }

            if (chain.has_allocations) {
                chain_list["allocation_samples"] = chain.allocation_samples;
            }

            if (chain.has_inclusion_parameter) {
                chain_list["inclusion_parameter_samples"] = chain.inclusion_parameter_samples;
            }

            if (chain.has_nuts_diagnostics) {
                chain_list["treedepth"] = chain.treedepth_samples;
                chain_list["divergent"] = chain.divergent_samples;
                chain_list["non_reversible"] = chain.non_reversible_samples;
                chain_list["energy"] = chain.energy_samples;
                chain_list["accept_prob"] = chain.accept_prob_samples;
            }

            if (chain.has_am_diagnostics) {
                chain_list["am_accept_prob"] = chain.am_accept_prob_samples;
            }
        }

        output[i] = chain_list;
    }

    return output;
}

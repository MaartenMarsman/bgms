#pragma once

#include "mcmc/samplers/sampler_base.h"
#include "mcmc/execution/sampler_config.h"
#include "mcmc/execution/warmup_schedule.h"
#include "models/base_model.h"

/**
 * GibbsSampler - conjugate Gibbs within-step for the GGM
 *
 * The within-step is an exact conjugate row sweep. With edge selection the
 * between-step (driven by the chain runner) uses the full-conditional edge
 * birth/death proposal, which needs no tuning, so the Gibbs sweep runs from the
 * first iteration with no warmup staging.
 *
 * The sweep logic (do_one_gibbs_step) is model-specific; this is a thin wrapper
 * providing the uniform sampler interface. Only constructed for models whose
 * full-conditionals are conjugate (currently the GGM); create_sampler gates it.
 */
class GibbsSampler : public SamplerBase {
public:
    /**
     * @param config    Sampler configuration (unused; the draw is exact).
     * @param schedule  Shared warmup schedule (unused; no adaptation).
     */
    GibbsSampler(const SamplerConfig& config, WarmupSchedule& schedule) {
        (void)config;
        (void)schedule;
    }

    /** Switch the model's between-step to the full-conditional edge proposal. */
    void initialize(BaseModel& model) override {
        model.set_conjugate_edge_proposal(true);
    }

    /** One conjugate Gibbs sweep; the exact draw always "accepts". */
    StepResult step(BaseModel& model, int iteration) override {
        model.do_one_gibbs_step(iteration);

        StepResult result;
        result.accept_prob = 1.0;
        return result;
    }
};

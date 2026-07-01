#pragma once

#include "mcmc/samplers/sampler_base.h"
#include "mcmc/execution/sampler_config.h"
#include "mcmc/execution/warmup_schedule.h"
#include "models/base_model.h"

/**
 * GibbsSampler - exact conjugate Gibbs sampler
 *
 * Delegates to the model's full-conditional Gibbs sweep. Because each draw is
 * exact there is no proposal-SD adaptation and no warmup tuning; this is a
 * thin wrapper that provides a uniform interface consistent with the other
 * samplers (Metropolis, NUTS), while the actual sweep logic
 * (do_one_gibbs_step) is model-specific.
 *
 * Only constructed for models whose full-conditionals are conjugate (currently
 * the GGM row-block Gibbs within-step); create_sampler gates this.
 */
class GibbsSampler : public SamplerBase {
public:
    /**
     * Construct a Gibbs sampler.
     * @param config    Sampler configuration (unused; the draw is exact).
     * @param schedule  Shared warmup schedule (unused; no adaptation).
     */
    GibbsSampler(const SamplerConfig& config, WarmupSchedule& schedule) {
        (void)config;
        (void)schedule;
    }

    /**
     * Perform one Gibbs sweep.
     *
     * The draw is exact, so the reported acceptance probability is always 1.
     */
    StepResult step(BaseModel& model, int iteration) override {
        model.do_one_gibbs_step(iteration);

        StepResult result;
        result.accept_prob = 1.0;
        return result;
    }
};

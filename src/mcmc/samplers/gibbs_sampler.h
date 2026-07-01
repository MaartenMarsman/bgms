#pragma once

#include "mcmc/samplers/sampler_base.h"
#include "mcmc/execution/sampler_config.h"
#include "mcmc/execution/warmup_schedule.h"
#include "models/base_model.h"

/**
 * GibbsSampler - conjugate Gibbs within-step for the GGM
 *
 * On a fixed graph the within-step is an exact conjugate row sweep from the
 * first iteration. With edge selection the between-step is a Roverato birth /
 * death move (driven by the chain runner) whose proposal SD needs tuning, so
 * warmup runs the adaptive-Metropolis within-step to tune those proposal SDs,
 * and production switches to the conjugate Gibbs sweep with the SDs frozen.
 *
 * The actual sweep logic (do_one_metropolis_step, do_one_gibbs_step) is
 * model-specific; this is a thin wrapper providing the uniform sampler
 * interface. Only constructed for models whose full-conditionals are
 * conjugate (currently the GGM); create_sampler gates this.
 */
class GibbsSampler : public SamplerBase {
public:
    /**
     * Construct a Gibbs sampler.
     * @param config    Sampler configuration; supplies the warmup length.
     * @param schedule  Shared warmup schedule (used to set up adaptation).
     */
    GibbsSampler(const SamplerConfig& config, WarmupSchedule& schedule)
        : schedule_(schedule), warmup_(config.no_warmup), initialized_(false) {}

    /**
     * Set up Metropolis adaptation, used to tune the between-step proposal SDs
     * during warmup when edge selection is active.
     */
    void initialize(BaseModel& model) override {
        if (initialized_) return;
        model.init_metropolis_adaptation(schedule_);
        initialized_ = true;
    }

    /**
     * One step. With edge selection, warmup uses the adaptive-Metropolis
     * within-step (tuning the proposal SDs the between-step reuses) and
     * production uses the exact conjugate Gibbs sweep. On a fixed graph the
     * Gibbs sweep is used throughout.
     */
    StepResult step(BaseModel& model, int iteration) override {
        if (!initialized_) initialize(model);

        StepResult result;
        if (model.has_edge_selection() && iteration < warmup_) {
            model.do_one_metropolis_step(iteration);
            result.accept_prob = model.last_metropolis_mean_accept_prob();
        } else {
            model.do_one_gibbs_step(iteration);
            result.accept_prob = 1.0;
        }
        return result;
    }

private:
    WarmupSchedule& schedule_;
    int warmup_;
    bool initialized_;
};

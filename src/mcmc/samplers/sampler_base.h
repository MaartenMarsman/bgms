#pragma once

#include <limits>
#include <RcppArmadillo.h>
#include "mcmc/execution/step_result.h"
#include "models/base_model.h"

// ---------------------------------------------------------------------------
// SamplerBase — abstract interface for all MCMC samplers
// ---------------------------------------------------------------------------

/**
 * SamplerBase - Abstract base class for MCMC samplers
 *
 * Provides a unified interface for all MCMC sampling algorithms:
 * - MetropolisSampler (component-wise random-walk Metropolis)
 * - NUTSSampler (No-U-Turn Sampler)
 *
 * The sampler internally decides whether to adapt based on the iteration
 * number and its warmup schedule reference.
 */
class SamplerBase {
public:
    virtual ~SamplerBase() = default;

    /**
     * Perform one MCMC step
     *
     * The sampler internally decides whether to adapt based on the
     * iteration number and its warmup schedule reference.
     *
     * @param model      The model to sample from
     * @param iteration  Current iteration (0-based, spans warmup + sampling)
     * @return StepResult with new state and diagnostics
     */
    virtual StepResult step(BaseModel& model, int iteration) = 0;

    /**
     * Initialize the sampler before the MCMC loop.
     * For gradient-based samplers, runs the step-size heuristic. Default no-op.
     */
    virtual void initialize(BaseModel& /*model*/) {}

    /**
     * Check if this sampler produces NUTS-style diagnostics
     * (tree depth, divergences, energy)
     */
    virtual bool has_nuts_diagnostics() const { return false; }

    /**
     * Warm-start the leapfrog step size (NUTS only): use `eps` as the initial
     * step size instead of the heuristic, with dual-averaging still live during
     * warmup. Non-gradient samplers ignore it. Default no-op.
     */
    virtual void set_warm_step_size(double /*eps*/) {}

    /**
     * @return The final (adaptation-averaged) leapfrog step size for a NUTS
     *         run, or NaN for samplers without one. Used to warm-start refits.
     */
    virtual double get_final_step_size() const {
        return std::numeric_limits<double>::quiet_NaN();
    }

    /**
     * Warm-start the diagonal inverse mass matrix (NUTS only): use `inv_mass`
     * (per-parameter variances, full theta layout) as the fixed metric for a
     * refit, with step-size dual averaging still live during the short warmup
     * and no windowed mass re-adaptation. Non-gradient samplers ignore it, and
     * an empty vector is a no-op. Default no-op.
     */
    virtual void set_warm_inv_mass(const arma::vec& /*inv_mass*/) {}

    /**
     * @return The final adapted diagonal inverse mass matrix for a NUTS run, or
     *         an empty vector for samplers without one. Used to warm-start refits.
     */
    virtual arma::vec get_final_inv_mass() const {
        return arma::vec();
    }
};

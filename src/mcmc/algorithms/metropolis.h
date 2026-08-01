#pragma once

#include <RcppArmadillo.h>
#include <functional>
#include "mcmc/execution/step_result.h"
struct SafeRNG;


/**
 * Performs one step of Random Walk Metropolis sampling for a scalar parameter
 *
 * Proposes from a symmetric normal distribution centered at the current state,
 * then accepts or rejects via the Metropolis-Hastings criterion.
 *
 * @param current_state  Current scalar parameter value
 * @param step_size      Standard deviation of the Gaussian proposal
 * @param log_post       Log-posterior function (scalar to scalar)
 * @param rng            Thread-safe random number generator
 * @return StepResult with accepted state (1-element vector) and acceptance probability
 */
StepResult metropolis_step(
    double current_state,
    double step_size,
    const std::function<double(double)>& log_post,
    SafeRNG& rng
);

/**
 * Performs one Random Walk Metropolis step with a cached current-state value
 *
 * Variant of metropolis_step() for callers that maintain the current-state
 * log-posterior across proposals: only the proposed state is evaluated. The
 * random-number stream matches metropolis_step() (one normal draw, one
 * uniform draw).
 *
 * @param current_state     Current scalar parameter value
 * @param step_size         Standard deviation of the Gaussian proposal
 * @param log_post_current  Log-posterior value at current_state
 * @param log_post_proposed Log-posterior function evaluated at the proposal
 * @param rng               Thread-safe random number generator
 * @return StepResult with accepted state (1-element vector) and acceptance probability
 */
StepResult metropolis_step_cached(
    double current_state,
    double step_size,
    double log_post_current,
    const std::function<double(double)>& log_post_proposed,
    SafeRNG& rng
);

#pragma once

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <memory>
#include "math/explog_macros.h"


/**
 * DiagnosticsBase - Polymorphic base for per-step diagnostics
 *
 * StepResult holds a shared_ptr<DiagnosticsBase> so the same return type
 * works for every sampler. Samplers that collect diagnostics (currently only
 * NUTS) define a derived struct; consumers downcast with dynamic_pointer_cast.
 * Samplers without diagnostics (Metropolis, HMC) leave the pointer null.
 */
struct DiagnosticsBase {
  virtual ~DiagnosticsBase() = default;
};



/**
 * NUTSDiagnostics - Per-iteration NUTS diagnostics (derives from DiagnosticsBase)
 */
struct NUTSDiagnostics : public DiagnosticsBase {
  int tree_depth;        ///< Depth of the trajectory tree
  bool divergent;        ///< Whether a divergence occurred
  double energy;         ///< Final Hamiltonian (-log posterior + kinetic energy)
  double accept_prob;    ///< Mean Metropolis acceptance over the trajectory's
                         ///  leapfrog steps (Stan's `accept_stat__`)
};



/**
 * StepResult - Outcome of one MCMC step
 *
 * Returned by SamplerBase::step(). The diagnostics pointer is non-null only
 * for NUTS (holds NUTSDiagnostics); Metropolis and HMC leave it null.
 */
struct StepResult {
  arma::vec state;                              ///< Accepted parameter vector
  double accept_prob;                           ///< Acceptance probability
  std::shared_ptr<DiagnosticsBase> diagnostics; ///< NUTS diagnostics, or null
};



/**
 * Robbins-Monro update for Metropolis proposal standard deviations
 *
 * Adjusts the proposal SD toward a target acceptance rate:
 *   sd += (observed_acceptance - target) * weight
 * The result is clamped to [0.001, 2.0]. NaN values are reset to 1.0.
 *
 * @param current_sd                          Current proposal standard deviation
 * @param observed_log_acceptance_probability Log acceptance probability from Metropolis step
 * @param rm_weight                           Robbins-Monro weight (e.g. iteration^{-0.75})
 * @param target_acceptance                   Target acceptance rate
 * @return Updated proposal SD, clamped to [0.001, 2.0]
 */
inline double update_proposal_sd_with_robbins_monro (
    const double current_sd,
    const double observed_log_acceptance_probability,
    const double rm_weight,
    const double target_acceptance
) {
  constexpr double rm_lower_bound = 0.001;
  constexpr double rm_upper_bound = 2.0;

  double observed_acceptance_probability = 1.0;
  if (observed_log_acceptance_probability < 0.0) {
    observed_acceptance_probability = MY_EXP (observed_log_acceptance_probability);
  }

  double updated_sd = current_sd +
    (observed_acceptance_probability - target_acceptance) * rm_weight;

  if (std::isnan (updated_sd)) {
    updated_sd = 1.0;
  }

  return std::clamp (updated_sd, rm_lower_bound, rm_upper_bound);
}

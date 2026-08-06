#ifndef BGMS_VARIABLE_HELPERS_H
#define BGMS_VARIABLE_HELPERS_H

#include <RcppArmadillo.h>
#include "math/explog_macros.h"


/**
 * Holds both log-normalizer and probabilities from joint computation.
 *
 * Used by logp_and_gradient to avoid duplicate probability/denominator
 * calculations.
 */
struct LogZAndProbs {
  /// Log-normalizer for each person.
  arma::vec log_Z;
  /// Probability matrix (num_persons x num_cats+1).
  arma::mat probs;
};


/**
 * Reusable scratch buffers for compute_logZ_and_probs_*_into.
 *
 * Owned by the caller (typically as a member of a model class, one per
 * MCMC chain). All buffers are resized via set_size() on each call,
 * which is a no-op once they've grown to the maximum needed size, so
 * after a few iterations no further heap allocations occur.
 */
struct LogZScratch {
  // Working buffers used by both ordinal and Blume-Capel paths.
  arma::vec eR, eB, pow_buf, den, term;
  arma::vec eM;          // ordinal: exp(main_param)
  arma::vec b_clamp;     // ordinal safe block: clamped bound
  arma::vec ex;          // ordinal/BC safe block: per-category exponent
  // Blume-Capel only
  arma::vec cat_vec, centered, theta, exp_theta;
  arma::vec pow_bound_low, pow_bound_high, pow_bound;
};


/**
 * Compute a numerically stable sum of the form:
 *
 *   denom = exp(-bound) + sum_{cat=0}^{K-1} exp(main_effect_param(cat)
 *                 + (cat + 1) * residual_score - bound)
 *
 * but evaluated efficiently using precomputed exponentials:
 *
 *   exp_r = exp(residual_score)
 *   exp_m = exp(main_effect_param)
 *   denom = exp(-bound) * ( 1 + sum_c exp_m[c] * exp_r^(c+1) )
 *
 * If non-finite values arise (overflow, underflow, NaN), a safe fallback
 * recomputes the naive version using direct exponentials.
 */
arma::vec compute_denom_ordinal(
    const arma::vec& residual,
    const arma::vec& main_eff,
    const arma::vec& bound
);

/**
 * Compute denom = Sigma_c exp( theta(c) + (c-ref)*r - b ), with
 *    theta(c) = lin_eff*(c-ref) + quad_eff*(c-ref)^2
 *    b    = max_c( theta(c) + (c-ref)*r )   (vectorized)
 *
 * Two modes:
 *
 * FAST (preexp + power-chain):
 *    denom = Sigma_c exp_theta[c] * exp(-b) * exp(r)^(c-ref)
 * Used only when all exponent terms are safe:
 *    |b| <= EXP_BOUND,
 *    underflow_bound >= -EXP_BOUND,
 *    num_cats*r - b <= EXP_BOUND.
 * This guarantees the recursive pow-chain stays finite.
 *
 * SAFE (direct evaluation):
 *    denom = Sigma_c exp(theta(c) + (c-ref)*r - b)
 * Used whenever any FAST-condition fails. Slower but always stable.
 *
 * FAST gives identical results when safe, otherwise SAFE is used.
 */
arma::vec compute_denom_blume_capel(
    const arma::vec& residual,
    const double lin_eff,
    const double quad_eff,
    const int ref,
    const int num_cats,
    arma::vec& b
);

/**
 * Compute category probabilities in a numerically stable manner.
 *
 * Uses pre-exp or bounded formulations depending on the magnitude of `bound`.
 *  - If |bound| <= 709 (EXP_BOUND): uses cheaper direct pre-exp computation
 *  - Else: clips bound at zero and applies stabilized scaling
 *
 * Clipping is required for bound < -709; the bounds improve stability
 * when large.
 *
 * Returns:
 *   probs: num_persons × (num_cats + 1) matrix of probabilities (row-normalized)
 */
arma::mat compute_probs_ordinal(
    const arma::vec& main_param,
    const arma::vec& residual_score,
    const arma::vec& bound,
    int num_cats
);

/**
 * Blume-Capel probabilities, numerically stable via FAST/SAFE split.
 *
 * Model:
 *   theta(c) = lin_eff * (c - ref) + quad_eff * (c - ref)^2,  c = 0..num_cats
 *   exps_i(c) = theta(c) + (c-ref) * r_i
 *   b_i       = max_c exps_i(c)
 *
 * Probabilities:
 *   p_i(c) proportional to exp( exps_i(c) - b_i )
 *
 * FAST (preexp + power-chain, same bounds as compute_denom_blume_capel):
 *   used when |b_i| <= EXP_BOUND and pow_bound_i = num_cats * r_i - b_i <= EXP_BOUND
 *
 * SAFE (direct):
 *   used otherwise: direct exp(theta(c) + (c-ref) * r_i - b_i)
 *
 * Under these conditions, denom is finite and > 0, so no one-hot fallback.
 */
arma::mat compute_probs_blume_capel(
    const arma::vec& residual,
    const double lin_eff,
    const double quad_eff,
    const int ref,
    const int num_cats,
    arma::vec& b
);

/**
 * Joint computation of log-normalizer and probabilities for ordinal variables.
 *
 * Computes both in a single pass, writing into caller-provided `out` and
 * reusing caller-provided `scratch`, which eliminates per-call heap
 * allocations after the first warm-up call.
 */
void compute_logZ_and_probs_ordinal_into(
    const arma::vec& main_param,
    const arma::vec& residual_score,
    const arma::vec& bound,
    int num_cats,
    LogZAndProbs& out,
    LogZScratch& scratch
);

/**
 * Joint computation of log-normalizer and probabilities for Blume-Capel variables.
 *
 * Computes both in a single pass, writing into caller-provided `out` and
 * reusing caller-provided `scratch`.
 */
void compute_logZ_and_probs_blume_capel_into(
    const arma::vec& residual,
    const double lin_eff,
    const double quad_eff,
    const int ref,
    const int num_cats,
    arma::vec& b,
    LogZAndProbs& out,
    LogZScratch& scratch
);

#endif // BGMS_VARIABLE_HELPERS_H
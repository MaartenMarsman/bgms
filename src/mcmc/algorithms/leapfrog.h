#pragma once

#include <RcppArmadillo.h>
#include <cstring>
#include <functional>
#include <utility>

// ---------------------------------------------------------------------------
// Memoizer — single-entry cache for joint log-posterior + gradient evaluations
// ---------------------------------------------------------------------------

/**
 * Memoizer - Single-entry cache for joint log-posterior and gradient evaluations.
 *
 * In NUTS, the typical access pattern within a leapfrog step is:
 *   1. cached_grad(theta) — compute gradient (and cache logp as side-effect)
 *   2. cached_log_post(theta) — retrieve the already-cached logp
 *
 * A single-entry cache is optimal here because each leapfrog step produces
 * a new unique theta: hash-map lookups would almost never hit, and hashing
 * an arma::vec element-by-element is expensive.
 *
 * The joint evaluation function computes both logp and gradient together
 * (since models often share most of the computation between the two).
 */
class Memoizer {
public:
  using JointFn = std::function<std::pair<double, arma::vec>(const arma::vec&)>;

  JointFn joint_fn;

  // Single-entry cache
  arma::vec cached_theta;
  double    cached_logp_val;
  arma::vec cached_grad_val;
  bool      has_cache = false;

  /**
   * Construct from separate log_post and grad functions.
   * Calls them independently (backward-compatible).
   */
  Memoizer(
    const std::function<double(const arma::vec&)>& lp,
    const std::function<arma::vec(const arma::vec&)>& gr
  ) : joint_fn([lp, gr](const arma::vec& theta) -> std::pair<double, arma::vec> {
        arma::vec g = gr(theta);
        double v = lp(theta);
        return {v, std::move(g)};
      }) {}

  /**
   * Construct from a joint function that computes both at once.
   */
  explicit Memoizer(JointFn jf) : joint_fn(std::move(jf)) {}

  double cached_log_post(const arma::vec& theta) {
    ensure_cached(theta);
    return cached_logp_val;
  }

  const arma::vec& cached_grad(const arma::vec& theta) {
    ensure_cached(theta);
    return cached_grad_val;
  }

  /**
   * Invalidate the cached gradient (e.g. after in-place position projection).
   *
   * After RATTLE projection modifies x, the cached gradient is for the
   * pre-projection position and must be recomputed.
   */
  void invalidate() { has_cache = false; }

  /**
   * Copy another Memoizer's cached entry (position, logp, gradient).
   *
   * Lets a secondary Memoizer start from an evaluation the primary one
   * already holds instead of recomputing it.
   */
  void seed_from(const Memoizer& other) {
    cached_theta = other.cached_theta;
    cached_logp_val = other.cached_logp_val;
    cached_grad_val = other.cached_grad_val;
    has_cache = other.has_cache;
  }

private:
  void ensure_cached(const arma::vec& theta) {
    if (has_cache &&
        theta.n_elem == cached_theta.n_elem &&
        std::memcmp(theta.memptr(), cached_theta.memptr(),
                    theta.n_elem * sizeof(double)) == 0) {
      return;
    }
    auto [lp, gr] = joint_fn(theta);
    cached_theta = theta;
    cached_logp_val = lp;
    cached_grad_val = std::move(gr);
    has_cache = true;
  }
};

// ---------------------------------------------------------------------------
// Leapfrog integrator — two variants (memoized for NUTS, joint for HMC)
// ---------------------------------------------------------------------------

/**
 * Performs a single leapfrog step with memoized gradient evaluation.
 * Used by NUTS tree-building.
 *
 * @param theta          Current position (parameter vector)
 * @param r              Current momentum vector
 * @param eps            Step size for integration
 * @param memo           Memoizer caching gradient evaluations
 * @param inv_mass_diag  Diagonal of the inverse mass matrix
 * @return Pair of (updated position, updated momentum)
 */
std::pair<arma::vec, arma::vec> leapfrog_memo(
    const arma::vec& theta,
    const arma::vec& r,
    double eps,
    Memoizer& memo,
    const arma::vec& inv_mass_diag
);

/**
 * Projection callback for SHAKE position constraint.
 * Projects position onto the constraint manifold c(q) = 0.
 */
using ProjectPositionFn = std::function<void(arma::vec& x)>;

/**
 * Projection callback for RATTLE velocity constraint.
 * Projects momentum onto the cotangent space: J M^{-1} r = 0.
 */
using ProjectMomentumFn = std::function<void(arma::vec& r, const arma::vec& x)>;

/**
 * Legacy projection callback combining position + momentum projection.
 * Retained for the test interface (ggm_gradient_interface.cpp).
 */
using ProjectFn = std::function<void(arma::vec& x, arma::vec& r)>;

/**
 * Performs a single constrained leapfrog step (RATTLE scheme).
 *
 * Structure follows Mici / Reich (1996):
 *   1. Half-step momentum
 *   2. Project momentum onto cotangent space
 *   3. Full-step position
 *   4. SHAKE: project position onto constraint manifold
 *   5. Momentum correction for constraint forces
 *   6. Second half-step momentum
 *   7. Project momentum onto cotangent space
 *
 * Position and momentum projections are separate, eliminating the
 * wasted PCG solve in the old bundled-projection implementation.
 *
 * @param theta            Current position (parameter vector)
 * @param r                Current momentum vector
 * @param eps              Step size for integration
 * @param memo             Memoizer caching gradient evaluations
 * @param inv_mass_diag    Diagonal of the inverse mass matrix
 * @param project_position SHAKE position projection callback
 * @param project_momentum RATTLE momentum projection callback
 * @return Pair of (updated position, updated momentum)
 */
std::pair<arma::vec, arma::vec> leapfrog_constrained(
    const arma::vec& theta,
    const arma::vec& r,
    double eps,
    Memoizer& memo,
    const arma::vec& inv_mass_diag,
    const ProjectPositionFn& project_position,
    const ProjectMomentumFn& project_momentum
);


// ---------------------------------------------------------------------------
// Constrained leapfrog with runtime reversibility check
// ---------------------------------------------------------------------------

/**
 * Result of a constrained leapfrog step with reversibility information.
 */
struct ConstrainedLeapfrogResult {
  arma::vec theta;    ///< Updated position
  arma::vec r;        ///< Updated momentum
  bool reversible;    ///< Whether the forward-backward check passed
};


/**
 * Constrained leapfrog step with a runtime reversibility check.
 *
 * Performs a forward constrained leapfrog step, then a backward step
 * (negate momentum, step forward, negate again).  If the round-trip
 * position differs from the original by more than factor * eps^2
 * in max-norm, the step is flagged as non-reversible.
 *
 * This is the runtime analogue of Mici's ConstrainedLeapfrogIntegrator
 * reverse check (Zappa et al., 2018; Lelièvre et al., 2019), adapted
 * to use eps^2-scaled tolerance matching the O(eps^2) column-coupling
 * error of the direct SHAKE solver.
 *
 * @param theta                Current position (parameter vector)
 * @param r                    Current momentum vector
 * @param eps                  Step size for integration
 * @param memo                 Memoizer caching gradient evaluations
 * @param inv_mass_diag        Diagonal of the inverse mass matrix
 * @param project_position     SHAKE position projection callback
 * @param project_momentum     RATTLE momentum projection callback
 * @param reverse_check_tol Factor for eps^2-scaled tolerance
 * @return ConstrainedLeapfrogResult with position, momentum, and reversibility flag
 */
ConstrainedLeapfrogResult leapfrog_constrained_checked(
    const arma::vec& theta,
    const arma::vec& r,
    double eps,
    Memoizer& memo,
    const arma::vec& inv_mass_diag,
    const ProjectPositionFn& project_position,
    const ProjectMomentumFn& project_momentum,
    double reverse_check_tol
);


/**
 * LeapfrogJointResult - Return type for multi-step leapfrog integration.
 *
 * Contains the final position, momentum, log-posterior, and gradient.
 */
struct LeapfrogJointResult {
  arma::vec theta;      ///< Final position
  arma::vec r;          ///< Final momentum
  double log_post;      ///< Log-posterior at final position
  arma::vec grad;       ///< Gradient at final position
};


/**
 * Multi-step leapfrog integration using a joint log_post+gradient function.
 * Used by HMC.
 *
 * Uses grad-only at intermediate steps; joint function at the final
 * position for both log_post and gradient. Accepts optional pre-computed
 * initial gradient to avoid recomputation.
 *
 * @param theta          Initial position
 * @param r              Initial momentum
 * @param eps            Step size
 * @param grad           Gradient-only function (for intermediate steps)
 * @param joint          Joint function returning (log_post, grad) pair
 * @param num_leapfrogs  Number of leapfrog steps
 * @param inv_mass_diag  Diagonal inverse mass matrix
 * @param init_grad      Optional pre-computed gradient at theta (nullptr to compute)
 * @return LeapfrogJointResult with final position, momentum, log_post, and gradient
 */
LeapfrogJointResult leapfrog(
    const arma::vec& theta,
    const arma::vec& r,
    double eps,
    const std::function<arma::vec(const arma::vec&)>& grad,
    const std::function<std::pair<double, arma::vec>(const arma::vec&)>& joint,
    int num_leapfrogs,
    const arma::vec& inv_mass_diag,
    const arma::vec* init_grad = nullptr
);

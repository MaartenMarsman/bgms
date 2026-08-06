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
 * LeapfrogJointResult - Return type for the joint leapfrog step.
 *
 * Contains the final position, momentum, and log-posterior.
 */
struct LeapfrogJointResult {
  arma::vec theta;      ///< Final position
  arma::vec r;          ///< Final momentum
  double log_post;      ///< Log-posterior at final position
};


/**
 * Single leapfrog step using a joint log_post+gradient function.
 * Used by the step-size heuristic in hamiltonian_utils.cpp.
 *
 * Evaluates the joint function at the new position for both log_post and
 * gradient. The gradient at the initial position is supplied by the caller,
 * which already holds it from its own joint evaluation.
 *
 * @param theta          Initial position
 * @param r              Initial momentum
 * @param eps            Step size
 * @param joint          Joint function returning (log_post, grad) pair
 * @param inv_mass_diag  Diagonal inverse mass matrix
 * @param init_grad      Pre-computed gradient at theta
 * @return LeapfrogJointResult with final position, momentum, and log_post
 */
LeapfrogJointResult leapfrog(
    const arma::vec& theta,
    const arma::vec& r,
    double eps,
    const std::function<std::pair<double, arma::vec>(const arma::vec&)>& joint,
    const arma::vec& inv_mass_diag,
    const arma::vec& init_grad
);

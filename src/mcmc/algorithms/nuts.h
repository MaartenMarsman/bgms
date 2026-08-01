#pragma once

#include <RcppArmadillo.h>
#include <functional>
#include <utility>
#include "mcmc/execution/step_result.h"
#include "mcmc/algorithms/leapfrog.h"
struct SafeRNG;


/**
 * BuildTreeResult - Return values of the recursive NUTS tree expansion
 *
 * Each call to build_tree expands the sampling path and may return new
 * candidate samples or indicate when the trajectory should terminate.
 */
struct BuildTreeResult {
  arma::vec theta_min;     ///< Leftmost position in the trajectory
  arma::vec r_min;         ///< Corresponding momentum at theta_min
  arma::vec theta_plus;    ///< Rightmost position in the trajectory
  arma::vec r_plus;        ///< Corresponding momentum at theta_plus
  arma::vec theta_prime;   ///< Current proposed sample (to possibly accept)
  arma::vec r_prime;       ///< Momentum at theta_prime (for energy diagnostics)
  double logp_prime;       ///< Log-posterior at theta_prime (avoids a re-eval in nuts_step)
  arma::vec rho;           ///< Sum of momenta along the subtree (for U-turn criterion)
  arma::vec p_sharp_beg;   ///< Sharp momentum (M^{-1} p) at subtree beginning
  arma::vec p_sharp_end;   ///< Sharp momentum (M^{-1} p) at subtree end
  arma::vec p_beg;         ///< Momentum at subtree beginning
  arma::vec p_end;         ///< Momentum at subtree end
  double log_sum_weight;   ///< log sum_i exp(H0 - h_i) across this subtree
  int s_prime;             ///< Stop flag (1 = continue, 0 = stop expansion)
  double alpha;            ///< Sum of min(1, exp(H0 - h)) across leapfrog steps
  int n_leapfrog;          ///< Number of leapfrog steps contributing to alpha
  bool divergent;          ///< Whether this subtree diverged
};



/**
 * Executes the No-U-Turn Sampler algorithm (NUTS)
 *
 * Takes a joint log_post+gradient function for efficient memoization.
 * The joint function computes both values together, which is more efficient
 * when they share common computations (e.g., normalization constants).
 *
 * @param init_theta        Initial position (parameter vector)
 * @param step_size         Step size for leapfrog integration
 * @param joint             Function returning (log_post, gradient) pair
 * @param inv_mass_diag     Diagonal of the inverse mass matrix
 * @param rng               Thread-safe random number generator
 * @param max_depth         Maximum tree depth (default = 10)
 * @return StepResult with position, acceptance probability, and NUTS diagnostics
 */
StepResult nuts_step(
    const arma::vec& init_theta,
    double step_size,
    const std::function<std::pair<double, arma::vec>(const arma::vec&)>& joint,
    const arma::vec& inv_mass_diag,
    SafeRNG& rng,
    int max_depth = 10
);

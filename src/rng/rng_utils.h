/**
 * @file rng_utils.h
 * @brief Thread-safe random number generation for MCMC sampling.
 *
 * Provides SafeRNG, a seedable wrapper around dqrng's xoshiro256++ engine,
 * and a set of distribution helpers that mirror R's `runif`, `rnorm`, etc.
 * Each model instance owns its own SafeRNG so parallel chains do not share
 * state. Boost.Random distributions are used instead of `<random>` to
 * guarantee identical output across compilers and platforms.
 *
 * Two groups of helpers:
 *   - **Scalar** (`runif`, `rnorm`, `rbern`, `rbeta`, `rexp`) — return a
 *     single draw.
 *   - **Armadillo** (`arma_rnorm_vec`, `arma_rnorm_mat`, `arma_runif_vec`,
 *     `arma_runif_mat`, `arma_randperm`) — fill vectors/matrices element-
 *     wise through the same scalar primitives.
 */

// [[Rcpp::depends(BH)]]
#pragma once

// the order of these two is mandatory, RcppArmadillo must come before dqrng
#include <RcppArmadillo.h>
// Only include xoshiro.h from dqrng - avoid dqrng_distribution.h which pulls
// in dqrng_generator.h -> convert_seed.h that has GCC 14 compatibility issues
#include <xoshiro.h>
#include <random>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/beta_distribution.hpp>
#include <boost/random/gamma_distribution.hpp>

// [[Rcpp::depends(dqrng, BH)]]

/**
 * Thread-safe random number generator wrapper.
 *
 * Wraps dqrng::xoshiro256plusplus with the `result_type` / `min()` / `max()`
 * interface required by `std::` and `boost::random` distributions.
 * Each MCMC chain owns one instance, seeded deterministically from the
 * user-supplied seed.
 */
struct SafeRNG {
  dqrng::xoshiro256plusplus eng;  ///< Underlying xoshiro256++ engine.

  // Default constructor // TODO: perhaps delete this to require a seed
  SafeRNG() : eng(1) {}

  /**
   * Construct with a deterministic seed.
   * @param seed  Integer seed forwarded to xoshiro256++
   */
  SafeRNG(const int seed) : eng(seed) {}

  /// Type alias required by distribution adaptors.
  using result_type = uint64_t;
  /// @return Minimum value the engine can produce.
  static constexpr result_type min() { return dqrng::xoshiro256plusplus::min(); }
  /// @return Maximum value the engine can produce.
  static constexpr result_type max() { return dqrng::xoshiro256plusplus::max(); }
  /// Advance the engine and return the next pseudorandom integer.
  result_type operator()() { return eng(); }
};


// ============================================================
// Scalar RNG helpers
// ============================================================

/**
 * Draw from Uniform(0, 1).
 * @param rng  Random number generator
 * @return A uniform random variate in [0, 1)
 */
inline double runif(SafeRNG& rng) {
  return boost::random::uniform_real_distribution<double>(0.0, 1.0)(rng.eng);
}

/**
 * Draw from Normal(mu, sigma).
 * @param rng    Random number generator
 * @param mu     Mean (default 0)
 * @param sigma  Standard deviation (default 1)
 * @return A normally distributed variate
 */
inline double rnorm(SafeRNG& rng, double mu = 0.0, double sigma = 1.0) {
  return boost::random::normal_distribution<double>(mu, sigma)(rng.eng);
}

/**
 * Draw from Bernoulli(p).
 * @param rng  Random number generator
 * @param p    Success probability
 * @return 1 with probability p, 0 otherwise
 */
inline int rbern(SafeRNG& rng, double p) {
  return (runif(rng) < p) ? 1 : 0;
}

/**
 * Draw from Beta(a, b).
 * @param rng  Random number generator
 * @param a    First shape parameter (alpha)
 * @param b    Second shape parameter (beta)
 * @return A beta-distributed variate in (0, 1)
 */
inline double rbeta(SafeRNG& rng, double a, double b) {
  return boost::random::beta_distribution<double>(a, b)(rng.eng);
}

/**
 * Draw from Exponential(lambda).
 * @param rng     Random number generator
 * @param lambda  Rate parameter (1 / mean)
 * @return An exponentially distributed variate
 */
inline double rexp(SafeRNG& rng, double lambda) {
  return boost::random::exponential_distribution<double>(lambda)(rng.eng);
}

/**
 * Draw from Gamma(shape, rate).
 * @param rng    Random number generator
 * @param shape  Shape parameter (alpha)
 * @param rate   Rate parameter (beta); mean = shape / rate
 * @return A gamma-distributed variate. boost parameterizes by scale = 1/rate.
 */
inline double rgamma(SafeRNG& rng, double shape, double rate) {
  return boost::random::gamma_distribution<double>(shape, 1.0 / rate)(rng.eng);
}

// ============================================================
// Armadillo RNG helpers
// ============================================================

/**
 * Fill a vector with Normal(mu, sigma) draws.
 * @param rng    Random number generator
 * @param n      Number of elements
 * @param mu     Mean (default 0)
 * @param sigma  Standard deviation (default 1)
 * @return Column vector of n independent normal variates
 */
inline arma::vec arma_rnorm_vec(SafeRNG& rng,
                                arma::uword n,
                                double mu = 0.0, double sigma = 1.0) {
  arma::vec out(n);
  for (arma::uword i = 0; i < n; ++i)
    out[i] = rnorm(rng, mu, sigma);
  return out;
}

/**
 * Fill a matrix with Normal(mu, sigma) draws (column-major order).
 * @param rng    Random number generator
 * @param nrow   Number of rows
 * @param ncol   Number of columns
 * @param mu     Mean (default 0)
 * @param sigma  Standard deviation (default 1)
 * @return Matrix of independent normal variates
 */
inline arma::mat arma_rnorm_mat(SafeRNG& rng,
                                arma::uword nrow, arma::uword ncol,
                                double mu = 0.0, double sigma = 1.0) {
  arma::mat out(nrow, ncol);
  for (arma::uword j = 0; j < ncol; ++j) {
    for (arma::uword i = 0; i < nrow; ++i) {
      out(i, j) = rnorm(rng, mu, sigma);
    }
  }
  return out;
}

/**
 * Fill a vector with Uniform(0, 1) draws.
 * @param rng  Random number generator
 * @param n    Number of elements
 * @return Column vector of n independent uniform variates
 */
inline arma::vec arma_runif_vec(SafeRNG& rng, arma::uword n) {
  arma::vec out(n);
  for (arma::uword i = 0; i < n; ++i)
    out[i] = runif(rng);
  return out;
}

/**
 * Fill a matrix with Uniform(0, 1) draws (column-major order).
 * @param rng   Random number generator
 * @param nrow  Number of rows
 * @param ncol  Number of columns
 * @return Matrix of independent uniform variates
 */
inline arma::mat arma_runif_mat(SafeRNG& rng,
                                arma::uword nrow, arma::uword ncol) {
  arma::mat out(nrow, ncol);
  for (arma::uword j = 0; j < ncol; ++j) {
    for (arma::uword i = 0; i < nrow; ++i) {
      out(i, j) = runif(rng);
    }
  }
  return out;
}

/**
 * Random permutation of 0, 1, ..., n-1 (like arma::randperm).
 * @param rng  Random number generator
 * @param n    Number of elements
 * @return Permuted index vector
 */
inline arma::uvec arma_randperm(SafeRNG& rng, arma::uword n) {
  arma::uvec out(n);
  std::iota(out.begin(), out.end(), 0);
  std::shuffle(out.begin(), out.end(), rng.eng);
  return out;
}


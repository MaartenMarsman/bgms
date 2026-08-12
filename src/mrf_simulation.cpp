// [[Rcpp::depends(RcppParallel, RcppArmadillo, dqrng, BH)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <tbb/global_control.h>
#include "math/explog_macros.h"
#include "rng/rng_utils.h"
#include "utils/progress_manager.h"
#include <vector>
#include <string>

using namespace Rcpp;
using namespace RcppParallel;


// ============================================================================
//   MRF Simulation Core Functions (Thread-Safe)
// ============================================================================

// Function: simulate_mrf
//
// Simulates observations from a Markov Random Field using Gibbs sampling.
// Supports both ordinal and Blume-Capel variable types.
//
// Inputs:
//  - num_states: Number of observations to simulate.
//  - num_variables: Number of variables in the MRF.
//  - num_categories: Number of categories per variable (on top of baseline 0).
//  - pairwise: Symmetric pairwise interaction matrix (diagonal ignored).
//  - main: Main effect parameters (variables x max_categories).
//          For ordinal: threshold parameters for categories 1..K.
//          For Blume-Capel: column 0 = linear (alpha), column 1 = quadratic (beta).
//  - variable_type: Type of each variable ("ordinal" or "blume-capel").
//  - baseline_category: Baseline category for Blume-Capel variables (0 for ordinal).
//  - iter: Number of Gibbs sampling iterations.
//  - rng: Thread-safe random number generator.
//
// Returns:
//  - Integer matrix of simulated observations (num_states x num_variables).
//
// Notes:
//  - Diagonal of pairwise matrix is explicitly ignored (set to zero internally).
//  - For ordinal variables, baseline_category should be 0.
arma::imat simulate_mrf(
    int num_states,
    int num_variables,
    const arma::ivec& num_categories,
    const arma::mat& pairwise,
    const arma::mat& main,
    const std::vector<std::string>& variable_type,
    const arma::ivec& baseline_category,
    int iter,
    SafeRNG& rng) {

  arma::imat observations(num_states, num_variables);
  int max_num_categories = arma::max(num_categories);
  arma::vec probabilities(max_num_categories + 1);
  // Category exponents of the full conditional, held so the largest can be
  // subtracted before exponentiating (as compute_probs_ordinal() and
  // compute_probs_blume_capel() already do on the prediction side). Without it
  // a large rest score overflows exp() to Inf and the draw degenerates.
  arma::vec exponents(max_num_categories + 1);
  double exponent = 0.0;
  double rest_score = 0.0;
  double cumsum = 0.0;
  double u = 0.0;
  int score = 0;

  // Copy pairwise and zero diagonal to prevent accidental self-interactions
  arma::mat pairwise_safe = pairwise;
  pairwise_safe.diag().zeros();

  // Random (uniform) starting values
  for(int variable = 0; variable < num_variables; variable++) {
    for(int person = 0; person < num_states; person++) {
      cumsum = 1.0;
      probabilities[0] = 1.0;
      for(int category = 0; category < num_categories[variable]; category++) {
        cumsum += 1;
        probabilities[category + 1] = cumsum;
      }

      u = cumsum * runif(rng);

      score = 0;
      while (score < num_categories[variable] && u > probabilities[score]) {
        score++;
      }
      observations(person, variable) = score;
    }
  }

  // Gibbs sampling iterations
  for(int iteration = 0; iteration < iter; iteration++) {
    for(int variable = 0; variable < num_variables; variable++) {
      for(int person = 0; person < num_states; person++) {
        // Compute rest score using centered parameterization
        // For ordinal variables with baseline_category=0, this is equivalent to obs * pairwise
        rest_score = 0.0;
        for(int vertex = 0; vertex < num_variables; vertex++) {
          int obs = observations(person, vertex);
          int ref = baseline_category[vertex];
          rest_score += 2.0 * (obs - ref) * pairwise_safe(vertex, variable);
        }

        // Both branches collect the category exponents first, subtract their
        // maximum, and only then exponentiate: the cumulative weights are
        // sampled from up to a common factor, so the shift leaves the
        // distribution untouched while keeping every term in range.
        double max_exponent = 0.0;
        if(variable_type[variable] == "blume-capel") {
          int ref = baseline_category[variable];
          for(int category = 0; category <= num_categories[variable]; category++) {
            const int s = category - ref;
            // Linear term
            exponent = main(variable, 0) * s;
            // Quadratic term
            exponent += main(variable, 1) * s * s;
            // Pairwise effects
            exponent += rest_score * s;
            exponents[category] = exponent;
            if(category == 0 || exponent > max_exponent) max_exponent = exponent;
          }
          cumsum = 0.0;
          for(int category = 0; category <= num_categories[variable]; category++) {
            cumsum += MY_EXP(exponents[category] - max_exponent);
            probabilities[category] = cumsum;
          }
        } else {
          // Ordinal: baseline category 0 carries exponent 0
          exponents[0] = 0.0;
          max_exponent = 0.0;
          for(int category = 0; category < num_categories[variable]; category++) {
            exponent = main(variable, category);
            exponent += (category + 1) * rest_score;
            exponents[category + 1] = exponent;
            if(exponent > max_exponent) max_exponent = exponent;
          }
          cumsum = 0.0;
          for(int category = 0; category <= num_categories[variable]; category++) {
            cumsum += MY_EXP(exponents[category] - max_exponent);
            probabilities[category] = cumsum;
          }
        }

        u = cumsum * runif(rng);

        // Sample category with bounds protection
        score = 0;
        int max_score = num_categories[variable];
        while (score < max_score && u > probabilities[score]) {
          score++;
        }
        observations(person, variable) = score;
      }
    }
  }

  return observations;
}


// ============================================================================
//   R Interface for mrfSampler()
// ============================================================================

// [[Rcpp::export]]
IntegerMatrix sample_omrf_gibbs(int num_states,
                                int num_variables,
                                IntegerVector num_categories,
                                NumericMatrix pairwise,
                                NumericMatrix main,
                                int iter,
                                int seed) {

  SafeRNG rng(seed);

  // Convert inputs to arma types
  arma::ivec num_categories_arma = Rcpp::as<arma::ivec>(num_categories);
  arma::mat pairwise_arma = Rcpp::as<arma::mat>(pairwise);
  arma::mat main_arma = Rcpp::as<arma::mat>(main);

  // Create ordinal defaults: all variables are "ordinal" with baseline_category = 0
  std::vector<std::string> variable_type(num_variables, "ordinal");
  arma::ivec baseline_category_arma(num_variables, arma::fill::zeros);

  // Simulate observations
  arma::imat result = simulate_mrf(
    num_states,
    num_variables,
    num_categories_arma,
    pairwise_arma,
    main_arma,
    variable_type,
    baseline_category_arma,
    iter,
    rng
  );

  // Check for user interrupt periodically (only in non-parallel context)
  Rcpp::checkUserInterrupt();

  return Rcpp::wrap(result);
}

// [[Rcpp::export]]
IntegerMatrix sample_bcomrf_gibbs(int num_states,
                                  int num_variables,
                                  IntegerVector num_categories,
                                  NumericMatrix pairwise,
                                  NumericMatrix main,
                                  StringVector variable_type_r,
                                  IntegerVector baseline_category,
                                  int iter,
                                  int seed) {

  SafeRNG rng(seed);

  // Convert inputs to arma/std types
  arma::ivec num_categories_arma = Rcpp::as<arma::ivec>(num_categories);
  arma::mat pairwise_arma = Rcpp::as<arma::mat>(pairwise);
  arma::mat main_arma = Rcpp::as<arma::mat>(main);
  arma::ivec baseline_category_arma = Rcpp::as<arma::ivec>(baseline_category);

  std::vector<std::string> variable_type(num_variables);
  for (int i = 0; i < num_variables; i++) {
    variable_type[i] = Rcpp::as<std::string>(variable_type_r[i]);
    // Ordinal variables must use baseline_category = 0 (category 0 is the reference)
    if (variable_type[i] != "blume-capel") {
      baseline_category_arma[i] = 0;
    }
  }

  // Simulate observations
  arma::imat result = simulate_mrf(
    num_states,
    num_variables,
    num_categories_arma,
    pairwise_arma,
    main_arma,
    variable_type,
    baseline_category_arma,
    iter,
    rng
  );

  // Check for user interrupt (only in non-parallel context)
  Rcpp::checkUserInterrupt();

  return Rcpp::wrap(result);
}


// ============================================================================
//   GGM Simulation (Direct Multivariate Normal Sampling)
// ============================================================================

// Simulate observations from a Gaussian Graphical Model.
//
// Given a precision matrix Omega, draws num_states observations from
// N(means, Omega^{-1}) using the Cholesky factorization of the covariance.
//
// Algorithm:
//   1. Compute Sigma = Omega^{-1} via arma::inv_sympd.
//   2. Cholesky decompose: L = chol(Sigma, "lower") so Sigma = L L'.
//   3. Draw Z ~ N(0, I) of size (num_states x p).
//   4. Return X = ones * means' + Z * L'.
//
// Takes the number of observations, a p x p positive-definite precision
// matrix Omega, a p-vector of means (can be all zeros), and a thread-safe
// RNG. Returns a num_states x p matrix of simulated continuous observations.
arma::mat simulate_ggm(
    int num_states,
    const arma::mat& precision,
    const arma::vec& means,
    SafeRNG& rng) {

  int p = precision.n_cols;

  // Covariance = inverse of precision
  arma::mat sigma = arma::inv_sympd(precision);

  // Lower Cholesky: sigma = L * L'
  arma::mat L = arma::chol(sigma, "lower");

  // Draw standard normal matrix: num_states x p
  arma::mat Z = arma_rnorm_mat(rng, num_states, p);

  // Transform: X = Z * L' (each row is one observation)
  arma::mat X = Z * L.t();

  // Add means
  X.each_row() += means.t();

  return X;
}


// ============================================================================
//   R Interface for GGM Simulation (standalone simulate_ggm)
// ============================================================================

// R-callable wrapper for single GGM simulation. Takes the number of
// observations, a p x p precision matrix Omega, a p-vector of means, and a
// random seed. Returns a num_states x p numeric matrix.
// [[Rcpp::export]]
NumericMatrix sample_ggm_direct(int num_states,
                                NumericMatrix precision,
                                NumericVector means,
                                int seed) {

  SafeRNG rng(seed);

  arma::mat precision_arma = Rcpp::as<arma::mat>(precision);
  arma::vec means_arma = Rcpp::as<arma::vec>(means);

  arma::mat result = simulate_ggm(
    num_states,
    precision_arma,
    means_arma,
    rng
  );

  Rcpp::checkUserInterrupt();

  return Rcpp::wrap(result);
}


// ============================================================================
//   Parallel Simulation for simulate.bgms() with Posterior Draws (Ordinal)
// ============================================================================

// Structure to hold individual simulation results
struct SimulationResult {
  arma::imat observations;
  int draw_index;
  bool error;
  std::string error_msg;
};


// Worker class for parallel simulation across posterior draws
class SimulationWorker : public RcppParallel::Worker {
public:
  const arma::mat& pairwise_samples;
  const arma::mat& main_samples;
  const arma::ivec& draw_indices;
  const int num_states;
  const int num_variables;
  const arma::ivec& num_categories;
  const std::vector<std::string>& variable_type;
  const arma::ivec& baseline_category;
  const int iter;
  const arma::ivec& main_param_counts;
  const std::vector<SafeRNG>& draw_rngs;
  ProgressManager& pm;
  std::vector<SimulationResult>& results;

  SimulationWorker(
    const arma::mat& pairwise_samples,
    const arma::mat& main_samples,
    const arma::ivec& draw_indices,
    int num_states,
    int num_variables,
    const arma::ivec& num_categories,
    const std::vector<std::string>& variable_type,
    const arma::ivec& baseline_category,
    int iter,
    const arma::ivec& main_param_counts,
    const std::vector<SafeRNG>& draw_rngs,
    ProgressManager& pm,
    std::vector<SimulationResult>& results
  ) :
    pairwise_samples(pairwise_samples),
    main_samples(main_samples),
    draw_indices(draw_indices),
    num_states(num_states),
    num_variables(num_variables),
    num_categories(num_categories),
    variable_type(variable_type),
    baseline_category(baseline_category),
    iter(iter),
    main_param_counts(main_param_counts),
    draw_rngs(draw_rngs),
    pm(pm),
    results(results)
  {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      if (pm.shouldExit()) return;

      SimulationResult result;
      result.draw_index = draw_indices[i];
      result.error = false;

      try {
        SafeRNG rng = draw_rngs[i];

        // Reconstruct pairwise matrix from flat vector
        arma::mat pairwise(num_variables, num_variables, arma::fill::zeros);
        int idx = 0;
        for (int col = 0; col < num_variables; col++) {
          for (int row = col + 1; row < num_variables; row++) {
            pairwise(row, col) = pairwise_samples(draw_indices[i] - 1, idx);
            pairwise(col, row) = pairwise_samples(draw_indices[i] - 1, idx);
            idx++;
          }
        }

        // Reconstruct main effect matrix
        int max_main = arma::max(main_param_counts);
        arma::mat main(num_variables, max_main, arma::fill::zeros);
        idx = 0;
        for (int v = 0; v < num_variables; v++) {
          for (int t = 0; t < main_param_counts[v]; t++) {
            main(v, t) = main_samples(draw_indices[i] - 1, idx);
            idx++;
          }
        }

        result.observations = simulate_mrf(
          num_states,
          num_variables,
          num_categories,
          pairwise,
          main,
          variable_type,
          baseline_category,
          iter,
          rng
        );

      } catch (const std::exception& e) {
        result.error = true;
        result.error_msg = e.what();
      } catch (...) {
        result.error = true;
        result.error_msg = "Unknown error";
      }

      results[i] = result;
      pm.update(0);
    }
  }
};


// Run parallel simulations across posterior draws. Takes pairwise samples
// (ndraws x n_pairwise), main/threshold samples (ndraws x n_main), 1-based
// draw indices, observations per draw, variable metadata (count, categories
// per variable, "ordinal" or "blume-capel" type, baseline categories), Gibbs
// iterations per simulation, thread count, seed, and progress bar type
// (0 = none, 1 = total, 2 = per-chain). Returns a list of simulation results
// (each an integer matrix).
// [[Rcpp::export]]
Rcpp::List run_simulation_parallel(
    const arma::mat& pairwise_samples,
    const arma::mat& main_samples,
    const arma::ivec& draw_indices,
    int num_states,
    int num_variables,
    const arma::ivec& num_categories,
    const Rcpp::StringVector& variable_type_r,
    const arma::ivec& baseline_category,
    int iter,
    int nThreads,
    int seed,
    int progress_type) {

  int ndraws = draw_indices.n_elem;

  // Convert variable_type to std::vector<std::string>
  // and enforce baseline_category = 0 for ordinal variables
  std::vector<std::string> variable_type(num_variables);
  arma::ivec baseline_category_safe = baseline_category;
  for (int i = 0; i < num_variables; i++) {
    variable_type[i] = Rcpp::as<std::string>(variable_type_r[i]);
    // Ordinal variables must use baseline_category = 0 (category 0 is the reference)
    if (variable_type[i] != "blume-capel") {
      baseline_category_safe[i] = 0;
    }
  }

  // Compute number of main parameters per variable
  arma::ivec main_param_counts(num_variables);
  for (int v = 0; v < num_variables; v++) {
    if (variable_type[v] == "blume-capel") {
      main_param_counts[v] = 2;  // linear and quadratic
    } else {
      main_param_counts[v] = num_categories[v];  // K thresholds for K+1 response options
    }
  }

  // Prepare one independent RNG per draw
  std::vector<SafeRNG> draw_rngs(ndraws);
  for (int d = 0; d < ndraws; d++) {
    draw_rngs[d] = SafeRNG(seed + d);
  }

  // Prepare results storage
  std::vector<SimulationResult> results(ndraws);

  ProgressManager pm(1, ndraws, 0, 50, progress_type);

  SimulationWorker worker(
    pairwise_samples,
    main_samples,
    draw_indices,
    num_states,
    num_variables,
    num_categories,
    variable_type,
    baseline_category_safe,
    iter,
    main_param_counts,
    draw_rngs,
    pm,
    results
  );

  {
    tbb::global_control control(
      tbb::global_control::max_allowed_parallelism, nThreads);
    parallelFor(0, ndraws, worker);
  }
  pm.finish();

  // Every worker returns as soon as the interrupt flag is set, so the draws it
  // had not reached hold default-constructed matrices. Those are not a shorter
  // simulation, they are empty ones, so the interrupt is raised rather than
  // handed back as data.
  if (pm.shouldExit()) {
    Rcpp::stop("Simulation interrupted by user.");
  }

  // Convert results to R list
  Rcpp::List output(ndraws);
  for (int i = 0; i < ndraws; i++) {
    if (results[i].error) {
      Rcpp::stop("Error in simulation draw %d: %s",
                 results[i].draw_index, results[i].error_msg.c_str());
    }
    output[i] = Rcpp::wrap(results[i].observations);
  }

  return output;
}


// ============================================================================
//   Parallel GGM Simulation for simulate.bgms() with Posterior Draws
// ============================================================================

// Structure to hold GGM simulation results
struct GGMSimulationResult {
  arma::mat observations;
  int draw_index;
  bool error;
  std::string error_msg;
};


// Worker class for parallel GGM simulation across posterior draws
class GGMSimulationWorker : public RcppParallel::Worker {
public:
  const arma::mat& pairwise_samples;
  const arma::mat& main_samples;
  const arma::ivec& draw_indices;
  const int num_states;
  const int num_variables;
  const arma::vec& means;
  const std::vector<SafeRNG>& draw_rngs;
  ProgressManager& pm;
  std::vector<GGMSimulationResult>& results;

  GGMSimulationWorker(
    const arma::mat& pairwise_samples,
    const arma::mat& main_samples,
    const arma::ivec& draw_indices,
    int num_states,
    int num_variables,
    const arma::vec& means,
    const std::vector<SafeRNG>& draw_rngs,
    ProgressManager& pm,
    std::vector<GGMSimulationResult>& results
  ) :
    pairwise_samples(pairwise_samples),
    main_samples(main_samples),
    draw_indices(draw_indices),
    num_states(num_states),
    num_variables(num_variables),
    means(means),
    draw_rngs(draw_rngs),
    pm(pm),
    results(results)
  {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      if (pm.shouldExit()) return;

      GGMSimulationResult result;
      result.draw_index = draw_indices[i];
      result.error = false;

      try {
        SafeRNG rng = draw_rngs[i];

        // Reconstruct precision matrix from flat vectors
        arma::mat precision(num_variables, num_variables, arma::fill::zeros);
        int idx = 0;
        for (int col = 0; col < num_variables; col++) {
          for (int row = col + 1; row < num_variables; row++) {
            precision(row, col) = pairwise_samples(draw_indices[i] - 1, idx);
            precision(col, row) = pairwise_samples(draw_indices[i] - 1, idx);
            idx++;
          }
        }
        for (int v = 0; v < num_variables; v++) {
          precision(v, v) = main_samples(draw_indices[i] - 1, v);
        }

        result.observations = simulate_ggm(
          num_states,
          precision,
          means,
          rng
        );

      } catch (const std::exception& e) {
        result.error = true;
        result.error_msg = e.what();
      } catch (...) {
        result.error = true;
        result.error_msg = "Unknown error";
      }

      results[i] = result;
      pm.update(0);
    }
  }
};


// Run parallel GGM simulations across posterior draws. Takes off-diagonal
// precision samples (ndraws x p*(p-1)/2), diagonal precision samples
// (ndraws x p), 1-based draw indices, observations per draw, the number of
// variables, a p-vector of means, thread count, seed, and progress bar type
// (0 = none, 1 = total, 2 = per-chain). Returns a list of simulation results
// (each a numeric matrix n x p).
// [[Rcpp::export]]
Rcpp::List run_ggm_simulation_parallel(
    const arma::mat& pairwise_samples,
    const arma::mat& main_samples,
    const arma::ivec& draw_indices,
    int num_states,
    int num_variables,
    const arma::vec& means,
    int nThreads,
    int seed,
    int progress_type) {

  int ndraws = draw_indices.n_elem;

  // Prepare one independent RNG per draw
  std::vector<SafeRNG> draw_rngs(ndraws);
  for (int d = 0; d < ndraws; d++) {
    draw_rngs[d] = SafeRNG(seed + d);
  }

  std::vector<GGMSimulationResult> results(ndraws);
  ProgressManager pm(1, ndraws, 0, 50, progress_type);

  GGMSimulationWorker worker(
    pairwise_samples,
    main_samples,
    draw_indices,
    num_states,
    num_variables,
    means,
    draw_rngs,
    pm,
    results
  );

  {
    tbb::global_control control(
      tbb::global_control::max_allowed_parallelism, nThreads);
    parallelFor(0, ndraws, worker);
  }
  pm.finish();

  // See the OMRF runner: an interrupted worker leaves its remaining draws
  // empty, which must not reach the caller as a result.
  if (pm.shouldExit()) {
    Rcpp::stop("Simulation interrupted by user.");
  }

  Rcpp::List output(ndraws);
  for (int i = 0; i < ndraws; i++) {
    if (results[i].error) {
      Rcpp::stop("Error in GGM simulation draw %d: %s",
                 results[i].draw_index, results[i].error_msg.c_str());
    }
    output[i] = Rcpp::wrap(results[i].observations);
  }

  return output;
}


// ============================================================================
//   Mixed MRF Simulation (Block Gibbs: Discrete + Continuous)
// ============================================================================

// Simulate observations from a mixed MRF using block Gibbs sampling.
//
// Each iteration updates all discrete variables from their full conditional
// given (x_{-s}, y), then updates all continuous variables from the
// conditional Gaussian.
//
// Takes the number of observations, the interaction blocks (pairwise_disc:
// p x p with zero diagonal, pairwise_cross: p x q, pairwise_cont: q x q
// negative-definite), mux (p x max_cats threshold / Blume-Capel parameters),
// muy (q-vector of continuous means), per-variable metadata (num_categories:
// counts on top of baseline 0; variable_type: "ordinal" or "blume-capel";
// baseline_category: reference categories), the number of Gibbs burn-in
// iterations, and a thread-safe RNG. Writes the n x p integer discrete
// observations to x_out and the n x q numeric continuous observations to
// y_out.
void simulate_mixed_mrf(
    int num_states,
    const arma::mat& pairwise_disc,
    const arma::mat& pairwise_cross,
    const arma::mat& pairwise_cont,
    const arma::mat& mux,
    const arma::vec& muy,
    const arma::ivec& num_categories,
    const std::vector<std::string>& variable_type,
    const arma::ivec& baseline_category,
    int iter,
    SafeRNG& rng,
    arma::imat& x_out,
    arma::mat& y_out) {

  int p = pairwise_disc.n_rows;
  int q = pairwise_cont.n_rows;

  // Precompute precision decomposition: precision = -2 * pairwise_cont
  arma::mat precision_y = -2.0 * pairwise_cont;
  arma::mat Sigma_y = arma::inv_sympd(precision_y);
  arma::mat L_Sigma = arma::chol(Sigma_y, "lower");

  // Precompute cross * covariance for conditional mean
  arma::mat cross_sigma = pairwise_cross * Sigma_y;  // p x q

  // Copy discrete interactions with zeroed diagonal
  arma::mat disc_safe = pairwise_disc;
  disc_safe.diag().zeros();

  // Generate each observation independently
  for (int obs = 0; obs < num_states; obs++) {
    // Initialize discrete variables uniformly
    arma::ivec x_current(p);
    for (int s = 0; s < p; s++) {
      int max_cat = num_categories(s);
      x_current(s) = static_cast<int>(runif(rng) * (max_cat + 1));
      if (x_current(s) > max_cat) x_current(s) = max_cat;
    }

    // Initialize continuous from marginal N(muy, Sigma_y)
    arma::vec z = arma_rnorm_vec(rng, q);
    arma::vec y_current = muy + L_Sigma * z;

    // Gibbs iterations
    for (int it = 0; it < iter; it++) {

      // --- Update discrete variables from f(x_s | x_{-s}, y) ---
      for (int s = 0; s < p; s++) {
        int Cs = num_categories(s);

        // Rest score from discrete neighbours
        double rest_discrete = 0.0;
        for (int k = 0; k < p; k++) {
          if (k != s) {
            int obs_k = x_current(k);
            int ref_k = baseline_category(k);
            rest_discrete += 2.0 * (obs_k - ref_k) * disc_safe(k, s);
          }
        }

        // Rest score from continuous (factor of 2)
        double rest_continuous = 0.0;
        for (int j = 0; j < q; j++) {
          rest_continuous += 2.0 * pairwise_cross(s, j) * y_current(j);
        }

        double rest = rest_discrete + rest_continuous;

        if (variable_type[s] == "blume-capel") {
          int ref = baseline_category(s);
          double alpha = mux(s, 0);
          double beta  = mux(s, 1);

          // Compute log-probabilities for categories 0..Cs
          arma::vec log_probs(Cs + 1);
          double max_lp = -std::numeric_limits<double>::infinity();
          for (int c = 0; c <= Cs; c++) {
            int d = c - ref;
            log_probs(c) = alpha * d + beta * d * d + d * rest;
            if (log_probs(c) > max_lp) max_lp = log_probs(c);
          }

          // Stabilize and convert to cumulative probabilities
          double cumsum = 0.0;
          arma::vec cum_probs(Cs + 1);
          for (int c = 0; c <= Cs; c++) {
            cumsum += MY_EXP(log_probs(c) - max_lp);
            cum_probs(c) = cumsum;
          }

          // Sample
          double u = runif(rng) * cumsum;
          int sampled = 0;
          while (sampled < Cs && u > cum_probs(sampled)) sampled++;
          x_current(s) = sampled;

        } else {
          // Ordinal: category 0 is reference with log-prob = 0
          arma::vec log_probs(Cs + 1);
          log_probs(0) = 0.0;
          double max_lp = 0.0;
          for (int c = 1; c <= Cs; c++) {
            log_probs(c) = mux(s, c - 1) + c * rest;
            if (log_probs(c) > max_lp) max_lp = log_probs(c);
          }

          double cumsum = 0.0;
          arma::vec cum_probs(Cs + 1);
          for (int c = 0; c <= Cs; c++) {
            cumsum += MY_EXP(log_probs(c) - max_lp);
            cum_probs(c) = cumsum;
          }

          double u = runif(rng) * cumsum;
          int sampled = 0;
          while (sampled < Cs && u > cum_probs(sampled)) sampled++;
          x_current(s) = sampled;
        }
      }

      // --- Update continuous variables from y | x ---
      // y | x ~ N(muy + 2 * cross_sigma^T * x_centered, Sigma_y)
      // Compute centered discrete observations
      arma::vec x_centered(p);
      for (int s = 0; s < p; s++) {
        x_centered(s) = static_cast<double>(x_current(s) - baseline_category(s));
      }

      // Conditional mean: muy + 2 * Sigma_y * cross^T * x_centered
      //   = muy + 2 * cross_sigma^T * x_centered
      arma::vec cond_mean = muy + 2.0 * cross_sigma.t() * x_centered;

      // Sample y ~ N(cond_mean, Sigma_y)
      arma::vec z2 = arma_rnorm_vec(rng, q);
      y_current = cond_mean + L_Sigma * z2;
    }

    // Store final state
    x_out.row(obs) = x_current.t();
    y_out.row(obs) = y_current.t();
  }
}


// ============================================================================
//   R Interface for Mixed MRF Simulation (standalone)
// ============================================================================

// [[Rcpp::export]]
Rcpp::List sample_mixed_mrf_gibbs(
    int num_states,
    NumericMatrix pairwise_disc_r,
    NumericMatrix pairwise_cross_r,
    NumericMatrix pairwise_cont_r,
    NumericMatrix mux_r,
    NumericVector muy_r,
    IntegerVector num_categories_r,
    Rcpp::StringVector variable_type_r,
    IntegerVector baseline_category_r,
    int iter,
    int seed) {

  SafeRNG rng(seed);

  int p = pairwise_disc_r.nrow();
  int q = pairwise_cont_r.nrow();

  arma::mat pairwise_disc = Rcpp::as<arma::mat>(pairwise_disc_r);
  arma::mat pairwise_cross = Rcpp::as<arma::mat>(pairwise_cross_r);
  arma::mat pairwise_cont = Rcpp::as<arma::mat>(pairwise_cont_r);
  arma::mat mux = Rcpp::as<arma::mat>(mux_r);
  arma::vec muy = Rcpp::as<arma::vec>(muy_r);
  arma::ivec num_categories = Rcpp::as<arma::ivec>(num_categories_r);
  arma::ivec baseline_category = Rcpp::as<arma::ivec>(baseline_category_r);

  // Convert variable_type
  std::vector<std::string> variable_type(p);
  for (int i = 0; i < p; i++) {
    variable_type[i] = Rcpp::as<std::string>(variable_type_r[i]);
    if (variable_type[i] != "blume-capel") {
      baseline_category[i] = 0;
    }
  }

  arma::imat x_out(num_states, p);
  arma::mat y_out(num_states, q);

  simulate_mixed_mrf(
    num_states, pairwise_disc, pairwise_cross, pairwise_cont, mux, muy,
    num_categories, variable_type, baseline_category,
    iter, rng, x_out, y_out
  );

  Rcpp::checkUserInterrupt();

  return Rcpp::List::create(
    Rcpp::Named("x") = Rcpp::wrap(x_out),
    Rcpp::Named("y") = Rcpp::wrap(y_out)
  );
}


// ============================================================================
//   Parallel Mixed MRF Simulation for simulate.bgms() with Posterior Draws
// ============================================================================

struct MixedSimulationResult {
  arma::imat x_observations;
  arma::mat y_observations;
  int draw_index;
  bool error;
  std::string error_msg;
};


class MixedSimulationWorker : public RcppParallel::Worker {
public:
  // Input: posterior samples as flat vectors (one row per draw)
  const arma::mat& mux_samples;       // ndraws x n_mux
  const arma::mat& disc_samples;       // ndraws x p*(p-1)/2
  const arma::mat& muy_samples;       // ndraws x q
  const arma::mat& cont_samples;       // ndraws x q*(q+1)/2
  const arma::mat& cross_samples;       // ndraws x p*q

  const arma::ivec& draw_indices;
  const int num_states;
  const int p, q;
  const arma::ivec& num_categories;
  const std::vector<std::string>& variable_type;
  const arma::ivec& baseline_category;
  const int iter;
  const arma::ivec& mux_param_counts;

  const std::vector<SafeRNG>& draw_rngs;
  ProgressManager& pm;
  std::vector<MixedSimulationResult>& results;

  MixedSimulationWorker(
    const arma::mat& mux_samples,
    const arma::mat& disc_samples,
    const arma::mat& muy_samples,
    const arma::mat& cont_samples,
    const arma::mat& cross_samples,
    const arma::ivec& draw_indices,
    int num_states,
    int p, int q,
    const arma::ivec& num_categories,
    const std::vector<std::string>& variable_type,
    const arma::ivec& baseline_category,
    int iter,
    const arma::ivec& mux_param_counts,
    const std::vector<SafeRNG>& draw_rngs,
    ProgressManager& pm,
    std::vector<MixedSimulationResult>& results
  ) :
    mux_samples(mux_samples),
    disc_samples(disc_samples),
    muy_samples(muy_samples),
    cont_samples(cont_samples),
    cross_samples(cross_samples),
    draw_indices(draw_indices),
    num_states(num_states),
    p(p), q(q),
    num_categories(num_categories),
    variable_type(variable_type),
    baseline_category(baseline_category),
    iter(iter),
    mux_param_counts(mux_param_counts),
    draw_rngs(draw_rngs),
    pm(pm),
    results(results)
  {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      if (pm.shouldExit()) return;

      MixedSimulationResult result;
      result.draw_index = draw_indices[i];
      result.error = false;

      try {
        SafeRNG rng = draw_rngs[i];
        int draw_row = draw_indices[i] - 1;  // 1-based to 0-based

        // Reconstruct mux (p x max_cats)
        int max_cats = arma::max(mux_param_counts);
        arma::mat mux(p, max_cats, arma::fill::zeros);
        int idx = 0;
        for (int s = 0; s < p; s++) {
          for (int c = 0; c < mux_param_counts[s]; c++) {
            mux(s, c) = mux_samples(draw_row, idx);
            idx++;
          }
        }

        // Reconstruct discrete interactions (p x p symmetric, zero diagonal)
        arma::mat pairwise_disc(p, p, arma::fill::zeros);
        idx = 0;
        for (int col = 0; col < p; col++) {
          for (int row = col + 1; row < p; row++) {
            pairwise_disc(row, col) = disc_samples(draw_row, idx);
            pairwise_disc(col, row) = disc_samples(draw_row, idx);
            idx++;
          }
        }

        // Reconstruct muy (q-vector)
        arma::vec muy(q);
        for (int j = 0; j < q; j++) {
          muy(j) = muy_samples(draw_row, j);
        }

        // Reconstruct continuous interactions (q x q symmetric, including diagonal)
        arma::mat pairwise_cont(q, q, arma::fill::zeros);
        idx = 0;
        for (int col = 0; col < q; col++) {
          for (int row = col; row < q; row++) {
            if (row == col) {
              pairwise_cont(row, col) = cont_samples(draw_row, idx);
            } else {
              pairwise_cont(row, col) = cont_samples(draw_row, idx);
              pairwise_cont(col, row) = cont_samples(draw_row, idx);
            }
            idx++;
          }
        }

        // Reconstruct cross interactions (p x q, row-major)
        arma::mat pairwise_cross(p, q, arma::fill::zeros);
        idx = 0;
        for (int s = 0; s < p; s++) {
          for (int j = 0; j < q; j++) {
            pairwise_cross(s, j) = cross_samples(draw_row, idx);
            idx++;
          }
        }

        // Simulate
        result.x_observations.set_size(num_states, p);
        result.y_observations.set_size(num_states, q);

        simulate_mixed_mrf(
          num_states, pairwise_disc, pairwise_cross, pairwise_cont, mux, muy,
          num_categories, variable_type, baseline_category,
          iter, rng,
          result.x_observations, result.y_observations
        );

      } catch (const std::exception& e) {
        result.error = true;
        result.error_msg = e.what();
      } catch (...) {
        result.error = true;
        result.error_msg = "Unknown error";
      }

      results[i] = result;
      pm.update(0);
    }
  }
};


// Run parallel mixed MRF simulations across posterior draws. The R layer
// splits the flat parameter vector into 5 sample matrices: mux (ndraws x
// n_mux), disc (ndraws x p*(p-1)/2), muy (ndraws x q), cont (ndraws x
// q*(q+1)/2), cross (ndraws x p*q). Also takes 1-based draw indices,
// observations per draw, the discrete/continuous variable counts p and q,
// per-variable metadata (num_categories, "ordinal" or "blume-capel" type,
// baseline categories), Gibbs burn-in iterations, thread count, seed, and
// progress bar type. Returns a list of lists, each containing "x" (integer
// matrix) and "y" (numeric matrix).
// [[Rcpp::export]]
Rcpp::List run_mixed_simulation_parallel(
    const arma::mat& mux_samples,
    const arma::mat& disc_samples,
    const arma::mat& muy_samples,
    const arma::mat& cont_samples,
    const arma::mat& cross_samples,
    const arma::ivec& draw_indices,
    int num_states,
    int p,
    int q,
    const arma::ivec& num_categories,
    const Rcpp::StringVector& variable_type_r,
    const arma::ivec& baseline_category,
    int iter,
    int nThreads,
    int seed,
    int progress_type) {

  int ndraws = draw_indices.n_elem;

  std::vector<std::string> variable_type(p);
  arma::ivec baseline_category_safe = baseline_category;
  for (int i = 0; i < p; i++) {
    variable_type[i] = Rcpp::as<std::string>(variable_type_r[i]);
    if (variable_type[i] != "blume-capel") {
      baseline_category_safe[i] = 0;
    }
  }

  // Compute mux param counts per discrete variable
  arma::ivec mux_param_counts(p);
  for (int s = 0; s < p; s++) {
    if (variable_type[s] == "blume-capel") {
      mux_param_counts[s] = 2;
    } else {
      mux_param_counts[s] = num_categories[s];
    }
  }

  std::vector<SafeRNG> draw_rngs(ndraws);
  for (int d = 0; d < ndraws; d++) {
    draw_rngs[d] = SafeRNG(seed + d);
  }

  std::vector<MixedSimulationResult> results(ndraws);
  ProgressManager pm(1, ndraws, 0, 50, progress_type);

  MixedSimulationWorker worker(
    mux_samples, disc_samples, muy_samples, cont_samples, cross_samples,
    draw_indices, num_states, p, q,
    num_categories, variable_type, baseline_category_safe,
    iter, mux_param_counts, draw_rngs, pm, results
  );

  {
    tbb::global_control control(
      tbb::global_control::max_allowed_parallelism, nThreads);
    parallelFor(0, ndraws, worker);
  }
  pm.finish();

  // See the OMRF runner: an interrupted worker leaves its remaining draws
  // empty, which must not reach the caller as a result.
  if (pm.shouldExit()) {
    Rcpp::stop("Simulation interrupted by user.");
  }

  Rcpp::List output(ndraws);
  for (int i = 0; i < ndraws; i++) {
    if (results[i].error) {
      Rcpp::stop("Error in mixed MRF simulation draw %d: %s",
                 results[i].draw_index, results[i].error_msg.c_str());
    }
    output[i] = Rcpp::List::create(
      Rcpp::Named("x") = Rcpp::wrap(results[i].x_observations),
      Rcpp::Named("y") = Rcpp::wrap(results[i].y_observations)
    );
  }

  return output;
}
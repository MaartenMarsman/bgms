#pragma once

#include <RcppArmadillo.h>
#include <cmath>
#include <algorithm>
#include "mcmc/execution/warmup_schedule.h"
#include "math/explog_macros.h"


/**
 * DualAveraging - Step-size adaptation via Nesterov's dual averaging
 *
 * Tracks a running estimate of the optimal log step size by integrating
 * acceptance-probability error against a target rate. Used by
 * NUTSAdaptationController during warmup stages.
 */
class DualAveraging {
public:
  /// Current log step size.
  double log_step_size;
  /// Smoothed log step size (final estimate).
  double log_step_size_avg;
  /// Running error statistic.
  double hbar;
  /// Bias term: log(10 * initial_step_size).
  double mu;
  /// Shrinkage parameter (default 0.05).
  double gamma;
  /// Stabilisation offset (default 10).
  double t0;
  /// Decay exponent for averaging weights (default 0.75).
  double kappa;
  /// Iteration counter.
  int t;

  DualAveraging(double initial_step_size)
    : log_step_size(MY_LOG(initial_step_size)),
      log_step_size_avg(MY_LOG(initial_step_size)),
      hbar(0.0),
      mu(MY_LOG(10.0 * initial_step_size)),
      gamma(0.05),
      t0(10.0),
      kappa(0.75),
      t(1) {}

  /**
   * Integrate one acceptance-probability observation and update the
   * current and averaged log step sizes.
   *
   * @param accept_prob   Observed acceptance probability
   * @param target_accept Target acceptance rate
   */
  void update(double accept_prob, double target_accept) {
    double eta = 1.0 / (t + t0);
    double error = target_accept - accept_prob;
    hbar = (1 - eta) * hbar + eta * error;
    log_step_size = mu - std::sqrt(t) / gamma * hbar;

    double weight = std::pow(t, -kappa);
    log_step_size_avg = weight * log_step_size + (1.0 - weight) * log_step_size_avg;
    t++;
  }

  /** Reset the adaptation state around a new step size. */
  void restart(double new_step_size) {
    log_step_size = MY_LOG(new_step_size);
    log_step_size_avg = MY_LOG(new_step_size);
    mu = MY_LOG(10.0 * new_step_size);
    hbar = 0.0;
    t = 1;
  }

  /** Current step size. */
  double current() const { return MY_EXP(log_step_size); }
  /** Smoothed step size (final estimate). */
  double averaged() const { return MY_EXP(log_step_size_avg); }
};


/**
 * DiagMassMatrixAccumulator - Online diagonal mass matrix estimator
 *
 * Accumulates Welford running variance of parameter samples, blended
 * with a weak prior to prevent degenerate estimates. Used by
 * NUTSAdaptationController within Stage-2 windows.
 */
class DiagMassMatrixAccumulator {
public:
  /// Number of samples accumulated.
  int count;
  /// Running mean of parameter samples.
  arma::vec mean;
  /// Running sum of squared deviations (Welford M2 statistic).
  arma::vec m2;

  DiagMassMatrixAccumulator(int dim)
    : count(0), mean(arma::zeros(dim)), m2(arma::zeros(dim)) {}

  /** Add one parameter sample to the Welford accumulator. */
  void update(const arma::vec& sample) {
    count++;
    arma::vec delta = sample - mean;
    mean += delta / count;
    arma::vec delta2 = sample - mean;
    m2 += delta % delta2;
  }

  /**
   * Per-dimension variance estimate, blended with a weak prior
   * (weight 5, variance 1e-3) to prevent degenerate values.
   */
  arma::vec variance() const {
    static constexpr double prior_weight = 5.0;
    static constexpr double prior_variance = 1e-3;
    double n = static_cast<double>(count);

    arma::vec empirical = m2 / std::max(1.0, n - 1.0);
    arma::vec prior = arma::ones(empirical.n_elem) * prior_variance;
    arma::vec var = (n / (n + prior_weight)) * empirical
    + (prior_weight / (n + prior_weight)) * prior;
    return var;
  }

  /** Clear the accumulator for the next adaptation window. */
  void reset() {
    count = 0;
    mean.zeros();
    m2.zeros();
  }
};

/**
 * NUTSAdaptationController - Warmup adaptation for NUTS
 *
 * Coordinates step-size dual averaging (Stages 1, 2, 3a, 3c) and
 * mass-matrix estimation in doubling windows (Stage 2). Step size
 * is frozen at the Stage 3b boundary.
 */
class NUTSAdaptationController {
public:
  NUTSAdaptationController(int dim,
                          double initial_step_size,
                          double target_accept,
                          WarmupSchedule& schedule_ref,
                          bool learn_mass_matrix = true)
    : schedule(schedule_ref),
      learn_mass_matrix_(learn_mass_matrix),
      mass_accumulator(dim),
      step_adapter(initial_step_size),
      inv_mass_(arma::ones<arma::vec>(dim)),
      step_size_(initial_step_size),
      target_accept_(target_accept),
      mass_matrix_updated_(false) {}

  /**
   * Advance adaptation one iteration: dual-averaging step-size update in
   * Stages 1/2/3a/3c, mass-matrix accumulation within Stage-2 windows,
   * and step-size freeze at the Stage-3b boundary.
   *
   * @param theta       Current parameter vector
   * @param accept_prob Acceptance probability of the last NUTS transition
   * @param iteration   Current warmup iteration
   */
  void update(const arma::vec& theta,
              double accept_prob,
              int iteration) {
    /* ---------------------------------------------------------
     * 1. STEP-SIZE ADAPTATION
     *    – runs in Stage-1, Stage-2, Stage-3a, and Stage-3c
     *    – Stage-3c continues adaptation because selection changes
     *      the model structure (including/excluding parameters)
     *    – pauses in Stage-3b (proposal SD learning phase)
     * --------------------------------------------------------- */
    if (schedule.in_stage1(iteration)  || schedule.in_stage2(iteration)  ||
    schedule.in_stage3a(iteration) || schedule.in_stage3c(iteration) )
    {
      step_adapter.update(accept_prob, target_accept_);
      step_size_ = step_adapter.current();
    }

    /* ---------------------------------------------------------
     * 2. MASS-MATRIX ADAPTATION
     *    – only while we are inside Stage-2
     * --------------------------------------------------------- */
    if (schedule.in_stage2(iteration) && learn_mass_matrix_) {
      mass_accumulator.update(theta);
      int w = schedule.current_window(iteration);
      if (iteration + 1 == schedule.window_ends[w]) {
        // inv_mass = variance (not 1/variance)
        // Higher variance → higher inverse mass → parameter moves more freely
        inv_mass_ = mass_accumulator.variance();
        mass_accumulator.reset();
        // Signal that mass matrix was updated - caller should run heuristic
        // and call reinit_stepsize() with the new step size
        mass_matrix_updated_ = true;
      }
    }

    /* ---------------------------------------------------------
     * 3. FREEZE ε AS SOON AS WE ENTER STAGE-3b or SAMPLING
     * --------------------------------------------------------- */
    if (iteration == schedule.stage3b_start || schedule.sampling(iteration)) {
      step_size_ = step_adapter.averaged();
    }
  }

  /** Step size in effect for the current iteration. */
  double current_step_size() const { return step_size_; }
  /** Smoothed step size used after warmup. */
  double final_step_size() const { return step_adapter.averaged(); }
  /** Diagonal of the inverse mass matrix (per-parameter variances). */
  const arma::vec& inv_mass_diag() const { return inv_mass_; }

  /**
   * Check if the mass matrix was just updated and needs step size re-initialization.
   * After calling this, call reinit_stepsize() with the result of the heuristic.
   */
  bool mass_matrix_just_updated() const { return mass_matrix_updated_; }

  /** Target acceptance probability for dual averaging and heuristic. */
  double target_acceptance() const { return target_accept_; }

  /**
   * Reinitialize step size adaptation after mass matrix update.
   * This should be called after running heuristic_initial_step_size() with
   * the new mass matrix to find an appropriate starting step size.
   *
   * - Set the new step size
   * - Set mu = log(10 * new_step_size) as the adaptation target
   * - Restart the dual averaging counters
   */
  void reinit_stepsize(double new_step_size) {
    step_size_ = new_step_size;
    step_adapter.restart(new_step_size);
    // Set mu to log(10 * epsilon) for dual averaging
    step_adapter.mu = MY_LOG(10.0 * new_step_size);
    mass_matrix_updated_ = false;
  }

private:
  WarmupSchedule& schedule;
  bool learn_mass_matrix_;
  DiagMassMatrixAccumulator mass_accumulator;
  DualAveraging step_adapter;
  arma::vec inv_mass_;
  double step_size_;
  double target_accept_;
  bool mass_matrix_updated_;
};

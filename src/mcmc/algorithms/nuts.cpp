#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <utility>
#include "math/log_sum_exp.h"
#include "mcmc/algorithms/leapfrog.h"
#include "mcmc/algorithms/nuts.h"
#include "mcmc/algorithms/hmc.h"
#include "rng/rng_utils.h"


// Multinomial No-U-Turn Sampler. Candidates are weighted by
// log w_i = H0 - h_i (the Hamiltonian offset of each state from the initial
// state), and selection within and between subtrees uses log-sum-exp
// progressive / biased sampling as in Stan's base_nuts.hpp. The generalized
// U-turn criterion follows Betancourt (2017).
//
// References:
//   Betancourt, M. (2017). A Conceptual Introduction to Hamiltonian Monte
//     Carlo. arXiv preprint arXiv:1701.02434. (Section A.3.2 covers
//     multinomial expansion and biased progressive sampling.)
//   Stan Development Team. base_nuts.hpp (BSD-3-Clause).
//     https://github.com/stan-dev/stan/blob/develop/src/stan/mcmc/hmc/nuts/base_nuts.hpp


// Computes the generalized U-turn criterion for the NUTS algorithm from the
// sharp momenta (M^{-1} p) at the backward and forward trajectory ends and
// rho, the sum of momenta along the trajectory. Returns true if the
// criterion is satisfied (continue expanding), false if a U-turn is detected
// (stop).
bool compute_criterion(const arma::vec& p_sharp_minus,
                       const arma::vec& p_sharp_plus,
                       const arma::vec& rho) {
  return arma::dot(p_sharp_plus, rho) > 0 && arma::dot(p_sharp_minus, rho) > 0;
}



// Recursively builds a binary tree of leapfrog steps in the NUTS algorithm.
//
// Each base-case leaf contributes log-weight H0 - h to the subtree's
// log_sum_weight. Sibling subtrees combine by symmetric multinomial
// sampling in log-space (Stan's base_nuts.hpp).
//
// Takes the position theta and momentum r at the base of the tree, the
// expansion direction v (-1 backward, +1 forward), the tree depth j, the
// leapfrog step size, the Hamiltonian H0 at the start of the trajectory, a
// memoizer for caching evaluations, the diagonal of the inverse mass matrix,
// and the RNG for the log-sum-exp progressive draws. Returns a
// BuildTreeResult with updated endpoints, candidate sample, log weight, and
// diagnostics.
BuildTreeResult build_tree(
    const arma::vec& theta,
    const arma::vec& r,
    int v,
    int j,
    double step_size,
    double H0,
    Memoizer& memo,
    const arma::vec& inv_mass_diag,
    SafeRNG& rng
) {
  constexpr double Delta_max = 1000.0;
  constexpr double neg_inf = -std::numeric_limits<double>::infinity();

  if (j == 0) {
    // ---- Base case: a single leapfrog step --------------------------------
    arma::vec theta_new, r_new;
    std::tie(theta_new, r_new) = leapfrog_memo(
      theta, r, v * step_size, memo, inv_mass_diag
    );

    // Evaluate posterior, kinetic energy, Hamiltonian.
    double logp = memo.cached_log_post(theta_new);
    double kin = kinetic_energy(r_new, inv_mass_diag);
    double h = -logp + kin;
    // Guard NaN and inf symmetrically: either flags a broken state and will
    // trigger the divergence branch below.
    if (!std::isfinite(h)) {
      h = std::numeric_limits<double>::infinity();
    }
    double alpha = std::min(1.0, MY_EXP(H0 - h));  // 0 when h = +inf

    arma::vec p_sharp = inv_mass_diag % r_new;

    BuildTreeResult result;
    result.theta_min = theta_new;
    result.theta_plus = theta_new;
    result.r_min = r_new;
    result.r_plus = r_new;
    result.rho = r_new;
    result.p_beg = r_new;
    result.p_end = r_new;
    result.r_prime = std::move(r_new);
    result.theta_prime = std::move(theta_new);
    result.logp_prime = logp;
    result.p_sharp_beg = p_sharp;
    result.p_sharp_end = std::move(p_sharp);
    result.alpha = alpha;
    result.n_leapfrog = 1;

    // Branch 3: divergence check. A divergent leaf contributes no weight to
    // candidate selection (log_sum_weight = -inf) but its alpha is still
    // summed into the trajectory-level Metropolis diagnostic.
    if ((h - H0) > Delta_max) {
      result.log_sum_weight = neg_inf;
      result.s_prime = 0;
      result.divergent = true;
    } else {
      result.log_sum_weight = H0 - h;
      result.s_prime = 1;
      result.divergent = false;
    }
    return result;
  }

  // ---- Recursive case: build first subtree, then second -------------------
  BuildTreeResult init_result = build_tree(
    theta, r, v, j - 1, step_size, H0, memo, inv_mass_diag, rng
  );

  if (init_result.s_prime == 0) {
    // First subtree is invalid; propagate as-is (log_sum_weight already set).
    return init_result;
  }

  bool divergent = init_result.divergent;

  arma::vec theta_min = std::move(init_result.theta_min);
  arma::vec r_min = std::move(init_result.r_min);
  arma::vec theta_plus = std::move(init_result.theta_plus);
  arma::vec r_plus = std::move(init_result.r_plus);
  arma::vec theta_prime = std::move(init_result.theta_prime);
  arma::vec r_prime = std::move(init_result.r_prime);
  double logp_prime = init_result.logp_prime;
  arma::vec rho_init = std::move(init_result.rho);
  arma::vec p_sharp_init_beg = std::move(init_result.p_sharp_beg);
  arma::vec p_sharp_init_end = std::move(init_result.p_sharp_end);
  arma::vec p_init_beg = std::move(init_result.p_beg);
  arma::vec p_init_end = std::move(init_result.p_end);
  double log_sum_weight_init = init_result.log_sum_weight;
  double alpha_prime = init_result.alpha;
  int n_leapfrog_prime = init_result.n_leapfrog;

  // Build the second subtree in the same direction.
  BuildTreeResult final_result;
  if (v == -1) {
    final_result = build_tree(
      theta_min, r_min, v, j - 1, step_size, H0, memo, inv_mass_diag, rng
    );
    theta_min = std::move(final_result.theta_min);
    r_min = std::move(final_result.r_min);
  } else {
    final_result = build_tree(
      theta_plus, r_plus, v, j - 1, step_size, H0, memo, inv_mass_diag, rng
    );
    theta_plus = std::move(final_result.theta_plus);
    r_plus = std::move(final_result.r_plus);
  }

  // Accumulate Metropolis contributions regardless of subtree validity: every
  // leapfrog step happened and contributes to the trajectory-level diagnostic.
  alpha_prime += final_result.alpha;
  n_leapfrog_prime += final_result.n_leapfrog;
  divergent = divergent || final_result.divergent;

  if (final_result.s_prime == 0) {
    // Second subtree invalid — return early with s_prime=0. The returned
    // log_sum_weight sums both halves' contributions; the outer level will
    // use this value in its own biased-progression step. The existing
    // theta_prime (from the init subtree) is preserved because no valid
    // candidate from the final subtree can be combined here (it would be
    // dominated by the -inf weight anyway).
    BuildTreeResult result;
    result.theta_min = std::move(theta_min);
    result.r_min = std::move(r_min);
    result.theta_plus = std::move(theta_plus);
    result.r_plus = std::move(r_plus);
    result.theta_prime = std::move(theta_prime);
    result.r_prime = std::move(r_prime);
    result.logp_prime = logp_prime;
    rho_init += final_result.rho;
    result.rho = std::move(rho_init);
    result.p_sharp_beg = std::move(p_sharp_init_beg);
    result.p_sharp_end = std::move(final_result.p_sharp_end);
    result.p_beg = std::move(p_init_beg);
    result.p_end = std::move(final_result.p_end);
    result.log_sum_weight =
      log_sum_exp(log_sum_weight_init, final_result.log_sum_weight);
    result.s_prime = 0;
    result.alpha = alpha_prime;
    result.n_leapfrog = n_leapfrog_prime;
    result.divergent = divergent;
    return result;
  }

  arma::vec rho_final = std::move(final_result.rho);
  arma::vec p_sharp_final_beg = std::move(final_result.p_sharp_beg);
  arma::vec p_sharp_final_end = std::move(final_result.p_sharp_end);
  arma::vec p_final_beg = std::move(final_result.p_beg);
  arma::vec p_final_end = std::move(final_result.p_end);
  double log_sum_weight_final = final_result.log_sum_weight;

  // Symmetric multinomial sample between the two sibling subtrees.
  // log_sum_weight_subtree = log_sum_exp(init, final) >= max(init, final),
  // so the first branch can only fire when init = -inf (every leaf in the
  // init subtree diverged). Stan keeps the two branches for clarity; we
  // match that layout.
  double log_sum_weight_subtree =
    log_sum_exp(log_sum_weight_init, log_sum_weight_final);
  if (log_sum_weight_final > log_sum_weight_subtree) {
    theta_prime = std::move(final_result.theta_prime);
    r_prime = std::move(final_result.r_prime);
    logp_prime = final_result.logp_prime;
  } else if (MY_LOG(runif(rng)) <
             log_sum_weight_final - log_sum_weight_subtree) {
    theta_prime = std::move(final_result.theta_prime);
    r_prime = std::move(final_result.r_prime);
    logp_prime = final_result.logp_prime;
  }

  arma::vec rho_subtree = rho_init + rho_final;

  // Generalized U-turn criterion (three checks, Betancourt section A.4).
  bool persist_criterion =
    compute_criterion(p_sharp_init_beg, p_sharp_final_end, rho_subtree);
  arma::vec rho_extended = rho_init + p_final_beg;
  persist_criterion = persist_criterion &&
    compute_criterion(p_sharp_init_beg, p_sharp_final_beg, rho_extended);
  rho_extended = rho_final + p_init_end;
  persist_criterion = persist_criterion &&
    compute_criterion(p_sharp_init_end, p_sharp_final_end, rho_extended);

  BuildTreeResult result;
  result.theta_min = std::move(theta_min);
  result.r_min = std::move(r_min);
  result.theta_plus = std::move(theta_plus);
  result.r_plus = std::move(r_plus);
  result.theta_prime = std::move(theta_prime);
  result.r_prime = std::move(r_prime);
  result.logp_prime = logp_prime;
  result.rho = std::move(rho_subtree);
  result.p_sharp_beg = std::move(p_sharp_init_beg);
  result.p_sharp_end = std::move(p_sharp_final_end);
  result.p_beg = std::move(p_init_beg);
  result.p_end = std::move(p_final_end);
  result.log_sum_weight = log_sum_weight_subtree;
  result.s_prime = persist_criterion ? 1 : 0;
  result.alpha = alpha_prime;
  result.n_leapfrog = n_leapfrog_prime;
  result.divergent = divergent;
  return result;
}


StepResult nuts_step(
    const arma::vec& init_theta,
    double step_size,
    const std::function<std::pair<double, arma::vec>(const arma::vec&)>& joint,
    const arma::vec& inv_mass_diag,
    SafeRNG& rng,
    int max_depth
) {
  // Create Memoizer with joint function
  Memoizer memo(joint);
  bool any_divergence = false;

  arma::vec r0 = arma::sqrt(1.0 / inv_mass_diag) % arma_rnorm_vec(rng, init_theta.n_elem);

  double logp0 = memo.cached_log_post(init_theta);
  double kin0 = kinetic_energy(r0, inv_mass_diag);
  double H0 = -logp0 + kin0;
  // Log-posterior of the currently selected sample, kept in sync with theta.
  double logp_selected = logp0;

  arma::vec theta_min = init_theta, r_min = r0;
  arma::vec theta_plus = init_theta, r_plus = r0;
  arma::vec theta = init_theta;
  arma::vec r = r0;

  arma::vec p_sharp_bck_bck = inv_mass_diag % r0;
  arma::vec p_sharp_fwd_fwd = p_sharp_bck_bck;
  arma::vec p_fwd_bck = r0;
  arma::vec p_sharp_fwd_bck = p_sharp_bck_bck;
  arma::vec p_bck_fwd = r0;
  arma::vec p_sharp_bck_fwd = p_sharp_bck_bck;
  // Regular (non-sharp) boundary momenta for the extreme ends of the
  // trajectory. Needed for the Stan-style "save old end as new inner
  // boundary" maintenance before each extension.
  arma::vec p_bck_bck = r0;
  arma::vec p_fwd_fwd = r0;
  arma::vec rho = r0;

  int j = 0;
  int s = 1;
  // Log-sum of state weights exp(H0 - h_i) accumulated across the trajectory.
  // Initialized to 0 = log exp(H0 - H0) so the initial position carries unit
  // weight; the sampler stays put if every subtree diverges.
  double log_sum_weight_traj = 0.0;
  // Accumulate Metropolis contributions across every leapfrog step (including
  // inside subtrees that terminated). Matches Stan's by-reference accumulation
  // in base_nuts.hpp.
  double sum_metro_prob = 0.0;
  int n_leapfrog_total = 0;

  while (s == 1 && j < max_depth) {
    int v = runif(rng) < 0.5 ? -1 : 1;
    arma::vec rho_fwd, rho_bck;

    BuildTreeResult result;
    if (v == -1) {
      rho_fwd = rho;
      // Save old backward-end as the forward half's new backward boundary.
      // Matches Stan base_nuts.hpp: p_fwd_bck = p_bck_bck before building.
      p_fwd_bck = p_bck_bck;
      p_sharp_fwd_bck = p_sharp_bck_bck;
      result = build_tree(
        theta_min, r_min, v, j, step_size, H0, memo, inv_mass_diag, rng
      );
      theta_min = std::move(result.theta_min);
      r_min = std::move(result.r_min);
      rho_bck = std::move(result.rho);
      // For a backward subtree, p_beg = first leaf built = interior (closest
      // to origin), p_end = last leaf built = outer-leftmost. Map accordingly.
      p_sharp_bck_bck = std::move(result.p_sharp_end);
      p_bck_bck = std::move(result.p_end);
      p_sharp_bck_fwd = std::move(result.p_sharp_beg);
      p_bck_fwd = std::move(result.p_beg);
    } else {
      rho_bck = rho;
      // Save old forward-end as the backward half's new forward boundary.
      // Matches Stan base_nuts.hpp: p_bck_fwd = p_fwd_fwd before building.
      p_bck_fwd = p_fwd_fwd;
      p_sharp_bck_fwd = p_sharp_fwd_fwd;
      result = build_tree(
        theta_plus, r_plus, v, j, step_size, H0, memo, inv_mass_diag, rng
      );
      theta_plus = std::move(result.theta_plus);
      r_plus = std::move(result.r_plus);
      rho_fwd = std::move(result.rho);
      p_sharp_fwd_fwd = std::move(result.p_sharp_end);
      p_fwd_fwd = std::move(result.p_end);
      p_fwd_bck = std::move(result.p_beg);
      p_sharp_fwd_bck = std::move(result.p_sharp_beg);
    }

    // Trajectory-level bookkeeping: always accumulate (subtrees that
    // terminated still consumed leapfrog steps).
    any_divergence = any_divergence || result.divergent;
    sum_metro_prob += result.alpha;
    n_leapfrog_total += result.n_leapfrog;

    // Biased progressive sampling: prefer the newly built subtree when its
    // weight exceeds the running trajectory weight (accept with probability
    // min(1, w_subtree / w_traj) in the multinomial formulation). Only eligible
    // when the subtree itself is valid.
    if (result.s_prime == 1) {
      if (result.log_sum_weight > log_sum_weight_traj ||
          MY_LOG(runif(rng)) <
            result.log_sum_weight - log_sum_weight_traj) {
        theta = std::move(result.theta_prime);
        r = std::move(result.r_prime);
        logp_selected = result.logp_prime;
      }
    }
    log_sum_weight_traj =
      log_sum_exp(log_sum_weight_traj, result.log_sum_weight);

    rho = rho_bck + rho_fwd;
    bool persist_criterion = true;
    if (result.s_prime == 1) {
      persist_criterion = compute_criterion(p_sharp_bck_bck, p_sharp_fwd_fwd, rho);
      arma::vec rho_extended = rho_bck + p_fwd_bck;
      persist_criterion = persist_criterion &&
        compute_criterion(p_sharp_bck_bck, p_sharp_fwd_bck, rho_extended);
      rho_extended = rho_fwd + p_bck_fwd;
      persist_criterion = persist_criterion &&
        compute_criterion(p_sharp_bck_fwd, p_sharp_fwd_fwd, rho_extended);
    }

    s = result.s_prime * (persist_criterion ? 1 : 0);
    j++;
  }

  double accept_prob = sum_metro_prob / static_cast<double>(n_leapfrog_total);
  double kin_final = kinetic_energy(r, inv_mass_diag);
  double energy = -logp_selected + kin_final;

  auto diag = std::make_shared<NUTSDiagnostics>();
  diag->tree_depth = j;
  diag->divergent = any_divergence;
  diag->energy = energy;
  diag->accept_prob = accept_prob;

  return {theta, accept_prob, diag};
}

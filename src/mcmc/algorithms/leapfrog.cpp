#include <RcppArmadillo.h>
#include <functional>
#include <utility>
#include "mcmc/algorithms/leapfrog.h"


std::pair<arma::vec, arma::vec> leapfrog_memo(
    const arma::vec& theta,
    const arma::vec& r,
    double eps,
    Memoizer& memo,
    const arma::vec& inv_mass_diag
) {
  arma::vec r_half = r;
  arma::vec theta_new = theta;

  const arma::vec& grad1 = memo.cached_grad(theta_new);

  r_half += 0.5 * eps * grad1;
  theta_new += eps * (inv_mass_diag % r_half);

  const arma::vec& grad2 = memo.cached_grad(theta_new);

  r_half += 0.5 * eps * grad2;

  return {theta_new, r_half};
}


std::pair<arma::vec, arma::vec> leapfrog_constrained(
    const arma::vec& theta,
    const arma::vec& r,
    double eps,
    Memoizer& memo,
    const arma::vec& inv_mass_diag,
    const ProjectPositionFn& project_position,
    const ProjectMomentumFn& project_momentum
) {
  arma::vec r_half = r;
  arma::vec theta_new = theta;

  // --- Step 1: Half-step momentum ---
  const arma::vec& grad1 = memo.cached_grad(theta_new);
  r_half += 0.5 * eps * grad1;

  // --- Step 2: Project momentum onto cotangent space ---
  project_momentum(r_half, theta_new);

  // --- Step 3: Full-step position ---
  theta_new += eps * (inv_mass_diag % r_half);

  // --- Step 4: SHAKE — position-only projection ---
  arma::vec theta_pre = theta_new;
  project_position(theta_new);

  // --- Step 5: Momentum correction for constraint forces ---
  arma::vec delta_x = theta_new - theta_pre;
  r_half += delta_x / (eps * inv_mass_diag);

  memo.invalidate();

  // --- Step 6: Second half-step momentum ---
  const arma::vec& grad2 = memo.cached_grad(theta_new);
  r_half += 0.5 * eps * grad2;

  // --- Step 7: Project momentum onto cotangent space ---
  project_momentum(r_half, theta_new);

  return {theta_new, r_half};
}


ConstrainedLeapfrogResult leapfrog_constrained_checked(
    const arma::vec& theta,
    const arma::vec& r,
    double eps,
    Memoizer& memo,
    const arma::vec& inv_mass_diag,
    const ProjectPositionFn& project_position,
    const ProjectMomentumFn& project_momentum,
    double reverse_check_tol
) {
  // --- Forward step ---
  auto [theta_new, r_new] = leapfrog_constrained(
    theta, r, eps, memo, inv_mass_diag,
    project_position, project_momentum
  );

  // --- Backward step (negate momentum, step forward) ---
  // Use a separate Memoizer so the caller's cache stays at theta_new, but
  // seed it with the caller's cached evaluation: the forward step left memo
  // at theta_new, which is exactly where the backward step starts.
  Memoizer back_memo(memo.joint_fn);
  back_memo.seed_from(memo);
  arma::vec r_back = -r_new;
  auto [theta_back, r_back_out] = leapfrog_constrained(
    theta_new, r_back, eps, back_memo, inv_mass_diag,
    project_position, project_momentum
  );

  // --- Reversibility check (eps^2-scaled max-norm) ---
  double max_diff = arma::max(arma::abs(theta_back - theta));
  double tol = reverse_check_tol * eps * eps;
  bool reversible = (max_diff <= tol);

  return {std::move(theta_new), std::move(r_new), reversible};
}


LeapfrogJointResult leapfrog(
    const arma::vec& theta_init,
    const arma::vec& r_init,
    double eps,
    const std::function<arma::vec(const arma::vec&)>& grad,
    const std::function<std::pair<double, arma::vec>(const arma::vec&)>& joint,
    int num_leapfrogs,
    const arma::vec& inv_mass_diag,
    const arma::vec* init_grad
) {
  arma::vec r = r_init;
  arma::vec theta = theta_init;

  // Use provided initial gradient or compute it
  arma::vec grad_theta = init_grad ? *init_grad : grad(theta_init);

  // All steps except the last one
  for (int step = 0; step < num_leapfrogs - 1; step++) {
    // Half-step momentum
    r += 0.5 * eps * grad_theta;

    // Full step position
    theta += eps * (inv_mass_diag % r);

    // Update gradient (intermediate position - only need grad)
    grad_theta = grad(theta);

    // Final half-step momentum
    r += 0.5 * eps * grad_theta;
  }

  // Final step: use joint to get both log_post and gradient
  if (num_leapfrogs >= 1) {
    // Half-step momentum
    r += 0.5 * eps * grad_theta;

    // Full step position
    theta += eps * (inv_mass_diag % r);

    // Use joint at final position
    auto [log_post_final, grad_final] = joint(theta);

    // Final half-step momentum
    r += 0.5 * eps * grad_final;

    return {theta, r, log_post_final, grad_final};
  }

  // Edge case: num_leapfrogs == 0 (shouldn't happen in practice)
  auto [log_post, grad_vec] = joint(theta);
  return {theta, r, log_post, grad_vec};
}

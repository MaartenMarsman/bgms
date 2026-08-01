#include <RcppArmadillo.h>
#include <algorithm>
#include <functional>
#include "mcmc/algorithms/hamiltonian_utils.h"
#include "mcmc/algorithms/leapfrog.h"
#include "math/explog_macros.h"
#include "rng/rng_utils.h"


double kinetic_energy(const arma::vec& r, const arma::vec& inv_mass_diag) {
  return 0.5 * arma::dot(r % inv_mass_diag, r);
}


double heuristic_initial_step_size(
    const arma::vec& theta,
    const std::function<arma::vec(const arma::vec&)>& grad,
    const std::function<std::pair<double, arma::vec>(const arma::vec&)>& joint,
    SafeRNG& rng,
    double target_acceptance,
    double init_step,
    int max_attempts
) {
  arma::vec inv_mass_diag = arma::ones<arma::vec>(theta.n_elem);
  return heuristic_initial_step_size(
    theta, grad, joint, inv_mass_diag, rng, target_acceptance, init_step, max_attempts
  );
}


double heuristic_initial_step_size(
    const arma::vec& theta,
    const std::function<arma::vec(const arma::vec&)>& grad,
    const std::function<std::pair<double, arma::vec>(const arma::vec&)>& joint,
    const arma::vec& inv_mass_diag,
    SafeRNG& rng,
    double target_acceptance,
    double init_step,
    int max_attempts
) {
  double eps = init_step;

  // Use joint at initial position (compute once - position doesn't change)
  auto [logp0, grad0] = joint(theta);

  // Sample initial momentum and evaluate
  arma::vec r = arma::sqrt(1.0 / inv_mass_diag) % arma_rnorm_vec(rng, theta.n_elem);
  double kin0 = kinetic_energy(r, inv_mass_diag);
  double H0 = logp0 - kin0;

  // One leapfrog step using leapfrog
  LeapfrogJointResult result = leapfrog(
    theta, r, eps, grad, joint, 1, inv_mass_diag, &grad0
  );

  double kin1 = kinetic_energy(result.r, inv_mass_diag);
  double H1 = result.log_post - kin1;

  double log_target = MY_LOG(target_acceptance);
  int direction = 2 * (H1 - H0 > log_target) - 1;  // +1 or -1

  int attempts = 0;
  while (direction * (H1 - H0) > direction * log_target && attempts < max_attempts) {
    eps = (direction == 1) ? 2.0 * eps : 0.5 * eps;

    // Resample momentum on each iteration for step size search
    r = arma::sqrt(1.0 / inv_mass_diag) % arma_rnorm_vec(rng, theta.n_elem);
    kin0 = kinetic_energy(r, inv_mass_diag);
    H0 = logp0 - kin0;

    // One leapfrog step from original position with new momentum
    result = leapfrog(theta, r, eps, grad, joint, 1, inv_mass_diag, &grad0);

    // Evaluate Hamiltonian
    kin1 = kinetic_energy(result.r, inv_mass_diag);
    H1 = result.log_post - kin1;

    attempts++;
  }

  return eps;
}

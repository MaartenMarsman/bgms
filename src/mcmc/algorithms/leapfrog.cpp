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


LeapfrogJointResult leapfrog(
    const arma::vec& theta_init,
    const arma::vec& r_init,
    double eps,
    const std::function<arma::vec(const arma::vec&)>& grad,
    const std::function<std::pair<double, arma::vec>(const arma::vec&)>& joint,
    const arma::vec& inv_mass_diag,
    const arma::vec* init_grad
) {
  arma::vec r = r_init;
  arma::vec theta = theta_init;

  // Use provided initial gradient or compute it
  arma::vec grad_theta = init_grad ? *init_grad : grad(theta_init);

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

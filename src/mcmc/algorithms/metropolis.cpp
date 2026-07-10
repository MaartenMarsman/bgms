#include <RcppArmadillo.h>
#include <algorithm>
#include <functional>
#include "mcmc/execution/step_result.h"
#include "mcmc/algorithms/metropolis.h"
#include "rng/rng_utils.h"


StepResult metropolis_step(
    double current_state,
    double step_size,
    const std::function<double(double)>& log_post,
    SafeRNG& rng
) {
  double proposed_state = rnorm(rng, current_state, step_size);
  double log_accept = log_post(proposed_state) - log_post(current_state);
  double accept_prob = std::min(1.0, MY_EXP(log_accept));

  double state = (runif(rng) < accept_prob) ? proposed_state : current_state;

  arma::vec State(1);
  State[0] = state;

  return {State, accept_prob};
}


StepResult metropolis_step_cached(
    double current_state,
    double step_size,
    double log_post_current,
    const std::function<double(double)>& log_post_proposed,
    SafeRNG& rng
) {
  double proposed_state = rnorm(rng, current_state, step_size);
  double log_accept = log_post_proposed(proposed_state) - log_post_current;
  double accept_prob = std::min(1.0, MY_EXP(log_accept));

  double state = (runif(rng) < accept_prob) ? proposed_state : current_state;

  arma::vec State(1);
  State[0] = state;

  return {State, accept_prob};
}

#pragma once

#include <RcppArmadillo.h>
#include <cmath>
#include "mcmc/execution/step_result.h"
#include "mcmc/execution/warmup_schedule.h"
#include "math/explog_macros.h"


/**
 * MetropolisAdaptationController - Robbins-Monro adaptation for Metropolis proposals
 *
 * Adjusts proposal standard deviations during warmup using acceptance
 * probabilities. Applied element-wise to a matrix of proposal SDs.
 */
class MetropolisAdaptationController {
public:
  /// Reference to the proposal standard deviation matrix (modified in place).
  arma::mat& proposal_sd;
  /// Total number of warmup iterations.
  const int total_warmup;
  /// Target Metropolis acceptance rate.
  const double target_accept;

  MetropolisAdaptationController(arma::mat& proposal_sd_matrix,
                          const WarmupSchedule& warmup,
                          double target_accept_rate = 0.44)
    : proposal_sd(proposal_sd_matrix),
      total_warmup(warmup.total_warmup),
      target_accept(target_accept_rate) {}

  /**
   * Robbins-Monro update of the proposal SDs for the masked entries.
   * A no-op outside the warmup phase.
   *
   * @param index_mask         Entries equal to 1 are updated
   * @param accept_prob_matrix Acceptance probabilities per entry
   * @param iteration          Current iteration (1-based within warmup)
   */
  void update(const arma::umat& index_mask,
              const arma::mat& accept_prob_matrix,
              int iteration) {

    if (iteration >= total_warmup || iteration < 1)
      return;

    const double rm_decay_rate = 0.75;
    const double rm_weight = std::pow(iteration, -rm_decay_rate);

    for (arma::uword i = 0; i < proposal_sd.n_rows; ++i) {
      for (arma::uword j = 0; j < proposal_sd.n_cols; ++j) {
        if (index_mask(i, j) == 1) {
          const double accept_log = MY_LOG(accept_prob_matrix(i, j));
          proposal_sd(i, j) = update_proposal_sd_with_robbins_monro(
            proposal_sd(i, j),
            accept_log,
            rm_weight,
            target_accept
          );
        }
      }
    }
  }
};

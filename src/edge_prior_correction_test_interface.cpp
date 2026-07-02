// -----------------------------------------------------------------------------
// edge_prior_correction_test_interface.cpp
//
// R-facing test entries for the normalizing-constant correction on the
// beta-bernoulli edge prior: the log C(theta) interpolator and the corrected
// inverse-CDF theta draw, so both can be tested against closed forms and
// quadrature without running a chain.
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>

#include "priors/edge_prior_correction.h"
#include "rng/rng_utils.h"

// -----------------------------------------------------------------------------
// test_correction_logC_interp:
//   Evaluate the table interpolator (linear inside the grid, linear-slope
//   extension beyond it) at the supplied points.
// -----------------------------------------------------------------------------

// [[Rcpp::export(name = "test_correction_logC_interp")]]
arma::vec test_correction_logC_interp(
    const arma::vec& theta_grid,
    const arma::vec& logC,
    const arma::vec& theta_eval
) {
    EdgePriorCorrection correction(theta_grid, logC);
    arma::vec out(theta_eval.n_elem);
    for (arma::uword k = 0; k < theta_eval.n_elem; ++k) {
        out[k] = correction.log_C(theta_eval[k]);
    }
    return out;
}

// -----------------------------------------------------------------------------
// test_corrected_bb_theta_draw:
//   Run the corrected theta draw as a chain of n_draws updates at fixed
//   conjugate shapes (a_post, b_post), starting from theta_init. The draw
//   window is centered at the current value, so consecutive draws form a
//   Markov chain targeting the corrected conditional.
// -----------------------------------------------------------------------------

// [[Rcpp::export(name = "test_corrected_bb_theta_draw")]]
arma::vec test_corrected_bb_theta_draw(
    double a_post,
    double b_post,
    const arma::vec& theta_grid,
    const arma::vec& logC,
    int n_draws,
    int seed,
    double theta_init = 0.5
) {
    EdgePriorCorrection correction(theta_grid, logC);
    SafeRNG rng(seed);
    arma::vec draws(n_draws);
    double current = theta_init;
    for (int k = 0; k < n_draws; ++k) {
        current = correction.draw_theta(rng, a_post, b_post, current);
        draws[k] = current;
    }
    return draws;
}

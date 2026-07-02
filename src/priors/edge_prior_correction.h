#pragma once

#include <RcppArmadillo.h>
#include "rng/rng_utils.h"

/**
 * Normalizing-constant correction for hierarchical edge priors on the GGM.
 *
 * Under the determinant-tilted joint spike-and-slab prior, the graph-prior
 * mean of the per-graph normalizer, C(theta), enters the full conditional of
 * the edge-prior hyperparameters, so the plain conjugate updates are biased.
 * This class holds the whole-graph log C(theta) curve (built at fit setup
 * from the tilted prior sampler; see R/correction_tables.R) and draws from
 * the corrected conditional
 *
 *   p(theta | Gamma) ∝ theta^(a-1) (1-theta)^(b-1) exp(-log C(theta)),
 *
 * with (a, b) the conjugate posterior shape parameters. The draw is
 * inverse-CDF on a fine grid (N = 401) centered at the current theta with
 * half-width 10 conjugate posterior standard deviations: the window scales
 * with the posterior, so resolution-in-sd is constant as theta approaches
 * the boundaries. log C enters only as differences, so the curve's additive
 * constant is immaterial. Outside the tabulated range the curve is extended
 * linearly with the boundary slope.
 */
class EdgePriorCorrection {
public:
    EdgePriorCorrection() = default;

    EdgePriorCorrection(const arma::vec& theta, const arma::vec& log_C)
        : theta_(theta), log_C_(log_C)
    {
        if (theta_.n_elem != log_C_.n_elem || theta_.n_elem < 2) {
            Rcpp::stop("EdgePriorCorrection: theta and log_C must have equal length >= 2.");
        }
    }

    bool active() const { return theta_.n_elem >= 2; }

    /** Linear interpolation of log C(theta), linear extension beyond the grid. */
    double log_C(double th) const {
        const arma::uword n = theta_.n_elem;
        if (th <= theta_[0]) {
            return log_C_[0] +
                (log_C_[1] - log_C_[0]) / (theta_[1] - theta_[0]) * (th - theta_[0]);
        }
        if (th >= theta_[n - 1]) {
            return log_C_[n - 1] +
                (log_C_[n - 1] - log_C_[n - 2]) / (theta_[n - 1] - theta_[n - 2]) *
                (th - theta_[n - 1]);
        }
        arma::uword lo = 0, hi = n - 1;
        while (hi - lo > 1) {
            arma::uword mid = (lo + hi) / 2;
            if (theta_[mid] <= th) lo = mid; else hi = mid;
        }
        double t = (th - theta_[lo]) / (theta_[hi] - theta_[lo]);
        return log_C_[lo] + t * (log_C_[hi] - log_C_[lo]);
    }

    /**
     * Draw theta from the corrected conditional.
     *
     * @param rng            Chain RNG.
     * @param a_post         Conjugate shape alpha + #included.
     * @param b_post         Conjugate shape beta + #excluded.
     * @param theta_current  Current theta (window center).
     */
    double draw_theta(SafeRNG& rng, double a_post, double b_post,
                      double theta_current) const {
        const int N = 401;
        double mu = a_post / (a_post + b_post);
        double sd = std::sqrt(mu * (1.0 - mu) / (a_post + b_post + 1.0));
        double half_width = std::min(0.49, std::max(1e-4, 10.0 * sd));
        double lo = std::max(1e-7, theta_current - half_width);
        double hi = std::min(1.0 - 1e-7, theta_current + half_width);

        double grid[N];
        double weight[N];
        double step = (hi - lo) / (N - 1);
        double max_logp = -std::numeric_limits<double>::infinity();
        for (int k = 0; k < N; k++) {
            double th = lo + step * k;
            grid[k] = th;
            weight[k] = (a_post - 1.0) * std::log(th) +
                (b_post - 1.0) * std::log1p(-th) - log_C(th);
            if (weight[k] > max_logp) max_logp = weight[k];
        }
        double total = 0.0;
        for (int k = 0; k < N; k++) {
            weight[k] = std::exp(weight[k] - max_logp);
            total += weight[k];
        }
        double u = runif(rng) * total;
        double cumulative = 0.0;
        int k = 0;
        for (; k < N - 1; k++) {
            cumulative += weight[k];
            if (cumulative >= u) break;
        }
        return grid[k];
    }

private:
    arma::vec theta_;
    arma::vec log_C_;
};

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


/**
 * Normalizing-constant correction data for the stochastic block edge prior.
 *
 * The block-model updates read the correction locally: the per-edge tilt
 * slope is f'(d) evaluated at the edge's min-endpoint expected degree
 * density, and the new-cluster collapsed marginal integrates the per-pair
 * curve f(theta) over a uniform quadrature grid. Holds
 *
 *   - (fprime_density, fprime): the slope f' tabulated against local
 *     density, extended as constants beyond the tabulated range;
 *   - (quad_theta, quad_f): the per-pair f(theta) curve pre-interpolated
 *     onto a uniform quadrature grid.
 *
 * On the mixed-MRF path the determinant tilt acts on the continuous
 * precision block only, so the tilt terms read continuous-continuous pairs
 * exclusively; is_continuous carries the per-node mask in model order
 * (discrete block first). An empty mask means all nodes are continuous
 * (the GGM path).
 *
 * Built from the same correction table as the beta-bernoulli update; see
 * R/correction_tables.R.
 */
class SBMCorrection {
public:
    SBMCorrection() = default;

    SBMCorrection(const arma::vec& fprime_density, const arma::vec& fprime,
                  const arma::vec& quad_theta, const arma::vec& quad_f,
                  const arma::uvec& is_continuous = arma::uvec())
        : fprime_density_(fprime_density), fprime_(fprime),
          quad_theta_(quad_theta), quad_f_(quad_f),
          is_continuous_(is_continuous)
    {
        if (fprime_density_.n_elem != fprime_.n_elem ||
            fprime_density_.n_elem < 2 ||
            quad_theta_.n_elem != quad_f_.n_elem || quad_theta_.n_elem < 2) {
            Rcpp::stop("SBMCorrection: curve inputs must have equal length >= 2.");
        }
    }

    bool active() const { return fprime_.n_elem >= 2; }

    /** Whether node i sits in the continuous block (empty mask: all do). */
    bool node_continuous(arma::uword i) const {
        return is_continuous_.n_elem == 0 || is_continuous_(i) != 0;
    }

    /** Number of continuous nodes among no_variables total. */
    arma::uword num_continuous(arma::uword no_variables) const {
        if (is_continuous_.n_elem == 0) return no_variables;
        return arma::accu(is_continuous_ != 0);
    }

    const arma::uvec& is_continuous() const { return is_continuous_; }

    /** Slope f' at local density d: linear interpolation, constant extension. */
    double fprime_at(double d) const {
        const arma::uword n = fprime_density_.n_elem;
        if (d <= fprime_density_[0]) return fprime_[0];
        if (d >= fprime_density_[n - 1]) return fprime_[n - 1];
        arma::uword lo = 0, hi = n - 1;
        while (hi - lo > 1) {
            arma::uword mid = (lo + hi) / 2;
            if (fprime_density_[mid] <= d) lo = mid; else hi = mid;
        }
        double t = (d - fprime_density_[lo]) /
            (fprime_density_[hi] - fprime_density_[lo]);
        return fprime_[lo] + t * (fprime_[hi] - fprime_[lo]);
    }

    const arma::vec& quad_theta() const { return quad_theta_; }
    const arma::vec& quad_f() const { return quad_f_; }

private:
    arma::vec fprime_density_;
    arma::vec fprime_;
    arma::vec quad_theta_;
    arma::vec quad_f_;
    arma::uvec is_continuous_;
};

#pragma once

#include <string>
#include <vector>
#include <RcppArmadillo.h>

/**
 * ChainResult - Storage for a single MCMC chain's output
 *
 * Holds samples, diagnostics, and error state for one chain.
 * Designed for use with both Metropolis and NUTS samplers.
 */
class ChainResult {

public:
    ChainResult() = default;

    /// True if the chain terminated with an error.
    bool        error = false;
    /// True if the chain was interrupted by the user.
    bool        userInterrupt = false;
    /// Error message (empty if none).
    std::string error_msg;

    /// Integer identifier for the chain (1-based).
    int         chain_id = 0;

    /// Parameter samples (param_dim x n_iter).
    arma::mat   samples;

    /// Edge indicator samples (n_edges x n_iter), only if edge_selection = true.
    arma::imat  indicator_samples;
    /// Whether indicator samples are stored.
    bool        has_indicators = false;

    /// SBM allocation samples (n_variables x n_iter), only if SBM edge prior.
    arma::imat  allocation_samples;
    /// Whether allocation samples are stored.
    bool        has_allocations = false;

    /// Sampled inclusion parameter (Beta-Bernoulli theta; n_iter).
    arma::vec   inclusion_parameter_samples;
    /// Whether inclusion-parameter samples are stored.
    bool        has_inclusion_parameter = false;

    /// NUTS tree depth diagnostics (n_iter).
    arma::ivec  treedepth_samples;
    /// NUTS divergent transition flags (n_iter).
    arma::ivec  divergent_samples;
    /// NUTS energy diagnostic (n_iter).
    arma::vec   energy_samples;
    /// NUTS mean per-trajectory Metropolis acceptance (n_iter).
    arma::vec   accept_prob_samples;
    /// Whether NUTS diagnostics are stored.
    bool        has_nuts_diagnostics = false;

    /// Adaptive-Metropolis mean per-iteration acceptance probability across
    /// all updated components (n_iter).
    arma::vec   am_accept_prob_samples;
    /// Whether AM diagnostics are stored.
    bool        has_am_diagnostics = false;

    /// Z-ratio engine constant block at end of run (frozen fit + hull when
    /// the calibrator produced one; the 6-slot additive block otherwise).
    arma::vec   zratio_addc;
    /// Z-ratio calibration anchor design rows (1, bre, m, cne, maxbd, dens).
    arma::mat   zratio_anchors_x;
    /// Z-ratio calibration anchor targets log(oracle) - log(additive).
    arma::vec   zratio_anchors_y;
    /// Z-ratio engine counters: n_hit, n_miss, n_pred, n_add, n_clamp,
    /// n_oracle, n_anchors, cache_size, frozen (0/1).
    arma::vec   zratio_counters;
    /// Graph density trace over the selection-enabled warmup iterations
    /// (recorded only when a Z-ratio calibration window is configured).
    std::vector<double> zratio_warmup_density;
    /// Edge-prior inclusion parameter trace over the same iterations.
    std::vector<double> zratio_warmup_theta;
    /// Whether Z-ratio engine diagnostics are stored.
    bool        has_zratio_diagnostics = false;

    /// In-chain trust gauge (hierarchical spec): D = fraction of edge moves
    /// that would flip under the exact ratio, its reference-noise floor, the
    /// signed-mean and spread of the log-ratio error s_e, and the pair
    /// counts (non-trivial seen / referenced / cap hits). Populated only when
    /// the gauge ran (gauge_ran = true).
    double      zratio_gauge_D = 0.0;
    double      zratio_gauge_noise_floor = 0.0;
    double      zratio_gauge_se_mean = 0.0;
    double      zratio_gauge_se_sd = 0.0;
    double      zratio_gauge_se_mcse = 0.0;
    long        zratio_gauge_n_ent = 0;
    long        zratio_gauge_n_ref = 0;
    long        zratio_gauge_n_capped = 0;
    bool        zratio_gauge_ran = false;
    /// Per-referenced-pair audit stream (edge endpoints, signed log-ratio
    /// error, reference MCSE); bounded by the per-sweep cap times the number
    /// of gauge sweeps.
    arma::ivec  zratio_gauge_pair_i;
    arma::ivec  zratio_gauge_pair_j;
    arma::vec   zratio_gauge_pair_se;
    arma::vec   zratio_gauge_pair_mcse;

    /**
     * Reserve storage for samples
     * @param param_dim  Number of parameters per sample
     * @param n_iter     Number of sampling iterations
     */
    void reserve(const size_t param_dim, const size_t n_iter) {
        samples.set_size(param_dim, n_iter);
        samples.fill(arma::datum::nan);
    }

    /**
     * Reserve storage for edge indicator samples
     * @param n_edges  Number of edges (p * (p - 1) / 2)
     * @param n_iter   Number of sampling iterations
     */
    void reserve_indicators(const size_t n_edges, const size_t n_iter) {
        indicator_samples.set_size(n_edges, n_iter);
        indicator_samples.fill(-1);
        has_indicators = true;
    }

    /**
     * Reserve storage for SBM allocation samples
     * @param n_variables  Number of variables
     * @param n_iter       Number of sampling iterations
     */
    void reserve_allocations(const size_t n_variables, const size_t n_iter) {
        allocation_samples.set_size(n_variables, n_iter);
        allocation_samples.fill(-1);
        has_allocations = true;
    }

    /**
     * Reserve storage for inclusion-parameter samples
     * @param n_iter  Number of sampling iterations
     */
    void reserve_inclusion_parameter(const size_t n_iter) {
        inclusion_parameter_samples.set_size(n_iter);
        inclusion_parameter_samples.fill(arma::datum::nan);
        has_inclusion_parameter = true;
    }

    /**
     * Reserve storage for NUTS diagnostics
     * @param n_iter  Number of sampling iterations
     */
    void reserve_nuts_diagnostics(const size_t n_iter) {
        // Integer diagnostics use -1 as a "not sampled" sentinel; the
        // floating ones use NaN. This keeps interrupted runs from returning
        // uninitialized values for iterations that never ran.
        treedepth_samples.set_size(n_iter);
        treedepth_samples.fill(-1);
        divergent_samples.set_size(n_iter);
        divergent_samples.fill(-1);
        energy_samples.set_size(n_iter);
        energy_samples.fill(arma::datum::nan);
        accept_prob_samples.set_size(n_iter);
        accept_prob_samples.fill(arma::datum::nan);
        has_nuts_diagnostics = true;
    }

    /**
     * Reserve storage for adaptive-Metropolis diagnostics
     * @param n_iter  Number of sampling iterations
     */
    void reserve_am_diagnostics(const size_t n_iter) {
        am_accept_prob_samples.set_size(n_iter);
        am_accept_prob_samples.fill(arma::datum::nan);
        has_am_diagnostics = true;
    }

    /**
     * Store a parameter sample
     * @param iter    Iteration index (0-based)
     * @param sample  Parameter vector
     */
    void store_sample(const size_t iter, const arma::vec& sample) {
        samples.col(iter) = sample;
    }

    /**
     * Store edge indicator sample
     * @param iter        Iteration index (0-based)
     * @param indicators  Edge indicator vector
     */
    void store_indicators(const size_t iter, const arma::ivec& indicators) {
        indicator_samples.col(iter) = indicators;
    }

    /**
     * Store SBM allocation sample
     * @param iter         Iteration index (0-based)
     * @param allocations  Allocation vector (1-based cluster labels)
     */
    void store_allocations(const size_t iter, const arma::ivec& allocations) {
        allocation_samples.col(iter) = allocations;
    }

    /**
     * Store the inclusion-parameter sample for one iteration
     * @param iter   Iteration index (0-based)
     * @param value  Sampled inclusion parameter (Beta-Bernoulli theta)
     */
    void store_inclusion_parameter(const size_t iter, const double value) {
        inclusion_parameter_samples(iter) = value;
    }

    /**
     * Store NUTS diagnostics for one iteration
     * @param iter         Iteration index (0-based)
     * @param tree_depth   Tree depth from NUTS
     * @param divergent    Whether a divergence occurred
     * @param energy       Final Hamiltonian energy
     * @param accept_prob  Mean Metropolis acceptance over the trajectory
     */
    void store_nuts_diagnostics(const size_t iter, int tree_depth, bool divergent, double energy, double accept_prob) {
        treedepth_samples(iter) = tree_depth;
        divergent_samples(iter) = divergent ? 1 : 0;
        energy_samples(iter) = energy;
        accept_prob_samples(iter) = accept_prob;
    }

    /**
     * Store adaptive-Metropolis diagnostics for one iteration
     * @param iter         Iteration index (0-based)
     * @param accept_prob  Mean acceptance probability across the sweep
     */
    void store_am_diagnostics(const size_t iter, double accept_prob) {
        am_accept_prob_samples(iter) = accept_prob;
    }
};

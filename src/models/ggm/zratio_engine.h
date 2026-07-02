#pragma once

#include <RcppArmadillo.h>
#include <unordered_map>
#include <string>

#include "rng/rng_utils.h"

/**
 * Deterministic per-edge normalizing-constant ratio for the hierarchical
 * prior specification p(K | Gamma) = rho_Gamma(K) / Z(Gamma).
 *
 * Every between-edge move carries J = Z(Gamma-)/Z(Gamma+). The engine
 * evaluates J from three integer counts of the toggled edge's mediating
 * neighbourhood (common-neighbour nodes, edges among them, bridge edges
 * between the exclusive neighbour sets) through a two-moment saddle over
 * cosine-transform tables of the tilted prior's pair integrals. No sampling
 * runs inside the acceptance step; distinct count tuples are evaluated once
 * and served from a persistent cache.
 *
 * Constant block `addc` (0-based):
 *   [0..5]   per-channel moment constants (CN node, CN-CN edge, bridge),
 *            built at fit time (R/zratio_tables.R).
 *   [6..11]  optional OLS correction (intercept, bre, m, cne, maxbdeg,
 *            dens), fit by the warm-up calibrator; engaged when the edge
 *            has bridge multiplicity >= 2.
 *   [12]     > 0.5 selects the direct ratio-scale correction (log J += fc
 *            after the cached saddle); the cache then keys on the counts
 *            only.
 *   [13..22] optional hull box (per-feature min/max of the calibration
 *            design, order bre, m, cne, maxbd, dens); outside the box the
 *            correction is zeroed so the frozen kernel never extrapolates.
 *
 * Conventions (bare scale): K_ii ~ Exp(beta), slab K_ij ~ N(0, sigma^2),
 * tilt |K|^delta; sigma = 2 * pairwise_scale, beta = scale_rate / 2 in bgms
 * parameter units. Reference: SV/Z sbc_prior_chain_exact.cpp and the
 * z_graph_prior deployed kernel (hier_chain_data.cpp).
 */
class ZRatioEngine {
public:
    ZRatioEngine(const arma::vec& addc, const arma::vec& tg,
                 const arma::vec& ihat, const arma::vec& ghat,
                 const arma::vec& wt, double psi0)
        : addc_(addc), tg_(tg), ihat_(ihat), ghat_(ghat), wt_(wt),
          psi0_(psi0) {}

    /**
     * log J = log( Z(Gamma-)/Z(Gamma+) ) for the edge (i, j) on graph G.
     * Independent of the toggled edge's own state, so one value serves the
     * add move and (negated) the delete move.
     */
    double log_zratio(const arma::imat& G, int i, int j);

    /**
     * Pre-evaluate the saddle for every count tuple in the bounded box
     * (nCN <= ncn_max, cne <= nCN(nCN-1)/2, bre <= bre_max) into the cache.
     */
    void precompute_table(int ncn_max, int bre_max);

    /** Two-moment saddle map over the cosine-transform grid. */
    double saddle_ratio(double s1, double s2) const;

    /**
     * Enable online calibration of the OLS correction during warm-up.
     *
     * While unfrozen, coupled-bridge blocks (maxbd >= 2) route through the
     * calibrator: identical block signatures are served from a cache,
     * blocks inside the anchor cloud's Mahalanobis hull use the current
     * fit, and uncovered blocks call the block-Gibbs local oracle, add an
     * anchor, and refit. freeze_calibration() packs the fit and its hull
     * box into the addc layout, after which the engine behaves exactly
     * like one constructed with a full 23-slot constant block.
     *
     * (delta, sigma, beta) are the bare-scale prior constants the oracle
     * samples under; rng must outlive the engine (the model's chain RNG).
     */
    void enable_calibration(double delta, double sigma, double beta,
                            SafeRNG* rng, int n_sweep = 300, int burn = 30,
                            double maha_thresh = 9.0, int min_anchors = 6);

    /** Refit and freeze: pack coefficients + hull box into addc[6..22]. */
    void freeze_calibration();

    bool calibrating() const { return calibration_enabled_ && !frozen_; }
    void set_rng(SafeRNG* rng) { rng_ = rng; }

    /**
     * Weighted block moments (S1, S2) from a block-Gibbs run on the
     * mediating block: n_sweep sweeps of the row-wise conjugate sampler
     * targeting the tilted block prior, Rao-Blackwellized through the
     * cross-resolvent trace moments. Returns false when no sweep yields a
     * finite moment pair.
     */
    bool block_oracle_moments(const arma::imat& a_blk, const arma::uvec& si,
                              const arma::uvec& sj, double& s1_out,
                              double& s2_out);

    long cache_size() const { return static_cast<long>(cache_.size()); }
    long n_hit() const { return n_hit_; }
    long n_miss() const { return n_miss_; }
    long n_pred() const { return n_pred_; }
    long n_add() const { return n_add_; }
    long n_clamp() const { return n_clamp_; }
    long n_oracle() const { return n_oracle_; }
    long n_anchors() const { return static_cast<long>(ay_.n_elem); }
    bool frozen() const { return frozen_; }
    const arma::vec& addc() const { return addc_; }

private:
    void gibbs_sweep_(arma::mat& k_blk,
                      const std::vector<arma::uvec>& nbr) const;
    bool inner_moments_(const arma::mat& k_blk, const arma::uvec& si,
                        const arma::uvec& sj, double& w, double& p1,
                        double& p2) const;
    void refit_();

    arma::vec addc_;
    arma::vec tg_, ihat_, ghat_, wt_;
    double psi0_;
    std::unordered_map<std::string, double> cache_;
    long n_hit_ = 0, n_miss_ = 0;
    long n_pred_ = 0, n_add_ = 0, n_clamp_ = 0;

    // Online-calibration state (inert unless enable_calibration ran).
    bool calibration_enabled_ = false;
    bool frozen_ = true;
    double delta_ = 0.0, sigma_ = 1.0, beta_ = 0.5;
    SafeRNG* rng_ = nullptr;
    int n_sweep_ = 300, burn_ = 30;
    double maha_thresh_ = 9.0;
    int min_anchors_ = 6;
    arma::mat ax_;                 // anchors: rows (1,bre,m,cne,maxbd,dens)
    arma::vec ay_;                 // anchors: log(oracle) - log(additive)
    arma::vec coef_;               // current OLS fit (empty before min_anchors)
    arma::vec maha_mu_;
    arma::mat maha_sinv_;
    std::unordered_map<std::string, double> corr_cache_;
    long n_oracle_ = 0;
};

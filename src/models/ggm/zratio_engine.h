#pragma once

#include <RcppArmadillo.h>
#include <unordered_map>
#include <string>
#include <cstdint>

#include "rng/rng_utils.h"

/**
 * Mediating-block descriptors for one candidate edge (i, j) on a graph.
 *
 * The block collects the common neighbours of the endpoints plus the
 * endpoints of 2-hop bridges between the exclusive neighbour sets; the
 * integer counts drive the additive saddle and the OLS correction, and
 * the adjacency + side memberships drive the block-Gibbs oracle.
 */
struct ZRatioBlock {
    bool valid = false;   ///< false: one side empty (isolated-edge ratio)
    int m = 0;            ///< number of mediating nodes
    int ncn = 0;          ///< common-neighbour nodes
    int cne = 0;          ///< edges among the common neighbours
    int bre = 0;          ///< bridge edges between the exclusive sides
    int maxbd = 0;        ///< maximum bridge degree over both sides
    double dens = 0.0;    ///< block edge density
    arma::imat a_blk;     ///< block adjacency (m x m)
    arma::uvec si;        ///< block rows adjacent to endpoint i
    arma::uvec sj;        ///< block rows adjacent to endpoint j
};

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
 * Conventions (standardized cell): K_ii ~ Exp(beta), slab K_ij ~ N(0,
 * sigma^2), tilt |K|^delta. The between-graph ratio is invariant under the
 * diagonal congruence Theta = A K A, so the constants are built at
 * sigma = 1, beta = eta = pairwise_scale * scale_rate in bgms parameter
 * units (R/zratio_tables.R, zratio_cell_constants), and the same cell
 * serves every user scale choice. Reference: SV/Z
 * sbc_prior_chain_exact.cpp and the z_graph_prior deployed kernel
 * (hier_chain_data.cpp).
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
     * Extract the mediating block of the edge (i, j): common neighbours,
     * 2-hop bridge endpoints, adjacency, side memberships, and the integer
     * counts (nCN, cne, bre, maxbd) with the block density. The toggled
     * edge's own state never enters.
     */
    ZRatioBlock extract_block(const arma::imat& G, int i, int j) const;

    /**
     * Deployed OLS correction for one block under the current constant
     * block: 0 for maxbd < 2 or when no fit is packed (addc[12] <= 0.5),
     * the OLS value inside the hull box, 0 outside it (clamped = true).
     * Reads state only; no counters move.
     */
    double deployed_correction(const ZRatioBlock& bl, bool& clamped) const;

    /**
     * Measurement-only audit of the edge (i, j): the deployed correction
     * versus the block-Gibbs local oracle on the same block.
     *
     * On success fills pred_out (deployed correction), oracle_out
     * (log(saddle on oracle moments) - log(additive saddle)) and bl_out,
     * and returns true. Returns false for one-sided blocks, non-positive
     * additive moments, or an oracle with no finite sweep. Requires
     * set_oracle_params (or enable_calibration) to have run; feeds
     * nothing back into the fit, caches, or counters.
     */
    bool audit_edge(const arma::imat& G, int i, int j, double& pred_out,
                    double& oracle_out, ZRatioBlock& bl_out);

    /**
     * Set the standardized-cell prior constants and RNG the block-Gibbs
     * oracle samples under, without entering calibration mode. rng must
     * outlive the engine.
     */
    void set_oracle_params(double delta, double sigma, double beta,
                           SafeRNG* rng, int n_sweep = 300, int burn = 30) {
        delta_ = delta;
        sigma_ = sigma;
        beta_ = beta;
        rng_ = rng;
        n_sweep_ = n_sweep;
        burn_ = burn;
    }

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
     * (delta, sigma, beta) are the standardized-cell prior constants the
     * oracle samples under; rng must outlive the engine (the model's
     * chain RNG).
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
    /// Calibration anchor design rows (1, bre, m, cne, maxbd, dens).
    const arma::mat& anchors_x() const { return ax_; }
    /// Calibration anchor targets log(oracle) - log(additive).
    const arma::vec& anchors_y() const { return ay_; }

private:
    void gibbs_sweep_(arma::mat& k_blk,
                      const std::vector<arma::uvec>& nbr) const;
    bool inner_moments_(const arma::mat& k_blk, const arma::uvec& si,
                        const arma::uvec& sj, double& w, double& p1,
                        double& p2) const;
    void refit_();

    /** Pack the (nCN, cne, bre) additive-cache counts into one integer key. */
    static std::uint64_t pack_count_key(int ncn, int cne, int bre) {
        return (static_cast<std::uint64_t>(ncn) << 42) |
               (static_cast<std::uint64_t>(cne) << 21) |
               static_cast<std::uint64_t>(bre);
    }

    arma::vec addc_;
    arma::vec tg_, ihat_, ghat_, wt_;
    double psi0_;
    std::unordered_map<std::uint64_t, double> cache_;
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

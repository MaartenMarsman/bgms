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
 * the adjacency + side memberships drive the exact Monte-Carlo evaluation.
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
 * One family's (CN or bipartite) theta-independent absolute-moment surface
 * (Option B). log(S) is a raw bivariate quadratic in (L, d), L = log(size),
 * d = component density, over the 9 monomials
 *   [1, L, L^2, d, d^2, L*d, L^2*d, L*d^2, L^2*d^2]
 * fit once per analysis at the deployment (eta, delta) to block-Gibbs anchors
 * (R build_surfaces). c1 predicts log-S1, c2 log-S2. Predictions clamp (size,
 * density) to the trained hull [size_lo, size_hi] x [dens_lo, dens_hi] and the
 * log-moment to its trained range +/- 0.1. Components smaller than size_min
 * fall back to the additive per-component moment (exact through pairwise
 * overlap below the smallest trained size); with size_min = 3 the size-1/2
 * (single-bridge) trivial components land there and additive == exact for them,
 * so there is no separate exact branch.
 */
struct SurfaceFamily {
    bool valid = false;
    arma::vec c1;   ///< 9 raw-poly coefficients for log-S1
    arma::vec c2;   ///< 9 raw-poly coefficients for log-S2
    double size_lo = 0.0, size_hi = 0.0;   ///< trained size hull (raw, pre-log)
    double dens_lo = 0.0, dens_hi = 0.0;   ///< trained density hull
    double l1_lo = 0.0, l1_hi = 0.0;       ///< trained log-S1 range (clamp +/-0.1)
    double l2_lo = 0.0, l2_hi = 0.0;       ///< trained log-S2 range
    double size_min = 0.0;                 ///< below this: additive fallback
};

/**
 * One decomposed mediating-block component, for surface deployment and for the
 * Stage-3 test harness. family 0 = common-neighbour cluster, 1 = bipartite
 * bridge. (size, dens) index the family surface; (e, na, nb) drive the additive
 * fallback and diagnostics. used_surface records which path served it.
 */
struct SurfaceComp {
    int family = 0;
    int size = 0;
    int e = 0;
    int na = 0, nb = 0;
    double dens = 0.0;
    double s1 = 0.0, s2 = 0.0;
    bool used_surface = false;
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
 *            only. Applied for every block with bridge multiplicity >= 2
 *            (the single deployment gate); never re-gated on the block's
 *            location in feature space.
 *   [13..22] hull box of the calibration design (per-feature min/max,
 *            order bre, m, cne, maxbd, dens), retained for diagnostics only.
 *            The frozen kernel no longer gates on it: reverting to the biased
 *            additive saddle outside a prior/size-dependent box is worse than
 *            extrapolating the smooth ratio-scale surface, so the correction
 *            extends past the training cloud.
 *
 * Conventions (standardized cell): K_ii ~ Gamma(alpha, beta) (alpha = 1 is
 * the exponential default), slab K_ij ~ N(0, sigma^2), tilt |K|^delta. The
 * between-graph ratio is invariant under the diagonal congruence
 * Theta = A K A, so the constants are built at sigma = 1,
 * beta = eta = pairwise_scale * scale_rate in bgms parameter units
 * (R/zratio_tables.R, zratio_cell_constants), and the same
 * (delta, eta, alpha) cell serves every user scale choice.
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
     * Attach the Option-B absolute-moment surfaces. Once set (and the cell is
     * the validated alpha = 1 Normal-slab family), log_zratio decomposes the
     * mediating block into disjoint components and sums the per-component
     * surface moments into (S1, S2) instead of the additive-counts saddle. The
     * additive + OLS path stays in place and is served whenever the surface is
     * absent or the cell is fenced (alpha != 1 / Cauchy), so no OLS machinery is
     * removed here.
     */
    void set_surface(const SurfaceFamily& cn, const SurfaceFamily& bip) {
        surf_cn_ = cn;
        surf_bip_ = bip;
        has_surface_ = cn.valid && bip.valid;
    }
    bool has_surface() const { return has_surface_; }

    /**
     * Surface-path moments for the edge (i, j): decompose the mediating block,
     * accumulate per-component (S1, S2), and return logR = log saddle_ratio.
     * Fills `comps` with the decomposition for inspection. Returns false when
     * the block is invalid (isolated-edge ratio) — logR is then log(psi0).
     * Drives the Stage-3 test harness; the hot path uses the lighter branch in
     * log_zratio.
     */
    bool surface_moments(const arma::imat& G, int i, int j, double& s1_out,
                         double& s2_out, double& logr_out,
                         std::vector<SurfaceComp>& comps);

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
     * otherwise the OLS value (applied everywhere, including past the
     * calibration cloud; clamped is always false, kept for interface
     * stability). Reads state only; no counters move.
     */
    double deployed_correction(const ZRatioBlock& bl, bool& clamped) const;

    /**
     * Set the standardized-cell prior constants and RNG the block-Gibbs
     * oracle samples under, without entering calibration mode. The frame is
     * standardized (unit slab scale), so only the diagonal rate eta is free;
     * sigma is fixed to 1. rng must outlive the engine. slab_cauchy selects
     * the Cauchy slab family: the block couplings then run omega-augmented
     * (scale-mixture of normals) and the endpoint legs mix per sweep,
     * matching the marginal-Cauchy normalizer the tables integrate.
     */
    void set_oracle_params(double delta, double eta, SafeRNG* rng,
                           int n_sweep = 300, int burn = 30,
                           bool slab_cauchy = false, double alpha = 1.0) {
        delta_ = delta;
        sigma_ = 1.0;
        beta_ = eta;
        alpha_ = alpha;
        rng_ = rng;
        n_sweep_ = n_sweep;
        burn_ = burn;
        oracle_slab_cauchy_ = slab_cauchy;
    }

    /**
     * Enable online calibration of the OLS correction during warm-up.
     *
     * While unfrozen, coupled-bridge blocks (maxbd >= 2) route through the
     * calibrator: identical block signatures are served from a cache;
     * otherwise a leverage gate decides between the fit and the oracle.
     * With x the block's design row and X the anchor design, the gate
     * serves the current fit when the standardized prediction variance
     * v(x) = x' (X'X + ridge)^{-1} x is at most gate_kappa, and calls the
     * block-Gibbs local oracle (adding an anchor and refitting) when it
     * exceeds it — so oracle runs concentrate where the surface is still
     * uncertain and stop once the visited feature space is spanned, drift
     * or no drift. max_anchors is a hard backstop on oracle runs.
     * freeze_calibration() packs the fit and the anchor feature box into
     * the addc layout, after which the engine behaves exactly like one
     * constructed with a full 23-slot constant block.
     *
     * (delta, eta, alpha) are the standardized-cell prior constants the
     * oracle samples under (unit slab scale, diagonal rate eta, diagonal
     * Gamma shape alpha); rng must outlive the engine (the model's chain
     * RNG). slab_cauchy selects the Cauchy slab family for the oracle (see
     * set_oracle_params).
     */
    void enable_calibration(double delta, double eta, SafeRNG* rng,
                            int n_sweep = 100, int burn = 30,
                            double gate_kappa = 1.0, int min_anchors = 6,
                            bool slab_cauchy = false, double alpha = 1.0,
                            int max_anchors = 100);

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

    /**
     * Block-local EXACT reference for the trust gauge: log R_e where
     * R_e = mean_H{ W(H) <phi, I_N> } / mean_H{ W(H) <phi, I_G> } over
     * n_draws rest-block precision draws H = K_R^{-1} of the mediating
     * block. Unlike block_oracle_moments (two-moment saddle collapse) the
     * endpoint couplings are integrated with the full product transform
     * phi(t|H) = prod_k (1 + u_k t^2)^{-1/2}, u_k = sigma^4 s_k^2 over ALL
     * singular values s_k of the cross-resolvent. Common random numbers
     * across the two averages, endpoints analytic. Fills logR_out and
     * mcse_out (batch-means MC standard error on the log scale). Returns
     * false when no draw yields a finite pair. Draws from the live rng_;
     * requires set_oracle_params / enable_calibration to have set the
     * standardized-cell prior constants.
     */
    bool block_reference_logR(const arma::imat& a_blk, const arma::uvec& si,
                              const arma::uvec& sj, int n_draws,
                              double& logR_out, double& mcse_out);

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
    /**
     * Shared worker behind extract_block. With counts_only the scalar
     * descriptors (valid, m, ncn, cne, bre, maxbd, dens) are read off G
     * directly and the block adjacency and side-membership vectors are
     * left empty; the values are identical to the full extraction.
     */
    void extract_block_(const arma::imat& G, int i, int j, bool counts_only,
                        ZRatioBlock& bl) const;
    /**
     * One row-wise sweep of the block sampler. When sigma_out is non-null
     * it receives the end-of-sweep block covariance k_blk^{-1} (the
     * SMW-maintained Sigma, exact up to within-sweep drift); returns
     * whether sigma_out holds a valid inverse (always true when sigma_out
     * is null).
     */
    bool gibbs_sweep_(arma::mat& k_blk, arma::mat& omega_blk,
                      const std::vector<arma::uvec>& nbr,
                      arma::mat* sigma_out = nullptr) const;
    /** Build neighbour lists, seed k_blk (+omega_blk under Cauchy), burn. */
    void init_block_(const arma::imat& a_blk, std::vector<arma::uvec>& nbr,
                     arma::mat& k_blk, arma::mat& omega_blk) const;
    bool inner_moments_(const arma::mat& r_inv, const arma::uvec& si,
                        const arma::uvec& sj, const arma::vec& wsi,
                        const arma::vec& wsj, double& w, double& p1,
                        double& p2) const;
    /** Per-draw full-product endpoint integrals for block_reference_logR. */
    bool inner_reference_(const arma::mat& r_inv, const arma::uvec& si,
                          const arma::uvec& sj, const arma::vec& wsi,
                          const arma::vec& wsj, double& w, double& fN,
                          double& gG, double& kappa2) const;
    void refit_();

    /**
     * exp(clamped raw-quadratic prediction) of a family surface at (size,
     * dens). s2 selects the log-S2 coefficients/range, else log-S1. Clamps
     * (size, dens) to the trained hull and the log-moment to its range +/- 0.1,
     * matching the R deploy predictor exactly.
     */
    double surface_eval_(const SurfaceFamily& f, bool s2, double size,
                         double dens) const;
    /**
     * Decompose the block into disjoint CN clusters and bipartite bridge
     * structures, accumulate per-component (S1, S2) — surface above size_min,
     * additive below — and (when comps != nullptr) record each component.
     * Requires the full block extraction (a_blk, si, sj).
     */
    void accumulate_surface_moments_(const ZRatioBlock& bl, double& s1_out,
                                     double& s2_out,
                                     std::vector<SurfaceComp>* comps) const;

    /** Pack the (nCN, cne, bre) additive-cache counts into one integer key. */
    static std::uint64_t pack_count_key(int ncn, int cne, int bre) {
        return (static_cast<std::uint64_t>(ncn) << 42) |
               (static_cast<std::uint64_t>(cne) << 21) |
               static_cast<std::uint64_t>(bre);
    }

    arma::vec addc_;
    arma::vec tg_, ihat_, ghat_, wt_;
    double psi0_;

    // Option-B absolute-moment surfaces (inert unless set_surface ran).
    bool has_surface_ = false;
    SurfaceFamily surf_cn_, surf_bip_;
    std::unordered_map<std::uint64_t, double> cache_;
    long n_hit_ = 0, n_miss_ = 0;
    long n_pred_ = 0, n_add_ = 0, n_clamp_ = 0;

    // Online-calibration state (inert unless enable_calibration ran).
    bool calibration_enabled_ = false;
    bool frozen_ = true;
    bool oracle_slab_cauchy_ = false;
    double delta_ = 0.0, sigma_ = 1.0, beta_ = 0.5, alpha_ = 1.0;
    SafeRNG* rng_ = nullptr;
    int n_sweep_ = 300, burn_ = 30;
    double gate_kappa_ = 1.0;
    int min_anchors_ = 6;
    int max_anchors_ = 100;
    arma::mat ax_;                 // anchors: rows (1,bre,m,cne,maxbd,dens)
    arma::vec ay_;                 // anchors: log(oracle) - log(additive)
    arma::vec coef_;               // current OLS fit (empty before min_anchors)
    arma::mat xtx_inv_;            // (X'X + ridge)^{-1} for the leverage gate
    std::unordered_map<std::string, double> corr_cache_;
    long n_oracle_ = 0;

    // Reused scratch for extract_block_. One engine serves one chain, and
    // the counts pass runs once per edge proposal, so per-call heap
    // allocation dominates the extraction cost without these.
    mutable std::vector<char> xb_in_r_;
    mutable std::vector<int> xb_excl_i_, xb_excl_j_, xb_rv_, xb_cn_, xb_sio_,
        xb_sjo_;
    mutable std::vector<unsigned char> xb_side_;
};

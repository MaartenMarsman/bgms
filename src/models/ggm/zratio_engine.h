#pragma once

#include <RcppArmadillo.h>
#include <unordered_map>
#include <array>
#include <utility>
#include <string>
#include <cstdint>

#include "rng/rng_utils.h"

/**
 * Mediating-block descriptors for one candidate edge (i, j) on a graph.
 *
 * The block collects the common neighbours of the endpoints plus the
 * endpoints of 2-hop bridges between the exclusive neighbour sets; the
 * integer counts drive the additive saddle, and the adjacency + side
 * memberships drive the exact Monte-Carlo evaluation and the surface deploy.
 */
struct ZRatioBlock {
    bool valid = false;   ///< false: one side empty (isolated-edge ratio)
    int m = 0;            ///< number of mediating nodes
    int ncn = 0;          ///< common-neighbour nodes
    int cne = 0;          ///< edges among the common neighbours
    int bre = 0;          ///< bridge edges between the exclusive sides
    int maxbd = 0;        ///< maximum bridge degree over both sides
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
 * One decomposed component in block-position coordinates, shared by the surface
 * deploy and the gold (per-component oracle) reference so both score the
 * identical decomposition. family 0 = CN cluster, 1 = bipartite bridge; nodes
 * are block-adjacency row indices; aside (bipartite only) marks the A-side per
 * node; (na, nb, e) are the side sizes and edge count.
 */
struct BlockComponent {
    int family = 0;
    std::vector<int> nodes;
    std::vector<char> aside;
    int na = 0, nb = 0, e = 0;
};

/** Build a SurfaceFamily from its R list (9 coeffs, hulls, ranges, size_min). */
inline SurfaceFamily surface_family_from_list(const Rcpp::List& s) {
    SurfaceFamily f;
    f.c1 = Rcpp::as<arma::vec>(s["c1"]);
    f.c2 = Rcpp::as<arma::vec>(s["c2"]);
    f.size_lo = Rcpp::as<double>(s["size_lo"]);
    f.size_hi = Rcpp::as<double>(s["size_hi"]);
    f.dens_lo = Rcpp::as<double>(s["dens_lo"]);
    f.dens_hi = Rcpp::as<double>(s["dens_hi"]);
    f.l1_lo = Rcpp::as<double>(s["l1_lo"]);
    f.l1_hi = Rcpp::as<double>(s["l1_hi"]);
    f.l2_lo = Rcpp::as<double>(s["l2_lo"]);
    f.l2_hi = Rcpp::as<double>(s["l2_hi"]);
    f.size_min = Rcpp::as<double>(s["size_min"]);
    f.valid = f.c1.n_elem == 9 && f.c2.n_elem == 9;
    return f;
}

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
 * Constant block `addc` (0-based, 6 slots): per-channel moment constants
 * [0,1] CN node (S1, S2), [2,3] CN-CN edge, [4,5] bridge, built at fit time
 * (R/zratio_tables.R). The alpha = 1 cell is corrected by the theta-independent
 * absolute-moment surface (set_surface); the alpha != 1 (Gamma-shape) cell
 * falls back to the additive saddle over these constants.
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
     * Attach the Option-B absolute-moment surfaces. Once set (at the validated
     * alpha = 1 diagonal, Normal or Cauchy slab), log_zratio decomposes the
     * mediating block into disjoint components and sums the per-component
     * surface moments into (S1, S2) instead of the additive-counts saddle. The
     * additive-counts saddle is served whenever the surface is absent or the
     * cell is fenced (a non-unit Gamma diagonal shape, alpha != 1).
     */
    void set_surface(const SurfaceFamily& cn, const SurfaceFamily& bip) {
        surf_cn_ = cn;
        surf_bip_ = bip;
        has_surface_ = cn.valid && bip.valid;
        // Attaching a new surface invalidates the deploy-time caches, whose
        // values are functions of the currently-attached surfaces.
        surf_cache_.clear();
        comp_cache_.clear();
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
     * Gold reference for the edge (i, j) on the SAME decomposition the surface
     * scores: each non-trivial component's moments come from the block-Gibbs
     * oracle (block_oracle_moments, at the set_oracle_params sweep count),
     * trivial components (CN size <= 2, single bridge) from the exact additive
     * kernel; the sum feeds saddle_ratio. Isolates the surface's moment
     * prediction from the closure, matching the companion's gold. Returns false
     * for an invalid (isolated-edge) block.
     */
    bool gold_moments(const arma::imat& G, int i, int j, double& s1_out,
                      double& s2_out, double& logr_out);

    /**
     * Extract the mediating block of the edge (i, j): common neighbours,
     * 2-hop bridge endpoints, adjacency, side memberships, and the integer
     * counts (nCN, cne, bre, maxbd) with the block density. The toggled
     * edge's own state never enters.
     */
    ZRatioBlock extract_block(const arma::imat& G, int i, int j) const;

    /**
     * Set the standardized-cell prior constants and RNG the block-Gibbs
     * oracle samples under. The frame is standardized (unit slab scale), so
     * only the diagonal rate eta is free; sigma is fixed to 1. rng must
     * outlive the engine. slab_cauchy selects the Cauchy slab family: the
     * block couplings then run omega-augmented (scale-mixture of normals) and
     * the endpoint legs mix per sweep, matching the marginal-Cauchy normalizer
     * the tables integrate. Used by the offline surface build, the gold
     * reference, and the trust gauge.
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
     * requires set_oracle_params to have set the
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
    /// Deploy-time extrapolation accounting: blocks with a component larger than
    /// the trained hull (clamped at deploy), and the largest such size seen.
    long n_extrap() const { return n_extrap_; }
    int max_extrap_size() const { return max_extrap_size_; }
    const arma::vec& addc() const { return addc_; }

private:
    /**
     * Shared worker behind extract_block. need_counts fills the scalar
     * descriptors (cne, bre, maxbd, dens) for the additive saddle;
     * need_adjacency fills the block adjacency and side-membership vectors for
     * the surface deploy and the oracle. The two are independent: the surface
     * hot path takes adjacency without the counts, the additive path the
     * reverse. valid, m, and ncn are always set.
     */
    void extract_block_(const arma::imat& G, int i, int j, bool need_adjacency,
                        bool need_counts, ZRatioBlock& bl) const;
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
                          double& gG) const;

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
     * structures (block-position coordinates). Shared by the surface deploy and
     * the gold reference so both score the identical decomposition. Requires
     * the full block extraction (a_blk, si, sj).
     */
    void decompose_(const ZRatioBlock& bl,
                    std::vector<BlockComponent>& out) const;
    /**
     * Accumulate per-component (S1, S2) — surface above size_min, additive
     * below — from the shared decomposition, recording each component when
     * comps != nullptr.
     */
    void accumulate_surface_moments_(const ZRatioBlock& bl, double& s1_out,
                                     double& s2_out,
                                     std::vector<SurfaceComp>* comps) const;

    /**
     * Hot-path surface logR for the deploy branch of log_zratio. Reads the
     * mediating-block node sets left in the extract scratch (xb_rv_, xb_cn_,
     * xb_sio_, xb_sjo_) after an extract_block_ pass with neither adjacency nor
     * counts, finds the CN clusters and bipartite bridges by union-find over
     * reused position-indexed scratch (no arma adjacency, no per-component heap
     * allocation), forms the canonical component-descriptor multiset, and serves
     * the saddle from surf_cache_ keyed on it. Option B makes a component's
     * moments a function of (family, size, density) alone, so a block's logR
     * depends only on this multiset; the cache therefore adds no approximation
     * beyond what the surface already assumes, and accumulating in canonical
     * order makes a hit bit-identical to a fresh evaluation. Matches the
     * decompose_-based reference path (surface_moments) up to floating-point
     * summation order (canonical vs DFS component accumulation).
     */
    double surface_logr_(const arma::imat& G);
    /**
     * Per-component (S1, S2) served from comp_cache_ keyed on the component
     * descriptor (family, size, e, na, nb): surface_eval above size_min, else
     * the additive per-component kernel. Mirrors accumulate_surface_moments_.
     */
    void comp_moments_(int family, int sz, int e, int na, int nb, double& c1,
                       double& c2);

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
    // Deploy-time surface caches. surf_cache_: block component-descriptor
    // multiset -> logR (skips the tilt-grid saddle on recurring blocks, the
    // surface analogue of the additive count-key cache_). comp_cache_: single
    // component descriptor -> (S1, S2) (skips the log/exp per recurring
    // component on cache misses). Both cleared when a new surface is attached.
    // Hashed on the raw key ints (FNV-1a); the keys are exact, so the lookup
    // structure has no effect on the returned values.
    struct IntSeqHash {
        static std::size_t mix(const int* p, std::size_t n) {
            std::uint64_t h = 1469598103934665603ULL;
            for (std::size_t k = 0; k < n; ++k) {
                h ^= static_cast<std::uint32_t>(p[k]);
                h *= 1099511628211ULL;
            }
            return static_cast<std::size_t>(h);
        }
        std::size_t operator()(const std::vector<int>& v) const {
            return mix(v.data(), v.size());
        }
        std::size_t operator()(const std::array<int, 5>& a) const {
            return mix(a.data(), a.size());
        }
    };
    std::unordered_map<std::vector<int>, double, IntSeqHash> surf_cache_;
    std::unordered_map<std::array<int, 5>, std::pair<double, double>, IntSeqHash>
        comp_cache_;
    long n_hit_ = 0, n_miss_ = 0;
    long n_pred_ = 0, n_add_ = 0;
    long n_extrap_ = 0;
    int max_extrap_size_ = 0;

    // Block-Gibbs oracle state (set by set_oracle_params; used by the surface
    // build, the gold reference, and the trust gauge).
    bool oracle_slab_cauchy_ = false;
    double delta_ = 0.0, sigma_ = 1.0, beta_ = 0.5, alpha_ = 1.0;
    SafeRNG* rng_ = nullptr;
    int n_sweep_ = 300, burn_ = 30;

    // Reused scratch for extract_block_. One engine serves one chain, and
    // the counts pass runs once per edge proposal, so per-call heap
    // allocation dominates the extraction cost without these.
    mutable std::vector<char> xb_in_r_;
    mutable std::vector<int> xb_excl_i_, xb_excl_j_, xb_rv_, xb_cn_, xb_sio_,
        xb_sjo_;
    mutable std::vector<unsigned char> xb_side_;

    // Reused scratch for surface_logr_ (one edge proposal at a time). The
    // union-find and accumulator arrays are indexed by block position; once
    // grown to the largest block seen they need no per-edge heap allocation.
    std::vector<std::array<int, 5>> sl_sig_;
    std::vector<int> sl_key_;
    std::vector<int> sl_uf_, sl_sz_, sl_e_, sl_na_, sl_nb_;

    int uf_find_(int x) {
        while (sl_uf_[x] != x) {
            sl_uf_[x] = sl_uf_[sl_uf_[x]];
            x = sl_uf_[x];
        }
        return x;
    }
    void uf_union_(int a, int b) {
        const int ra = uf_find_(a), rb = uf_find_(b);
        if (ra != rb) sl_uf_[ra] = rb;
    }
};

#pragma once

#include <RcppArmadillo.h>
#include <unordered_map>
#include <string>

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

    long cache_size() const { return static_cast<long>(cache_.size()); }
    long n_hit() const { return n_hit_; }
    long n_miss() const { return n_miss_; }
    long n_pred() const { return n_pred_; }
    long n_add() const { return n_add_; }
    long n_clamp() const { return n_clamp_; }

private:
    arma::vec addc_;
    arma::vec tg_, ihat_, ghat_, wt_;
    double psi0_;
    std::unordered_map<std::string, double> cache_;
    long n_hit_ = 0, n_miss_ = 0;
    long n_pred_ = 0, n_add_ = 0, n_clamp_ = 0;
};

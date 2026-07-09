#pragma once

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>

#include "models/ggm/zratio_engine.h"

/**
 * In-chain trust gauge accumulator for the hierarchical-spec per-edge
 * Z-ratio kernel.
 *
 * During K assessment sweeps (deployed selection passes) the gauge compares,
 * at every non-trivial edge move, the deployed log J against the block-local
 * exact reference log R_e and records D = the fraction of accept/reject
 * decisions that would flip under the exact ratio.
 *
 * Per referenced pair, with the actual acceptance log-ratio ln_alpha and
 * sign (+1 add, -1 delete):
 *   s_e      = log J - log R_e                       (log-ratio error, nats)
 *   la_ref   = ln_alpha - sign * s_e                 (exact acceptance ratio)
 *   dalpha   = | min(1,e^ln_alpha) - min(1,e^la_ref) |
 * Per sweep D_all = (n_ent / nE) * mean(dalpha over referenced pairs); the
 * chain D pools by averaging D_all over the K sweeps.
 * Non-trivial = m >= 2 (covers the additive zone, not only the
 * corrected zone); trivial blocks have an exact ratio and are skipped free.
 * Referenced pairs are capped per sweep; the uncapped remainder is counted
 * (n_ent) and never silently dropped.
 */
struct ZRatioGauge {
    bool active = false;
    int cap = 25;         ///< max referenced pairs per sweep
    int n_draws = 120;    ///< reference block-Gibbs draws per pair
    double nE = 0.0;      ///< q(q-1)/2, set at activation

    // Per-sweep scratch (reset by begin_sweep).
    long ent_sweep = 0;
    long ref_sweep = 0;
    double dal_sum_sweep = 0.0;

    // Pooled across sweeps.
    double sum_D = 0.0;
    int n_sweeps = 0;
    long n_ent = 0;       ///< non-trivial pairs seen
    long n_ref = 0;       ///< pairs actually referenced
    long n_capped = 0;    ///< cap hits (non-trivial pairs not referenced)
    double se_sum = 0.0, se_sum2 = 0.0;
    double noise_sum = 0.0;   ///< sum of per-pair dalpha noise (reference MCSE)

    void reset() {
        ent_sweep = ref_sweep = 0;
        dal_sum_sweep = 0.0;
        sum_D = 0.0;
        n_sweeps = 0;
        n_ent = n_ref = n_capped = 0;
        se_sum = se_sum2 = noise_sum = 0.0;
    }
    void begin_sweep() { ent_sweep = 0; ref_sweep = 0; dal_sum_sweep = 0.0; }
    void end_sweep() {
        const double d_all = (nE > 0.0 && ref_sweep > 0)
            ? (static_cast<double>(ent_sweep) / nE) *
                  (dal_sum_sweep / static_cast<double>(ref_sweep))
            : 0.0;
        sum_D += d_all;
        ++n_sweeps;
    }

    double D() const { return n_sweeps > 0 ? sum_D / n_sweeps : 0.0; }
    double se_mean() const { return n_ref > 0 ? se_sum / n_ref : NA_REAL; }
    double se_sd() const {
        if (n_ref < 2) return NA_REAL;
        const double m = se_sum / n_ref;
        const double v = se_sum2 / n_ref - m * m;
        return v > 0.0 ? std::sqrt(v) : 0.0;
    }
    /// Reference-noise floor on the D scale (same scaling as D()).
    double noise_floor() const {
        return (n_sweeps > 0 && nE > 0.0 && n_ref > 0)
            ? (static_cast<double>(n_ent) /
               (static_cast<double>(n_sweeps) * nE)) *
                  (noise_sum / static_cast<double>(n_ref))
            : 0.0;
    }
};

/**
 * Reference one non-trivial edge move and fold it into the gauge. No-op when
 * the gauge is inactive, the engine is absent, or the block is trivial
 * (m < 2, exact ratio). Reads the current graph G (the model's live
 * indicators or continuous subgraph).
 */
inline void zratio_gauge_record(ZRatioGauge& g, ZRatioEngine* engine,
                                const arma::imat& G, int i, int j,
                                double ln_alpha, double log_j, int sign) {
    if (!g.active || engine == nullptr) return;
    ZRatioBlock bl = engine->extract_block(G, i, j);
    if (!bl.valid || bl.m < 2) return;   // trivial: exact ratio, s_e = 0
    ++g.ent_sweep;
    ++g.n_ent;
    if (g.ref_sweep >= g.cap) {
        ++g.n_capped;
        return;
    }
    double log_r = 0.0, mcse = 0.0;
    if (!engine->block_reference_logR(bl.a_blk, bl.si, bl.sj, g.n_draws, log_r,
                                      mcse)) {
        return;
    }
    ++g.ref_sweep;
    ++g.n_ref;
    const double se = log_j - log_r;
    const double la_ref = ln_alpha - sign * se;
    const double a_hat = std::min(1.0, std::exp(ln_alpha));
    const double a_ref = std::min(1.0, std::exp(la_ref));
    const double dalpha = std::abs(a_hat - a_ref);
    const double dnoise =
        std::abs(std::min(1.0, std::exp(ln_alpha + mcse)) -
                 std::min(1.0, std::exp(ln_alpha - mcse))) /
        2.0;
    g.dal_sum_sweep += dalpha;
    g.se_sum += se;
    g.se_sum2 += se * se;
    g.noise_sum += dnoise;
}

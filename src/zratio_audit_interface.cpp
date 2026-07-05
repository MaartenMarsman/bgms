// -----------------------------------------------------------------------------
// zratio_audit_interface.cpp
//
// R-facing entries for the post-sampling Z-ratio alarm suite (hierarchical
// spec). The scan enumerates mediating-block descriptors and the deployed
// correction over one visited graph; the audit scores picked edges against
// the block-Gibbs local oracle. Both are measurement-only: nothing feeds
// back into a fit, a cache, or a chain. Reference: SV/z_graph_prior
// w5c0b_alarm_settlement.R (alarm scan v2, A3 audit) and
// w5c2_failure_corpus.R (audit_v3 incl. the additive-zone channel).
// -----------------------------------------------------------------------------

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "models/ggm/zratio_engine.h"
#include "rng/rng_utils.h"

// -----------------------------------------------------------------------------
// zratio_scan_graph:
//   Mediating-block descriptors for every upper-triangle pair of G whose
//   block is two-sided. Returns a matrix with columns (i, j, ncn, cne,
//   bre, maxbd, m, dens, pred, clamped); i, j are 1-based. `pred` is the
//   deployed correction under `addc` (OLS value inside the hull, 0 outside
//   it or without a packed fit) and `clamped` flags the hull fallback.
// -----------------------------------------------------------------------------

// [[Rcpp::export(name = "zratio_scan_graph")]]
arma::mat zratio_scan_graph(
    arma::imat G,
    arma::vec addc,
    arma::vec tg,
    arma::vec ihat,
    arma::vec ghat,
    arma::vec wt,
    double psi0
) {
    ZRatioEngine engine(addc, tg, ihat, ghat, wt, psi0);
    const int q = static_cast<int>(G.n_rows);
    std::vector<double> rows;
    rows.reserve(static_cast<size_t>(q) * (q - 1) * 5);
    for (int i = 0; i < q - 1; i++) {
        for (int j = i + 1; j < q; j++) {
            ZRatioBlock bl = engine.extract_block(G, i, j);
            if (!bl.valid) continue;
            bool clamped = false;
            double pred = engine.deployed_correction(bl, clamped);
            rows.push_back(i + 1.0);
            rows.push_back(j + 1.0);
            rows.push_back(bl.ncn);
            rows.push_back(bl.cne);
            rows.push_back(bl.bre);
            rows.push_back(bl.maxbd);
            rows.push_back(bl.m);
            rows.push_back(bl.dens);
            rows.push_back(pred);
            rows.push_back(clamped ? 1.0 : 0.0);
        }
    }
    const arma::uword ncol = 10;
    arma::mat out(rows.size() / ncol, ncol);
    for (arma::uword r = 0; r < out.n_rows; r++) {
        for (arma::uword c = 0; c < ncol; c++) {
            out(r, c) = rows[r * ncol + c];
        }
    }
    return out;
}

// -----------------------------------------------------------------------------
// zratio_audit_edges:
//   Measurement-only oracle audit of picked edges on one graph. Per edge:
//   the deployed correction under `addc` versus the block-Gibbs local
//   oracle discrepancy log(saddle on oracle moments) - log(additive), and
//   the audit score |pred - oracle|. (delta, sigma, beta) are the
//   bare-scale prior constants. Rows with ok = 0 had a one-sided block,
//   non-positive additive moments, or an oracle with no finite sweep;
//   their scores are NA.
// -----------------------------------------------------------------------------

// [[Rcpp::export(name = "zratio_audit_edges")]]
Rcpp::List zratio_audit_edges(
    arma::imat G,
    arma::imat edges,
    arma::vec addc,
    arma::vec tg,
    arma::vec ihat,
    arma::vec ghat,
    arma::vec wt,
    double psi0,
    double delta,
    double sigma,
    double beta,
    int n_sweep,
    int burn,
    int seed
) {
    ZRatioEngine engine(addc, tg, ihat, ghat, wt, psi0);
    SafeRNG rng(seed);
    engine.set_oracle_params(delta, sigma, beta, &rng, n_sweep, burn);
    const arma::uword n = edges.n_rows;
    arma::vec err(n), pred(n), oracle(n);
    arma::ivec ok(n);
    err.fill(arma::datum::nan);
    pred.fill(arma::datum::nan);
    oracle.fill(arma::datum::nan);
    for (arma::uword e = 0; e < n; e++) {
        double p_out = 0, o_out = 0;
        ZRatioBlock bl;
        bool good = engine.audit_edge(G, edges(e, 0) - 1, edges(e, 1) - 1,
                                      p_out, o_out, bl);
        ok[e] = good ? 1 : 0;
        if (good) {
            pred[e] = p_out;
            oracle[e] = o_out;
            err[e] = std::abs(p_out - o_out);
        }
    }
    return Rcpp::List::create(
        Rcpp::_["err"] = err,
        Rcpp::_["pred"] = pred,
        Rcpp::_["oracle"] = oracle,
        Rcpp::_["ok"] = ok
    );
}

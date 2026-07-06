// -----------------------------------------------------------------------------
// zratio_test_interface.cpp
//
// R-facing test entries for the hierarchical-spec per-edge Z-ratio engine.
// Kept separate from the sampler dispatch: these drive ZRatioEngine in
// isolation so its saddle map, neighbourhood counts, correction, clamp, and
// cache can be validated against the reference implementation (SV/Z
// sbc_prior_chain_exact.cpp) without any chain machinery.
// -----------------------------------------------------------------------------

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "models/ggm/zratio_engine.h"

// -----------------------------------------------------------------------------
// zratio_test_eval:
//   Evaluate log J = log( Z(G-)/Z(G+) ) for one or more edges on a fixed
//   graph through a shared engine (and hence a shared cache). `edges` is an
//   n x 2 matrix of 1-based (i, j) pairs. Returns the vector of log ratios
//   plus cache/counter diagnostics.
// -----------------------------------------------------------------------------

// [[Rcpp::export(name = "zratio_test_eval")]]
Rcpp::List zratio_test_eval(
    arma::imat G,
    arma::imat edges,
    arma::vec addc,
    arma::vec tg,
    arma::vec ihat,
    arma::vec ghat,
    arma::vec wt,
    double psi0
) {
    ZRatioEngine engine(addc, tg, ihat, ghat, wt, psi0);
    arma::vec out(edges.n_rows);
    for (arma::uword e = 0; e < edges.n_rows; ++e) {
        out[e] = engine.log_zratio(G, edges(e, 0) - 1, edges(e, 1) - 1);
    }
    return Rcpp::List::create(
        Rcpp::_["log_zratio"] = out,
        Rcpp::_["cache_size"] = engine.cache_size(),
        Rcpp::_["n_hit"] = engine.n_hit(),
        Rcpp::_["n_miss"] = engine.n_miss(),
        Rcpp::_["n_pred"] = engine.n_pred(),
        Rcpp::_["n_add"] = engine.n_add(),
        Rcpp::_["n_clamp"] = engine.n_clamp()
    );
}

// -----------------------------------------------------------------------------
// zratio_test_saddle:
//   The bare two-moment saddle map at (s1, s2), for closed-form checks.
// -----------------------------------------------------------------------------

// [[Rcpp::export(name = "zratio_test_saddle")]]
double zratio_test_saddle(
    double s1,
    double s2,
    arma::vec addc,
    arma::vec tg,
    arma::vec ihat,
    arma::vec ghat,
    arma::vec wt,
    double psi0
) {
    ZRatioEngine engine(addc, tg, ihat, ghat, wt, psi0);
    return engine.saddle_ratio(s1, s2);
}

// -----------------------------------------------------------------------------
// zratio_test_calibrated_eval:
//   Drive the online calibrator in isolation: evaluate the (graph, edge)
//   stream in order with calibration enabled, freezing after
//   `freeze_after` evaluations (0 = never freeze). `graphs` is a list of
//   q x q integer matrices, one per row of `edges`. Returns the log
//   ratios, the calibration counters, and the packed post-freeze addc.
//   With a fresh engine, a large n_sweep, and a single edge this doubles
//   as the block-oracle "truth" for that edge.
// -----------------------------------------------------------------------------

// [[Rcpp::export(name = "zratio_test_calibrated_eval")]]
Rcpp::List zratio_test_calibrated_eval(
    Rcpp::List graphs,
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
    int seed,
    int n_sweep,
    int burn,
    int freeze_after,
    bool slab_cauchy = false
) {
    ZRatioEngine engine(addc, tg, ihat, ghat, wt, psi0);
    SafeRNG rng(seed);
    engine.enable_calibration(delta, sigma, beta, &rng, n_sweep, burn, 9.0, 6,
                              slab_cauchy);
    arma::vec out(edges.n_rows);
    for (arma::uword e = 0; e < edges.n_rows; ++e) {
        if (freeze_after > 0 &&
            e == static_cast<arma::uword>(freeze_after)) {
            engine.freeze_calibration();
        }
        arma::imat G = Rcpp::as<arma::imat>(graphs[e]);
        out[e] = engine.log_zratio(G, edges(e, 0) - 1, edges(e, 1) - 1);
    }
    return Rcpp::List::create(
        Rcpp::_["log_zratio"] = out,
        Rcpp::_["n_oracle"] = engine.n_oracle(),
        Rcpp::_["n_anchors"] = engine.n_anchors(),
        Rcpp::_["n_pred"] = engine.n_pred(),
        Rcpp::_["n_add"] = engine.n_add(),
        Rcpp::_["n_clamp"] = engine.n_clamp(),
        Rcpp::_["frozen"] = engine.frozen(),
        Rcpp::_["addc"] = engine.addc()
    );
}

// -----------------------------------------------------------------------------
// zratio_test_precompute:
//   Fill the cache over the bounded count box, then evaluate the requested
//   edges; asserts (in tests) that precomputed and lazily computed values
//   agree and that no cache miss occurs after preloading.
// -----------------------------------------------------------------------------

// [[Rcpp::export(name = "zratio_test_precompute")]]
Rcpp::List zratio_test_precompute(
    arma::imat G,
    arma::imat edges,
    arma::vec addc,
    arma::vec tg,
    arma::vec ihat,
    arma::vec ghat,
    arma::vec wt,
    double psi0,
    int ncn_max,
    int bre_max
) {
    ZRatioEngine engine(addc, tg, ihat, ghat, wt, psi0);
    engine.precompute_table(ncn_max, bre_max);
    long preload = engine.cache_size();
    arma::vec out(edges.n_rows);
    for (arma::uword e = 0; e < edges.n_rows; ++e) {
        out[e] = engine.log_zratio(G, edges(e, 0) - 1, edges(e, 1) - 1);
    }
    return Rcpp::List::create(
        Rcpp::_["log_zratio"] = out,
        Rcpp::_["preload"] = preload,
        Rcpp::_["n_miss"] = engine.n_miss()
    );
}

// -----------------------------------------------------------------------------
// zratio_test_interface.cpp
//
// R-facing test entries for the hierarchical-spec per-edge Z-ratio engine.
// Kept separate from the sampler dispatch: these drive ZRatioEngine in
// isolation so its saddle map, neighbourhood counts, correction, clamp, and
// cache can be tested without any chain machinery.
// -----------------------------------------------------------------------------

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "models/ggm/zratio_engine.h"
#include "models/ggm/zratio_law.h"

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
        Rcpp::_["n_add"] = engine.n_add()
    );
}

// -----------------------------------------------------------------------------
// zratio_test_reference:
//   Block-local EXACT reference log R_e for the edge (i, j) (1-based) on G,
//   via block_reference_logR (full-product endpoint transform). Returns the
//   reference, its batch-means MCSE, the deployed log J for the same edge,
//   and block descriptors. Drives the trust-gauge reference in isolation.
// -----------------------------------------------------------------------------

// [[Rcpp::export(name = "zratio_test_reference")]]
Rcpp::List zratio_test_reference(
    arma::imat G,
    int i,
    int j,
    arma::vec addc,
    arma::vec tg,
    arma::vec ihat,
    arma::vec ghat,
    arma::vec wt,
    double psi0,
    double delta,
    double eta,
    int n_draws,
    int burn,
    int seed,
    bool slab_cauchy = false,
    double alpha = 1.0
) {
    ZRatioEngine engine(addc, tg, ihat, ghat, wt, psi0);
    SafeRNG rng(seed);
    engine.set_oracle_params(delta, eta, &rng, n_draws, burn, slab_cauchy,
                             alpha);
    ZRatioBlock bl = engine.extract_block(G, i - 1, j - 1);
    if (!bl.valid) {
        return Rcpp::List::create(Rcpp::_["valid"] = false);
    }
    double logR = NA_REAL, mcse = NA_REAL;
    bool ok =
        engine.block_reference_logR(bl.a_blk, bl.si, bl.sj, n_draws, logR, mcse);
    return Rcpp::List::create(
        Rcpp::_["valid"] = true,
        Rcpp::_["ok"] = ok,
        Rcpp::_["logR"] = logR,
        Rcpp::_["mcse"] = mcse,
        Rcpp::_["m"] = bl.m,
        Rcpp::_["maxbd"] = bl.maxbd,
        Rcpp::_["ncn"] = bl.ncn,
        Rcpp::_["log_zratio"] = engine.log_zratio(G, i - 1, j - 1)
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
// zratio_test_precompute:
//   Fill the cache over the bounded count box, then evaluate the requested
//   edges; asserts (in tests) that precomputed and lazily computed values
//   agree and that no cache miss occurs after preloading.
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// zratio_test_surface_eval:
//   Drive the Option-B surface deploy for one edge (i, j) (1-based) on G. The
//   `surface` list carries cn/bip sublists (c1, c2 = 9 raw-poly coeffs; the
//   size/dens/log-moment hulls; size_min). Returns the summed (S1, S2), the
//   surface logR, the hot-path log_zratio (which routes through the surface
//   branch under the alpha = 1 Normal cell), and the per-component decomposition
//   for validation against the R deploy_surface.
// -----------------------------------------------------------------------------

// [[Rcpp::export(name = "zratio_test_surface_eval")]]
Rcpp::List zratio_test_surface_eval(
    arma::imat G,
    int i,
    int j,
    arma::vec addc,
    arma::vec tg,
    arma::vec ihat,
    arma::vec ghat,
    arma::vec wt,
    double psi0,
    Rcpp::List surface,
    double delta,
    double eta,
    bool slab_cauchy = false,
    double alpha = 1.0
) {
    ZRatioEngine engine(addc, tg, ihat, ghat, wt, psi0);
    engine.set_oracle_params(delta, eta, nullptr, 300, 30, slab_cauchy, alpha);
    engine.set_surface(surface_family_from_list(surface["cn"]),
                       surface_family_from_list(surface["bip"]));
    double s1 = NA_REAL, s2 = NA_REAL, logr = NA_REAL;
    std::vector<SurfaceComp> comps;
    bool valid = engine.surface_moments(G, i - 1, j - 1, s1, s2, logr, comps);
    double lz = engine.log_zratio(G, i - 1, j - 1);
    const int nc = static_cast<int>(comps.size());
    Rcpp::IntegerVector fam(nc), sz(nc), ee(nc), na(nc), nb(nc), used(nc);
    Rcpp::NumericVector dens(nc), cs1(nc), cs2(nc);
    for (int k = 0; k < nc; ++k) {
        fam[k] = comps[k].family;
        sz[k] = comps[k].size;
        ee[k] = comps[k].e;
        na[k] = comps[k].na;
        nb[k] = comps[k].nb;
        used[k] = comps[k].used_surface ? 1 : 0;
        dens[k] = comps[k].dens;
        cs1[k] = comps[k].s1;
        cs2[k] = comps[k].s2;
    }
    return Rcpp::List::create(
        Rcpp::_["valid"] = valid,
        Rcpp::_["has_surface"] = engine.has_surface(),
        Rcpp::_["log_zratio"] = lz,
        Rcpp::_["logR"] = logr,
        Rcpp::_["S1"] = s1,
        Rcpp::_["S2"] = s2,
        Rcpp::_["comp"] = Rcpp::DataFrame::create(
            Rcpp::_["family"] = fam, Rcpp::_["size"] = sz, Rcpp::_["e"] = ee,
            Rcpp::_["na"] = na, Rcpp::_["nb"] = nb, Rcpp::_["dens"] = dens,
            Rcpp::_["s1"] = cs1, Rcpp::_["s2"] = cs2,
            Rcpp::_["used_surface"] = used));
}

// -----------------------------------------------------------------------------
// zratio_test_surface_batch:
//   Set the Option-B surface once, then evaluate log_zratio for every edge in
//   `edges` (1-based) on a SHARED engine, so the deploy-time surface cache
//   persists across the sweep exactly as in the sampler. Returns the logR
//   vector and the prediction/cache diagnostics, for parity and timing checks.
// -----------------------------------------------------------------------------

// [[Rcpp::export(name = "zratio_test_surface_batch")]]
Rcpp::List zratio_test_surface_batch(
    arma::imat G,
    arma::imat edges,
    arma::vec addc,
    arma::vec tg,
    arma::vec ihat,
    arma::vec ghat,
    arma::vec wt,
    double psi0,
    Rcpp::List surface,
    double delta,
    double eta,
    bool slab_cauchy = false,
    double alpha = 1.0
) {
    ZRatioEngine engine(addc, tg, ihat, ghat, wt, psi0);
    engine.set_oracle_params(delta, eta, nullptr, 300, 30, slab_cauchy, alpha);
    engine.set_surface(surface_family_from_list(surface["cn"]),
                       surface_family_from_list(surface["bip"]));
    arma::vec out(edges.n_rows);
    for (arma::uword e = 0; e < edges.n_rows; ++e) {
        out[e] = engine.log_zratio(G, edges(e, 0) - 1, edges(e, 1) - 1);
    }
    return Rcpp::List::create(
        Rcpp::_["log_zratio"] = out,
        Rcpp::_["n_pred"] = engine.n_pred(),
        Rcpp::_["n_add"] = engine.n_add(),
        Rcpp::_["n_extrap"] = engine.n_extrap(),
        Rcpp::_["max_extrap_size"] = engine.max_extrap_size()
    );
}

// -----------------------------------------------------------------------------
// zratio_test_gold_moments:
//   Per-component gold reference logR for the edge (i, j) (1-based) on G: the
//   block-Gibbs oracle on each non-trivial component's own sub-adjacency,
//   summed through the same saddle closure the surface uses. Drives the Stage-5
//   gold edge-toggle comparison against the deployed surface and the additive
//   baseline.
// -----------------------------------------------------------------------------

// [[Rcpp::export(name = "zratio_test_gold_moments")]]
Rcpp::List zratio_test_gold_moments(
    arma::imat G,
    int i,
    int j,
    arma::vec addc,
    arma::vec tg,
    arma::vec ihat,
    arma::vec ghat,
    arma::vec wt,
    double psi0,
    double delta,
    double eta,
    int n_sweep,
    int burn,
    int seed,
    bool slab_cauchy = false,
    double alpha = 1.0
) {
    ZRatioEngine engine(addc, tg, ihat, ghat, wt, psi0);
    SafeRNG rng(seed);
    engine.set_oracle_params(delta, eta, &rng, n_sweep, burn, slab_cauchy,
                             alpha);
    double s1 = NA_REAL, s2 = NA_REAL, logr = NA_REAL;
    bool valid = engine.gold_moments(G, i - 1, j - 1, s1, s2, logr);
    return Rcpp::List::create(
        Rcpp::_["valid"] = valid,
        Rcpp::_["logR"] = logr,
        Rcpp::_["S1"] = s1,
        Rcpp::_["S2"] = s2
    );
}

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

// -----------------------------------------------------------------------------
// zratio_law_moments:
//   Test-only entry for the DORMANT mu-first CPA analytic CN law (zratio_law.h,
//   not wired into the default build; see that header). Builds the (eta, delta)
//   grids, runs the self-consistent solve for a component of `n` common
//   neighbours at mean degree D, gates on the psi residual, and dresses to
//   (S1, S2). `certified` is FALSE when the solve does not converge within
//   `budget` or the dressing is non-positive. Deterministic (no RNG). Drives
//   the law-engine validation in test-zratio-law.R.
// -----------------------------------------------------------------------------

// [[Rcpp::export(name = "zratio_law_moments")]]
Rcpp::List zratio_law_moments(
    double eta,
    double delta,
    double D,
    double n,
    int budget = 3000
) {
    double s1 = NA_REAL, s2 = NA_REAL, psi_gap = NA_REAL;
    bool certified =
        zratio_law::law_moments(eta, delta, D, n, budget, s1, s2, psi_gap);
    return Rcpp::List::create(
        Rcpp::_["certified"] = certified,
        Rcpp::_["S1"] = s1,
        Rcpp::_["S2"] = s2,
        Rcpp::_["psi_gap"] = psi_gap
    );
}

// -----------------------------------------------------------------------------
// zratio_interface.cpp
//
// Production R-facing entry for the hierarchical-spec per-edge Z-ratio engine.
// Kept apart from the zratio_test_* entries (zratio_test_interface.cpp), which
// drive the engine in isolation for tests: this file holds the entry the fit
// itself calls.
// -----------------------------------------------------------------------------

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "models/ggm/zratio_engine.h"

// -----------------------------------------------------------------------------
// zratio_block_oracle_moments:
//   Weighted block moments (S1, S2) from the block-Gibbs oracle on a BARE
//   mediating component, standing alone (no host graph). This is the offline
//   anchor generator for the Option-B absolute-moment surfaces (called on every
//   alpha = 1 hierarchical fit by R/zratio_surfaces.R): for a CN cluster every
//   node is a common neighbour, so si = sj = all block rows; for a bipartite
//   bridge structure si = side-A rows, sj = side-B rows (both 0-based into
//   `a_blk`). Reuses the bgms-aligned Schur+rank-2-SMW kernel, so the anchors
//   are drawn from the same law the sampler uses. addc/tg/ihat/ghat/wt/psi0 are
//   unused by the moment recipe (it reads only delta, eta, sigma = 1) but the
//   engine constructor requires them; pass the fit-time constants for the cell.
// -----------------------------------------------------------------------------

// [[Rcpp::export(name = "zratio_block_oracle_moments")]]
Rcpp::List zratio_block_oracle_moments(
    arma::imat a_blk,
    arma::uvec si,
    arma::uvec sj,
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
    double s1 = NA_REAL, s2 = NA_REAL;
    bool ok = engine.block_oracle_moments(a_blk, si, sj, s1, s2);
    return Rcpp::List::create(
        Rcpp::_["ok"] = ok,
        Rcpp::_["S1"] = s1,
        Rcpp::_["S2"] = s2
    );
}

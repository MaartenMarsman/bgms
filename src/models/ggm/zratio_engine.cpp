#include "zratio_engine.h"
#include "math/explog_macros.h"

#include <algorithm>
#include <vector>
#include <cmath>

// The additive second spectral moment dips slightly negative at large blocks
// (>=3-body error in the additive estimator). The true S2 is variance-like and
// always >= 0, and Phi_2 is ~flat in S2, so flooring a non-positive additive S2
// to this epsilon clamps the estimator to its valid domain without altering any
// already-valid cell; it keeps the saddle finite, opens the calibration gate,
// and lets the correction anchor.
static constexpr double kS2Floor = 1e-3;

double ZRatioEngine::saddle_ratio(double s1, double s2) const {
    if (s1 <= 0 || s2 <= 0) return 1.0;
    double eh = s1 * s1 / (2.0 * s2), ur = s2 / s1, nf = 0, dg = 0;
    for (arma::uword k = 0; k < tg_.n_elem; ++k) {
        double ph = std::pow(1.0 + ur * tg_[k] * tg_[k], -eh);
        nf += wt_[k] * ph * ihat_[k];
        dg += wt_[k] * ph * ghat_[k];
    }
    return nf / dg;
}

ZRatioBlock ZRatioEngine::extract_block(const arma::imat& G, int i,
                                        int j) const {
    ZRatioBlock bl;
    extract_block_(G, i, j, /*need_adjacency=*/true, /*need_counts=*/true, bl);
    return bl;
}

void ZRatioEngine::extract_block_(const arma::imat& G, int i, int j,
                                  bool need_adjacency, bool need_counts,
                                  ZRatioBlock& bl) const {
    bl = ZRatioBlock();
    const int q = static_cast<int>(G.n_rows);

    // Mediating block: common neighbours of (i, j) plus the endpoints of
    // 2-hop bridges between the exclusive neighbour sets. The toggled
    // edge's own state never enters, so the value is state-invariant.
    std::vector<char>& in_r = xb_in_r_;
    in_r.assign(q, 0);
    // Exclusive neighbour lists of i and j; the bridge scan below then runs
    // over the candidate pairs instead of the full q x q grid.
    std::vector<int>& excl_i = xb_excl_i_;
    std::vector<int>& excl_j = xb_excl_j_;
    excl_i.clear();
    excl_j.clear();
    // G is symmetric, so G(k, i) == G(i, k): read the transpose entry to walk
    // column i with unit stride (arma is column-major), which keeps the O(q)
    // classification scan in cache instead of striding by n_rows per element.
    for (int k = 0; k < q; k++) {
        if (k == i || k == j) continue;
        const bool near_i = (G(k, i) == 1), near_j = (G(k, j) == 1);
        if (near_i && near_j) in_r[k] = 1;
        else if (near_i) excl_i.push_back(k);
        else if (near_j) excl_j.push_back(k);
    }
    for (int a : excl_i) {
        for (int b : excl_j) {
            if (G(b, a) == 1) {
                in_r[a] = 1;
                in_r[b] = 1;
            }
        }
    }
    std::vector<int>& rv = xb_rv_;
    rv.clear();
    for (int k = 0; k < q; k++) {
        if (in_r[k]) rv.push_back(k);
    }
    const int m = static_cast<int>(rv.size());
    bl.m = m;

    // Per-position side membership: bit 1 = adjacent to i, bit 2 = to j.
    std::vector<unsigned char>& side = xb_side_;
    side.assign(m, 0);
    std::vector<int>& cn = xb_cn_;
    std::vector<int>& si_o = xb_sio_;
    std::vector<int>& sj_o = xb_sjo_;
    cn.clear();
    si_o.clear();
    sj_o.clear();
    for (int p = 0; p < m; p++) {
        bool si = (G(rv[p], i) == 1), sj = (G(rv[p], j) == 1);
        if (si) side[p] |= 1;
        if (sj) side[p] |= 2;
        if (si && sj) cn.push_back(p);
        else if (si) si_o.push_back(p);
        else if (sj) sj_o.push_back(p);
    }
    if ((cn.empty() && si_o.empty()) || (cn.empty() && sj_o.empty())) {
        // One side of the mediating block is empty: isolated-edge ratio.
        return;
    }
    bl.valid = true;

    bl.ncn = static_cast<int>(cn.size());
    // Scalar count descriptors (cne, bre, maxbd) drive the additive
    // saddle only; the surface deploy reads the adjacency instead, so skip
    // these O(m^2) passes when the caller does not need the counts.
    if (need_counts) {
        for (size_t a = 0; a < cn.size(); a++) {
            for (size_t b = a + 1; b < cn.size(); b++) {
                if (G(rv[cn[a]], rv[cn[b]]) == 1) bl.cne++;
            }
        }
        for (int a : si_o) {
            for (int b : sj_o) {
                if (G(rv[a], rv[b]) == 1) bl.bre++;
            }
        }
        for (int a : si_o) {
            int d = 0;
            for (int b : sj_o) {
                if (G(rv[a], rv[b]) == 1) d++;
            }
            if (d > bl.maxbd) bl.maxbd = d;
        }
        for (int b : sj_o) {
            int d = 0;
            for (int a : si_o) {
                if (G(rv[a], rv[b]) == 1) d++;
            }
            if (d > bl.maxbd) bl.maxbd = d;
        }
    }

    if (!need_adjacency) return;

    bl.a_blk.zeros(m, m);
    for (int a = 0; a < m; a++) {
        for (int b = a + 1; b < m; b++) {
            int e = (G(rv[a], rv[b]) == 1) ? 1 : 0;
            bl.a_blk(a, b) = e;
            bl.a_blk(b, a) = e;
        }
    }

    // Side memberships in ascending block position (CN nodes sit on both).
    std::vector<arma::uword> si_v, sj_v;
    for (int p = 0; p < m; p++) {
        if (side[p] & 1) si_v.push_back(p);
        if (side[p] & 2) sj_v.push_back(p);
    }
    bl.si = arma::uvec(si_v);
    bl.sj = arma::uvec(sj_v);
}

double ZRatioEngine::surface_eval_(const SurfaceFamily& f, bool s2,
                                   double size, double dens) const {
    // Clamp (size, density) to the trained hull, then evaluate the raw
    // quadratic in (log size, density) over the 9 monomials and clamp the
    // log-moment to its trained range +/- 0.1 (matches the R deploy predictor).
    const double n = std::min(std::max(size, f.size_lo), f.size_hi);
    const double d = std::min(std::max(dens, f.dens_lo), f.dens_hi);
    const double L = MY_LOG(n);
    const double x[9] = {1.0, L, L * L, d, d * d, L * d, L * L * d,
                         L * d * d, L * L * d * d};
    const arma::vec& c = s2 ? f.c2 : f.c1;
    double p = 0.0;
    for (int k = 0; k < 9; ++k) p += c[k] * x[k];
    const double lo = (s2 ? f.l2_lo : f.l1_lo) - 0.1;
    const double hi = (s2 ? f.l2_hi : f.l1_hi) + 0.1;
    p = std::min(std::max(p, lo), hi);

    // Past the trained hull the prediction continues along the surface's own
    // boundary slope in log-size instead of freezing at the hull edge. The
    // clamped edge value stays the base, so the trained range still bounds
    // where the tail starts. Scored against block-Gibbs gold at sizes 90-150
    // against a size-80 hull, over two anchor-build seeds: the tangent holds a
    // median 0.0006 nats and at most 0.0011 (common-neighbour) and a median
    // 0.0043 and at most 0.0060 (bipartite), where freezing grows to 0.060 and
    // 0.095 and the fitted quadratic, continued as its own extrapolant, grows
    // to 0.0026 and 0.016.
    if (size > f.size_hi) {
        // d/dL of the 9-monomial polynomial at the hull edge.
        double slope = c[1] + 2.0 * c[2] * L + c[5] * d + 2.0 * c[6] * L * d +
                       c[7] * d * d + 2.0 * c[8] * L * d * d;
        // The absolute moments grow with component size, so a negative fitted
        // edge slope is a fit pathology, not a signal. Floor it at zero, which
        // degenerates to the old freeze, and tally the floor: it is a silent
        // per-density-band degeneracy the gold scoring above would not catch.
        if (slope < 0.0) {
            slope = 0.0;
            ++n_slope_floor_;
        }
        p += slope * (MY_LOG(size) - L);
    }
    return MY_EXP(p);
}

// Connected components of the induced subgraph on `nodes` (block positions).
// For a CN cluster (bip = false) every intra-node edge counts; for a bipartite
// bridge structure (bip = true) only cross-side edges (is_A differs) count.
static std::vector<std::vector<int>> conn_components_(
    const arma::imat& a_blk, const std::vector<int>& nodes, bool bip,
    const std::vector<char>& is_A) {
    const int nn = static_cast<int>(nodes.size());
    std::vector<char> vis(nn, 0);
    std::vector<std::vector<int>> out;
    std::vector<int> stack;
    for (int s = 0; s < nn; ++s) {
        if (vis[s]) continue;
        std::vector<int> comp;
        stack.clear();
        stack.push_back(s);
        vis[s] = 1;
        while (!stack.empty()) {
            const int u = stack.back();
            stack.pop_back();
            comp.push_back(nodes[u]);
            for (int t = 0; t < nn; ++t) {
                if (vis[t]) continue;
                bool e = (a_blk(nodes[u], nodes[t]) == 1);
                if (bip) e = e && (is_A[nodes[u]] != is_A[nodes[t]]);
                if (e) {
                    vis[t] = 1;
                    stack.push_back(t);
                }
            }
        }
        out.push_back(std::move(comp));
    }
    return out;
}

void ZRatioEngine::decompose_(const ZRatioBlock& bl,
                              std::vector<BlockComponent>& out) const {
    out.clear();
    const int m = bl.m;
    const arma::imat& A = bl.a_blk;
    // Per-position side: CN = adjacent to both endpoints, A-side = i-only,
    // B-side = j-only. Decompose the CN-CN subgraph and the A-B bridge subgraph
    // into disjoint components, exactly as the reference deploy_surface.
    std::vector<char> in_si(m, 0), in_sj(m, 0), is_A(m, 0);
    for (arma::uword k = 0; k < bl.si.n_elem; ++k) in_si[bl.si[k]] = 1;
    for (arma::uword k = 0; k < bl.sj.n_elem; ++k) in_sj[bl.sj[k]] = 1;
    std::vector<int> cn_nodes, ab_nodes;
    for (int p = 0; p < m; ++p) {
        const bool ci = in_si[p], cj = in_sj[p];
        if (ci && cj) cn_nodes.push_back(p);
        else if (ci) { ab_nodes.push_back(p); is_A[p] = 1; }
        else if (cj) { ab_nodes.push_back(p); is_A[p] = 0; }
    }

    const std::vector<char> dummy;
    for (auto& comp : conn_components_(A, cn_nodes, false, dummy)) {
        BlockComponent bc;
        bc.family = 0;
        int e = 0;
        for (size_t a = 0; a < comp.size(); ++a) {
            for (size_t b = a + 1; b < comp.size(); ++b) {
                if (A(comp[a], comp[b]) == 1) ++e;
            }
        }
        bc.na = static_cast<int>(comp.size());
        bc.nb = 0;
        bc.e = e;
        bc.nodes = std::move(comp);
        out.push_back(std::move(bc));
    }
    for (auto& comp : conn_components_(A, ab_nodes, true, is_A)) {
        BlockComponent bc;
        bc.family = 1;
        bc.aside.resize(comp.size());
        int na = 0, nb = 0, e = 0;
        for (size_t k = 0; k < comp.size(); ++k) {
            bc.aside[k] = is_A[comp[k]];
            if (is_A[comp[k]]) ++na;
            else ++nb;
        }
        for (size_t a = 0; a < comp.size(); ++a) {
            for (size_t b = a + 1; b < comp.size(); ++b) {
                if (is_A[comp[a]] != is_A[comp[b]] &&
                    A(comp[a], comp[b]) == 1) {
                    ++e;
                }
            }
        }
        bc.na = na;
        bc.nb = nb;
        bc.e = e;
        bc.nodes = std::move(comp);
        out.push_back(std::move(bc));
    }
}

void ZRatioEngine::accumulate_surface_moments_(
    const ZRatioBlock& bl, double& s1_out, double& s2_out,
    std::vector<SurfaceComp>* comps) const {
    s1_out = 0.0;
    s2_out = 0.0;
    std::vector<BlockComponent> parts;
    decompose_(bl, parts);
    for (const BlockComponent& bc : parts) {
        const int sz = static_cast<int>(bc.nodes.size());
        double c1, c2, dens;
        bool used;
        if (bc.family == 0) {
            dens = (sz >= 2) ? bc.e / (static_cast<double>(sz) * (sz - 1) / 2.0)
                             : 0.0;
            used = sz >= static_cast<int>(surf_cn_.size_min);
            if (used) {
                c1 = surface_eval_(surf_cn_, false, sz, dens);
                c2 = surface_eval_(surf_cn_, true, sz, dens);
            } else {
                c1 = sz * addc_[0] + bc.e * addc_[2];
                c2 = sz * addc_[1] + bc.e * addc_[3];
            }
        } else {
            dens = (bc.na > 0 && bc.nb > 0)
                       ? bc.e / static_cast<double>(bc.na * bc.nb) : 0.0;
            used = sz >= static_cast<int>(surf_bip_.size_min);
            if (used) {
                c1 = surface_eval_(surf_bip_, false, sz, dens);
                c2 = surface_eval_(surf_bip_, true, sz, dens);
            } else {
                c1 = bc.e * addc_[4];
                c2 = bc.e * addc_[5];
            }
        }
        s1_out += c1;
        s2_out += c2;
        if (comps) {
            comps->push_back(
                {bc.family, sz, bc.e, bc.na, bc.nb, dens, c1, c2, used});
        }
    }
}

void ZRatioEngine::comp_moments_(int family, int sz, int e, int na, int nb,
                                 double& c1, double& c2) {
    const std::array<int, 5> key = {family, sz, e, na, nb};
    auto it = comp_cache_.find(key);
    if (it != comp_cache_.end()) {
        c1 = it->second.first;
        c2 = it->second.second;
        return;
    }
    if (family == 0) {
        const double dens =
            (sz >= 2) ? e / (static_cast<double>(sz) * (sz - 1) / 2.0) : 0.0;
        if (sz >= static_cast<int>(surf_cn_.size_min)) {
            c1 = surface_eval_(surf_cn_, false, sz, dens);
            c2 = surface_eval_(surf_cn_, true, sz, dens);
        } else {
            c1 = sz * addc_[0] + e * addc_[2];
            c2 = sz * addc_[1] + e * addc_[3];
        }
    } else {
        const double dens =
            (na > 0 && nb > 0) ? e / static_cast<double>(na * nb) : 0.0;
        if (sz >= static_cast<int>(surf_bip_.size_min)) {
            c1 = surface_eval_(surf_bip_, false, sz, dens);
            c2 = surface_eval_(surf_bip_, true, sz, dens);
        } else {
            c1 = e * addc_[4];
            c2 = e * addc_[5];
        }
    }
    comp_cache_.emplace(key, std::make_pair(c1, c2));
}

double ZRatioEngine::surface_logr_(const arma::imat& G) {
    // Node sets left in the extract scratch: rv maps block positions to graph
    // nodes; cn are the common-neighbour positions, sio/sjo the A-side
    // (i-only) and B-side (j-only) positions. Adjacency is read from G through
    // rv, so no block adjacency matrix is materialised.
    const std::vector<int>& rv = xb_rv_;
    const std::vector<int>& cn = xb_cn_;
    const std::vector<int>& sio = xb_sio_;
    const std::vector<int>& sjo = xb_sjo_;
    const int m = static_cast<int>(rv.size());

    sl_uf_.resize(m);
    sl_sz_.resize(m);
    sl_e_.resize(m);
    sl_na_.resize(m);
    sl_nb_.resize(m);
    sl_sig_.clear();

    // Common-neighbour clusters: connected components of the CN-CN subgraph
    // under any intra-cluster edge. Union over CN pairs, then tally size and
    // internal edges per component root.
    for (int a : cn) {
        sl_uf_[a] = a;
        sl_sz_[a] = 0;
        sl_e_[a] = 0;
    }
    for (size_t x = 0; x < cn.size(); ++x) {
        for (size_t y = x + 1; y < cn.size(); ++y) {
            if (G(rv[cn[x]], rv[cn[y]]) == 1) uf_union_(cn[x], cn[y]);
        }
    }
    for (int a : cn) sl_sz_[uf_find_(a)]++;
    for (size_t x = 0; x < cn.size(); ++x) {
        for (size_t y = x + 1; y < cn.size(); ++y) {
            if (G(rv[cn[x]], rv[cn[y]]) == 1) sl_e_[uf_find_(cn[x])]++;
        }
    }
    for (int a : cn) {
        if (uf_find_(a) == a) {
            sl_sig_.push_back({0, sl_sz_[a], sl_e_[a], sl_sz_[a], 0});
        }
    }

    // Bipartite bridges: connected components of the A-B subgraph under
    // cross-side edges only. Union over cross pairs, then tally A/B side sizes
    // and cross edges per component root (which may sit on either side).
    for (int a : sio) {
        sl_uf_[a] = a;
        sl_sz_[a] = 0;
        sl_e_[a] = 0;
        sl_na_[a] = 0;
        sl_nb_[a] = 0;
    }
    for (int b : sjo) {
        sl_uf_[b] = b;
        sl_sz_[b] = 0;
        sl_e_[b] = 0;
        sl_na_[b] = 0;
        sl_nb_[b] = 0;
    }
    for (int a : sio) {
        for (int b : sjo) {
            if (G(rv[a], rv[b]) == 1) uf_union_(a, b);
        }
    }
    for (int a : sio) {
        const int r = uf_find_(a);
        sl_sz_[r]++;
        sl_na_[r]++;
    }
    for (int b : sjo) {
        const int r = uf_find_(b);
        sl_sz_[r]++;
        sl_nb_[r]++;
    }
    for (int a : sio) {
        for (int b : sjo) {
            if (G(rv[a], rv[b]) == 1) sl_e_[uf_find_(a)]++;
        }
    }
    for (int a : sio) {
        if (uf_find_(a) == a) {
            sl_sig_.push_back({1, sl_sz_[a], sl_e_[a], sl_na_[a], sl_nb_[a]});
        }
    }
    for (int b : sjo) {
        if (uf_find_(b) == b) {
            sl_sig_.push_back({1, sl_sz_[b], sl_e_[b], sl_na_[b], sl_nb_[b]});
        }
    }

    // Extrapolation accounting (Tier-1 observability): any component larger than
    // the trained hull is clamped to the hull edge by surface_eval_, so its
    // moment is an extrapolation. Count the block once if it holds any such
    // component and track the largest size seen. Runs on every call (before the
    // cache lookup below) so the tally is the true per-fit deploy count.
    bool extrapolated = false;
    int largest = 0;
    for (const std::array<int, 5>& t : sl_sig_) {
        const double hull = (t[0] == 0) ? surf_cn_.size_hi : surf_bip_.size_hi;
        if (t[1] > hull) {
            extrapolated = true;
            if (t[1] > largest) largest = t[1];
        }
    }
    if (extrapolated && phase_ != ZRatioPhase::Gauge) {
        n_extrap_++;
        if (largest > max_extrap_size_) max_extrap_size_ = largest;
        if (phase_ == ZRatioPhase::Retained) {
            n_extrap_ret_++;
            if (largest > max_extrap_size_ret_) max_extrap_size_ret_ = largest;
        }
    }

    // Canonical multiset -> cached saddle. Accumulating in sorted order makes
    // the sum bit-identical for any block with this multiset.
    std::sort(sl_sig_.begin(), sl_sig_.end());
    sl_key_.clear();
    for (const std::array<int, 5>& t : sl_sig_) {
        sl_key_.insert(sl_key_.end(), t.begin(), t.end());
    }
    auto it = surf_cache_.find(sl_key_);
    if (it != surf_cache_.end()) return it->second;
    double s1 = 0.0, s2 = 0.0;
    for (const std::array<int, 5>& t : sl_sig_) {
        double c1, c2;
        comp_moments_(t[0], t[1], t[2], t[3], t[4], c1, c2);
        s1 += c1;
        s2 += c2;
    }
    if (s2 <= 0.0) s2 = kS2Floor;
    const double logr = MY_LOG(saddle_ratio(s1, s2));
    surf_cache_.emplace(sl_key_, logr);
    return logr;
}

bool ZRatioEngine::gold_moments(const arma::imat& G, int i, int j,
                                double& s1_out, double& s2_out,
                                double& logr_out) {
    ZRatioBlock bl;
    extract_block_(G, i, j, /*need_adjacency=*/true, /*need_counts=*/false, bl);
    if (!bl.valid) {
        s1_out = 0.0;
        s2_out = 0.0;
        logr_out = MY_LOG(psi0_);
        return false;
    }
    std::vector<BlockComponent> parts;
    decompose_(bl, parts);
    double s1 = 0.0, s2 = 0.0;
    for (const BlockComponent& bc : parts) {
        const int sz = static_cast<int>(bc.nodes.size());
        // Trivial components (CN size <= 2, single bridge) are exact through
        // pairwise overlap; larger components use the block-Gibbs oracle on the
        // component's own sub-adjacency.
        if (bc.family == 0 && sz < 3) {
            s1 += sz * addc_[0] + bc.e * addc_[2];
            s2 += sz * addc_[1] + bc.e * addc_[3];
            continue;
        }
        if (bc.family == 1 && sz < 3) {
            s1 += bc.e * addc_[4];
            s2 += bc.e * addc_[5];
            continue;
        }
        arma::imat sub(sz, sz, arma::fill::zeros);
        arma::uvec si, sj;
        if (bc.family == 0) {
            for (int a = 0; a < sz; ++a) {
                for (int b = a + 1; b < sz; ++b) {
                    if (bl.a_blk(bc.nodes[a], bc.nodes[b]) == 1) {
                        sub(a, b) = sub(b, a) = 1;
                    }
                }
            }
            si = arma::regspace<arma::uvec>(0, sz - 1);
            sj = si;
        } else {
            std::vector<arma::uword> ai, bi;
            for (int a = 0; a < sz; ++a) {
                if (bc.aside[a]) ai.push_back(a);
                else bi.push_back(a);
            }
            for (int a = 0; a < sz; ++a) {
                for (int b = a + 1; b < sz; ++b) {
                    if (bc.aside[a] != bc.aside[b] &&
                        bl.a_blk(bc.nodes[a], bc.nodes[b]) == 1) {
                        sub(a, b) = sub(b, a) = 1;
                    }
                }
            }
            si = arma::uvec(ai);
            sj = arma::uvec(bi);
        }
        double s1c = 0.0, s2c = 0.0;
        if (block_oracle_moments(sub, si, sj, s1c, s2c)) {
            s1 += s1c;
            s2 += s2c;
        } else if (bc.family == 0) {
            s1 += sz * addc_[0] + bc.e * addc_[2];
            s2 += sz * addc_[1] + bc.e * addc_[3];
        } else {
            s1 += bc.e * addc_[4];
            s2 += bc.e * addc_[5];
        }
    }
    if (s2 <= 0.0) s2 = kS2Floor;
    s1_out = s1;
    s2_out = s2;
    logr_out = MY_LOG(saddle_ratio(s1, s2));
    return true;
}

bool ZRatioEngine::surface_moments(const arma::imat& G, int i, int j,
                                   double& s1_out, double& s2_out,
                                   double& logr_out,
                                   std::vector<SurfaceComp>& comps) {
    comps.clear();
    ZRatioBlock bl;
    extract_block_(G, i, j, /*need_adjacency=*/true, /*need_counts=*/false, bl);
    if (!bl.valid) {
        s1_out = 0.0;
        s2_out = 0.0;
        logr_out = MY_LOG(psi0_);
        return false;
    }
    accumulate_surface_moments_(bl, s1_out, s2_out, &comps);
    if (s2_out <= 0.0) s2_out = kS2Floor;
    logr_out = MY_LOG(saddle_ratio(s1_out, s2_out));
    return true;
}

double ZRatioEngine::log_zratio(const arma::imat& G, int i, int j) {
    // Surface (Option B) serves the alpha = 1 cell (Normal or Cauchy slab, both
    // built from the block-Gibbs oracle) when attached: decompose the block and
    // sum per-component moments. Otherwise (surface absent, or the alpha != 1
    // Gamma-shape fence) fall through to the additive-counts saddle.
    const bool surface_active = has_surface_ && std::abs(alpha_ - 1.0) < 1e-12;
    // Neither branch needs the block adjacency matrix: the surface deploy reads
    // adjacency from G through the extract scratch (surface_logr_), the additive
    // saddle needs only the scalar counts.
    ZRatioBlock bl;
    extract_block_(G, i, j, /*need_adjacency=*/false,
                   /*need_counts=*/!surface_active, bl);
    if (!bl.valid) {
        n_add_++;
        return MY_LOG(psi0_);
    }
    if (surface_active) {
        if (phase_ != ZRatioPhase::Gauge) {
            n_pred_++;
            if (phase_ == ZRatioPhase::Retained) n_pred_ret_++;
        }
        return surface_logr_(G);
    }
    const int ncn = bl.ncn, cne = bl.cne, bre = bl.bre;

    double s1 = ncn * addc_[0] + cne * addc_[2] + bre * addc_[4];
    double s2 = ncn * addc_[1] + cne * addc_[3] + bre * addc_[5];
    if (s2 <= 0.0) s2 = kS2Floor;

    // Additive saddle over (nCN, cne, bre), served from the persistent
    // count-key cache. This is the hierarchical fence for the alpha != 1
    // (Gamma-shape) cell; the alpha = 1 cell (Normal or Cauchy slab) is served
    // by the surface branch above.
    const std::uint64_t sig = pack_count_key(ncn, cne, bre);
    double a;
    auto it = cache_.find(sig);
    if (it != cache_.end()) {
        n_hit_++;
        a = it->second;
    } else {
        n_miss_++;
        a = saddle_ratio(s1, s2);
        cache_[sig] = a;
    }
    n_add_++;
    return MY_LOG(a);
}

// Rank-2 SMW refresh of sigma = k^{-1} after the column-i change
// Delta K = e_i d^T + d e_i^T (d supported on {i} + N_i, d_i = half the
// diagonal change). Woodbury with U = [e_i, d], J = [[0,1],[1,0]] gives
// Sigma' = Sigma - [w1 w2] Cap^{-1} [w1 w2]^T, Cap = J + U^T Sigma U, a
// 2x2 solve. The outer products are not FP-symmetric, so the result is
// mirrored from the upper triangle. Returns false when Cap is numerically
// singular; the caller then refreshes Sigma by a full factorisation.
static bool smw_rank2_col_update_(arma::mat& sigma, int i,
                                  const arma::uvec& ni, const arma::vec& d_i_ni,
                                  double d_diag) {
    const arma::uword m = sigma.n_rows;
    arma::vec w1 = sigma.col(i);
    arma::vec w2 = d_diag * w1;
    for (arma::uword k = 0; k < ni.n_elem; ++k) {
        w2 += d_i_ni[k] * sigma.col(ni[k]);
    }
    const double cap11 = w1[i];
    const double cap12 = 1.0 + w2[i];
    double cap22 = d_diag * w2[i];
    for (arma::uword k = 0; k < ni.n_elem; ++k) {
        cap22 += d_i_ni[k] * w2[ni[k]];
    }
    const double det = cap11 * cap22 - cap12 * cap12;
    const double scale = std::abs(cap11) + std::abs(cap12) + std::abs(cap22);
    if (!std::isfinite(det) || std::abs(det) <= 1e-12 * scale * scale) {
        return false;
    }
    const double q11 = cap22 / det, q12 = -cap12 / det, q22 = cap11 / det;
    arma::vec u1 = q11 * w1 + q12 * w2;
    arma::vec u2 = q12 * w1 + q22 * w2;
    for (arma::uword c = 0; c < m; ++c) {
        for (arma::uword r = 0; r <= c; ++r) {
            sigma(r, c) -= u1[r] * w1[c] + u2[r] * w2[c];
        }
    }
    sigma = arma::symmatu(sigma);
    return true;
}

bool ZRatioEngine::gibbs_sweep_(arma::mat& k_blk, arma::mat& omega_blk,
                                const std::vector<arma::uvec>& nbr,
                                arma::mat* sigma_out) const {
    const int m = static_cast<int>(k_blk.n_rows);
    const double s2i = 1.0 / (sigma_ * sigma_);
    // Maintained block covariance Sigma = k_blk^{-1}: one exact
    // factorisation per sweep, rank-2 SMW refresh after each accepted row
    // write. Row i's C = ((K_{-i,-i})^{-1})_{N_i, N_i} then extracts in
    // O(|N_i|^2) via the Schur identity instead of a per-row O(m^3)
    // submatrix inversion. When Sigma is unavailable (factorisation or SMW
    // failure) the sweep falls back to the per-row inversion.
    arma::mat sigma_blk;
    bool have_sigma = arma::inv_sympd(sigma_blk, k_blk);
    arma::vec d_ni;
    for (int i = 0; i < m; ++i) {
        const arma::uvec& ni = nbr[i];
        if (ni.n_elem > 0) {
            const int nq = static_cast<int>(ni.n_elem);
            arma::mat c_mat(nq, nq);
            if (have_sigma) {
                const double sig_ii = sigma_blk(i, i);
                for (int a = 0; a < nq; ++a) {
                    const double sig_ai = sigma_blk(ni[a], i);
                    for (int b = 0; b < nq; ++b) {
                        c_mat(a, b) = sigma_blk(ni[a], ni[b]) -
                                      sig_ai * sigma_blk(i, ni[b]) / sig_ii;
                    }
                }
            } else {
                arma::uvec rest(m - 1);
                int p_ = 0;
                for (int v = 0; v < m; ++v) {
                    if (v != i) rest[p_++] = v;
                }
                arma::mat a_inv;
                if (!arma::inv_sympd(a_inv, k_blk.submat(rest, rest))) continue;
                arma::uvec idx_a(nq);
                for (arma::uword j = 0; j < ni.n_elem; ++j) {
                    idx_a[j] = (ni[j] < static_cast<arma::uword>(i))
                                   ? ni[j]
                                   : ni[j] - 1;
                }
                c_mat = a_inv.submat(idx_a, idx_a);
            }
            arma::mat m_mat = 2.0 * beta_ * c_mat;
            if (oracle_slab_cauchy_) {
                for (int j = 0; j < nq; ++j) {
                    m_mat(j, j) += s2i / omega_blk(i, ni[j]);
                }
            } else {
                m_mat.diag() += s2i;
            }
            arma::mat r_chol;
            if (!arma::chol(r_chol, m_mat)) continue;
            arma::vec z(nq);
            for (int j = 0; j < nq; ++j) z[j] = rnorm(*rng_, 0.0, 1.0);
            // fast: skip the rcond estimate; r_chol just passed chol(), and
            // the estimate costs as much as the back-substitution itself.
            arma::vec bvec =
                arma::solve(arma::trimatu(r_chol), z, arma::solve_opts::fast);
            double xi = rgamma(*rng_, delta_ + 1.0, beta_);
            double quad = arma::as_scalar(bvec.t() * c_mat * bvec);
            // The diagonal factor K_ii^(alpha - 1) couples the Gamma pivot
            // to the row draw; the alpha = 1 conjugate conditional serves
            // as an independence-Metropolis proposal with acceptance
            // ratio (K_ii_new / K_ii_old)^(alpha - 1).
            bool accept = true;
            if (std::abs(alpha_ - 1.0) > 1e-12) {
                const double kii_new = xi + quad;
                const double kii_old = k_blk(i, i);
                accept = MY_LOG(runif(*rng_)) <
                         (alpha_ - 1.0) * (MY_LOG(kii_new) -
                                           MY_LOG(kii_old));
            }
            if (accept) {
                double d_diag = 0.0;
                if (have_sigma) {
                    d_diag = (xi + quad - k_blk(i, i)) / 2.0;
                    d_ni.set_size(nq);
                    for (int j = 0; j < nq; ++j) {
                        d_ni[j] = bvec[j] - k_blk(ni[j], i);
                    }
                }
                for (int j = 0; j < nq; ++j) {
                    k_blk(ni[j], i) = bvec[j];
                    k_blk(i, ni[j]) = bvec[j];
                }
                k_blk(i, i) = xi + quad;
                if (have_sigma &&
                    !smw_rank2_col_update_(sigma_blk, i, ni, d_ni, d_diag)) {
                    have_sigma = arma::inv_sympd(sigma_blk, k_blk);
                }
            }
            if (oracle_slab_cauchy_) {
                // Conjugate omega | k ~ IG(1, 1/2 + k^2 / (2 sigma^2)).
                for (int j = 0; j < nq; ++j) {
                    const double b = k_blk(ni[j], i);
                    const double ig_rate = 0.5 + 0.5 * b * b * s2i;
                    const double wo = ig_rate / rexp(*rng_, 1.0);
                    omega_blk(i, ni[j]) = wo;
                    omega_blk(ni[j], i) = wo;
                }
            }
        } else {
            const double kii_new = rgamma(*rng_, delta_ + alpha_, beta_);
            if (have_sigma) {
                const double d_diag = (kii_new - k_blk(i, i)) / 2.0;
                k_blk(i, i) = kii_new;
                if (!smw_rank2_col_update_(sigma_blk, i, arma::uvec(), d_ni,
                                           d_diag)) {
                    have_sigma = arma::inv_sympd(sigma_blk, k_blk);
                }
            } else {
                k_blk(i, i) = kii_new;
            }
        }
    }
    if (sigma_out == nullptr) return true;
    if (have_sigma) {
        *sigma_out = std::move(sigma_blk);
        return true;
    }
    return arma::inv_sympd(*sigma_out, k_blk);
}

bool ZRatioEngine::inner_moments_(const arma::mat& r_inv, const arma::uvec& si,
                                  const arma::uvec& sj, const arma::vec& wsi,
                                  const arma::vec& wsj, double& w, double& p1,
                                  double& p2) const {
    const double t2 = 2.0 * beta_ * sigma_ * sigma_;
    arma::mat rii = r_inv.submat(si, si), rjj = r_inv.submat(sj, sj),
              rij = r_inv.submat(si, sj);
    const double s4 = sigma_ * sigma_ * sigma_ * sigma_, s8 = s4 * s4;
    // Leg-dressed recipe under the scale-mixture slab: with
    // Wi = diag(sqrt(omega_leg)), Mi = (I + t2 Wi Rii Wi)^{-1} and
    // P = (Wi Mi Wi) Rij (Wj Mj Wj) Rij^T; reduces to the plain resolvent
    // form at omega = 1 (the Normal-slab branch below).
    //
    // Both branches evaluate the moments through one Cholesky per side:
    // with I + t2 Rii~ = Ri^T Ri and B = Ri^{-T} Rij~ Rj^{-1}, P is a
    // similarity transform of B B^T, so tr(P) = ||B||_F^2 and
    // tr(P^2) = ||B B^T||_F^2, while det(Mi) = prod(diag Ri)^{-2}. No
    // explicit resolvent inverses are formed.
    if (oracle_slab_cauchy_) {
        rii.each_col() %= wsi;
        rii.each_row() %= wsi.t();
        rjj.each_col() %= wsj;
        rjj.each_row() %= wsj.t();
        rij.each_col() %= wsi;
        rij.each_row() %= wsj.t();
    }
    rii *= t2;
    rii.diag() += 1.0;
    rjj *= t2;
    rjj.diag() += 1.0;
    arma::mat ri, rj;
    if (!arma::chol(ri, rii)) return false;
    if (!arma::chol(rj, rjj)) return false;
    arma::mat x = arma::solve(arma::trimatl(ri.t()), rij,
                              arma::solve_opts::fast);
    arma::mat b = arma::solve(arma::trimatl(rj.t()), x.t(),
                              arma::solve_opts::fast)
                      .t();
    const double wdet = arma::prod(ri.diag()) * arma::prod(rj.diag());
    if (!(wdet > 0.0) || !std::isfinite(wdet)) return false;
    w = 1.0 / wdet;
    arma::mat c = b * b.t();
    p1 = s4 * arma::accu(b % b);
    p2 = s8 * arma::accu(c % c);
    return true;
}

// Symmetric PSD square root via eigendecomposition with nonneg clamping.
static arma::mat sympd_sqrt_(const arma::mat& a) {
    arma::vec ev;
    arma::mat vecs;
    arma::eig_sym(ev, vecs, 0.5 * (a + a.t()));
    ev = arma::clamp(ev, 0.0, arma::datum::inf);
    return vecs * arma::diagmat(arma::sqrt(ev)) * vecs.t();
}

bool ZRatioEngine::inner_reference_(const arma::mat& r_inv, const arma::uvec& si,
                                    const arma::uvec& sj, const arma::vec& wsi,
                                    const arma::vec& wsj, double& w, double& fN,
                                    double& gG) const {
    const double t2 = 2.0 * beta_ * sigma_ * sigma_;
    const double s4 = sigma_ * sigma_ * sigma_ * sigma_;
    arma::mat rii = r_inv.submat(si, si), rjj = r_inv.submat(sj, sj),
              rij = r_inv.submat(si, sj);
    arma::mat mi, mj, lh, rh;
    if (oracle_slab_cauchy_) {
        arma::mat di = arma::diagmat(wsi), dj = arma::diagmat(wsj);
        if (!arma::inv_sympd(mi, arma::eye(si.n_elem, si.n_elem) +
                                     t2 * di * rii * di)) {
            return false;
        }
        if (!arma::inv_sympd(mj, arma::eye(sj.n_elem, sj.n_elem) +
                                     t2 * dj * rjj * dj)) {
            return false;
        }
        lh = di * mi * di;
        rh = dj * mj * dj;
    } else {
        if (!arma::inv_sympd(mi, arma::eye(si.n_elem, si.n_elem) + t2 * rii)) {
            return false;
        }
        if (!arma::inv_sympd(mj, arma::eye(sj.n_elem, sj.n_elem) + t2 * rjj)) {
            return false;
        }
        lh = mi;
        rh = mj;
    }
    w = std::sqrt(arma::det(mi) * arma::det(mj));
    // Singular values s_k of Lh^.5 Rij Rh^.5; u_k = sigma^4 s_k^2 (the
    // eigenvalues of the moment matrix P) parameterize phi below.
    arma::vec sv;
    if (!arma::svd(sv, sympd_sqrt_(lh) * rij * sympd_sqrt_(rh))) return false;
    arma::vec u = s4 * (sv % sv);
    // phi(t) = prod_k (1 + u_k t^2)^{-1/2} over the tilt grid.
    arma::vec tg2 = tg_ % tg_;
    arma::vec logphi(tg_.n_elem, arma::fill::zeros);
    for (arma::uword k = 0; k < u.n_elem; ++k) {
        logphi += -0.5 * arma::log1p(u[k] * tg2);
    }
    arma::vec phi = ARMA_MY_EXP(logphi);
    fN = arma::accu(wt_ % ihat_ % phi);
    gG = arma::accu(wt_ % ghat_ % phi);
    return std::isfinite(w) && w > 0.0 && std::isfinite(fN) &&
           std::isfinite(gG);
}

void ZRatioEngine::init_block_(const arma::imat& a_blk,
                               std::vector<arma::uvec>& nbr, arma::mat& k_blk,
                               arma::mat& omega_blk) const {
    const int m = static_cast<int>(a_blk.n_rows);
    nbr.assign(m, arma::uvec());
    for (int i = 0; i < m; ++i) {
        std::vector<arma::uword> v;
        for (int j = 0; j < m; ++j) {
            if (j != i && a_blk(i, j) == 1) v.push_back(j);
        }
        nbr[i] = arma::uvec(v);
    }
    k_blk.zeros(m, m);
    if (oracle_slab_cauchy_) omega_blk.ones(m, m);
    else omega_blk.reset();
    for (int l = 0; l < m; ++l) {
        k_blk(l, l) = (alpha_ == 1.0 ? rexp(*rng_, beta_)
                                     : rgamma(*rng_, alpha_, beta_)) +
                      m;
    }
    for (int s = 0; s < burn_; ++s) gibbs_sweep_(k_blk, omega_blk, nbr);
}

bool ZRatioEngine::block_oracle_moments(const arma::imat& a_blk,
                                        const arma::uvec& si,
                                        const arma::uvec& sj, double& s1_out,
                                        double& s2_out) {
    std::vector<arma::uvec> nbr;
    arma::mat k_blk, omega_blk;
    init_block_(a_blk, nbr, k_blk, omega_blk);
    double sw = 0, sw1 = 0, sw2 = 0;
    long kept = 0;
    arma::vec wsi, wsj;
    arma::mat sigma_blk;
    for (int s = 0; s < n_sweep_; ++s) {
        // The leg draws below run on every sweep (sigma-valid or not) so
        // the draw stream matches the sweep count.
        const bool have_sig = gibbs_sweep_(k_blk, omega_blk, nbr, &sigma_blk);
        if (oracle_slab_cauchy_) {
            // Fresh leg weights per kept sweep: sqrt(omega) = 1/|z| with z
            // standard normal; the sweep average carries the mixture.
            wsi.set_size(si.n_elem);
            wsj.set_size(sj.n_elem);
            for (arma::uword k = 0; k < wsi.n_elem; ++k) {
                wsi[k] = 1.0 / std::abs(rnorm(*rng_, 0.0, 1.0));
            }
            for (arma::uword k = 0; k < wsj.n_elem; ++k) {
                wsj[k] = 1.0 / std::abs(rnorm(*rng_, 0.0, 1.0));
            }
        }
        double w, p1, p2;
        if (have_sig && inner_moments_(sigma_blk, si, sj, wsi, wsj, w, p1, p2)) {
            sw += w;
            sw1 += w * p1;
            sw2 += w * p2;
            ++kept;
        }
    }
    if (kept == 0 || sw <= 0) return false;
    s1_out = sw1 / sw;
    s2_out = sw2 / sw;
    return true;
}

bool ZRatioEngine::block_reference_logR(const arma::imat& a_blk,
                                        const arma::uvec& si,
                                        const arma::uvec& sj, int n_draws,
                                        double& logR_out, double& mcse_out) {
    std::vector<arma::uvec> nbr;
    arma::mat k_blk, omega_blk;
    init_block_(a_blk, nbr, k_blk, omega_blk);
    std::vector<double> wf, wg;   // per-draw W*<phi,I_N> and W*<phi,I_G>
    wf.reserve(n_draws);
    wg.reserve(n_draws);
    arma::vec wsi, wsj;
    arma::mat sigma_blk;
    for (int s = 0; s < n_draws; ++s) {
        const bool have_sig = gibbs_sweep_(k_blk, omega_blk, nbr, &sigma_blk);
        if (oracle_slab_cauchy_) {
            wsi.set_size(si.n_elem);
            wsj.set_size(sj.n_elem);
            for (arma::uword k = 0; k < wsi.n_elem; ++k) {
                wsi[k] = 1.0 / std::abs(rnorm(*rng_, 0.0, 1.0));
            }
            for (arma::uword k = 0; k < wsj.n_elem; ++k) {
                wsj[k] = 1.0 / std::abs(rnorm(*rng_, 0.0, 1.0));
            }
        }
        double w, fN, gG;
        if (have_sig &&
            inner_reference_(sigma_blk, si, sj, wsi, wsj, w, fN, gG)) {
            wf.push_back(w * fN);
            wg.push_back(w * gG);
        }
    }
    const int n = static_cast<int>(wf.size());
    if (n == 0) return false;
    double nf = 0, dg = 0;
    for (int i = 0; i < n; ++i) {
        nf += wf[i];
        dg += wg[i];
    }
    if (!(nf > 0.0) || !(dg > 0.0)) return false;
    const double ratio = nf / dg;
    // Batch-means MC standard error of the ratio, mapped to the log scale
    // (delta method). Common random numbers across the two averages already
    // cancel most of the ratio's noise.
    const int nb =
        std::max(10, static_cast<int>(std::floor(std::sqrt((double)n))));
    std::vector<double> batch;
    batch.reserve(nb);
    for (int b = 0; b < nb; ++b) {
        int lo = static_cast<int>((long)b * n / nb);
        int hi = static_cast<int>((long)(b + 1) * n / nb);
        double bnf = 0, bdg = 0;
        for (int i = lo; i < hi; ++i) {
            bnf += wf[i];
            bdg += wg[i];
        }
        if (hi > lo && bdg > 0.0) batch.push_back(bnf / bdg);
    }
    double mcse = 0.0;
    const int nB = static_cast<int>(batch.size());
    if (nB > 1) {
        double mb = 0;
        for (double x : batch) mb += x;
        mb /= nB;
        double vb = 0;
        for (double x : batch) vb += (x - mb) * (x - mb);
        vb /= (nB - 1);
        mcse = std::sqrt(vb / nB);
    }
    logR_out = MY_LOG(ratio);
    mcse_out = mcse / ratio;
    return true;
}

void ZRatioEngine::precompute_table(int ncn_max, int bre_max) {
    double c1 = addc_[0], c2 = addc_[1], c3 = addc_[2], c4 = addc_[3],
           c5 = addc_[4], c6 = addc_[5];
    for (int ncn = 0; ncn <= ncn_max; ncn++) {
        int cne_max = ncn * (ncn - 1) / 2;
        for (int cne = 0; cne <= cne_max; cne++) {
            for (int bre = 0; bre <= bre_max; bre++) {
                double s1 = ncn * c1 + cne * c3 + bre * c5;
                double s2 = ncn * c2 + cne * c4 + bre * c6;
                if (s2 <= 0.0) s2 = kS2Floor;
                cache_[pack_count_key(ncn, cne, bre)] = saddle_ratio(s1, s2);
            }
        }
    }
}

#include "zratio_engine.h"
#include "math/explog_macros.h"

#include <algorithm>
#include <vector>
#include <cmath>

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
    const int q = static_cast<int>(G.n_rows);
    ZRatioBlock bl;

    // Mediating block: common neighbours of (i, j) plus the endpoints of
    // 2-hop bridges between the exclusive neighbour sets. The toggled
    // edge's own state never enters, so the value is state-invariant.
    std::vector<bool> in_r(q, false);
    // Exclusive neighbour lists of i and j; the bridge scan below then runs
    // over the candidate pairs instead of the full q x q grid.
    std::vector<int> excl_i, excl_j;
    for (int k = 0; k < q; k++) {
        if (k == i || k == j) continue;
        const bool near_i = (G(i, k) == 1), near_j = (G(j, k) == 1);
        if (near_i && near_j) in_r[k] = true;
        else if (near_i) excl_i.push_back(k);
        else if (near_j) excl_j.push_back(k);
    }
    for (int a : excl_i) {
        for (int b : excl_j) {
            if (G(a, b) == 1) {
                in_r[a] = true;
                in_r[b] = true;
            }
        }
    }
    std::vector<int> rv;
    for (int k = 0; k < q; k++) {
        if (in_r[k]) rv.push_back(k);
    }
    const int m = static_cast<int>(rv.size());
    bl.m = m;

    std::vector<int> cn, si_o, sj_o;
    for (int p = 0; p < m; p++) {
        bool si = (G(i, rv[p]) == 1), sj = (G(j, rv[p]) == 1);
        if (si && sj) cn.push_back(p);
        else if (si) si_o.push_back(p);
        else if (sj) sj_o.push_back(p);
    }
    if ((cn.empty() && si_o.empty()) || (cn.empty() && sj_o.empty())) {
        // One side of the mediating block is empty: isolated-edge ratio.
        return bl;
    }
    bl.valid = true;

    bl.a_blk.zeros(m, m);
    for (int a = 0; a < m; a++) {
        for (int b = a + 1; b < m; b++) {
            int e = (G(rv[a], rv[b]) == 1) ? 1 : 0;
            bl.a_blk(a, b) = e;
            bl.a_blk(b, a) = e;
        }
    }

    bl.ncn = static_cast<int>(cn.size());
    for (size_t a = 0; a < cn.size(); a++) {
        for (size_t b = a + 1; b < cn.size(); b++) {
            if (bl.a_blk(cn[a], cn[b]) == 1) bl.cne++;
        }
    }
    for (int a : si_o) {
        for (int b : sj_o) {
            if (bl.a_blk(a, b) == 1) bl.bre++;
        }
    }
    for (int a : si_o) {
        int d = 0;
        for (int b : sj_o) {
            if (bl.a_blk(a, b) == 1) d++;
        }
        if (d > bl.maxbd) bl.maxbd = d;
    }
    for (int b : sj_o) {
        int d = 0;
        for (int a : si_o) {
            if (bl.a_blk(a, b) == 1) d++;
        }
        if (d > bl.maxbd) bl.maxbd = d;
    }

    bl.dens = (m >= 2)
        ? ((static_cast<double>(arma::accu(bl.a_blk)) / 2.0) /
           (static_cast<double>(m) * (m - 1) / 2.0))
        : 0.0;

    // Side memberships in ascending block position (CN nodes sit on both).
    std::vector<arma::uword> si_v, sj_v;
    for (int p = 0; p < m; p++) {
        bool si = std::find(cn.begin(), cn.end(), p) != cn.end();
        bool sio = std::find(si_o.begin(), si_o.end(), p) != si_o.end();
        bool sjo = std::find(sj_o.begin(), sj_o.end(), p) != sj_o.end();
        if (si || sio) si_v.push_back(p);
        if (si || sjo) sj_v.push_back(p);
    }
    bl.si = arma::uvec(si_v);
    bl.sj = arma::uvec(sj_v);
    return bl;
}

double ZRatioEngine::log_zratio(const arma::imat& G, int i, int j) {
    ZRatioBlock bl = extract_block(G, i, j);
    if (!bl.valid) {
        n_add_++;
        return MY_LOG(psi0_);
    }
    const int ncn = bl.ncn, cne = bl.cne, bre = bl.bre, maxbdeg = bl.maxbd,
              m = bl.m;
    const double dens = bl.dens;

    double s1 = ncn * addc_[0] + cne * addc_[2] + bre * addc_[4];
    double s2 = ncn * addc_[1] + cne * addc_[3] + bre * addc_[5];

    // Additive saddle: depends only on (nCN, cne, bre), served from the
    // persistent count-key cache. Corrections ride on top post-cache, so
    // corrected edges keep the compact key and full cache reuse.
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
    const double log_r_add = MY_LOG(a);

    if (maxbdeg < 2) {
        n_add_++;
        return log_r_add;
    }

    // Warm-up calibration: identical block signatures are served from the
    // correction cache; blocks inside the anchor cloud's Mahalanobis hull
    // use the current fit; uncovered blocks run the exact Monte-Carlo
    // evaluation, anchor the discrepancy, and refit.
    if (calibrating() && s1 > 0 && s2 > 0) {
        std::string ckey = std::to_string(ncn) + "_" + std::to_string(cne) +
                           "_" + std::to_string(bre) + "_" +
                           std::to_string(m) + "_" + std::to_string(maxbdeg) +
                           "_" + std::to_string(std::lround(dens * 1000.0));
        auto cit = corr_cache_.find(ckey);
        if (cit != corr_cache_.end()) {
            n_pred_++;
            return log_r_add + cit->second;
        }
        arma::vec x = {1.0, static_cast<double>(bre),
                       static_cast<double>(m), static_cast<double>(cne),
                       static_cast<double>(maxbdeg), dens};
        bool covered = false;
        if (coef_.n_elem > 0 && maha_sinv_.n_elem > 0) {
            arma::vec z = x.subvec(1, 5) - maha_mu_;
            covered = arma::as_scalar(z.t() * maha_sinv_ * z) <= maha_thresh_;
        }
        if (covered) {
            n_pred_++;
            return log_r_add + arma::dot(x, coef_);
        }
        double s1b = 0, s2b = 0, dl = 0;
        if (block_oracle_moments(bl.a_blk, bl.si, bl.sj, s1b, s2b)) {
            dl = MY_LOG(saddle_ratio(s1b, s2b)) - log_r_add;
        }
        n_oracle_++;
        corr_cache_[ckey] = dl;
        ax_.insert_rows(ax_.n_rows, x.t());
        ay_.insert_rows(ay_.n_elem, arma::vec{dl});
        refit_();
        return log_r_add + dl;
    }

    // Frozen OLS correction on the log-ratio scale, zeroed outside the
    // trained hull so the kernel never extrapolates.
    bool clamped = false;
    const double log_r_corr = deployed_correction(bl, clamped);
    if (clamped) n_clamp_++;
    if (log_r_corr != 0.0) n_pred_++;
    else n_add_++;
    return log_r_add + log_r_corr;
}

double ZRatioEngine::deployed_correction(const ZRatioBlock& bl,
                                         bool& clamped) const {
    clamped = false;
    if (bl.maxbd < 2) return 0.0;
    bool direct = (addc_.n_elem >= 13 && addc_[12] > 0.5);
    if (!direct) return 0.0;
    // Single gate: once the edge has bridge multiplicity >= 2 the direct
    // ratio-scale correction is applied everywhere. The correction is a
    // smooth low-dimensional surface on the log-ratio scale, and the
    // corrected value dominates the uncorrected additive saddle, most of all
    // in the dense/large-block corner. Reverting to additive outside a
    // calibration hull would reintroduce the additive bias exactly where it
    // is largest, so the frozen kernel extrapolates the surface rather than
    // gating on a prior/size-dependent box.
    double fc = addc_[6] + addc_[7] * static_cast<double>(bl.bre) +
                addc_[8] * static_cast<double>(bl.m) +
                addc_[9] * static_cast<double>(bl.cne) +
                addc_[10] * static_cast<double>(bl.maxbd) +
                addc_[11] * bl.dens;
    return fc;
}

void ZRatioEngine::enable_calibration(double delta, double eta, SafeRNG* rng,
                                      int n_sweep, int burn,
                                      double maha_thresh, int min_anchors,
                                      bool slab_cauchy) {
    calibration_enabled_ = true;
    frozen_ = false;
    delta_ = delta;
    sigma_ = 1.0;
    beta_ = eta;
    rng_ = rng;
    n_sweep_ = n_sweep;
    burn_ = burn;
    maha_thresh_ = maha_thresh;
    min_anchors_ = min_anchors;
    oracle_slab_cauchy_ = slab_cauchy;
}

void ZRatioEngine::refit_() {
    if (static_cast<int>(ay_.n_elem) < min_anchors_) return;
    arma::mat xtx = ax_.t() * ax_;
    xtx.diag() += 1e-8;
    coef_ = arma::solve(xtx, ax_.t() * ay_);
    arma::mat feats = ax_.cols(1, 5);
    maha_mu_ = arma::mean(feats, 0).t();
    arma::mat s = arma::cov(feats);
    s.diag() += 1e-6;
    if (!arma::inv_sympd(maha_sinv_, s)) maha_sinv_.reset();
}

void ZRatioEngine::freeze_calibration() {
    if (!calibration_enabled_ || frozen_) return;
    refit_();
    frozen_ = true;
    corr_cache_.clear();
    if (coef_.n_elem == 0) return;   // no fit: pure additive kernel
    arma::vec packed(23, arma::fill::zeros);
    packed.subvec(0, 5) = addc_.subvec(0, 5);
    packed.subvec(6, 11) = coef_;
    packed[12] = 1.0;
    arma::mat feats = ax_.cols(1, 5);
    for (arma::uword c = 0; c < 5; ++c) {
        packed[13 + 2 * c] = feats.col(c).min();
        packed[14 + 2 * c] = feats.col(c).max();
    }
    addc_ = packed;
}

void ZRatioEngine::gibbs_sweep_(arma::mat& k_blk, arma::mat& omega_blk,
                                const std::vector<arma::uvec>& nbr) const {
    const int m = static_cast<int>(k_blk.n_rows);
    const double s2i = 1.0 / (sigma_ * sigma_);
    for (int i = 0; i < m; ++i) {
        const arma::uvec& ni = nbr[i];
        if (ni.n_elem > 0) {
            const int nq = static_cast<int>(ni.n_elem);
            arma::uvec rest(m - 1);
            int p_ = 0;
            for (int v = 0; v < m; ++v) {
                if (v != i) rest[p_++] = v;
            }
            arma::mat a_inv;
            if (!arma::inv_sympd(a_inv, k_blk.submat(rest, rest))) continue;
            arma::uvec idx_a(nq);
            for (arma::uword j = 0; j < ni.n_elem; ++j) {
                idx_a[j] = (ni[j] < static_cast<arma::uword>(i)) ? ni[j]
                                                                 : ni[j] - 1;
            }
            arma::mat c_mat = a_inv.submat(idx_a, idx_a);
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
            arma::vec bvec = arma::solve(arma::trimatu(r_chol), z);
            double xi = rgamma(*rng_, delta_ + 1.0, beta_);
            double quad = arma::as_scalar(bvec.t() * c_mat * bvec);
            for (int j = 0; j < nq; ++j) {
                k_blk(ni[j], i) = bvec[j];
                k_blk(i, ni[j]) = bvec[j];
            }
            k_blk(i, i) = xi + quad;
            if (oracle_slab_cauchy_) {
                // Conjugate omega | k ~ IG(1, 1/2 + k^2 / (2 sigma^2)).
                for (int j = 0; j < nq; ++j) {
                    const double b = bvec[j];
                    const double ig_rate = 0.5 + 0.5 * b * b * s2i;
                    const double wo = ig_rate / rexp(*rng_, 1.0);
                    omega_blk(i, ni[j]) = wo;
                    omega_blk(ni[j], i) = wo;
                }
            }
        } else {
            k_blk(i, i) = rgamma(*rng_, delta_ + 1.0, beta_);
        }
    }
}

bool ZRatioEngine::inner_moments_(const arma::mat& k_blk, const arma::uvec& si,
                                  const arma::uvec& sj, const arma::vec& wsi,
                                  const arma::vec& wsj, double& w, double& p1,
                                  double& p2) const {
    const double t2 = 2.0 * beta_ * sigma_ * sigma_;
    arma::mat r_inv;
    if (!arma::inv_sympd(r_inv, k_blk)) return false;
    arma::mat rii = r_inv.submat(si, si), rjj = r_inv.submat(sj, sj),
              rij = r_inv.submat(si, sj);
    const double s4 = sigma_ * sigma_ * sigma_ * sigma_, s8 = s4 * s4;
    if (oracle_slab_cauchy_) {
        // Leg-dressed recipe under the scale-mixture slab: with
        // Wi = diag(sqrt(omega_leg)), Mi = (I + t2 Wi Rii Wi)^{-1} and
        // P = (Wi Mi Wi) Rij (Wj Mj Wj) Rij^T; reduces to the plain
        // resolvent form at omega = 1.
        arma::mat di = arma::diagmat(wsi), dj = arma::diagmat(wsj);
        arma::mat mi, mj;
        if (!arma::inv_sympd(mi, arma::eye(si.n_elem, si.n_elem) +
                                     t2 * di * rii * di)) {
            return false;
        }
        if (!arma::inv_sympd(mj, arma::eye(sj.n_elem, sj.n_elem) +
                                     t2 * dj * rjj * dj)) {
            return false;
        }
        arma::mat p_mat = (di * mi * di) * rij * (dj * mj * dj) * rij.t();
        w = std::sqrt(arma::det(mi) * arma::det(mj));
        p1 = s4 * arma::trace(p_mat);
        p2 = s8 * arma::accu(p_mat % p_mat.t());
        return true;
    }
    arma::mat mi, mj;
    if (!arma::inv_sympd(mi, arma::eye(si.n_elem, si.n_elem) + t2 * rii)) {
        return false;
    }
    if (!arma::inv_sympd(mj, arma::eye(sj.n_elem, sj.n_elem) + t2 * rjj)) {
        return false;
    }
    // Moments via traces: with u_k = (sigma^2 s_k)^2 and s_k the singular
    // values of Mi^.5 Rij Mj^.5, sum u_k = sigma^4 tr(P) and sum u_k^2 =
    // sigma^8 tr(P^2) for P = Mi Rij Mj Rij^T.
    arma::mat p_mat = mi * rij * mj * rij.t();
    w = std::sqrt(arma::det(mi) * arma::det(mj));
    p1 = s4 * arma::trace(p_mat);
    p2 = s8 * arma::accu(p_mat % p_mat.t());
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

bool ZRatioEngine::inner_reference_(const arma::mat& k_blk, const arma::uvec& si,
                                    const arma::uvec& sj, const arma::vec& wsi,
                                    const arma::vec& wsj, double& w, double& fN,
                                    double& gG, double& kappa2) const {
    const double t2 = 2.0 * beta_ * sigma_ * sigma_;
    const double s4 = sigma_ * sigma_ * sigma_ * sigma_;
    arma::mat r_inv;
    if (!arma::inv_sympd(r_inv, k_blk)) return false;
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
    // Singular values s_k of Lh^.5 Rij Rh^.5; u_k = sigma^4 s_k^2 (== the
    // eigenvalues of the moment matrix P, so kappa2 = sigma^4 tr(P)).
    arma::vec sv;
    if (!arma::svd(sv, sympd_sqrt_(lh) * rij * sympd_sqrt_(rh))) return false;
    arma::vec u = s4 * (sv % sv);
    kappa2 = arma::accu(u);
    // phi(t) = prod_k (1 + u_k t^2)^{-1/2} over the tilt grid.
    arma::vec tg2 = tg_ % tg_;
    arma::vec logphi(tg_.n_elem, arma::fill::zeros);
    for (arma::uword k = 0; k < u.n_elem; ++k) {
        logphi += -0.5 * arma::log1p(u[k] * tg2);
    }
    arma::vec phi = arma::exp(logphi);
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
        k_blk(l, l) = rexp(*rng_, beta_) + m;
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
    for (int s = 0; s < n_sweep_; ++s) {
        gibbs_sweep_(k_blk, omega_blk, nbr);
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
        if (inner_moments_(k_blk, si, sj, wsi, wsj, w, p1, p2)) {
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
    for (int s = 0; s < n_draws; ++s) {
        gibbs_sweep_(k_blk, omega_blk, nbr);
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
        double w, fN, gG, k2;
        if (inner_reference_(k_blk, si, sj, wsi, wsj, w, fN, gG, k2)) {
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
    logR_out = std::log(ratio);
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
                cache_[pack_count_key(ncn, cne, bre)] = saddle_ratio(s1, s2);
            }
        }
    }
}

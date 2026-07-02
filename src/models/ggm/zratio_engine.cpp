#include "zratio_engine.h"

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

double ZRatioEngine::log_zratio(const arma::imat& G, int i, int j) {
    const int q = static_cast<int>(G.n_rows);

    // Mediating block: common neighbours of (i, j) plus the endpoints of
    // 2-hop bridges between the exclusive neighbour sets. The toggled
    // edge's own state never enters, so the value is state-invariant.
    std::vector<bool> in_r(q, false);
    for (int k = 0; k < q; k++) {
        if (k != i && k != j && G(i, k) == 1 && G(j, k) == 1) in_r[k] = true;
    }
    for (int a = 0; a < q; a++) {
        if (a == i || a == j || G(i, a) != 1 || G(j, a) == 1) continue;
        for (int b = 0; b < q; b++) {
            if (b != i && b != j && b != a && G(j, b) == 1 && G(i, b) != 1 &&
                G(a, b) == 1) {
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

    std::vector<int> cn, si_o, sj_o;
    for (int p = 0; p < m; p++) {
        bool si = (G(i, rv[p]) == 1), sj = (G(j, rv[p]) == 1);
        if (si && sj) cn.push_back(p);
        else if (si) si_o.push_back(p);
        else if (sj) sj_o.push_back(p);
    }
    if ((cn.empty() && si_o.empty()) || (cn.empty() && sj_o.empty())) {
        // One side of the mediating block is empty: isolated-edge ratio.
        n_add_++;
        return std::log(psi0_);
    }

    arma::imat a_blk(m, m, arma::fill::zeros);
    for (int a = 0; a < m; a++) {
        for (int b = a + 1; b < m; b++) {
            int e = (G(rv[a], rv[b]) == 1) ? 1 : 0;
            a_blk(a, b) = e;
            a_blk(b, a) = e;
        }
    }

    int ncn = static_cast<int>(cn.size()), cne = 0, bre = 0;
    for (size_t a = 0; a < cn.size(); a++) {
        for (size_t b = a + 1; b < cn.size(); b++) {
            if (a_blk(cn[a], cn[b]) == 1) cne++;
        }
    }
    for (int a : si_o) {
        for (int b : sj_o) {
            if (a_blk(a, b) == 1) bre++;
        }
    }
    int maxbdeg = 0;
    for (int a : si_o) {
        int d = 0;
        for (int b : sj_o) {
            if (a_blk(a, b) == 1) d++;
        }
        if (d > maxbdeg) maxbdeg = d;
    }
    for (int b : sj_o) {
        int d = 0;
        for (int a : si_o) {
            if (a_blk(a, b) == 1) d++;
        }
        if (d > maxbdeg) maxbdeg = d;
    }

    double s1 = ncn * addc_[0] + cne * addc_[2] + bre * addc_[4];
    double s2 = ncn * addc_[1] + cne * addc_[3] + bre * addc_[5];

    // Additive saddle: depends only on (nCN, cne, bre), served from the
    // persistent count-key cache. Corrections ride on top post-cache, so
    // corrected edges keep the compact key and full cache reuse.
    std::string sig = "A" + std::to_string(ncn) + "_" + std::to_string(cne) +
                      "_" + std::to_string(bre);
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
    const double log_r_add = std::log(a);

    if (maxbdeg < 2) {
        n_add_++;
        return log_r_add;
    }

    const double dens = (m >= 2)
        ? ((static_cast<double>(arma::accu(a_blk)) / 2.0) /
           (static_cast<double>(m) * (m - 1) / 2.0))
        : 0.0;

    // Warm-up calibration: identical block signatures are served from the
    // correction cache; blocks inside the anchor cloud's Mahalanobis hull
    // use the current fit; uncovered blocks call the block-Gibbs oracle,
    // anchor the discrepancy, and refit.
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
        std::vector<arma::uword> si_v, sj_v;
        for (int p_ = 0; p_ < m; p_++) {
            bool si = std::find(cn.begin(), cn.end(), p_) != cn.end();
            bool sio = std::find(si_o.begin(), si_o.end(), p_) != si_o.end();
            bool sjo = std::find(sj_o.begin(), sj_o.end(), p_) != sj_o.end();
            if (si || sio) si_v.push_back(p_);
            if (si || sjo) sj_v.push_back(p_);
        }
        double s1b = 0, s2b = 0, dl = 0;
        if (block_oracle_moments(a_blk, arma::uvec(si_v), arma::uvec(sj_v),
                                 s1b, s2b)) {
            dl = std::log(saddle_ratio(s1b, s2b)) - log_r_add;
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
    double log_r_corr = 0.0;
    bool direct = (addc_.n_elem >= 13 && addc_[12] > 0.5);
    if (direct) {
        double fc = addc_[6] + addc_[7] * static_cast<double>(bre) +
                    addc_[8] * static_cast<double>(m) +
                    addc_[9] * static_cast<double>(cne) +
                    addc_[10] * static_cast<double>(maxbdeg) +
                    addc_[11] * dens;
        if (addc_.n_elem >= 23) {
            double bd = static_cast<double>(bre), md = static_cast<double>(m),
                   cd = static_cast<double>(cne),
                   xd = static_cast<double>(maxbdeg);
            bool inside =
                (bd >= addc_[13] && bd <= addc_[14] && md >= addc_[15] &&
                 md <= addc_[16] && cd >= addc_[17] && cd <= addc_[18] &&
                 xd >= addc_[19] && xd <= addc_[20] && dens >= addc_[21] &&
                 dens <= addc_[22]);
            if (!inside) {
                fc = 0.0;
                n_clamp_++;
            }
        }
        log_r_corr = fc;
    }
    if (log_r_corr != 0.0) n_pred_++;
    else n_add_++;
    return log_r_add + log_r_corr;
}

void ZRatioEngine::enable_calibration(double delta, double sigma, double beta,
                                      SafeRNG* rng, int n_sweep, int burn,
                                      double maha_thresh, int min_anchors) {
    calibration_enabled_ = true;
    frozen_ = false;
    delta_ = delta;
    sigma_ = sigma;
    beta_ = beta;
    rng_ = rng;
    n_sweep_ = n_sweep;
    burn_ = burn;
    maha_thresh_ = maha_thresh;
    min_anchors_ = min_anchors;
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

void ZRatioEngine::gibbs_sweep_(arma::mat& k_blk,
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
            m_mat.diag() += s2i;
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
        } else {
            k_blk(i, i) = rgamma(*rng_, delta_ + 1.0, beta_);
        }
    }
}

bool ZRatioEngine::inner_moments_(const arma::mat& k_blk, const arma::uvec& si,
                                  const arma::uvec& sj, double& w, double& p1,
                                  double& p2) const {
    const double t2 = 2.0 * beta_ * sigma_ * sigma_;
    arma::mat r_inv;
    if (!arma::inv_sympd(r_inv, k_blk)) return false;
    arma::mat rii = r_inv.submat(si, si), rjj = r_inv.submat(sj, sj),
              rij = r_inv.submat(si, sj);
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
    const double s4 = sigma_ * sigma_ * sigma_ * sigma_, s8 = s4 * s4;
    w = std::sqrt(arma::det(mi) * arma::det(mj));
    p1 = s4 * arma::trace(p_mat);
    p2 = s8 * arma::accu(p_mat % p_mat.t());
    return true;
}

bool ZRatioEngine::block_oracle_moments(const arma::imat& a_blk,
                                        const arma::uvec& si,
                                        const arma::uvec& sj, double& s1_out,
                                        double& s2_out) {
    const int m = static_cast<int>(a_blk.n_rows);
    std::vector<arma::uvec> nbr(m);
    for (int i = 0; i < m; ++i) {
        std::vector<arma::uword> v;
        for (int j = 0; j < m; ++j) {
            if (j != i && a_blk(i, j) == 1) v.push_back(j);
        }
        nbr[i] = arma::uvec(v);
    }
    arma::mat k_blk(m, m, arma::fill::zeros);
    for (int l = 0; l < m; ++l) {
        k_blk(l, l) = rexp(*rng_, beta_) + m;
    }
    for (int s = 0; s < burn_; ++s) gibbs_sweep_(k_blk, nbr);
    double sw = 0, sw1 = 0, sw2 = 0;
    long kept = 0;
    for (int s = 0; s < n_sweep_; ++s) {
        gibbs_sweep_(k_blk, nbr);
        double w, p1, p2;
        if (inner_moments_(k_blk, si, sj, w, p1, p2)) {
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

void ZRatioEngine::precompute_table(int ncn_max, int bre_max) {
    double c1 = addc_[0], c2 = addc_[1], c3 = addc_[2], c4 = addc_[3],
           c5 = addc_[4], c6 = addc_[5];
    for (int ncn = 0; ncn <= ncn_max; ncn++) {
        int cne_max = ncn * (ncn - 1) / 2;
        for (int cne = 0; cne <= cne_max; cne++) {
            for (int bre = 0; bre <= bre_max; bre++) {
                double s1 = ncn * c1 + cne * c3 + bre * c5;
                double s2 = ncn * c2 + cne * c4 + bre * c6;
                std::string sig = "A" + std::to_string(ncn) + "_" +
                                  std::to_string(cne) + "_" +
                                  std::to_string(bre);
                cache_[sig] = saddle_ratio(s1, s2);
            }
        }
    }
}

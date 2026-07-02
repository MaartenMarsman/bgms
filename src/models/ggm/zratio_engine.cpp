#include "zratio_engine.h"

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

    // OLS correction on the log-ratio scale, engaged for coupled-bridge
    // blocks (maxbdeg >= 2) inside the trained hull; zeroed outside it so
    // the frozen kernel falls back to the additive baseline.
    double log_r_corr = 0.0;
    bool direct = (addc_.n_elem >= 13 && addc_[12] > 0.5);
    if (direct && addc_.n_elem >= 12 && maxbdeg >= 2) {
        double dens = (m >= 2)
            ? ((static_cast<double>(arma::accu(a_blk)) / 2.0) /
               (static_cast<double>(m) * (m - 1) / 2.0))
            : 0.0;
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
        if (log_r_corr != 0.0) n_pred_++;
    }
    if (log_r_corr == 0.0) n_add_++;

    // The uncorrected saddle depends only on (nCN, cne, bre); the direct
    // correction rides on top post-cache, so corrected edges keep the
    // compact count key and full cache reuse.
    std::string sig = "A" + std::to_string(ncn) + "_" + std::to_string(cne) +
                      "_" + std::to_string(bre);
    auto it = cache_.find(sig);
    if (it != cache_.end()) {
        n_hit_++;
        return std::log(it->second) + log_r_corr;
    }
    n_miss_++;
    double a = saddle_ratio(s1, s2);
    cache_[sig] = a;
    return std::log(a) + log_r_corr;
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

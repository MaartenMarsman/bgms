#pragma once

#include <Rcpp.h>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>

// -----------------------------------------------------------------------------
// zratio_law.h
//
// DORMANT large-q insurance -- NOT wired into the default surface build. The
// deployed CN anchor source is the Monte-Carlo block oracle
// (zratio_block_oracle_moments): it is cheaper at the sizes bgms reaches
// (block-Gibbs anchor ~40-55 ms vs a law solve ~120-180 ms at q=50/size<=44)
// and ties it on gold accuracy, so the all-MC surface wins on both build cost
// and precision. This analytic law is kept only because its cost is ~flat in
// giant component size while the oracle grows ~cubically, so a crossover exists
// somewhere above the current reach; it is validated (test-zratio-law.R) and
// reachable via the test-only export zratio_law_moments should a large-q regime
// ever make it load-bearing. It is not used in the paper.
//
// The mu-first CPA analytic common-neighbour (CN) law: a deterministic,
// tableless C++ port of the companion's build-time CN anchor engine
// (mu-law-solve.cpp + eval_mu_law + re_param2). One self-consistent spectral
// solve per (size, density) cell yields the absolute per-component moments
// (S1, S2).
//
// Three pieces, all deterministic (no RNG), pure STL + R::dgamma:
//   * build_law_grids  -- the (eta, delta) grid setup (port of re_param2): the
//     coarse gamma quadrature and spectral grid for the solve, plus the fine
//     gamma grid and node law nu0 for the dressing.
//   * mu_law_solve     -- the complex-plane damped-Picard / theta-bisection /
//     psi-secant self-consistent solve (port of mu_law_solve_full), bounded by
//     a total inner-iteration budget; reports the residual psi_gap.
//   * eval_mu_law_dress-- the plain-Gamma dressing (port of eval_mu_law, tilt
//     branch): a 7-point stencil of the dressed Stieltjes profile at z = -t2
//     whose finite differences give S1 = s^4 n g'(-t2), S2 = s^8 n g'''(-t2)/6.
//
// Certification gate: dress only when abs(psi_gap) <= 1e-4; otherwise the cell
// is uncertified. The law covers the Normal slab (alpha = 1) CN family only.
// -----------------------------------------------------------------------------

namespace zratio_law {

// Trapezoid weights on a non-uniform grid: w[i] = (dx_left + dx_right) / 2,
// matching the companion's c(diff(x), 0) / 2 + c(0, diff(x)) / 2.
inline void trapezoid_weights(const std::vector<double>& x,
                              std::vector<double>& w) {
    const int n = static_cast<int>(x.size());
    w.assign(n, 0.0);
    for (int i = 0; i < n; ++i) {
        const double left = (i > 0) ? (x[i] - x[i - 1]) : 0.0;
        const double right = (i < n - 1) ? (x[i + 1] - x[i]) : 0.0;
        w[i] = 0.5 * (left + right);
    }
}

// Grids for one (eta, delta) cell. sigma = 1 in the standardized frame; the
// engine's eta is the companion's beta directly, and t2 = 2 beta sigma^2.
struct LawGrids {
    std::vector<double> gamc, gamw;  // coarse gamma quadrature (800) for solve
    std::vector<double> xs, dxw;     // spectral grid (370) for solve + dressing
    std::vector<double> gam, nu0;    // fine gamma grid (4000) + node law nu0
    double eps, t2, sigma2, delta, beta;
    int ng() const { return static_cast<int>(gamc.size()); }
    int nx() const { return static_cast<int>(xs.size()); }
};

// Port of re_param2: grids scale with the density via sc = max(1, 2 / beta);
// t2 = 2 beta sigma^2. Deterministic, built once per (eta, delta).
inline LawGrids build_law_grids(double eta, double delta) {
    LawGrids G;
    const double beta = eta, sigma = 1.0;
    G.beta = beta;
    G.sigma2 = sigma * sigma;
    G.delta = delta;
    G.t2 = 2.0 * beta * sigma * sigma;
    G.eps = 2e-4;
    const double sc = std::max(1.0, 2.0 / beta);

    // coarse gamma quadrature: (seq(sqrt(1e-6), sqrt(30 sc), 800))^2
    const int NG = 800;
    G.gamc.resize(NG);
    {
        const double lo = std::sqrt(1e-6), hi = std::sqrt(30.0 * sc);
        for (int j = 0; j < NG; ++j) {
            const double u = lo + (hi - lo) * j / (NG - 1);
            G.gamc[j] = u * u;
        }
    }
    trapezoid_weights(G.gamc, G.gamw);

    // spectral grid: c(seq(-1.5, -0.004, 70), (seq(1e-4, sqrt(60 sc), 300))^2)
    const int NXn = 70, NXp = 300;
    G.xs.resize(NXn + NXp);
    for (int i = 0; i < NXn; ++i) {
        G.xs[i] = -1.5 + (-0.004 - (-1.5)) * i / (NXn - 1);
    }
    {
        const double lo = 1e-4, hi = std::sqrt(60.0 * sc);
        for (int i = 0; i < NXp; ++i) {
            const double u = lo + (hi - lo) * i / (NXp - 1);
            G.xs[NXn + i] = u * u;
        }
    }
    trapezoid_weights(G.xs, G.dxw);

    // fine gamma grid + node law: nu0 prop exp(delta log gam - beta gam)
    const int NF = 4000;
    G.gam.resize(NF);
    G.nu0.resize(NF);
    {
        const double lo = 1e-6, hi = 30.0 * sc;
        double mx = -1e300;
        for (int j = 0; j < NF; ++j) {
            G.gam[j] = lo + (hi - lo) * j / (NF - 1);
            const double ld = delta * std::log(G.gam[j]) - beta * G.gam[j];
            G.nu0[j] = ld;
            if (ld > mx) mx = ld;
        }
        double s = 0.0;
        for (int j = 0; j < NF; ++j) {
            G.nu0[j] = std::exp(G.nu0[j] - mx);
            s += G.nu0[j];
        }
        for (int j = 0; j < NF; ++j) G.nu0[j] /= s;
    }
    return G;
}

// Converged density state of the solve (mu on the spectral grid + CPA scalars).
struct LawSol {
    std::vector<double> mu;
    double Zn, th, ga, cf, dg;
    double psi, psi_gap;
    int n_solve;
};

namespace detail {

// nuc: tilted gamma weights on the coarse grid for a given psi (log-sum-exp).
inline void make_nuc(const LawGrids& G, double psi, std::vector<double>& nuc) {
    double mx = -1e300;
    for (int j = 0; j < G.ng(); ++j) {
        const double ld = G.delta * std::log(G.gamc[j]) - G.beta * G.gamc[j]
                        + std::log(G.gamc[j])
                        - std::log(G.gamc[j] + G.t2 + psi);
        nuc[j] = ld;
        if (ld > mx) mx = ld;
    }
    double s = 0.0;
    for (int j = 0; j < G.ng(); ++j) {
        nuc[j] = std::exp(nuc[j] - mx) * G.gamw[j];
        s += nuc[j];
    }
    for (int j = 0; j < G.ng(); ++j) nuc[j] /= s;
}

struct Inner {
    std::vector<std::complex<double>> g;
    std::vector<double> mu;
    double Zn, th, ga, cf, dg;
    long bursts;
};

inline double gy_of(const LawGrids& G, const std::vector<double>& mu, double Zn,
                    double y) {
    double s = 0.0;
    for (int i = 0; i < G.nx(); ++i) s += G.dxw[i] * mu[i] / (G.xs[i] - y);
    return s / Zn;
}

// inner solve: damped Picard bursts + CPA scalar refresh (solve_mu_cpa logic).
inline void solve_inner(const LawGrids& G, const std::vector<double>& nuc,
                        double D, double n, double dmp, double tol, int maxit,
                        Inner& S) {
    const double rl = std::min(1.0, D / (n - 1.0));
    const int nx = G.nx(), ng = G.ng();
    if (static_cast<int>(S.g.size()) != nx) {
        S.g.assign(nx, std::complex<double>(0.0, 0.05));
    }
    std::vector<std::complex<double>> z(nx), kf(nx);
    for (int i = 0; i < nx; ++i) z[i] = std::complex<double>(G.xs[i], G.eps);
    double th = 1.0, cf = 1.0, ga = 0.0;
    for (int j = 0; j < ng; ++j) ga += nuc[j] / (G.gamc[j] + G.t2);
    S.mu.assign(nx, 0.0);
    double dg = 1e300, dmp_eff = dmp;
    long bursts = 0;
    int restarts = 0;
    const int nblk = (maxit + 9) / 10;
    for (int blk = 0; blk < nblk; ++blk) {
        for (int i = 0; i < nx; ++i) {
            kf[i] = D * cf * G.sigma2 * z[i] / (z[i] + th * G.t2);
        }
        for (int it = 0; it < 10; ++it) {  // one burst
            dg = 0.0;
            for (int i = 0; i < nx; ++i) {
                const std::complex<double> w = z[i] - kf[i] * (ga - S.g[i]);
                const double wr = w.real(), wi = w.imag();
                double s0 = 0.0, s1 = 0.0;
                for (int j = 0; j < ng; ++j) {
                    const double A = G.gamc[j] - wr;
                    const double inv = nuc[j] / (A * A + wi * wi);
                    s0 += inv;
                    s1 += A * inv;
                }
                const std::complex<double> gn(s1, wi * s0);
                const double d = std::abs(gn - S.g[i]);
                if (d > dg) dg = d;
                S.g[i] = (1.0 - dmp_eff) * S.g[i] + dmp_eff * gn;
            }
            ++bursts;
            if (dg < tol) break;
        }
        if (!std::isfinite(dg)) {                        // divergence guard
            if (++restarts > 3) { dg = 1e300; break; }
            S.g.assign(nx, std::complex<double>(0.0, 0.05));
            dmp_eff = std::max(0.1, 0.5 * dmp_eff);
            th = 1.0; cf = 1.0;
            ga = 0.0;
            for (int j = 0; j < ng; ++j) ga += nuc[j] / (G.gamc[j] + G.t2);
            dg = 1e300;
            continue;
        }
        double Zn = 0.0;
        for (int i = 0; i < nx; ++i) {
            S.mu[i] = S.g[i].imag() / M_PI;
            Zn += G.dxw[i] * S.mu[i];
        }
        S.Zn = std::max(Zn, 1e-12);
        double th_new = 1.0;
        if (rl < 1.0 - 1e-9) {
            auto cond = [&](double t) {
                const double Lb = G.t2 * gy_of(G, S.mu, S.Zn, -t * G.t2);
                return rl * (1.0 - t) / (1.0 + (1.0 - t) * Lb)
                     - (1.0 - rl) * t / (1.0 - t * Lb);
            };
            double lo = 1e-6, hi = 1.0 - 1e-9;
            double flo = cond(lo), fhi = cond(hi);
            (void) fhi;
            if (flo * fhi <= 0.0) {
                for (int b = 0; b < 80; ++b) {
                    const double mid = 0.5 * (lo + hi), fm = cond(mid);
                    if (flo * fm <= 0.0) { hi = mid; fhi = fm; }
                    else { lo = mid; flo = fm; }
                    if (hi - lo < 1e-13) break;
                }
                th_new = 0.5 * (lo + hi);
            } else {
                th_new = th;
            }
        }
        const double ga_new = gy_of(G, S.mu, S.Zn, -th_new * G.t2);
        const double dsc =
            std::max(std::abs(th_new - th), std::abs(ga_new - ga));
        th = th_new; ga = ga_new;
        const double Lb = G.t2 * ga;
        cf = 1.0 - Lb * (1.0 - th) / (1.0 + (1.0 - th) * Lb);
        if (dg < tol && dsc < 1e-9) break;
    }
    S.th = th; S.ga = ga; S.cf = cf; S.dg = dg; S.bursts = bursts;
}

inline double Q0_of_inner(const LawGrids& G, const Inner& S, double D) {
    const double a = S.th * G.t2;
    if (std::abs(a - G.t2) > 1e-6) {
        const double gy_t2 = gy_of(G, S.mu, S.Zn, -G.t2);
        return D * S.cf * G.sigma2 * (-G.t2) * (S.ga - gy_t2) / (a - G.t2);
    }
    double s = 0.0;
    for (int i = 0; i < G.nx(); ++i) {
        s += G.dxw[i] * S.mu[i] / std::pow(G.xs[i] + G.t2, 2.0);
    }
    return D * S.cf * G.sigma2 * G.t2 * s / S.Zn;
}

// one full psi pass: secant with staged tolerances (law_mu_fast logic), all
// inner solves drawing on a shared total-iteration budget.
inline void psi_pass(const LawGrids& G, double D, double n, double dmp,
                     int maxit_srch, double tol_psi, int max_eval, Inner& S,
                     std::vector<double>& nuc, double& best_p, double& best_f,
                     int& n_solve, long budget, long& used) {
    auto Fmap = [&](double p, double tol_in, int maxit_in) -> double {
        const long rem = budget - used;
        if (rem < 10L) return NAN;
        make_nuc(G, p, nuc);
        solve_inner(G, nuc, D, n, dmp, tol_in,
                    static_cast<int>(std::min(static_cast<long>(maxit_in), rem)),
                    S);
        used += S.bursts;
        ++n_solve;
        return Q0_of_inner(G, S, D);
    };
    double p0 = best_p;
    double q0 = Fmap(p0, 1e-7, maxit_srch);
    if (!std::isfinite(q0)) return;
    double f0 = q0 - p0;
    if (std::abs(f0) < best_f) { best_p = p0; best_f = std::abs(f0); }
    if (best_f < tol_psi) return;
    double p1 = 0.5 * (p0 + q0);
    for (int it = 0; it < max_eval; ++it) {
        if (used >= budget) return;
        const double tol_in =
            std::min(1e-7, std::max(1e-10, 1e-2 * std::abs(f0)));
        const double q1 = Fmap(p1, tol_in, maxit_srch);
        if (!std::isfinite(q1)) return;
        const double f1 = q1 - p1;
        if (std::abs(f1) < best_f) { best_p = p1; best_f = std::abs(f1); }
        if (std::abs(f1) < tol_psi) return;
        const double step = f1 * (p1 - p0) / (f1 - f0);
        double p2 = p1 - step;
        if (!std::isfinite(p2) || p2 < 0.0 ||
            std::abs(step) > 5.0 * std::abs(p1 - p0) + 1.0) {
            p2 = 0.5 * (p1 + q1);
        }
        p0 = p1; f0 = f1; p1 = p2;
    }
}

}  // namespace detail

// Full self-consistent solve (port of mu_law_solve_full): cold start, staged
// psi secant + hard-corner retry + final polish. Fills sol; sol.psi_gap is the
// convergence residual the caller gates on.
inline void mu_law_solve(const LawGrids& G, double D, double n, double tol_psi,
                         int max_eval, int budget_iters, LawSol& sol) {
    detail::Inner S;
    std::vector<double> nuc(G.ng());
    double best_p = 0.0, best_f = 1e300;
    int n_solve = 0;
    const long budget = std::max(static_cast<long>(budget_iters), 100L);
    long used = 0;

    detail::psi_pass(G, D, n, 0.5, 2000, tol_psi, max_eval, S, nuc, best_p,
                     best_f, n_solve, budget, used);
    if (best_f > 1e-6 && used < budget) {  // hard-corner retry
        detail::psi_pass(G, D, n, 0.25, 8000, tol_psi, 25, S, nuc, best_p,
                         best_f, n_solve, budget, used);
    }
    // final polish at the best psi so the dressing sees a converged mu.
    detail::make_nuc(G, best_p, nuc);
    detail::solve_inner(G, nuc, D, n, 0.5, 1e-10, best_f > 1e-4 ? 1500 : 6000, S);
    ++n_solve;
    sol.mu = S.mu;
    sol.Zn = S.Zn;
    sol.th = S.th;
    sol.ga = S.ga;
    sol.cf = S.cf;
    sol.dg = S.dg;
    sol.psi = best_p;
    sol.psi_gap = detail::Q0_of_inner(G, S, D) - best_p;
    sol.n_solve = n_solve;
}

// Plain-Gamma dressing (port of eval_mu_law, tilt = TRUE branch): a 7-point
// stencil of the dressed profile at z = -t2 whose finite differences give
//   S1 = s^4 n g'(-t2),  S2 = s^8 n g'''(-t2) / 6.
// Returns true and fills s1/s2 on success; false if the dressing produces a
// non-finite or non-positive moment (the caller then falls back to MC).
inline bool eval_mu_law_dress(const LawGrids& G, const LawSol& sol, double D,
                              double n, double& s1_out, double& s2_out) {
    const double h = 0.1;
    const double t2 = G.t2, sigma2 = G.sigma2;
    const double rl = std::min(1.0, D / (n - 1.0));
    const double a = sol.th * t2;
    const double cf = sol.cf;
    const int nx = G.nx(), nf = static_cast<int>(G.gam.size());
    const double Zn = sol.Zn;
    const std::vector<double>& mu = sol.mu;
    const std::vector<double>& XS = G.xs;
    const std::vector<double>& DXW = G.dxw;

    auto gy = [&](double y) {
        double s = 0.0;
        for (int i = 0; i < nx; ++i) s += DXW[i] * mu[i] / (XS[i] - y);
        return s / Zn;
    };
    auto Emu = [&](auto f) {
        double s = 0.0;
        for (int i = 0; i < nx; ++i) s += DXW[i] * mu[i] * f(XS[i]);
        return s / Zn;
    };

    const double Lb = t2 * sol.ga;
    const double gpa = Emu([&](double x) { return 1.0 / ((x + a) * (x + a)); });
    const double Qa = D * cf * sigma2 * a * gpa;
    auto Qof = [&](double y) {
        if (std::abs(y + a) < 1e-9) return Qa;
        return D * cf * sigma2 * y * (sol.ga - gy(y)) / (y + a);
    };

    double zs[7];
    for (int j = 0; j < 7; ++j) zs[j] = -t2 + (j - 3) * h;
    const double m1 = n - 1.0;
    const double tt2v = -t2 * (1.0 - sol.th) / (1.0 + (1.0 - sol.th) * Lb);
    const double t0v = sol.th * t2 / (1.0 - sol.th * Lb);
    const double Lam = rl * tt2v * tt2v + (1.0 - rl) * t0v * t0v;
    const double chi1 = Emu([&](double x) { return 1.0 / ((x + a) * (x + a)); });

    double Vj[7];
    for (int j = 0; j < 7; ++j) {
        const double zj = zs[j];
        auto A = [&](double x) { return 1.0 / x - 1.0 / (x - zj); };
        const double e_lad =
            Emu([&](double x) { return A(x) * x * x / ((x + a) * (x + a)); });
        const double lad = Lam * e_lad * e_lad / (1.0 - Lam * chi1);
        const double e2x2 = Emu([&](double x) { return A(x) * A(x) * x * x; });
        const double e2x3 =
            Emu([&](double x) { return A(x) * A(x) * x * x * x / (x + a); });
        const double e2x4 = Emu([&](double x) {
            return A(x) * A(x) * x * x * x * x / ((x + a) * (x + a));
        });
        Vj[j] = (sigma2 * sigma2 / (t2 * t2)) * m1 *
                (2.0 * (e2x2 - 2.0 * e2x3 + e2x4) + 3.0 * lad);
        if (!(Vj[j] > 0.0) || !std::isfinite(Vj[j])) return false;
    }

    const double Q0v = Qof(-t2);
    double gd[7];
    const int NQ = 400;
    std::vector<double> qs(NQ), wq(NQ);
    for (int j = 0; j < 7; ++j) {
        const double zj = zs[j];
        const double Qj = Qof(zj);
        const double qhi = Qj + 14.0 * std::sqrt(Vj[j]);
        const double shape = Qj * Qj / Vj[j], scale = Vj[j] / Qj;
        double wsum = 0.0;
        for (int k = 0; k < NQ; ++k) {
            qs[k] = 1e-8 + (qhi - 1e-8) * k / (NQ - 1);
            wq[k] = R::dgamma(qs[k], shape, scale, 0);
            wsum += wq[k];
        }
        for (int k = 0; k < NQ; ++k) wq[k] /= wsum;

        double num = 0.0, nb = 0.0;
        for (int k = 0; k < NQ; ++k) {
            const double q0 = qs[k] * Q0v / Qj;
            const double denom_shift = qs[k] - zj;  // gam + (qs - zj) > 0
            double sn = 0.0, sb = 0.0;
            for (int gi = 0; gi < nf; ++gi) {
                const double gam = G.gam[gi];
                const double tl = gam / (gam + t2 + q0);
                const double nug = G.nu0[gi];
                sb += nug * tl;
                sn += nug * tl / (gam + denom_shift);
            }
            num += wq[k] * sn;
            nb += wq[k] * sb;
        }
        gd[j] = num / nb;
        if (!std::isfinite(gd[j])) return false;
    }

    const double d1 =
        (-gd[0] + 9.0 * gd[1] - 45.0 * gd[2] + 45.0 * gd[4] - 9.0 * gd[5] +
         gd[6]) / (60.0 * h);
    const double d3 =
        (gd[0] - 8.0 * gd[1] + 13.0 * gd[2] - 13.0 * gd[4] + 8.0 * gd[5] -
         gd[6]) / (8.0 * h * h * h);
    s1_out = sigma2 * sigma2 * n * d1;
    s2_out = sigma2 * sigma2 * sigma2 * sigma2 * n * d3 / 6.0;
    return std::isfinite(s1_out) && std::isfinite(s2_out) &&
           s1_out > 0.0 && s2_out > 0.0;
}

// One-shot certified CN moment: build grids, solve, gate on psi_gap, dress.
// Returns true with (s1, s2) on a certified, positive solve; false otherwise
// (uncertified psi, or a non-finite/non-positive dressing). psi_gap is always
// written for diagnostics.
inline bool law_moments(double eta, double delta, double D, double n,
                        int budget_iters, double& s1_out, double& s2_out,
                        double& psi_gap_out) {
    const LawGrids G = build_law_grids(eta, delta);
    LawSol sol;
    mu_law_solve(G, D, n, 1e-7, 15, budget_iters, sol);
    psi_gap_out = sol.psi_gap;
    if (!std::isfinite(sol.psi_gap) || std::abs(sol.psi_gap) > 1e-4) {
        return false;
    }
    return eval_mu_law_dress(G, sol, D, n, s1_out, s2_out);
}

}  // namespace zratio_law

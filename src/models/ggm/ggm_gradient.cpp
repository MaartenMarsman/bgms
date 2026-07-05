#include "models/ggm/ggm_gradient.h"

#include <cmath>
#include <limits>
#include <stdexcept>

// =====================================================================
// rebuild
// =====================================================================

void GGMGradientEngine::rebuild(
    const GraphConstraintStructure& structure,
    size_t n,
    const arma::mat& suf_stat,
    const BaseParameterPrior& interaction_prior,
    const BaseParameterPrior& diagonal_prior,
    double determinant_tilt)
{
    structure_ = &structure;
    n_ = n;
    p_ = structure.p;
    suf_stat_ = &suf_stat;
    interaction_prior_ = &interaction_prior;
    diagonal_prior_ = &diagonal_prior;
    delta_ = determinant_tilt;
}

void GGMGradientEngine::rebuild(const GraphConstraintStructure& structure) {
    structure_ = &structure;
    p_ = structure.p;
    n_ = 0;
    suf_stat_ = nullptr;
    interaction_prior_ = nullptr;
    diagonal_prior_ = nullptr;
    delta_ = 0.0;
}

// =====================================================================
// build_Aq
// =====================================================================
// Build constraint matrix A_q for column q.
//
// For each excluded index i (rows where edge (i,q) is absent),
// row r of A_q is: A_q[r, 0:i] = Phi[0:i, i].
// This encodes the zero-constraint on K[i,q] = sum_l Phi[l,i]*Phi[l,q].

void GGMGradientEngine::build_Aq(
    const arma::mat& Phi,
    const ColumnConstraints& col,
    size_t q,
    arma::mat& Aq)
{
    Aq.zeros(col.m_q, q);
    for (size_t r = 0; r < col.m_q; ++r) {
        size_t i = col.excluded_indices[r];
        for (size_t l = 0; l <= i; ++l) {
            Aq(r, l) = Phi(l, i);
        }
    }
}

// =====================================================================
// givens_qr
// =====================================================================
// Compute null-space basis N_q and QR R-diagonal from A_q.
//
// Uses QR of A_q^T (which is q x m_q). The last (q - rank) columns
// of Q form the null-space basis. The diagonal of R gives the
// Givens QR decomposition of n x m matrix M (n >= m).
//
// Computes M = Q R via bottom-to-top Givens rotations. Stores the
// rotation sequence for reverse-mode differentiation in the backward
// pass. The convention is:
//   G = [[c, s], [-s, c]]  applied to rows (r1, r2) of W from the left.
//   Q accumulates G^T on columns.
// R diagonal entries are always positive by construction.

void GGMGradientEngine::givens_qr(
    const arma::mat& M,
    arma::mat& Q,
    arma::mat& R,
    arma::vec& R_diag,
    std::vector<GivensRotation>& rots)
{
    size_t n = M.n_rows;
    size_t m = M.n_cols;
    R = M;            // working copy, will become R
    Q.eye(n, n);      // accumulate Q
    rots.clear();

    for (size_t j = 0; j < m; ++j) {
        if (j + 1 < n) {
            for (size_t i = n - 1; i > j; --i) {
                size_t r1 = i - 1;
                size_t r2 = i;
                double a = R(r1, j);
                double b = R(r2, j);
                double r_val = std::sqrt(a * a + b * b);

                if (r_val < 1e-300) continue;  // skip identity rotation

                double c = a / r_val;
                double s = b / r_val;

                rots.push_back({c, s, r1, r2, j});

                // G * R: rotate rows (r1, r2) of R
                for (size_t k = j; k < m; ++k) {
                    double w1 = R(r1, k), w2 = R(r2, k);
                    R(r1, k) =  c * w1 + s * w2;
                    R(r2, k) = -s * w1 + c * w2;
                }
                // Q * G^T: rotate columns (r1, r2) of Q
                for (size_t l = 0; l < n; ++l) {
                    double q1 = Q(l, r1), q2 = Q(l, r2);
                    Q(l, r1) =  c * q1 + s * q2;
                    Q(l, r2) = -s * q1 + c * q2;
                }
            }
        }
    }

    // Extract positive R diagonal
    size_t rank = std::min(n, m);
    R_diag.set_size(rank);
    for (size_t j = 0; j < rank; ++j) {
        R_diag(j) = std::abs(R(j, j));
    }
}

// =====================================================================
// forward_map
// =====================================================================

ForwardMapResult GGMGradientEngine::forward_map(const arma::vec& theta) const {
    ForwardMapResult result;
    result.Phi.zeros(p_, p_);
    result.psi.set_size(p_);
    result.Nq.resize(p_);
    result.R_diag.resize(p_);
    result.givens_rotations.resize(p_);
    result.Q_full.resize(p_);
    result.R_full.resize(p_);

    arma::mat Aq_buf;  // reusable buffer for A_q

    for (size_t q = 0; q < p_; ++q) {
        const auto& col = structure_->columns[q];
        size_t offset = structure_->theta_offsets[q];

        if (q == 0) {
            // Column 0: just the diagonal
            double psi_q = theta(offset);
            result.psi(q) = psi_q;
            result.Phi(0, 0) = std::exp(psi_q);
            continue;
        }

        // Extract f_q from theta
        size_t d_q = col.d_q;
        arma::vec f_q;
        if (d_q > 0) {
            f_q = theta.subvec(offset, offset + d_q - 1);
        }

        // psi_q is after f_q
        double psi_q = theta(offset + d_q);
        result.psi(q) = psi_q;
        result.Phi(q, q) = std::exp(psi_q);

        // Build constraint matrix A_q
        size_t m_q = col.m_q;

        if (m_q == 0 && d_q == q) {
            // No constraints: x_q = f_q directly. N_q is the identity and is
            // never materialized; the backward pass copies x_bar straight
            // into f_bar for m_q == 0 columns.
            result.Nq[q].reset();
            result.R_diag[q].reset();
            result.givens_rotations[q].clear();
            for (size_t k = 0; k < d_q; ++k) {
                result.Phi(k, q) = f_q(k);
            }
        } else if (d_q == 0) {
            // Fully constrained: x_q = 0
            // Givens QR for R_diag (Jacobian) and stored rotations
            if (m_q > 0) {
                build_Aq(result.Phi, col, q, Aq_buf);
                givens_qr(Aq_buf.t(),
                          result.Q_full[q], result.R_full[q],
                          result.R_diag[q], result.givens_rotations[q]);
            }
            result.Nq[q].reset();
            // x_q stays zero (already zeroed)
        } else {
            // General case: build A_q, Givens QR, null space
            build_Aq(result.Phi, col, q, Aq_buf);
            givens_qr(Aq_buf.t(),
                      result.Q_full[q], result.R_full[q],
                      result.R_diag[q], result.givens_rotations[q]);

            // N_q = last d_q columns of Q
            result.Nq[q] = result.Q_full[q].cols(m_q, q - 1);

            // x_q = N_q f_q
            arma::vec x_q = result.Nq[q] * f_q;
            for (size_t k = 0; k < q; ++k) {
                result.Phi(k, q) = x_q(k);
            }
        }
    }

    // Jacobian: log|det J| = p*log(2) + 2*sum(psi) + sum_{i<p}(p-i)*psi_i
    //                       - sum_q sum_j log|R_{q,jj}|
    double ldj = static_cast<double>(p_) * std::log(2.0);
    for (size_t q = 0; q < p_; ++q) {
        ldj += 2.0 * result.psi(q);
    }
    for (size_t i = 0; i + 1 < p_; ++i) {
        ldj += static_cast<double>(p_ - 1 - i) * result.psi(i);
    }
    for (size_t q = 1; q < p_; ++q) {
        const auto& rd = result.R_diag[q];
        for (size_t j = 0; j < rd.n_elem; ++j) {
            ldj -= std::log(rd(j));
        }
    }
    result.log_det_jacobian = ldj;

    // K = Phi^T Phi
    result.K = result.Phi.t() * result.Phi;

    return result;
}

// =====================================================================
// logp_and_gradient
// =====================================================================
// Null-space gradient using the Phi-space adjoint for the initial
// Phi_bar computation, eliminating the K^{-1} triangular solve.
//
// The backward pass has two phases:
//   1. Phi_bar from data + priors via direct Phi-space formulation:
//        Phi_bar = -(Phi*S + 2*Phi) + Cauchy adjoint.
//      This replaces the K^{-1} = Phi^{-1} Phi^{-T} solve and the
//      K_bar_sym matrix multiply from the original approach.
//   2. Reverse-Givens cross-column adjoint (unchanged): differentiates
//      through the Givens QR to capture both the Jacobian gradient
//      and the null-space basis rotation.
//
// Since log|K| = 2*sum(psi) depends only on the diagonal of Phi,
// the off-diagonal entries of Phi_bar.col(q) are identical to those
// from the K^{-1} approach. The log-det contribution to psi_bar is
// added explicitly (+n) instead of flowing through K^{-1}.

std::pair<double, arma::vec> GGMGradientEngine::logp_and_gradient(
    const arma::vec& theta) const
{
    // --- Forward pass: theta -> Phi, K via null-space constraints ---
    ForwardMapResult fm = forward_map(theta);
    const arma::mat& Phi = fm.Phi;
    const arma::mat& K = fm.K;
    const arma::mat& S = *suf_stat_;
    double n = static_cast<double>(n_);


    // Check for degenerate Phi (extreme theta pushed by leapfrog).
    double min_diag = Phi.diag().min();
    if (!std::isfinite(min_diag) || min_diag < 1e-15) {
        return {-std::numeric_limits<double>::infinity(),
                arma::vec(theta.n_elem, arma::fill::zeros)};
    }

    // P = Phi * S — reused for value and gradient.
    // Phi is upper triangular; trimatu dispatches to BLAS dtrmm,
    // halving the FLOP count vs dense gemm.
    arma::mat P = arma::trimatu(Phi) * S;

    // --- Log-posterior value ---
    double log_det_K = 2.0 * arma::accu(fm.psi);
    double log_lik = (n / 2.0) * log_det_K - 0.5 * arma::accu(Phi % P);

    // Interaction prior on included off-diagonal K entries.
    // Applied on the K_yy = -K/2 scale to match the mixed-MRF convention,
    // so that `interaction_prior` refers to the same quantity across all
    // bgms model families (pure GGM and mixed MRF continuous blocks).
    double log_slab = 0.0;
    for (size_t q = 1; q < p_; ++q) {
        for (size_t i : structure_->columns[q].included_indices) {
            log_slab += interaction_prior_->logp(-0.5 * K(i, q));
        }
    }

    // Prior on the partial-association diagonal: -K_yy_{ii} = K_{ii}/2.
    // Evaluated at 0.5 * K_{ii} so `diagonal_prior_` refers to the same
    // quantity (the partial association) as `interaction_prior_`.
    double log_diag_prior = 0.0;
    for (size_t i = 0; i < p_; ++i) {
        log_diag_prior += diagonal_prior_->logp(0.5 * K(i, i));
    }

    // Determinant tilt: adds delta_ * log|K| = 2*delta_ * sum(psi) to the
    // log-prior. Pushes the chain away from the PD-cone boundary.
    double log_tilt = delta_ * log_det_K;

    double lp = log_lik + log_slab + log_diag_prior + fm.log_det_jacobian + log_tilt;

    if (!std::isfinite(lp)) {
        return {lp, arma::vec(theta.n_elem, arma::fill::zeros)};
    }

    // --- Backward pass ---

    // Phase 1: Phi_bar from data + priors (direct Phi-space, no K^{-1})
    // Data:  d/dPhi [-0.5 tr(Phi^T Phi S)] = -Phi S = -P
    // Diagonal prior: d/dPhi [log p(K_ii/2)]
    //   = (1/2) * grad(K_ii/2) * dK_ii/dPhi
    //   = (1/2) * grad(K_ii/2) * 2 Phi(:,i) = grad(K_ii/2) * Phi(:,i).
    arma::mat Phi_bar = -P;
    for (size_t i = 0; i < p_; ++i) {
        double dg = diagonal_prior_->grad(0.5 * K(i, i));
        Phi_bar.col(i).head(i + 1) += dg * Phi.col(i).head(i + 1);
    }

    // Interaction prior adjoint on included edges.
    // Prior is on K_yy_{ij} = -0.5 * K_{ij}, so the chain rule gives
    //   d/dK_{ij} log p(-0.5 K_{ij}) = -0.5 * prior'(-0.5 K_{ij})
    for (size_t q = 1; q < p_; ++q) {
        for (size_t i : structure_->columns[q].included_indices) {
            double kyy_ij = -0.5 * K(i, q);
            double d = -0.5 * interaction_prior_->grad(kyy_ij);
            Phi_bar.col(q).head(i + 1) += d * Phi.col(i).head(i + 1);
            Phi_bar.col(i).head(i + 1) += d * Phi.col(q).head(i + 1);
        }
    }

    // Phase 2: reverse-Givens extraction of the theta gradient from Phi_bar.
    // The log-det data term contributes +n to every psi (log|K| = 2*sum(psi));
    // the determinant tilt contributes +2*delta.
    arma::vec gradient(theta.n_elem, arma::fill::zeros);
    theta_gradient_from_phi_bar(theta, 0, fm, Phi_bar, n + 2.0 * delta_,
                                gradient);

    return {lp, gradient};
}

// =====================================================================
// theta_gradient_from_phi_bar
// =====================================================================
// Reverse-mode extraction of the (f_q, psi_q) gradient from a caller-
// seeded Phi-space adjoint. Processes columns right-to-left, extracting
// the theta gradient and accumulating the cross-column adjoint into
// Phi_bar (which is consumed as workspace).

void GGMGradientEngine::theta_gradient_from_phi_bar(
    const arma::vec& theta,
    size_t theta_offset,
    const ForwardMapResult& fm,
    arma::mat& Phi_bar,
    double psi_extra,
    arma::vec& gradient) const
{
    const arma::mat& Phi = fm.Phi;

    for (size_t q = p_; q-- > 0; ) {
        const auto& col = structure_->columns[q];
        size_t offset = theta_offset + structure_->theta_offsets[q];
        size_t d_q = col.d_q;

        // --- psi_q gradient ---
        // Phi_bar(q,q) * Phi(q,q) = chain rule through exp(psi_q)
        // +2 from Jacobian: d/dpsi [2*psi] = 2
        // +(p-1-q) from Jacobian: d/dpsi [(p-1-q)*psi] for q < p-1
        // +psi_extra from the caller (log-det data term, determinant tilt)
        double psi_bar = Phi_bar(q, q) * Phi(q, q);
        psi_bar += psi_extra + 2.0;
        if (q + 1 < p_) {
            psi_bar += static_cast<double>(p_ - 1 - q);
        }
        gradient(offset + d_q) = psi_bar;

        if (q == 0) continue;

        // --- f_q gradient via N_q ---
        size_t m_q = col.m_q;

        if (d_q > 0 && m_q == 0) {
            // Unconstrained column: N_q = I, so f_bar = x_bar directly.
            for (size_t k = 0; k < d_q; ++k) {
                gradient(offset + k) = Phi_bar(k, q);
            }
            continue;
        }

        arma::vec x_bar = Phi_bar.col(q).head(q);

        if (d_q > 0) {
            const arma::mat& Nq = fm.Nq[q];
            arma::vec f_bar = Nq.t() * x_bar;
            for (size_t k = 0; k < d_q; ++k) {
                gradient(offset + k) = f_bar(k);
            }
        }

        // --- Cross-column adjoint (reverse-Givens) ---
        if (m_q == 0) continue;

        const auto& rotations = fm.givens_rotations[q];
        size_t n_rot = rotations.size();

        // Initialize W_bar (R_bar) and Q_bar with seed adjoints.
        arma::mat W_bar(q, m_q, arma::fill::zeros);
        arma::mat Q_bar(q, q, arma::fill::zeros);

        size_t rank = std::min(m_q, q);
        const arma::mat& R = fm.R_full[q];
        for (size_t j = 0; j < rank; ++j) {
            W_bar(j, j) = -1.0 / R(j, j);
        }

        if (d_q > 0) {
            arma::vec f_q = theta.subvec(offset, offset + d_q - 1);
            for (size_t k = 0; k < d_q; ++k) {
                for (size_t l = 0; l < q; ++l) {
                    Q_bar(l, m_q + k) = x_bar(l) * f_q(k);
                }
            }
        }

        arma::mat Q_work = fm.Q_full[q];
        arma::mat W_work = fm.R_full[q];

        for (size_t k = n_rot; k-- > 0; ) {
            const auto& rot = rotations[k];
            double cc = rot.c, ss = rot.s;
            size_t r1 = rot.r1, r2 = rot.r2;

            for (size_t cj = 0; cj < m_q; ++cj) {
                double w1 = W_work(r1, cj), w2 = W_work(r2, cj);
                W_work(r1, cj) =  cc * w1 - ss * w2;
                W_work(r2, cj) =  ss * w1 + cc * w2;
            }
            for (size_t ri = 0; ri < q; ++ri) {
                double q1 = Q_work(ri, r1), q2 = Q_work(ri, r2);
                Q_work(ri, r1) =  cc * q1 - ss * q2;
                Q_work(ri, r2) =  ss * q1 + cc * q2;
            }

            double c_bar = 0.0, s_bar = 0.0;
            for (size_t cj = 0; cj < m_q; ++cj) {
                c_bar += W_bar(r1, cj) * W_work(r1, cj) +
                         W_bar(r2, cj) * W_work(r2, cj);
                s_bar += W_bar(r1, cj) * W_work(r2, cj) -
                         W_bar(r2, cj) * W_work(r1, cj);
            }
            for (size_t ri = 0; ri < q; ++ri) {
                c_bar += Q_bar(ri, r1) * Q_work(ri, r1) +
                         Q_bar(ri, r2) * Q_work(ri, r2);
                s_bar += Q_bar(ri, r1) * Q_work(ri, r2) -
                         Q_bar(ri, r2) * Q_work(ri, r1);
            }

            for (size_t cj = 0; cj < m_q; ++cj) {
                double wb1 = W_bar(r1, cj), wb2 = W_bar(r2, cj);
                W_bar(r1, cj) = cc * wb1 - ss * wb2;
                W_bar(r2, cj) = ss * wb1 + cc * wb2;
            }

            for (size_t ri = 0; ri < q; ++ri) {
                double qb1 = Q_bar(ri, r1), qb2 = Q_bar(ri, r2);
                Q_bar(ri, r1) = cc * qb1 - ss * qb2;
                Q_bar(ri, r2) = ss * qb1 + cc * qb2;
            }

            size_t j_col = rot.col;
            double a = W_work(r1, j_col);
            double b = W_work(r2, j_col);
            double r_val = std::sqrt(a * a + b * b);
            if (r_val > 1e-300) {
                double r3 = r_val * r_val * r_val;
                W_bar(r1, j_col) += (c_bar * b * b - s_bar * a * b) / r3;
                W_bar(r2, j_col) += (-c_bar * a * b + s_bar * a * a) / r3;
            }
        }

        for (size_t r = 0; r < m_q; ++r) {
            size_t i = col.excluded_indices[r];
            for (size_t l = 0; l <= i; ++l) {
                Phi_bar(l, i) += W_bar(l, r);
            }
        }
    }
}


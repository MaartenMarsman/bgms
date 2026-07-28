#include "models/ggm/ggm_model.h"
#include "rng/rng_utils.h"
#include "math/explog_macros.h"
#include "math/cholupdate.h"
#include "mcmc/execution/chain_result.h"
#include "mcmc/execution/step_result.h"
#include "mcmc/execution/warmup_schedule.h"

void GGMModel::collect_chain_diagnostics(ChainResult& chain_result) const {
    if (!zratio_engine_) return;
    const ZRatioEngine& engine = *zratio_engine_;
    chain_result.has_zratio_diagnostics = true;
    chain_result.zratio_addc = engine.addc();
    chain_result.zratio_counters = {
        static_cast<double>(engine.n_hit()),
        static_cast<double>(engine.n_miss()),
        static_cast<double>(engine.n_pred()),
        static_cast<double>(engine.n_add()),
        static_cast<double>(engine.cache_size()),
        static_cast<double>(engine.n_extrap()),
        static_cast<double>(engine.max_extrap_size())
    };
    if (zratio_gauge_.n_sweeps > 0) {
        chain_result.zratio_gauge_ran = true;
        chain_result.zratio_gauge_D = zratio_gauge_.D();
        chain_result.zratio_gauge_noise_floor = zratio_gauge_.noise_floor();
        chain_result.zratio_gauge_se_mean = zratio_gauge_.se_mean();
        chain_result.zratio_gauge_se_sd = zratio_gauge_.se_sd();
        chain_result.zratio_gauge_se_mcse = zratio_gauge_.se_mcse();
        chain_result.zratio_gauge_n_ent = zratio_gauge_.n_ent;
        chain_result.zratio_gauge_n_ref = zratio_gauge_.n_ref;
        chain_result.zratio_gauge_n_capped = zratio_gauge_.n_capped;
        chain_result.zratio_gauge_pair_i =
            arma::conv_to<arma::ivec>::from(zratio_gauge_.rec_i);
        chain_result.zratio_gauge_pair_j =
            arma::conv_to<arma::ivec>::from(zratio_gauge_.rec_j);
        chain_result.zratio_gauge_pair_se = arma::vec(zratio_gauge_.rec_se);
        chain_result.zratio_gauge_pair_mcse =
            arma::vec(zratio_gauge_.rec_mcse);
    }
}

// =====================================================================
// NUTS gradient support
// =====================================================================

void GGMModel::ensure_constraint_structure() {
    if (!constraint_dirty_) return;
    constraint_structure_.build(edge_indicators_);
    gradient_engine_.rebuild(constraint_structure_, n_, suf_stat_, *interaction_prior_, *diagonal_prior_, determinant_tilt_);
    constraint_dirty_ = false;
    theta_valid_ = false;
}

void GGMModel::recompute_theta() const {
    if (theta_valid_) return;

    // Build constraint structure (const-safe: structure is already built
    // by ensure_constraint_structure before any gradient call)
    const auto& cs = constraint_structure_;
    theta_.set_size(cs.active_dim);

    arma::mat Aq_buf;

    for (size_t q = 0; q < p_; ++q) {
        const auto& col = cs.columns[q];
        size_t offset = cs.theta_offsets[q];

        // psi_q = log(phi_qq)
        double psi_q = MY_LOG(cholesky_of_precision_(q, q));
        theta_(offset + col.d_q) = psi_q;

        if (q == 0 || col.d_q == 0) continue;

        // Build A_q, compute null-space basis N_q via Givens QR
        arma::mat Q_tmp, R_tmp;
        arma::vec R_diag;
        std::vector<GivensRotation> rots_tmp;
        GGMGradientEngine::build_Aq(cholesky_of_precision_, col, q, Aq_buf);
        GGMGradientEngine::givens_qr(Aq_buf.t(), Q_tmp, R_tmp, R_diag, rots_tmp);
        arma::mat Nq = Q_tmp.cols(col.m_q, q - 1);

        // f_q = N_q^T x_q
        arma::vec x_q = cholesky_of_precision_.col(q).head(q);
        arma::vec f_q = Nq.t() * x_q;

        for (size_t k = 0; k < col.d_q; ++k) {
            theta_(offset + k) = f_q(k);
        }
    }

    theta_valid_ = true;
}

size_t GGMModel::parameter_dimension() const {
    // Lazy: if constraint structure hasn't been built, use full dimension
    if (constraint_dirty_) {
        return p_ + p_ * (p_ - 1) / 2;
    }
    return constraint_structure_.active_dim;
}

size_t GGMModel::full_parameter_dimension() const {
    return p_ + p_ * (p_ - 1) / 2;
}

arma::vec GGMModel::get_vectorized_parameters() const {
    // Ensure the constraint structure is built so we can compute theta
    if (constraint_dirty_) {
        // const_cast is safe: ensure_constraint_structure only modifies
        // the constraint cache, not the model state
        const_cast<GGMModel*>(this)->ensure_constraint_structure();
    }
    recompute_theta();
    return theta_;
}

arma::vec GGMModel::get_full_vectorized_parameters() const {
    if (constraint_dirty_) {
        const_cast<GGMModel*>(this)->ensure_constraint_structure();
    }
    recompute_theta();

    const auto& cs = constraint_structure_;
    arma::vec full(cs.full_dim, arma::fill::zeros);

    // Scatter each active theta entry into its slot in the full (zero-padded)
    // vector; excluded off-diagonal slots stay zero.
    cs.for_each_active_full_pair([&](size_t active_idx, size_t full_idx) {
        full(full_idx) = theta_(active_idx);
    });

    return full;
}

void GGMModel::set_vectorized_parameters(const arma::vec& parameters) {
    ensure_constraint_structure();

    // Run forward map: theta -> Phi -> K
    const ForwardMapResult& fm = gradient_engine_.forward_map(parameters);

    // Update internal state
    precision_matrix_ = fm.K;
    cholesky_of_precision_ = fm.Phi;
    bool ok = arma::solve(inv_cholesky_of_precision_, arma::trimatu(cholesky_of_precision_),
                          arma::eye(p_, p_), arma::solve_opts::fast);
    if (!ok) {
        refresh_cholesky();
    } else {
        covariance_matrix_ = inv_cholesky_of_precision_ * inv_cholesky_of_precision_.t();
        log_det_precision_ = cholesky_helpers::get_log_det(cholesky_of_precision_);
    }

    // Cache theta
    theta_ = parameters;
    theta_valid_ = true;
}

std::pair<double, arma::vec> GGMModel::logp_and_gradient(
    const arma::vec& parameters)
{
    ensure_constraint_structure();
    return gradient_engine_.logp_and_gradient(parameters);
}


arma::vec GGMModel::get_active_inv_mass() const {
    if (constraint_dirty_) {
        const_cast<GGMModel*>(this)->ensure_constraint_structure();
    }

    const auto& cs = constraint_structure_;

    if (inv_mass_.n_elem == 0) {
        return arma::ones<arma::vec>(cs.active_dim);
    }

    // The theta-space (unconstrained) NUTS path is the only caller of this
    // function. get_full_vectorized_parameters() scatters each active theta
    // entry into a definite slot in the full-dim vector (included off-diag
    // slot for f_q[k]; full_psi_offset for psi_q; excluded slots stay 0).
    // Welford in NUTSAdaptationController therefore tracks the variance of
    // the theta entry at that slot directly, so the active inverse mass is
    // simply the included-slot subset — no N_q rotation needed.
    if (inv_mass_.n_elem == cs.full_dim) {
        // Gather the full-dim inverse mass back into the active layout — the
        // exact inverse of the get_full_vectorized_parameters scatter.
        arma::vec active(cs.active_dim);
        cs.for_each_active_full_pair([&](size_t active_idx, size_t full_idx) {
            active(active_idx) = inv_mass_(full_idx);
        });
        return active;
    }

    // Fallback: return inv_mass_ as-is (dimensions should match active_dim)
    return inv_mass_;
}


void GGMModel::get_constants(size_t i, size_t j) {
    // GGM stores K directly, so the precision entries are precision_matrix_.
    constants_ = cholesky_helpers::precision_proposal_constants(
        log_det_precision_, covariance_matrix_, i, j,
        precision_matrix_(i, j), precision_matrix_(j, j));
}

double GGMModel::constrained_diagonal(const double x) const {
    if (x == 0) {
        return constants_[5];
    } else {
        return constants_[4] + std::pow((x - constants_[2]) / constants_[3], 2);
    }
}

double GGMModel::log_density_impl(const arma::mat& omega, const arma::mat& phi) const {

    double logdet_omega = cholesky_helpers::get_log_det(phi);
    double trace_prod = arma::accu(omega % suf_stat_);

    double log_likelihood = n_ * (p_ * MY_LOG(2 * arma::datum::pi) / 2 + logdet_omega / 2) - trace_prod / 2;

    return log_likelihood;
}

double GGMModel::log_det_ratio_edge(size_t i, size_t j) const {
    // Rank-2 matrix-determinant lemma: log|K_prop| - log|K_curr| where K_prop
    // differs from K_curr at entries (i,j), (j,i), and (j,j). GGM stores K
    // directly, so the precision scalars are precision_matrix_ entries.
    return cholesky_helpers::log_det_ratio_edge_kernel(
        covariance_matrix_, i, j,
        precision_matrix_(i, j), precision_proposal_(i, j),
        precision_matrix_(j, j), precision_proposal_(j, j));
}

double GGMModel::log_det_ratio_diag(size_t j) const {
    return cholesky_helpers::log_det_ratio_diag_kernel(
        covariance_matrix_, j,
        precision_matrix_(j, j), precision_proposal_(j, j));
}

double GGMModel::log_density_impl_edge(size_t i, size_t j) const {
    // Log-likelihood ratio (not the full log-likelihood).
    double Ui2 = precision_matrix_(i, j) - precision_proposal_(i, j);
    double Uj2 = (precision_matrix_(j, j) - precision_proposal_(j, j)) / 2;
    double logdet = log_det_ratio_edge(i, j);
    double trace_prod = -2 * (suf_stat_(j, j) * Uj2 + suf_stat_(i, j) * Ui2);
    return (n_ * logdet - trace_prod) / 2;
}

double GGMModel::log_density_impl_diag(size_t j) const {
    double Uj2 = (precision_matrix_(j, j) - precision_proposal_(j, j)) / 2;
    double logdet = log_det_ratio_diag(j);
    double trace_prod = -2 * suf_stat_(j, j) * Uj2;
    return (n_ * logdet - trace_prod) / 2;
}

bool GGMModel::proposal_is_positive_definite_() const {
    arma::mat R_chk;
    return arma::chol(R_chk, precision_proposal_);
}

double GGMModel::ggm_edge_move(size_t i, size_t j) {
    get_constants(i, j);
    double Phi_q1q  = constants_[0];
    (void)constants_[1]; // Phi_q1q1 computed in get_constants but unused here

    size_t e = j * (j + 1) / 2 + i; // parameter index in vectorized form (column-major upper triangle)
    double proposal_sd = proposal_sds_(e);

    double phi_prop       = rnorm(rng_, Phi_q1q, proposal_sd);
    double omega_prop_q1q = constants_[2] + constants_[3] * phi_prop;
    double omega_prop_qq  = constrained_diagonal(omega_prop_q1q);

    // Only the (i,j), (j,i), (j,j) entries of the proposal are read on the
    // data path; the full-matrix copy is needed only for the prior-only PD
    // check below.
    if (n_ == 0) precision_proposal_ = precision_matrix_;
    precision_proposal_(i, j) = omega_prop_q1q;
    precision_proposal_(j, i) = omega_prop_q1q;
    precision_proposal_(j, j) = omega_prop_qq;

    // Prior-only chains have no likelihood anchor vetoing non-PD proposals;
    // reject them explicitly (see proposal_is_positive_definite_).
    if (n_ == 0 && !proposal_is_positive_definite_()) {
        return -arma::datum::inf;
    }

    double ln_alpha = log_density_impl_edge(i, j);

    // Determinant-tilt prior: |K|^delta contributes
    //   delta * (log|K_prop| - log|K_curr|)
    // to the MH ratio. log_det_ratio_edge uses the rank-2 matrix-determinant
    // lemma in O(p) via the cached covariance, so this is essentially free.
    if (determinant_tilt_ != 0.0) {
        ln_alpha += determinant_tilt_ * log_det_ratio_edge(i, j);
    }

    // Interaction prior on K_yy_{ij} = -0.5 * Omega_{ij}
    ln_alpha += interaction_prior_->logp(-0.5 * precision_proposal_(i, j));
    ln_alpha -= interaction_prior_->logp(-0.5 * precision_matrix_(i, j));

    // Gamma(shape, rate) prior on changed diagonal K_jj. The Roverato move
    // slaves K_jj = c_3 + phi_{q-1,q}^2, so K_jj moves with the off-diagonal
    // and its prior must be re-evaluated.
    ln_alpha += diagonal_prior_->logp(0.5 * precision_proposal_(j, j));
    ln_alpha -= diagonal_prior_->logp(0.5 * precision_matrix_(j, j));

    if (MY_LOG(runif(rng_)) < ln_alpha) {
        double omega_ij_old = precision_matrix_(i, j);
        double omega_jj_old = precision_matrix_(j, j);

        precision_matrix_(i, j) = omega_prop_q1q;
        precision_matrix_(j, i) = omega_prop_q1q;
        precision_matrix_(j, j) = omega_prop_qq;

        cholesky_update_after_edge(omega_ij_old, omega_jj_old, i, j);
    }

    return ln_alpha;
}

double GGMModel::update_edge_parameter(size_t i, size_t j) {
    if (edge_indicators_(i, j) == 0) {
        return 0.0; // Edge is not included; skip update (AR irrelevant, masked out)
    }
    return std::min(1.0, MY_EXP(ggm_edge_move(i, j)));
}

void GGMModel::cholesky_update_after_edge(double omega_ij_old, double omega_jj_old, size_t i, size_t j)
{

    v2_[0] = omega_ij_old - precision_proposal_(i, j);
    v2_[1] = (omega_jj_old - precision_proposal_(j, j)) / 2;

    vf1_[i] = v1_[0];
    vf1_[j] = v1_[1];
    vf2_[i] = v2_[0];
    vf2_[j] = v2_[1];

    // log|K| advances by the rank-2 matrix-determinant-lemma ratio -- the
    // same quantity the MH ratio uses -- against the pre-update covariance.
    // chol(K) itself is not maintained per accept; every sweep ends with a
    // refresh_cholesky() that recomputes the factor and an exact log-det.
    log_det_precision_ += cholesky_helpers::log_det_ratio_edge_kernel(
        covariance_matrix_, i, j,
        omega_ij_old, precision_proposal_(i, j),
        omega_jj_old, precision_proposal_(j, j));

    const arma::uvec support =
        {static_cast<arma::uword>(i), static_cast<arma::uword>(j)};
    apply_rank2_chol_smw_update_(support);

    // reset for next iteration
    vf1_[i] = 0.0;
    vf1_[j] = 0.0;
    vf2_[i] = 0.0;
    vf2_[j] = 0.0;

}

void GGMModel::apply_rank2_chol_smw_update_(const arma::uvec& support)
{
    // K_new = K_old + vf1 vf2^T + vf2 vf1^T = K_old + u1 u1^T - u2 u2^T,
    // where u1 = (vf1 + vf2) / sqrt(2), u2 = (vf1 - vf2) / sqrt(2).
    //
    // `support` lists the nonzero indices of vf1/vf2 (hence of u1/u2): {i,j}
    // for the edge accept, {i} + N_i for a row-Gibbs row. Sigma * u then
    // touches only those columns -- O(p |support|) instead of the dense
    // O(p^2) gemv.
    //
    // chol(K) is not maintained here: nothing reads the factor between
    // accepts within a sweep (the moves read Sigma and the incrementally
    // maintained log-det), and every sweep ends with refresh_cholesky().
    u1_ = (vf1_ + vf2_) / sqrt(2);
    u2_ = (vf1_ - vf2_) / sqrt(2);

    // Sherman-Morrison-Woodbury rank-2 update of covariance_matrix_ = inv(K),
    // O(p^2) total:
    //   K_new = K_old + M D M^T, M = [u1, u2], D = diag(+1, -1)
    //   inv(K_new) = inv(K_old) - A C^{-1} A^T,
    //     A = inv(K_old) M = [a1, a2],
    //     C = D^{-1} + M^T inv(K_old) M = diag(+1, -1) + symmetric 2x2.
    // Capacitance singularity (|det C| ~ 0) falls back to refresh_cholesky().
    // inv_cholesky_of_precision_ is not maintained here: only
    // refresh_cholesky() and set_vectorized_parameters() write it, and
    // nothing reads it between accepts.
    // a1 = Sigma u1, a2 = Sigma u2 (both dense length p -- the outer-product
    // Sigma updates below stay O(p^2)). u1/u2 are zero outside `support`, so
    // when the support is small we gather just those columns once and matvec
    // against them (O(p |support|)); when it is near-full the gather copy
    // costs more than the dense gemv, so fall back. Gather wins while
    // ~3 |support| < 2p.
    arma::vec a1, a2;
    if (3 * support.n_elem < 2 * p_) {
        const arma::mat Scols = covariance_matrix_.cols(support);
        a1 = Scols * u1_.elem(support);
        a2 = Scols * u2_.elem(support);
    } else {
        a1 = covariance_matrix_ * u1_;
        a2 = covariance_matrix_ * u2_;
    }
    // u1/u2 vanish off `support`, so the capacitance dots restrict to it.
    const arma::vec u1s = u1_.elem(support);
    const arma::vec u2s = u2_.elem(support);
    double c11 =  1.0 + arma::dot(u1s, a1.elem(support));
    double c12 =        arma::dot(u1s, a2.elem(support));
    double c22 = -1.0 + arma::dot(u2s, a2.elem(support));
    double det = c11 * c22 - c12 * c12;
    if (!std::isfinite(det) || std::abs(det) < 1e-14) {
        refresh_cholesky();
    } else {
        const double inv_c00 =  c22 / det;
        const double inv_c11 =  c11 / det;
        const double inv_c01 = -c12 / det;
        // Delta Sigma = -inv_c00 a1 a1^T - inv_c11 a2 a2^T
        //               - inv_c01 (a1 a2^T + a2 a1^T)
        // regrouped as two rank-1 outer products with bundled weights:
        //   Delta Sigma = -(a1 b1^T + a2 b2^T),
        //     b1 = inv_c00 a1 + inv_c01 a2,
        //     b2 = inv_c01 a1 + inv_c11 a2.
        const arma::vec b1 = inv_c00 * a1 + inv_c01 * a2;
        const arma::vec b2 = inv_c01 * a1 + inv_c11 * a2;
        covariance_matrix_ -= a1 * b1.t();
        covariance_matrix_ -= a2 * b2.t();
        // The outer products are symmetric in exact arithmetic but not in
        // floating point (b1(j) rounds once, so a1(i) b1(j) != a1(j) b1(i));
        // downstream chol() calls on Sigma-derived submatrices require exact
        // symmetry. Mirror the upper triangle in place -- symmatu on a
        // self-assignment would materialise a p x p temporary per accept.
        for (arma::uword c = 1; c < p_; ++c) {
            for (arma::uword r = 0; r < c; ++r) {
                covariance_matrix_(c, r) = covariance_matrix_(r, c);
            }
        }
        // Validate the touched rows of Sigma K = I. When K passes near a
        // singular state (a tiny row-Gibbs xi draw at n = 0, delta = 0 makes
        // this legitimate, not exceptional), Sigma legitimately blows up to
        // ~1/xi; the SMW update that moves K away from that state then
        // subtracts two huge outer products and cancellation destroys Sigma
        // in absolute terms. Everything downstream (proposal constants,
        // the row-Gibbs Schur matrix) reads Sigma and would write a
        // non-positive-definite K from the garbage. The probe costs
        // O(p |support|) -- the same order as the matvec above -- and
        // repairs the cache from K the moment accuracy is lost.
        if (!sigma_rows_consistent_(support)) {
            refresh_cholesky();
        }
    }
}

bool GGMModel::row_block_gibbs_eligible() {
    // Scope: a Normal or Cauchy slab on the off-diagonals and any Gamma(., .)
    // on K_ii/2. The determinant tilt is a shift in the xi Gamma shape; a
    // Gamma shape != 1 is an independent-MH correction; a Cauchy slab is the
    // scale-mixture-of-normals with a per-edge weight omega. None of these
    // gate eligibility.
    const bool normal = dynamic_cast<const NormalPrior*>(interaction_prior_.get()) != nullptr;
    const bool cauchy = dynamic_cast<const CauchyPrior*>(interaction_prior_.get()) != nullptr;
    if (!normal && !cauchy) return false;
    if (dynamic_cast<const GammaScalePrior*>(diagonal_prior_.get()) == nullptr) return false;
    return true;
}

double GGMModel::slab_scale_() const {
    if (const auto* n = dynamic_cast<const NormalPrior*>(interaction_prior_.get()))
        return n->scale();
    if (const auto* c = dynamic_cast<const CauchyPrior*>(interaction_prior_.get()))
        return c->scale();
    return 1.0;  // unreachable once row_block_gibbs_eligible() holds
}

bool GGMModel::slab_is_cauchy_() const {
    return dynamic_cast<const CauchyPrior*>(interaction_prior_.get()) != nullptr;
}

void GGMModel::refresh_cauchy_omega_() {
    if (!slab_is_cauchy_()) return;
    const double sigma = slab_scale_();
    const double two_sig2 = 2.0 * sigma * sigma;
    for (size_t i = 0; i + 1 < p_; ++i) {
        for (size_t j = i + 1; j < p_; ++j) {
            double w;
            if (edge_indicators_(i, j) == 1) {
                // Active: omega | K ~ InvGamma(1, 1/2 + kyy^2/(2 sigma^2)),
                // drawn as 1 / Gamma(1, rate) since InvGamma(a, b) = 1/Gamma(a, b).
                const double kyy = -0.5 * precision_matrix_(i, j);
                w = 1.0 / rgamma(rng_, 1.0, 0.5 + kyy * kyy / two_sig2);
            } else {
                // Inactive (K_ij = 0): omega ~ prior InvGamma(1/2, 1/2). Kept
                // fresh so an edge-selection add can condition on a valid weight.
                w = 1.0 / rgamma(rng_, 0.5, 0.5);
            }
            omega_(i, j) = w;
            omega_(j, i) = w;
        }
    }
}

void GGMModel::update_row_block_gibbs(size_t i) {
    // Gaussian-Gamma draw of (beta = K_{N_i, i}, kii = K_{i,i}) given
    // A = K_{-i, -i}, S, and the slab x Gamma prior. The determinant tilt
    // enters as a shift in the xi shape; a Gamma shape alpha != 1 is corrected
    // by an independent-MH accept below; a Cauchy slab enters through the
    // per-edge weights omega_ (fixed at 1 for a Normal slab).
    //
    // Slab sigma (std on K_yy_ij = -K_ij/2), diagonal Gamma rate beta0, and
    // shape alpha are read from the (already eligibility-checked) priors.
    const auto* diag  = static_cast<const GammaScalePrior*>(diagonal_prior_.get());
    const double sigma = slab_scale_();
    const double beta0 = diag->rate();
    const double alpha = diag->shape();
    const bool   cauchy = slab_is_cauchy_();
    const double s_ii  = suf_stat_(i, i);

    // Active neighbour set N_i (in row order). Reuse gibbs_Ni_ -- clearing a
    // std::vector keeps its capacity, so this avoids a per-row reserve(p-1).
    gibbs_Ni_.clear();
    for (size_t k = 0; k < p_; ++k) {
        if (k != i && edge_indicators_(i, k) == 1)
            gibbs_Ni_.push_back(static_cast<arma::uword>(k));
    }
    const std::vector<arma::uword>& Ni = gibbs_Ni_;
    const size_t q = Ni.size();

    // Stash old K column entries so the rank-2 update can encode the delta.
    // beta_old and the other per-row buffers below alias reused members:
    // set_size keeps the existing allocation when q fits.
    const double kii_old = precision_matrix_(i, i);
    arma::vec& beta_old = gibbs_beta_old_;
    beta_old.set_size(q);
    for (size_t k = 0; k < q; ++k) beta_old(k) = precision_matrix_(i, Ni[k]);

    // xi shape: n/2 + delta + 1 (alpha = 1). The determinant tilt |K|^delta
    // contributes delta * log(xi) to the K_ii log-kernel (|K| = |A| * xi), so
    // it shifts the Gamma shape and leaves the rate unchanged.
    const double xi_shape = static_cast<double>(n_) / 2.0 + determinant_tilt_ + 1.0;
    const double xi_rate  = (beta0 + s_ii) / 2.0;

    arma::vec beta_new(q, arma::fill::zeros);
    double kii_new;

    if (q == 0) {
        // No active neighbours: K_{i,i} = xi, no beta draw.
        kii_new = rgamma(rng_, xi_shape, xi_rate);
    } else {
        // C = (A^{-1})_{N_i, N_i} via Schur on Sigma:
        //   C_{kl} = Sigma_{N_i[k], N_i[l]} - Sigma_{N_i[k], i} Sigma_{i, N_i[l]} / Sigma_{ii}
        const double sigma_ii = covariance_matrix_(i, i);
        arma::vec& sigma_iNi = gibbs_sigma_iNi_;
        sigma_iNi.set_size(q);
        for (size_t k = 0; k < q; ++k) sigma_iNi(k) = covariance_matrix_(i, Ni[k]);
        arma::mat& C = gibbs_C_;
        C.set_size(q, q);
        for (size_t k = 0; k < q; ++k) {
            for (size_t l = 0; l < q; ++l) {
                C(k, l) = covariance_matrix_(Ni[k], Ni[l])
                          - sigma_iNi(k) * sigma_iNi(l) / sigma_ii;
            }
        }

        // M = (beta0 + S_ii) C + diag(1/(4 sigma^2 omega_k)) -- symmetric PD.
        // The slab precision on K_ij is 1/(4 sigma^2) since K_yy_ij = -K_ij/2;
        // the Cauchy scale-mixture divides it by the per-edge weight omega_k
        // (omega = 1 recovers the Normal slab).
        const double inv_4sig2 = 1.0 / (4.0 * sigma * sigma);
        arma::mat M = (beta0 + s_ii) * C;
        for (size_t k = 0; k < q; ++k) {
            const double w_k = cauchy ? omega_(i, Ni[k]) : 1.0;
            M(k, k) += inv_4sig2 / w_k;
        }

        arma::mat L_M;
        if (!arma::chol(L_M, M, "lower")) {
            // M is PD by construction (C PD as submatrix of A^{-1}, prior
            // precision > 0). A failure here means numerical trouble; skip
            // this row's update rather than corrupt K.
            return;
        }

        // S_{N_i, i} vector.
        arma::vec& s_Ni_i = gibbs_s_Ni_i_;
        s_Ni_i.set_size(q);
        for (size_t k = 0; k < q; ++k) s_Ni_i(k) = suf_stat_(Ni[k], i);

        // Mean mu = -M^{-1} S_{N_i, i}. Two triangular solves: L y = -S, L^T mu = y.
        arma::vec y  = arma::solve(arma::trimatl(L_M), -s_Ni_i);
        arma::vec mu = arma::solve(arma::trimatu(L_M.t()), y);

        // beta = mu + L^{-T} z, z ~ N(0, I_q): one triangular solve.
        arma::vec z = arma_rnorm_vec(rng_, q);
        arma::vec w = arma::solve(arma::trimatu(L_M.t()), z);
        beta_new = mu + w;

        // xi ~ Gamma; K_{i,i} = xi + beta^T C beta.
        const double xi = rgamma(rng_, xi_shape, xi_rate);
        kii_new = xi + arma::as_scalar(beta_new.t() * C * beta_new);
    }

    // Gamma shape alpha != 1: the prior carries an extra (kii/2)^(alpha-1)
    // factor that does not factorize across (beta, xi). Treat the alpha = 1
    // joint draw as an independent-MH proposal; the beta and S_ii factors
    // cancel, leaving the ratio (kii_new/kii_old)^(alpha-1). The guard keeps
    // the alpha = 1 path free of the runif draw.
    if (std::abs(alpha - 1.0) > 1e-12) {
        const double log_mh = (alpha - 1.0) * (MY_LOG(kii_new) - MY_LOG(kii_old));
        if (MY_LOG(runif(rng_)) >= log_mh) {
            return;  // reject: leave precision_matrix_, chol(K), Sigma unchanged
        }
    }

    // Write the new K column / row, then apply the symmetric rank-2 update to
    // chol(K) and Sigma. The rank-2 decomposition is
    //   Delta K = e_i vf2^T + vf2 e_i^T,
    //   (vf2)_i      = (kii_new - kii_old) / 2,
    //   (vf2)_{N_i} = beta_new - beta_old,
    //   (vf2)_k      = 0  otherwise
    // which reproduces the sparse column-change at row/col i without touching
    // the (k, l) entries off row i.
    precision_matrix_(i, i) = kii_new;
    for (size_t k = 0; k < q; ++k) {
        precision_matrix_(i, Ni[k]) = beta_new(k);
        precision_matrix_(Ni[k], i) = beta_new(k);
    }

    vf1_[i] = 1.0;
    vf2_[i] = (kii_new - kii_old) / 2.0;
    for (size_t k = 0; k < q; ++k) vf2_[Ni[k]] = beta_new(k) - beta_old(k);

    // Support of the rank-2 update: {i} + N_i. The SMW matvec gathers only
    // these columns of Sigma (O(p q) vs dense O(p^2)). Reuse gibbs_support_.
    gibbs_support_.set_size(q + 1);
    gibbs_support_[0] = static_cast<arma::uword>(i);
    for (size_t k = 0; k < q; ++k) gibbs_support_[k + 1] = Ni[k];

    // Sigma is refreshed so the next row's Schur extraction is exact;
    // chol(K) is rebuilt once in do_one_gibbs_step after the sweep.
    apply_rank2_chol_smw_update_(gibbs_support_);

    vf1_[i] = 0.0;
    vf2_[i] = 0.0;
    for (size_t k = 0; k < q; ++k) vf2_[Ni[k]] = 0.0;
}

void GGMModel::do_one_gibbs_step(int /*iteration*/) {
    // One full sweep: each column of K drawn from its full-conditional. For a
    // Cauchy slab, alternate the row sweep (conditional on omega) with a
    // refresh of the scale-mixture weights omega (conditional on K).
    for (size_t i = 0; i < p_; ++i) {
        update_row_block_gibbs(i);
    }
    // The row draws read only Sigma. One O(p^3) factorisation here restores
    // chol(K), inv chol(K), the log-det, and an exact Sigma (clearing the
    // SMW drift the sweep accumulates), ready for the edge-indicator
    // between-step that reads them next.
    refresh_cholesky();
    refresh_cauchy_omega_();
}

double GGMModel::ggm_diag_move(size_t i) {
    double logdet_omega = log_det_precision_;
    double logdet_omega_sub_ii = logdet_omega + MY_LOG(covariance_matrix_(i, i));

    size_t e = i * (i + 3) / 2; // parameter index in vectorized form (column-major upper triangle, i==j)
    double proposal_sd = proposal_sds_(e);

    double theta_curr = (logdet_omega - logdet_omega_sub_ii) / 2;
    double theta_prop = rnorm(rng_, theta_curr, proposal_sd);

    // Only the (i,i) entry of the proposal is read on the data path; the
    // full-matrix copy is needed only for the prior-only PD check below.
    if (n_ == 0) precision_proposal_ = precision_matrix_;
    precision_proposal_(i, i) = precision_matrix_(i, i) - MY_EXP(theta_curr) * MY_EXP(theta_curr) + MY_EXP(theta_prop) * MY_EXP(theta_prop);

    // Prior-only chains have no likelihood anchor vetoing non-PD proposals;
    // reject them explicitly (see proposal_is_positive_definite_). K_ii > 0
    // by construction is not sufficient: a small K_ii relative to the
    // off-diagonals can still violate PD.
    if (n_ == 0 && !proposal_is_positive_definite_()) {
        return -arma::datum::inf;
    }

    double ln_alpha = log_density_impl_diag(i);

    // Determinant-tilt prior: |K|^delta contributes delta * log_det_ratio
    // to the MH ratio. Rank-1 update => O(1) via the cached covariance.
    if (determinant_tilt_ != 0.0) {
        ln_alpha += determinant_tilt_ * log_det_ratio_diag(i);
    }

    ln_alpha += diagonal_prior_->logp(0.5 * precision_proposal_(i, i));
    ln_alpha -= diagonal_prior_->logp(0.5 * precision_matrix_(i, i));
    ln_alpha += 2.0 * (theta_prop - theta_curr); // Jacobian: dK_ii/dtheta = 2*exp(2*theta)

    if (MY_LOG(runif(rng_)) < ln_alpha) {
        double omega_ii = precision_matrix_(i, i);
        precision_matrix_(i, i) = precision_proposal_(i, i);
        cholesky_update_after_diag(omega_ii, i);
    }

    return ln_alpha;
}

double GGMModel::update_diagonal_parameter(size_t i) {
    return std::min(1.0, MY_EXP(ggm_diag_move(i)));
}

void GGMModel::cholesky_update_after_diag(double omega_ii_old, size_t i)
{
    // log|K| advances by the rank-1 determinant-lemma ratio against the
    // pre-update covariance; chol(K) is not maintained per accept (rebuilt
    // once per sweep by refresh_cholesky(), which also resets the log-det).
    log_det_precision_ += cholesky_helpers::log_det_ratio_diag_kernel(
        covariance_matrix_, i, omega_ii_old, precision_proposal_(i, i));

    // SMW rank-1 update of covariance_matrix_ = inv(K), O(p^2):
    //   K_new = K_old + alpha e_i e_i^T, alpha = K_new(i,i) - K_old(i,i)
    //   inv(K_new) = inv(K_old) - alpha / (1 + alpha Sigma_ii) c_i c_i^T,
    //     c_i = Sigma.col(i).
    // Near-singular denominator falls back to refresh_cholesky().
    double alpha = precision_proposal_(i, i) - omega_ii_old;
    arma::vec ci = covariance_matrix_.col(i);
    double denom = 1.0 + alpha * ci(i);
    if (!std::isfinite(denom) || std::abs(denom) < 1e-14) {
        refresh_cholesky();
    } else {
        // Exactly symmetric: entry (i,j) is coeff * (ci(i) * ci(j)) and IEEE
        // multiplication commutes, so no symmatu reflection is needed here
        // (unlike the rank-2 update, where b1/b2 round independently).
        covariance_matrix_ -= (alpha / denom) * (ci * ci.t());
        // Same near-singular-passage guard as the rank-2 update.
        const arma::uvec row_i = {static_cast<arma::uword>(i)};
        if (!sigma_rows_consistent_(row_i)) {
            refresh_cholesky();
        }
    }
}


double GGMModel::update_edge_indicator_parameter_pair(size_t i, size_t j) {

    size_t e = j * (j + 1) / 2 + i; // parameter index in vectorized form (column-major upper triangle)
    double proposal_sd = proposal_sds_(e);

    // Acceptance probability (raw alpha) of the birth/death proposal, set in
    // each branch. The caller derives the RB draw J and the odds accumulators.
    double alpha = 0.0;

    if (edge_indicators_(i, j) == 1) {
        // Propose to turn OFF the edge. Only the (i,j), (j,i), (j,j)
        // entries of the proposal are read below.
        precision_proposal_(i, j) = 0.0;
        precision_proposal_(j, i) = 0.0;

        // Update diagonal to preserve positive-definiteness
        get_constants(i, j);
        precision_proposal_(j, j) = constrained_diagonal(0.0);

        double ln_alpha = log_density_impl_edge(i, j);

        // Determinant-tilt prior: |K|^delta contributes delta * log_det_ratio
        // to the MH ratio. The rank-2 update at (i,j),(j,j) makes this O(p).
        if (determinant_tilt_ != 0.0) {
            ln_alpha += determinant_tilt_ * log_det_ratio_edge(i, j);
        }

        ln_alpha += MY_LOG(1.0 - inclusion_probability_(i, j)) - MY_LOG(inclusion_probability_(i, j));

        // Hierarchical spec: the delete ratio carries -log J with
        // J = Z(Gamma-)/Z(Gamma+) (the add ratio carries +log J below).
        double log_j_del = 0.0;
        if (zratio_engine_) {
            log_j_del = zratio_engine_->log_zratio(edge_indicators_,
                                                   static_cast<int>(i),
                                                   static_cast<int>(j));
            ln_alpha -= log_j_del;
        }

        ln_alpha += R::dnorm(precision_matrix_(i, j) / constants_[3], 0.0, proposal_sd, true) - MY_LOG(constants_[3]);
        // Slab in K_yy coords; proposal in K_ij coords. Jacobian |dK_yy/dK_ij| = 1/2.
        ln_alpha -= interaction_prior_->logp(-0.5 * precision_matrix_(i, j)) - MY_LOG(2.0);

        // Gamma(shape, rate) prior on changed diagonal K_jj. The Roverato move
        // slaves K_jj = c_3 + phi_{q-1,q}^2, so K_jj moves with the off-diagonal
        // and its prior must be re-evaluated.
        ln_alpha += diagonal_prior_->logp(0.5 * precision_proposal_(j, j));
        ln_alpha -= diagonal_prior_->logp(0.5 * precision_matrix_(j, j));

        if (zratio_gauge_.active) {
            zratio_gauge_record(zratio_gauge_, zratio_engine_.get(),
                                edge_indicators_, static_cast<int>(i),
                                static_cast<int>(j), ln_alpha, log_j_del, -1);
        }

        alpha = MY_EXP(std::min(0.0, ln_alpha));   // death proposal (gamma = 1)

        if (MY_LOG(runif(rng_)) < ln_alpha) {

            // Store old values for Cholesky update
            double omega_ij_old = precision_matrix_(i, j);
            double omega_jj_old = precision_matrix_(j, j);

            // Update omega
            precision_matrix_(i, j) = 0.0;
            precision_matrix_(j, i) = 0.0;
            precision_matrix_(j, j) = precision_proposal_(j, j);

            // Update edge indicator
            edge_indicators_(i, j) = 0;
            edge_indicators_(j, i) = 0;

            cholesky_update_after_edge(omega_ij_old, omega_jj_old, i, j);

            invalidate_gradient_cache();
        }

    } else {
        // Propose to turn ON the edge
        double epsilon = rnorm(rng_, 0.0, proposal_sd);

        // Get constants for current state (with edge OFF)
        get_constants(i, j);
        double omega_prop_ij = constants_[3] * epsilon;
        double omega_prop_jj = constrained_diagonal(omega_prop_ij);

        // Only the (i,j), (j,i), (j,j) entries of the proposal are read below.
        precision_proposal_(i, j) = omega_prop_ij;
        precision_proposal_(j, i) = omega_prop_ij;
        precision_proposal_(j, j) = omega_prop_jj;

        double ln_alpha = log_density_impl_edge(i, j);

        // Determinant-tilt prior: |K|^delta contributes delta * log_det_ratio
        // to the MH ratio.
        if (determinant_tilt_ != 0.0) {
            ln_alpha += determinant_tilt_ * log_det_ratio_edge(i, j);
        }

        ln_alpha += MY_LOG(inclusion_probability_(i, j)) - MY_LOG(1.0 - inclusion_probability_(i, j));

        // Hierarchical spec: the add ratio carries +log J.
        double log_j_add = 0.0;
        if (zratio_engine_) {
            log_j_add = zratio_engine_->log_zratio(edge_indicators_,
                                                   static_cast<int>(i),
                                                   static_cast<int>(j));
            ln_alpha += log_j_add;
        }

        // Slab in K_yy coords; proposal in K_ij coords. Jacobian |dK_yy/dK_ij| = 1/2.
        ln_alpha += interaction_prior_->logp(-0.5 * omega_prop_ij) - MY_LOG(2.0);

        // Gamma(shape, rate) prior on changed diagonal K_jj. The Roverato move
        // slaves K_jj = c_3 + phi_{q-1,q}^2, so K_jj moves with the off-diagonal
        // and its prior must be re-evaluated.
        ln_alpha += diagonal_prior_->logp(0.5 * precision_proposal_(j, j));
        ln_alpha -= diagonal_prior_->logp(0.5 * precision_matrix_(j, j));

        // Proposal term: proposed edge value given it was generated from truncated normal
        ln_alpha -= R::dnorm(omega_prop_ij / constants_[3], 0.0, proposal_sd, true) - MY_LOG(constants_[3]);

        if (zratio_gauge_.active) {
            zratio_gauge_record(zratio_gauge_, zratio_engine_.get(),
                                edge_indicators_, static_cast<int>(i),
                                static_cast<int>(j), ln_alpha, log_j_add, 1);
        }

        alpha = MY_EXP(std::min(0.0, ln_alpha));   // birth proposal (gamma = 0)

        if (MY_LOG(runif(rng_)) < ln_alpha) {
            // Accept: turn ON the edge
            // Store old values for Cholesky update
            double omega_ij_old = precision_matrix_(i, j);
            double omega_jj_old = precision_matrix_(j, j);

            // Update omega
            precision_matrix_(i, j) = omega_prop_ij;
            precision_matrix_(j, i) = omega_prop_ij;
            precision_matrix_(j, j) = omega_prop_jj;

            // Update edge indicator
            edge_indicators_(i, j) = 1;
            edge_indicators_(j, i) = 1;

            cholesky_update_after_edge(omega_ij_old, omega_jj_old, i, j);

            invalidate_gradient_cache();
        }
    }

    return alpha;
}

void GGMModel::do_one_metropolis_step(int iteration) {
    // Collect per-slot accept probabilities for the Robbins-Monro adapter.
    // proposal_sds_ is stored as a flat dim_-length vec indexed by the
    // upper-triangle scheme `e = j * (j + 1) / 2 + i`; we mirror that here
    // as a dim_ x 1 matrix.
    arma::mat accept_prob(dim_, 1, arma::fill::zeros);
    arma::umat index_mask(dim_, 1, arma::fill::zeros);

    // Update off-diagonals (upper triangle)
    for (size_t i = 0; i < p_ - 1; ++i) {
        for (size_t j = i + 1; j < p_; ++j) {
            double ap = update_edge_parameter(i, j);
            if (edge_indicators_(i, j) == 1) {
                size_t e = j * (j + 1) / 2 + i;
                accept_prob(e, 0) = ap;
                index_mask(e, 0) = 1;
            }
        }
    }

    // Update diagonals
    for (size_t i = 0; i < p_; ++i) {
        double ap = update_diagonal_parameter(i);
        size_t e = i * (i + 3) / 2;
        accept_prob(e, 0) = ap;
        index_mask(e, 0) = 1;
    }

    if (metropolis_adapter_) {
        metropolis_adapter_->update(index_mask, accept_prob, iteration);
    }

    // chol(K) and the log-det were deferred during the sweep (the moves read
    // Sigma and the incrementally maintained log-det only). One factorisation
    // here restores the factor, an exact Sigma, and an exact log-det for the
    // between-step and sample recording that follow.
    refresh_cholesky();
}

void GGMModel::init_metropolis_adaptation(const WarmupSchedule& schedule) {
    metropolis_adapter_ = std::make_unique<MetropolisAdaptationController>(
        proposal_sds_, schedule, target_accept_);
}

void GGMModel::prepare_iteration() {
    // Shuffle edge visit order for random-scan edge selection.
    // Called unconditionally to keep RNG state consistent.
    shuffled_edge_order_ = arma_randperm(rng_, num_pairwise_);
}

void GGMModel::update_edge_indicators() {
    for (size_t idx = 0; idx < num_pairwise_; ++idx) {
        size_t flat = shuffled_edge_order_(idx);
        size_t i = edge_pairs_(flat, 0);
        size_t j = edge_pairs_(flat, 1);
        // Capture the pre-move state before the toggle, then the acceptance
        // probability, at the row-major upper-triangle index (i = 0..p-1,
        // j = i..p-1, diagonal included) so both are aligned with
        // get_vectorized_indicator_parameters().
        const int pre = edge_indicators_(i, j);
        double alpha = use_conjugate_edge_proposal_
            ? update_edge_indicator_conjugate(i, j)
            : update_edge_indicator_parameter_pair(i, j);
        if (rb_alpha_.n_elem > 0) {
            size_t e = i * p_ - i * (i - 1) / 2 + (j - i);
            rb_alpha_(e) = alpha;
            rb_pregamma_(e) = pre;
        }
    }
    // Same rationale as the end-of-Metropolis-step refresh.
    refresh_cholesky();
}

arma::vec GGMModel::get_vectorized_rb_inclusion() {
    // J = alpha for a birth (gamma = 0), 1 - alpha for a death (gamma = 1);
    // diagonal slots (pregamma = -1) are never proposed and stay at 0.
    arma::vec j = rb_alpha_;
    for (arma::uword e = 0; e < j.n_elem; ++e) {
        if (rb_pregamma_(e) == 1) j(e) = 1.0 - rb_alpha_(e);
    }
    return j;
}

arma::vec GGMModel::get_vectorized_rb_alpha() {
    return rb_alpha_;
}

arma::ivec GGMModel::get_vectorized_rb_pregamma() {
    return rb_pregamma_;
}

double GGMModel::update_edge_indicator_conjugate(size_t i, size_t j) {
    // Full-conditional (MoMS) edge birth/death for the joint spec, Normal slab,
    // alpha = 1. The cofactor move preserves |K| (the determinant tilt and the
    // likelihood determinant cancel), so F(phi) is Gaussian in the cofactor
    // coordinate and the proposal is its exact conditional. The acceptance then
    // reduces to the inclusion odds times p_slab(0)/q(0) and does not depend on
    // the proposed value.
    get_constants(i, j);
    const double c1 = constants_[2];
    const double c2 = constants_[3];               // > 0
    const double sigma  = slab_scale_();
    const double inv_4s2 = 1.0 / (4.0 * sigma * sigma);
    const auto* diagp = static_cast<const GammaScalePrior*>(diagonal_prior_.get());
    const double beta0 = diagp->rate();
    const double alpha = diagp->shape();

    // Conditional on the Cauchy scale-mixture weight omega (= 1 for a Normal
    // slab), k_ij has slab variance 4 sigma^2 omega, so the slab precision is
    // inv_4s2 / omega. Everything else is the Normal-slab move.
    const double w = slab_is_cauchy_() ? omega_(i, j) : 1.0;
    const double inv_4s2w = inv_4s2 / w;

    // Gaussian full conditional of phi (rest of K fixed, edge present):
    //   Q  = c2^2/(4 sigma^2 omega) + beta0 + S_jj
    //   mu = -c2 (S_ij + c1/(4 sigma^2 omega)) / Q
    const double Q  = c2 * c2 * inv_4s2w + beta0 + suf_stat_(j, j);
    const double mu = -c2 * (suf_stat_(i, j) + c1 * inv_4s2w) / Q;

    // Proposal ordinate at the spike phi_0 = -c1/c2, on the k_ij scale
    // (q(0) = q_phi(phi_0)/c2), and the (conditional Gaussian) slab density of
    // k_ij at 0: N(0; 0, 4 sigma^2 omega).
    const double phi0 = -c1 / c2;
    const double d = phi0 - mu;
    const double log_q0 = 0.5 * (MY_LOG(Q) - MY_LOG(2.0 * arma::datum::pi))
                          - 0.5 * Q * d * d - MY_LOG(c2);
    const double log_pslab0 =
        -0.5 * MY_LOG(2.0 * arma::datum::pi * 4.0 * sigma * sigma * w);
    const double log_odds = MY_LOG(inclusion_probability_(i, j))
                            - MY_LOG(1.0 - inclusion_probability_(i, j));

    // A_add = odds * p_slab(0) / q(0); delete accepts with the reciprocal.
    // For a Gamma shape alpha != 1 the alpha = 1 Gaussian is an independence
    // proposal; the target's extra (k_jj/2)^(alpha-1) factor adds the
    // correction (k_jj/k_jj(0))^(alpha-1) (the 1/2 cancels in the ratio),
    // where k_jj(0) = constants_[5] is the spike diagonal.
    // Hierarchical spec: the add ratio also carries the per-graph normalizer
    // ratio J = Z(Gamma-)/Z(Gamma+), state-invariant for the toggled edge,
    // so the delete reciprocal handles it with the same value.
    double log_A_add = log_odds + log_pslab0 - log_q0;
    double log_j_conj = 0.0;
    if (zratio_engine_) {
        log_j_conj = zratio_engine_->log_zratio(edge_indicators_,
                                                static_cast<int>(i),
                                                static_cast<int>(j));
        log_A_add += log_j_conj;
    }
    const bool alpha_ne_1 = std::abs(alpha - 1.0) > 1e-12;

    // Acceptance probability (raw) of the birth/death proposal; the caller
    // derives the RB draw J and the odds accumulators. (Named accept_prob to
    // avoid the Gamma-shape `alpha` above.)
    double accept_prob = 0.0;

    if (edge_indicators_(i, j) == 0) {
        // Add: draw phi* from the full conditional, then accept.
        const double phi_star = rnorm(rng_, mu, 1.0 / std::sqrt(Q));
        const double kij = c1 + c2 * phi_star;
        const double kjj = constrained_diagonal(kij);
        double log_A = log_A_add;
        if (alpha_ne_1) {
            log_A += (alpha - 1.0) * (MY_LOG(kjj) - MY_LOG(constants_[5]));
        }
        if (zratio_gauge_.active) {
            zratio_gauge_record(zratio_gauge_, zratio_engine_.get(),
                                edge_indicators_, static_cast<int>(i),
                                static_cast<int>(j), log_A, log_j_conj, 1);
        }
        accept_prob = MY_EXP(std::min(0.0, log_A));   // birth proposal (gamma = 0)
        if (MY_LOG(runif(rng_)) < log_A) {
            const double omega_ij_old = precision_matrix_(i, j);
            const double omega_jj_old = precision_matrix_(j, j);
            precision_proposal_(i, j) = kij;
            precision_proposal_(j, j) = kjj;
            precision_matrix_(i, j) = kij;
            precision_matrix_(j, i) = kij;
            precision_matrix_(j, j) = kjj;
            edge_indicators_(i, j) = 1;
            edge_indicators_(j, i) = 1;
            cholesky_update_after_edge(omega_ij_old, omega_jj_old, i, j);
            invalidate_gradient_cache();
        }
    } else {
        // Delete: deterministic to the spike; accept with 1 / A_add evaluated
        // at the current slab state, so the alpha correction uses the current
        // diagonal against the spike diagonal.
        double log_A = -log_A_add;
        if (alpha_ne_1) {
            log_A -= (alpha - 1.0)
                     * (MY_LOG(precision_matrix_(j, j)) - MY_LOG(constants_[5]));
        }
        if (zratio_gauge_.active) {
            zratio_gauge_record(zratio_gauge_, zratio_engine_.get(),
                                edge_indicators_, static_cast<int>(i),
                                static_cast<int>(j), log_A, log_j_conj, -1);
        }
        accept_prob = MY_EXP(std::min(0.0, log_A));   // death proposal (gamma = 1)
        if (MY_LOG(runif(rng_)) < log_A) {
            const double kjj = constants_[5];      // constrained_diagonal(0)
            const double omega_ij_old = precision_matrix_(i, j);
            const double omega_jj_old = precision_matrix_(j, j);
            precision_proposal_(i, j) = 0.0;
            precision_proposal_(j, j) = kjj;
            precision_matrix_(i, j) = 0.0;
            precision_matrix_(j, i) = 0.0;
            precision_matrix_(j, j) = kjj;
            edge_indicators_(i, j) = 0;
            edge_indicators_(j, i) = 0;
            cholesky_update_after_edge(omega_ij_old, omega_jj_old, i, j);
            invalidate_gradient_cache();
        }
    }

    return accept_prob;
}

void GGMModel::tune_proposal_sd(int iteration, const WarmupSchedule& schedule) {
    auto rm_weight_opt = schedule.rm_weight_for_proposal_sd(iteration);
    if (!rm_weight_opt) return;
    const double rm_weight = *rm_weight_opt;
    const double target_accept = target_accept_;

    // Off-diagonal sweeps
    for (size_t i = 0; i < p_ - 1; ++i) {
        for (size_t j = i + 1; j < p_; ++j) {
            if (edge_indicators_(i, j) == 0) continue;

            // Same proposal/accept/update as the sampling path; only the
            // Robbins-Monro adaptation of proposal_sds_ differs (it consumes
            // the raw ln_alpha rather than the min(1,exp()) accept prob).
            size_t e = j * (j + 1) / 2 + i;
            double ln_alpha = ggm_edge_move(i, j);

            proposal_sds_(e) = update_proposal_sd_with_robbins_monro(
                proposal_sds_(e), ln_alpha, rm_weight, target_accept);
        }
    }

    // Diagonal sweeps
    for (size_t i = 0; i < p_; ++i) {
        size_t e = i * (i + 3) / 2;
        double ln_alpha = ggm_diag_move(i);

        proposal_sds_(e) = update_proposal_sd_with_robbins_monro(
            proposal_sds_(e), ln_alpha, rm_weight, target_accept);
    }

    // Invalidate gradient cache after MH updates
    invalidate_gradient_cache();

    // Same rationale as the end-of-Metropolis-step refresh.
    refresh_cholesky();
}

bool GGMModel::sigma_rows_consistent_(const arma::uvec& rows) const {
    // (Sigma K)(r, r) = 1 exactly; both matrices are symmetric, so the
    // check reads two contiguous columns per row. A violation means the
    // SMW-maintained Sigma has lost absolute accuracy (near-singular
    // passage), not that K is wrong -- K is the source of truth.
    for (arma::uword k = 0; k < rows.n_elem; ++k) {
        const arma::uword r = rows[k];
        double d = arma::dot(covariance_matrix_.col(r), precision_matrix_.col(r));
        if (!std::isfinite(d) || std::abs(d - 1.0) > kSigmaProbeTol_) {
            return false;
        }
    }
    return true;
}

void GGMModel::refresh_cholesky() {
    cholesky_of_precision_ = arma::chol(precision_matrix_, "upper");
    arma::solve(inv_cholesky_of_precision_, arma::trimatu(cholesky_of_precision_),
                arma::eye(p_, p_), arma::solve_opts::fast);
    covariance_matrix_ = inv_cholesky_of_precision_ * inv_cholesky_of_precision_.t();
    log_det_precision_ = cholesky_helpers::get_log_det(cholesky_of_precision_);
}


void GGMModel::initialize_precision_from_mle() {
    // With n=0 there is no data; keep the identity initialization.
    if (n_ == 0) return;

    // Regularized MLE: K = n * inv(S + delta * I).
    // delta = trace(S) / (p * n) gives scale-appropriate shrinkage toward I.
    double trace_s = arma::trace(suf_stat_);
    double delta = trace_s / static_cast<double>(p_ * n_);
    arma::mat S_reg = suf_stat_ + delta * arma::eye(p_, p_);
    arma::mat K_init;
    if (arma::inv_sympd(K_init, S_reg)) {
        precision_matrix_ = static_cast<double>(n_) * K_init;

        // The samplers maintain the invariant that an excluded edge has a
        // zero precision entry, so the initial state must satisfy it too:
        // zero the excluded entries whether the graph is fixed sparse or a
        // sparse initial state under edge selection.
        bool any_excluded = false;
        for (size_t i = 0; i < p_ - 1; ++i) {
            for (size_t j = i + 1; j < p_; ++j) {
                if (edge_indicators_(i, j) == 0) {
                    precision_matrix_(i, j) = 0.0;
                    precision_matrix_(j, i) = 0.0;
                    any_excluded = true;
                }
            }
        }
        if (any_excluded) {
            // Make diagonally dominant to ensure PD after zeroing.
            for (size_t i = 0; i < p_; ++i) {
                double row_sum = 0.0;
                for (size_t j = 0; j < p_; ++j) {
                    if (j != i) row_sum += std::abs(precision_matrix_(i, j));
                }
                if (precision_matrix_(i, i) <= row_sum) {
                    precision_matrix_(i, i) = row_sum + 0.1;
                }
            }
        }

        refresh_cholesky();
    }
    // If inv_sympd fails, keep the identity initialization.
}


// =============================================================================
// Missing data imputation
// =============================================================================

void GGMModel::update_suf_stat_for_imputation(int variable, int person, double delta) {
    // INVARIANT: observations_(person, variable) must still hold x_old when
    // this function is called. The loop adds 2 * delta * x_old to the (v,v)
    // entry; the delta^2 correction completes the diagonal update.
    for (size_t q = 0; q < p_; q++) {
        suf_stat_(variable, q) += delta * observations_(person, q);
        suf_stat_(q, variable) += delta * observations_(person, q);
    }
    suf_stat_(variable, variable) += delta * delta;
}

void GGMModel::impute_missing() {
    if (!has_missing_) return;

    const int num_missings = missing_index_.n_rows;

    for (int miss = 0; miss < num_missings; miss++) {
        const int person = missing_index_(miss, 0);
        const int variable = missing_index_(miss, 1);

        // Compute conditional mean: mu = -sum_{k != v} omega_{vk} * x_{ik} / omega_{vv}
        double conditional_mean = 0.0;
        for (size_t k = 0; k < p_; k++) {
            if (k != static_cast<size_t>(variable)) {
                conditional_mean += precision_matrix_(variable, k) * observations_(person, k);
            }
        }
        conditional_mean = -conditional_mean / precision_matrix_(variable, variable);

        // Conditional variance: 1 / omega_{vv}
        double conditional_sd = std::sqrt(1.0 / precision_matrix_(variable, variable));

        // Sample new value
        double x_new = rnorm(rng_, conditional_mean, conditional_sd);
        double x_old = observations_(person, variable);
        double delta = x_new - x_old;

        // Incrementally update suf_stat_ (observations_ still holds x_old)
        update_suf_stat_for_imputation(variable, person, delta);

        // Now update the observation
        observations_(person, variable) = x_new;
    }

    // Full recompute at end of sweep to eliminate floating-point drift
    // (matches OMRF pattern; cost is O(np^2), negligible for typical sizes)
    suf_stat_ = observations_.t() * observations_;
}


// =============================================================================
// Factory function
// =============================================================================

GGMModel createGGMModelFromR(
    const Rcpp::List& inputFromR,
    const arma::mat& prior_inclusion_prob,
    const arma::imat& initial_edge_indicators,
    const bool edge_selection,
    std::unique_ptr<BaseParameterPrior> interaction_prior,
    std::unique_ptr<BaseParameterPrior> diagonal_prior,
    const bool na_impute
) {

    if (inputFromR.containsElementNamed("n") && inputFromR.containsElementNamed("suf_stat")) {
        int n = Rcpp::as<int>(inputFromR["n"]);
        arma::mat suf_stat = Rcpp::as<arma::mat>(inputFromR["suf_stat"]);
        return GGMModel(
            n,
            suf_stat,
            prior_inclusion_prob,
            initial_edge_indicators,
            edge_selection,
            std::move(interaction_prior),
            std::move(diagonal_prior)
        );
    } else if (inputFromR.containsElementNamed("X")) {
        arma::mat X = Rcpp::as<arma::mat>(inputFromR["X"]);
        return GGMModel(
            X,
            prior_inclusion_prob,
            initial_edge_indicators,
            edge_selection,
            std::move(interaction_prior),
            std::move(diagonal_prior),
            na_impute
        );
    } else {
        throw std::invalid_argument("Input list must contain either 'X' or both 'n' and 'suf_stat'.");
    }

}

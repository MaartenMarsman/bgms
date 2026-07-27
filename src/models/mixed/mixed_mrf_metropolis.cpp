// mixed_mrf_metropolis.cpp — MixedMRFModel within-model MH updates.
//
// The per-parameter Metropolis-Hastings sweeps (main effects, discrete/cross/continuous
// pairwise effects, edge indicators, continuous mean) and the rank-2 precision math they
// rely on: get_precision_constants, log_det_ratio_yy_*, log_ggm_ratio_*, and the
// cholesky_update_after_precision_* helpers. The load-bearing std::abs sign convention
// in the det-ratio helpers must be preserved (see architecture.md, Numerical Considerations).
#include <RcppArmadillo.h>
#include <utility>
#include "models/mixed/mixed_mrf_model.h"
#include "rng/rng_utils.h"
#include "mcmc/execution/step_result.h"
#include "math/explog_macros.h"


// =============================================================================
// update_main_effect
// =============================================================================
// MH update for one main-effect parameter.
//   Ordinal: main_effects_discrete_(s, c) = threshold for category c+1  (c in [0, C_s-1])
//   Blume-Capel: main_effects_discrete_(s, 0) = linear α, main_effects_discrete_(s, 1) = quadratic β
//                (c indexes 0 or 1 for BC)
//
// The accept/reject uses log_marginal_omrf(s) + beta-type prior.
// =============================================================================

double MixedMRFModel::update_main_effect(int s, int c, std::optional<double> rm_weight) {
    double& current = main_effects_discrete_(s, c);
    double proposal_sd = proposal_sd_main_discrete_(s, c);

    double current_val = current;
    double proposed = rnorm(rng_, current_val, proposal_sd);

    // Current log-posterior (likelihood term from the sweep cache)
    double ll_curr = ll_marginal_cache_(s)
                   + threshold_prior_->logp(current_val);

    // Proposed log-posterior
    current = proposed;
    double marginal_prop = log_marginal_omrf_cached(s);
    double ll_prop = marginal_prop + threshold_prior_->logp(proposed);

    double ln_alpha = ll_prop - ll_curr;

    if(MY_LOG(runif(rng_)) >= ln_alpha) {
        current = current_val;  // reject
    } else {
        ll_marginal_cache_(s) = marginal_prop;
    }

    if (rm_weight) {
        proposal_sd_main_discrete_(s, c) = update_proposal_sd_with_robbins_monro(
            proposal_sd_main_discrete_(s, c), ln_alpha, *rm_weight, target_accept_);
    }
    return std::min(1.0, MY_EXP(ln_alpha));
}


// =============================================================================
// update_continuous_mean
// =============================================================================
// MH update for one continuous mean parameter main_effects_continuous_(j).
// The accept/reject uses the cached GGM value plus a rank-1 quadratic-form
// delta, and a Normal(0, 1) prior. conditional_mean_ is only touched on
// accept (column-j shift by delta).
// =============================================================================

double MixedMRFModel::update_continuous_mean(int j, std::optional<double> rm_weight) {
    double current_val = main_effects_continuous_(j);
    double proposed = rnorm(rng_, current_val, proposal_sd_main_continuous_(j));
    double delta = proposed - current_val;

    // Current log-posterior (likelihood terms from the sweep caches)
    double ll_curr = ll_ggm_cache_ + arma::accu(ll_marginal_cache_)
                   + means_prior_->logp(current_val);

    // Proposed state: μ_j moves column j of the conditional mean and the
    // rest-score offsets 2 A_xy μ; everything else is unchanged.
    //
    // GGM part: ΔM = delta · 1 e_j' with K unchanged, so
    //   quad_prop - quad_curr = -2 delta (1'D) K[:,j] + n delta² K_jj,
    // read off the residual column sums (D = Y - conditional mean).
    arma::rowvec resid_colsum = arma::sum(continuous_observations_, 0)
                              - arma::sum(conditional_mean_, 0);
    double quad_delta =
        -2.0 * delta * arma::dot(resid_colsum,
                                 -2.0 * pairwise_effects_continuous_.col(j))
        + static_cast<double>(n_) * delta * delta
              * (-2.0 * pairwise_effects_continuous_(j, j));
    double ggm_prop = ll_ggm_cache_ - quad_delta / 2.0;

    main_effects_continuous_(j) = proposed;
    cross_bias_prop_ = cross_bias_ + (2.0 * delta) * pairwise_effects_cross_.col(j);

    for(size_t s = 0; s < p_; ++s)
        ll_marginal_prop_(s) = log_marginal_omrf_from(
            s, marginal_matvec_, marginal_interactions_(s, s), cross_bias_prop_(s));

    double ll_prop = ggm_prop + arma::accu(ll_marginal_prop_)
                   + means_prior_->logp(proposed);

    double ln_alpha = ll_prop - ll_curr;

    if(MY_LOG(runif(rng_)) >= ln_alpha) {
        main_effects_continuous_(j) = current_val;  // reject
    } else {
        conditional_mean_.col(j) += delta;
        ll_ggm_cache_ = ggm_prop;
        ll_marginal_cache_ = ll_marginal_prop_;
        cross_bias_ = cross_bias_prop_;
    }

    if (rm_weight) {
        proposal_sd_main_continuous_(j) = update_proposal_sd_with_robbins_monro(
            proposal_sd_main_continuous_(j), ln_alpha, *rm_weight, target_accept_);
    }
    return std::min(1.0, MY_EXP(ln_alpha));
}


// =============================================================================
// update_pairwise_discrete
// =============================================================================
// MH update for one discrete-discrete interaction pairwise_effects_discrete_(i, j).
// Symmetric: sets both (i,j) and (j,i).
// Acceptance: log_marginal_omrf(i) + log_marginal_omrf(j) + Cauchy prior.
// =============================================================================

double MixedMRFModel::update_pairwise_discrete(int i, int j, std::optional<double> rm_weight) {
    double current_val = pairwise_effects_discrete_(i, j);
    double proposed = rnorm(rng_, current_val, proposal_sd_pairwise_discrete_(i, j));
    double delta = proposed - current_val;

    // Current log-posterior (likelihood terms from the sweep cache)
    double ll_curr = ll_marginal_cache_(i) + ll_marginal_cache_(j)
                   + interaction_prior_->logp(current_val);

    // An A_xx(i,j) change moves one entry of M, so the cached X · M shifts
    // by delta in columns i and j only.
    pairwise_effects_discrete_(i, j) = proposed;
    pairwise_effects_discrete_(j, i) = proposed;
    refresh_marginal_interactions_entry(i, j);
    matvec_col_i_scratch_ = marginal_matvec_.col(i);
    matvec_col_j_scratch_ = marginal_matvec_.col(j);
    marginal_matvec_.col(i) += delta * discrete_observations_dbl_.col(j);
    marginal_matvec_.col(j) += delta * discrete_observations_dbl_.col(i);

    double marginal_prop_i = log_marginal_omrf_cached(i);
    double marginal_prop_j = log_marginal_omrf_cached(j);
    double ll_prop = marginal_prop_i + marginal_prop_j
                   + interaction_prior_->logp(proposed);

    double ln_alpha = ll_prop - ll_curr;

    if(MY_LOG(runif(rng_)) >= ln_alpha) {
        pairwise_effects_discrete_(i, j) = current_val;  // reject
        pairwise_effects_discrete_(j, i) = current_val;
        refresh_marginal_interactions_entry(i, j);
        marginal_matvec_.col(i) = matvec_col_i_scratch_;
        marginal_matvec_.col(j) = matvec_col_j_scratch_;
    } else {
        ll_marginal_cache_(i) = marginal_prop_i;
        ll_marginal_cache_(j) = marginal_prop_j;
    }

    if (rm_weight) {
        proposal_sd_pairwise_discrete_(i, j) = update_proposal_sd_with_robbins_monro(
            proposal_sd_pairwise_discrete_(i, j), ln_alpha, *rm_weight, target_accept_);
    }
    return std::min(1.0, MY_EXP(ln_alpha));
}


// =============================================================================
// Rank-1 precision proposal helpers (permutation-free)
// =============================================================================
// Direct analogs of GGMModel::get_constants / constrained_diagonal,
// operating on precision = -2 * pairwise_effects_continuous_,
// cholesky_of_precision_, and covariance_continuous_.
//
// All constants and proposals live in precision space. Conversion to/from
// pairwise_effects_continuous_ happens at the outer call sites.
// =============================================================================

void MixedMRFModel::get_precision_constants(int i, int j) {
    // Kyy = -2 * pairwise_effects_continuous_; delegate to the shared kernel,
    // passing the resolved precision entries.
    size_t ui = static_cast<size_t>(i);
    size_t uj = static_cast<size_t>(j);
    cont_constants_ = cholesky_helpers::precision_proposal_constants(
        log_det_precision_, covariance_continuous_, ui, uj,
        -2.0 * pairwise_effects_continuous_(ui, uj),
        -2.0 * pairwise_effects_continuous_(uj, uj));
}

double MixedMRFModel::precision_constrained_diagonal(double x) const {
    if(x == 0.0) {
        return cont_constants_[5];
    } else {
        double t = (x - cont_constants_[2]) / cont_constants_[3];
        return cont_constants_[4] + t * t;
    }
}


// =============================================================================
// log_ggm_ratio_edge
// =============================================================================
// Log-likelihood ratio for a rank-2 off-diagonal precision change: matrix
// determinant lemma for the log-det part, Woodbury for the proposed
// covariance, and a rank-2 expansion of the quadratic-form difference.
// With ΔK = vf1 vf2' + vf2 vf1' and ΔΣ = -(w1 s2' + w2 s1'), the conditional
// mean moves by ΔM = a1 s2' + a2 s1' (a_k = -2 X A_xy w_k) and
//   quad_prop - quad_curr = -2 tr(ΔM' D K) + tr(K ΔM' ΔM) + tr(ΔK D_p' D_p),
// D = Y - conditional mean, D_p = D - ΔM: O(nq + q²) dot products instead of
// the O(nq²) from-scratch quadratic forms. Assumes precision_proposal_ is
// filled.
// =============================================================================

double MixedMRFModel::log_det_ratio_yy_edge(int i, int j) const {
    // Kyy = -2 * pairwise_effects_continuous_; delegate to the shared kernel.
    size_t ui = static_cast<size_t>(i);
    size_t uj = static_cast<size_t>(j);
    return cholesky_helpers::log_det_ratio_edge_kernel(
        covariance_continuous_, ui, uj,
        -2.0 * pairwise_effects_continuous_(ui, uj), precision_proposal_(ui, uj),
        -2.0 * pairwise_effects_continuous_(uj, uj), precision_proposal_(uj, uj));
}

double MixedMRFModel::log_det_ratio_yy_diag(int i) const {
    size_t ui = static_cast<size_t>(i);
    return cholesky_helpers::log_det_ratio_diag_kernel(
        covariance_continuous_, ui,
        -2.0 * pairwise_effects_continuous_(ui, ui), precision_proposal_(ui, ui));
}


double MixedMRFModel::log_ggm_ratio_edge(int i, int j, arma::mat& cov_prop_out) const {
    size_t ui = static_cast<size_t>(i);
    size_t uj = static_cast<size_t>(j);

    // --- Log-determinant ratio via matrix determinant lemma ---
    // ΔΩ has 3 nonzero entries: (i,j), (j,i), (j,j).
    // Ui = old - new off-diag, Uj = (old - new diag) / 2. The same Ui/Uj also
    // drive the Woodbury covariance update below, so they are kept here; the
    // log-det ratio itself is the canonical rank-2 det-lemma in
    // log_det_ratio_yy_edge (recomputes the identical Ui/Uj internally).
    double Ui = -2.0 * pairwise_effects_continuous_(ui, uj) - precision_proposal_(ui, uj);
    double Uj = (-2.0 * pairwise_effects_continuous_(uj, uj) - precision_proposal_(uj, uj)) / 2.0;

    double logdet_ratio = log_det_ratio_yy_edge(i, j);

    // --- Proposed covariance via Woodbury ---
    // ΔΩ = vf1 vf2' + vf2 vf1' where vf1 = [0,...,-1,...] (j-th),
    //   vf2 = [0,...,Ui,...,Uj,...] (i-th and j-th).
    // s1 = Σ vf1 = -Σ[:,j], s2 = Σ vf2 = Ui*Σ[:,i] + Uj*Σ[:,j]
    cont_s1_ = -covariance_continuous_.col(uj);
    cont_s2_ = Ui * covariance_continuous_.col(ui) + Uj * covariance_continuous_.col(uj);
    const arma::vec& s1 = cont_s1_;
    const arma::vec& s2 = cont_s2_;

    // 2×2 core matrix T = I + [vf2,vf1]' [s1,s2]
    // T = [1 + vf2's1,  vf2's2;  vf1's1,  1 + vf1's2]
    double t11 = 1.0 + Ui * s1(ui) + Uj * s1(uj);     // 1 + vf2' s1
    double t12 = Ui * s2(ui) + Uj * s2(uj);            // vf2' s2
    double t21 = -s1(uj);                              // vf1' s1 = Σ(j,j)
    double t22 = 1.0 - s2(uj);                         // 1 + vf1' s2

    double det_T = t11 * t22 - t12 * t21;

    // T^{-1}
    double inv_t11 =  t22 / det_T;
    double inv_t12 = -t12 / det_T;
    double inv_t21 = -t21 / det_T;
    double inv_t22 =  t11 / det_T;

    // Σ' = Σ - [s1,s2] T^{-1} [s2',s1']
    //     = Σ - (inv_t11*s1 + inv_t21*s2)*s2' - (inv_t12*s1 + inv_t22*s2)*s1'
    arma::vec w1 = inv_t11 * s1 + inv_t21 * s2;  // coefficient for s2' row
    arma::vec w2 = inv_t12 * s1 + inv_t22 * s2;  // coefficient for s1' row
    arma::mat cov_prop = covariance_continuous_ - w1 * s2.t() - w2 * s1.t();

    // --- Quadratic-form difference through the rank-2 structure ---
    // ΔΣ = -(w1 s2' + w2 s1'), so ΔM = 2 X A_xy ΔΣ = a1 s2' + a2 s1' with
    // a_k = -2 (X A_xy) w_k read off the sweep cache.
    cont_a1_ = -2.0 * (cross_matvec_ * w1);
    cont_a2_ = -2.0 * (cross_matvec_ * w2);
    resid_scratch_ = continuous_observations_ - conditional_mean_;
    const arma::mat& D = resid_scratch_;

    // K s and D' a contractions (K = -2 A_yy)
    arma::vec k1 = -2.0 * (pairwise_effects_continuous_ * s1);
    arma::vec k2 = -2.0 * (pairwise_effects_continuous_ * s2);
    arma::vec g1 = D.t() * cont_a1_;
    arma::vec g2 = D.t() * cont_a2_;

    // -2 tr(ΔM' D K)
    double quad_lin = -2.0 * (arma::dot(g1, k2) + arma::dot(g2, k1));

    // tr(K ΔM' ΔM)
    double a11 = arma::dot(cont_a1_, cont_a1_);
    double a12 = arma::dot(cont_a1_, cont_a2_);
    double a22 = arma::dot(cont_a2_, cont_a2_);
    double quad_sq = a11 * arma::dot(s2, k2)
                   + a12 * (arma::dot(s1, k2) + arma::dot(s2, k1))
                   + a22 * arma::dot(s1, k1);

    // tr(ΔΩ D_p' D_p) = 2 (D_p vf2)'(D_p vf1) with D_p = D - ΔM
    double s2_vf2 = Ui * s2(ui) + Uj * s2(uj);
    double s1_vf2 = Ui * s1(ui) + Uj * s1(uj);
    arma::vec dp_vf1 = -D.col(uj) + s2(uj) * cont_a1_ + s1(uj) * cont_a2_;
    arma::vec dp_vf2 = Ui * D.col(ui) + Uj * D.col(uj)
                     - s2_vf2 * cont_a1_ - s1_vf2 * cont_a2_;
    double quad_dk = 2.0 * arma::dot(dp_vf2, dp_vf1);

    double quad_delta = quad_lin + quad_sq + quad_dk;

    double n = static_cast<double>(n_);
    cov_prop_out = std::move(cov_prop);
    return n / 2.0 * logdet_ratio - quad_delta / 2.0;
}


// =============================================================================
// log_ggm_ratio_diag
// =============================================================================
// Log-likelihood ratio for a rank-1 diagonal precision change.
// Same structure as log_ggm_ratio_edge but simpler (Ui = 0): the
// quadratic-form difference contracts through ΔΣ = c s s' in O(nq).
// =============================================================================

double MixedMRFModel::log_ggm_ratio_diag(int i, arma::mat& cov_prop_out) const {
    size_t ui = static_cast<size_t>(i);

    // Current precision diagonal
    double precision_ii = -2.0 * pairwise_effects_continuous_(ui, ui);

    // --- Log-determinant ratio (rank-1) ---
    // Uj also drives the Sherman-Morrison covariance update below, so it is
    // kept here; the log-det ratio itself is the canonical rank-1 det-lemma in
    // log_det_ratio_yy_diag (recomputes the identical Uj internally).
    double Uj = (precision_ii - precision_proposal_(ui, ui)) / 2.0;

    double logdet_ratio = log_det_ratio_yy_diag(i);

    // --- Proposed covariance via Sherman-Morrison (rank-1 special case) ---
    // ΔΩ = -2Uj * e_i e_i', so Σ' = Σ + 2Uj * Σ[:,i] Σ[i,:]' / (1 - 2Uj * Σ(i,i))
    cont_s1_ = covariance_continuous_.col(ui);
    const arma::vec& s = cont_s1_;
    double denom = 1.0 - 2.0 * Uj * covariance_continuous_(ui, ui);
    double coef = 2.0 * Uj / denom;
    arma::mat cov_prop = covariance_continuous_ + coef * s * s.t();

    // --- Quadratic-form difference through the rank-1 structure ---
    // ΔΣ = coef * s s', so ΔM = 2 X A_xy ΔΣ = a s' with a = 2 coef (X A_xy) s;
    //   quad_prop - quad_curr
    //     = -2 (D'a)·(K s) + (a·a) (s'K s) - 2Uj ||D_p[:,i]||²,
    // D = Y - conditional mean, D_p[:,i] = D[:,i] - s(i) a.
    cont_a1_ = (2.0 * coef) * (cross_matvec_ * s);
    resid_scratch_ = continuous_observations_ - conditional_mean_;
    const arma::mat& D = resid_scratch_;

    arma::vec k = -2.0 * (pairwise_effects_continuous_ * s);
    arma::vec g = D.t() * cont_a1_;

    double quad_lin = -2.0 * arma::dot(g, k);
    double quad_sq = arma::dot(cont_a1_, cont_a1_) * arma::dot(s, k);
    arma::vec dp_col = D.col(ui) - s(ui) * cont_a1_;
    double quad_dk = -2.0 * Uj * arma::dot(dp_col, dp_col);

    double quad_delta = quad_lin + quad_sq + quad_dk;

    double n = static_cast<double>(n_);
    cov_prop_out = std::move(cov_prop);
    return n / 2.0 * logdet_ratio - quad_delta / 2.0;
}


// =============================================================================
// cholesky_update_after_precision_edge
// =============================================================================
// Rank-2 Cholesky update after accepting an off-diagonal precision change.
// Decomposes ΔΩ = vf1*vf2' + vf2*vf1' into two rank-1 ops.
// Then recomputes inv_cholesky_of_precision_ and covariance_continuous_.
// =============================================================================

bool MixedMRFModel::cholesky_update_after_precision_edge(
    double old_ij, double old_jj, int i, int j)
{
    cont_v2_[0] = old_ij - precision_proposal_(i, j);
    cont_v2_[1] = (old_jj - precision_proposal_(j, j)) / 2.0;

    cont_vf1_[i] = cont_v1_[0];   // 0
    cont_vf1_[j] = cont_v1_[1];   // -1
    cont_vf2_[i] = cont_v2_[0];
    cont_vf2_[j] = cont_v2_[1];

    cont_u1_ = (cont_vf1_ + cont_vf2_) / std::sqrt(2.0);
    cont_u2_ = (cont_vf1_ - cont_vf2_) / std::sqrt(2.0);

    cholesky_update(cholesky_of_precision_, cont_u1_);
    bool down_ok = cholesky_downdate(cholesky_of_precision_, cont_u2_);

    // Update the inverse Cholesky; if the downdate lost positive definiteness
    // or the rank-2 update has drifted into ill-conditioning, rebuild the
    // decomposition from scratch (mirrors GGMModel's drift-guard).
    // pairwise_effects_continuous_ already holds the accepted value here, so
    // the rebuild reconstructs the accepted state.
    bool incremental = down_ok &&
        arma::inv(inv_cholesky_of_precision_, arma::trimatu(cholesky_of_precision_));
    if (incremental) {
        covariance_continuous_ = inv_cholesky_of_precision_ * inv_cholesky_of_precision_.t();
        log_det_precision_ = cholesky_helpers::get_log_det(cholesky_of_precision_);
    } else {
        recompute_pairwise_effects_continuous_decomposition();
    }

    cont_vf1_[i] = 0.0;
    cont_vf1_[j] = 0.0;
    cont_vf2_[i] = 0.0;
    cont_vf2_[j] = 0.0;

    return incremental;
}


// =============================================================================
// cholesky_update_after_precision_diag
// =============================================================================
// Rank-1 Cholesky update after accepting a diagonal precision change.
// =============================================================================

bool MixedMRFModel::cholesky_update_after_precision_diag(double old_ii, int i) {
    double delta = old_ii - precision_proposal_(i, i);
    bool downdate = delta > 0.0;

    cont_vf1_[i] = std::sqrt(std::abs(delta));

    bool down_ok = true;
    if(downdate)
        down_ok = cholesky_downdate(cholesky_of_precision_, cont_vf1_);
    else
        cholesky_update(cholesky_of_precision_, cont_vf1_);

    // Update the inverse Cholesky; fall back to a full rebuild if the
    // downdate lost positive definiteness or on drift (mirrors GGMModel's
    // drift-guard).
    bool incremental = down_ok &&
        arma::inv(inv_cholesky_of_precision_, arma::trimatu(cholesky_of_precision_));
    if (incremental) {
        covariance_continuous_ = inv_cholesky_of_precision_ * inv_cholesky_of_precision_.t();
        log_det_precision_ = cholesky_helpers::get_log_det(cholesky_of_precision_);
    } else {
        recompute_pairwise_effects_continuous_decomposition();
    }

    cont_vf1_[i] = 0.0;

    return incremental;
}


// =============================================================================
// update_pairwise_effects_continuous_offdiag
// =============================================================================
// MH update for one off-diagonal element of the precision matrix pairwise_effects_continuous_(i, j).
// Uses rank-1 Cholesky infrastructure (GGM-style, no permutation):
//   1. Extract constants from covariance_continuous_ and cholesky_of_precision_
//   2. Propose on the unconstrained Cholesky scale
//   3. Map to precision space with constrained diagonal
//   4. Evaluate rank-2 log-likelihood ratio
//   5. On accept: rank-1 Cholesky update
//
// Constants and proposals live in precision space.
// Prior: interaction_prior_ on off-diagonal,
//        diagonal_prior_ on negative diagonal.
// Storage: pairwise_effects_continuous_ = -1/2 * precision.
// =============================================================================

double MixedMRFModel::update_pairwise_effects_continuous_offdiag(int i, int j, std::optional<double> rm_weight) {
    get_precision_constants(i, j);

    double phi_curr = cont_constants_[0];  // Phi_q1q
    double phi_prop = rnorm(rng_, phi_curr, proposal_sd_pairwise_continuous_(i, j));

    // Propose in precision space
    double theta_prop_ij = cont_constants_[2] + cont_constants_[3] * phi_prop;
    double theta_prop_jj = precision_constrained_diagonal(theta_prop_ij);

    // Current precision values
    double theta_curr_ij = -2.0 * pairwise_effects_continuous_(i, j);
    double theta_curr_jj = -2.0 * pairwise_effects_continuous_(j, j);

    // Fill proposal matrix in precision space
    precision_proposal_ = -2.0 * pairwise_effects_continuous_;
    precision_proposal_(i, j) = theta_prop_ij;
    precision_proposal_(j, i) = theta_prop_ij;
    precision_proposal_(j, j) = theta_prop_jj;

    arma::mat cov_prop;
    double ggm_ratio = log_ggm_ratio_edge(i, j, cov_prop);
    double ln_alpha = ggm_ratio;

    // Determinant-tilt prior on the Kyy block: |Kyy|^delta contributes
    //   delta * (log|Kyy_prop| - log|Kyy_curr|)
    // to the MH ratio. Rank-2 lemma, O(q) via the cached covariance.
    if (determinant_tilt_yy_ != 0.0) {
        ln_alpha += determinant_tilt_yy_ * log_det_ratio_yy_edge(i, j);
    }

    // OMRF ratio with proposed continuous interactions. The marginal
    // interactions M = A_xx + 2 A_xy Σ A_xy^T depend on Σ; when Kyy changes,
    // Σ changes too, so the proposed-state marginals use Σ' (cov_prop).
    ln_alpha += omrf_ratio_for_covariance_change(cov_prop);

    // Prior ratio:
    //   Cauchy(0, scale) on off-diagonal = -1/2 * precision_ij
    //   Gamma(1, 1) on negative diagonal = 1/2 * precision_jj
    double cont_prop_ij = -0.5 * theta_prop_ij;
    double cont_curr_ij = pairwise_effects_continuous_(i, j);

    ln_alpha += interaction_prior_->logp(cont_prop_ij);
    ln_alpha -= interaction_prior_->logp(cont_curr_ij);

    // Gamma(1,1) prior on changed diagonal K_jj
    ln_alpha += diagonal_prior_->logp(0.5 * theta_prop_jj);
    ln_alpha -= diagonal_prior_->logp(0.5 * (-2.0 * pairwise_effects_continuous_(j, j)));

    if(MY_LOG(runif(rng_)) < ln_alpha) {
        // Pass old precision values to Cholesky update
        double old_theta_ij = theta_curr_ij;
        double old_theta_jj = theta_curr_jj;

        // Store: pairwise_effects_continuous_ = -1/2 * precision
        pairwise_effects_continuous_(i, j) = -0.5 * theta_prop_ij;
        pairwise_effects_continuous_(j, i) = -0.5 * theta_prop_ij;
        pairwise_effects_continuous_(j, j) = -0.5 * theta_prop_jj;

        if (cholesky_update_after_precision_edge(old_theta_ij, old_theta_jj, i, j)) {
            adopt_kyy_proposal_caches(ggm_ratio, /*rank2=*/true);
        } else {
            recompute_am_caches();
        }
    }

    if (rm_weight) {
        proposal_sd_pairwise_continuous_(i, j) = update_proposal_sd_with_robbins_monro(
            proposal_sd_pairwise_continuous_(i, j), ln_alpha, *rm_weight, target_accept_);
    }
    return std::min(1.0, MY_EXP(ln_alpha));
}


// =============================================================================
// update_pairwise_effects_continuous_diag
// =============================================================================
// MH update for one diagonal element of the precision matrix.
// Proposes on the log-Cholesky scale to ensure positivity of precision.
// Uses rank-1 Cholesky update on accept.
// Prior: Gamma(1, 1) on negative diagonal + Jacobian for log-scale proposal.
// =============================================================================

double MixedMRFModel::update_pairwise_effects_continuous_diag(int i, std::optional<double> rm_weight) {
    double logdet = log_det_precision_;
    double logdet_sub_ii = logdet + MY_LOG(covariance_continuous_(i, i));

    double theta_curr = (logdet - logdet_sub_ii) / 2.0;
    double theta_prop = rnorm(rng_, theta_curr, proposal_sd_pairwise_continuous_(i, i));

    // Current precision diagonal
    double theta_ii_curr = -2.0 * pairwise_effects_continuous_(i, i);
    double theta_ii_prop = theta_ii_curr
        - MY_EXP(theta_curr) * MY_EXP(theta_curr)
        + MY_EXP(theta_prop) * MY_EXP(theta_prop);

    // Fill proposal in precision space
    precision_proposal_ = -2.0 * pairwise_effects_continuous_;
    precision_proposal_(i, i) = theta_ii_prop;

    arma::mat cov_prop;
    double ggm_ratio = log_ggm_ratio_diag(i, cov_prop);
    double ln_alpha = ggm_ratio;

    // Determinant-tilt prior: rank-1 lemma, O(1) via the cached covariance.
    if (determinant_tilt_yy_ != 0.0) {
        ln_alpha += determinant_tilt_yy_ * log_det_ratio_yy_diag(i);
    }

    // OMRF ratio with proposed continuous interactions. Use Σ' (cov_prop) for
    // the proposed-state marginals, not the cached Σ.
    ln_alpha += omrf_ratio_for_covariance_change(cov_prop);

    // Prior ratio: Gamma(1,1) on K_ii (precision diagonal)
    ln_alpha += diagonal_prior_->logp(0.5 * theta_ii_prop);
    ln_alpha -= diagonal_prior_->logp(0.5 * theta_ii_curr);

    // Jacobian: dK_ii/dtheta = 2*exp(2*theta)
    ln_alpha += 2.0 * (theta_prop - theta_curr);

    if(MY_LOG(runif(rng_)) < ln_alpha) {
        // Pass old precision value to Cholesky update
        double old_theta_ii = theta_ii_curr;

        // Store: pairwise_effects_continuous_ = -1/2 * precision
        pairwise_effects_continuous_(i, i) = -0.5 * theta_ii_prop;

        if (cholesky_update_after_precision_diag(old_theta_ii, i)) {
            adopt_kyy_proposal_caches(ggm_ratio, /*rank2=*/false);
        } else {
            recompute_am_caches();
        }
    }

    if (rm_weight) {
        proposal_sd_pairwise_continuous_(i, i) = update_proposal_sd_with_robbins_monro(
            proposal_sd_pairwise_continuous_(i, i), ln_alpha, *rm_weight, target_accept_);
    }
    return std::min(1.0, MY_EXP(ln_alpha));
}


// =============================================================================
// update_pairwise_cross
// =============================================================================
// MH update for one cross-type interaction pairwise_effects_cross_(i, j).
// Acceptance: sum_s log_marginal_omrf(s) + a rank-1 GGM quadratic-form delta
// + Cauchy prior. Parameters and conditional_mean_ are only touched on accept.
// =============================================================================

double MixedMRFModel::update_pairwise_cross(int i, int j, std::optional<double> rm_weight) {
    double current_val = pairwise_effects_cross_(i, j);
    double proposed = rnorm(rng_, current_val, proposal_sd_pairwise_cross_(i, j));
    double delta = proposed - current_val;

    // Current log-posterior (likelihood terms from the sweep caches)
    double ll_curr = ll_ggm_cache_ + arma::accu(ll_marginal_cache_)
                   + interaction_prior_->logp(current_val);

    // An A_xy(i,j) change of delta moves M by the rank-2 update
    //   ΔM = 2 delta (e_i u' + u e_i') + 2 delta² Σ_jj e_i e_i',  u = A_xy Σ.col(j)
    // so X · M' follows from the cached matvecs without rebuilding M, and the
    // conditional mean shifts by the rank-1 term 2 delta x_i Σ.row(j).
    arma::vec u = pairwise_effects_cross_ * covariance_continuous_.col(j);
    marginal_matvec_prop_ = marginal_matvec_
        + (2.0 * delta) * discrete_observations_dbl_.col(i) * u.t();
    marginal_matvec_prop_.col(i) +=
        (2.0 * delta) * (cross_matvec_ * covariance_continuous_.col(j))
        + (2.0 * delta * delta * covariance_continuous_(j, j))
              * discrete_observations_dbl_.col(i);
    mdiag_prop_ = marginal_interactions_.diag();
    mdiag_prop_(i) += 4.0 * delta * u(i)
                    + 2.0 * delta * delta * covariance_continuous_(j, j);
    cross_bias_prop_ = cross_bias_;
    cross_bias_prop_(i) += 2.0 * delta * main_effects_continuous_(j);

    // GGM part: the conditional-mean shift is ΔM_y = a v' with
    // a = 2 delta x_i, v = Σ[:,j]; K is unchanged, so
    //   quad_prop - quad_curr = -2 (D'a)·(K v) + (a·a) (v'K v),
    // D = Y - conditional mean.
    arma::vec a = (2.0 * delta) * discrete_observations_dbl_.col(i);
    arma::vec kv = -2.0 * (pairwise_effects_continuous_ * covariance_continuous_.col(j));
    resid_scratch_ = continuous_observations_ - conditional_mean_;
    arma::vec g = resid_scratch_.t() * a;
    double quad_delta = -2.0 * arma::dot(g, kv)
        + arma::dot(a, a) * arma::dot(covariance_continuous_.col(j), kv);
    double ggm_prop = ll_ggm_cache_ - quad_delta / 2.0;

    for(size_t s = 0; s < p_; ++s)
        ll_marginal_prop_(s) = log_marginal_omrf_from(
            s, marginal_matvec_prop_, mdiag_prop_(s), cross_bias_prop_(s));

    double ll_prop = ggm_prop + arma::accu(ll_marginal_prop_)
                   + interaction_prior_->logp(proposed);

    double ln_alpha = ll_prop - ll_curr;

    if(MY_LOG(runif(rng_)) < ln_alpha) {
        pairwise_effects_cross_(i, j) = proposed;
        adopt_cross_proposal_caches(i, j, delta, u, ggm_prop);
    }

    if (rm_weight) {
        proposal_sd_pairwise_cross_(i, j) = update_proposal_sd_with_robbins_monro(
            proposal_sd_pairwise_cross_(i, j), ln_alpha, *rm_weight, target_accept_);
    }
    return std::min(1.0, MY_EXP(ln_alpha));
}


// =============================================================================
// update_edge_indicator_discrete
// =============================================================================
// Metropolis-Hastings add-delete move for a discrete-discrete edge (i, j).
//   Add (G=0→1): propose k ~ N(0, σ), accept with slab + Hastings.
//   Delete (G=1→0): set k = 0, accept with reverse terms.
// =============================================================================

void MixedMRFModel::update_edge_indicator_discrete(int i, int j) {
    double k_curr = pairwise_effects_discrete_(i, j);
    double prop_sd = proposal_sd_pairwise_discrete_(i, j);

    int g_curr = gxx(i, j);
    int g_prop = 1 - g_curr;

    double k_prop;
    if(g_prop == 1) {
        k_prop = rnorm(rng_, k_curr, prop_sd);  // k_curr = 0 on a true add
    } else {
        k_prop = 0.0;
    }

    // --- Likelihood ratio ---
    double delta = k_prop - k_curr;
    double ll_curr = ll_marginal_cache_(i) + ll_marginal_cache_(j);

    pairwise_effects_discrete_(i, j) = k_prop;
    pairwise_effects_discrete_(j, i) = k_prop;
    refresh_marginal_interactions_entry(i, j);
    matvec_col_i_scratch_ = marginal_matvec_.col(i);
    matvec_col_j_scratch_ = marginal_matvec_.col(j);
    marginal_matvec_.col(i) += delta * discrete_observations_dbl_.col(j);
    marginal_matvec_.col(j) += delta * discrete_observations_dbl_.col(i);

    double marginal_prop_i = log_marginal_omrf_cached(i);
    double marginal_prop_j = log_marginal_omrf_cached(j);
    double ll_prop = marginal_prop_i + marginal_prop_j;

    // Restore; the accept branch re-applies
    pairwise_effects_discrete_(i, j) = k_curr;
    pairwise_effects_discrete_(j, i) = k_curr;
    refresh_marginal_interactions_entry(i, j);
    marginal_matvec_.col(i) = matvec_col_i_scratch_;
    marginal_matvec_.col(j) = matvec_col_j_scratch_;

    double ln_alpha = ll_prop - ll_curr;

    // Discrete slab prior: interaction_prior_
    // Must match the prior used in logp_and_gradient.
    if(g_prop == 1) {
        // Add: slab prior, subtract proposal density, inclusion prior
        ln_alpha += interaction_prior_->logp(k_prop);
        ln_alpha -= R::dnorm(k_prop, k_curr, prop_sd, true);
        ln_alpha += MY_LOG(inclusion_probability_(i, j))
                  - MY_LOG(1.0 - inclusion_probability_(i, j));
    } else {
        // Delete: subtract slab prior, add reverse proposal density, inclusion prior
        ln_alpha -= interaction_prior_->logp(k_curr);
        ln_alpha += R::dnorm(k_curr, k_prop, prop_sd, true);
        ln_alpha -= MY_LOG(inclusion_probability_(i, j))
                  - MY_LOG(1.0 - inclusion_probability_(i, j));
    }

    // Rao-Blackwellized inclusion draw (g_prop == 1 is a birth, gamma = 0).
    rb_edge_(i, j) = (g_prop == 1)
        ? MY_EXP(std::min(0.0, ln_alpha))
        : 1.0 - MY_EXP(std::min(0.0, ln_alpha));

    if(MY_LOG(runif(rng_)) < ln_alpha) {
        pairwise_effects_discrete_(i, j) = k_prop;
        pairwise_effects_discrete_(j, i) = k_prop;
        set_gxx(i, j, g_prop);
        constraint_dirty_ = true;
        refresh_marginal_interactions_entry(i, j);
        marginal_matvec_.col(i) += delta * discrete_observations_dbl_.col(j);
        marginal_matvec_.col(j) += delta * discrete_observations_dbl_.col(i);
        ll_marginal_cache_(i) = marginal_prop_i;
        ll_marginal_cache_(j) = marginal_prop_j;
    }
}


// =============================================================================
// update_edge_indicator_continuous
// =============================================================================
// Metropolis-Hastings add-delete move for a continuous-continuous edge (i, j).
// Uses Cholesky reparameterization (permute-free constants extraction).
// All proposals and constants live in precision space.
//   Add (G=0→1): propose ε ~ N(0, σ), precision_ij = C[3]*ε, constrain diagonal.
//   Delete (G=1→0): set precision_ij = 0, constrain diagonal.
// =============================================================================

void MixedMRFModel::update_edge_indicator_continuous(int i, int j) {
    get_precision_constants(i, j);

    int g_curr = gyy(i, j);
    int g_prop = 1 - g_curr;

    double theta_prop_ij, theta_prop_jj;

    if(g_prop == 1) {
        // Add: propose from N(0, σ) on reparameterized scale
        double epsilon = rnorm(rng_, 0.0, proposal_sd_pairwise_continuous_(i, j));
        theta_prop_ij = cont_constants_[3] * epsilon;
        theta_prop_jj = precision_constrained_diagonal(theta_prop_ij);
    } else {
        // Delete: set off-diagonal to 0 in precision space
        theta_prop_ij = 0.0;
        theta_prop_jj = precision_constrained_diagonal(0.0);
    }

    // Fill proposal in precision space
    precision_proposal_ = -2.0 * pairwise_effects_continuous_;
    precision_proposal_(i, j) = theta_prop_ij;
    precision_proposal_(j, i) = theta_prop_ij;
    precision_proposal_(j, j) = theta_prop_jj;

    // --- Likelihood ratio ---
    arma::mat cov_prop;
    double ggm_ratio = log_ggm_ratio_edge(i, j, cov_prop);
    double ln_alpha = ggm_ratio;

    // Determinant-tilt prior: see update_pairwise_effects_continuous_offdiag.
    if (determinant_tilt_yy_ != 0.0) {
        ln_alpha += determinant_tilt_yy_ * log_det_ratio_yy_edge(i, j);
    }

    // OMRF ratio with proposed continuous interactions. Use Σ' (cov_prop) for
    // the proposed-state marginals, not the cached Σ.
    ln_alpha += omrf_ratio_for_covariance_change(cov_prop);

    // --- Spike-and-slab terms ---
    // off-diagonal = -1/2 * precision_ij
    double cont_prop_ij = -0.5 * theta_prop_ij;
    double cont_curr_ij = pairwise_effects_continuous_(i, j);

    // Gamma(1,1) prior on changed diagonal K_jj (always present, not part of spike-and-slab)
    ln_alpha += diagonal_prior_->logp(0.5 * theta_prop_jj);
    ln_alpha -= diagonal_prior_->logp(0.5 * (-2.0 * pairwise_effects_continuous_(j, j)));

    // Slab in K_yy coords; proposal in K_ij coords. Jacobian |dK_yy/dK_ij| = 1/2.
    if(g_prop == 1) {
        // Add: slab prior on proposed off-diagonal
        ln_alpha += interaction_prior_->logp(cont_prop_ij) - MY_LOG(2.0);
        // Subtract proposal density (in K_ij coords)
        ln_alpha -= R::dnorm(theta_prop_ij / cont_constants_[3], 0.0,
                             proposal_sd_pairwise_continuous_(i, j), true)
                  - MY_LOG(cont_constants_[3]);
        // Inclusion prior: log(π / (1-π))
        ln_alpha += MY_LOG(inclusion_probability_(p_ + i, p_ + j))
                  - MY_LOG(1.0 - inclusion_probability_(p_ + i, p_ + j));
    } else {
        // Delete: subtract slab prior on current off-diagonal
        ln_alpha -= interaction_prior_->logp(cont_curr_ij) - MY_LOG(2.0);
        // Add reverse proposal density
        double theta_curr_ij = -2.0 * cont_curr_ij;
        ln_alpha += R::dnorm(theta_curr_ij / cont_constants_[3], 0.0,
                             proposal_sd_pairwise_continuous_(i, j), true)
                  - MY_LOG(cont_constants_[3]);
        // Inclusion prior: log((1-π) / π)
        ln_alpha -= MY_LOG(inclusion_probability_(p_ + i, p_ + j))
                  - MY_LOG(1.0 - inclusion_probability_(p_ + i, p_ + j));
    }

    // Hierarchical spec on the continuous block: the add ratio carries
    // +log J with J = Z(Gamma_yy-)/Z(Gamma_yy+) and the delete ratio -log J.
    // Mediating-block counts are read off the continuous subgraph; the
    // toggled edge's own state never enters.
    double log_j_cont = 0.0;
    if (zratio_engine_) {
        log_j_cont = zratio_engine_->log_zratio(continuous_subgraph(), i, j);
        ln_alpha += (g_prop == 1) ? log_j_cont : -log_j_cont;
    }

    if (zratio_gauge_.active) {
        zratio_gauge_record(zratio_gauge_, zratio_engine_.get(),
                            continuous_subgraph(), i, j, ln_alpha, log_j_cont,
                            g_prop == 1 ? 1 : -1);
    }

    // Rao-Blackwellized inclusion draw (g_prop == 1 is a birth, gamma = 0).
    rb_edge_(p_ + i, p_ + j) = (g_prop == 1)
        ? MY_EXP(std::min(0.0, ln_alpha))
        : 1.0 - MY_EXP(std::min(0.0, ln_alpha));

    if(MY_LOG(runif(rng_)) < ln_alpha) {
        // Pass old precision values to Cholesky update
        double old_theta_ij = -2.0 * pairwise_effects_continuous_(i, j);
        double old_theta_jj = -2.0 * pairwise_effects_continuous_(j, j);

        // Store: pairwise_effects_continuous_ = -1/2 * precision
        pairwise_effects_continuous_(i, j) = -0.5 * theta_prop_ij;
        pairwise_effects_continuous_(j, i) = -0.5 * theta_prop_ij;
        pairwise_effects_continuous_(j, j) = -0.5 * theta_prop_jj;

        set_gyy(i, j, g_prop);
        constraint_dirty_ = true;
        if (cholesky_update_after_precision_edge(old_theta_ij, old_theta_jj, i, j)) {
            adopt_kyy_proposal_caches(ggm_ratio, /*rank2=*/true);
        } else {
            recompute_am_caches();
        }
    }
}


// =============================================================================
// update_edge_indicator_cross
// =============================================================================
// Metropolis-Hastings add-delete move for a cross-type edge (i, j).
//   Add (G=0→1): propose k ~ N(0, σ).
//   Delete (G=1→0): set k = 0.
// =============================================================================

void MixedMRFModel::update_edge_indicator_cross(int i, int j) {
    double k_curr = pairwise_effects_cross_(i, j);
    double prop_sd = proposal_sd_pairwise_cross_(i, j);

    int g_curr = gxy(i, j);
    int g_prop = 1 - g_curr;

    double k_prop;
    if(g_prop == 1) {
        k_prop = rnorm(rng_, k_curr, prop_sd);  // k_curr = 0 on a true add
    } else {
        k_prop = 0.0;
    }

    // --- Likelihood ratio ---
    // Same rank-2 M update and rank-1 GGM quadratic-form delta as
    // update_pairwise_cross, with delta = k_prop - k_curr. Parameters and
    // conditional_mean_ are only touched on accept.
    double delta = k_prop - k_curr;
    double ll_curr = ll_ggm_cache_ + arma::accu(ll_marginal_cache_);

    arma::vec u = pairwise_effects_cross_ * covariance_continuous_.col(j);
    marginal_matvec_prop_ = marginal_matvec_
        + (2.0 * delta) * discrete_observations_dbl_.col(i) * u.t();
    marginal_matvec_prop_.col(i) +=
        (2.0 * delta) * (cross_matvec_ * covariance_continuous_.col(j))
        + (2.0 * delta * delta * covariance_continuous_(j, j))
              * discrete_observations_dbl_.col(i);
    mdiag_prop_ = marginal_interactions_.diag();
    mdiag_prop_(i) += 4.0 * delta * u(i)
                    + 2.0 * delta * delta * covariance_continuous_(j, j);
    cross_bias_prop_ = cross_bias_;
    cross_bias_prop_(i) += 2.0 * delta * main_effects_continuous_(j);

    arma::vec a = (2.0 * delta) * discrete_observations_dbl_.col(i);
    arma::vec kv = -2.0 * (pairwise_effects_continuous_ * covariance_continuous_.col(j));
    resid_scratch_ = continuous_observations_ - conditional_mean_;
    arma::vec g = resid_scratch_.t() * a;
    double quad_delta = -2.0 * arma::dot(g, kv)
        + arma::dot(a, a) * arma::dot(covariance_continuous_.col(j), kv);
    double ggm_prop = ll_ggm_cache_ - quad_delta / 2.0;

    for(size_t s = 0; s < p_; ++s)
        ll_marginal_prop_(s) = log_marginal_omrf_from(
            s, marginal_matvec_prop_, mdiag_prop_(s), cross_bias_prop_(s));

    double ll_prop = ggm_prop + arma::accu(ll_marginal_prop_);

    double ln_alpha = ll_prop - ll_curr;

    if(g_prop == 1) {
        // Add
        ln_alpha += interaction_prior_->logp(k_prop);
        ln_alpha -= R::dnorm(k_prop, k_curr, prop_sd, true);
        ln_alpha += MY_LOG(inclusion_probability_(i, p_ + j))
                  - MY_LOG(1.0 - inclusion_probability_(i, p_ + j));
    } else {
        // Delete
        ln_alpha -= interaction_prior_->logp(k_curr);
        ln_alpha += R::dnorm(k_curr, k_prop, prop_sd, true);
        ln_alpha -= MY_LOG(inclusion_probability_(i, p_ + j))
                  - MY_LOG(1.0 - inclusion_probability_(i, p_ + j));
    }

    // Rao-Blackwellized inclusion draw (g_prop == 1 is a birth, gamma = 0).
    rb_edge_(i, p_ + j) = (g_prop == 1)
        ? MY_EXP(std::min(0.0, ln_alpha))
        : 1.0 - MY_EXP(std::min(0.0, ln_alpha));

    if(MY_LOG(runif(rng_)) < ln_alpha) {
        pairwise_effects_cross_(i, j) = k_prop;
        set_gxy(i, j, g_prop);
        constraint_dirty_ = true;
        adopt_cross_proposal_caches(i, j, delta, u, ggm_prop);
    }
}

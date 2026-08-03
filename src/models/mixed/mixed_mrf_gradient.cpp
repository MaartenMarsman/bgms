// mixed_mrf_gradient.cpp — MixedMRFModel gradient engine (NUTS support).
//
// gradient, logp_and_gradient, and the full-space variant, plus the gradient-cache
// bookkeeping (ensure_gradient_cache, invalidate_gradient_cache) and the
// NUTS-vector -> working-temps unpacking. Pairs with the constrained-NUTS projection
// cache in mixed_mrf_model.cpp.
#include <RcppArmadillo.h>
#include "models/mixed/mixed_mrf_model.h"
#include "utils/variable_helpers.h"
#include "math/explog_macros.h"


// =============================================================================
// Gradient cache
// =============================================================================
// The gradient cache stores precomputed index mappings and observed-statistic
// contributions that do not change during a leapfrog trajectory.  It is
// invalidated whenever edge indicators change (same pattern as the OMRF).
// =============================================================================

void MixedMRFModel::ensure_gradient_cache() {
    if(gradient_cache_valid_) return;

    // --- Build index matrix for pairwise_effects_discrete_ upper-triangular entries ---
    // Maps (i, j) to a position in the flat gradient vector (offset from
    // the start of pairwise_discrete entries, which sits at num_main_).
    disc_index_cache_.set_size(p_, p_);
    disc_index_cache_.zeros();

    int num_active_disc = 0;
    for(size_t i = 0; i + 1 < p_; ++i) {
        for(size_t j = i + 1; j < p_; ++j) {
            if(edge_indicators_(i, j) == 1) {
                disc_index_cache_(i, j) = num_main_ + num_active_disc;
                disc_index_cache_(j, i) = disc_index_cache_(i, j);
                num_active_disc++;
            }
        }
    }

    // --- Build index matrix for pairwise_effects_cross_ entries ---
    // Maps (i, j) to a position in the flat gradient vector (offset from
    // the start of pairwise_cross entries, which sits at num_main_ + active_kxx + q).
    cross_index_cache_.set_size(p_, q_);
    cross_index_cache_.zeros();

    int cross_offset = num_main_ + num_active_disc + static_cast<int>(q_);
    int num_active_cross = 0;
    for(size_t i = 0; i < p_; ++i) {
        for(size_t j = 0; j < q_; ++j) {
            if(edge_indicators_(i, p_ + j) == 1) {
                cross_index_cache_(i, j) = cross_offset + num_active_cross;
                num_active_cross++;
            }
        }
    }

    // main_effects_continuous_ offset in gradient vector
    main_effects_continuous_grad_offset_ = num_main_ + num_active_disc;

    // --- Precompute observed statistics portion of the gradient ---
    size_t active_dim = num_main_ + num_active_disc + q_ + num_active_cross
                      + chol_constraint_structure_.active_dim;
    grad_obs_cache_.set_size(active_dim);
    grad_obs_cache_.zeros();

    // Cholesky block offset in gradient vector
    chol_grad_offset_ = num_main_ + num_active_disc + q_ + num_active_cross;

    // Observed statistics for discrete main effects
    int offset = 0;
    for(size_t s = 0; s < p_; ++s) {
        if(is_ordinal_variable_(s)) {
            int C_s = num_categories_(s);
            for(int c = 0; c < C_s; ++c) {
                grad_obs_cache_(offset + c) = counts_per_category_(c + 1, s);
            }
            offset += C_s;
        } else {
            grad_obs_cache_(offset)     = blume_capel_stats_(0, s);
            grad_obs_cache_(offset + 1) = blume_capel_stats_(1, s);
            offset += 2;
        }
    }

    // Observed statistics for pairwise_effects_discrete_ edges
    for(size_t i = 0; i + 1 < p_; ++i) {
        for(size_t j = i + 1; j < p_; ++j) {
            if(edge_indicators_(i, j) == 0) continue;
            int loc = disc_index_cache_(i, j);
            // Factor 4: K = σ, and the log-PL has edge (i,j) in two conditionals,
            // giving d/dK [4K·(x^Tx)] = 4·(x^Tx)
            grad_obs_cache_(loc) = 4.0 * arma::dot(
                discrete_observations_dbl_.col(i),
                discrete_observations_dbl_.col(j)
            );
        }
    }

    // No precomputed observed stats for means or cross effects — those depend on
    // continuous_observations_ combined with current parameters, so they
    // are computed fresh each logp_and_gradient call.

    // Cache transpose of discrete observations for vectorized pairwise gradient
    discrete_observations_dbl_t_ = discrete_observations_dbl_.t();

    gradient_cache_valid_ = true;
}


void MixedMRFModel::invalidate_gradient_cache() {
    gradient_cache_valid_ = false;
    theta_yy_valid_ = false;
}


// =============================================================================
// Unvectorize NUTS parameters into temporaries
// =============================================================================
// Unpacks a NUTS-dimension parameter vector into temporary matrices without
// mutating model state.  Used during leapfrog trajectory evaluation.
// =============================================================================

void MixedMRFModel::unvectorize_nuts_to_temps(
    const arma::vec& params,
    arma::mat& temp_main_discrete,
    arma::mat& temp_pairwise_discrete,
    arma::vec& temp_main_continuous,
    arma::mat& temp_pairwise_cross
) const {
    size_t idx = 0;

    // 1. main_effects_discrete_
    for(size_t s = 0; s < p_; ++s) {
        if(is_ordinal_variable_(s)) {
            for(int c = 0; c < num_categories_(s); ++c) {
                temp_main_discrete(s, c) = params(idx++);
            }
        } else {
            temp_main_discrete(s, 0) = params(idx++);
            temp_main_discrete(s, 1) = params(idx++);
        }
    }

    // 2. pairwise_effects_discrete_ upper-triangular (active only)
    for(size_t i = 0; i + 1 < p_; ++i) {
        for(size_t j = i + 1; j < p_; ++j) {
            if(edge_indicators_(i, j) == 1) {
                temp_pairwise_discrete(i, j) = params(idx++);
                temp_pairwise_discrete(j, i) = temp_pairwise_discrete(i, j);
            }
        }
    }

    // 3. main_effects_continuous_
    for(size_t j = 0; j < q_; ++j) {
        temp_main_continuous(j) = params(idx++);
    }

    // 4. pairwise_effects_cross_ row-major (active only)
    for(size_t i = 0; i < p_; ++i) {
        for(size_t j = 0; j < q_; ++j) {
            if(edge_indicators_(i, p_ + j) == 1) {
                temp_pairwise_cross(i, j) = params(idx++);
            }
        }
    }
}


// =============================================================================
// logp_and_gradient — marginal pseudo-likelihood
// =============================================================================
// Computes the log pseudo-posterior and its gradient with respect to the
// NUTS parameters (μ_x, A_xx, μ_y, A_xy, R) where R is the Cholesky
// factor of the continuous precision matrix Ω = R^T R.
//
// The pseudo-log-posterior is:
//   l(θ) = sum_s log p(x_s | x_{-s})       [OMRF conditionals — marginal]
//            + log p(y | x)                 [GGM conditional]
//            + log π(θ)                     [priors on all params]
//            + log |det J|                  [Cholesky Jacobian]
//
// The OMRF conditionals use the marginal effective-interaction matrix
// Θ = A_xx + 2 A_xy Σ_yy A_xy' (continuous block integrated out).
// =============================================================================

std::pair<double, arma::vec> MixedMRFModel::logp_and_gradient(
    const arma::vec& parameters)
{
    ensure_constraint_structure();
    ensure_gradient_cache();

    // --- Unvectorize into temporaries (blocks 1–4) ---
    arma::mat temp_main_discrete = main_effects_discrete_;
    arma::mat temp_pairwise_discrete = pairwise_effects_discrete_;
    arma::vec temp_main_continuous = main_effects_continuous_;
    arma::mat temp_pairwise_cross = pairwise_effects_cross_;
    unvectorize_nuts_to_temps(parameters, temp_main_discrete, temp_pairwise_discrete, temp_main_continuous, temp_pairwise_cross);

    // --- Unpack block 5: Kyy theta (f_q, psi_q) via the engine forward map,
    //     which enforces excluded-edge zeros through the null-space bases ---
    size_t chol_offset = static_cast<size_t>(chol_grad_offset_);
    size_t chol_dim = chol_constraint_structure_.active_dim;
    const ForwardMapResult& fm = yy_engine_.forward_map(
        arma::vec(parameters.subvec(chol_offset, chol_offset + chol_dim - 1)));
    const arma::mat& temp_cholesky = fm.Phi;

    // Guard against degenerate Cholesky (extreme theta pushed by leapfrog)
    double min_diag = temp_cholesky.diag().min();
    if(!std::isfinite(min_diag) || min_diag < 1e-15) {
        return {-std::numeric_limits<double>::infinity(),
                arma::vec(parameters.n_elem, arma::fill::zeros)};
    }

    const arma::mat& temp_precision = fm.K;
    arma::mat temp_inv_chol;
    bool solve_ok = arma::solve(temp_inv_chol, arma::trimatu(temp_cholesky),
                                arma::eye(q_, q_), arma::solve_opts::fast);
    if(!solve_ok) {
        return {-std::numeric_limits<double>::infinity(),
                arma::vec(parameters.n_elem, arma::fill::zeros)};
    }
    arma::mat temp_covariance = temp_inv_chol * temp_inv_chol.t();
    double temp_log_det = 2.0 * arma::accu(fm.psi);

    // --- Derived quantities ---
    // Conditional mean: M_i = μ_y' + 2 x_i' A_xy Σ_yy  (n x q)
    arma::mat temp_cond_mean = 2.0 * discrete_observations_dbl_ * temp_pairwise_cross * temp_covariance;
    temp_cond_mean.each_row() += temp_main_continuous.t();

    // Residual: D = Y - M  (n x q)
    arma::mat D = continuous_observations_ - temp_cond_mean;

    // Marginal PL effective discrete interaction matrix
    //   M = A_xx + 2 A_xy Σ_yy A_xy'
    // (see recompute_marginal_interactions in mixed_mrf_model.cpp).
    arma::mat temp_marginal;
    temp_marginal = temp_pairwise_discrete + 2.0 * temp_pairwise_cross * temp_covariance * temp_pairwise_cross.t();

    // Rest-score ingredients shared by every variable: one n x p GEMM instead
    // of p per-variable GEMVs, and the p cross-bias scalars in one GEMV.
    arma::mat X_marginal = discrete_observations_dbl_ * temp_marginal;
    arma::vec cross_bias = 2.0 * (temp_pairwise_cross * temp_main_continuous);

    // Start gradient from observed-statistics cache
    arma::vec grad = grad_obs_cache_;

    double logp = 0.0;

    // For marginal PL: precompute A_xy Σ_yy (used in cross-contributions)
    arma::mat cross_times_cov;  // p x q
    arma::mat Theta_bar;        // p x p marginal-PL coupling for precision gradient
    cross_times_cov = temp_pairwise_cross * temp_covariance;

    // =========================================================================
    // Part 1: OMRF conditionals
    // =========================================================================
    // The loop collects per-variable expected scores and residual scalars;
    // the pairwise, cross, and Θ̄ contributions are assembled after the loop
    // from batched GEMMs instead of per-variable GEMVs.

    arma::mat E_all(n_, p_, arma::fill::none);        // expected score per obs
    arma::vec diff_diag_all(p_, arma::fill::none);    // dot(x_s, x_s) - sum(E_sq)
    arma::vec sum_obs_minus_E_all(p_, arma::fill::none);

    int main_effects_discrete_offset = 0;
    for(size_t s = 0; s < p_; ++s) {
        int C_s = num_categories_(s);

        // --- Rest score for variable s ---
        arma::vec rest;
        // Marginal: Θ-based rest + A_xy μ_y bias
        double precision_ss = temp_marginal(s, s);
        rest = 2.0 * (X_marginal.col(s)
                    - discrete_observations_dbl_.col(s) * precision_ss)
             + cross_bias(s);

        // Numerator contribution to logp from discrete sufficient statistics
        // (already in grad_obs_cache_ as counts, but logp needs the actual
        // dot-products). Marginal self-interaction quadratic contribution,
        // dot(x_s, rest), and the main-effect sums.
        logp += precision_ss * arma::dot(
            discrete_observations_dbl_.col(s),
            discrete_observations_dbl_.col(s));
        logp += arma::dot(discrete_observations_dbl_.col(s), rest);
        if(is_ordinal_variable_(s)) {
            for(int c = 1; c <= C_s; ++c) {
                logp += static_cast<double>(counts_per_category_(c, s)) * temp_main_discrete(s, c - 1);
            }
        } else {
            logp += temp_main_discrete(s, 0) * static_cast<double>(blume_capel_stats_(0, s))
                  + temp_main_discrete(s, 1) * static_cast<double>(blume_capel_stats_(1, s));
        }

        if(is_ordinal_variable_(s)) {
            arma::vec main_param = temp_main_discrete.row(s).cols(0, C_s - 1).t();

            // Marginal PL: absorb marginal self-interaction into main_param
            double precision_ss = temp_marginal(s, s);
            for(int c = 0; c < C_s; ++c) {
                main_param(c) += static_cast<double>((c + 1) * (c + 1)) * precision_ss;
            }

            // bound = per-observation upper bound on log-scores for numerical
            // stability. Must cover max_c(main_param(c) + (c+1)*rest(i)).
            // The highest-category term main_param(C_s-1) + C_s*rest dominates
            // when rest > 0; category 0 (score = 0) dominates when rest << 0.
            arma::vec bound = main_param(C_s - 1) + static_cast<double>(C_s) * rest;
            bound = arma::max(bound, arma::zeros<arma::vec>(bound.n_elem));

            // Fill-in-place using persistent per-chain scratch.
            compute_logZ_and_probs_ordinal_into(
                main_param, rest, bound, C_s, logz_out_, logz_scratch_
            );

            // log pseudo-posterior contribution
            logp -= arma::accu(logz_out_.log_Z);

            // Main-effect gradient: ∂/∂main_effects_discrete_{s,c} = count_c - sum_i prob(c)
            for(int c = 0; c < C_s; ++c) {
                grad(main_effects_discrete_offset + c) -= arma::accu(logz_out_.probs.col(c + 1));
            }

            // Expected value E_s[c+1|rest] per observation
            arma::vec weights = arma::regspace<arma::vec>(1, C_s);
            E_all.col(s) = logz_out_.probs.cols(1, C_s) * weights;

            arma::vec weights_sq = arma::square(weights);
            arma::vec E_sq = logz_out_.probs.cols(1, C_s) * weights_sq;

            diff_diag_all(s) = arma::dot(
                discrete_observations_dbl_.col(s),
                discrete_observations_dbl_.col(s)) - arma::accu(E_sq);

            sum_obs_minus_E_all(s) = arma::accu(discrete_observations_dbl_.col(s))
                                   - arma::accu(E_all.col(s));

            main_effects_discrete_offset += C_s;
        } else {
            // --- Blume-Capel variable ---
            int ref = baseline_category_(s);
            double lin_eff = temp_main_discrete(s, 0);
            double quad_eff = temp_main_discrete(s, 1);

            // Marginal PL: absorb marginal self-interaction into quadratic effect
            double effective_quad = quad_eff;
            effective_quad += temp_marginal(s, s);

            arma::vec bc_bound;
            compute_logZ_and_probs_blume_capel_into(
                rest, lin_eff, effective_quad, ref, C_s, bc_bound,
                logz_out_, logz_scratch_
            );

            logp -= arma::accu(logz_out_.log_Z);

            arma::vec score = arma::regspace<arma::vec>(0, C_s) - static_cast<double>(ref);
            arma::vec sq_score = arma::square(score);

            // Main-effect gradient
            grad(main_effects_discrete_offset)     -= arma::accu(logz_out_.probs * score);
            grad(main_effects_discrete_offset + 1) -= arma::accu(logz_out_.probs * sq_score);

            // Expected score per person
            E_all.col(s) = logz_out_.probs * score;

            arma::vec E_sq = logz_out_.probs * sq_score;

            diff_diag_all(s) = arma::dot(
                discrete_observations_dbl_.col(s),
                discrete_observations_dbl_.col(s)) - arma::accu(E_sq);

            sum_obs_minus_E_all(s) = arma::accu(discrete_observations_dbl_.col(s))
                                   - arma::accu(E_all.col(s));

            main_effects_discrete_offset += 2;
        }
    }

    // --- Batched pairwise/cross gradient assembly ---
    // Two X^T * (...) GEMMs replace the 2p per-variable GEMVs:
    //   pw_grad_all.col(s)  = X^T E_s
    //   diff_pw_all.col(s)  = X^T (x_s - E_s), with the (s, s) entry zeroed
    arma::mat pw_grad_all = discrete_observations_dbl_t_ * E_all;
    arma::mat diff_pw_all = discrete_observations_dbl_t_
                          * (discrete_observations_dbl_ - E_all);
    diff_pw_all.diag().zeros();

    // Θ̄ coupling for the precision gradient: Θ̄_{s,t} = 2 diff_pw_s(t),
    // diagonal = dot(x_s, x_s) - sum(E_sq)
    Theta_bar = 2.0 * diff_pw_all.t();
    Theta_bar.diag() = diff_diag_all;

    // Pairwise discrete gradient: sum_i x_{i,t} * (x_{i,s}+1 - E_s)
    // Factor 2: chain rule d/dK = 2 × d/dσ; edge (i, j) collects the
    // variable-i and variable-j conditional contributions.
    for(size_t i = 0; i + 1 < p_; ++i) {
        for(size_t j = i + 1; j < p_; ++j) {
            if(edge_indicators_(i, j) == 0) continue;
            int loc = disc_index_cache_(i, j);
            grad(loc) -= 2.0 * pw_grad_all(j, i);
            grad(loc) -= 2.0 * pw_grad_all(i, j);
        }
    }

    // Pairwise_cross gradient from marginal OMRF (through Θ):
    // ∂marginal_{st}/∂pairwise_effects_cross_{a,j} has two terms:
    //   = 2 [Σyy pairwise_effects_cross_t']_j δ_{as} + 2 [pairwise_effects_cross_s Σyy]_j δ_{at}
    // Self-contribution (a=s): 4 diff_pw_s^T (A_xy Σyy) row s → diff_pw_all^T term
    // Cross-contribution (a=t): 4 diff_pw_s(t) (A_xy Σyy) row s → diff_pw_all term
    // Diagonal effective interaction: 4 diff_diag_s [A_xy Σyy] row s
    // Rest-score bias: 2 sum_i(x_{is} - E_s) μy^T
    arma::mat cross_grad_all =
        4.0 * ((diff_pw_all + diff_pw_all.t()) * cross_times_cov)
      + 4.0 * (cross_times_cov.each_col() % diff_diag_all)
      + 2.0 * (sum_obs_minus_E_all * temp_main_continuous.t());

    for(size_t s = 0; s < p_; ++s) {
        for(size_t j = 0; j < q_; ++j) {
            if(edge_indicators_(s, p_ + j) == 0) continue;
            grad(cross_index_cache_(s, j)) += cross_grad_all(s, j);
        }
    }

    // Continuous mean gradient from marginal OMRF:
    // ∂l_s/∂main_effects_continuous_j = 2 pairwise_effects_cross_{sj} * sum_i (x_{is} - E_s)
    arma::vec mean_grad_omrf = 2.0 * (temp_pairwise_cross.t() * sum_obs_minus_E_all);
    for(size_t j = 0; j < q_; ++j) {
        grad(main_effects_continuous_grad_offset_ + j) += mean_grad_omrf(j);
    }

    // =========================================================================
    // Part 2: GGM conditional log-likelihood and gradients
    // =========================================================================
    // log p(y | x) = n/2 (log|Ω| - q log(2π)) - ½ trace(Ω D'D)
    // where Ω = R'R (precision), D = Y - M

    double quad_sum = arma::accu((D * temp_precision) % D);
    logp += static_cast<double>(n_) / 2.0 *
            (-static_cast<double>(q_) * MY_LOG(2.0 * arma::datum::pi)
             + temp_log_det)
          - quad_sum / 2.0;

    // ∂/∂μ_y: Ω * sum_over_rows(D)
    arma::vec D_colsums = arma::sum(D, 0).t();  // q-vector
    arma::vec grad_main_effects_continuous_ggm = temp_precision * D_colsums;

    for(size_t j = 0; j < q_; ++j) {
        grad(main_effects_continuous_grad_offset_ + j) += grad_main_effects_continuous_ggm(j);
    }

    // ∂/∂A_xy: The GGM conditional depends on A_xy through M.
    // ∂M/∂pairwise_effects_cross_{s,j} = 2 x_s [Σ_yy]_{j,:}
    // ∂logp_ggm/∂A_xy = 2 X' D  (shortcut: Θ Σ_yy = I eliminates Θ)
    //
    // Correctly: ∂(−½ trace(Θ D'D))/∂pairwise_effects_cross_{s,j}
    //   = trace(Θ D' ∂M/∂pairwise_effects_cross_{s,j})
    //   = trace(Θ D' · 2 x_s [Σ_yy]_{j,:})
    //   = 2 [x_s' D Θ Σ_yy]_j
    //   = 2 [x_s' D]_j    (since Θ Σ_yy = I)
    arma::mat grad_pairwise_effects_cross_ggm = 2.0 * discrete_observations_dbl_t_ * D;  // p x q

    for(size_t i = 0; i < p_; ++i) {
        for(size_t j = 0; j < q_; ++j) {
            if(edge_indicators_(i, p_ + j) == 0) continue;
            int loc = cross_index_cache_(i, j);
            grad(loc) += grad_pairwise_effects_cross_ggm(i, j);
        }
    }

    // =========================================================================
    // Part 3: Prior log-densities and gradient contributions
    // =========================================================================

    // --- main_effects_discrete_ priors ---
    main_effects_discrete_offset = 0;
    for(size_t s = 0; s < p_; ++s) {
        if(is_ordinal_variable_(s)) {
            int C_s = num_categories_(s);
            for(int c = 0; c < C_s; ++c) {
                double val = temp_main_discrete(s, c);
                logp += threshold_prior_->logp(val);
                grad(main_effects_discrete_offset + c) += threshold_prior_->grad(val);
            }
            main_effects_discrete_offset += C_s;
        } else {
            for(int k = 0; k < 2; ++k) {
                double val = temp_main_discrete(s, k);
                logp += threshold_prior_->logp(val);
                grad(main_effects_discrete_offset + k) += threshold_prior_->grad(val);
            }
            main_effects_discrete_offset += 2;
        }
    }

    // --- pairwise_effects_discrete_ priors ---
    for(size_t i = 0; i + 1 < p_; ++i) {
        for(size_t j = i + 1; j < p_; ++j) {
            if(edge_indicators_(i, j) == 0) continue;
            int loc = disc_index_cache_(i, j);
            double val = temp_pairwise_discrete(i, j);
            logp += interaction_prior_->logp(val);
            grad(loc) += interaction_prior_->grad(val);
        }
    }

    // --- main_effects_continuous_ priors: Normal(0, 1) ---
    for(size_t j = 0; j < q_; ++j) {
        double val = temp_main_continuous(j);
        logp += means_prior_->logp(val);
        grad(main_effects_continuous_grad_offset_ + j) += means_prior_->grad(val);
    }

    // --- pairwise_effects_cross_ priors ---
    for(size_t i = 0; i < p_; ++i) {
        for(size_t j = 0; j < q_; ++j) {
            if(edge_indicators_(i, p_ + j) == 0) continue;
            int loc = cross_index_cache_(i, j);
            double val = temp_pairwise_cross(i, j);
            logp += interaction_prior_->logp(val);
            grad(loc) += interaction_prior_->grad(val);
        }
    }

    // =========================================================================
    // Part 4: Precision gradient via Cholesky parameterization
    // =========================================================================
    // Compute Ω̄ = ∂ℓ/∂Ω, then map to R̄ = ∂ℓ/∂R, then to position gradient.
    //
    // Ω̄ = (n/2) Σ − ½ D^T D − 2 Σ A_xy^T X^T D + priors on Ω
    //     + [marginal PL coupling through Θ]

    // --- Phase 1: GGM conditional contribution ---
    arma::mat Omega_bar = 0.5 * static_cast<double>(n_) * temp_covariance
                        - 0.5 * D.t() * D;

    // --- Phase 2: Conditional-mean coupling ---
    // M_i = μ_y + 2 Σ A_xy^T x_i depends on Σ = Ω^{-1}
    // ∂ℓ/∂Σ_{ab} from GGM conditional = 2 [A_xy^T X^T D Ω]_{ab}
    // Mapping: ∂ℓ/∂Ω += −Σ (∂ℓ/∂Σ) Σ = −2 Σ A_xy^T X^T D
    Omega_bar -= 2.0 * temp_covariance * temp_pairwise_cross.t()
               * discrete_observations_dbl_t_ * D;

    // --- Phase 2b: Marginal PL coupling through M ---
    // M = A_xx + 2 A_xy Σ A_xy^T depends on Σ
    // ∂M_{st}/∂Σ_{ab} = 2 A_xy_{s,a} A_xy_{t,b}
    //   →  ∂ℓ/∂Σ = 2 A_xy^T Θ̄ A_xy,  Θ̄_{s,t} = ∂l_s/∂M_{s,t}
    // ∂ℓ/∂Ω += −Σ (∂ℓ/∂Σ) Σ = −2 Σ A_xy^T Θ̄ A_xy Σ
    Omega_bar -= 2.0 * temp_covariance * temp_pairwise_cross.t()
               * Theta_bar * temp_pairwise_cross * temp_covariance;

    // --- Phase 3: Priors on precision entries ---
    // Prior on the partial-association diagonal: -K_yy_{jj} = Theta_{jj}/2.
    //   logp uses 0.5 * Theta_jj.
    //   d/dTheta_jj log p(Theta_jj/2) = 0.5 * grad(Theta_jj/2).
    for(size_t j = 0; j < q_; ++j) {
        double half_kjj = 0.5 * temp_precision(j, j);
        logp += diagonal_prior_->logp(half_kjj);
        Omega_bar(j, j) += 0.5 * diagonal_prior_->grad(half_kjj);
    }
    // Interaction prior on off-diagonal Kyy_{ij} = -Ω_{ij}/2 (upper triangle only).
    // Gated on edge_indicators_: inactive edges have K_yy_{ij} = 0 (point mass)
    // and contribute no slab density, matching the GGM convention at
    // ggm_gradient.cpp where the slab is summed over included_indices only.
    // Only add to Omega_bar(i,j), not (j,i): the symmetrization
    // Ω̄ + Ω̄ᵀ in Phase 4 handles the lower triangle automatically.
    // The prior is on Kyy_{ij}, so we evaluate at -Ω_{ij}/2 and apply
    // chain rule: ∂logπ/∂Ω_{ij} = ∂logπ/∂Kyy_{ij} · (-1/2).
    for(size_t i = 0; i + 1 < q_; ++i) {
        for(size_t j = i + 1; j < q_; ++j) {
            if(edge_indicators_(p_ + i, p_ + j) == 0) continue;
            double kyy_val = -0.5 * temp_precision(i, j);
            logp += interaction_prior_->logp(kyy_val);
            Omega_bar(i, j) += -0.5 * interaction_prior_->grad(kyy_val);
        }
    }

    // --- Phase 4: Map Ω̄ → R̄ → position gradient ---
    // R̄ = R (Ω̄ + Ω̄^T)
    arma::mat Omega_bar_sym = Omega_bar + Omega_bar.t();
    arma::mat R_bar = temp_cholesky * Omega_bar_sym;

    // Parameterization Jacobian for the (f_q, psi_q) null-space coordinates:
    //   log|det J| = q*log(2) + 2*sum(psi) + sum_i (q-1-i)*psi_i
    //              - sum_q sum_j log R_diag_{q,j}
    // (the R_diag term replaces the RATTLE-era Pfaffian; the two agree at
    // identity mass since A_q^T = Q_q R_q gives A_q A_q^T = R_q^T R_q).
    logp += fm.log_det_jacobian;

    // Determinant tilt on the Kyy block: adds delta * log|Kyy| = 2*delta * sum(psi)
    // to the log-prior. Pushes the continuous-block precision matrix away from
    // the PD-cone boundary. delta = 0 recovers the untilted target.
    logp += determinant_tilt_yy_ * temp_log_det;

    if(!std::isfinite(logp)) {
        return {logp, arma::vec(grad.n_elem, arma::fill::zeros)};
    }

    // Extract the (f_q, psi_q) gradient from the Phi-space adjoint via the
    // engine's reverse-Givens pass. The determinant tilt contributes
    // +2*delta to every psi gradient; the Jacobian terms are handled inside.
    yy_engine_.theta_gradient_from_phi_bar(
        parameters, static_cast<size_t>(chol_grad_offset_), fm, R_bar,
        2.0 * determinant_tilt_yy_, grad);

    return {logp, grad};
}

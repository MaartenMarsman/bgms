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
    for(size_t i = 0; i < p_ - 1; ++i) {
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
                      + num_cholesky_;
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
    for(size_t i = 0; i < p_ - 1; ++i) {
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
    for(size_t i = 0; i < p_ - 1; ++i) {
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
    ensure_gradient_cache();
    ensure_constraint_structure();

    // --- Unvectorize into temporaries (blocks 1–4) ---
    arma::mat temp_main_discrete = main_effects_discrete_;
    arma::mat temp_pairwise_discrete = pairwise_effects_discrete_;
    arma::vec temp_main_continuous = main_effects_continuous_;
    arma::mat temp_pairwise_cross = pairwise_effects_cross_;
    unvectorize_nuts_to_temps(parameters, temp_main_discrete, temp_pairwise_discrete, temp_main_continuous, temp_pairwise_cross);

    // --- Unpack block 5: Cholesky of precision ---
    arma::mat temp_cholesky(q_, q_, arma::fill::zeros);
    size_t chol_idx = static_cast<size_t>(chol_grad_offset_);
    for(size_t j = 0; j < q_; ++j) {
        for(size_t i = 0; i < j; ++i) {
            temp_cholesky(i, j) = parameters(chol_idx++);
        }
        temp_cholesky(j, j) = std::exp(parameters(chol_idx++));
    }

    // Guard against degenerate Cholesky (extreme theta pushed by leapfrog)
    double min_diag = temp_cholesky.diag().min();
    if(!std::isfinite(min_diag) || min_diag < 1e-15) {
        return {-std::numeric_limits<double>::infinity(),
                arma::vec(parameters.n_elem, arma::fill::zeros)};
    }

    arma::mat temp_precision = temp_cholesky.t() * temp_cholesky;
    arma::mat temp_inv_chol;
    bool solve_ok = arma::solve(temp_inv_chol, arma::trimatu(temp_cholesky),
                                arma::eye(q_, q_), arma::solve_opts::fast);
    if(!solve_ok) {
        return {-std::numeric_limits<double>::infinity(),
                arma::vec(parameters.n_elem, arma::fill::zeros)};
    }
    arma::mat temp_covariance = temp_inv_chol * temp_inv_chol.t();
    double temp_log_det = 2.0 * arma::sum(arma::log(temp_cholesky.diag()));

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
    Theta_bar = arma::zeros<arma::mat>(p_, p_);

    // =========================================================================
    // Part 1: OMRF conditionals
    // =========================================================================

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
            arma::vec E = logz_out_.probs.cols(1, C_s) * weights;

            // Pairwise discrete gradient: sum_i x_{i,t} * (x_{i,s}+1 - E_s)
            // (uses pre-transposed discrete observations for BLAS efficiency)
            // Factor 2: chain rule d/dK = 2 × d/dσ
            arma::vec pw_grad = discrete_observations_dbl_t_ * E;
            for(size_t t = 0; t < p_; ++t) {
                if(edge_indicators_(s, t) == 0 || s == t) continue;
                int loc = (s < t) ? disc_index_cache_(s, t) : disc_index_cache_(t, s);
                grad(loc) -= 2.0 * pw_grad(t);
            }

            // Additional pairwise_discrete gradient from Θ_ss in denominator:
            // ∂/∂pairwise_effects_discrete_{st} through Θ_ss: zero (∂Θ_ss/∂pairwise_effects_discrete_st = δ_{st})
            // So pairwise_discrete gradient from Θ rest scores is already handled above.

            // Pairwise_cross gradient from marginal OMRF (through Θ):
            // ∂marginal_{st}/∂pairwise_effects_cross_{a,j} has two terms:
            //   = 2 [Σyy pairwise_effects_cross_t']_j δ_{as} + 2 [pairwise_effects_cross_s Σyy]_j δ_{at}
            // Self-contribution (a=s): from rest_s → pairwise_effects_cross_s
            // Cross-contribution (a=t): from rest_s → pairwise_effects_cross_t for each t≠s

            arma::vec weights_sq = arma::square(weights);
            arma::vec E_sq = logz_out_.probs.cols(1, C_s) * weights_sq;

            arma::vec diff_pw = discrete_observations_dbl_t_ *
                (discrete_observations_dbl_.col(s) - E);
            diff_pw(s) = 0.0;

            double diff_diag = arma::dot(
                discrete_observations_dbl_.col(s),
                discrete_observations_dbl_.col(s)) - arma::accu(E_sq);

            double sum_obs_minus_E = arma::accu(discrete_observations_dbl_.col(s)) - arma::accu(E);

            // Accumulate Θ̄ for precision gradient coupling
            for(size_t t = 0; t < p_; ++t) {
                if(t != s) Theta_bar(s, t) += 2.0 * diff_pw(t);
            }
            Theta_bar(s, s) += diff_diag;

            // Self-contribution: a = s
            // Off-diagonal effective interaction: ∂Θ_{st}/∂pairwise_effects_cross_{s,j} = 2 [Σyy pairwise_effects_cross_t']_j
            // Diagonal effective interaction: ∂Θ_{ss}/∂pairwise_effects_cross_{s,j} = 4 [Σyy pairwise_effects_cross_s']_j
            // Rest-score bias: ∂(2 pairwise_effects_cross_s μy)/∂pairwise_effects_cross_{s,j} = 2 μy_j
            arma::rowvec cross_self = 4.0 * (diff_pw.t() * temp_pairwise_cross) * temp_covariance
                                  + 4.0 * diff_diag * cross_times_cov.row(s)
                                  + 2.0 * sum_obs_minus_E * temp_main_continuous.t();

            for(size_t j = 0; j < q_; ++j) {
                if(edge_indicators_(s, p_ + j) == 0) continue;
                int loc = cross_index_cache_(s, j);
                grad(loc) += cross_self(j);
            }

            // Cross-contribution: a = t, for each t ≠ s
            // ∂l_s/∂pairwise_effects_cross_{t,:} = diff_pw(t) * 2 * pairwise_effects_cross_s * Σyy
            arma::rowvec V_s = 4.0 * cross_times_cov.row(s);
            for(size_t t = 0; t < p_; ++t) {
                if(t == s || std::abs(diff_pw(t)) < 1e-300) continue;
                for(size_t j = 0; j < q_; ++j) {
                    if(edge_indicators_(t, p_ + j) == 0) continue;
                    int loc = cross_index_cache_(t, j);
                    grad(loc) += diff_pw(t) * V_s(j);
                }
            }

            // Continuous mean gradient from marginal OMRF:
            // ∂l_s/∂main_effects_continuous_j = 2 pairwise_effects_cross_{sj} * sum_i (x_{is} - E_s)
            for(size_t j = 0; j < q_; ++j) {
                grad(main_effects_continuous_grad_offset_ + j) += 2.0 * temp_pairwise_cross(s, j) * sum_obs_minus_E;
            }

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
            arma::vec E = logz_out_.probs * score;

            // Pairwise discrete gradient
            // Factor 2: chain rule d/dK = 2 × d/dσ
            arma::vec pw_grad = discrete_observations_dbl_t_ * E;
            for(size_t t = 0; t < p_; ++t) {
                if(edge_indicators_(s, t) == 0 || s == t) continue;
                int loc = (s < t) ? disc_index_cache_(s, t) : disc_index_cache_(t, s);
                grad(loc) -= 2.0 * pw_grad(t);
            }

            // Pairwise_cross gradient from marginal OMRF (same structure as ordinal)
            arma::vec E_sq = logz_out_.probs * sq_score;

            arma::vec diff_pw = discrete_observations_dbl_t_ *
                (discrete_observations_dbl_.col(s) - E);
            diff_pw(s) = 0.0;

            double diff_diag = arma::dot(
                discrete_observations_dbl_.col(s),
                discrete_observations_dbl_.col(s)) - arma::accu(E_sq);

            double sum_obs_minus_E = arma::accu(discrete_observations_dbl_.col(s)) - arma::accu(E);

            // Accumulate Θ̄ for precision gradient coupling
            for(size_t t = 0; t < p_; ++t) {
                if(t != s) Theta_bar(s, t) += 2.0 * diff_pw(t);
            }
            Theta_bar(s, s) += diff_diag;

            // Self-contribution: a = s
            arma::rowvec cross_self = 4.0 * (diff_pw.t() * temp_pairwise_cross) * temp_covariance
                                  + 4.0 * diff_diag * cross_times_cov.row(s)
                                  + 2.0 * sum_obs_minus_E * temp_main_continuous.t();

            for(size_t j = 0; j < q_; ++j) {
                if(edge_indicators_(s, p_ + j) == 0) continue;
                int loc = cross_index_cache_(s, j);
                grad(loc) += cross_self(j);
            }

            // Cross-contribution: a = t, for each t ≠ s
            arma::rowvec V_s = 4.0 * cross_times_cov.row(s);
            for(size_t t = 0; t < p_; ++t) {
                if(t == s || std::abs(diff_pw(t)) < 1e-300) continue;
                for(size_t j = 0; j < q_; ++j) {
                    if(edge_indicators_(t, p_ + j) == 0) continue;
                    int loc = cross_index_cache_(t, j);
                    grad(loc) += diff_pw(t) * V_s(j);
                }
            }

            // Continuous mean gradient from marginal OMRF
            for(size_t j = 0; j < q_; ++j) {
                grad(main_effects_continuous_grad_offset_ + j) += 2.0 * temp_pairwise_cross(s, j) * sum_obs_minus_E;
            }

            main_effects_discrete_offset += 2;
        }
    }

    // Add numerator contribution to logp from discrete sufficient statistics
    // (already in grad_obs_cache_ as counts, but logp needs the actual dot-products)
    main_effects_discrete_offset = 0;
    for(size_t s = 0; s < p_; ++s) {
        int C_s = num_categories_(s);
        arma::vec rest;
        double precision_ss = temp_marginal(s, s);
        rest = 2.0 * (X_marginal.col(s)
                    - discrete_observations_dbl_.col(s) * precision_ss)
             + cross_bias(s);
        // Marginal self-interaction quadratic contribution
        logp += precision_ss * arma::dot(
            discrete_observations_dbl_.col(s),
            discrete_observations_dbl_.col(s));
        // Numerator: dot(x_s, rest) + main-effect sums
        logp += arma::dot(discrete_observations_dbl_.col(s), rest);

        if(is_ordinal_variable_(s)) {
            for(int c = 1; c <= C_s; ++c) {
                logp += static_cast<double>(counts_per_category_(c, s)) * temp_main_discrete(s, c - 1);
            }
        } else {
            logp += temp_main_discrete(s, 0) * static_cast<double>(blume_capel_stats_(0, s))
                  + temp_main_discrete(s, 1) * static_cast<double>(blume_capel_stats_(1, s));
        }
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
    for(size_t i = 0; i < p_ - 1; ++i) {
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
    for(size_t i = 0; i < q_ - 1; ++i) {
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

    // Cholesky-to-K Jacobian (graph-agnostic) + per-column Pfaffian correction.
    // Mirrors GGMGradientEngine::logp_and_gradient_full:
    //   ldj      = q*log(2) + Σ_j (q+1-j) ψ_j
    //   pfaffian = 0.5 * Σ_qq log det(A_qq A_qq^T)        (identity mass here)
    //   logp    += ldj - pfaffian
    // Theta-space integration uses identity mass on the manifold; for the
    // full-space (RATTLE) path we plug through the integrator's inverse-mass
    // diagonal. For the full Kyy graph, all A_qq are empty, so the Pfaffian
    // collapses to 0 and the formula reduces to the original q-j+1 weight.
    logp += static_cast<double>(q_) * std::log(2.0);
    for(size_t j = 0; j < q_; ++j) {
        logp += static_cast<double>(q_ + 1 - j) * std::log(temp_cholesky(j, j));
    }

    // Determinant tilt on the Kyy block: adds delta * log|Kyy| = 2*delta * sum(psi)
    // to the log-prior. Pushes the continuous-block precision matrix away from
    // the PD-cone boundary. delta = 0 recovers the untilted target.
    logp += determinant_tilt_yy_ * temp_log_det;

    const auto& cs = chol_constraint_structure_;
    arma::mat Aq_buf;
    std::vector<arma::mat> G_chol(q_);
    std::vector<arma::mat> Aq_cache(q_);
    double pfaffian = 0.0;
    for(size_t col = 1; col < q_; ++col) {
        const auto& cc = cs.columns[col];
        if(cc.m_q == 0) continue;

        GGMGradientEngine::build_Aq(temp_cholesky, cc, col, Aq_buf);
        Aq_cache[col] = Aq_buf;

        // Identity mass (theta-space): G_q = A_q A_q^T.
        arma::mat G_q = Aq_buf * Aq_buf.t();

        arma::mat L_q;
        bool chol_ok = arma::chol(L_q, G_q, "lower");
        if(!chol_ok) {
            double ridge = 1e-12 * (arma::trace(G_q) /
                                    static_cast<double>(cc.m_q) + 1.0);
            chol_ok = arma::chol(L_q, G_q + ridge * arma::eye(cc.m_q, cc.m_q),
                                 "lower");
            if(!chol_ok) {
                return {-std::numeric_limits<double>::infinity(),
                        arma::vec(grad.n_elem, arma::fill::zeros)};
            }
        }
        G_chol[col] = L_q;
        pfaffian += arma::accu(arma::log(arma::diagvec(L_q)));
    }
    logp -= pfaffian;

    // Pfaffian adjoint: d/dA_q [-0.5 log det(A_q A_q^T)] = -G_q^{-1} A_q.
    // Each A_q(r, l) = R(l, i_r) for l <= i_r, so the adjoint flows back to
    // column i_r of R_bar at rows l = 0..i_r.
    for(size_t col = 1; col < q_; ++col) {
        const auto& cc = cs.columns[col];
        if(cc.m_q == 0) continue;

        const arma::mat& L_q = G_chol[col];
        const arma::mat& Aq = Aq_cache[col];

        arma::mat Z = arma::solve(arma::trimatl(L_q), Aq,
                                  arma::solve_opts::fast);
        Z = arma::solve(arma::trimatu(L_q.t()), Z,
                        arma::solve_opts::fast);

        for(size_t r = 0; r < cc.m_q; ++r) {
            size_t i_r = cc.excluded_indices[r];
            for(size_t l = 0; l <= i_r; ++l) {
                R_bar(l, i_r) -= Z(r, l);
            }
        }
    }

    // Extract position gradient from R_bar with the unified weight (q+1-j)
    // on the diagonal-psi entries. The determinant tilt adds +2*delta to the
    // psi-bar entries (d/dpsi [delta * 2 * sum(psi)] = 2 * delta).
    size_t gidx = static_cast<size_t>(chol_grad_offset_);
    for(size_t j = 0; j < q_; ++j) {
        double w_j = static_cast<double>(q_ + 1 - j);
        for(size_t i = 0; i < j; ++i) {
            grad(gidx++) = R_bar(i, j);
        }
        grad(gidx++) = R_bar(j, j) * temp_cholesky(j, j) + w_j
                       + 2.0 * determinant_tilt_yy_;
    }

    return {logp, grad};
}


// =============================================================================
// logp_and_gradient_full — full-space version for RATTLE
// =============================================================================
// Same computation as logp_and_gradient() but with fixed full-dimension
// indexing: all Kxx, Kxy slots present (zeros for excluded edges).
// No edge gating — gradient is computed for every parameter.
// The Cholesky Jacobian is included (parameterization, not constraint).
// =============================================================================

std::pair<double, arma::vec> MixedMRFModel::logp_and_gradient_full(
    const arma::vec& x)
{
    ensure_constraint_structure();
    const size_t full_dim = full_parameter_dimension();

    // --- Unpack all 5 blocks from full-space vector ---
    size_t idx = 0;

    // Block 1: main_effects_discrete_
    arma::mat temp_main_discrete(p_, max_cats_, arma::fill::zeros);
    for(size_t s = 0; s < p_; ++s) {
        if(is_ordinal_variable_(s)) {
            for(int c = 0; c < num_categories_(s); ++c) {
                temp_main_discrete(s, c) = x(idx++);
            }
        } else {
            temp_main_discrete(s, 0) = x(idx++);
            temp_main_discrete(s, 1) = x(idx++);
        }
    }

    // Block 2: ALL pairwise_effects_discrete_ upper-triangular
    arma::mat temp_pairwise_discrete(p_, p_, arma::fill::zeros);
    for(size_t i = 0; i < p_ - 1; ++i) {
        for(size_t j = i + 1; j < p_; ++j) {
            temp_pairwise_discrete(i, j) = x(idx);
            temp_pairwise_discrete(j, i) = x(idx);
            idx++;
        }
    }

    // Block 3: main_effects_continuous_
    arma::vec temp_main_continuous(q_);
    for(size_t j = 0; j < q_; ++j) {
        temp_main_continuous(j) = x(idx++);
    }

    // Block 4: ALL pairwise_effects_cross_ row-major
    arma::mat temp_pairwise_cross(p_, q_, arma::fill::zeros);
    for(size_t i = 0; i < p_; ++i) {
        for(size_t j = 0; j < q_; ++j) {
            temp_pairwise_cross(i, j) = x(idx++);
        }
    }

    // Block 5: Cholesky of precision
    const size_t chol_offset = idx;
    arma::mat temp_cholesky(q_, q_, arma::fill::zeros);
    for(size_t j = 0; j < q_; ++j) {
        for(size_t i = 0; i < j; ++i) {
            temp_cholesky(i, j) = x(idx++);
        }
        temp_cholesky(j, j) = std::exp(x(idx++));
    }

    // Guard against degenerate Cholesky
    double min_diag = temp_cholesky.diag().min();
    if(!std::isfinite(min_diag) || min_diag < 1e-15) {
        return {-std::numeric_limits<double>::infinity(),
                arma::vec(full_dim, arma::fill::zeros)};
    }

    arma::mat temp_precision = temp_cholesky.t() * temp_cholesky;
    arma::mat temp_inv_chol;
    bool solve_ok = arma::solve(temp_inv_chol, arma::trimatu(temp_cholesky),
                                arma::eye(q_, q_), arma::solve_opts::fast);
    if(!solve_ok) {
        return {-std::numeric_limits<double>::infinity(),
                arma::vec(full_dim, arma::fill::zeros)};
    }
    arma::mat temp_covariance = temp_inv_chol * temp_inv_chol.t();
    double temp_log_det = 2.0 * arma::sum(arma::log(temp_cholesky.diag()));

    // --- Derived quantities ---
    arma::mat temp_cond_mean = 2.0 * discrete_observations_dbl_ * temp_pairwise_cross * temp_covariance;
    temp_cond_mean.each_row() += temp_main_continuous.t();
    arma::mat D = continuous_observations_ - temp_cond_mean;

    arma::mat temp_marginal;
    temp_marginal = temp_pairwise_discrete + 2.0 * temp_pairwise_cross * temp_covariance * temp_pairwise_cross.t();

    // Rest-score ingredients shared by every variable: one n x p GEMM instead
    // of p per-variable GEMVs, and the p cross-bias scalars in one GEMV.
    arma::mat X_marginal = discrete_observations_dbl_ * temp_marginal;
    arma::vec cross_bias = 2.0 * (temp_pairwise_cross * temp_main_continuous);

    // Initialize gradient (full dimension)
    arma::vec grad(full_dim, arma::fill::zeros);

    double logp = 0.0;

    // Full-space index offsets (fixed layout)
    const size_t kxx_offset = num_main_;
    const size_t mean_offset = num_main_ + num_pairwise_xx_;
    const size_t kxy_offset = num_main_ + num_pairwise_xx_ + q_;

    // For marginal PL: precompute helpers
    arma::mat cross_times_cov;
    arma::mat Theta_bar;
    cross_times_cov = temp_pairwise_cross * temp_covariance;
    Theta_bar = arma::zeros<arma::mat>(p_, p_);

    // =========================================================================
    // Part 1: OMRF conditionals (same as active-space but full indexing)
    // =========================================================================

    // Helper: flat index for Kxx(i,j) in full vector (i < j)
    auto kxx_idx = [&](size_t i, size_t j) -> size_t {
        // Row-major upper triangle: offset = sum_{r=0}^{i-1} (p-1-r) + (j-i-1)
        return kxx_offset + i * (2 * p_ - 3 - i) / 2 + (j - 1);
    };

    // Helper: flat index for Kxy(i,j) in full vector
    auto kxy_idx = [&](size_t i, size_t j) -> size_t {
        return kxy_offset + i * q_ + j;
    };

    // Precompute observed statistics for Kxx gradient (same as grad_obs_cache_)
    // Main effect observed stats
    int main_offset = 0;
    for(size_t s = 0; s < p_; ++s) {
        if(is_ordinal_variable_(s)) {
            int C_s = num_categories_(s);
            for(int c = 0; c < C_s; ++c) {
                grad(main_offset + c) = counts_per_category_(c + 1, s);
            }
            main_offset += C_s;
        } else {
            grad(main_offset)     = blume_capel_stats_(0, s);
            grad(main_offset + 1) = blume_capel_stats_(1, s);
            main_offset += 2;
        }
    }

    // Kxx observed stats
    for(size_t i = 0; i < p_ - 1; ++i) {
        for(size_t j = i + 1; j < p_; ++j) {
            grad(kxx_idx(i, j)) = 4.0 * arma::dot(
                discrete_observations_dbl_.col(i),
                discrete_observations_dbl_.col(j));
        }
    }

    // OMRF conditionals loop
    int main_effects_discrete_offset = 0;
    for(size_t s = 0; s < p_; ++s) {
        int C_s = num_categories_(s);

        // Rest score
        arma::vec rest;
        double precision_ss = temp_marginal(s, s);
        rest = 2.0 * (X_marginal.col(s)
                    - discrete_observations_dbl_.col(s) * precision_ss)
             + cross_bias(s);

        if(is_ordinal_variable_(s)) {
            arma::vec main_param = temp_main_discrete.row(s).cols(0, C_s - 1).t();

            double precision_ss = temp_marginal(s, s);
            for(int c = 0; c < C_s; ++c) {
                main_param(c) += static_cast<double>((c + 1) * (c + 1)) * precision_ss;
            }

            arma::vec bound = main_param(C_s - 1) + static_cast<double>(C_s) * rest;
            bound = arma::max(bound, arma::zeros<arma::vec>(bound.n_elem));

            compute_logZ_and_probs_ordinal_into(
                main_param, rest, bound, C_s, logz_out_, logz_scratch_
            );

            logp -= arma::accu(logz_out_.log_Z);

            // Main-effect gradient
            for(int c = 0; c < C_s; ++c) {
                grad(main_effects_discrete_offset + c) -= arma::accu(logz_out_.probs.col(c + 1));
            }

            arma::vec weights = arma::regspace<arma::vec>(1, C_s);
            arma::vec E = logz_out_.probs.cols(1, C_s) * weights;

            // Pairwise discrete gradient (ALL edges, no gating)
            arma::vec pw_grad = discrete_observations_dbl_t_ * E;
            for(size_t t = 0; t < p_; ++t) {
                if(s == t) continue;
                size_t loc = (s < t) ? kxx_idx(s, t) : kxx_idx(t, s);
                grad(loc) -= 2.0 * pw_grad(t);
            }

            arma::vec weights_sq = arma::square(weights);
            arma::vec E_sq = logz_out_.probs.cols(1, C_s) * weights_sq;

            arma::vec diff_pw = discrete_observations_dbl_t_ *
                (discrete_observations_dbl_.col(s) - E);
            diff_pw(s) = 0.0;

            double diff_diag = arma::dot(
                discrete_observations_dbl_.col(s),
                discrete_observations_dbl_.col(s)) - arma::accu(E_sq);

            double sum_obs_minus_E = arma::accu(discrete_observations_dbl_.col(s)) - arma::accu(E);

            for(size_t t = 0; t < p_; ++t) {
                if(t != s) Theta_bar(s, t) += 2.0 * diff_pw(t);
            }
            Theta_bar(s, s) += diff_diag;

            arma::rowvec cross_self = 4.0 * (diff_pw.t() * temp_pairwise_cross) * temp_covariance
                                  + 4.0 * diff_diag * cross_times_cov.row(s)
                                  + 2.0 * sum_obs_minus_E * temp_main_continuous.t();

            for(size_t j = 0; j < q_; ++j) {
                grad(kxy_idx(s, j)) += cross_self(j);
            }

            arma::rowvec V_s = 4.0 * cross_times_cov.row(s);
            for(size_t t = 0; t < p_; ++t) {
                if(t == s || std::abs(diff_pw(t)) < 1e-300) continue;
                for(size_t j = 0; j < q_; ++j) {
                    grad(kxy_idx(t, j)) += diff_pw(t) * V_s(j);
                }
            }

            for(size_t j = 0; j < q_; ++j) {
                grad(mean_offset + j) += 2.0 * temp_pairwise_cross(s, j) * sum_obs_minus_E;
            }

            main_effects_discrete_offset += C_s;
        } else {
            // Blume-Capel variable
            int ref = baseline_category_(s);
            double lin_eff = temp_main_discrete(s, 0);
            double quad_eff = temp_main_discrete(s, 1);

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

            grad(main_effects_discrete_offset)     -= arma::accu(logz_out_.probs * score);
            grad(main_effects_discrete_offset + 1) -= arma::accu(logz_out_.probs * sq_score);

            arma::vec E = logz_out_.probs * score;

            arma::vec pw_grad = discrete_observations_dbl_t_ * E;
            for(size_t t = 0; t < p_; ++t) {
                if(s == t) continue;
                size_t loc = (s < t) ? kxx_idx(s, t) : kxx_idx(t, s);
                grad(loc) -= 2.0 * pw_grad(t);
            }

            arma::vec E_sq = logz_out_.probs * sq_score;

            arma::vec diff_pw = discrete_observations_dbl_t_ *
                (discrete_observations_dbl_.col(s) - E);
            diff_pw(s) = 0.0;

            double diff_diag = arma::dot(
                discrete_observations_dbl_.col(s),
                discrete_observations_dbl_.col(s)) - arma::accu(E_sq);

            double sum_obs_minus_E = arma::accu(discrete_observations_dbl_.col(s)) - arma::accu(E);

            for(size_t t = 0; t < p_; ++t) {
                if(t != s) Theta_bar(s, t) += 2.0 * diff_pw(t);
            }
            Theta_bar(s, s) += diff_diag;

            arma::rowvec cross_self = 4.0 * (diff_pw.t() * temp_pairwise_cross) * temp_covariance
                                  + 4.0 * diff_diag * cross_times_cov.row(s)
                                  + 2.0 * sum_obs_minus_E * temp_main_continuous.t();

            for(size_t j = 0; j < q_; ++j) {
                grad(kxy_idx(s, j)) += cross_self(j);
            }

            arma::rowvec V_s = 4.0 * cross_times_cov.row(s);
            for(size_t t = 0; t < p_; ++t) {
                if(t == s || std::abs(diff_pw(t)) < 1e-300) continue;
                for(size_t j = 0; j < q_; ++j) {
                    grad(kxy_idx(t, j)) += diff_pw(t) * V_s(j);
                }
            }

            for(size_t j = 0; j < q_; ++j) {
                grad(mean_offset + j) += 2.0 * temp_pairwise_cross(s, j) * sum_obs_minus_E;
            }

            main_effects_discrete_offset += 2;
        }
    }

    // Numerator contribution to logp
    main_effects_discrete_offset = 0;
    for(size_t s = 0; s < p_; ++s) {
        int C_s = num_categories_(s);
        arma::vec rest;
        double precision_ss = temp_marginal(s, s);
        rest = 2.0 * (X_marginal.col(s)
                    - discrete_observations_dbl_.col(s) * precision_ss)
             + cross_bias(s);
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
    }

    // =========================================================================
    // Part 2: GGM conditional
    // =========================================================================

    double quad_sum = arma::accu((D * temp_precision) % D);
    logp += static_cast<double>(n_) / 2.0 *
            (-static_cast<double>(q_) * MY_LOG(2.0 * arma::datum::pi)
             + temp_log_det)
          - quad_sum / 2.0;

    arma::vec D_colsums = arma::sum(D, 0).t();
    arma::vec grad_mean_ggm = temp_precision * D_colsums;
    for(size_t j = 0; j < q_; ++j) {
        grad(mean_offset + j) += grad_mean_ggm(j);
    }

    arma::mat grad_cross_ggm = 2.0 * discrete_observations_dbl_t_ * D;
    for(size_t i = 0; i < p_; ++i) {
        for(size_t j = 0; j < q_; ++j) {
            grad(kxy_idx(i, j)) += grad_cross_ggm(i, j);
        }
    }

    // =========================================================================
    // Part 3: Priors
    // =========================================================================

    // Main effects priors
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

    // Kxx priors. Gated on edge_indicators_: inactive edges have K_xx_{ij} = 0
    // (point mass) and contribute no slab density. Matches the GGM convention.
    for(size_t i = 0; i < p_ - 1; ++i) {
        for(size_t j = i + 1; j < p_; ++j) {
            if(edge_indicators_(i, j) == 0) continue;
            double val = temp_pairwise_discrete(i, j);
            logp += interaction_prior_->logp(val);
            grad(kxx_idx(i, j)) += interaction_prior_->grad(val);
        }
    }

    // Continuous mean priors: Normal(0, 1)
    for(size_t j = 0; j < q_; ++j) {
        double val = temp_main_continuous(j);
        logp += means_prior_->logp(val);
        grad(mean_offset + j) += means_prior_->grad(val);
    }

    // Kxy priors. Gated on edge_indicators_: inactive edges have K_xy_{ij} = 0
    // (point mass) and contribute no slab density.
    for(size_t i = 0; i < p_; ++i) {
        for(size_t j = 0; j < q_; ++j) {
            if(edge_indicators_(i, p_ + j) == 0) continue;
            double val = temp_pairwise_cross(i, j);
            logp += interaction_prior_->logp(val);
            grad(kxy_idx(i, j)) += interaction_prior_->grad(val);
        }
    }

    // =========================================================================
    // Part 4: Precision gradient via Cholesky parameterization
    // =========================================================================

    arma::mat Omega_bar = 0.5 * static_cast<double>(n_) * temp_covariance
                        - 0.5 * D.t() * D;

    Omega_bar -= 2.0 * temp_covariance * temp_pairwise_cross.t()
               * discrete_observations_dbl_t_ * D;

    Omega_bar -= 2.0 * temp_covariance * temp_pairwise_cross.t()
               * Theta_bar * temp_pairwise_cross * temp_covariance;

    // Prior on the partial-association diagonal: -K_yy_{jj} = Theta_{jj}/2.
    //   d/dTheta_jj log p(Theta_jj/2) = 0.5 * grad(Theta_jj/2).
    for(size_t j = 0; j < q_; ++j) {
        double half_kjj = 0.5 * temp_precision(j, j);
        logp += diagonal_prior_->logp(half_kjj);
        Omega_bar(j, j) += 0.5 * diagonal_prior_->grad(half_kjj);
    }
    // Interaction prior on off-diagonal Kyy_{ij} = -Ω_{ij}/2.
    // Gated on edge_indicators_: inactive edges have K_yy_{ij} = 0 (point mass)
    // and contribute no slab density, matching the GGM convention at
    // ggm_gradient.cpp where the slab is summed over included_indices only.
    // Chain rule: ∂logπ/∂Ω_{ij} = ∂logπ/∂Kyy_{ij} · (-1/2)
    for(size_t i = 0; i < q_ - 1; ++i) {
        for(size_t j = i + 1; j < q_; ++j) {
            if(edge_indicators_(p_ + i, p_ + j) == 0) continue;
            double kyy_val = -0.5 * temp_precision(i, j);
            logp += interaction_prior_->logp(kyy_val);
            Omega_bar(i, j) += -0.5 * interaction_prior_->grad(kyy_val);
        }
    }

    // R̄ = R (Ω̄ + Ω̄ᵀ)
    arma::mat Omega_bar_sym = Omega_bar + Omega_bar.t();
    arma::mat R_bar = temp_cholesky * Omega_bar_sym;

    // Cholesky-to-K Jacobian (graph-agnostic) + mass-weighted per-column
    // Pfaffian correction for the RATTLE manifold marginal:
    //   ldj      = q*log(2) + Σ_j (q+1-j) ψ_j
    //   pfaffian = 0.5 * Σ_qq log det(A_qq diag(M_qq^{-1}) A_qq^T)
    //   logp    += ldj - pfaffian
    // Mass diagonal is plumbed from this->inv_mass_ (empty ⇒ identity).
    // For the full Kyy graph, all A_qq are empty and Pfaffian = 0.
    logp += static_cast<double>(q_) * std::log(2.0);
    for(size_t j = 0; j < q_; ++j) {
        logp += static_cast<double>(q_ + 1 - j) * std::log(temp_cholesky(j, j));
    }

    // Determinant tilt on the Kyy block: adds delta * log|Kyy| = 2*delta * sum(psi)
    // to the log-prior. Pushes the continuous-block precision matrix away from
    // the PD-cone boundary. delta = 0 recovers the untilted target.
    logp += determinant_tilt_yy_ * temp_log_det;

    const auto& cs = chol_constraint_structure_;
    const bool identity_mass = inv_mass_.is_empty();
    arma::mat Aq_buf;
    std::vector<arma::mat> G_chol(q_);
    std::vector<arma::mat> Aq_cache(q_);
    std::vector<arma::vec> inv_mass_q_cache(q_);
    double pfaffian = 0.0;
    for(size_t col = 1; col < q_; ++col) {
        const auto& cc = cs.columns[col];
        if(cc.m_q == 0) continue;

        GGMGradientEngine::build_Aq(temp_cholesky, cc, col, Aq_buf);
        Aq_cache[col] = Aq_buf;

        arma::vec inv_mass_q(col);
        if(identity_mass) {
            inv_mass_q.ones();
        } else {
            size_t off_q = chol_block_offset_ + cs.full_theta_offsets[col];
            for(size_t l = 0; l < col; ++l) {
                inv_mass_q(l) = inv_mass_(off_q + l);
            }
        }
        inv_mass_q_cache[col] = inv_mass_q;

        // G_q = A_q diag(inv_mass_q) A_q^T
        arma::mat Aq_scaled = Aq_buf;
        Aq_scaled.each_row() %= inv_mass_q.t();
        arma::mat G_q = Aq_scaled * Aq_buf.t();

        arma::mat L_q;
        bool chol_ok = arma::chol(L_q, G_q, "lower");
        if(!chol_ok) {
            double ridge = 1e-12 * (arma::trace(G_q) /
                                    static_cast<double>(cc.m_q) + 1.0);
            chol_ok = arma::chol(L_q, G_q + ridge * arma::eye(cc.m_q, cc.m_q),
                                 "lower");
            if(!chol_ok) {
                return {-std::numeric_limits<double>::infinity(),
                        arma::vec(full_dim, arma::fill::zeros)};
            }
        }
        G_chol[col] = L_q;
        pfaffian += arma::accu(arma::log(arma::diagvec(L_q)));
    }
    logp -= pfaffian;

    // Pfaffian adjoint: d/dA_q [-0.5 log det(A_q M_q^{-1} A_q^T)] flows back
    // to R_bar at the excluded-edge columns. dA_q = G_q^{-1} A_q · diag(M_q^{-1}).
    for(size_t col = 1; col < q_; ++col) {
        const auto& cc = cs.columns[col];
        if(cc.m_q == 0) continue;

        const arma::mat& L_q = G_chol[col];
        const arma::mat& Aq = Aq_cache[col];
        const arma::vec& inv_mass_q = inv_mass_q_cache[col];

        arma::mat Z = arma::solve(arma::trimatl(L_q), Aq,
                                  arma::solve_opts::fast);
        Z = arma::solve(arma::trimatu(L_q.t()), Z,
                        arma::solve_opts::fast);

        arma::mat dAq = Z;
        dAq.each_row() %= inv_mass_q.t();

        for(size_t r = 0; r < cc.m_q; ++r) {
            size_t i_r = cc.excluded_indices[r];
            for(size_t l = 0; l <= i_r; ++l) {
                R_bar(l, i_r) -= dAq(r, l);
            }
        }
    }

    // Extract position gradient from R_bar with the unified weight (q+1-j)
    // on the diagonal-psi entries. The determinant tilt adds +2*delta to the
    // psi-bar entries (d/dpsi [delta * 2 * sum(psi)] = 2 * delta).
    size_t gidx = chol_offset;
    for(size_t j = 0; j < q_; ++j) {
        double w_j = static_cast<double>(q_ + 1 - j);
        for(size_t i = 0; i < j; ++i) {
            grad(gidx++) = R_bar(i, j);
        }
        grad(gidx++) = R_bar(j, j) * temp_cholesky(j, j) + w_j
                       + 2.0 * determinant_tilt_yy_;
    }

    return {logp, grad};
}

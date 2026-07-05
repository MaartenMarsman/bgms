// mixed_mrf_likelihoods.cpp — MixedMRFModel likelihood terms.
//
// The two log-likelihood contributions of the mixed model: log_conditional_ggm
// (continuous block conditional on the discrete state) and log_marginal_omrf
// (discrete block). Consumed by the gradient and Metropolis update bodies in the
// sibling translation units.
#include <RcppArmadillo.h>
#include "models/mixed/mixed_mrf_model.h"
#include "utils/variable_helpers.h"
#include "math/explog_macros.h"


// =============================================================================
// log_marginal_omrf
// =============================================================================
// Marginal OMRF pseudolikelihood for discrete variable s:
//   log f(x_s | x_{-s}) using M = A_xx + 2 A_xy Σ_yy A_xy'
//
// After integrating y out of the joint density L(x,y) = m_x(x) + x'A_xx x
// + 2 x'A_xy y + b_y' y + y'A_yy y with A_yy = -Λ/2, the marginal log-density
// is m_x(x) + x'M x + 2 (A_xy μ_y)' x + const.  Reading off the x_s-conditional:
//
//   log p(x_s=c | x_{-s}) ∝ main_x_s(c) + c² M_ss + c · rest_s,
//   rest_s = 2 Σ_{j≠s} M_{sj} x_j + 2 (A_xy μ_y)_s.
// =============================================================================

double MixedMRFModel::log_marginal_omrf(int s) const {
    // Rest score: 2 · M · x minus self-interaction, plus cross-bias.
    // Factor 2 from x'Mx derivative.
    double precision_ss = marginal_interactions_(s, s);
    arma::vec rest = 2.0 * (discrete_observations_dbl_ * marginal_interactions_.col(s)
                          - discrete_observations_dbl_.col(s) * precision_ss)
                   + 2.0 * arma::dot(pairwise_effects_cross_.row(s), main_effects_continuous_);
    return log_marginal_omrf_given_rest(s, rest, precision_ss);
}

double MixedMRFModel::log_marginal_omrf_from(int s, const arma::mat& matvec,
                                             double precision_ss, double bias_s) const {
    // Same rest score as log_marginal_omrf, with the O(np) matrix-vector
    // product replaced by a cached X · M column.
    arma::vec rest = 2.0 * (matvec.col(s) - discrete_observations_dbl_.col(s) * precision_ss)
                   + bias_s;
    return log_marginal_omrf_given_rest(s, rest, precision_ss);
}

double MixedMRFModel::log_marginal_omrf_cached(int s) const {
    return log_marginal_omrf_from(s, marginal_matvec_,
                                  marginal_interactions_(s, s), cross_bias_(s));
}

double MixedMRFModel::log_marginal_omrf_given_rest(
    int s, const arma::vec& rest, double precision_ss) const {
    int C_s = num_categories_(s);

    // Numerator: dot(x_s, rest) + precision_ss * dot(x_s, x_s) + main effects
    double numer = arma::dot(discrete_observations_dbl_.col(s), rest)
                 + precision_ss * arma::dot(discrete_observations_dbl_.col(s),
                                        discrete_observations_dbl_.col(s));

    if(is_ordinal_variable_(s)) {
        for(int c = 1; c <= C_s; ++c) {
            numer += static_cast<double>(counts_per_category_(c, s)) * main_effects_discrete_(s, c - 1);
        }

        // Denominator: main_param(c) = μ_x(s,c) + (c+1)^2 Θ_ss
        arma::vec main_param(C_s);
        for(int c = 0; c < C_s; ++c) {
            main_param(c) = main_effects_discrete_(s, c) + static_cast<double>((c + 1) * (c + 1)) * precision_ss;
        }

        arma::vec bound = arma::clamp(
            static_cast<double>(C_s) * rest, 0.0, arma::datum::inf);
        arma::vec denom = compute_denom_ordinal(rest, main_param, bound);

        return numer - arma::accu(bound + ARMA_MY_LOG(denom));
    } else {
        // Blume-Capel: alpha * sum(x) + beta * sum(x^2)
        double alpha = main_effects_discrete_(s, 0);
        double beta = main_effects_discrete_(s, 1);
        numer += alpha * static_cast<double>(blume_capel_stats_(0, s))
               + beta * static_cast<double>(blume_capel_stats_(1, s));

        // Denominator: theta_c includes marginal_interactions_(s,s) * (c - ref)^2
        int ref = baseline_category_(s);
        double effective_beta = beta + precision_ss;

        arma::vec bound;
        arma::vec denom = compute_denom_blume_capel(
            rest, alpha, effective_beta, ref, C_s, bound
        );

        return numer - arma::accu(bound + ARMA_MY_LOG(denom));
    }
}


// =============================================================================
// log_conditional_ggm
// =============================================================================
// Conditional GGM log-likelihood: log f(y | x)
//   y | x ~ N(conditional_mean_, covariance_continuous_)
//
// Uses cached covariance_continuous_, log_det_precision_, and conditional_mean_.
// The quadratic form uses precision = -2 * pairwise_effects_continuous_.
// =============================================================================

double MixedMRFModel::log_conditional_ggm() const {
    arma::mat D = continuous_observations_ - conditional_mean_;

    // Quadratic form: trace(Precision D'D)
    double quad_sum = arma::accu((D * (-2.0 * pairwise_effects_continuous_)) % D);

    return static_cast<double>(n_) / 2.0 *
           (-static_cast<double>(q_) * MY_LOG(2.0 * arma::datum::pi)
            + log_det_precision_)
         - quad_sum / 2.0;
}


// =============================================================================
// omrf_ratio_for_covariance_change
// =============================================================================
// A Kyy proposal changes Σ, which moves the whole marginal interaction matrix
// M = A_xx + 2 A_xy Σ A_xy'. The change factors through ΔΣ:
//   ΔM = 2 A_xy ΔΣ A_xy',   X · ΔM = (X A_xy) (ΔΣ A_xy') = cross_matvec_ · E
// so the proposed rest scores follow from the cached matvecs in O(nq(q+p))
// instead of rebuilding M and its p rest-score products.
// =============================================================================

double MixedMRFModel::omrf_ratio_for_covariance_change(const arma::mat& cov_prop) {
    arma::mat delta_sigma = cov_prop - covariance_continuous_;
    arma::mat E = delta_sigma * pairwise_effects_cross_.t();  // q x p

    marginal_matvec_prop_ = marginal_matvec_ + 2.0 * (cross_matvec_ * E);
    for(size_t s = 0; s < p_; ++s) {
        mdiag_prop_(s) = marginal_interactions_(s, s)
                       + 2.0 * arma::dot(pairwise_effects_cross_.row(s), E.col(s));
        ll_marginal_prop_(s) = log_marginal_omrf_from(
            s, marginal_matvec_prop_, mdiag_prop_(s), cross_bias_(s));
    }
    return arma::accu(ll_marginal_prop_) - arma::accu(ll_marginal_cache_);
}

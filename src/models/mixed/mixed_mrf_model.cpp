// mixed_mrf_model.cpp — MixedMRFModel core.
//
// Construction, sufficient statistics, parameter (de)vectorization, the
// do_one_metropolis_step driver, missing-data imputation, proposal-SD tuning,
// and the position/momentum projection cache for constrained NUTS. The likelihood,
// gradient, and within-model MH update bodies live in the sibling translation units
// (mixed_mrf_likelihoods.cpp, mixed_mrf_gradient.cpp, mixed_mrf_metropolis.cpp).
#include <RcppArmadillo.h>
#include <utility>
#include "models/mixed/mixed_mrf_model.h"
#include "mcmc/execution/chain_result.h"
#include "math/explog_macros.h"
#include "rng/rng_utils.h"
#include "mcmc/execution/warmup_schedule.h"


// =============================================================================
// Constructor
// =============================================================================

MixedMRFModel::MixedMRFModel(
    const arma::imat& discrete_observations,
    const arma::mat& continuous_observations,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const arma::mat& inclusion_probability,
    const arma::imat& initial_edge_indicators,
    bool edge_selection,
    std::unique_ptr<BaseParameterPrior> interaction_prior,
    std::unique_ptr<BaseParameterPrior> threshold_prior,
    std::unique_ptr<BaseParameterPrior> means_prior,
    std::unique_ptr<BaseParameterPrior> diagonal_prior,
    int seed
) :
    n_(discrete_observations.n_rows),
    p_(discrete_observations.n_cols),
    q_(continuous_observations.n_cols),
    discrete_observations_(discrete_observations),
    continuous_observations_(continuous_observations),
    num_categories_(num_categories),
    is_ordinal_variable_(is_ordinal_variable),
    baseline_category_(baseline_category),
    edge_indicators_(initial_edge_indicators),
    inclusion_probability_(inclusion_probability),
    edge_selection_(edge_selection),
    edge_selection_active_(false),
    interaction_prior_(std::move(interaction_prior)),
    threshold_prior_(std::move(threshold_prior)),
    means_prior_(std::move(means_prior)),
    diagonal_prior_(std::move(diagonal_prior)),
    rng_(seed)
{
    // Dimension counts
    num_main_ = count_num_main_effects();
    num_pairwise_xx_ = (p_ * (p_ - 1)) / 2;
    num_pairwise_yy_ = (q_ * (q_ - 1)) / 2;
    num_cross_ = p_ * q_;
    num_cholesky_ = (q_ * (q_ + 1)) / 2;

    max_cats_ = num_categories_.max();

    // Center Blume-Capel observations at baseline category so that all
    // downstream code operates in a shifted coordinate system where the
    // reference corresponds to zero (same convention as OMRFModel).
    for(size_t s = 0; s < p_; ++s) {
        if(!is_ordinal_variable_(s)) {
            discrete_observations_.col(s) -= baseline_category_(s);
        }
    }
    discrete_observations_dbl_ = arma::conv_to<arma::mat>::from(discrete_observations_);
    discrete_observations_dbl_t_ = discrete_observations_dbl_.t();

    // Compute sufficient statistics
    compute_sufficient_statistics();

    // Initialize parameters to zero
    main_effects_discrete_ = arma::zeros<arma::mat>(p_, max_cats_);
    main_effects_continuous_ = arma::zeros<arma::vec>(q_);
    pairwise_effects_discrete_ = arma::zeros<arma::mat>(p_, p_);
    pairwise_effects_continuous_ = -0.5 * arma::eye<arma::mat>(q_, q_);
    pairwise_effects_cross_ = arma::zeros<arma::mat>(p_, q_);

    // Initialize proposal SDs
    proposal_sd_main_discrete_ = arma::ones<arma::mat>(p_, max_cats_);
    proposal_sd_main_continuous_ = arma::ones<arma::mat>(q_, 1);
    proposal_sd_pairwise_discrete_ = arma::ones<arma::mat>(p_, p_);
    proposal_sd_pairwise_continuous_ = arma::ones<arma::mat>(q_, q_);
    proposal_sd_pairwise_cross_ = arma::ones<arma::mat>(p_, q_);

    // Initialize precision caches (precision starts as identity)
    cholesky_of_precision_ = arma::eye<arma::mat>(q_, q_);
    inv_cholesky_of_precision_ = arma::eye<arma::mat>(q_, q_);
    covariance_continuous_ = arma::eye<arma::mat>(q_, q_);
    log_det_precision_ = 0.0;

    // Rank-1 Cholesky update workspace
    precision_proposal_ = arma::mat(q_, q_, arma::fill::none);
    cont_vf1_ = arma::zeros<arma::vec>(q_);
    cont_vf2_ = arma::zeros<arma::vec>(q_);
    cont_u1_ = arma::zeros<arma::vec>(q_);
    cont_u2_ = arma::zeros<arma::vec>(q_);

    // Initialize conditional mean: M = μ_y' + 2 X cross_int Sigma_yy
    //   With cross_int = 0 and precision = I, this reduces to 0.
    conditional_mean_ = arma::zeros<arma::mat>(n_, q_);

    // Initialize marginal interactions: disc_int + 2 * cross_int * Sigma_yy * cross_int'
    //   With cross_int = 0, this is zero.
    marginal_interactions_ = arma::zeros<arma::mat>(p_, p_);
    cross_term_ = arma::zeros<arma::mat>(p_, p_);

    // Initialize edge-order permutation vectors. The counts are size_t, so a
    // zero count (one continuous or one discrete variable) must not reach
    // regspace(0, count - 1), which would underflow to SIZE_MAX.
    edge_order_xx_ = num_pairwise_xx_
        ? arma::regspace<arma::uvec>(0, num_pairwise_xx_ - 1) : arma::uvec();
    edge_order_yy_ = num_pairwise_yy_
        ? arma::regspace<arma::uvec>(0, num_pairwise_yy_ - 1) : arma::uvec();
    edge_order_xy_ = num_cross_
        ? arma::regspace<arma::uvec>(0, num_cross_ - 1) : arma::uvec();

    // Flat-index -> (i, j) lookup tables for the upper-triangle edge loops
    // (row-major: (0,1),(0,2),...,(1,2),...).
    auto build_edge_pairs = [](size_t dim, size_t count) {
        arma::umat pairs(count, 2);
        size_t flat = 0;
        for (size_t i = 0; i + 1 < dim; ++i) {
            for (size_t j = i + 1; j < dim; ++j, ++flat) {
                pairs(flat, 0) = i;
                pairs(flat, 1) = j;
            }
        }
        return pairs;
    };
    edge_pairs_xx_ = build_edge_pairs(p_, num_pairwise_xx_);
    edge_pairs_yy_ = build_edge_pairs(q_, num_pairwise_yy_);

    // Detect sparse initial graph (constraints without edge selection)
    if(!edge_selection_) {
        size_t max_edges = num_pairwise_xx_ + num_pairwise_yy_ + num_cross_;
        size_t num_edges = 0;
        for(size_t i = 0; i < p_ + q_; ++i)
            for(size_t j = i + 1; j < p_ + q_; ++j)
                if(edge_indicators_(i, j) == 1) num_edges++;
        has_sparse_graph_ = (num_edges < max_edges);
    }
}


// =============================================================================
// Copy constructor
// =============================================================================

MixedMRFModel::MixedMRFModel(const MixedMRFModel& other)
    : BaseModel(other),
      target_accept_(other.target_accept_),
      determinant_tilt_yy_(other.determinant_tilt_yy_),
      n_(other.n_),
      p_(other.p_),
      q_(other.q_),
      num_main_(other.num_main_),
      num_pairwise_xx_(other.num_pairwise_xx_),
      num_pairwise_yy_(other.num_pairwise_yy_),
      num_cross_(other.num_cross_),
      num_cholesky_(other.num_cholesky_),
      discrete_observations_(other.discrete_observations_),
      discrete_observations_dbl_(other.discrete_observations_dbl_),
      continuous_observations_(other.continuous_observations_),
      num_categories_(other.num_categories_),
      max_cats_(other.max_cats_),
      is_ordinal_variable_(other.is_ordinal_variable_),
      baseline_category_(other.baseline_category_),
      missing_index_discrete_(other.missing_index_discrete_),
      missing_index_continuous_(other.missing_index_continuous_),
      has_missing_(other.has_missing_),
      counts_per_category_(other.counts_per_category_),
      blume_capel_stats_(other.blume_capel_stats_),
      main_effects_discrete_(other.main_effects_discrete_),
      main_effects_continuous_(other.main_effects_continuous_),
      pairwise_effects_discrete_(other.pairwise_effects_discrete_),
      pairwise_effects_continuous_(other.pairwise_effects_continuous_),
      pairwise_effects_cross_(other.pairwise_effects_cross_),
      edge_indicators_(other.edge_indicators_),
      inclusion_probability_(other.inclusion_probability_),
      edge_selection_(other.edge_selection_),
      edge_selection_active_(other.edge_selection_active_),
      interaction_prior_(other.interaction_prior_->clone()),
      threshold_prior_(other.threshold_prior_->clone()),
      means_prior_(other.means_prior_->clone()),
      diagonal_prior_(other.diagonal_prior_->clone()),
      proposal_sd_main_discrete_(other.proposal_sd_main_discrete_),
      proposal_sd_main_continuous_(other.proposal_sd_main_continuous_),
      proposal_sd_pairwise_discrete_(other.proposal_sd_pairwise_discrete_),
      proposal_sd_pairwise_continuous_(other.proposal_sd_pairwise_continuous_),
      proposal_sd_pairwise_cross_(other.proposal_sd_pairwise_cross_),
      cholesky_of_precision_(other.cholesky_of_precision_),
      inv_cholesky_of_precision_(other.inv_cholesky_of_precision_),
      covariance_continuous_(other.covariance_continuous_),
      log_det_precision_(other.log_det_precision_),
      marginal_interactions_(other.marginal_interactions_),
      cross_term_(other.cross_term_),
      conditional_mean_(other.conditional_mean_),
      cont_constants_(other.cont_constants_),
      precision_proposal_(other.precision_proposal_),
      cont_v1_(other.cont_v1_),
      cont_v2_(other.cont_v2_),
      cont_vf1_(other.cont_vf1_),
      cont_vf2_(other.cont_vf2_),
      cont_u1_(other.cont_u1_),
      cont_u2_(other.cont_u2_),
      discrete_observations_dbl_t_(other.discrete_observations_dbl_t_),
      gradient_cache_valid_(false),
      chol_constraint_structure_(other.chol_constraint_structure_),
      chol_block_offset_(other.chol_block_offset_),
      constraint_dirty_(other.constraint_dirty_),
      has_sparse_graph_(other.has_sparse_graph_),
      rng_(other.rng_),
      edge_order_xx_(other.edge_order_xx_),
      edge_order_yy_(other.edge_order_yy_),
      edge_order_xy_(other.edge_order_xy_),
      edge_pairs_xx_(other.edge_pairs_xx_),
      edge_pairs_yy_(other.edge_pairs_yy_)
{
    // The gradient engine holds a pointer to the source's constraint
    // structure; force a rebuild so the clone binds to its own copy.
    constraint_dirty_ = true;
    theta_yy_valid_ = false;

    // Deep-copy the Z-ratio engine (per-chain caches never cross threads)
    // and rebind its RNG to this clone's stream.
    if (other.zratio_engine_) {
        zratio_engine_ = std::make_shared<ZRatioEngine>(*other.zratio_engine_);
        zratio_engine_->set_rng(&rng_);
    }
}


void MixedMRFModel::collect_chain_diagnostics(ChainResult& chain_result) const {
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


// =============================================================================
// Sufficient statistics
// =============================================================================

void MixedMRFModel::compute_sufficient_statistics() {
    // Category counts for ordinal variables
    counts_per_category_ = arma::zeros<arma::imat>(max_cats_ + 1, p_);
    for(size_t s = 0; s < p_; ++s) {
        if(is_ordinal_variable_(s)) {
            for(size_t i = 0; i < n_; ++i) {
                int cat = discrete_observations_(i, s);
                if(cat >= 0 && cat <= num_categories_(s)) {
                    counts_per_category_(cat, s)++;
                }
            }
        }
    }

    // Blume-Capel statistics (linear and quadratic sums of centered obs)
    blume_capel_stats_ = arma::zeros<arma::imat>(2, p_);
    for(size_t s = 0; s < p_; ++s) {
        if(!is_ordinal_variable_(s)) {
            for(size_t i = 0; i < n_; ++i) {
                int val = discrete_observations_(i, s);  // already centered
                blume_capel_stats_(0, s) += val;
                blume_capel_stats_(1, s) += val * val;
            }
        }
    }
}


size_t MixedMRFModel::count_num_main_effects() const {
    size_t count = 0;
    for(size_t s = 0; s < p_; ++s) {
        if(is_ordinal_variable_(s)) {
            count += num_categories_(s);
        } else {
            count += 2;  // linear α and quadratic β
        }
    }
    return count;
}


// =============================================================================
// Cache maintenance
// =============================================================================

void MixedMRFModel::recompute_conditional_mean() {
    // M = μ_y' + 2 X A_xy Σ_yy
    conditional_mean_ = 2.0 * discrete_observations_dbl_ * pairwise_effects_cross_ * covariance_continuous_;
    conditional_mean_.each_row() += main_effects_continuous_.t();
}

void MixedMRFModel::recompute_pairwise_effects_continuous_decomposition() {
    // Cholesky on precision = -2 * pairwise_effects_continuous_
    arma::mat precision = -2.0 * pairwise_effects_continuous_;
    cholesky_of_precision_ = arma::chol(precision, "upper");
    arma::inv(inv_cholesky_of_precision_, arma::trimatu(cholesky_of_precision_));
    covariance_continuous_ = inv_cholesky_of_precision_ * inv_cholesky_of_precision_.t();
    log_det_precision_ = cholesky_helpers::get_log_det(cholesky_of_precision_);
}

void MixedMRFModel::recompute_marginal_interactions() {
    // Marginal PL effective interaction after integrating y from the joint:
    //   M = A_xx + 2 A_xy Σ_yy A_xy'
    // The log-marginal has x'Mx, so the x_s-conditional rest score carries
    // a factor 2 on M.col(s) (consumed at the call site in log_marginal_omrf).
    cross_term_ = 2.0 * pairwise_effects_cross_ * covariance_continuous_ * pairwise_effects_cross_.t();
    marginal_interactions_ = pairwise_effects_discrete_ + cross_term_;
}

void MixedMRFModel::refresh_marginal_interactions_entry(int i, int j) {
    // A_xx-only change: the cross term is untouched, so only the (i,j)/(j,i)
    // entries of M move.
    double value = pairwise_effects_discrete_(i, j) + cross_term_(i, j);
    marginal_interactions_(i, j) = value;
    marginal_interactions_(j, i) = value;
}

void MixedMRFModel::recompute_am_caches() {
    // Full refresh: the accept paths maintain marginal_interactions_,
    // cross_term_, the matvecs, conditional_mean_, and the log-likelihood
    // values by low-rank updates within a sweep; recomputing everything from
    // the parameters here bounds the accumulated floating-point drift to one
    // sweep.
    recompute_marginal_interactions();
    marginal_matvec_ = discrete_observations_dbl_ * marginal_interactions_;
    cross_matvec_ = discrete_observations_dbl_ * pairwise_effects_cross_;
    cross_bias_ = 2.0 * pairwise_effects_cross_ * main_effects_continuous_;
    recompute_conditional_mean_from_cross_matvec();

    if(ll_marginal_cache_.n_elem != p_) {
        ll_marginal_cache_.set_size(p_);
        marginal_matvec_prop_.set_size(n_, p_);
        mdiag_prop_.set_size(p_);
        ll_marginal_prop_.set_size(p_);
        cross_bias_prop_.set_size(p_);
        matvec_col_i_scratch_.set_size(n_);
        matvec_col_j_scratch_.set_size(n_);
    }
    for(size_t s = 0; s < p_; ++s)
        ll_marginal_cache_(s) = log_marginal_omrf_cached(s);
    ll_ggm_cache_ = log_conditional_ggm();
}

void MixedMRFModel::recompute_conditional_mean_from_cross_matvec() {
    // M = μ_y' + 2 X A_xy Σ_yy with X A_xy read from the sweep cache.
    conditional_mean_ = 2.0 * cross_matvec_ * covariance_continuous_;
    conditional_mean_.each_row() += main_effects_continuous_.t();
}

void MixedMRFModel::adopt_kyy_proposal_caches(double ggm_ratio, bool rank2) {
    // The proposal scratch holds the accepted per-variable state.
    std::swap(marginal_matvec_, marginal_matvec_prop_);
    std::swap(ll_marginal_cache_, ll_marginal_prop_);

    // M and the cross term move by 2 A_xy ΔΣ A_xy' = 2 A_xy E, with
    // E = ΔΣ A_xy' cached by omrf_ratio_for_covariance_change. The diagonal
    // is overwritten with mdiag_prop_ so it matches the values the adopted
    // per-variable marginals were computed with.
    arma::mat m_delta = 2.0 * (pairwise_effects_cross_ * cross_delta_scratch_);
    cross_term_ += m_delta;
    marginal_interactions_ += m_delta;
    marginal_interactions_.diag() = mdiag_prop_;
    cross_term_.diag() = mdiag_prop_ - pairwise_effects_discrete_.diag();

    // Conditional mean: ΔM = a1 s2' + a2 s1' (rank 1: a1 s1'), the same
    // factors the accepted GGM ratio was evaluated with.
    if(rank2) {
        conditional_mean_ += cont_a1_ * cont_s2_.t() + cont_a2_ * cont_s1_.t();
    } else {
        conditional_mean_ += cont_a1_ * cont_s1_.t();
    }
    ll_ggm_cache_ += ggm_ratio;
}

void MixedMRFModel::adopt_cross_proposal_caches(
    int i, int j, double delta, const arma::vec& u, double ggm_prop)
{
    std::swap(marginal_matvec_, marginal_matvec_prop_);
    std::swap(ll_marginal_cache_, ll_marginal_prop_);
    cross_bias_ = cross_bias_prop_;
    ll_ggm_cache_ = ggm_prop;

    // ΔM = 2δ (e_i u' + u e_i') + 2δ² Σ_jj e_i e_i' touches row and column i
    // of M and the cross term; the (i, i) entry is overwritten with
    // mdiag_prop_(i), matching the adopted per-variable marginals.
    marginal_interactions_.row(i) += (2.0 * delta) * u.t();
    marginal_interactions_.col(i) += (2.0 * delta) * u;
    marginal_interactions_(i, i) = mdiag_prop_(i);
    cross_term_.row(i) += (2.0 * delta) * u.t();
    cross_term_.col(i) += (2.0 * delta) * u;
    cross_term_(i, i) = mdiag_prop_(i) - pairwise_effects_discrete_(i, i);

    // X · A_xy moves in column j only; the conditional mean by the rank-1
    // shift the accepted GGM ratio was evaluated with.
    cross_matvec_.col(j) += delta * discrete_observations_dbl_.col(i);
    conditional_mean_ += (2.0 * delta) * discrete_observations_dbl_.col(i)
                       * covariance_continuous_.row(j);
}


// =============================================================================
// Constraint structure (RATTLE)
// =============================================================================

void MixedMRFModel::ensure_constraint_structure() {
    if(!constraint_dirty_) return;

    // --- Cholesky constraints (Gyy block) ---
    // Extract q x q sub-block of edge_indicators_ for the continuous-continuous edges
    arma::imat gyy_indicators(q_, q_, arma::fill::ones);
    for(size_t i = 0; i < q_; ++i) {
        for(size_t j = i + 1; j < q_; ++j) {
            int val = edge_indicators_(p_ + i, p_ + j);
            gyy_indicators(i, j) = val;
            gyy_indicators(j, i) = val;
        }
    }
    chol_constraint_structure_.build(gyy_indicators);

    // --- Cholesky block offset (full-layout offset of the Kyy theta block) ---
    chol_block_offset_ = num_main_ + num_pairwise_xx_ + q_ + num_cross_;

    yy_engine_.rebuild(chol_constraint_structure_);
    theta_yy_valid_ = false;

    constraint_dirty_ = false;
}

void MixedMRFModel::recompute_theta_yy() const {
    if (theta_yy_valid_) return;

    // Inverse of the engine forward map: psi_q = log(phi_qq) and
    // f_q = N_q^T x_q with N_q from the Givens QR of A_q^T. The constraint
    // structure is already built by ensure_constraint_structure.
    const auto& cs = chol_constraint_structure_;
    theta_yy_.set_size(cs.active_dim);

    arma::mat Aq_buf;

    for (size_t q = 0; q < q_; ++q) {
        const auto& col = cs.columns[q];
        size_t offset = cs.theta_offsets[q];

        theta_yy_(offset + col.d_q) = MY_LOG(cholesky_of_precision_(q, q));

        if (q == 0 || col.d_q == 0) continue;

        if (col.m_q == 0) {
            // Unconstrained column: f_q = x_q directly.
            for (size_t k = 0; k < col.d_q; ++k) {
                theta_yy_(offset + k) = cholesky_of_precision_(k, q);
            }
            continue;
        }

        arma::mat Q_tmp, R_tmp;
        arma::vec R_diag;
        std::vector<GivensRotation> rots_tmp;
        GGMGradientEngine::build_Aq(cholesky_of_precision_, col, q, Aq_buf);
        GGMGradientEngine::givens_qr(Aq_buf.t(), Q_tmp, R_tmp, R_diag, rots_tmp);
        arma::mat Nq = Q_tmp.cols(col.m_q, q - 1);

        arma::vec x_q = cholesky_of_precision_.col(q).head(q);
        arma::vec f_q = Nq.t() * x_q;
        for (size_t k = 0; k < col.d_q; ++k) {
            theta_yy_(offset + k) = f_q(k);
        }
    }

    theta_yy_valid_ = true;
}



// =============================================================================
// Parameter vectorization
// =============================================================================

// NUTS vectorization order (includes Cholesky of precision):
//   1. main_effects_discrete_: per-variable (ordinal: C_s thresholds; BC: 2 coefficients)
//   2. pairwise_effects_discrete_: upper-triangular, row-major  — p(p-1)/2
//   3. main_effects_continuous_: all q means
//   4. pairwise_effects_cross_: all p*q entries, row-major
//   5. Cholesky of precision: column-by-column, each column j has j off-diagonal
//      entries R_{0j},...,R_{(j-1)j} followed by ψ_j = log(R_{jj}) — q(q+1)/2
//
// Storage vectorization order (stores pairwise_effects_continuous_ = -Ω/2):
//   1–4. Same as NUTS order (Cholesky block NOT stored — A_yy entries stored instead)
//   5. pairwise_effects_continuous_: upper-triangle including diagonal — q(q+1)/2

size_t MixedMRFModel::parameter_dimension() const {
    if(constraint_dirty_) {
        const_cast<MixedMRFModel*>(this)->ensure_constraint_structure();
    }
    // Active NUTS parameters + Kyy theta block (f_q, psi_q per column)
    size_t dim = num_main_ + q_ + chol_constraint_structure_.active_dim;

    // Active pairwise_effects_discrete_ edges
    for(size_t i = 0; i < p_ - 1; ++i) {
        for(size_t j = i + 1; j < p_; ++j) {
            if(gxx(i, j)) dim++;
        }
    }

    // Active pairwise_effects_cross_ edges
    for(size_t i = 0; i < p_; ++i) {
        for(size_t j = 0; j < q_; ++j) {
            if(gxy(i, j)) dim++;
        }
    }

    return dim;
}

size_t MixedMRFModel::full_parameter_dimension() const {
    // All NUTS params + Cholesky block
    return num_main_ + num_pairwise_xx_ + q_ + num_cross_ + (q_ * (q_ + 1)) / 2;
}

size_t MixedMRFModel::storage_dimension() const {
    // All parameters including pairwise_effects_continuous_
    return num_main_ + num_pairwise_xx_ + q_ +
           (q_ * (q_ + 1)) / 2 + num_cross_;
}

arma::vec MixedMRFModel::get_vectorized_parameters() const {
    if(constraint_dirty_) {
        const_cast<MixedMRFModel*>(this)->ensure_constraint_structure();
    }
    arma::vec out(parameter_dimension());
    size_t idx = 0;

    // 1. main_effects_discrete_
    for(size_t s = 0; s < p_; ++s) {
        if(is_ordinal_variable_(s)) {
            for(int c = 0; c < num_categories_(s); ++c) {
                out(idx++) = main_effects_discrete_(s, c);
            }
        } else {
            out(idx++) = main_effects_discrete_(s, 0);
            out(idx++) = main_effects_discrete_(s, 1);
        }
    }

    // 2. pairwise_effects_discrete_ upper-triangular (included edges only)
    for(size_t i = 0; i < p_ - 1; ++i) {
        for(size_t j = i + 1; j < p_; ++j) {
            if(gxx(i, j) == 1) {
                out(idx++) = pairwise_effects_discrete_(i, j);
            }
        }
    }

    // 3. main_effects_continuous_
    for(size_t j = 0; j < q_; ++j) {
        out(idx++) = main_effects_continuous_(j);
    }

    // 4. pairwise_effects_cross_ row-major (included edges only)
    for(size_t i = 0; i < p_; ++i) {
        for(size_t j = 0; j < q_; ++j) {
            if(gxy(i, j) == 1) {
                out(idx++) = pairwise_effects_cross_(i, j);
            }
        }
    }

    // 5. Kyy theta block: (f_q, psi_q) per column via the null-space map
    recompute_theta_yy();
    out.subvec(idx, idx + chol_constraint_structure_.active_dim - 1) = theta_yy_;

    return out;
}

arma::vec MixedMRFModel::get_full_vectorized_parameters() const {
    if(constraint_dirty_) {
        const_cast<MixedMRFModel*>(this)->ensure_constraint_structure();
    }
    // All NUTS parameters + Cholesky, fixed size (inactive edges zeroed)
    arma::vec out(full_parameter_dimension(), arma::fill::zeros);
    size_t idx = 0;

    // 1. main_effects_discrete_
    for(size_t s = 0; s < p_; ++s) {
        if(is_ordinal_variable_(s)) {
            for(int c = 0; c < num_categories_(s); ++c) {
                out(idx++) = main_effects_discrete_(s, c);
            }
        } else {
            out(idx++) = main_effects_discrete_(s, 0);
            out(idx++) = main_effects_discrete_(s, 1);
        }
    }

    // 2. pairwise_effects_discrete_ upper-triangular (all entries, zeros for inactive)
    for(size_t i = 0; i < p_ - 1; ++i) {
        for(size_t j = i + 1; j < p_; ++j) {
            out(idx++) = pairwise_effects_discrete_(i, j);
        }
    }

    // 3. main_effects_continuous_
    for(size_t j = 0; j < q_; ++j) {
        out(idx++) = main_effects_continuous_(j);
    }

    // 4. pairwise_effects_cross_ row-major (all entries, zeros for inactive)
    for(size_t i = 0; i < p_; ++i) {
        for(size_t j = 0; j < q_; ++j) {
            out(idx++) = pairwise_effects_cross_(i, j);
        }
    }

    // 5. Kyy theta block scattered into the full Cholesky-slot layout
    //    (excluded off-diagonal slots stay zero)
    recompute_theta_yy();
    const auto& cs = chol_constraint_structure_;
    cs.for_each_active_full_pair([&](size_t active_idx, size_t full_idx) {
        out(idx + full_idx) = theta_yy_(active_idx);
    });

    return out;
}

arma::vec MixedMRFModel::get_storage_vectorized_parameters() const {
    // All parameters including pairwise_effects_continuous_, fixed size
    arma::vec out(storage_dimension(), arma::fill::zeros);
    size_t idx = 0;

    // 1. main_effects_discrete_
    for(size_t s = 0; s < p_; ++s) {
        if(is_ordinal_variable_(s)) {
            for(int c = 0; c < num_categories_(s); ++c) {
                out(idx++) = main_effects_discrete_(s, c);
            }
        } else {
            out(idx++) = main_effects_discrete_(s, 0);
            out(idx++) = main_effects_discrete_(s, 1);
        }
    }

    // 2. pairwise_effects_discrete_ upper-triangular
    for(size_t i = 0; i < p_ - 1; ++i) {
        for(size_t j = i + 1; j < p_; ++j) {
            out(idx++) = pairwise_effects_discrete_(i, j);
        }
    }

    // 3. main_effects_continuous_
    for(size_t j = 0; j < q_; ++j) {
        out(idx++) = main_effects_continuous_(j);
    }

    // 4. pairwise_effects_cross_ row-major
    for(size_t i = 0; i < p_; ++i) {
        for(size_t j = 0; j < q_; ++j) {
            out(idx++) = pairwise_effects_cross_(i, j);
        }
    }

    // 5. pairwise_effects_continuous_ upper-triangle including diagonal
    for(size_t i = 0; i < q_; ++i) {
        for(size_t j = i; j < q_; ++j) {
            out(idx++) = pairwise_effects_continuous_(i, j);
        }
    }

    return out;
}

void MixedMRFModel::set_vectorized_parameters(const arma::vec& params) {
    ensure_constraint_structure();
    size_t idx = 0;

    // 1. main_effects_discrete_
    for(size_t s = 0; s < p_; ++s) {
        if(is_ordinal_variable_(s)) {
            for(int c = 0; c < num_categories_(s); ++c) {
                main_effects_discrete_(s, c) = params(idx++);
            }
        } else {
            main_effects_discrete_(s, 0) = params(idx++);
            main_effects_discrete_(s, 1) = params(idx++);
        }
    }

    // 2. pairwise_effects_discrete_ upper-triangular (included edges only)
    for(size_t i = 0; i < p_ - 1; ++i) {
        for(size_t j = i + 1; j < p_; ++j) {
            if(gxx(i, j) == 1) {
                pairwise_effects_discrete_(i, j) = params(idx);
                pairwise_effects_discrete_(j, i) = params(idx);
                idx++;
            }
        }
    }

    // 3. main_effects_continuous_
    for(size_t j = 0; j < q_; ++j) {
        main_effects_continuous_(j) = params(idx++);
    }

    // 4. pairwise_effects_cross_ row-major (included edges only)
    for(size_t i = 0; i < p_; ++i) {
        for(size_t j = 0; j < q_; ++j) {
            if(gxy(i, j) == 1) {
                pairwise_effects_cross_(i, j) = params(idx++);
            }
        }
    }

    // 5. Kyy theta block: forward map (f_q, psi_q) -> (Phi, K) with
    //    excluded-edge zeros enforced by the null-space parameterization
    size_t chol_dim = chol_constraint_structure_.active_dim;
    arma::vec theta_yy = params.subvec(idx, idx + chol_dim - 1);
    const ForwardMapResult& fm = yy_engine_.forward_map(theta_yy);
    cholesky_of_precision_ = fm.Phi;
    pairwise_effects_continuous_ = -0.5 * fm.K;
    bool ok = arma::solve(inv_cholesky_of_precision_,
                          arma::trimatu(cholesky_of_precision_),
                          arma::eye(q_, q_), arma::solve_opts::fast);
    if(!ok) {
        // Fallback: recompute from scratch
        recompute_pairwise_effects_continuous_decomposition();
    } else {
        covariance_continuous_ = inv_cholesky_of_precision_ *
                                 inv_cholesky_of_precision_.t();
        log_det_precision_ = 2.0 * arma::accu(fm.psi);
    }
    theta_yy_ = std::move(theta_yy);
    theta_yy_valid_ = true;

    // Refresh caches
    recompute_conditional_mean();
    recompute_marginal_interactions();
}

arma::vec MixedMRFModel::get_active_inv_mass() const {
    if(constraint_dirty_) {
        const_cast<MixedMRFModel*>(this)->ensure_constraint_structure();
    }

    arma::vec active(parameter_dimension());
    // Main effects: always active
    active.head(num_main_) = inv_mass_.head(num_main_);

    size_t offset_full = num_main_;
    size_t offset_active = num_main_;

    // pairwise_effects_discrete_ included edges
    for(size_t i = 0; i < p_ - 1; ++i) {
        for(size_t j = i + 1; j < p_; ++j) {
            if(gxx(i, j) == 1) {
                active(offset_active++) = inv_mass_(offset_full);
            }
            offset_full++;
        }
    }

    // main_effects_continuous_: always active
    for(size_t j = 0; j < q_; ++j) {
        active(offset_active++) = inv_mass_(offset_full++);
    }

    // pairwise_effects_cross_ included edges
    for(size_t i = 0; i < p_; ++i) {
        for(size_t j = 0; j < q_; ++j) {
            if(gxy(i, j) == 1) {
                active(offset_active++) = inv_mass_(offset_full);
            }
            offset_full++;
        }
    }

    // Kyy theta block: each (f_q, psi_q) entry carries the mass of its
    // full Cholesky slot (excluded slots are skipped by the traversal)
    const auto& cs = chol_constraint_structure_;
    cs.for_each_active_full_pair([&](size_t active_idx, size_t full_idx) {
        active(offset_active + active_idx) = inv_mass_(offset_full + full_idx);
    });

    return active;
}

arma::ivec MixedMRFModel::get_vectorized_indicator_parameters() {
    size_t total = num_pairwise_xx_ + num_pairwise_yy_ + num_cross_;
    arma::ivec out(total);
    size_t idx = 0;

    // 1. Upper-triangle of Gxx
    for(size_t i = 0; i < p_ - 1; ++i) {
        for(size_t j = i + 1; j < p_; ++j) {
            out(idx++) = gxx(i, j);
        }
    }

    // 2. Upper-triangle of Gyy
    for(size_t i = 0; i < q_ - 1; ++i) {
        for(size_t j = i + 1; j < q_; ++j) {
            out(idx++) = gyy(i, j);
        }
    }

    // 3. Full Gxy block row-major
    for(size_t i = 0; i < p_; ++i) {
        for(size_t j = 0; j < q_; ++j) {
            out(idx++) = gxy(i, j);
        }
    }

    return out;
}


// =============================================================================
// Infrastructure
// =============================================================================

void MixedMRFModel::set_seed(int seed) {
    rng_ = SafeRNG(seed);
}

std::unique_ptr<BaseModel> MixedMRFModel::clone() const {
    return std::make_unique<MixedMRFModel>(*this);
}


// =============================================================================
// Missing data imputation
// =============================================================================

void MixedMRFModel::set_missing_data(const arma::imat& missing_discrete,
                                      const arma::imat& missing_continuous) {
    missing_index_discrete_ = missing_discrete;
    missing_index_continuous_ = missing_continuous;
    has_missing_ = (missing_index_discrete_.n_rows > 0 ||
                    missing_index_continuous_.n_rows > 0);
}

void MixedMRFModel::impute_missing() {
    if(!has_missing_) return;

    // --- Phase 1: Impute discrete entries ---
    const int num_disc_missing = missing_index_discrete_.n_rows;
    if(num_disc_missing > 0) {
        arma::vec category_probabilities(max_cats_ + 1);

        for(int miss = 0; miss < num_disc_missing; miss++) {
            const int person = missing_index_discrete_(miss, 0);
            const int variable = missing_index_discrete_(miss, 1);
            const int num_cats = num_categories_(variable);

            // Rest score: 2 * sum_t x_vt A_xx(t,s) + 2 sum_j y_vj A_xy(s,j)
            // A_xx diagonal is zero, so no self-interaction subtraction needed
            double rest_v = 0.0;
            for(size_t t = 0; t < p_; t++) {
                rest_v += 2.0 * discrete_observations_dbl_(person, t) * pairwise_effects_discrete_(t, variable);
            }
            for(size_t j = 0; j < q_; j++) {
                rest_v += 2.0 * continuous_observations_(person, j) * pairwise_effects_cross_(variable, j);
            }

            double cumsum = 0.0;

            // Max-shift the category exponents so exp() cannot overflow; the
            // shift cancels in the normalized inverse-transform draw below.
            if(is_ordinal_variable_(variable)) {
                // P(x=0) ∝ exp(0), P(x=c) ∝ exp(c · rest + μ_x(s, c-1))
                double max_exp = 0.0;
                for(int c = 1; c <= num_cats; c++) {
                    double e = static_cast<double>(c) * rest_v +
                               main_effects_discrete_(variable, c - 1);
                    if(e > max_exp) max_exp = e;
                }
                cumsum = MY_EXP(-max_exp);
                category_probabilities(0) = cumsum;
                for(int c = 1; c <= num_cats; c++) {
                    double exponent = static_cast<double>(c) * rest_v +
                                      main_effects_discrete_(variable, c - 1) - max_exp;
                    cumsum += MY_EXP(exponent);
                    category_probabilities(c) = cumsum;
                }
            } else {
                // Blume-Capel: categories centered at baseline
                const int ref = baseline_category_(variable);
                double alpha = main_effects_discrete_(variable, 0);
                double beta = main_effects_discrete_(variable, 1);
                double max_exp = -arma::datum::inf;
                for(int cat = 0; cat <= num_cats; cat++) {
                    const int score = cat - ref;
                    double e = alpha * score + beta * score * score + score * rest_v;
                    if(e > max_exp) max_exp = e;
                }
                cumsum = 0.0;
                for(int cat = 0; cat <= num_cats; cat++) {
                    const int score = cat - ref;
                    double exponent = alpha * score +
                                      beta * score * score +
                                      score * rest_v - max_exp;
                    cumsum += MY_EXP(exponent);
                    category_probabilities(cat) = cumsum;
                }
            }

            // Sample via inverse-transform
            double u = runif(rng_) * cumsum;
            int sampled = 0;
            while(u > category_probabilities(sampled)) {
                sampled++;
            }

            int new_value = sampled;
            if(!is_ordinal_variable_(variable)) {
                new_value -= baseline_category_(variable);
            }
            const int old_value = discrete_observations_(person, variable);

            if(new_value != old_value) {
                discrete_observations_(person, variable) = new_value;
                discrete_observations_dbl_(person, variable) =
                    static_cast<double>(new_value);

                if(is_ordinal_variable_(variable)) {
                    counts_per_category_(old_value, variable)--;
                    counts_per_category_(new_value, variable)++;
                } else {
                    blume_capel_stats_(0, variable) += (new_value - old_value);
                    blume_capel_stats_(1, variable) +=
                        (new_value * new_value - old_value * old_value);
                }
            }
        }
    }

    // --- Phase 2: Refresh caches that depend on the discrete data ---
    if(num_disc_missing > 0) {
        // The gradient reads discrete_observations_dbl_t_ directly, not via
        // ensure_gradient_cache; hold it equal to the transpose of the
        // discrete data.
        discrete_observations_dbl_t_ = discrete_observations_dbl_.t();
        // conditional_mean_ is a function of the discrete data and is read by
        // every MH acceptance ratio (log_conditional_ggm).
        recompute_conditional_mean();
    }

    // --- Phase 3: Impute continuous entries ---
    const int num_cont_missing = missing_index_continuous_.n_rows;
    if(num_cont_missing > 0) {
        for(int miss = 0; miss < num_cont_missing; miss++) {
            const int person = missing_index_continuous_(miss, 0);
            const int variable = missing_index_continuous_(miss, 1);

            // Conditional: y_vj | y_{v,-j}, x ~ N(mu*, 1/precision_jj)
            // mu* = M_vj - sum_{k!=j} (interaction_jk / interaction_jj) * (y_vk - M_vk)
            double precision_jj = -2.0 * pairwise_effects_continuous_(variable, variable);
            double cond_mean = conditional_mean_(person, variable);
            for(size_t k = 0; k < q_; k++) {
                if(k != static_cast<size_t>(variable)) {
                    cond_mean -= (pairwise_effects_continuous_(variable, k) / pairwise_effects_continuous_(variable, variable)) *
                        (continuous_observations_(person, k) -
                         conditional_mean_(person, k));
                }
            }
            double cond_sd = std::sqrt(1.0 / precision_jj);

            continuous_observations_(person, variable) =
                rnorm(rng_, cond_mean, cond_sd);
        }
    }

    // Invalidate gradient cache (observations changed)
    invalidate_gradient_cache();
}


// =============================================================================
// Stubs (to be implemented in later phases)
// =============================================================================

void MixedMRFModel::do_one_metropolis_step(int iteration) {
    recompute_am_caches();

    // Per-slot accept-probability and visit-mask matrices for the five
    // proposal-SD storages. Only entries we actually visit get mask=1; the
    // adapter only RM-updates those slots.
    arma::mat  ar_main_disc  = arma::zeros<arma::mat >(p_, max_cats_);
    arma::umat mask_main_disc= arma::zeros<arma::umat>(p_, max_cats_);
    arma::mat  ar_main_cont  = arma::zeros<arma::mat >(q_, 1);
    arma::umat mask_main_cont= arma::zeros<arma::umat>(q_, 1);
    arma::mat  ar_pair_disc  = arma::zeros<arma::mat >(p_, p_);
    arma::umat mask_pair_disc= arma::zeros<arma::umat>(p_, p_);
    arma::mat  ar_pair_cont  = arma::zeros<arma::mat >(q_, q_);
    arma::umat mask_pair_cont= arma::zeros<arma::umat>(q_, q_);
    arma::mat  ar_pair_cross = arma::zeros<arma::mat >(p_, q_);
    arma::umat mask_pair_cross= arma::zeros<arma::umat>(p_, q_);

    // Step 1: main effects (ordinal thresholds or BC α/β)
    for(size_t s = 0; s < p_; ++s) {
        if(is_ordinal_variable_(s)) {
            for(int c = 0; c < num_categories_(s); ++c) {
                ar_main_disc(s, c) = update_main_effect(s, c, std::nullopt);
                mask_main_disc(s, c) = 1;
            }
        } else {
            ar_main_disc(s, 0) = update_main_effect(s, 0, std::nullopt);
            ar_main_disc(s, 1) = update_main_effect(s, 1, std::nullopt);
            mask_main_disc(s, 0) = 1;
            mask_main_disc(s, 1) = 1;
        }
    }

    // Step 2: continuous means
    for(size_t j = 0; j < q_; ++j) {
        ar_main_cont(j, 0) = update_continuous_mean(j, std::nullopt);
        mask_main_cont(j, 0) = 1;
    }

    // Step 3: pairwise_effects_discrete_ (upper triangle, edge-gated)
    for(size_t i = 0; i < p_ - 1; ++i)
        for(size_t j = i + 1; j < p_; ++j)
            if(!edge_selection_active_ || gxx(i, j) == 1) {
                ar_pair_disc(i, j) = update_pairwise_discrete(i, j, std::nullopt);
                mask_pair_disc(i, j) = 1;
            }

    // Step 4: pairwise_effects_continuous_ (off-diag + diagonal, edge-gated)
    if(q_ >= 2) {
        for(size_t i = 0; i < q_ - 1; ++i)
            for(size_t j = i + 1; j < q_; ++j)
                if(!edge_selection_active_ || gyy(i, j) == 1) {
                    ar_pair_cont(i, j) = update_pairwise_effects_continuous_offdiag(i, j, std::nullopt);
                    mask_pair_cont(i, j) = 1;
                }
    }
    for(size_t i = 0; i < q_; ++i) {
        ar_pair_cont(i, i) = update_pairwise_effects_continuous_diag(i, std::nullopt);
        mask_pair_cont(i, i) = 1;
    }

    // Step 5: pairwise_effects_cross_ (edge-gated)
    for(size_t i = 0; i < p_; ++i)
        for(size_t j = 0; j < q_; ++j)
            if(!edge_selection_active_ || gxy(i, j) == 1) {
                ar_pair_cross(i, j) = update_pairwise_cross(i, j, std::nullopt);
                mask_pair_cross(i, j) = 1;
            }

    // Robbins-Monro batch update on each storage's adapter (MH mode only).
    if (mh_adapter_main_discrete_)
        mh_adapter_main_discrete_->update(mask_main_disc, ar_main_disc, iteration);
    if (mh_adapter_main_continuous_)
        mh_adapter_main_continuous_->update(mask_main_cont, ar_main_cont, iteration);
    if (mh_adapter_pairwise_discrete_)
        mh_adapter_pairwise_discrete_->update(mask_pair_disc, ar_pair_disc, iteration);
    if (mh_adapter_pairwise_continuous_)
        mh_adapter_pairwise_continuous_->update(mask_pair_cont, ar_pair_cont, iteration);
    if (mh_adapter_pairwise_cross_)
        mh_adapter_pairwise_cross_->update(mask_pair_cross, ar_pair_cross, iteration);
}

void MixedMRFModel::init_metropolis_adaptation(const WarmupSchedule& schedule) {
    mh_adapter_main_discrete_ = std::make_unique<MetropolisAdaptationController>(
        proposal_sd_main_discrete_, schedule, target_accept_);
    mh_adapter_main_continuous_ = std::make_unique<MetropolisAdaptationController>(
        proposal_sd_main_continuous_, schedule, target_accept_);
    mh_adapter_pairwise_discrete_ = std::make_unique<MetropolisAdaptationController>(
        proposal_sd_pairwise_discrete_, schedule, target_accept_);
    mh_adapter_pairwise_continuous_ = std::make_unique<MetropolisAdaptationController>(
        proposal_sd_pairwise_continuous_, schedule, target_accept_);
    mh_adapter_pairwise_cross_ = std::make_unique<MetropolisAdaptationController>(
        proposal_sd_pairwise_cross_, schedule, target_accept_);
}

void MixedMRFModel::sweep_within_model_mh(std::optional<double> rm_weight) {
    recompute_am_caches();

    // Step 1: Update all main effects (ordinal thresholds or BC α/β)
    for(size_t s = 0; s < p_; ++s) {
        if(is_ordinal_variable_(s)) {
            for(int c = 0; c < num_categories_(s); ++c)
                update_main_effect(s, c, rm_weight);
        } else {
            update_main_effect(s, 0, rm_weight);  // linear α
            update_main_effect(s, 1, rm_weight);  // quadratic β
        }
    }

    // Step 2: Update all continuous means
    for(size_t j = 0; j < q_; ++j)
        update_continuous_mean(j, rm_weight);

    // Step 3: Update pairwise_effects_discrete_ (upper triangle, edge-gated)
    for(size_t i = 0; i < p_ - 1; ++i)
        for(size_t j = i + 1; j < p_; ++j)
            if(!edge_selection_active_ || gxx(i, j) == 1)
                update_pairwise_discrete(i, j, rm_weight);

    // Step 4: Update pairwise_effects_continuous_ (off-diag + diagonal, edge-gated)
    if(q_ >= 2) {
        for(size_t i = 0; i < q_ - 1; ++i)
            for(size_t j = i + 1; j < q_; ++j)
                if(!edge_selection_active_ || gyy(i, j) == 1)
                    update_pairwise_effects_continuous_offdiag(i, j, rm_weight);
    }
    for(size_t i = 0; i < q_; ++i)
        update_pairwise_effects_continuous_diag(i, rm_weight);

    // Step 5: Update pairwise_effects_cross_ (edge-gated)
    for(size_t i = 0; i < p_; ++i)
        for(size_t j = 0; j < q_; ++j)
            if(!edge_selection_active_ || gxy(i, j) == 1)
                update_pairwise_cross(i, j, rm_weight);

    // Edge-indicator updates are handled by ChainRunner, not here.
    // (Matches the OMRF pattern; avoids double-counting indicator proposals.)
}

void MixedMRFModel::update_edge_indicators() {
    if(!edge_selection_active_) return;

    invalidate_gradient_cache();
    recompute_am_caches();

    // Discrete-discrete edges (shuffled order)
    for(size_t e = 0; e < num_pairwise_xx_; ++e) {
        size_t idx = edge_order_xx_(e);
        update_edge_indicator_discrete(edge_pairs_xx_(idx, 0), edge_pairs_xx_(idx, 1));
    }

    // Continuous-continuous edges (shuffled order)
    for(size_t e = 0; e < num_pairwise_yy_; ++e) {
        size_t idx = edge_order_yy_(e);
        update_edge_indicator_continuous(edge_pairs_yy_(idx, 0), edge_pairs_yy_(idx, 1));
    }

    // Cross edges (shuffled order)
    for(size_t e = 0; e < num_cross_; ++e) {
        size_t idx = edge_order_xy_(e);
        size_t i = idx / q_;
        size_t j = idx % q_;
        update_edge_indicator_cross(i, j);
    }
}

void MixedMRFModel::prepare_iteration() {
    // Shuffle edge-update order to avoid order bias.
    // Always called, even when edge selection is off, to keep RNG consistent.
    edge_order_xx_ = arma_randperm(rng_, num_pairwise_xx_);
    edge_order_yy_ = arma_randperm(rng_, num_pairwise_yy_);
    edge_order_xy_ = arma_randperm(rng_, num_cross_);
}

void MixedMRFModel::tune_proposal_sd(int iteration, const WarmupSchedule& schedule) {
    auto rm_weight_opt = schedule.rm_weight_for_proposal_sd(iteration);
    if (!rm_weight_opt) return;
    // Stage-3b sweep: re-run every within-model MH proposal with RM
    // applied to its proposal-SD slot via *rm_weight_opt. Sampler-agnostic.
    sweep_within_model_mh(rm_weight_opt);
}

#pragma once

#include <array>
#include <memory>
#include "models/base_model.h"
#include "math/cholesky_helpers.h"
#include "rng/rng_utils.h"
#include "models/ggm/graph_constraint_structure.h"
#include "models/ggm/ggm_gradient.h"
#include "models/ggm/zratio_engine.h"
#include "models/ggm/zratio_gauge.h"
#include "priors/parameter_prior.h"
#include "mcmc/samplers/metropolis_adaptation.h"


/**
 * GGMModel - Gaussian Graphical Model
 *
 * Bayesian inference on the precision matrix (inverse covariance) of a
 * multivariate Gaussian via element-wise Metropolis-Hastings. Edge
 * selection uses a spike-and-slab prior with Cauchy slab.
 *
 * The Cholesky factor of the precision matrix is maintained incrementally
 * through rank-1 updates/downdates after each element change.
 *
 * Inherits from BaseModel for compatibility with the generic MCMC framework.
 */
class GGMModel : public BaseModel {
public:

    /**
     * Construct from raw observations.
     *
     * Computes the sufficient-statistic matrix S = X'X from the raw data.
     * When na_impute is true, the observation matrix is retained for
     * full-conditional imputation of missing entries.
     *
     * @param observations          Raw data matrix (n x p)
     * @param inclusion_probability Prior inclusion probabilities for each edge
     * @param initial_edge_indicators Initial edge inclusion indicators
     * @param edge_selection        Enable edge selection (spike-and-slab)
     * @param interaction_prior      Polymorphic prior on the off-diagonal interaction/pairwise parameters
     * @param diagonal_prior         Polymorphic prior on the diagonal precision parameters
     * @param na_impute             Retain observations for missing-data imputation
     */
    GGMModel(
            const arma::mat& observations,
            const arma::mat& inclusion_probability,
            const arma::imat& initial_edge_indicators,
            const bool edge_selection,
            std::unique_ptr<BaseParameterPrior> interaction_prior,
            std::unique_ptr<BaseParameterPrior> diagonal_prior,
            const bool na_impute = false
    )   // Delegate to the sufficient-statistics constructor: the raw-data path
        // differs only in deriving n (centered data has n-1 effective df) and
        // S = X'X from the observations, and in retaining the observations for
        // missing-data imputation. All member init and the MLE warm-start live
        // in the delegated constructor.
        : GGMModel(observations.n_rows - 1,
                   observations.t() * observations,
                   inclusion_probability,
                   initial_edge_indicators,
                   edge_selection,
                   std::move(interaction_prior),
                   std::move(diagonal_prior))
    {
        if (na_impute) {
            observations_ = observations;
        }
    }

    /**
     * Construct from sufficient statistics.
     *
     * Bypasses raw data storage; useful when only X'X and n are available.
     * Missing-data imputation is not supported with this constructor.
     *
     * @param n                     Number of observations
     * @param suf_stat              Sufficient-statistic matrix X'X (p x p)
     * @param inclusion_probability Prior inclusion probabilities for each edge
     * @param initial_edge_indicators Initial edge inclusion indicators
     * @param edge_selection        Enable edge selection (spike-and-slab)
     * @param interaction_prior      Polymorphic prior on the off-diagonal interaction/pairwise parameters
     * @param diagonal_prior         Polymorphic prior on the diagonal precision parameters
     */
    GGMModel(
            const int n,
            const arma::mat& suf_stat,
            const arma::mat& inclusion_probability,
            const arma::imat& initial_edge_indicators,
            const bool edge_selection,
            std::unique_ptr<BaseParameterPrior> interaction_prior,
            std::unique_ptr<BaseParameterPrior> diagonal_prior
    ) : n_(n),
        p_(suf_stat.n_cols),
        dim_((p_ * (p_ + 1)) / 2),
        suf_stat_(suf_stat),
        inclusion_probability_(inclusion_probability),
        edge_selection_(edge_selection),
        interaction_prior_(std::move(interaction_prior)),
        diagonal_prior_(std::move(diagonal_prior)),
        precision_matrix_(arma::eye<arma::mat>(p_, p_)),
        cholesky_of_precision_(arma::eye<arma::mat>(p_, p_)),
        inv_cholesky_of_precision_(arma::eye<arma::mat>(p_, p_)),
        covariance_matrix_(arma::eye<arma::mat>(p_, p_)),
        omega_(arma::ones<arma::mat>(p_, p_)),
        edge_indicators_(initial_edge_indicators),
        vectorized_indicator_parameters_(edge_selection_ ? dim_ : 0),
        rb_alpha_(edge_selection_ ? dim_ : 0, arma::fill::zeros),
        rb_pregamma_(edge_selection_ ? dim_ : 0),
        proposal_sds_(arma::mat(dim_, 1, arma::fill::ones) * 0.25),
        num_pairwise_(p_ * (p_ - 1) / 2),
        precision_proposal_(arma::mat(p_, p_, arma::fill::none))
    {
        edge_pairs_.set_size(num_pairwise_, 2);
        size_t flat = 0;
        for (size_t i = 0; i + 1 < p_; ++i) {
            for (size_t j = i + 1; j < p_; ++j, ++flat) {
                edge_pairs_(flat, 0) = i;
                edge_pairs_(flat, 1) = j;
            }
        }
        initialize_precision_from_mle();
        // Pre-move state sentinel: -1 marks an edge not yet proposed (and the
        // never-proposed diagonal slots), which the odds accumulator skips.
        if (rb_pregamma_.n_elem > 0) rb_pregamma_.fill(-1);
    }

    /** Copy constructor for cloning (required for parallel chains). */
    GGMModel(const GGMModel& other)
        : BaseModel(other),
          target_accept_(other.target_accept_),
          determinant_tilt_(other.determinant_tilt_),
          n_(other.n_),
          p_(other.p_),
          dim_(other.dim_),
          suf_stat_(other.suf_stat_),
          inclusion_probability_(other.inclusion_probability_),
          edge_selection_(other.edge_selection_),
          interaction_prior_(other.interaction_prior_->clone()),
          diagonal_prior_(other.diagonal_prior_->clone()),
          precision_matrix_(other.precision_matrix_),
          cholesky_of_precision_(other.cholesky_of_precision_),
          inv_cholesky_of_precision_(other.inv_cholesky_of_precision_),
          covariance_matrix_(other.covariance_matrix_),
          log_det_precision_(other.log_det_precision_),
          omega_(other.omega_),
          edge_indicators_(other.edge_indicators_),
          vectorized_indicator_parameters_(other.vectorized_indicator_parameters_),
          rb_alpha_(other.rb_alpha_),
          rb_pregamma_(other.rb_pregamma_),
          proposal_sds_(other.proposal_sds_),
          shuffled_edge_order_(other.shuffled_edge_order_),
          edge_pairs_(other.edge_pairs_),
          num_pairwise_(other.num_pairwise_),
          rng_(other.rng_),
          observations_(other.observations_),
          has_missing_(other.has_missing_),
          missing_index_(other.missing_index_),
          precision_proposal_(other.precision_proposal_),
          constraint_structure_(other.constraint_structure_),
          gradient_engine_(other.gradient_engine_),
          constraint_dirty_(other.constraint_dirty_),
          theta_valid_(other.theta_valid_),
          theta_(other.theta_),
          zratio_engine_(other.zratio_engine_
                             ? std::make_shared<ZRatioEngine>(*other.zratio_engine_)
                             : nullptr)
    {
        // The engine's oracle draws from this clone's chain RNG.
        if (zratio_engine_) zratio_engine_->set_rng(&rng_);
    }

    /**
     * Attach the per-edge Z-ratio engine, switching the between-edge moves
     * to the hierarchical prior specification p(K | Gamma) = rho/Z(Gamma):
     * the add acceptance gains log J = log(Z(Gamma-)/Z(Gamma+)) and the
     * delete acceptance its negation. Each chain clone deep-copies the
     * engine, so per-chain caches never cross threads.
     */
    void set_zratio_engine(std::shared_ptr<ZRatioEngine> engine) {
        zratio_engine_ = std::move(engine);
        if (zratio_engine_) zratio_engine_->set_rng(&rng_);
    }

    /** Trust gauge available iff the hierarchical Z-ratio engine is attached. */
    bool gauge_available() const override { return zratio_engine_ != nullptr; }

    void set_gauge_active(bool on, int n_draws, int cap) override {
        if (on) {
            zratio_gauge_.reset();
            zratio_gauge_.nE = static_cast<double>(get_num_pairwise());
            zratio_gauge_.n_draws = n_draws;
            zratio_gauge_.cap = cap;
        }
        zratio_gauge_.active = on;
    }
    void set_zratio_phase(ZRatioPhase phase) override {
        if (zratio_engine_) zratio_engine_->set_phase(phase);
    }
    void gauge_begin_sweep() override { zratio_gauge_.begin_sweep(); }
    void gauge_end_sweep() override { zratio_gauge_.end_sweep(); }

    /**
     * Copy the Z-ratio engine's end-of-run state (cache/hit counters and the
     * trust-gauge block) into the chain result. No-op without an engine.
     */
    void collect_chain_diagnostics(ChainResult& chain_result) const override;

    /** @return true when edge selection is enabled. */
    bool has_edge_selection()  const override { return edge_selection_; }
    /** @return true when missing-data imputation is active. */
    bool has_missing_data()    const override { return has_missing_; }

    /** Impute missing entries from full-conditional normal distributions. */
    void impute_missing() override;

    /**
     * Register missing-data locations.
     *
     * @param missing_index  M x 2 matrix of 0-based (row, col) indices
     * @throws std::logic_error if the model was constructed without na_impute
     */
    void set_missing_data(const arma::imat& missing_index) {
        if (observations_.n_elem == 0) {
            throw std::logic_error(
                "set_missing_data() called but observations_ is empty. "
                "The model must be constructed with na_impute=true to retain observations.");
        }
        missing_index_ = missing_index;
        has_missing_ = (missing_index.n_rows > 0 && missing_index.n_cols == 2);
    }

    /**
     * Set the Robbins-Monro target acceptance rate used by the
     * adaptive-Metropolis updates of this GGM. Honoured by all
     * Metropolis sweeps (off-diagonal and diagonal).
     */
    void set_metropolis_target_accept(double target) override {
        target_accept_ = target;
    }

    /**
     * Construct Robbins-Monro adaptation controller for the per-iteration
     * MH proposal SDs. Called once by MetropolisSampler before warmup;
     * under NUTS this is never called and the controller stays null.
     */
    void init_metropolis_adaptation(const WarmupSchedule& schedule) override;

    /**
     * @return true iff the row-block Gibbs within-step supports the current
     * priors: a Normal or Cauchy slab on the K_yy off-diagonals and a Gamma
     * prior on K_ii/2. The Gamma shape (alpha != 1) and the determinant tilt
     * (delta != 0) do NOT gate eligibility -- they are handled by the
     * independent-MH correction and the xi Gamma-shape shift, respectively.
     */
    bool row_block_gibbs_eligible();

    /**
     * Set the determinant-tilt exponent delta. Adds delta * log|K| to the
     * log-prior, pushing the chain away from the PD-cone boundary. delta = 0
     * (default) recovers the untilted target. Currently consumed only by
     * the NUTS gradient engine; the MH path is unchanged. Triggers an engine
     * rebuild on the next gradient call.
     */
    void set_determinant_tilt(double delta) {
        determinant_tilt_ = delta;
        constraint_dirty_ = true;
    }

    /** Shuffle edge visit order (random scan). */
    void prepare_iteration() override;

    /** Sweep over edges in shuffled order, proposing add/remove moves. */
    void update_edge_indicators() override;

    /**
     * Element-wise MH updates for proposal-SD tuning during stage 3b.
     *
     * Runs off-diagonal and diagonal Metropolis updates with
     * Robbins-Monro adaptation, following the OMRF pattern.
     */
    void tune_proposal_sd(int iteration, const WarmupSchedule& schedule) override;

    /**
     * Combined log-posterior and gradient for NUTS.
     *
     * Uses the free-element Cholesky parameterization:
     * theta = (psi_1, f_2, psi_2, ..., f_p, psi_p) where psi_q = log(phi_qq)
     * and x_q = N_q f_q gives the off-diagonal Cholesky entries.
     *
     * @param parameters  Active theta vector (dimension = p + |E|)
     * @return (log-posterior, gradient) pair
     */
    std::pair<double, arma::vec> logp_and_gradient(
        const arma::vec& parameters) override;

    /**
     * Set model state from a theta vector (inverse of get_vectorized_parameters).
     *
     * Runs the forward map theta -> Phi -> K and updates all internal
     * matrices (precision, Cholesky, inverse Cholesky, covariance).
     *
     * @param parameters  Active theta vector (dimension = p + |E|)
     */
    void set_vectorized_parameters(const arma::vec& parameters) override;

    /**
     * Perform one full Metropolis sweep.
     *
     * Iterates over all off-diagonal entries (edge parameter updates) and
     * all diagonal entries; edge-indicator add-delete moves are handled
     * separately in update_edge_indicators().
     *
     * @param iteration  Current iteration index (for Robbins-Monro adaptation)
     */
    void do_one_metropolis_step(int iteration = -1) override;

    /**
     * Perform one full row-block Gibbs sweep.
     *
     * Iterates over rows i = 0..p-1, drawing each column of K from its exact
     * conjugate full-conditional via update_row_block_gibbs(i). No proposal
     * adaptation (the draw is exact); edge-indicator add-delete moves are
     * handled separately in update_edge_indicators().
     *
     * Precondition: row_block_gibbs_eligible() == true.
     *
     * @param iteration  Current iteration index (unused; kept for interface
     *                   symmetry with do_one_metropolis_step).
     */
    void do_one_gibbs_step(int iteration = -1) override;

    /** Enable the full-conditional edge birth/death proposal (Gibbs path). */
    void set_conjugate_edge_proposal(bool enable) override {
        use_conjugate_edge_proposal_ = enable;
    }

    /**
     * @return Active theta dimension: p + |E| (diagonals + included edges).
     *
     * Changes when edge indicators toggle. Used by NUTS for leapfrog
     * integration.
     */
    size_t parameter_dimension() const override;

    /**
     * @return Full theta dimension: p + p(p-1)/2 (all possible off-diag slots).
     *
     * Fixed across all graphs. Used by the adaptation controller for
     * mass-matrix sizing.
     */
    size_t full_parameter_dimension() const override;

    /**
     * @return Storage dimension for sample output: p(p+1)/2 (upper triangle of K).
     *
     * Preserves the existing output contract: downstream R code expects
     * the upper triangle of the precision matrix.
     */
    size_t storage_dimension() const override { return dim_; }

    /**
     * Set random seed for reproducibility.
     * @param seed  Integer seed value
     */
    void set_seed(int seed) override {
        rng_ = SafeRNG(seed);
    }

    /**
     * @return Active theta vector: (psi_1, f_2, psi_2, ..., f_p, psi_p).
     *
     * Dimension = parameter_dimension() = p + |E|. Used by NUTS as the
     * current state. Recomputed lazily from Phi when stale.
     */
    arma::vec get_vectorized_parameters() const override;

    /**
     * @return Full (zero-padded) theta vector for mass-matrix adaptation.
     *
     * Dimension = full_parameter_dimension() = p + p(p-1)/2. Inactive
     * edges have their f_q slots set to zero.
     */
    arma::vec get_full_vectorized_parameters() const override;

    /**
     * @return Upper triangle of the precision matrix for sample storage.
     *
     * Preserves the existing output contract.
     */
    arma::vec get_storage_vectorized_parameters() const override {
        return extract_upper_triangle();
    }

    /** @return Upper triangle of the edge-indicator matrix as an integer vector. */
    arma::ivec get_vectorized_indicator_parameters() override {
        size_t e = 0;
        for (size_t i = 0; i < p_; ++i) {
            for (size_t j = i; j < p_; ++j) {
                vectorized_indicator_parameters_(e) = edge_indicators_(i, j);
                ++e;
            }
        }
        return vectorized_indicator_parameters_;
    }

    /**
     * Get per-edge Rao-Blackwellized inclusion draws from the most recent
     * update_edge_indicators() sweep, ordered to match
     * get_vectorized_indicator_parameters().
     */
    arma::vec get_vectorized_rb_inclusion() override;

    /** Per-edge acceptance probability from the last sweep (raw alpha). */
    arma::vec get_vectorized_rb_alpha() override;

    /** Per-edge pre-move indicator state from the last sweep (0/1; -1 = none). */
    arma::ivec get_vectorized_rb_pregamma() override;

    /** @return Reference to the model's random number generator. */
    SafeRNG& get_rng() override { return rng_; }

    /** @return Current edge-indicator matrix. */
    const arma::imat& get_edge_indicators() const override {
        return edge_indicators_;
    }

    /** @return Mutable reference to the prior inclusion-probability matrix. */
    arma::mat& get_inclusion_probability() override {
        return inclusion_probability_;
    }

    /** @return const reference to the current precision matrix K. */
    const arma::mat& get_precision_matrix() const {
        return precision_matrix_;
    }

    /** @return Number of variables (p). */
    int get_num_variables() const override {
        return static_cast<int>(p_);
    }

    /** @return Number of unique off-diagonal pairs p(p-1)/2. */
    int get_num_pairwise() const override {
        return static_cast<int>(p_ * (p_ - 1) / 2);
    }

    /**
     * @return Active subset of the inverse mass diagonal.
     *
     * Filters the full inv_mass_ (dimension p + p(p-1)/2) to active
     * parameters only (dimension p + |E|). For columns where N_q != I,
     * rotates the per-Cholesky-entry variances into f_q coordinates.
     */
    arma::vec get_active_inv_mass() const override;

    /** @return Deep copy of this model. */
    std::unique_ptr<BaseModel> clone() const override {
        return std::make_unique<GGMModel>(*this);
    }

private:

    // Robbins-Monro target acceptance rate for adaptive-Metropolis
    // proposal-SD tuning. Set via set_metropolis_target_accept(); defaults
    // to 0.44 (componentwise random-walk Metropolis optimum).
    double target_accept_ = 0.44;

    /// Per-iteration adaptation controller (MH mode only — under NUTS this
    /// stays null and the stage-3b path in tune_proposal_sd is used instead).
    std::unique_ptr<MetropolisAdaptationController> metropolis_adapter_;

    // Determinant-tilt exponent (see set_determinant_tilt). Forwarded to
    // GGMGradientEngine on every rebuild.
    double determinant_tilt_ = 0.0;

    /** Extract upper triangle of the precision matrix into a vector. */
    arma::vec extract_upper_triangle() const {
        arma::vec result(dim_);
        size_t e = 0;
        for (size_t i = 0; i < p_; ++i) {
            for (size_t j = i; j < p_; ++j) {
                result(e) = precision_matrix_(i, j);
                ++e;
            }
        }
        return result;
    }

    /// Number of observations.
    size_t n_;
    /// Number of variables.
    size_t p_;
    /// Number of upper-triangle elements: p(p+1)/2.
    size_t dim_;
    /// Sufficient-statistic matrix X'X (p x p).
    arma::mat suf_stat_;
    /// Prior inclusion probabilities (p x p, symmetric).
    arma::mat inclusion_probability_;
    /// Whether the model was constructed with edge selection.
    bool edge_selection_;
    /// Use the full-conditional edge birth/death proposal (Gibbs sampler) in
    /// place of the random-walk Roverato proposal. Set by the GibbsSampler.
    bool use_conjugate_edge_proposal_ = false;
    /// Prior on off-diagonal precision elements (interactions).
    std::unique_ptr<BaseParameterPrior> interaction_prior_;
    /// Prior on diagonal precision elements (scale).
    std::unique_ptr<BaseParameterPrior> diagonal_prior_;

    /// Precision matrix Omega, its Cholesky factor R (Omega = R'R),
    /// inverse Cholesky factor, and covariance matrix.
    arma::mat precision_matrix_, cholesky_of_precision_, inv_cholesky_of_precision_, covariance_matrix_;
    double log_det_precision_ = 0.0;    ///< Cached log|K|, refreshed with the Cholesky caches
    /// Per-edge scale-mixture weight for the Cauchy slab: K_yy_ij | omega_ij ~
    /// N(0, sigma^2 omega_ij), omega_ij ~ InvGamma(1/2, 1/2). Fixed at 1 for a
    /// Normal slab (the mixture collapses to the plain Normal). Used only by
    /// the row-block Gibbs within-step.
    arma::mat omega_;
    /// Current edge-indicator matrix (p x p, symmetric, 0/1).
    arma::imat edge_indicators_;
    /// Pre-allocated storage returned by get_vectorized_indicator_parameters().
    arma::ivec vectorized_indicator_parameters_;
    /// Per-edge acceptance probability (raw alpha) and pre-move indicator
    /// state from the last update_edge_indicators() sweep. Same indexing
    /// (row-major upper triangle of size dim_) as
    /// vectorized_indicator_parameters_; diagonal entries are never proposed
    /// (pregamma stays -1). The RB draw J and the odds accumulators derive
    /// from these two.
    arma::vec rb_alpha_;
    arma::ivec rb_pregamma_;

    /// Proposal standard deviations for Metropolis updates (one per element,
    /// stored as a (dim_, 1) matrix so it can be wrapped by
    /// MetropolisAdaptationController).
    arma::mat proposal_sds_;

    /// Shuffled edge visit order for random-scan edge selection.
    arma::uvec shuffled_edge_order_;
    arma::umat edge_pairs_;             ///< num_pairwise x 2 flat-index -> (i, j) table (row-major upper triangle)
    /// Number of unique off-diagonal pairs: p(p-1)/2.
    size_t num_pairwise_ = 0;
    /// Random number generator.
    SafeRNG rng_;

    /// Raw observation matrix (n x p), only populated when na_impute=true.
    arma::mat observations_;
    /// Whether missing-data imputation is active.
    bool has_missing_ = false;
    /// M x 2 matrix of 0-based (row, col) indices of missing entries.
    arma::imat missing_index_;

    /**
     * Incrementally adjust S = X'X after replacing one observation value.
     *
     * @param variable  Column index of the changed variable
     * @param person    Row index of the changed observation
     * @param delta     Change in value (new - old)
     */
    void update_suf_stat_for_imputation(int variable, int person, double delta);

    /// Scratch matrix for proposed precision values.
    arma::mat precision_proposal_;

    /**
     * Workspace for conditional precision reparameterization.
     *
     * - [0] Phi_q1q
     * - [1] Phi_q1q1
     * - [2] omega_ij - Phi_q1q * Phi_q1q1
     * - [3] Phi_q1q1
     * - [4] omega_jj - Phi_q1q^2
     * - [5] constrained diagonal at x = 0
     */
    std::array<double, 6> constants_{};

    /**
     * Work vectors for rank-2 Cholesky update.
     *
     * A symmetric rank-2 update  A + vf1*vf2' + vf2*vf1'  is decomposed
     * into two rank-1 updates via  u1 = (vf1+vf2)/sqrt(2),
     * u2 = (vf1-vf2)/sqrt(2).
     */
    arma::vec v1_ = {0, -1};
    arma::vec v2_ = {0, 0};
    arma::vec vf1_ = arma::zeros<arma::vec>(p_);
    arma::vec vf2_ = arma::zeros<arma::vec>(p_);
    arma::vec u1_ = arma::zeros<arma::vec>(p_);
    arma::vec u2_ = arma::zeros<arma::vec>(p_);

    /**
     * Row-block Gibbs per-row scratch, reused across the p rows of a sweep so
     * the within-step does not churn the allocator. These are the buffers
     * filled element-by-element (loops) or via set_size, which Armadillo
     * resizes by reusing the existing allocation when it fits -- unlike the
     * solve/chol/matmul results, which materialise a temporary and are left
     * as locals. gibbs_Ni_ in particular replaces a per-row std::vector with
     * a reserve(p-1) that allocated ~p on every call.
     */
    std::vector<arma::uword> gibbs_Ni_;
    arma::uvec gibbs_support_;
    arma::vec gibbs_beta_old_;
    arma::vec gibbs_sigma_iNi_;
    arma::vec gibbs_s_Ni_i_;
    arma::mat gibbs_C_;

    /**
     * Propose a new off-diagonal precision entry via a normal perturbation
     * on an unconstrained reparameterization. Accepts or rejects with a
     * Metropolis ratio using the Gaussian likelihood and Cauchy prior.
     *
     * @param i  Row index (i < j)
     * @param j  Column index
     * @return   Metropolis acceptance probability min(1, exp(ln_alpha)),
     *           or 0.0 if the edge is inactive (caller masks it out).
     */
    double update_edge_parameter(size_t i, size_t j);

    /**
     * Closed-form Gibbs draw for row i of K given the rest of K and the graph.
     * The per-row primitive driven by do_one_gibbs_step().
     *
     * Origin: Wang's row-by-row block Gibbs update for Gaussian graphical
     * models (each column of the precision matrix sampled from its exact
     * full-conditional).
     *
     * Partitions K with A = K_{-i,-i}, beta = K_{N_i, i}, kii = K_{i,i}, where
     * N_i is the active neighbour set of i. With Normal slab N(0, sigma^2) on
     * K_yy_ij = -K_ij/2 and Gamma(alpha = 1, beta0) prior on K_ii/2, the
     * conditional (beta, xi = kii - beta^T C beta) is conjugate:
     *
     *   xi   | rest ~ Gamma(n/2 + delta + 1, (beta0 + S_ii)/2)
     *   beta | rest ~ N(-M^{-1} S_{N_i, i}, M^{-1}),
     *     M = (beta0 + S_ii) C + diag(1/(4 sigma^2 omega_k)),
     *     C = (A^{-1})_{N_i, N_i} = Sigma_{N_i, N_i} - Sigma_{N_i, i} Sigma_{i, N_i}/Sigma_ii.
     *
     * The determinant tilt delta enters as a shift in the xi Gamma shape. A
     * Gamma shape alpha != 1 on the diagonal is handled by an independent-MH
     * accept on (kii_new/kii_old)^(alpha-1) around the alpha = 1 proposal. A
     * Cauchy slab enters through the scale-mixture weights omega_k (fixed at 1
     * for a Normal slab), refreshed per sweep by do_one_gibbs_step().
     *
     * Precondition: row_block_gibbs_eligible() == true; covariance_matrix_
     * holds K^{-1} up to date with precision_matrix_.
     *
     * @param i  Row index.
     */
    void update_row_block_gibbs(size_t i);

    /** @return the slab scale sigma on K_yy (Normal or Cauchy interaction prior). */
    double slab_scale_() const;

    /** @return true when the interaction (slab) prior is Cauchy. */
    bool slab_is_cauchy_() const;

    /**
     * Refresh the Cauchy scale-mixture weights omega_ from the current K.
     * With all of K held fixed, each active edge weight is a closed-form draw
     *   omega_ij | K ~ InvGamma(1, 1/2 + K_yy_ij^2 / (2 sigma^2)),
     * K_yy_ij = -K_ij/2. No-op for a Normal slab.
     */
    void refresh_cauchy_omega_();

    /**
     * Propose a new diagonal precision entry on the log scale.
     * Accepts or rejects with a Metropolis ratio using the Gaussian
     * likelihood, a Gamma(1,1) prior, and a Jacobian correction.
     *
     * @param i  Diagonal index
     * @return   Metropolis acceptance probability min(1, exp(ln_alpha)).
     */
    double update_diagonal_parameter(size_t i);

    /**
     * Within-model MH move primitives. Each proposes, computes the full
     * acceptance ratio, and applies the move (with Cholesky update) on accept,
     * returning the RAW ln_alpha. Both the sampling path
     * (update_edge_parameter/update_diagonal_parameter) and the Robbins-Monro
     * tuning path (tune_proposal_sd) call these, so the proposal RNG and math
     * are defined in exactly one place. Callers own the edge-active guard, the
     * cache-flag invalidation, and the min(1,exp())/RM consumption of ln_alpha.
     *
     * ggm_edge_move assumes the edge (i,j) is active (caller short-circuits).
     */
    double ggm_edge_move(size_t i, size_t j);
    double ggm_diag_move(size_t i);

    /**
     * Positive-definiteness canary for prior-only chains (n == 0).
     *
     * With data, the likelihood term (n/2) * log|K| vetoes proposals that
     * leave the PD cone. Without data there is no such anchor: the
     * reparameterization constants are computed from the incrementally
     * maintained covariance, whose floating-point drift can place a proposal
     * outside the cone with finite acceptance probability. An accepted
     * non-PD state invalidates the Cholesky machinery and the next
     * refresh_cholesky() throws. ggm_edge_move, ggm_diag_move and
     * update_edge_indicator_conjugate reject such proposals explicitly when
     * n == 0.
     *
     * @return true if precision_proposal_ admits a Cholesky factorization
     */
    bool proposal_is_positive_definite_() const;

    /**
     * Positive-definiteness of an edge proposal given its three changed
     * entries, for callers that have not filled precision_proposal_.
     *
     * The conjugate birth/death move (update_edge_indicator_conjugate) forms
     * its proposal as scalars and never writes the full matrix on the data
     * path, so it fills precision_proposal_ here before deferring to
     * proposal_is_positive_definite_. Only the prior-only path calls this, so
     * the O(p^2) copy and O(p^3) factorization are off the data path.
     *
     * @param i, j       Edge indices (i < j).
     * @param kij, kjj   Proposed K(i,j) and K(j,j).
     * @return true if the proposed precision matrix admits a Cholesky
     *         factorization
     */
    bool edge_proposal_is_positive_definite_(size_t i, size_t j, double kij,
                                             double kjj);

    /**
     * Metropolis-Hastings add-delete move for an edge indicator.
     *
     * If the edge is on, proposes deletion; if off, proposes a new value
     * from a scaled normal. Acceptance combines the likelihood ratio,
     * Bernoulli prior odds, Cauchy slab, and proposal density.
     *
     * @param i  Row index (i < j)
     * @param j  Column index
     * @return Rao-Blackwellized inclusion draw J = gamma + (1 - 2 gamma) alpha
     */
    double update_edge_indicator_parameter_pair(size_t i, size_t j);

    /**
     * Full-conditional edge birth/death for the joint spec (Normal slab,
     * alpha = 1). The cofactor move preserves |K|, so F(phi) is Gaussian and
     * the proposal is the exact conditional of the toggled coordinate; the
     * acceptance reduces to the inclusion odds times p_slab(0)/q(0),
     * independent of the proposed value. No proposal-SD tuning. Used by the
     * Gibbs sampler in place of update_edge_indicator_parameter_pair.
     * @return Rao-Blackwellized inclusion draw J = gamma + (1 - 2 gamma) alpha
     */
    double update_edge_indicator_conjugate(size_t i, size_t j);

    /**
     * Precompute reparameterization constants for the (i, j) element.
     *
     * Derives six values from the cofactor structure of the inverse
     * precision matrix that allow off-diagonal proposals on an
     * unconstrained scale while deterministically satisfying the
     * positive-definiteness constraint on the diagonal.
     *
     * @param i  Row index
     * @param j  Column index
     */
    void get_constants(size_t i, size_t j);



    /**
     * Return the diagonal value omega_jj required to keep the precision
     * matrix positive definite after changing the off-diagonal element to x.
     *
     * @param x  Proposed off-diagonal value omega_ij
     * @return   Constrained diagonal value omega_jj
     */
    double constrained_diagonal(const double x) const;

    /**
     * Log-likelihood ratio for a proposed off-diagonal element change,
     * computed via the matrix-determinant lemma (rank-2 update).
     *
     * @param i  Row index of the changed element
     * @param j  Column index of the changed element
     */
    double log_density_impl_edge(size_t i, size_t j) const;

    /**
     * Log-likelihood ratio for a proposed diagonal element change,
     * computed via the matrix-determinant lemma (rank-1 update).
     *
     * @param j  Index of the changed diagonal element
     */
    double log_density_impl_diag(size_t j) const;

    /**
     * log|K_prop| - log|K_curr| for a rank-2 off-diagonal proposal at (i, j),
     * computed via the matrix-determinant lemma in O(p). Reads
     * precision_matrix_, precision_proposal_, and covariance_matrix_; assumes
     * precision_proposal_ has already been set up at (i, j), (j, i), (j, j).
     * Used to add the determinant-tilt term delta * (log|K_prop| - log|K_curr|)
     * to MH ratios.
     */
    double log_det_ratio_edge(size_t i, size_t j) const;

    /**
     * log|K_prop| - log|K_curr| for a rank-1 diagonal proposal at j.
     * Computed via the matrix-determinant lemma in O(1). Reads the same
     * cached state as log_det_ratio_edge.
     */
    double log_det_ratio_diag(size_t j) const;



    /**
     * Refresh Sigma and log|K| after changing an off-diagonal element.
     *
     * Advances log_det_precision_ by the rank-2 determinant-lemma ratio,
     * encodes the accepted (i,j)/(j,j) change as the vf1/vf2 rank-2 pair,
     * and delegates the covariance update to apply_rank2_chol_smw_update_.
     * chol(K) is not maintained per accept (rebuilt once per sweep).
     *
     * @param omega_ij_old  Previous value of omega(i,j)
     * @param omega_jj_old  Previous value of omega(j,j)
     * @param i             Row index
     * @param j             Column index
     */
    void cholesky_update_after_edge(double omega_ij_old, double omega_jj_old, size_t i, size_t j);

    /**
     * Apply a symmetric rank-2 update to Sigma for a K update already written.
     *
     * Given vf1_, vf2_ of length p, this refreshes
     *   Sigma <- Sherman-Morrison-Woodbury rank-2 update (O(p^2)) for
     *   K_new  = K_old + vf1 vf2^T + vf2 vf1^T,
     * with a fallback to refresh_cholesky() when the 2x2 capacitance is
     * near-singular.
     *
     * Inputs are taken from the model's vf1_, vf2_ scratch members so callers
     * can populate them in-place without an extra copy. The helper does not
     * touch precision_matrix_ -- the caller must already have written the
     * post-update entries it represents. Generic in vf1, vf2: the edge update
     * passes sparse 2-entry vectors; the row-block Gibbs sweep reuses it with
     * full-vector inputs.
     *
     * `support` lists the nonzero indices of vf1/vf2 ({i,j} for an edge
     * accept, {i} + N_i for a row-Gibbs row); the SMW matvec Sigma * u
     * touches only those columns, O(p |support|) instead of a dense O(p^2)
     * gemv, with a dense fallback when the support is near-full.
     *
     * chol(K) and log|K| are NOT maintained here. The accept paths advance
     * log_det_precision_ by the determinant-lemma ratio before calling this
     * helper, nothing reads chol(K) between accepts within a sweep, and every
     * sweep (Metropolis, tuning, edge-indicator, row-Gibbs) ends with
     * refresh_cholesky(), which rebuilds the factor and resets Sigma and the
     * log-det exactly -- bounding the SMW floating-point drift per sweep.
     */
    void apply_rank2_chol_smw_update_(const arma::uvec& support);

    /**
     * Refresh Sigma and log|K| after changing a diagonal element.
     *
     * Advances log_det_precision_ by the rank-1 determinant-lemma ratio and
     * applies a Sherman-Morrison rank-1 update (O(p^2)) to the covariance
     * matrix, with a fallback to refresh_cholesky() when the scalar
     * capacitance is near-singular. chol(K) is not maintained per accept
     * (see apply_rank2_chol_smw_update_).
     *
     * @param omega_ii_old  Previous value of omega(i,i)
     * @param i             Diagonal index
     */
    void cholesky_update_after_diag(double omega_ii_old, size_t i);

    /**
     * Check (Sigma K)(r, r) = 1 on the given rows.
     *
     * Guards the SMW-maintained covariance against catastrophic cancellation
     * when K passes near a singular state (legitimate under the prior at
     * delta = 0, where a row-Gibbs xi draw can be arbitrarily small): the
     * update that moves K away from the near-singular state subtracts huge
     * outer products and Sigma loses absolute accuracy. Callers refresh from
     * K on violation before the drifted Sigma can corrupt proposal constants
     * or the row-Gibbs Schur extraction (which would write a
     * non-positive-definite K).
     *
     * @param rows  Row indices touched by the update just applied.
     * @return false when any |(Sigma K)(r, r) - 1| exceeds kSigmaProbeTol_
     *         or is non-finite.
     */
    bool sigma_rows_consistent_(const arma::uvec& rows) const;

    /** Tolerance on |(Sigma K)(r, r) - 1| in the per-accept probe. */
    static constexpr double kSigmaProbeTol_ = 1e-6;

    /**
     * Recompute chol(K), its inverse, Sigma, and log|K| from the precision
     * matrix. Called once at the end of every sweep (the accept paths
     * maintain only Sigma and the log-det incrementally) and when a
     * per-accept probe or SMW capacitance guard detects a drifted Sigma.
     */
    void refresh_cholesky();

    /**
     * Initialize precision matrix at the regularized MLE.
     *
     * Computes K = n * inv(S + delta * I) where delta provides
     * Ledoit-Wolf-style shrinkage toward identity. Gives NUTS a
     * starting point near the posterior mode, avoiding the step-size
     * instability that arises when starting from K = I far from the
     * mode.
     */
    void initialize_precision_from_mle();

    // =================================================================
    // NUTS gradient support
    // =================================================================

    /// Graph constraint structure (rebuilt when edge indicators change).
    GraphConstraintStructure constraint_structure_;
    /// Gradient engine for the free-element Cholesky parameterization.
    GGMGradientEngine gradient_engine_;
    /// Whether the constraint structure needs rebuilding.
    bool constraint_dirty_ = true;
    /// Whether theta_ is in sync with cholesky_of_precision_.
    mutable bool theta_valid_ = false;
    /// Cached theta vector (active parameterization).
    mutable arma::vec theta_;
    /// Per-edge Z-ratio engine (hierarchical prior spec); null on the
    /// joint spec. Deep-copied per chain clone (owns a mutable cache).
    std::shared_ptr<ZRatioEngine> zratio_engine_;
    /// In-chain trust-gauge accumulator; active only during assessment sweeps.
    ZRatioGauge zratio_gauge_;

public:
    /**
     * Mark the NUTS gradient caches stale after a change to the precision
     * matrix or the edge set. Both the constraint structure (active dimension
     * depends on the edge indicators) and the cached theta parameterization
     * become invalid, so the MH within-K, edge-selection, and proposal-SD
     * tuning paths all funnel their invalidation through here rather than
     * repeating the two flag writes. (set_determinant_tilt is intentionally
     * NOT routed through this: it only needs a gradient-engine rebuild, not a
     * theta invalidation.)
     */
    void invalidate_gradient_cache() {
        constraint_dirty_ = true;
        theta_valid_ = false;
    }

    /**
     * Rebuild the constraint structure and gradient engine from current
     * edge indicators. Called lazily before gradient evaluation.
     */
    void ensure_constraint_structure();

    /**
     * Convert the current Cholesky factor to the theta parameterization.
     *
     * For each column q, computes psi_q = log(phi_qq) and
     * f_q = N_q^T x_q where x_q = Phi[0:q-1, q].
     */
    void recompute_theta() const;
};

/**
 * Construct a GGMModel from an R list.
 *
 * Dispatches to the sufficient-statistics constructor (when the list
 * contains `n` and `suf_stat`) or the raw-data constructor (when the
 * list contains `X`).
 *
 * @param inputFromR              R list with data (either `X` or `n` + `suf_stat`)
 * @param inclusion_probability   Prior inclusion probabilities for each edge
 * @param initial_edge_indicators Initial edge inclusion indicators
 * @param edge_selection          Enable edge selection (spike-and-slab)
 * @param interaction_prior       Prior on pairwise interaction / off-diagonal precision parameters
 * @param diagonal_prior          Prior on diagonal precision parameters
 * @param na_impute               Retain observations for missing-data imputation
 * @return Fully constructed GGMModel
 */
GGMModel createGGMModelFromR(
    const Rcpp::List& inputFromR,
    const arma::mat& inclusion_probability,
    const arma::imat& initial_edge_indicators,
    const bool edge_selection,
    std::unique_ptr<BaseParameterPrior> interaction_prior,
    std::unique_ptr<BaseParameterPrior> diagonal_prior,
    const bool na_impute = false
);

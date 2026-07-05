#pragma once

#include <array>
#include <memory>
#include <optional>
#include "models/base_model.h"
#include "models/ggm/graph_constraint_structure.h"
#include "models/ggm/ggm_gradient.h"
#include "models/ggm/zratio_engine.h"
#include "math/cholesky_helpers.h"
#include "math/cholupdate.h"
#include "rng/rng_utils.h"
#include "priors/parameter_prior.h"
#include "mcmc/samplers/metropolis_adaptation.h"
#include "utils/variable_helpers.h"

/**
 * MixedMRFModel - Mixed Markov Random Field Model
 *
 * Joint model for p discrete (ordinal or Blume-Capel) variables x and
 * q continuous variables y.  The joint density is:
 *
 *   log f(x, y) ∝ Σ_s μ_{x,s}(x_s) + x' A_xx x + y' A_yy y + 2 x' A_xy y
 *
 * All three interaction blocks (A_xx, A_yy, A_xy) enter the density at
 * the same scale, so a Cauchy(0, scale) prior has the same meaning for
 * every block.
 *
 * A_yy is stored as a negative semi-definite matrix with negative diagonal.
 * The positive-definite precision matrix is Precision = -2 A_yy.
 * Internally the Cholesky decomposition and covariance
 * cache operate on Precision.
 *
 * Uses the marginal pseudo-likelihood, with and without edge selection
 * via spike-and-slab priors.
 *
 * Discrete variables are either ordinal (free category thresholds, category
 * 0 as reference) or Blume-Capel (linear α + quadratic β, user-specified
 * reference).  Blume-Capel observations are centered at their baseline
 * category in the constructor, matching the OMRFModel convention.
 *
 * Inherits from BaseModel for compatibility with the generic MCMC framework
 * (ChainRunner, MetropolisSampler, WarmupSchedule).
 */
class MixedMRFModel : public BaseModel {
public:

    // =========================================================================
    // Construction
    // =========================================================================

    /**
     * Construct from raw observations.
     *
     * @param discrete_observations   Integer matrix of discrete observations (n × p, 0-based)
     * @param continuous_observations  Continuous observations (n × q)
     * @param num_categories       Number of categories per discrete variable (p-vector)
     * @param is_ordinal_variable  1 = ordinal, 0 = Blume-Capel (p-vector)
     * @param baseline_category    Reference category per discrete variable (p-vector)
     * @param inclusion_probability Prior inclusion probabilities ((p+q) × (p+q))
     * @param initial_edge_indicators Initial edge inclusion matrix ((p+q) × (p+q))
     * @param edge_selection       Enable edge selection (spike-and-slab)
     * @param interaction_prior     Polymorphic prior on pairwise interactions
     * @param threshold_prior      Polymorphic prior on main effects / thresholds
     * @param means_prior          Polymorphic prior on continuous means
     * @param diagonal_prior       Polymorphic prior on precision diagonal
     * @param seed                 RNG seed for reproducibility
     */
    MixedMRFModel(
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
    );

    /** Copy constructor for cloning (required for parallel chains). */
    MixedMRFModel(const MixedMRFModel& other);

    // =========================================================================
    // Capability queries
    // =========================================================================

    /** @return true when edge selection is enabled. */
    bool has_edge_selection() const override { return edge_selection_; }
    /** @return true when missing-data imputation is active. */
    bool has_missing_data() const override { return has_missing_; }

    // =========================================================================
    // Core sampling methods
    // =========================================================================

    /**
     * Combined log pseudo-posterior and gradient evaluation.
     * @param parameters  NUTS-dimension parameter vector
     * @return Pair of (log-pseudo-posterior, gradient)
     */
    std::pair<double, arma::vec> logp_and_gradient(
        const arma::vec& parameters) override;

    /**
     * Perform one full Metropolis sweep over all parameter groups.
     * @param iteration  Current iteration (for Robbins-Monro adaptation)
     */
    void do_one_metropolis_step(int iteration = -1) override;

    /**
     * Set the Robbins-Monro target acceptance rate used by the
     * adaptive-Metropolis updates of this mixed model. Honoured by the
     * proposal-SD tuning of every Metropolis component (discrete main,
     * continuous main, pairwise discrete, pairwise continuous, pairwise
     * cross).
     */
    void set_metropolis_target_accept(double target) override {
        target_accept_ = target;
    }

    /**
     * Set the determinant-tilt exponent delta for the Kyy block. Adds
     * delta * log|Kyy| to the log-prior, pushing the continuous-block
     * precision matrix away from the PD-cone boundary. delta = 0
     * (default) recovers the untilted target. Currently consumed only
     * by the NUTS gradient paths; the MH path is unchanged.
     */
    void set_determinant_tilt_yy(double delta) {
        determinant_tilt_yy_ = delta;
    }

    /**
     * Attach the per-edge Z-ratio engine, switching the continuous-block
     * between-edge moves to the hierarchical prior specification
     * p(K_yy | Gamma_yy) = rho/Z(Gamma_yy): the add acceptance gains
     * log J = log(Z(Gamma-)/Z(Gamma+)) and the delete acceptance its
     * negation, with the mediating-block counts read off the continuous
     * subgraph. Each chain clone deep-copies the engine.
     */
    void set_zratio_engine(std::shared_ptr<ZRatioEngine> engine) {
        zratio_engine_ = std::move(engine);
        if (zratio_engine_) zratio_engine_->set_rng(&rng_);
    }

    /** Freeze the Z-ratio calibrator at the warmup/sampling boundary. */
    void on_warmup_end() override {
        if (zratio_engine_) zratio_engine_->freeze_calibration();
    }

    /**
     * Copy the Z-ratio engine's end-of-run state (counters, frozen
     * constant block, calibration anchors) into the chain result. No-op
     * without an engine.
     */
    void collect_chain_diagnostics(ChainResult& chain_result) const override;

    /**
     * Construct Robbins-Monro adaptation controllers for the per-iteration
     * MH proposal SDs. Called once by MetropolisSampler before warmup; under
     * NUTS this is never called and the controllers stay null. One adapter
     * per proposal-SD storage (5 in total).
     */
    void init_metropolis_adaptation(const WarmupSchedule& schedule) override;

    /**
     * Tune proposal SDs via Robbins-Monro (Stage 3b).
     *
     * Re-runs every MH sweep with the schedule-supplied RM weight applied to
     * each proposal-SD slot. Outside stage 3b the schedule returns nullopt
     * and this is a no-op. Mirrors `OMRFModel::tune_proposal_sd` /
     * `GGMModel::tune_proposal_sd`; sampler-agnostic by construction.
     */
    void tune_proposal_sd(int iteration, const WarmupSchedule& schedule) override;

    /**
     * Shuffle edge update order at the start of each iteration.
     * Advances the RNG state consistently even when edge selection is off.
     */
    void prepare_iteration() override;

    // =========================================================================
    // Edge selection
    // =========================================================================

    /** Perform one sweep of Metropolis-Hastings edge add-delete moves. */
    void update_edge_indicators() override;

    /**
     * Enable or disable edge-selection proposals.
     * @param active  true to enable edge add-delete moves
     */
    void set_edge_selection_active(bool active) override {
        edge_selection_active_ = active;
    }

    // =========================================================================
    // Parameter vectorization
    // =========================================================================

    /**
     * Dimensionality of the active NUTS parameter space.
     * Includes Cholesky entries for the continuous precision (always full q(q+1)/2).
     * When edge selection is active, excludes discrete/cross parameters for inactive edges.
     */
    size_t parameter_dimension() const override;

    /**
     * Full NUTS dimension (all params regardless of edge state).
     * Includes q(q+1)/2 Cholesky entries. Used for mass-matrix sizing.
     */
    size_t full_parameter_dimension() const override;

    /**
     * Storage dimension (all parameters including continuous precision
     * as A_yy entries, regardless of edge state). Used for sample storage.
     */
    size_t storage_dimension() const override;

    /** Get active NUTS parameters as a flat vector (includes Cholesky block). */
    arma::vec get_vectorized_parameters() const override;

    /** Get all NUTS parameters (inactive edges zeroed, includes Cholesky block). */
    arma::vec get_full_vectorized_parameters() const override;

    /** Get all parameters as A_yy entries for sample storage. */
    arma::vec get_storage_vectorized_parameters() const override;

    /** Set NUTS parameters from a flat vector (includes Cholesky block). */
    void set_vectorized_parameters(const arma::vec& params) override;

    /** Get vectorized edge indicators (Gxx upper-tri, Gyy upper-tri, Gxy full). */
    arma::ivec get_vectorized_indicator_parameters() override;

    /** Get active subset of inverse mass diagonal (includes Cholesky block). */
    arma::vec get_active_inv_mass() const override;

    // =========================================================================
    // Infrastructure
    // =========================================================================

    /** Set random seed for reproducibility. */
    void set_seed(int seed) override;

    /** Clone the model for parallel execution. */
    std::unique_ptr<BaseModel> clone() const override;

    /** @return Reference to the model's random number generator. */
    SafeRNG& get_rng() override { return rng_; }

    /** @return Current edge-indicator matrix ((p+q) × (p+q)). */
    const arma::imat& get_edge_indicators() const override {
        return edge_indicators_;
    }

    /** @return Mutable reference to the prior inclusion-probability matrix. */
    arma::mat& get_inclusion_probability() override {
        return inclusion_probability_;
    }

    /** @return Total number of variables (p + q). */
    int get_num_variables() const override {
        return static_cast<int>(p_ + q_);
    }

    /**
     * Number of unique off-diagonal pairs in the (p+q) × (p+q) indicator
     * matrix: p(p-1)/2 + q(q-1)/2 + p*q.
     */
    int get_num_pairwise() const override {
        return static_cast<int>(num_pairwise_xx_ + num_pairwise_yy_ + num_cross_);
    }

    // =========================================================================
    // Missing data
    // =========================================================================

    /** Impute missing entries from full-conditional distributions. */
    void impute_missing() override;

    /**
     * Register missing-data locations for discrete and continuous sub-matrices.
     *
     * @param missing_discrete  M_d x 2 matrix of 0-based (row, col) indices into discrete_observations_
     * @param missing_continuous M_c x 2 matrix of 0-based (row, col) indices into continuous_observations_
     */
    void set_missing_data(const arma::imat& missing_discrete,
                          const arma::imat& missing_continuous);

private:

    // Robbins-Monro target acceptance rate for adaptive-Metropolis
    // proposal-SD tuning. Set via set_metropolis_target_accept(); defaults
    // to 0.44 (componentwise random-walk Metropolis optimum).
    double target_accept_ = 0.44;

    /// Per-edge Z-ratio engine for the hierarchical spec on the continuous
    /// block (null under the joint spec). Deep-copied per chain clone.
    std::shared_ptr<ZRatioEngine> zratio_engine_;

    // Determinant-tilt exponent on the Kyy block (see set_determinant_tilt_yy).
    // Adds determinant_tilt_yy_ * log|Kyy| to the NUTS log-prior; MH ratios
    // are not yet adjusted.
    double determinant_tilt_yy_ = 0.0;

    /// Per-iteration adaptation controllers (MH mode only — under NUTS these
    /// stay null and the stage-3b path in tune_proposal_sd is used instead).
    /// One adapter per proposal-SD storage; off-diag and diag of the continuous
    /// pairwise share `proposal_sd_pairwise_continuous_` and thus one adapter.
    std::unique_ptr<MetropolisAdaptationController> mh_adapter_main_discrete_;
    std::unique_ptr<MetropolisAdaptationController> mh_adapter_main_continuous_;
    std::unique_ptr<MetropolisAdaptationController> mh_adapter_pairwise_discrete_;
    std::unique_ptr<MetropolisAdaptationController> mh_adapter_pairwise_continuous_;
    std::unique_ptr<MetropolisAdaptationController> mh_adapter_pairwise_cross_;

    // =========================================================================
    // Counts and dimensions
    // =========================================================================

    size_t n_;                          ///< Number of observations
    size_t p_;                          ///< Number of discrete variables
    size_t q_;                          ///< Number of continuous variables
    size_t num_main_;                   ///< Total main-effect params (sum C_s for ord + 2 per BC)
    size_t num_pairwise_xx_;            ///< p(p-1)/2
    size_t num_pairwise_yy_;            ///< q(q-1)/2
    size_t num_cross_;                  ///< p * q
    size_t num_cholesky_ = 0;           ///< q(q+1)/2 — number of Cholesky entries

    // =========================================================================
    // Data
    // =========================================================================

    arma::imat discrete_observations_;   ///< Discrete observations (n x p), BC columns centered
    arma::mat discrete_observations_dbl_; ///< Double version (post-centering)
    arma::mat continuous_observations_;  ///< Continuous observations (n x q)
    arma::ivec num_categories_;         ///< Categories per discrete variable (p-vector)
    int max_cats_;                      ///< max(num_categories)
    arma::uvec is_ordinal_variable_;    ///< 1 = ordinal, 0 = Blume-Capel (p-vector)
    arma::ivec baseline_category_;      ///< Reference category per discrete variable (p-vector)

    // =========================================================================
    // Missing data
    // =========================================================================

    arma::imat missing_index_discrete_;   ///< M_d x 2 (row, col) for missing discrete entries
    arma::imat missing_index_continuous_; ///< M_c x 2 (row, col) for missing continuous entries
    bool has_missing_ = false;            ///< Whether imputation is active

    // =========================================================================
    // Sufficient statistics
    // =========================================================================

    arma::imat counts_per_category_;    ///< (max_cats+1) x p category counts (ordinal only)
    arma::imat blume_capel_stats_;      ///< 2 x p linear/quadratic sums (BC only)

    // =========================================================================
    // Parameters
    // =========================================================================

    arma::mat main_effects_discrete_;                     ///< p x max_cats main effects (thresholds or alpha/beta)
    arma::vec main_effects_continuous_;                     ///< q-vector continuous means
    arma::mat pairwise_effects_discrete_;                     ///< p x p discrete interactions (symmetric, zero diag)
    arma::mat pairwise_effects_continuous_;                     ///< q x q continuous interaction matrix (negative-definite)
    arma::mat pairwise_effects_cross_;                     ///< p x q cross-type interactions

    // =========================================================================
    // Edge indicators
    // =========================================================================

    /// Combined (p+q) x (p+q) indicator matrix.
    /// Gxx block: rows [0,p), cols [0,p) -- symmetric, zero diag.
    /// Gyy block: rows [p,p+q), cols [p,p+q) -- symmetric, zero diag.
    /// Gxy block: rows [0,p), cols [p,p+q) -- full p x q rectangle.
    arma::imat edge_indicators_;
    arma::mat inclusion_probability_;   ///< Prior inclusion probabilities
    bool edge_selection_;               ///< Enable edge selection
    bool edge_selection_active_;        ///< Currently in edge selection phase

    // =========================================================================
    // Priors
    // =========================================================================

    std::unique_ptr<BaseParameterPrior> interaction_prior_;   ///< Prior on pairwise interactions
    std::unique_ptr<BaseParameterPrior> threshold_prior_;    ///< Prior on main effects / thresholds
    std::unique_ptr<BaseParameterPrior> means_prior_;        ///< Prior on continuous means
    std::unique_ptr<BaseParameterPrior> diagonal_prior_;     ///< Prior on precision diagonal

    // =========================================================================
    // Proposal SDs (Robbins-Monro adapted)
    // =========================================================================

    arma::mat proposal_sd_main_discrete_;             ///< p x max_cats
    arma::mat proposal_sd_main_continuous_;             ///< q x 1 (mat-shaped for MetropolisAdaptationController)
    arma::mat proposal_sd_pairwise_discrete_;             ///< p x p
    arma::mat proposal_sd_pairwise_continuous_;             ///< q x q
    arma::mat proposal_sd_pairwise_cross_;             ///< p x q

    // =========================================================================
    // Cached quantities
    // =========================================================================

    arma::mat cholesky_of_precision_;       ///< q x q upper Cholesky R (Precision = R'R)
    arma::mat inv_cholesky_of_precision_;   ///< q x q R^{-1} (upper triangular)
    arma::mat covariance_continuous_;       ///< q x q Σ = Precision^{-1}
    double log_det_precision_;              ///< log|Precision|
    arma::mat marginal_interactions_;                       ///< p x p marginal PL interaction matrix
    arma::mat cross_term_;                  ///< p x p cached 2 A_xy Σ A_xy' (marginal PL cross term)
    arma::mat conditional_mean_;            ///< n x q conditional mean

    // Rank-1 Cholesky update workspace
    std::array<double, 6> cont_constants_{};  ///< Reparameterization constants
    arma::mat precision_proposal_;        ///< q x q scratch for proposed precision
    // 2-element scratch vectors for the symmetric rank-2 precision update;
    // entries are set per edge before use (mirrors GGMModel::v1_/v2_).
    arma::vec cont_v1_ = {0, -1};
    arma::vec cont_v2_ = {0, 0};
    arma::vec cont_vf1_;                      ///< q-vector, zeroed between uses
    arma::vec cont_vf2_;                      ///< q-vector, zeroed between uses
    arma::vec cont_u1_;                       ///< q-vector workspace
    arma::vec cont_u2_;                       ///< q-vector workspace

    // =========================================================================
    // Gradient cache (populated by ensure_gradient_cache)
    // =========================================================================

    arma::mat discrete_observations_dbl_t_; ///< p x n transpose (BLAS gradient)
    arma::vec grad_obs_cache_;          ///< Cached observed-data gradient component
    arma::imat disc_index_cache_;        ///< p x p map from (i,j) to gradient index
    arma::imat cross_index_cache_;        ///< p x q map from (i,j) to gradient index
    int main_effects_continuous_grad_offset_ = 0;           ///< Offset of main_effects_continuous block in gradient vector
    int chol_grad_offset_ = 0;          ///< Offset of Cholesky block in gradient vector
    bool gradient_cache_valid_ = false; ///< Whether gradient cache is current

    // Per-chain scratch for compute_logZ_and_probs_*_into. Reused across
    // every call to logp_and_gradient and across every variable inside it.
    mutable LogZAndProbs logz_out_;
    mutable LogZScratch  logz_scratch_;

    // =========================================================================
    // RATTLE constraint structure
    // =========================================================================

    /// Cholesky constraint structure (per-column excluded/included for Gyy block).
    GraphConstraintStructure chol_constraint_structure_;
    /// Kyy-block theta-space engine (forward map + reverse-Givens adjoint).
    GGMGradientEngine yy_engine_;
    /// Cached Kyy theta block (f_q, psi_q per column), lazily recomputed.
    mutable arma::vec theta_yy_;
    /// Whether theta_yy_ matches the current cholesky_of_precision_ and graph.
    mutable bool theta_yy_valid_ = false;
    /// Offset of Cholesky block (Block 5) in the full-space vector.
    size_t chol_block_offset_ = 0;
    /// Whether constraint structure needs rebuilding.
    bool constraint_dirty_ = true;
    /// Whether initial graph is sparse (constraints without edge selection).
    bool has_sparse_graph_ = false;

    // =========================================================================
    // RNG and edge-update order
    // =========================================================================

    SafeRNG rng_;                       ///< Per-chain random number generator
    arma::uvec edge_order_xx_;          ///< Shuffled xx-edge pair indices
    arma::uvec edge_order_yy_;          ///< Shuffled yy-edge pair indices
    arma::uvec edge_order_xy_;          ///< Shuffled xy-edge pair indices
    arma::umat edge_pairs_xx_;          ///< num_pairwise_xx x 2 flat-index -> (i, j) table
    arma::umat edge_pairs_yy_;          ///< num_pairwise_yy x 2 flat-index -> (i, j) table

    // =========================================================================
    // Private helpers
    // =========================================================================

    /** Count total main-effect parameters across all discrete variables. */
    size_t count_num_main_effects() const;

    /** Compute category counts and BC sufficient statistics from discrete_observations_. */
    void compute_sufficient_statistics();

    /** Recompute conditional_mean_ from main_effects_continuous_, pairwise_effects_cross_, covariance_continuous_. */
    void recompute_conditional_mean();

    /** Recompute cholesky_of_precision_, inv_cholesky_of_precision_, covariance_continuous_, log_det_precision_ from pairwise_effects_continuous_. */
    void recompute_pairwise_effects_continuous_decomposition();

    /** Recompute marginal_interactions_ from pairwise_effects_discrete_, pairwise_effects_cross_, covariance_continuous_ (marginal PL only). Refreshes cross_term_. */
    void recompute_marginal_interactions();

    /** Refresh marginal_interactions_(i,j)/(j,i) from pairwise_effects_discrete_ and the cached cross_term_. Valid only while pairwise_effects_cross_ and covariance_continuous_ are unchanged since the last cross_term_ refresh. */
    void refresh_marginal_interactions_entry(int i, int j);

    /** Rebuild Cholesky constraint structure and excluded-edge index lists. */
    void ensure_constraint_structure();

    /** Recompute theta_yy_ from cholesky_of_precision_ (inverse of the engine forward map). */
    void recompute_theta_yy() const;

    // =========================================================================
    // Gradient helpers (implemented in mixed_mrf_gradient.cpp)
    // =========================================================================

    /** Rebuild gradient index maps after edge-indicator changes. */
    void ensure_gradient_cache();

    /** Mark gradient cache as stale (call after edge-indicator changes). */
    void invalidate_gradient_cache();

    /** Unpack NUTS-vector into temporary parameter matrices (no model mutation). */
    void unvectorize_nuts_to_temps(
        const arma::vec& params,
        arma::mat& temp_main_discrete,
        arma::mat& temp_pairwise_discrete,
        arma::vec& temp_main_continuous,
        arma::mat& temp_pairwise_cross
    ) const;

    // =========================================================================
    // Likelihood functions (implemented in mixed_mrf_likelihoods.cpp)
    // =========================================================================

    /** Marginal OMRF pseudolikelihood for discrete variable s, using marginal_interactions_. */
    double log_marginal_omrf(int s) const;

    /** Conditional GGM log-likelihood: log f(y | x), using cached decomposition. */
    double log_conditional_ggm() const;

    // =========================================================================
    // MH update functions (implemented in mixed_mrf_metropolis.cpp)
    // =========================================================================

    // --- Rank-1 precision proposal helpers (permutation-free) ---

    // Extract reparameterization constants for the (i,j) off-diagonal precision update.
    // Populates cont_constants_[0..5] from cholesky_of_precision_ and covariance_continuous_.
    void get_precision_constants(int i, int j);

    // Constrained diagonal value for a proposed off-diagonal precision element.
    double precision_constrained_diagonal(double x) const;

    // Log-likelihood ratio for a proposed off-diagonal precision change (rank-2).
    // Assumes precision_proposal_ is already filled by the caller. Writes the
    // proposed covariance Σ' (computed via Woodbury) to cov_prop_out so callers
    // can use it to recompute marginal_interactions_ for the OMRF likelihood
    // ratio at the proposed Kyy.
    double log_ggm_ratio_edge(int i, int j, arma::mat& cov_prop_out) const;

    // Log-likelihood ratio for a proposed diagonal precision change (rank-1).
    // Assumes precision_proposal_ is already filled by the caller. Writes the
    // proposed covariance Σ' (computed via Sherman-Morrison) to cov_prop_out.
    double log_ggm_ratio_diag(int i, arma::mat& cov_prop_out) const;

    // log|Kyy_prop| - log|Kyy_curr| for a rank-2 off-diagonal proposal at
    // (i, j), via the matrix-determinant lemma in O(q). Reads
    // pairwise_effects_continuous_, precision_proposal_, and
    // covariance_continuous_; assumes precision_proposal_ has the proposed
    // Kyy at (i, j), (j, i), (j, j) already filled. Used to add the
    // determinant-tilt term delta_yy * (log|Kyy_prop| - log|Kyy_curr|) to MH
    // ratios.
    double log_det_ratio_yy_edge(int i, int j) const;

    // log|Kyy_prop| - log|Kyy_curr| for a rank-1 diagonal proposal at i.
    // Computed via the matrix-determinant lemma in O(1).
    double log_det_ratio_yy_diag(int i) const;

    // Rank-1 Cholesky update after accepting an off-diagonal precision change.
    void cholesky_update_after_precision_edge(double old_ij, double old_jj, int i, int j);

    // Rank-1 Cholesky update after accepting a diagonal precision change.
    void cholesky_update_after_precision_diag(double old_ii, int i);

    // --- Parameter update sweeps ---

    /**
     * Run every within-model MH proposal once (main effects, continuous
     * means, pairwise discrete/continuous/cross). Edge-indicator updates
     * are handled separately by ChainRunner.
     *
     * If `rm_weight` is set, each proposal also Robbins-Monro-updates its
     * proposal-SD slot using `ln_alpha` and `*rm_weight`. If nullopt, only
     * the accept/reject step runs.
     *
     * Shared between `do_one_metropolis_step` (called by MetropolisSampler
     * every iteration, with nullopt) and `tune_proposal_sd` (called every
     * iteration but a no-op outside stage 3b; during 3b passes the
     * schedule-supplied weight).
     */
    void sweep_within_model_mh(std::optional<double> rm_weight);

    // All within-model MH update sweeps below take an optional `rm_weight`
    // and return the Metropolis acceptance probability for the proposal.
    // - From `do_one_metropolis_step`: called with std::nullopt; the returned
    //   AR is collected per slot and fed to the MetropolisAdaptationController
    //   batch update after the sweep.
    // - From `tune_proposal_sd` during stage 3b: called with the schedule's
    //   `rm_weight_for_proposal_sd(iter)`; Robbins-Monro updates the matching
    //   proposal-SD slot inline using `ln_alpha` and `*rm_weight`. Returned
    //   AR is discarded.

    /** Update one main-effect: main_effects_discrete_(s, c). Ordinal threshold or BC α/β. */
    double update_main_effect(int s, int c, std::optional<double> rm_weight);

    /** Update one continuous mean: main_effects_continuous_(j). */
    double update_continuous_mean(int j, std::optional<double> rm_weight);

    /** Update one discrete interaction: pairwise_effects_discrete_(i, j). Symmetric. */
    double update_pairwise_discrete(int i, int j, std::optional<double> rm_weight);

    /** Update one off-diagonal precision element. Cholesky-based. */
    double update_pairwise_effects_continuous_offdiag(int i, int j, std::optional<double> rm_weight);

    /** Update one diagonal precision element. Log-scale Cholesky. */
    double update_pairwise_effects_continuous_diag(int i, std::optional<double> rm_weight);

    /** Update one cross interaction: pairwise_effects_cross_(i, j). */
    double update_pairwise_cross(int i, int j, std::optional<double> rm_weight);

    // --- Edge-indicator update sweeps ---

    /** Metropolis-Hastings add-delete move for one discrete-discrete edge. */
    void update_edge_indicator_discrete(int i, int j);

    /** Metropolis-Hastings add-delete move for one continuous-continuous edge. */
    void update_edge_indicator_continuous(int i, int j);

    /** Metropolis-Hastings add-delete move for one cross-type edge. */
    void update_edge_indicator_cross(int i, int j);

    // =========================================================================
    // Edge-indicator accessor helpers
    // =========================================================================

    int gxx(int i, int j) const { return edge_indicators_(i, j); }
    int gyy(int i, int j) const { return edge_indicators_(p_ + i, p_ + j); }
    int gxy(int i, int j) const { return edge_indicators_(i, p_ + j); }

    void set_gxx(int i, int j, int val) {
        edge_indicators_(i, j) = val;
        edge_indicators_(j, i) = val;
    }
    void set_gyy(int i, int j, int val) {
        edge_indicators_(p_ + i, p_ + j) = val;
        edge_indicators_(p_ + j, p_ + i) = val;
    }
    void set_gxy(int i, int j, int val) {
        edge_indicators_(i, p_ + j) = val;
        edge_indicators_(p_ + j, i) = val;
    }

    /**
     * Continuous-block adjacency Gamma_yy (q x q, unit diagonal) for the
     * Z-ratio engine's mediating-block extraction.
     */
    arma::imat continuous_subgraph() const {
        arma::imat g = edge_indicators_.submat(p_, p_, p_ + q_ - 1,
                                               p_ + q_ - 1);
        g.diag().ones();
        return g;
    }
};

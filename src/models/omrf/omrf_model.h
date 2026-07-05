#pragma once

#include <memory>
#include <functional>
#include "models/base_model.h"
#include "mcmc/samplers/metropolis_adaptation.h"
#include "rng/rng_utils.h"
#include "mcmc/execution/step_result.h"
#include "utils/common_helpers.h"
#include "utils/variable_helpers.h"
#include "priors/parameter_prior.h"

/**
 * OMRFModel - Ordinal Markov Random Field Model
 *
 * A class-based implementation of the OMRF model for Bayesian inference on
 * ordinal and Blume-Capel variables. This class encapsulates:
 *   - Parameter storage (main effects, pairwise effects, edge indicators)
 *   - Sufficient statistics computation
 *   - Log-pseudoposterior and gradient evaluations
 *   - Adaptive Metropolis-Hastings updates for individual parameters
 *   - NUTS updates for joint parameter sampling
 *   - Edge selection (spike-and-slab) with asymmetric proposals
 *
 * Inherits from BaseModel for compatibility with the generic MCMC framework.
 */
class OMRFModel : public BaseModel {
public:

    /**
     * Constructor from raw observations
     *
     * @param observations        Integer matrix of categorical observations (persons × variables)
     * @param num_categories      Number of categories per variable
     * @param inclusion_probability Prior inclusion probabilities for edges
     * @param initial_edge_indicators Initial edge inclusion matrix
     * @param is_ordinal_variable Indicator (1 = ordinal, 0 = Blume-Capel)
     * @param baseline_category   Reference categories for Blume-Capel variables
     * @param main_alpha          Beta prior hyperparameter α for main effects
     * @param main_beta           Beta prior hyperparameter β for main effects
     * @param pairwise_scale      Scale parameter of Cauchy prior on interactions
     * @param edge_selection      Enable edge selection (spike-and-slab)
     */
    OMRFModel(
        const arma::imat& observations,
        const arma::ivec& num_categories,
        const arma::mat& inclusion_probability,
        const arma::imat& initial_edge_indicators,
        const arma::uvec& is_ordinal_variable,
        const arma::ivec& baseline_category,
        std::unique_ptr<BaseParameterPrior> interaction_prior,
        std::unique_ptr<BaseParameterPrior> threshold_prior,
        bool edge_selection
    );

    /**
     * Copy constructor for cloning (required for parallel chains)
     */
    OMRFModel(const OMRFModel& other);

    // =========================================================================
    // BaseModel interface implementation
    // =========================================================================

    /** @return true when edge selection is enabled. */
    bool has_edge_selection() const override { return edge_selection_; }
    /** @return true when missing-data imputation is active. */
    bool has_missing_data() const override { return has_missing_; }

    /**
     * Combined log-posterior and gradient evaluation (more efficient)
     */
    std::pair<double, arma::vec> logp_and_gradient(const arma::vec& parameters) override;

    /**
     * Perform one adaptive Metropolis step (updates all parameters)
     * @param iteration  Current iteration (for Robbins-Monro adaptation)
     */
    void do_one_metropolis_step(int iteration = -1) override;

    /**
     * @return Mean Metropolis acceptance probability over the most recent
     *         do_one_metropolis_step() sweep. NaN before the first step.
     */
    double last_metropolis_mean_accept_prob() const override {
        return last_mh_mean_accept_;
    }

    /**
     * Set the Robbins-Monro target acceptance rate used by the
     * adaptive-Metropolis sampler. Honoured by both
     * init_metropolis_adaptation() (forwarded to the controllers) and
     * tune_proposal_sd() (the post-step Robbins-Monro update).
     */
    void set_metropolis_target_accept(double target) override {
        target_accept_ = target;
    }

    /**
     * Initialize Metropolis adaptation controllers for proposal-SD tuning
     * Must be called before warmup begins (e.g., by MetropolisSampler on first step)
     */
    void init_metropolis_adaptation(const WarmupSchedule& schedule) override;

    /**
     * Stage 3b: tune pairwise proposal SDs via Robbins-Monro
     * Called every iteration from the runner; checks schedule internally
     */
    void tune_proposal_sd(int iteration, const WarmupSchedule& schedule) override;

    /**
     * Return dimensionality of active parameter space
     */
    size_t parameter_dimension() const override;

    /**
     * Set random seed for reproducibility
     */
    void set_seed(int seed) override;

    /**
     * Get vectorized parameters (main effects + active pairwise effects)
     */
    arma::vec get_vectorized_parameters() const override;

    /**
     * Set parameters from vectorized form
     */
    void set_vectorized_parameters(const arma::vec& parameters) override;

    /**
     * Get vectorized edge indicators
     */
    arma::ivec get_vectorized_indicator_parameters() override;

    /**
     * Clone the model for parallel execution.
     */
    std::unique_ptr<BaseModel> clone() const override;

    /** @return Reference to the model's random number generator. */
    SafeRNG& get_rng() override { return rng_; }

    // =========================================================================
    // OMRF-specific methods
    // =========================================================================



    /**
     * Shuffle edge update order at the start of each iteration.
     * The shuffled order is stored in shuffled_edge_order_ for use by
     * update_edge_indicators(). Called unconditionally to advance the
     * RNG state consistently.
     */
    void prepare_iteration() override;

    /**
     * Update edge indicators via Metropolis-Hastings
     */
    void update_edge_indicators() override;

    /**
     * Impute missing values (if any)
     */
    void impute_missing() override;

    /**
     * Set missing data information
     */
    void set_missing_data(const arma::imat& missing_index);

    // =========================================================================
    // Accessors
    // =========================================================================

    /** @return Current main-effect parameter matrix (p x max_cats). */
    const arma::mat& get_main_effects() const { return main_effects_; }
    /** @return Current pairwise interaction matrix (p x p, symmetric). */
    const arma::mat& get_pairwise_effects() const { return pairwise_effects_; }
    /** @return Current edge-indicator matrix (p x p, symmetric, 0/1). */
    const arma::imat& get_edge_indicators() const override { return edge_indicators_; }
    /** @return Mutable reference to the prior inclusion-probability matrix. */
    arma::mat& get_inclusion_probability() override { return inclusion_probability_; }
    /** @return Residual matrix X * pairwise_effects (n x p). */
    const arma::mat& get_residual_matrix() const { return residual_matrix_; }

    /**
     * Replace all main-effect parameters.
     * @param main_effects  New main-effect matrix (p x max_cats)
     */
    void set_main_effects(const arma::mat& main_effects) { main_effects_ = main_effects; }
    /**
     * Replace all pairwise effects and update residuals.
     * @param pairwise_effects  New pairwise interaction matrix (p x p)
     */
    void set_pairwise_effects(const arma::mat& pairwise_effects);
    /**
     * Replace all edge indicators.
     * @param edge_indicators  New edge-indicator matrix (p x p)
     */
    void set_edge_indicators(const arma::imat& edge_indicators) { edge_indicators_ = edge_indicators; }

    /** @return Number of variables (p) as int. */
    int get_num_variables() const override { return static_cast<int>(p_); }
    /** @return Number of unique off-diagonal pairs p(p-1)/2 as int. */
    int get_num_pairwise() const override { return static_cast<int>(num_pairwise_); }
    /** @return Number of variables (p). */
    size_t num_variables() const { return p_; }
    /** @return Number of observations (n). */
    size_t num_observations() const { return n_; }
    /** @return Total number of main-effect parameters across all variables. */
    size_t num_main_effects() const { return num_main_; }
    /** @return Number of unique pairwise interactions p(p-1)/2. */
    size_t num_pairwise_effects() const { return num_pairwise_; }

    /** @return Number of variables (shorthand for interface compatibility). */
    size_t get_p() const { return p_; }
    /** @return Number of observations (shorthand for interface compatibility). */
    size_t get_n() const { return n_; }

    /**
     * Set the inverse mass matrix diagonal for NUTS.
     * @param inv_mass  Diagonal elements of the inverse mass matrix
     */
    void set_inv_mass(const arma::vec& inv_mass) override { inv_mass_ = inv_mass; }
    /** @return Current inverse mass matrix diagonal. */
    const arma::vec& get_inv_mass() const override { return inv_mass_; }

    /**
     * Get full dimension (main + ALL pairwise, regardless of edge indicators)
     * Used for fixed-size sample storage
     */
    size_t full_parameter_dimension() const override { return num_main_ + num_pairwise_; }

    /**
     * Get all parameters in a fixed-size vector (inactive edges are 0)
     * Used for sample storage to avoid dimension changes
     */
    arma::vec get_full_vectorized_parameters() const override;

    /** @return Mutable reference to main-effect proposal SDs (for external adaptation). */
    arma::mat& get_proposal_sd_main() { return proposal_sd_main_; }
    /** @return Mutable reference to pairwise proposal SDs (for external adaptation). */
    arma::mat& get_proposal_sd_pairwise() { return proposal_sd_pairwise_; }

    /**
     * Enable or disable edge-selection proposals.
     * @param active  true to enable edge add-delete moves
     */
    void set_edge_selection_active(bool active) override { edge_selection_active_ = active; }
    /** @return true when edge-selection proposals are currently active. */
    bool is_edge_selection_active() const { return edge_selection_active_; }

private:
    // =========================================================================
    // Data members
    // =========================================================================

    // Most recent mean Metropolis acceptance probability over all updated
    // pairwise + main-effect components. Reset every do_one_metropolis_step().
    double last_mh_mean_accept_ = std::numeric_limits<double>::quiet_NaN();

    // Robbins-Monro target acceptance rate for adaptive-Metropolis
    // proposal-SD tuning. Set via set_metropolis_target_accept(), defaults
    // to 0.44 (the componentwise random-walk Metropolis optimum).
    double target_accept_ = 0.44;

    // Data
    size_t n_;                          ///< Number of observations
    size_t p_;                          ///< Number of variables
    arma::imat observations_;           ///< Categorical observations (n x p)
    arma::mat observations_double_;     ///< Observations as double (for efficient matrix ops)
    arma::mat observations_double_t_;   ///< Transposed observations (for BLAS pairwise gradient)
    arma::ivec num_categories_;         ///< Categories per variable
    arma::uvec is_ordinal_variable_;    ///< 1 = ordinal, 0 = Blume-Capel
    arma::ivec baseline_category_;      ///< Reference category for Blume-Capel

    // Sufficient statistics
    arma::imat counts_per_category_;    ///< Category counts (max_cats+1 x p)
    arma::imat blume_capel_stats_;      ///< [linear_sum, quadratic_sum] for BC vars (2 x p)
    arma::imat pairwise_stats_;         ///< X^T X
    arma::mat residual_matrix_;         ///< X * pairwise_effects (n x p)

    // Per-variable log normalizer sum_i (bound + log denom) at the current
    // state. Refreshed at the top of each MH/indicator sweep and maintained
    // on accept, so proposals stop re-evaluating the current state.
    arma::vec log_denominator_cache_;   ///< p

    // Parameters
    arma::mat main_effects_;            ///< Main effect parameters (p x max_cats)
    arma::mat pairwise_effects_;        ///< Pairwise interactions (p x p, symmetric)
    arma::imat edge_indicators_;        ///< Edge inclusion indicators (p x p, symmetric binary)

    // Priors
    arma::mat inclusion_probability_;   ///< Prior inclusion probabilities
    std::unique_ptr<BaseParameterPrior> interaction_prior_; ///< Prior on pairwise interactions
    std::unique_ptr<BaseParameterPrior> threshold_prior_;  ///< Prior on main effects / thresholds

    // Model configuration
    bool edge_selection_;               ///< Enable edge selection
    bool edge_selection_active_;        ///< Currently in edge selection phase

    // Dimension tracking
    size_t num_main_;                   ///< Total number of main effect parameters
    size_t num_pairwise_;               ///< Number of possible pairwise effects

    // Proposal SDs (adapted by MetropolisAdaptationController during warmup)
    arma::mat proposal_sd_main_;        ///< Proposal SD for main effects
    arma::mat proposal_sd_pairwise_;    ///< Proposal SD for pairwise effects

    // Metropolis adaptation controllers (created by init_metropolis_adaptation)
    std::unique_ptr<MetropolisAdaptationController> metropolis_main_adapter_;      ///< Main-effect adapter
    std::unique_ptr<MetropolisAdaptationController> metropolis_pairwise_adapter_;  ///< Pairwise-effect adapter

    // RNG
    SafeRNG rng_;                       ///< Per-chain random number generator

    // NUTS settings
    arma::vec inv_mass_;                ///< Inverse mass diagonal

    // Missing data handling
    bool has_missing_;                  ///< Whether the data contains missing values
    arma::imat missing_index_;          ///< (row, col) indices of missing entries

    // Cached gradient components
    arma::vec grad_obs_cache_;          ///< Cached observed-data gradient
    arma::imat index_matrix_cache_;     ///< Cached parameter index map
    bool gradient_cache_valid_;         ///< Whether the gradient cache is current

    // Per-chain scratch for compute_logZ_and_probs_*_into. Reused across
    // every call to logp_and_gradient and across every variable inside it.
    // After the first few calls the buffers stabilise at max size and no
    // further heap allocations happen on the hot path.
    mutable LogZAndProbs logz_out_;
    mutable LogZScratch  logz_scratch_;

    // Interaction indexing (for edge updates)
    arma::imat interaction_index_;      ///< Maps edge pair to index
    arma::uvec shuffled_edge_order_;    ///< Pre-shuffled order (set in prepare_iteration)

    // =========================================================================
    // Private helper methods
    // =========================================================================

    /**
     * Compute sufficient statistics from observations
     */
    void compute_sufficient_statistics();

    /**
     * Count total number of main effect parameters
     */
    size_t count_num_main_effects_internal() const;

    /**
     * Build interaction index matrix
     */
    void build_interaction_index();

    /**
     * Update residual matrix after pairwise effects change
     */
    void update_residual_matrix();

    /**
     * Incrementally update two residual columns after a single pairwise effect change
     */
    void update_residual_columns(int var1, int var2, double delta);

    /**
     * Invalidate gradient cache (call after parameter changes)
     */
    void invalidate_gradient_cache() { gradient_cache_valid_ = false; }

    /**
     * Ensure gradient cache is valid
     */
    void ensure_gradient_cache();

    // -------------------------------------------------------------------------
    // Log-posterior components
    // -------------------------------------------------------------------------

    /**
     * Log normalizer sum_i (bound + log denom) for one variable at the
     * current state; the expensive shared piece of every MH acceptance.
     */
    double compute_log_denominator(int variable) const;

    /**
     * Log normalizer for variable under a pairwise shift: the residual
     * column moves by 2 * obs_other * delta.
     */
    double compute_log_denominator_shifted(
        int variable, const arma::vec& obs_other, double delta) const;

    /** Refill log_denominator_cache_ for all variables. */
    void recompute_log_denominators();

    /**
     * One cached MH step on pairwise effect (var1, var2): shared by
     * update_pairwise_effect and the stage-3b tuner (which runs it without
     * the edge gate). Returns the acceptance probability.
     */
    double mh_pairwise_step(int var1, int var2);

    // -------------------------------------------------------------------------
    // Parameter vectorization
    // -------------------------------------------------------------------------

    /**
     * Flatten parameters to vector
     */
    arma::vec vectorize_parameters() const;

    /**
     * Flatten parameters into pre-allocated vector (avoids allocation)
     */
    void vectorize_parameters_into(arma::vec& param_vec) const;

    /**
     * Unflatten vector to parameter matrices
     */
    void unvectorize_parameters(const arma::vec& param_vec);

    /**
     * Unvectorize a parameter vector into temporary main/pairwise matrices,
     * then compute the corresponding residual matrix.
     */
    void unvectorize_to_temps(
        const arma::vec& parameters,
        arma::mat& temp_main,
        arma::mat& temp_pairwise,
        arma::mat& temp_residual
    ) const;

    /**
     * Extract active inverse mass (only for included edges)
     */
    arma::vec get_active_inv_mass() const override;

    /**
     * Extract active inverse mass into pre-allocated vector (avoids allocation)
     */
    void get_active_inv_mass_into(arma::vec& active_inv_mass) const;

    // -------------------------------------------------------------------------
    // Metropolis updates
    // -------------------------------------------------------------------------

    /**
     * Update single main effect parameter via Metropolis
     * @return acceptance probability (for Metropolis adaptation)
     */
    double update_main_effect_parameter(int variable, int category, int parameter);

    /**
     * Update single pairwise effect via Metropolis
     * @return acceptance probability (for Metropolis adaptation)
     */
    double update_pairwise_effect(int var1, int var2);

    /**
     * Update single edge indicator (spike-and-slab)
     */
    void update_edge_indicator(int var1, int var2);
};


/**
 * Factory function to create OMRFModel from R inputs
 */
OMRFModel createOMRFModelFromR(
    const Rcpp::List& inputFromR,
    const arma::mat& inclusion_probability,
    const arma::imat& initial_edge_indicators,
    std::unique_ptr<BaseParameterPrior> interaction_prior,
    std::unique_ptr<BaseParameterPrior> threshold_prior,
    bool edge_selection = true
);

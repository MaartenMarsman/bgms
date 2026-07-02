#pragma once

#include <RcppArmadillo.h>
#include <stdexcept>
#include <memory>
#include <limits>

// Forward declarations
struct StepResult;
struct SafeRNG;
struct WarmupSchedule;

/**
 * BaseModel — Abstract interface for all graphical models.
 *
 * Defines the virtual methods that the MCMC framework (MetropolisSampler,
 * NUTSSampler, ChainRunner) calls during sampling. Most methods are pure
 * virtual (`= 0`) so the compiler enforces implementation in every
 * subclass. The exceptions are gradient-only methods (logp_and_gradient,
 * set_vectorized_parameters) which throw at runtime for models that do not
 * support NUTS.
 *
 * Subclass hierarchy:
 *   - GGMModel      — Gaussian Graphical Model (precision matrix, Metropolis + NUTS)
 *   - OMRFModel     — Ordinal Markov Random Field (Metropolis + NUTS)
 *   - MixedMRFModel — Mixed discrete + continuous MRF (Metropolis + NUTS)
 *
 * Methods fall into several groups:
 *   - **Capability queries** (has_edge_selection, has_constraints, has_missing_data)
 *     — run/data-dependent toggles the runner and samplers branch on. (Sampler
 *     choice itself is config-driven, not capability-driven.)
 *   - **Sampling steps** (do_one_metropolis_step, logp_and_gradient)
 *     — called by the sampler each iteration.
 *   - **Edge selection** (update_edge_indicators)
 *     — spike-and-slab structure learning.
 *   - **Parameter access** (get/set_vectorized_parameters, get_full_vectorized_parameters)
 *     — used by NUTS for momentum-based proposals and by the runner for output.
 *   - **Adaptation** (set_inv_mass, tune_proposal_sd)
 *     — tuned during warmup.
 *   - **Missing data** (has_missing_data, impute_missing)
 *     — full-conditional imputation between iterations.
 */
class BaseModel {
public:
    virtual ~BaseModel() = default;

    // =========================================================================
    // Capability queries
    // =========================================================================

    /** @return true if the model supports edge selection (spike-and-slab). */
    virtual bool has_edge_selection() const { return false; }

    // =========================================================================
    // Core sampling methods
    // =========================================================================

    /**
     * Combined log-(pseudo)posterior and gradient evaluation.
     *
     * More efficient than calling log-posterior and gradient() separately
     * because intermediate computations can be shared.
     *
     * @param parameters  Vectorized parameters at which to evaluate
     * @return Pair of (log-posterior value, gradient vector)
     */
    virtual std::pair<double, arma::vec> logp_and_gradient(
        const arma::vec& parameters) {
        throw std::runtime_error("logp_and_gradient not implemented for this model");
    }

    /**
     * Perform one full Metropolis sweep over all parameters.
     *
     * The model handles its own parameter grouping (e.g. off-diagonal,
     * diagonal, edge indicators in GGM; main, pairwise in OMRF).
     *
     * @param iteration  Current iteration index (for Robbins-Monro adaptation)
     */
    virtual void do_one_metropolis_step(int iteration = -1) = 0;

    /**
     * Perform one full Gibbs sweep over all parameters.
     *
     * Default throws: only models with an exact conjugate full-conditional
     * sweep (currently GGM row-block Gibbs) override this. The GibbsSampler
     * wrapper is only constructed for models that support it.
     *
     * @param iteration  Current iteration index (unused by exact samplers;
     *                   kept for interface symmetry with the Metropolis step).
     */
    virtual void do_one_gibbs_step(int /*iteration*/ = -1) {
        throw std::runtime_error("do_one_gibbs_step not implemented for this model");
    }

    /**
     * Enable the full-conditional edge birth/death proposal for the between-model
     * step, used by the Gibbs sampler in place of the random-walk Roverato
     * proposal (which needs tuning). Default no-op; only the GGM overrides it.
     */
    virtual void set_conjugate_edge_proposal(bool /*enable*/) {}

    /**
     * Mean Metropolis acceptance probability across all components updated
     * in the most recent do_one_metropolis_step() call.
     *
     * Surfaced as the per-iteration `am_accept_prob__` trace for the
     * adaptive-Metropolis sampler. This is not `fit$am_diag$accept_prob`,
     * which summarize_am_diagnostics computes in R as a per-parameter
     * empirical move rate. Defaults to NaN for models that do not implement
     * Metropolis updates or have not yet taken a step.
     *
     * @return Mean acceptance probability over the last sweep, or NaN.
     */
    virtual double last_metropolis_mean_accept_prob() const {
        return std::numeric_limits<double>::quiet_NaN();
    }

    /**
     * Set the target Metropolis acceptance rate for Robbins-Monro proposal
     * adaptation. Called by the sampler entry points (sample_omrf,
     * sample_mixed_mrf, sample_ggm) before the MCMC loop, with the value
     * the user passed to bgm()'s `target_accept` argument.
     *
     * Default: no-op for models that do not perform Metropolis adaptation.
     */
    virtual void set_metropolis_target_accept(double /*target*/) {}

    /**
     * Initialize Metropolis adaptation controllers.
     *
     * Called once before the MCMC loop begins. Subclasses store the
     * warmup schedule for later use by Robbins-Monro adaptation.
     *
     * @param schedule  Warmup schedule defining adaptation phases
     */
    virtual void init_metropolis_adaptation(const WarmupSchedule& /*schedule*/) {}

    /**
     * Tune pairwise proposal SDs via Robbins-Monro (warmup Stage 3b).
     *
     * Called every iteration from the runner; the implementation checks
     * the schedule internally to decide whether to adapt.
     *
     * @param iteration  Current iteration index
     * @param schedule   Warmup schedule
     */
    virtual void tune_proposal_sd(int /*iteration*/, const WarmupSchedule& /*schedule*/) {}

    /**
     * Called at the start of every iteration before edge selection and sampling.
     *
     * Subclasses may use this to shuffle edge update order or advance the
     * RNG state consistently.
     */
    virtual void prepare_iteration() {}

    /**
     * Called once at the warmup/sampling boundary (before the first
     * retained iteration). Default no-op; models with warm-up-calibrated
     * state (e.g. the GGM Z-ratio engine) freeze it here.
     */
    virtual void on_warmup_end() {}

    // =========================================================================
    // Edge selection
    // =========================================================================

    /**
     * Update edge indicators via Metropolis-Hastings add-delete moves.
     *
     * Only meaningful when has_edge_selection() returns true. Each derived
     * model implements its own edge-indicator sweep here.
     */
    virtual void update_edge_indicators() = 0;

    // =========================================================================
    // Parameter vectorization
    // =========================================================================

    /** @return Active parameters as a flat vector (dimension may change with edge selection). */
    virtual arma::vec get_vectorized_parameters() const = 0;

    /**
     * Set parameters from a flat vector (inverse of get_vectorized_parameters).
     * @param parameters  Vectorized parameter values
     */
    virtual void set_vectorized_parameters(const arma::vec& parameters) {
        throw std::runtime_error("set_vectorized_parameters method must be implemented in derived class");
    }

    /** @return Edge indicators as a flat integer vector. */
    virtual arma::ivec get_vectorized_indicator_parameters() = 0;

    /**
     * @return Full parameter dimension (fixed size, includes inactive parameters).
     *
     * Used by NUTSSampler for mass-matrix sizing and adaptation.
     * For most models this equals the storage dimension. For models where
     * some parameters are not sampled by NUTS (e.g., MixedMRFModel's
     * continuous precision),
     * this returns the NUTS-block dimension.
     * Defaults to parameter_dimension().
     */
    virtual size_t full_parameter_dimension() const {
        return parameter_dimension();
    }

    /**
     * @return All parameters in a fixed-size vector (inactive edges are 0).
     *
     * Used by NUTSSampler for adaptation (online covariance).
     * Dimension must match full_parameter_dimension().
     */
    virtual arma::vec get_full_vectorized_parameters() const = 0;

    /** @return Dimensionality of the active parameter space. Pure virtual. */
    virtual size_t parameter_dimension() const = 0;

    /**
     * @return Dimension for sample storage (includes all parameters).
     *
     * For most models this equals full_parameter_dimension(). Override
     * when storage needs more entries than the NUTS block (e.g., continuous
     * precision parameters in MixedMRFModel).
     */
    virtual size_t storage_dimension() const {
        return full_parameter_dimension();
    }

    /**
     * @return All parameters in a fixed-size vector for sample storage.
     *
     * Dimension must match storage_dimension(). Default delegates to
     * get_full_vectorized_parameters().
     */
    virtual arma::vec get_storage_vectorized_parameters() const {
        return get_full_vectorized_parameters();
    }

    // =========================================================================
    // Infrastructure
    // =========================================================================

    /**
     * Set the random seed for reproducibility.
     * @param seed  Integer seed value
     */
    virtual void set_seed(int seed) = 0;

    /** @return Deep copy of this model (for parallel chains). */
    virtual std::unique_ptr<BaseModel> clone() const = 0;

    /** @return Reference to the model's random number generator. */
    virtual SafeRNG& get_rng() = 0;

    // =========================================================================
    // NUTS adaptation
    // =========================================================================

    /**
     * Set the inverse mass matrix diagonal for NUTS.
     * @param inv_mass  Diagonal elements of the inverse mass matrix
     */
    virtual void set_inv_mass(const arma::vec& inv_mass) { inv_mass_ = inv_mass; }
    /** @return Current inverse mass matrix diagonal. */
    virtual const arma::vec& get_inv_mass() const { return inv_mass_; }

    /**
     * @return Active subset of the inverse mass diagonal.
     *
     * For models with edge selection, this may return only the entries
     * corresponding to included edges. Default: returns the full diagonal.
     */
    virtual arma::vec get_active_inv_mass() const { return inv_mass_; }

    // =========================================================================
    // RATTLE constrained integration
    // =========================================================================

    /** @return true if the model has constraints requiring RATTLE projection. */
    virtual bool has_constraints() const { return false; }

    /**
     * Full-dimension position for RATTLE integration.
     *
     * Returns the current parameters in the full (fixed-dimension)
     * coordinate system used by RATTLE. Default: delegates to
     * get_vectorized_parameters().
     */
    virtual arma::vec get_full_position() const {
        return get_vectorized_parameters();
    }

    /**
     * Set model state from a full-dimension RATTLE position vector.
     *
     * Updates all derived matrices (precision, covariance, etc.) from
     * the full position vector. Default: delegates to
     * set_vectorized_parameters().
     *
     * @param x  Full-dimension position vector
     */
    virtual void set_full_position(const arma::vec& x) {
        set_vectorized_parameters(x);
    }

    /**
     * Full-space log-posterior and gradient for RATTLE integration.
     *
     * Evaluates the log-posterior and its gradient in the full
     * (fixed-dimension) parameter space. Default: delegates to
     * logp_and_gradient().
     *
     * @param x  Full-dimension position vector
     * @return (log-posterior value, gradient vector)
     */
    virtual std::pair<double, arma::vec> logp_and_gradient_full(
        const arma::vec& x) {
        return logp_and_gradient(x);
    }

    /**
     * Project position onto the constraint manifold (in-place).
     *
     * For models with linear constraints on the parameter space (e.g., zero
     * entries in the precision matrix), this modifies x so that all
     * constraints are satisfied. Default: no-op.
     *
     * @param x  Full-dimension position vector (modified in-place)
     */
    virtual void project_position(arma::vec& x) const { (void)x; }

    /**
     * Project position onto the constraint manifold (mass-weighted SHAKE).
     *
     * Uses the correction direction M^{-1} J^T for RATTLE-correct
     * symplecticity. Default: delegates to identity-mass overload.
     *
     * @param x              Full-dimension position vector (modified)
     * @param inv_mass_diag  Diagonal of the inverse mass matrix
     */
    virtual void project_position(arma::vec& x,
                                  const arma::vec& inv_mass_diag) const {
        (void)inv_mass_diag;
        project_position(x);
    }

    /**
     * Project momentum onto the cotangent space of the constraint manifold.
     *
     * Ensures the momentum vector lies in the tangent space of the
     * constraint surface at the current position. Default: no-op.
     *
     * @param r  Momentum vector (modified in-place)
     * @param x  Current position (after projection)
     */
    virtual void project_momentum(arma::vec& r, const arma::vec& x) const {
        (void)r; (void)x;
    }

    /**
     * Project momentum onto the cotangent space (mass-weighted RATTLE).
     *
     * Enforces J M^{-1} r = 0 for RATTLE-correct symplecticity.
     * Default: delegates to identity-mass overload.
     *
     * @param r              Momentum vector (modified in-place)
     * @param x              Current position (after projection)
     * @param inv_mass_diag  Diagonal of the inverse mass matrix
     */
    virtual void project_momentum(arma::vec& r, const arma::vec& x,
                                  const arma::vec& inv_mass_diag) const {
        (void)inv_mass_diag;
        project_momentum(r, x);
    }

    /** Reset any cached state used by projection solvers. */
    virtual void reset_projection_cache() {}

    // =========================================================================
    // Edge selection control
    // =========================================================================

    /**
     * Enable or disable edge-selection proposals.
     * @param active  true to enable edge add-delete moves
     */
    virtual void set_edge_selection_active(bool active) {
        (void)active;
    }

    // =========================================================================
    // Missing data
    // =========================================================================

    /** @return true when missing-data imputation is active. */
    virtual bool has_missing_data() const { return false; }

    /** Impute missing entries from full-conditional distributions. */
    virtual void impute_missing() {}

    // =========================================================================
    // Edge prior support
    // =========================================================================

    /** @return Current edge-indicator matrix. */
    virtual const arma::imat& get_edge_indicators() const = 0;

    /** @return Mutable reference to the prior inclusion-probability matrix. */
    virtual arma::mat& get_inclusion_probability() = 0;

    /** @return Number of variables (p). */
    virtual int get_num_variables() const = 0;

    /** @return Number of unique off-diagonal pairs p(p-1)/2. */
    virtual int get_num_pairwise() const = 0;

protected:
    BaseModel() = default;
    /// Inverse mass matrix diagonal for NUTS.
    arma::vec inv_mass_;
};

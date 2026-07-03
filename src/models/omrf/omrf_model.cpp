#include <RcppArmadillo.h>
#include "models/omrf/omrf_model.h"
#include "rng/rng_utils.h"
#include "mcmc/algorithms/hmc.h"
#include "mcmc/algorithms/nuts.h"
#include "mcmc/algorithms/metropolis.h"
#include "mcmc/execution/step_result.h"
#include "mcmc/samplers/metropolis_adaptation.h"
#include "mcmc/execution/chain_runner.h"
#include "math/explog_macros.h"
#include "utils/common_helpers.h"
#include "utils/variable_helpers.h"


// =============================================================================
// Constructor
// =============================================================================

OMRFModel::OMRFModel(
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::mat& inclusion_probability,
    const arma::imat& initial_edge_indicators,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    std::unique_ptr<BaseParameterPrior> interaction_prior,
    std::unique_ptr<BaseParameterPrior> threshold_prior,
    bool edge_selection
) :
    n_(observations.n_rows),
    p_(observations.n_cols),
    observations_(observations),
    num_categories_(num_categories),
    is_ordinal_variable_(is_ordinal_variable),
    baseline_category_(baseline_category),
    inclusion_probability_(inclusion_probability),
    interaction_prior_(std::move(interaction_prior)),
    threshold_prior_(std::move(threshold_prior)),
    edge_selection_(edge_selection),
    edge_selection_active_(false),
    has_missing_(false),
    gradient_cache_valid_(false)
{
    // Initialize parameter dimensions
    num_main_ = count_num_main_effects_internal();
    num_pairwise_ = (p_ * (p_ - 1)) / 2;

    // Initialize parameters
    int max_cats = num_categories_.max();
    main_effects_ = arma::zeros<arma::mat>(p_, max_cats);
    pairwise_effects_ = arma::zeros<arma::mat>(p_, p_);
    edge_indicators_ = initial_edge_indicators;

    // Initialize proposal SDs
    proposal_sd_main_ = arma::ones<arma::mat>(p_, max_cats);
    proposal_sd_pairwise_ = arma::ones<arma::mat>(p_, p_);

    // Initialize mass matrix
    inv_mass_ = arma::ones<arma::vec>(num_main_ + num_pairwise_);

    // Center observations for Blume-Capel variables (x - baseline) so that
    // ALL downstream code — sufficient statistics, residuals, gradients,
    // log-pseudoposterior, imputation — operates in the same coordinate
    // system. For ordinal variables baseline=0, so this is a no-op.
    for (size_t v = 0; v < p_; ++v) {
        if (!is_ordinal_variable_(v)) {
            observations_.col(v) -= baseline_category_(v);
        }
    }
    observations_double_ = arma::conv_to<arma::mat>::from(observations_);
    observations_double_t_ = observations_double_.t();

    // Compute sufficient statistics
    compute_sufficient_statistics();

    // Initialize residual matrix
    update_residual_matrix();

    // Build interaction index
    build_interaction_index();
}


// =============================================================================
// Copy constructor
// =============================================================================

OMRFModel::OMRFModel(const OMRFModel& other)
    : BaseModel(other),
      last_mh_mean_accept_(other.last_mh_mean_accept_),
      target_accept_(other.target_accept_),
      n_(other.n_),
      p_(other.p_),
      observations_(other.observations_),
      observations_double_(other.observations_double_),
      observations_double_t_(other.observations_double_t_),
      num_categories_(other.num_categories_),
      is_ordinal_variable_(other.is_ordinal_variable_),
      baseline_category_(other.baseline_category_),
      counts_per_category_(other.counts_per_category_),
      blume_capel_stats_(other.blume_capel_stats_),
      pairwise_stats_(other.pairwise_stats_),
      residual_matrix_(other.residual_matrix_),
      main_effects_(other.main_effects_),
      pairwise_effects_(other.pairwise_effects_),
      edge_indicators_(other.edge_indicators_),
      inclusion_probability_(other.inclusion_probability_),
      interaction_prior_(other.interaction_prior_->clone()),
      threshold_prior_(other.threshold_prior_->clone()),
      edge_selection_(other.edge_selection_),
      edge_selection_active_(other.edge_selection_active_),
      num_main_(other.num_main_),
      num_pairwise_(other.num_pairwise_),
      proposal_sd_main_(other.proposal_sd_main_),
      proposal_sd_pairwise_(other.proposal_sd_pairwise_),
      rng_(other.rng_),
      inv_mass_(other.inv_mass_),
      has_missing_(other.has_missing_),
      missing_index_(other.missing_index_),
      grad_obs_cache_(other.grad_obs_cache_),
      index_matrix_cache_(other.index_matrix_cache_),
      gradient_cache_valid_(other.gradient_cache_valid_),
      interaction_index_(other.interaction_index_),
      shuffled_edge_order_(other.shuffled_edge_order_)
{
}


// =============================================================================
// Sufficient statistics computation
// =============================================================================

void OMRFModel::compute_sufficient_statistics() {
    int max_cats = num_categories_.max();

    // Category counts for ordinal variables
    counts_per_category_ = arma::zeros<arma::imat>(max_cats + 1, p_);
    for (size_t v = 0; v < p_; ++v) {
        if (is_ordinal_variable_(v)) {
            for (size_t i = 0; i < n_; ++i) {
                int cat = observations_(i, v);
                if (cat >= 0 && cat <= num_categories_(v)) {
                    counts_per_category_(cat, v)++;
                }
            }
        }
    }

    // Blume-Capel statistics (linear and quadratic sums)
    blume_capel_stats_ = arma::zeros<arma::imat>(2, p_);
    for (size_t v = 0; v < p_; ++v) {
        if (!is_ordinal_variable_(v)) {
            for (size_t i = 0; i < n_; ++i) {
                int s = observations_(i, v);         // already centered
                blume_capel_stats_(0, v) += s;       // linear
                blume_capel_stats_(1, v) += s * s;   // quadratic
            }
        }
    }

    // Pairwise statistics (X^T X) - use pre-computed transformed observations
    arma::mat ps = observations_double_.t() * observations_double_;
    pairwise_stats_ = arma::conv_to<arma::imat>::from(ps);
}


size_t OMRFModel::count_num_main_effects_internal() const {
    size_t count = 0;
    for (size_t v = 0; v < p_; ++v) {
        if (is_ordinal_variable_(v)) {
            count += num_categories_(v);
        } else {
            count += 2;  // linear and quadratic for Blume-Capel
        }
    }
    return count;
}


void OMRFModel::build_interaction_index() {
    interaction_index_ = arma::zeros<arma::imat>(num_pairwise_, 3);
    int idx = 0;
    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            interaction_index_(idx, 0) = idx;
            interaction_index_(idx, 1) = v1;
            interaction_index_(idx, 2) = v2;
            idx++;
        }
    }
}


void OMRFModel::update_residual_matrix() {
    residual_matrix_ = 2.0 * observations_double_ * pairwise_effects_;
}


void OMRFModel::update_residual_columns(int var1, int var2, double delta) {
    residual_matrix_.col(var1) += 2.0 * delta * observations_double_.col(var2);
    residual_matrix_.col(var2) += 2.0 * delta * observations_double_.col(var1);
}


void OMRFModel::set_pairwise_effects(const arma::mat& pairwise_effects) {
    pairwise_effects_ = pairwise_effects;
    update_residual_matrix();
    invalidate_gradient_cache();
}


// =============================================================================
// BaseModel interface implementation
// =============================================================================

size_t OMRFModel::parameter_dimension() const {
    // Count active parameters: main effects + included pairwise effects
    size_t active = num_main_;
    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            if (edge_indicators_(v1, v2) == 1) {
                active++;
            }
        }
    }
    return active;
}


void OMRFModel::set_seed(int seed) {
    rng_ = SafeRNG(seed);
}


std::unique_ptr<BaseModel> OMRFModel::clone() const {
    return std::make_unique<OMRFModel>(*this);
}


void OMRFModel::init_metropolis_adaptation(const WarmupSchedule& schedule) {
    metropolis_main_adapter_ = std::make_unique<MetropolisAdaptationController>(
        proposal_sd_main_, schedule, target_accept_);
    metropolis_pairwise_adapter_ = std::make_unique<MetropolisAdaptationController>(
        proposal_sd_pairwise_, schedule, target_accept_);
}


void OMRFModel::tune_proposal_sd(int iteration, const WarmupSchedule& schedule) {
    auto rm_weight_opt = schedule.rm_weight_for_proposal_sd(iteration);
    if (!rm_weight_opt) return;
    const double rm_weight = *rm_weight_opt;
    const double target_accept = target_accept_;

    const int num_variables = static_cast<int>(p_);

    for (int variable1 = 0; variable1 < num_variables - 1; variable1++) {
        for (int variable2 = variable1 + 1; variable2 < num_variables; variable2++) {
            double current = pairwise_effects_(variable1, variable2);
            double proposal_sd = proposal_sd_pairwise_(variable1, variable2);

            auto log_post = [&](double theta) {
                double delta = theta - current;
                return log_pseudoposterior_pairwise_at_delta(variable1, variable2, delta);
            };

            StepResult result = metropolis_step(current, proposal_sd, log_post, rng_);

            double value = result.state[0];
            pairwise_effects_(variable1, variable2) = value;
            pairwise_effects_(variable2, variable1) = value;

            if (current != value) {
                double delta = value - current;
                residual_matrix_.col(variable1) += 2.0 * observations_double_.col(variable2) * delta;
                residual_matrix_.col(variable2) += 2.0 * observations_double_.col(variable1) * delta;
            }

            proposal_sd = update_proposal_sd_with_robbins_monro(
                proposal_sd, MY_LOG(result.accept_prob), rm_weight, target_accept);
            proposal_sd_pairwise_(variable1, variable2) = proposal_sd;
            proposal_sd_pairwise_(variable2, variable1) = proposal_sd;
        }
    }

    invalidate_gradient_cache();
}


// =============================================================================
// Parameter vectorization
// =============================================================================

arma::vec OMRFModel::vectorize_parameters() const {
    // Count active parameters
    int num_active = 0;
    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            if (edge_indicators_(v1, v2) == 1) {
                num_active++;
            }
        }
    }

    arma::vec param_vec(num_main_ + num_active);
    int offset = 0;

    // Main effects
    for (size_t v = 0; v < p_; ++v) {
        if (is_ordinal_variable_(v)) {
            int num_cats = num_categories_(v);
            for (int c = 0; c < num_cats; ++c) {
                param_vec(offset++) = main_effects_(v, c);
            }
        } else {
            param_vec(offset++) = main_effects_(v, 0);  // linear
            param_vec(offset++) = main_effects_(v, 1);  // quadratic
        }
    }

    // Active pairwise effects
    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            if (edge_indicators_(v1, v2) == 1) {
                param_vec(offset++) = pairwise_effects_(v1, v2);
            }
        }
    }

    return param_vec;
}


void OMRFModel::unvectorize_parameters(const arma::vec& param_vec) {
    int offset = 0;

    // Main effects
    for (size_t v = 0; v < p_; ++v) {
        if (is_ordinal_variable_(v)) {
            int num_cats = num_categories_(v);
            for (int c = 0; c < num_cats; ++c) {
                main_effects_(v, c) = param_vec(offset++);
            }
        } else {
            main_effects_(v, 0) = param_vec(offset++);  // linear
            main_effects_(v, 1) = param_vec(offset++);  // quadratic
        }
    }

    // Active pairwise effects
    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            if (edge_indicators_(v1, v2) == 1) {
                double val = param_vec(offset++);
                pairwise_effects_(v1, v2) = val;
                pairwise_effects_(v2, v1) = val;
            }
        }
    }

    update_residual_matrix();
    invalidate_gradient_cache();
}


void OMRFModel::unvectorize_to_temps(
    const arma::vec& parameters,
    arma::mat& temp_main,
    arma::mat& temp_pairwise,
    arma::mat& temp_residual
) const {
    int offset = 0;
    for (size_t v = 0; v < p_; ++v) {
        if (is_ordinal_variable_(v)) {
            int num_cats = num_categories_(v);
            for (int c = 0; c < num_cats; ++c) {
                temp_main(v, c) = parameters(offset++);
            }
        } else {
            temp_main(v, 0) = parameters(offset++);
            temp_main(v, 1) = parameters(offset++);
        }
    }

    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            if (edge_indicators_(v1, v2) == 1) {
                temp_pairwise(v1, v2) = parameters(offset++);
                temp_pairwise(v2, v1) = temp_pairwise(v1, v2);
            }
        }
    }

    temp_residual = 2.0 * observations_double_ * temp_pairwise;
}


arma::vec OMRFModel::get_vectorized_parameters() const {
    return vectorize_parameters();
}


void OMRFModel::set_vectorized_parameters(const arma::vec& parameters) {
    unvectorize_parameters(parameters);
}


arma::vec OMRFModel::get_full_vectorized_parameters() const {
    // Fixed-size vector: all main effects + ALL pairwise effects
    arma::vec param_vec(num_main_ + num_pairwise_);
    int offset = 0;

    // Main effects
    for (size_t v = 0; v < p_; ++v) {
        if (is_ordinal_variable_(v)) {
            int num_cats = num_categories_(v);
            for (int c = 0; c < num_cats; ++c) {
                param_vec(offset++) = main_effects_(v, c);
            }
        } else {
            param_vec(offset++) = main_effects_(v, 0);  // linear
            param_vec(offset++) = main_effects_(v, 1);  // quadratic
        }
    }

    // ALL pairwise effects (zeros for inactive edges)
    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            param_vec(offset++) = pairwise_effects_(v1, v2);
        }
    }

    return param_vec;
}


arma::ivec OMRFModel::get_vectorized_indicator_parameters() {
    arma::ivec indicators(num_pairwise_);
    int idx = 0;
    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            indicators(idx++) = edge_indicators_(v1, v2);
        }
    }
    return indicators;
}


arma::vec OMRFModel::get_active_inv_mass() const {
    if (!edge_selection_active_) {
        return inv_mass_;
    }

    // Count active parameters
    int num_active = 0;
    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            if (edge_indicators_(v1, v2) == 1) {
                num_active++;
            }
        }
    }

    arma::vec active_inv_mass(num_main_ + num_active);
    active_inv_mass.head(num_main_) = inv_mass_.head(num_main_);

    int offset_full = num_main_;
    int offset_active = num_main_;

    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            if (edge_indicators_(v1, v2) == 1) {
                active_inv_mass(offset_active) = inv_mass_(offset_full);
                offset_active++;
            }
            offset_full++;
        }
    }

    return active_inv_mass;
}


void OMRFModel::vectorize_parameters_into(arma::vec& param_vec) const {
    // Count active parameters
    int num_active = 0;
    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            if (edge_indicators_(v1, v2) == 1) {
                num_active++;
            }
        }
    }

    // Resize if needed (should rarely happen after first call)
    size_t needed_size = num_main_ + num_active;
    if (param_vec.n_elem != needed_size) {
        param_vec.set_size(needed_size);
    }

    int offset = 0;

    // Main effects
    for (size_t v = 0; v < p_; ++v) {
        if (is_ordinal_variable_(v)) {
            int num_cats = num_categories_(v);
            for (int c = 0; c < num_cats; ++c) {
                param_vec(offset++) = main_effects_(v, c);
            }
        } else {
            param_vec(offset++) = main_effects_(v, 0);  // linear
            param_vec(offset++) = main_effects_(v, 1);  // quadratic
        }
    }

    // Active pairwise effects
    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            if (edge_indicators_(v1, v2) == 1) {
                param_vec(offset++) = pairwise_effects_(v1, v2);
            }
        }
    }
}


void OMRFModel::get_active_inv_mass_into(arma::vec& active_inv_mass) const {
    if (!edge_selection_active_) {
        // No edge selection - just use full inv_mass
        if (active_inv_mass.n_elem != inv_mass_.n_elem) {
            active_inv_mass.set_size(inv_mass_.n_elem);
        }
        active_inv_mass = inv_mass_;
        return;
    }

    // Count active parameters
    int num_active = 0;
    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            if (edge_indicators_(v1, v2) == 1) {
                num_active++;
            }
        }
    }

    size_t needed_size = num_main_ + num_active;
    if (active_inv_mass.n_elem != needed_size) {
        active_inv_mass.set_size(needed_size);
    }

    active_inv_mass.head(num_main_) = inv_mass_.head(num_main_);

    int offset_full = num_main_;
    int offset_active = num_main_;

    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            if (edge_indicators_(v1, v2) == 1) {
                active_inv_mass(offset_active) = inv_mass_(offset_full);
                offset_active++;
            }
            offset_full++;
        }
    }
}


// =============================================================================
// Log-pseudoposterior computation
// =============================================================================

double OMRFModel::log_pseudoposterior_main_component(int variable, int category, int parameter) const {
    double log_posterior = 0.0;

    const int num_cats = num_categories_(variable);
    arma::vec bound = num_cats * residual_matrix_.col(variable);

    if (is_ordinal_variable_(variable)) {
        const double value = main_effects_(variable, category);
        log_posterior += value * counts_per_category_(category + 1, variable);
        log_posterior += threshold_prior_->logp(value);

        arma::vec residual_score = residual_matrix_.col(variable);
        arma::vec main_effect_param = main_effects_.row(variable).cols(0, num_cats - 1).t();

        arma::vec denom = compute_denom_ordinal(residual_score, main_effect_param, bound);
        log_posterior -= arma::accu(bound + ARMA_MY_LOG(denom));
    } else {
        const double value = main_effects_(variable, parameter);
        const double linear_main_effect = main_effects_(variable, 0);
        const double quadratic_main_effect = main_effects_(variable, 1);
        const int ref = baseline_category_(variable);

        log_posterior += value * blume_capel_stats_(parameter, variable);
        log_posterior += threshold_prior_->logp(value);

        arma::vec residual_score = residual_matrix_.col(variable);
        arma::vec denom(n_, arma::fill::zeros);

        denom = compute_denom_blume_capel(
            residual_score, linear_main_effect, quadratic_main_effect, ref, num_cats, bound
        );

        log_posterior -= arma::accu(bound + ARMA_MY_LOG(denom));
    }

    return log_posterior;
}


double OMRFModel::compute_log_likelihood_ratio_for_variable(
    int variable,
    const arma::ivec& interacting_score,
    double proposed_state,
    double current_state
) const {
    // Convert interaction score vector to double precision
    arma::vec interaction = arma::conv_to<arma::vec>::from(interacting_score);

    const int num_persons = static_cast<int>(n_);
    const int num_cats = num_categories_(variable);

    // Compute adjusted linear predictors without the current interaction
    arma::vec residual_score = residual_matrix_.col(variable) - 2.0 * interaction * current_state;
    arma::vec bounds = residual_score * num_cats;

    arma::vec denom_current = arma::zeros(num_persons);
    arma::vec denom_proposed = arma::zeros(num_persons);

    if (is_ordinal_variable_(variable)) {
        arma::vec main_param = main_effects_.row(variable).cols(0, num_cats - 1).t();

        denom_current += compute_denom_ordinal(
            residual_score + 2.0 * interaction * current_state, main_param, bounds
        );
        denom_proposed += compute_denom_ordinal(
            residual_score + 2.0 * interaction * proposed_state, main_param, bounds
        );
    } else {
        const int ref_cat = baseline_category_(variable);

        denom_current = compute_denom_blume_capel(
            residual_score + 2.0 * interaction * current_state, main_effects_(variable, 0),
            main_effects_(variable, 1), ref_cat, num_cats, bounds
        );
        double log_ratio = arma::accu(ARMA_MY_LOG(denom_current) + bounds);

        denom_proposed = compute_denom_blume_capel(
            residual_score + 2.0 * interaction * proposed_state, main_effects_(variable, 0),
            main_effects_(variable, 1), ref_cat, num_cats, bounds
        );
        log_ratio -= arma::accu(ARMA_MY_LOG(denom_proposed) + bounds);

        return log_ratio;
    }

    // Accumulated log-likelihood difference across persons
    return arma::accu(ARMA_MY_LOG(denom_current) - ARMA_MY_LOG(denom_proposed));
}


double OMRFModel::log_pseudolikelihood_ratio_interaction(
    int variable1,
    int variable2,
    double proposed_state,
    double current_state
) const {
    double log_ratio = 0.0;
    const double delta = proposed_state - current_state;

    arma::ivec score1 = observations_.col(variable1);
    arma::ivec score2 = observations_.col(variable2);

    log_ratio += 4.0 * pairwise_stats_(variable1, variable2) * delta;

    log_ratio += compute_log_likelihood_ratio_for_variable(
        variable1, score2, proposed_state, current_state
    );

    log_ratio += compute_log_likelihood_ratio_for_variable(
        variable2, score1, proposed_state, current_state
    );

    return log_ratio;
}


double OMRFModel::log_pseudoposterior_pairwise_at_delta(int var1, int var2, double delta) const {
    const int num_observations = static_cast<int>(n_);
    const double proposed_value = pairwise_effects_(var1, var2) + delta;

    double log_pseudo_posterior = 4.0 * proposed_value * pairwise_stats_(var1, var2);

    const arma::vec& obs_var1 = observations_double_.col(var1);
    const arma::vec& obs_var2 = observations_double_.col(var2);

    for (int var : {var1, var2}) {
        int num_cats = num_categories_(var);
        const arma::vec& obs_other = (var == var1) ? obs_var2 : obs_var1;

        arma::vec residual_score = residual_matrix_.col(var) + 2.0 * obs_other * delta;
        arma::vec denominator = arma::zeros(num_observations);
        arma::vec bound = num_cats * residual_score;

        if (is_ordinal_variable_(var)) {
            arma::vec main_effect_param = main_effects_.row(var).cols(0, num_cats - 1).t();
            denominator += compute_denom_ordinal(residual_score, main_effect_param, bound);
        } else {
            const int ref = baseline_category_(var);
            denominator = compute_denom_blume_capel(
                residual_score, main_effects_(var, 0), main_effects_(var, 1), ref, num_cats, bound
            );
        }

        log_pseudo_posterior -= arma::accu(ARMA_MY_LOG(denominator));
        log_pseudo_posterior -= arma::accu(bound);
    }

    if (edge_indicators_(var1, var2) == 1) {
        log_pseudo_posterior += interaction_prior_->logp(proposed_value);
    }

    return log_pseudo_posterior;
}


// =============================================================================
// Gradient computation
// =============================================================================

void OMRFModel::ensure_gradient_cache() {
    if (gradient_cache_valid_) return;

    const int num_variables = static_cast<int>(p_);
    const int num_main = static_cast<int>(num_main_);
    index_matrix_cache_.set_size(num_variables, num_variables);
    index_matrix_cache_.zeros();

    // Count active pairwise effects + build index map
    int num_active = 0;
    for (int i = 0; i < num_variables - 1; i++) {
        for (int j = i + 1; j < num_variables; j++) {
            if (edge_indicators_(i, j) == 1) {
                index_matrix_cache_(i, j) = num_main + num_active++;
                index_matrix_cache_(j, i) = index_matrix_cache_(i, j);
            }
        }
    }

    // Allocate gradient vector (main + active pairwise only)
    grad_obs_cache_.set_size(num_main + num_active);
    grad_obs_cache_.zeros();

    // Observed statistics for main effects
    int offset = 0;
    for (int variable = 0; variable < num_variables; variable++) {
        if (is_ordinal_variable_(variable)) {
            const int num_cats = num_categories_(variable);
            for (int cat = 0; cat < num_cats; cat++) {
                grad_obs_cache_(offset + cat) = counts_per_category_(cat + 1, variable);
            }
            offset += num_cats;
        } else {
            grad_obs_cache_(offset) = blume_capel_stats_(0, variable);
            grad_obs_cache_(offset + 1) = blume_capel_stats_(1, variable);
            offset += 2;
        }
    }

    // Observed statistics for pairwise effects
    for (int i = 0; i < num_variables - 1; i++) {
        for (int j = i + 1; j < num_variables; j++) {
            if (edge_indicators_(i, j) == 0) continue;
            int location = index_matrix_cache_(i, j);
            grad_obs_cache_(location) = 4.0 * pairwise_stats_(i, j);
        }
    }

    gradient_cache_valid_ = true;
}

std::pair<double, arma::vec> OMRFModel::logp_and_gradient(const arma::vec& parameters) {
    ensure_gradient_cache();

    arma::mat temp_main(main_effects_.n_rows, main_effects_.n_cols, arma::fill::none);
    arma::mat temp_pairwise(p_, p_, arma::fill::zeros);
    arma::mat temp_residual;
    unvectorize_to_temps(parameters, temp_main, temp_pairwise, temp_residual);

    const int num_variables = static_cast<int>(p_);

    double log_pp = 0.0;
    arma::vec gradient = grad_obs_cache_;

    // ---- Main effects: priors + sufficient statistics ----
    for (int variable = 0; variable < num_variables; variable++) {
        if (is_ordinal_variable_(variable)) {
            const int num_cats = num_categories_(variable);
            for (int cat = 0; cat < num_cats; cat++) {
                double value = temp_main(variable, cat);
                log_pp += counts_per_category_(cat + 1, variable) * value;
                log_pp += threshold_prior_->logp(value);
            }
        } else {
            double value = temp_main(variable, 0);
            log_pp += threshold_prior_->logp(value);
            log_pp += blume_capel_stats_(0, variable) * value;

            value = temp_main(variable, 1);
            log_pp += threshold_prior_->logp(value);
            log_pp += blume_capel_stats_(1, variable) * value;
        }
    }

    // ---- Pairwise effects: priors + sufficient statistics ----
    for (int var1 = 0; var1 < num_variables - 1; var1++) {
        for (int var2 = var1 + 1; var2 < num_variables; var2++) {
            if (edge_indicators_(var1, var2) == 0) continue;

            double value = temp_pairwise(var1, var2);
            log_pp += 4.0 * pairwise_stats_(var1, var2) * value;
            log_pp += interaction_prior_->logp(value);
        }
    }

    // ---- Per-variable: joint computation of log-normalizer and gradient ----
    int offset = 0;
    for (int variable = 0; variable < num_variables; variable++) {
        const int num_cats = num_categories_(variable);
        arma::vec residual_score = temp_residual.col(variable);
        arma::vec bound = num_cats * residual_score;

        if (is_ordinal_variable_(variable)) {
            arma::vec main_param = temp_main.row(variable).cols(0, num_cats - 1).t();

            // Fill-in-place: persistent per-chain scratch reused across variables
            // and iterations, eliminating per-call heap allocations.
            compute_logZ_and_probs_ordinal_into(
                main_param, residual_score, bound, num_cats,
                logz_out_, logz_scratch_
            );

            // Use log_Z for log-pseudoposterior
            log_pp -= arma::accu(logz_out_.log_Z);

            // Use probs for gradient
            for (int cat = 0; cat < num_cats; cat++) {
                gradient(offset + cat) -= arma::accu(logz_out_.probs.col(cat + 1));
            }

            // Pairwise gradient contributions (vectorized using BLAS)
            arma::vec weights = arma::regspace<arma::vec>(1, num_cats);
            arma::vec E = logz_out_.probs.cols(1, num_cats) * weights;
            arma::vec pw_grad = observations_double_t_ * E;
            for (int j = 0; j < num_variables; j++) {
                if (edge_indicators_(variable, j) == 0 || variable == j) continue;
                int location = (variable < j) ? index_matrix_cache_(variable, j) : index_matrix_cache_(j, variable);
                gradient(location) -= 2.0 * pw_grad(j);
            }
            offset += num_cats;
        } else {
            const int ref = baseline_category_(variable);
            const double lin_eff = temp_main(variable, 0);
            const double quad_eff = temp_main(variable, 1);

            // Fill-in-place Blume-Capel
            compute_logZ_and_probs_blume_capel_into(
                residual_score, lin_eff, quad_eff, ref, num_cats, bound,
                logz_out_, logz_scratch_
            );

            // Use log_Z for log-pseudoposterior
            log_pp -= arma::accu(logz_out_.log_Z);

            // Use probs for gradient
            arma::vec score = arma::regspace<arma::vec>(0, num_cats) - static_cast<double>(ref);
            arma::vec sq_score = arma::square(score);

            gradient(offset)     -= arma::accu(logz_out_.probs * score);
            gradient(offset + 1) -= arma::accu(logz_out_.probs * sq_score);

            // Pairwise gradient contributions (vectorized using BLAS)
            arma::vec E = logz_out_.probs * score;
            arma::vec pw_grad = observations_double_t_ * E;
            for (int j = 0; j < num_variables; j++) {
                if (edge_indicators_(variable, j) == 0 || variable == j) continue;
                int location = (variable < j) ? index_matrix_cache_(variable, j) : index_matrix_cache_(j, variable);
                gradient(location) -= 2.0 * pw_grad(j);
            }
            offset += 2;
        }
    }

    // ---- Priors: gradient contributions ----
    offset = 0;
    for (int variable = 0; variable < num_variables; variable++) {
        if (is_ordinal_variable_(variable)) {
            const int num_cats = num_categories_(variable);
            for (int cat = 0; cat < num_cats; cat++) {
                gradient(offset + cat) += threshold_prior_->grad(temp_main(variable, cat));
            }
            offset += num_cats;
        } else {
            for (int k = 0; k < 2; k++) {
                gradient(offset + k) += threshold_prior_->grad(temp_main(variable, k));
            }
            offset += 2;
        }
    }
    for (int i = 0; i < num_variables - 1; i++) {
        for (int j = i + 1; j < num_variables; j++) {
            if (edge_indicators_(i, j) == 0) continue;
            int location = index_matrix_cache_(i, j);
            const double effect = temp_pairwise(i, j);
            gradient(location) += interaction_prior_->grad(effect);
        }
    }

    return {log_pp, gradient};
}


// =============================================================================
// Metropolis-Hastings updates
// =============================================================================

double OMRFModel::update_main_effect_parameter(int variable, int category, int parameter) {
    double& current = is_ordinal_variable_(variable)
        ? main_effects_(variable, category)
        : main_effects_(variable, parameter);

    double proposal_sd = is_ordinal_variable_(variable)
        ? proposal_sd_main_(variable, category)
        : proposal_sd_main_(variable, parameter);

    auto log_post = [&](double theta) {
        current = theta;
        return log_pseudoposterior_main_component(variable, category, parameter);
    };

    StepResult result = metropolis_step(current, proposal_sd, log_post, rng_);
    current = result.state[0];
    return result.accept_prob;
}


double OMRFModel::update_pairwise_effect(int var1, int var2) {
    if (edge_indicators_(var1, var2) == 0) return 1.0;

    double current_value = pairwise_effects_(var1, var2);
    double proposal_sd = proposal_sd_pairwise_(var1, var2);

    auto log_post = [&](double theta) {
        double delta = theta - current_value;
        return log_pseudoposterior_pairwise_at_delta(var1, var2, delta);
    };

    StepResult result = metropolis_step(current_value, proposal_sd, log_post, rng_);

    double value = result.state[0];
    pairwise_effects_(var1, var2) = value;
    pairwise_effects_(var2, var1) = value;

    if (current_value != value) {
        double delta = value - current_value;
        residual_matrix_.col(var1) += 2.0 * observations_double_.col(var2) * delta;
        residual_matrix_.col(var2) += 2.0 * observations_double_.col(var1) * delta;
    }

    return result.accept_prob;
}


void OMRFModel::update_edge_indicator(int var1, int var2) {
    const double current_state = pairwise_effects_(var1, var2);

    const bool proposing_addition = (edge_indicators_(var1, var2) == 0);
    const double proposed_state = proposing_addition
        ? rnorm(rng_, current_state, proposal_sd_pairwise_(var1, var2))
        : 0.0;

    double log_accept = log_pseudolikelihood_ratio_interaction(
        var1, var2, proposed_state, current_state
    );

    const double inclusion_probability_ij = inclusion_probability_(var1, var2);
    const double sd = proposal_sd_pairwise_(var1, var2);

    if (proposing_addition) {
        log_accept += interaction_prior_->logp(proposed_state);
        log_accept -= R::dnorm(proposed_state, current_state, sd, true);
        log_accept += MY_LOG(inclusion_probability_ij) - MY_LOG(1.0 - inclusion_probability_ij);
    } else {
        log_accept -= interaction_prior_->logp(current_state);
        log_accept += R::dnorm(current_state, proposed_state, sd, true);
        log_accept -= MY_LOG(inclusion_probability_ij) - MY_LOG(1.0 - inclusion_probability_ij);
    }

    if (MY_LOG(runif(rng_)) < log_accept) {
        const int updated_indicator = 1 - edge_indicators_(var1, var2);
        edge_indicators_(var1, var2) = updated_indicator;
        edge_indicators_(var2, var1) = updated_indicator;

        pairwise_effects_(var1, var2) = proposed_state;
        pairwise_effects_(var2, var1) = proposed_state;

        const double delta = proposed_state - current_state;
        residual_matrix_.col(var1) += 2.0 * observations_double_.col(var2) * delta;
        residual_matrix_.col(var2) += 2.0 * observations_double_.col(var1) * delta;
    }
}


// =============================================================================
// Main update methods
// =============================================================================

void OMRFModel::do_one_metropolis_step(int iteration) {
    // Track running mean acceptance probability across all components
    // updated this iteration; exposed via last_metropolis_mean_accept_prob().
    double sum_accept = 0.0;
    int    n_accept   = 0;

    // --- Pairwise effects sweep ---
    arma::mat accept_prob_pairwise = arma::zeros<arma::mat>(p_, p_);
    arma::umat index_mask_pairwise = arma::zeros<arma::umat>(p_, p_);

    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            double ap = update_pairwise_effect(v1, v2);
            if (edge_indicators_(v1, v2) == 1) {
                accept_prob_pairwise(v1, v2) = ap;
                index_mask_pairwise(v1, v2) = 1;
                sum_accept += ap;
                ++n_accept;
            }
        }
    }

    if (metropolis_pairwise_adapter_) {
        metropolis_pairwise_adapter_->update(index_mask_pairwise, accept_prob_pairwise, iteration);
    }

    // --- Main effects sweep ---
    arma::umat index_mask_main = arma::ones<arma::umat>(
        proposal_sd_main_.n_rows, proposal_sd_main_.n_cols);
    arma::mat accept_prob_main = arma::ones<arma::mat>(
        proposal_sd_main_.n_rows, proposal_sd_main_.n_cols);

    for (size_t v = 0; v < p_; ++v) {
        if (is_ordinal_variable_(v)) {
            int num_cats = num_categories_(v);
            for (int c = 0; c < num_cats; ++c) {
                double ap = update_main_effect_parameter(v, c, -1);
                accept_prob_main(v, c) = ap;
                sum_accept += ap;
                ++n_accept;
            }
        } else {
            for (int p = 0; p < 2; ++p) {
                double ap = update_main_effect_parameter(v, -1, p);
                accept_prob_main(v, p) = ap;
                sum_accept += ap;
                ++n_accept;
            }
        }
    }

    if (metropolis_main_adapter_) {
        metropolis_main_adapter_->update(index_mask_main, accept_prob_main, iteration);
    }

    last_mh_mean_accept_ = (n_accept > 0)
        ? sum_accept / static_cast<double>(n_accept)
        : std::numeric_limits<double>::quiet_NaN();

    invalidate_gradient_cache();
}


void OMRFModel::prepare_iteration() {
    // Shuffle edge order unconditionally to advance the RNG state consistently.
    shuffled_edge_order_ = arma_randperm(rng_, num_pairwise_);
}


void OMRFModel::update_edge_indicators() {
    for (size_t i = 0; i < num_pairwise_; ++i) {
        int idx = shuffled_edge_order_(i);
        int var1 = interaction_index_(idx, 1);
        int var2 = interaction_index_(idx, 2);
        update_edge_indicator(var1, var2);
    }
}



void OMRFModel::impute_missing() {
    if (!has_missing_) return;

    const int num_variables = p_;
    const int num_missings = missing_index_.n_rows;
    const int max_num_categories = num_categories_.max();

    arma::vec category_probabilities(max_num_categories + 1);

    for (int miss = 0; miss < num_missings; miss++) {
        const int person = missing_index_(miss, 0);
        const int variable = missing_index_(miss, 1);

        const double residual_score = residual_matrix_(person, variable);
        const int num_cats = num_categories_(variable);
        const bool is_ordinal = is_ordinal_variable_(variable);

        double cumsum = 0.0;

        if (is_ordinal) {
            cumsum = 1.0;
            category_probabilities[0] = cumsum;
            for (int cat = 0; cat < num_cats; cat++) {
                const int score = cat + 1;
                const double exponent = main_effects_(variable, cat) + score * residual_score;
                cumsum += MY_EXP(exponent);
                category_probabilities[score] = cumsum;
            }
        } else {
            const int ref = baseline_category_(variable);
            cumsum = 0.0;

            for (int cat = 0; cat <= num_cats; cat++) {
                const int score = cat - ref;
                const double exponent =
                    main_effects_(variable, 0) * score +
                    main_effects_(variable, 1) * score * score +
                    score * residual_score;
                cumsum += MY_EXP(exponent);
                category_probabilities[cat] = cumsum;
            }
        }

        // Sample from categorical distribution via inverse transform
        const double u = runif(rng_) * cumsum;
        int sampled_score = 0;
        while (u > category_probabilities[sampled_score]) {
            sampled_score++;
        }

        int new_value = sampled_score;
        if (!is_ordinal)
            new_value -= baseline_category_(variable);
        const int old_value = observations_(person, variable);

        if (new_value != old_value) {
            observations_(person, variable) = new_value;
            observations_double_(person, variable) = static_cast<double>(new_value);

            if (is_ordinal) {
                counts_per_category_(old_value, variable)--;
                counts_per_category_(new_value, variable)++;
            } else {
                const int delta = new_value - old_value;
                const int delta_sq = new_value * new_value - old_value * old_value;
                blume_capel_stats_(0, variable) += delta;
                blume_capel_stats_(1, variable) += delta_sq;
            }

            // Incrementally update residuals across all variables
            for (int var = 0; var < num_variables; var++) {
                const double delta_score = 2.0 * (new_value - old_value) * pairwise_effects_(var, variable);
                residual_matrix_(person, var) += delta_score;
            }
        }
    }

    // Recompute pairwise sufficient statistics
    arma::mat ps = observations_double_.t() * observations_double_;
    pairwise_stats_ = arma::conv_to<arma::imat>::from(ps);

    // Update cached transpose so gradients use current imputed values
    observations_double_t_ = observations_double_.t();

    // Sufficient statistics changed; gradient cache must be rebuilt
    invalidate_gradient_cache();
}


void OMRFModel::set_missing_data(const arma::imat& missing_index) {
    missing_index_ = missing_index;
    has_missing_ = (missing_index.n_rows > 0 && missing_index.n_cols == 2);
}


// =============================================================================
// Factory function
// =============================================================================

OMRFModel createOMRFModelFromR(
    const Rcpp::List& inputFromR,
    const arma::mat& inclusion_probability,
    const arma::imat& initial_edge_indicators,
    std::unique_ptr<BaseParameterPrior> interaction_prior,
    std::unique_ptr<BaseParameterPrior> threshold_prior,
    bool edge_selection
) {
    arma::imat observations = Rcpp::as<arma::imat>(inputFromR["observations"]);
    arma::ivec num_categories = Rcpp::as<arma::ivec>(inputFromR["num_categories"]);
    arma::uvec is_ordinal_variable = Rcpp::as<arma::uvec>(inputFromR["is_ordinal_variable"]);
    arma::ivec baseline_category = Rcpp::as<arma::ivec>(inputFromR["baseline_category"]);

    return OMRFModel(
        observations,
        num_categories,
        inclusion_probability,
        initial_edge_indicators,
        is_ordinal_variable,
        baseline_category,
        std::move(interaction_prior),
        std::move(threshold_prior),
        edge_selection
    );
}



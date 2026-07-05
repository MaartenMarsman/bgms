#pragma once

#include <RcppArmadillo.h>
#include <memory>
#include <utility>
#include "mcmc/algorithms/hmc.h"
#include "mcmc/algorithms/leapfrog.h"
#include "mcmc/algorithms/nuts.h"
#include "mcmc/execution/sampler_config.h"
#include "mcmc/execution/step_result.h"
#include "mcmc/execution/warmup_schedule.h"
#include "mcmc/samplers/nuts_adaptation.h"
#include "mcmc/samplers/sampler_base.h"
#include "models/base_model.h"

/**
 * NUTSSampler - No-U-Turn Sampler with warmup adaptation
 *
 * Adaptive tree-depth leapfrog integration with step-size dual averaging
 * and diagonal mass-matrix estimation. Owns its NUTSAdaptationController
 * directly.
 *
 * Coordinates the multi-stage warmup:
 *  - Stage-3c boundary: edge selection just activated → restart dual averaging
 *    so adaptation can tune to the new geometry quickly.
 *  - Mass-matrix update: when the controller emits a new mass matrix, re-run
 *    the step-size heuristic with the new metric.
 *  - Phase-aware reverse check: observe during warmup, enforce during sampling
 *    (constrained integration only).
 *
 * For constrained models (edge selection or sparse graph), uses RATTLE
 * integration: full Cholesky space with position and momentum projection at
 * each leapfrog step.
 */
class NUTSSampler : public SamplerBase {
public:
    explicit NUTSSampler(const SamplerConfig& config, WarmupSchedule& schedule)
        : step_size_(config.initial_step_size),
          target_acceptance_(config.target_acceptance),
          schedule_(schedule),
          max_tree_depth_(config.max_tree_depth),
          reverse_check_(config.reverse_check),
          reverse_check_tol_(config.reverse_check_tol),
          learn_mass_matrix_(config.learn_mass_matrix),
          initialized_(false)
    {}

    bool has_nuts_diagnostics() const override { return true; }

    StepResult step(BaseModel& model, int iteration) override {
        // Stage 3c boundary: edge selection just activated.
        // Restart dual averaging so adaptation can tune to the new
        // geometry (changed active parameters) quickly.
        if (schedule_.in_stage3c(iteration) && !stage3c_initialized_) {
            stage3c_initialized_ = true;
            if (uses_constrained_integration(model)) {
                // Re-run step-size heuristic with constrained integrator.
                // The step size from stages 1-2 was tuned for unconstrained
                // leapfrog; constraints are now active, so re-tune.
                SafeRNG& rng = model.get_rng();
                arma::vec x = model.get_full_position();
                arma::vec inv_mass = model.get_inv_mass();
                auto joint_fn = [&model](const arma::vec& params)
                    -> std::pair<double, arma::vec> {
                    return model.logp_and_gradient_full(params);
                };
                ProjectPositionFn proj_pos = [&model, &inv_mass](arma::vec& pos) {
                    model.project_position(pos, inv_mass);
                };
                ProjectMomentumFn proj_mom = [&model, &inv_mass](arma::vec& mom, const arma::vec& pos) {
                    model.project_momentum(mom, pos, inv_mass);
                };
                double new_eps = heuristic_initial_step_size_constrained(
                    x, joint_fn, inv_mass, proj_pos, proj_mom, rng,
                    target_acceptance_, nuts_adapt_->current_step_size());
                nuts_adapt_->reinit_stepsize(new_eps);
            } else {
                nuts_adapt_->reinit_stepsize(nuts_adapt_->current_step_size());
            }
        }

        // Use adaptation controller's current step size for this iteration
        step_size_ = nuts_adapt_->current_step_size();

        // Phase-aware reverse check: observe during warmup, enforce during sampling.
        // The check always runs (recording non_reversible counts),
        // but only terminates trees / rejects steps when enforcing.
        enforce_reverse_check_ = schedule_.sampling(iteration);

        StepResult result = uses_constrained_integration(model)
            ? do_constrained_step(model)
            : do_unconstrained_step(model);

        // Let the adaptation controller handle step-size and mass-matrix logic.
        // For RATTLE (constrained) models, feed x-space samples so the mass
        // matrix is estimated in the same coordinate system NUTS operates in.
        arma::vec full_params = uses_constrained_integration(model)
            ? model.get_full_position()
            : model.get_full_vectorized_parameters();
        nuts_adapt_->update(full_params, result.accept_prob, iteration);

        // If mass matrix was just updated, apply it and re-run the step-size heuristic
        if (nuts_adapt_->mass_matrix_just_updated()) {
            arma::vec new_inv_mass = nuts_adapt_->inv_mass_diag();
            model.set_inv_mass(new_inv_mass);

            SafeRNG& rng = model.get_rng();

            if (uses_constrained_integration(model)) {
                arma::vec x = model.get_full_position();
                auto joint_fn = [&model](const arma::vec& params)
                    -> std::pair<double, arma::vec> {
                    return model.logp_and_gradient_full(params);
                };
                ProjectPositionFn proj_pos = [&model, &new_inv_mass](arma::vec& pos) {
                    model.project_position(pos, new_inv_mass);
                };
                ProjectMomentumFn proj_mom = [&model, &new_inv_mass](arma::vec& mom, const arma::vec& pos) {
                    model.project_momentum(mom, pos, new_inv_mass);
                };
                double new_eps = heuristic_initial_step_size_constrained(
                    x, joint_fn, new_inv_mass, proj_pos, proj_mom, rng,
                    target_acceptance_, nuts_adapt_->current_step_size());
                nuts_adapt_->reinit_stepsize(new_eps);
            } else {
                arma::vec theta = model.get_vectorized_parameters();
                auto grad_fn = [&model](const arma::vec& params) -> arma::vec {
                    return model.logp_and_gradient(params).second;
                };
                auto joint_fn = [&model](const arma::vec& params)
                    -> std::pair<double, arma::vec> {
                    return model.logp_and_gradient(params);
                };
                arma::vec active_inv_mass = model.get_active_inv_mass();
                double new_eps = heuristic_initial_step_size(
                    theta, grad_fn, joint_fn, active_inv_mass, rng,
                    target_acceptance_, nuts_adapt_->current_step_size());
                nuts_adapt_->reinit_stepsize(new_eps);
            }
        }

        // Update step_size_ from controller (may have changed due to mass update)
        step_size_ = nuts_adapt_->current_step_size();

        return result;
    }

    void initialize(BaseModel& model) override {
        if (initialized_) return;
        do_initialize(model);
        initialized_ = true;
    }

    double get_step_size() const { return step_size_; }
    double get_averaged_step_size() const {
        return nuts_adapt_ ? nuts_adapt_->final_step_size() : step_size_;
    }
    const arma::vec& get_inv_mass() const { return nuts_adapt_->inv_mass_diag(); }

private:
    bool uses_constrained_integration(const BaseModel& model) const {
        return model.has_constraints();
    }

    StepResult do_unconstrained_step(BaseModel& model) {
        arma::vec theta = model.get_vectorized_parameters();
        SafeRNG& rng = model.get_rng();

        auto joint_fn = [&model](const arma::vec& params)
            -> std::pair<double, arma::vec> {
            return model.logp_and_gradient(params);
        };

        arma::vec active_inv_mass = model.get_active_inv_mass();

        StepResult result = nuts_step(
            theta, step_size_, joint_fn,
            active_inv_mass, rng, max_tree_depth_
        );

        model.set_vectorized_parameters(result.state);
        return result;
    }

    StepResult do_constrained_step(BaseModel& model) {
        model.reset_projection_cache();
        arma::vec x = model.get_full_position();
        SafeRNG& rng = model.get_rng();

        auto joint_fn = [&model](const arma::vec& params)
            -> std::pair<double, arma::vec> {
            return model.logp_and_gradient_full(params);
        };

        arma::vec inv_mass = model.get_inv_mass();

        ProjectPositionFn proj_pos = [&model, &inv_mass](arma::vec& pos) {
            model.project_position(pos, inv_mass);
        };

        ProjectMomentumFn proj_mom = [&model, &inv_mass](arma::vec& mom, const arma::vec& pos) {
            model.project_momentum(mom, pos, inv_mass);
        };

        StepResult result = nuts_step(
            x, step_size_, joint_fn,
            inv_mass, rng, max_tree_depth_,
            &proj_pos, &proj_mom,
            reverse_check_ && enforce_reverse_check_,
            reverse_check_tol_
        );

        model.set_full_position(result.state);
        return result;
    }

    void do_initialize(BaseModel& model) {
        int dim = static_cast<int>(model.full_parameter_dimension());
        SafeRNG& rng = model.get_rng();

        // Initialize inverse mass to ones
        arma::vec init_inv_mass = arma::ones<arma::vec>(dim);
        model.set_inv_mass(init_inv_mass);

        double init_eps;

        if (uses_constrained_integration(model)) {
            // Project initial position onto constraint manifold before
            // computing step size. The MLE initialization may violate
            // K_ij = 0 constraints for excluded edges.
            arma::vec x = model.get_full_position();
            arma::vec r_dummy = arma::zeros<arma::vec>(x.n_elem);
            model.project_position(x);
            model.project_momentum(r_dummy, x);
            model.set_full_position(x);

            x = model.get_full_position();
            auto grad_fn = [&model](const arma::vec& params) -> arma::vec {
                return model.logp_and_gradient_full(params).second;
            };
            auto joint_fn = [&model](const arma::vec& params)
                -> std::pair<double, arma::vec> {
                return model.logp_and_gradient_full(params);
            };
            init_eps = heuristic_initial_step_size(
                x, grad_fn, joint_fn, rng, target_acceptance_);
        } else {
            arma::vec theta = model.get_vectorized_parameters();
            auto grad_fn = [&model](const arma::vec& params) -> arma::vec {
                return model.logp_and_gradient(params).second;
            };
            auto joint_fn = [&model](const arma::vec& params)
                -> std::pair<double, arma::vec> {
                return model.logp_and_gradient(params);
            };
            init_eps = heuristic_initial_step_size(
                theta, grad_fn, joint_fn, rng, target_acceptance_);
        }

        step_size_ = init_eps;

        // Construct the adaptation controller with the shared schedule
        nuts_adapt_ = std::make_unique<NUTSAdaptationController>(
            dim, init_eps, target_acceptance_, schedule_,
            learn_mass_matrix_);
    }

    // --- Configuration / state ---
    double step_size_;
    double target_acceptance_;
    WarmupSchedule& schedule_;
    int max_tree_depth_;
    bool reverse_check_;
    double reverse_check_tol_;
    bool learn_mass_matrix_;

    // --- Lifecycle flags ---
    bool initialized_;
    bool stage3c_initialized_ = false;

    /// Whether the reverse check should enforce (reject) this iteration.
    /// Set in step() before do_*_step(). Read by the constrained path.
    bool enforce_reverse_check_ = false;

    // --- Adaptation controller (owns step size + mass matrix) ---
    std::unique_ptr<NUTSAdaptationController> nuts_adapt_;
};

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
 *
 * Integration always runs in the active (theta-space) parameterization:
 * graph constraints are enforced by the models' null-space coordinates,
 * so edge selection and sparse graphs need no projected integrator.
 */
class NUTSSampler : public SamplerBase {
public:
    explicit NUTSSampler(const SamplerConfig& config, WarmupSchedule& schedule)
        : step_size_(config.initial_step_size),
          target_acceptance_(config.target_acceptance),
          schedule_(schedule),
          max_tree_depth_(config.max_tree_depth),
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
            nuts_adapt_->reinit_stepsize(nuts_adapt_->current_step_size());
        }

        // Use adaptation controller's current step size for this iteration
        step_size_ = nuts_adapt_->current_step_size();

        StepResult result = do_step(model);

        // Let the adaptation controller handle step-size and mass-matrix
        // logic. The mass matrix is estimated on the full (zero-padded)
        // theta layout so entries keep their slots across active-set changes.
        nuts_adapt_->update(model.get_full_vectorized_parameters(),
                            result.accept_prob, iteration);

        // If mass matrix was just updated, apply it and re-run the step-size heuristic
        if (nuts_adapt_->mass_matrix_just_updated()) {
            arma::vec new_inv_mass = nuts_adapt_->inv_mass_diag();
            model.set_inv_mass(new_inv_mass);

            SafeRNG& rng = model.get_rng();

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
    StepResult do_step(BaseModel& model) {
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

    void do_initialize(BaseModel& model) {
        int dim = static_cast<int>(model.full_parameter_dimension());
        SafeRNG& rng = model.get_rng();

        // Initialize inverse mass to ones
        arma::vec init_inv_mass = arma::ones<arma::vec>(dim);
        model.set_inv_mass(init_inv_mass);

        arma::vec theta = model.get_vectorized_parameters();
        auto grad_fn = [&model](const arma::vec& params) -> arma::vec {
            return model.logp_and_gradient(params).second;
        };
        auto joint_fn = [&model](const arma::vec& params)
            -> std::pair<double, arma::vec> {
            return model.logp_and_gradient(params);
        };
        double init_eps = heuristic_initial_step_size(
            theta, grad_fn, joint_fn, rng, target_acceptance_);

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
    bool learn_mass_matrix_;

    // --- Lifecycle flags ---
    bool initialized_;
    bool stage3c_initialized_ = false;

    // --- Adaptation controller (owns step size + mass matrix) ---
    std::unique_ptr<NUTSAdaptationController> nuts_adapt_;
};

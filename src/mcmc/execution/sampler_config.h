#pragma once

#include <string>

/**
 * SamplerConfig - Configuration for MCMC sampling
 *
 * Holds all settings for the generic MCMC runner, including:
 * - Sampler type selection (NUTS or adaptive-metropolis)
 * - Iteration counts
 * - NUTS-specific parameters
 * - Edge selection settings
 */
struct SamplerConfig {
    /// Sampler type: "nuts" or "adaptive-metropolis".
    std::string sampler_type = "adaptive-metropolis";

    /// Number of post-warmup iterations.
    int no_iter = 1000;
    /// Number of warmup iterations.
    int no_warmup = 500;

    /// Maximum NUTS tree depth.
    int max_tree_depth = 10;
    /// Initial step size for gradient-based samplers.
    double initial_step_size = 0.1;
    /// Target acceptance rate for dual-averaging adaptation.
    double target_acceptance = 0.8;

    /// Adapt the (diagonal) mass matrix during warmup.
    bool learn_mass_matrix = true;

    /// Enable spike-and-slab edge selection.
    bool edge_selection = false;

    /// Enable missing-data imputation during sampling.
    bool na_impute = false;

    /// Enable runtime reversibility check for constrained integration.
    bool reverse_check = true;

    /// Factor for the eps²-scaled reversibility tolerance (tol = factor * eps²).
    double reverse_check_tol = 0.5;

    /// Random seed.
    int seed = 42;

    /// Default constructor.
    SamplerConfig() = default;
};

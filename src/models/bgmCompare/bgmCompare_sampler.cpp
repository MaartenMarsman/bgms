#include <RcppArmadillo.h>
#include "models/bgmCompare/bgmCompare_helper.h"
#include "models/bgmCompare/bgmCompare_logp_and_grad.h"
#include "models/bgmCompare/bgmCompare_sampler.h"
#include "models/bgmCompare/bgmCompare_output.h"
#include "mcmc/samplers/metropolis_adaptation.h"
#include "mcmc/samplers/nuts_adaptation.h"
#include "mcmc/algorithms/hmc.h"
#include "mcmc/algorithms/leapfrog.h"
#include "mcmc/algorithms/nuts.h"
#include "mcmc/algorithms/metropolis.h"
#include "rng/rng_utils.h"
#include "math/explog_macros.h"
#include <string>
#include "utils/progress_manager.h"
#include "utils/common_helpers.h"
#include "priors/parameter_prior.h"



// Imputes missing observations for the bgmCompare model.
//
// This function performs single imputation of missing values during Gibbs sampling.
// Each missing entry is resampled from its conditional distribution given:
//   - the current main and pairwise effect parameters,
//   - the observed data for that individual,
//   - group-specific sufficient statistics.
//
// Workflow:
//  1. For each missing entry, identify its (person, variable, group).
//  2. Compute group-specific main and pairwise effects via projections.
//  3. Calculate unnormalized probabilities for all categories of the variable:
//     - Ordinal: softmax using category-specific thresholds.
//     - Blume–Capel: quadratic + linear score with baseline centering.
//  4. Sample a new category with inverse transform sampling.
//  5. If the imputed value differs from the old one, update:
//       - `observations` (raw data matrix),
//       - `counts_per_category` or `blume_capel_stats` (main-effect sufficient stats),
//       - `pairwise_stats` (pairwise sufficient stats).
//
// Inputs:
//  - main_effects, pairwise_effects: Current parameter matrices.
//  - main_effect_indices, pairwise_effect_indices: Lookup tables for variable/pair rows.
//  - inclusion_indicator: Indicates which differences/pairs are included.
//  - projection: Group projection matrix.
//  - observations: Data matrix [persons × variables]; updated in place.
//  - num_groups: Number of groups.
//  - group_membership: Group assignment for each person.
//  - group_indices: Row ranges [start,end] for each group.
//  - counts_per_category: Group-level sufficient statistics for ordinal variables.
//  - blume_capel_stats: Group-level sufficient statistics for Blume–Capel variables.
//  - pairwise_stats: Group-level sufficient statistics for pairwise interactions.
//  - num_categories: Number of categories for each variable in each group.
//  - missing_data_indices: Matrix of (person, variable) pairs with missing values.
//  - is_ordinal_variable: Indicator vector (1 = ordinal, 0 = Blume–Capel).
//  - baseline_category: Reference categories for Blume–Capel variables.
//  - rng: Random number generator.
//
// Notes:
//  - The function updates both raw data and sufficient statistics in-place.
//  - Group-specific pairwise effects are recomputed per missing entry.
//  - For efficiency, you may consider incremental updates to `pairwise_stats`
//    instead of full recomputation (`obs.t() * obs`) after each change.
void impute_missing_bgmcompare(
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::imat& inclusion_indicator,
    const arma::mat& projection,
    arma::imat& observations,
    const int num_groups,
    const arma::ivec& group_membership,
    const arma::imat& group_indices,
    std::vector<arma::imat>& counts_per_category,
    std::vector<arma::imat>& blume_capel_stats,
    std::vector<arma::mat>& pairwise_stats,
    const arma::ivec& num_categories,
    const arma::imat& missing_data_indices,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    SafeRNG& rng
) {
  const int num_variables = observations.n_cols;
  const int num_missings = missing_data_indices.n_rows;
  const int max_num_categories = arma::max(num_categories);

  arma::vec category_response_probabilities(max_num_categories + 1);
  double exponent, cumsum, u;
  int score, person, variable, new_value, old_value, group;

  //Impute missing data
  for(int missing = 0; missing < num_missings; missing++) {
    // Identify the observation to impute
    person = missing_data_indices(missing, 0);
    variable = missing_data_indices(missing, 1);
    group = group_membership[person];

    const arma::vec proj_g = projection.row(group).t();
    // Compute thresholds for the variable in the given group
    arma::vec group_main_effects = compute_group_main_effects(
      variable, num_groups, main_effects,  main_effect_indices, proj_g);

    // Generate a new observation based on the model
    arma::mat group_pairwise_effects(num_variables, num_variables, arma::fill::zeros);
    for(int v1 = 0; v1 < num_variables-1; v1++) {
      for(int v2 = v1 + 1; v2 < num_variables; v2++) {
        double w = compute_group_pairwise_effects(
            v1, v2, num_groups, pairwise_effects, pairwise_effect_indices,
            inclusion_indicator, proj_g
        );
        group_pairwise_effects(v1, v2) = w;
        group_pairwise_effects(v2, v1) = w;
      }
    }

    double rest_score =
      arma::as_scalar(observations.row(person) * group_pairwise_effects.col(variable));

    if(is_ordinal_variable[variable] == true) {
      // For regular binary or ordinal variables
      cumsum = 1.0;
      category_response_probabilities[0] = 1.0;
      for(int category = 1; category <= num_categories(variable); category++) {
        exponent = group_main_effects(category - 1);
        exponent += category * rest_score;
        cumsum += MY_EXP(exponent);
        category_response_probabilities[category] = cumsum;
      }
    } else {
      // For Blume-Capel variables
      cumsum = 0.0;
      const int ref = baseline_category[variable];
      for(int category = 0; category <= num_categories(variable); category++) {
        score = category - ref;
        exponent = group_main_effects[0] * score;
        exponent += group_main_effects[1] * score * score;
        exponent += rest_score * score;
        cumsum += MY_EXP(exponent);
        category_response_probabilities[category] = cumsum;
      }
    }

    // Sample a new value based on computed probabilities
    u = cumsum * runif(rng);
    score = 0;
    while (u > category_response_probabilities[score]) {
      score++;
    }

    new_value = score;
    if(!is_ordinal_variable[variable])
      new_value -= baseline_category[variable];
    old_value = observations(person, variable);

    if(old_value != new_value) {
      // Update raw observations
      observations(person, variable) = new_value;

      // Update sufficient statistics for main effects
      if(is_ordinal_variable[variable] == true) {
        arma::imat counts_per_category_group = counts_per_category[group];
        if(old_value > 0)
          counts_per_category_group(old_value-1, variable)--;
        if(new_value > 0)
          counts_per_category_group(new_value-1, variable)++;
        counts_per_category[group] = counts_per_category_group;
      } else {
        arma::imat blume_capel_stats_group = blume_capel_stats[group];
        blume_capel_stats_group(0, variable) -= old_value;
        blume_capel_stats_group(0, variable) += new_value;
        blume_capel_stats_group(1, variable) -= old_value * old_value;
        blume_capel_stats_group(1, variable) += new_value * new_value;
        blume_capel_stats[group] = blume_capel_stats_group;
      }

      // Update sufficient statistics for pairwise effects
      const int r0 = group_indices(group, 0);
      const int r1 = group_indices(group, 1);
      arma::mat obs = arma::conv_to<arma::mat>::from(observations.rows(r0, r1));
      arma::mat pairwise_stats_group = obs.t() * obs; // crossprod
      pairwise_stats[group] = pairwise_stats_group;
    }
  }
  return;
}



// Updates main effect parameters in bgmCompare using a random-walk Metropolis step.
//
// For each variable, the function proposes new parameter values for either:
//   - all categories (ordinal variables), or
//   - two parameters (linear and quadratic, Blume–Capel variables).
//
// If group-specific differences are enabled (`inclusion_indicator(v,v) == 1`),
// additional parameters (one per group contrast) are also updated.
//
// Each proposed parameter is evaluated via
// `log_pseudoposterior_main_component()`, and accepted/rejected using
// the Metropolis–Hastings rule. Proposal standard deviations are adapted
// online with `MetropolisAdaptationController`.
//
// Workflow:
//  1. Iterate over all variables.
//  2. For each category (ordinal) or parameter (Blume–Capel):
//      - Update the "overall" effect (h=0).
//      - Optionally update group-difference effects (h=1..G-1).
//  3. Record acceptance probabilities and update adaptation statistics.
//
// Inputs:
//  - main_effects: Matrix of main effect parameters [rows = effects, cols = groups];
//                  updated in place.
//  - pairwise_effects: Current pairwise effects (passed through to log posterior).
//  - main_effect_indices: Row index ranges for each variable’s main effects.
//  - pairwise_effect_indices: Index map for pairwise effects.
//  - inclusion_indicator: Indicator matrix; diagonal entries control group differences.
//  - projection: Group projection matrix.
//  - num_categories: Number of categories for each variable.
//  - observations: Data matrix [persons × variables].
//  - num_groups: Number of groups (G).
//  - group_indices: Row ranges per group in `observations`.
//  - counts_per_category, blume_capel_stats: Group-specific sufficient statistics.
//  - is_ordinal_variable: Indicator for ordinal vs. Blume–Capel.
//  - baseline_category: Reference categories (Blume–Capel only).
//  - difference_prior: Prior (BaseParameterPrior) on group differences.
//  - threshold_prior: Prior (BaseParameterPrior) on main effects (thresholds).
//  - iteration: Current iteration index (for adaptation).
//  - rwm_adapt: Adaptation controller for proposal SDs.
//  - rng: Random number generator.
//  - proposal_sd_main: Proposal standard deviations [same shape as `main_effects`];
//                      updated in place.
//
// Notes:
//  - Acceptance probabilities are stored per parameter and fed to `metropolis_adapt.update()`.
//  - This function does not alter pairwise effects, but passes them into
//    the posterior for likelihood consistency.
//  - The helper lambda `do_update` encapsulates the proposal/accept/revert loop
//    for a single parameter, improving readability.
void update_main_effects_metropolis_bgmcompare (
    arma::mat& main_effects,
    arma::mat& pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::imat& inclusion_indicator,
    const arma::mat& projection,
    const arma::ivec& num_categories,
    const arma::imat& observations,
    const int num_groups,
    const arma::imat& group_indices,
    const std::vector<arma::imat>& counts_per_category,
    const std::vector<arma::imat>& blume_capel_stats,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const int iteration,
    MetropolisAdaptationController& metropolis_adapt,
    SafeRNG& rng,
    arma::mat& proposal_sd_main,
    const BaseParameterPrior& difference_prior,
    const BaseParameterPrior& threshold_prior
) {
  const int num_vars = observations.n_cols;
  arma::umat index_mask_main = arma::zeros<arma::umat>(proposal_sd_main.n_rows,
                                                       proposal_sd_main.n_cols);
  arma::mat accept_prob_main = arma::zeros<arma::mat>(proposal_sd_main.n_rows,
                                                      proposal_sd_main.n_cols);

  // --- helper for one update ---
  auto do_update = [&](int row, int variable, int category, int par, int h) {
    index_mask_main(row, h) = 1;
    double& current = main_effects(row, h);
    double proposal_sd = proposal_sd_main(row, h);

    auto log_post = [&](double theta) {
      main_effects(row, h) = theta;
      return log_pseudoposterior_main_component(
        main_effects, pairwise_effects, main_effect_indices,
        pairwise_effect_indices, projection, observations,
        group_indices, num_categories, counts_per_category,
        blume_capel_stats, num_groups, inclusion_indicator,
        is_ordinal_variable, baseline_category,
        variable, category, par, h,
        difference_prior, threshold_prior
      );
    };

    StepResult result = metropolis_step(current, proposal_sd, log_post, rng);
    current = result.state[0];
    accept_prob_main(row, h) = result.accept_prob;
  };

  // --- loop over variables ---
  for (int variable = 0; variable < num_vars; ++variable) {
    int base_category_index = main_effect_indices(variable, 0);
    int num_cats = num_categories(variable);
    bool group_differences = (inclusion_indicator(variable, variable) == 1);

    if (is_ordinal_variable[variable]) {
      // ordinal: loop categories
      for (int category = 0; category < num_cats; ++category) {
        int row = base_category_index + category;
        int hmax = group_differences ? num_groups : 1;
        for (int h = 0; h < hmax; ++h)
          do_update(row, variable, category, -1, h);
      }
    } else {
      // non-ordinal: two parameters
      for (int par = 0; par < 2; ++par) {
        int row = base_category_index + par;
        int hmax = group_differences ? num_groups : 1;
        for (int h = 0; h < hmax; ++h)
          do_update(row, variable, -1, par, h);
      }
    }
  }

  metropolis_adapt.update(index_mask_main, accept_prob_main, iteration);
}




// Updates pairwise interaction parameters in bgmCompare using a random-walk
// Metropolis step.
//
// For each variable pair (var1,var2), the function proposes new parameter
// values for:
//   - the overall interaction (h=0), and
//   - optionally group-difference effects (h=1..G-1) if enabled in
//     `inclusion_indicator`.
//
// Each proposed parameter is evaluated via
// `log_pseudoposterior_pair_component()`, and accepted/rejected using the
// Metropolis–Hastings rule. Proposal standard deviations are adapted online
// through `MetropolisAdaptationController`.
//
// Workflow:
//  1. Iterate over all unique pairs of variables.
//  2. For each pair, update the overall effect (h=0).
//  3. If group differences are active, update group-specific difference
//     effects (h=1..G-1).
//  4. Record acceptance probabilities and update proposal SDs via `rwm_adapt`.
//
// Inputs:
//  - main_effects: Matrix of main effect parameters (passed through to log posterior).
//  - pairwise_effects: Matrix of pairwise interaction parameters [rows = pairs, cols = groups];
//                      updated in place.
//  - main_effect_indices: Index map for main effects (per variable).
//  - pairwise_effect_indices: Row index map for pairwise effects (per var1,var2).
//  - inclusion_indicator: Indicator matrix; off-diagonal entries control group differences.
//  - projection: Group projection matrix.
//  - num_categories: Number of categories per variable.
//  - observations: Data matrix [persons × variables].
//  - num_groups: Number of groups (G).
//  - group_indices: Row ranges per group in `observations`.
//  - pairwise_stats: Group-specific sufficient statistics for pairwise effects.
//  - is_ordinal_variable: Indicator for ordinal vs. Blume–Capel variables.
//  - baseline_category: Reference categories (Blume–Capel only).
//  - pairwise_scale: Scale parameter for overall interaction prior.
//  - difference_scale: Scale parameter for group difference priors.
//  - iteration: Current iteration index (for adaptation).
//  - rwm_adapt: Adaptation controller for proposal SDs.
//  - rng: Random number generator.
//  - proposal_sd_pair: Proposal standard deviations [same shape as `pairwise_effects`];
//                      updated in place.
//
// Notes:
//  - Acceptance probabilities are tracked per parameter and fed to
//    `metropolis_adapt.update()`.
//  - The helper lambda `do_update` encapsulates the proposal/accept/reject
//    logic for a single parameter.
void update_pairwise_effects_metropolis_bgmcompare (
    arma::mat& main_effects,
    arma::mat& pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::imat& inclusion_indicator,
    const arma::mat& projection,
    const arma::ivec& num_categories,
    const arma::imat& observations,
    const int num_groups,
    const arma::imat& group_indices,
    const std::vector<arma::mat>& pairwise_stats,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const int iteration,
    MetropolisAdaptationController& metropolis_adapt,
    SafeRNG& rng,
    arma::mat& proposal_sd_pair,
    const BaseParameterPrior& interaction_prior,
    const BaseParameterPrior& difference_prior
) {
  int num_variables = observations.n_cols;
  int num_pairs = num_variables * (num_variables - 1) / 2;
  arma::mat accept_prob_pair = arma::zeros<arma::mat>(num_pairs, num_groups);
  arma::umat index_mask_pair = arma::zeros<arma::umat>(num_pairs, num_groups);

  // --- Build group-specific pairwise matrices and residual matrices ---
  std::vector<arma::mat> pairwise_groups(num_groups);
  std::vector<arma::mat> residual_matrices(num_groups);
  std::vector<arma::mat> obs_double_groups(num_groups);

  for (int g = 0; g < num_groups; ++g) {
    const int r0 = group_indices(g, 0);
    const int r1 = group_indices(g, 1);
    const arma::vec proj_g = projection.row(g).t();

    // Build group-specific pairwise effects matrix
    arma::mat pairwise_group(num_variables, num_variables, arma::fill::zeros);
    for (int v = 0; v < num_variables - 1; ++v) {
      for (int u = v + 1; u < num_variables; ++u) {
        double w = compute_group_pairwise_effects(
          v, u, num_groups, pairwise_effects, pairwise_effect_indices,
          inclusion_indicator, proj_g
        );
        pairwise_group(v, u) = w;
        pairwise_group(u, v) = w;
      }
    }
    pairwise_groups[g] = pairwise_group;

    // Convert observations to double and compute residual matrix
    obs_double_groups[g] = arma::conv_to<arma::mat>::from(observations.rows(r0, r1));
    residual_matrices[g] = obs_double_groups[g] * pairwise_group;
  }

  // --- helper for one update using optimized residual-based function ---
  auto do_update = [&](int var1, int var2, int h) {
    int idx = pairwise_effect_indices(var1, var2);
    index_mask_pair(idx, h) = 1;
    double current = pairwise_effects(idx, h);
    double proposal_sd = proposal_sd_pair(idx, h);

    auto log_post = [&](double theta) {
      double delta = theta - current;
      return log_pseudoposterior_pair_component(
        main_effects, pairwise_effects, main_effect_indices,
        pairwise_effect_indices, projection, observations, group_indices,
        num_categories, pairwise_stats, residual_matrices, num_groups,
        inclusion_indicator, is_ordinal_variable, baseline_category,
        var1, var2, h, delta,
        interaction_prior, difference_prior
      );
    };

    StepResult result = metropolis_step(current, proposal_sd, log_post, rng);
    double value = result.state[0];

    // Update residual matrices if move was accepted
    if (current != value) {
      double delta = value - current;
      pairwise_effects(idx, h) = value;

      for (int g = 0; g < num_groups; ++g) {
        const arma::vec proj_g = projection.row(g).t();
        double delta_g = (h == 0) ? delta : delta * proj_g(h - 1);

        // Update pairwise_group for this group
        pairwise_groups[g](var1, var2) += delta_g;
        pairwise_groups[g](var2, var1) += delta_g;

        // Update residual matrix columns
        residual_matrices[g].col(var1) += obs_double_groups[g].col(var2) * delta_g;
        residual_matrices[g].col(var2) += obs_double_groups[g].col(var1) * delta_g;
      }
    }

    accept_prob_pair(idx, h) = result.accept_prob;
  };

  for (int var1 = 0; var1 < num_variables - 1; var1++) {
    for (int var2 = var1 + 1; var2 < num_variables; var2++) {
      bool group_differences = (inclusion_indicator(var1, var2) == 1);
      int hmax = group_differences ? num_groups : 1;
      for (int h = 0; h < hmax; h++) {
        do_update(var1, var2, h);
      }
    }
  }

  metropolis_adapt.update(index_mask_pair, accept_prob_pair, iteration);
}



// Heuristically determine an initial NUTS step size for bgmCompare.
//
// This function vectorizes the current model parameters, then repeatedly
// simulates short HMC trajectories to calibrate a stable starting step size
// that achieves a target acceptance rate.
//
// Workflow:
//  1. Vectorize current parameters into a single state vector.
//  2. Define closures for log-posterior evaluation and gradient computation:
//     - `log_post`: unpacks parameters and evaluates the log pseudoposterior.
//     - `grad`: unpacks parameters and evaluates the gradient of the
//       pseudoposterior.
//  3. Pass these to `heuristic_initial_step_size`, which runs the heuristic
//     tuning loop.
//
// Inputs:
//  - main_effects: Matrix of main-effect parameters [n_main_rows × G].
//  - pairwise_effects: Matrix of pairwise interaction parameters [n_pairs × G].
//  - main_effect_indices: Row index ranges for each variable’s main effects.
//  - pairwise_effect_indices: Row index map for pairwise effects.
//  - inclusion_indicator: Matrix marking which main and pairwise differences
//                         are active.
//  - projection: Group projection matrix (encodes contrasts).
//  - num_categories: Number of categories per variable [V].
//  - observations: Data matrix [persons × variables].
//  - num_groups: Number of groups (G).
//  - group_indices: Row ranges per group in `observations`.
//  - counts_per_category: Per-group sufficient statistics for ordinal variables.
//  - blume_capel_stats: Per-group sufficient statistics for Blume–Capel variables.
//  - pairwise_stats: Per-group sufficient statistics for pairwise effects.
//  - is_ordinal_variable: Indicator for ordinal vs. Blume–Capel variables [V].
//  - baseline_category: Reference categories for Blume–Capel variables [V].
//  - pairwise_scale: Scale parameter for overall pairwise priors.
//  - difference_scale: Scale parameter for group difference priors.
//  - main_alpha, main_beta: Hyperparameters for Beta prior on main effects.
//  - target_acceptance: Desired acceptance probability (e.g. 0.8).
//  - rng: Random number generator.
//
// Returns:
//  - A double value for the initial NUTS step size.
//
// Notes:
//  - This routine is only used during warmup to initialize
//    `NUTSAdaptationController`.
//  - Correct indexing of parameters relies on `build_index_maps` to ensure
//    consistency between vectorization and gradient computation.
double find_initial_stepsize_bgmcompare(
    arma::mat& main_effects,
    arma::mat& pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::imat& inclusion_indicator,
    const arma::mat& projection,
    const arma::ivec& num_categories,
    const arma::imat& observations,
    const int num_groups,
    const arma::imat& group_indices,
    const std::vector<arma::imat>& counts_per_category,
    const std::vector<arma::imat>& blume_capel_stats,
    const std::vector<arma::mat>& pairwise_stats,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const double target_acceptance,
    SafeRNG& rng,
    const BaseParameterPrior& interaction_prior,
    const BaseParameterPrior& difference_prior,
    const BaseParameterPrior& threshold_prior
) {
  arma::vec theta = vectorize_model_parameters_bgmcompare(
    main_effects, pairwise_effects, inclusion_indicator, main_effect_indices,
    pairwise_effect_indices, num_categories, is_ordinal_variable
  );
  arma::mat current_main = main_effects;
  arma::mat current_pair = pairwise_effects;

  // Pre-convert observations to double once (avoids repeated conversion in gradient evaluations)
  const arma::mat obs_double = arma::conv_to<arma::mat>::from(observations);

  auto index_maps = build_index_maps(
    main_effects, pairwise_effects,
    inclusion_indicator,
    main_effect_indices,
    pairwise_effect_indices, num_categories, is_ordinal_variable
  );
  auto& main_index = index_maps.first;
  auto& pair_index = index_maps.second;

  arma::vec grad_obs_act = gradient_observed_active(
    main_effect_indices, pairwise_effect_indices, projection,
    observations, group_indices, num_categories, inclusion_indicator,
    counts_per_category, blume_capel_stats, pairwise_stats,
    num_groups, is_ordinal_variable, baseline_category,
    main_index, pair_index
  );

  auto grad = [&](const arma::vec& theta_vec) {
    unvectorize_model_parameters_bgmcompare(
      theta_vec, current_main, current_pair, inclusion_indicator,
      main_effect_indices, pairwise_effect_indices, num_groups, num_categories,
      is_ordinal_variable
    );

    return gradient(
      current_main, current_pair, main_effect_indices, pairwise_effect_indices,
      projection, obs_double, group_indices, num_categories,
      counts_per_category, blume_capel_stats,
      pairwise_stats, num_groups, inclusion_indicator,
      is_ordinal_variable, baseline_category,
      main_index, pair_index,
      grad_obs_act,
      interaction_prior, difference_prior, threshold_prior
    );
  };

  auto joint = [&](const arma::vec& theta_vec) {
    unvectorize_model_parameters_bgmcompare(
      theta_vec, current_main, current_pair, inclusion_indicator,
      main_effect_indices, pairwise_effect_indices, num_groups, num_categories,
      is_ordinal_variable
    );

    return logp_and_gradient(
      current_main, current_pair, main_effect_indices, pairwise_effect_indices,
      projection, obs_double, group_indices, num_categories,
      counts_per_category, blume_capel_stats,
      pairwise_stats, num_groups, inclusion_indicator,
      is_ordinal_variable, baseline_category,
      main_index, pair_index, grad_obs_act,
      interaction_prior, difference_prior, threshold_prior
    );
  };

  return heuristic_initial_step_size(theta, grad, joint, rng, target_acceptance);
}



// Perform one No-U-Turn Sampler (NUTS) update step for the bgmCompare model.
//
// The function:
//  1. Vectorizes the current parameter state (main + pairwise effects).
//  2. Defines closures for log-posterior evaluation and gradient calculation.
//  3. Runs a NUTS trajectory (adaptive tree-based extension of HMC).
//  4. Unpacks the accepted state back into main and pairwise matrices.
//  5. Updates the adaptation controller with acceptance probability.
//
// Inputs:
//  - main_effects, pairwise_effects: Current parameter matrices, updated in place.
//  - main_effect_indices, pairwise_effect_indices: Index maps for parameters.
//  - inclusion_indicator: Indicates active main and pairwise differences.
//  - projection: Group projection matrix for contrasts.
//  - num_categories: Number of categories per variable [V].
//  - observations: Data matrix [N × V].
//  - num_groups: Number of groups.
//  - group_indices: Row ranges for each group in `observations`.
//  - counts_per_category, blume_capel_stats: Per-group sufficient statistics.
//  - pairwise_stats: Per-group pairwise sufficient statistics.
//  - is_ordinal_variable: Marks ordinal vs. Blume–Capel variables.
//  - baseline_category: Reference categories for Blume–Capel variables.
//  - pairwise_scale: Scale of overall pairwise prior.
//  - difference_scale: Scale of group-difference prior.
//  - main_alpha, main_beta: Hyperparameters for main-effect priors.
//  - nuts_max_depth: Maximum tree depth for NUTS doubling procedure.
//  - iteration: Current sampler iteration (for adaptation scheduling).
//  - nuts_adapt: Adaptation controller for step size and mass matrix.
//  - learn_mass_matrix: Whether to adapt the mass matrix (unused inside NUTS but relevant to controller).
//  - selection: If true, restrict mass matrix to active parameters only.
//  - rng: Random number generator.
//
// Returns:
//  - A `StepResult` containing the accepted state and diagnostics
//    (e.g. tree depth, divergences, energy).
//
// Notes:
//  - This variant is specific to bgmCompare, where parameters are stored in
//    row-wise structures with group-difference columns.
//  - Consistency between vectorization, unvectorization, and gradient
//    indexing is enforced via `build_index_maps`.
//  - Diagnostics from the returned `StepResult` can be used to monitor
//    sampler stability (e.g. divergences, tree depth).
StepResult update_nuts_bgmcompare(
    arma::mat& main_effects,
    arma::mat& pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::imat& inclusion_indicator,
    const arma::mat& projection,
    const arma::ivec& num_categories,
    const arma::imat& observations,
    const int num_groups,
    const arma::imat& group_indices,
    const std::vector<arma::imat>& counts_per_category,
    const std::vector<arma::imat>& blume_capel_stats,
    const std::vector<arma::mat>& pairwise_stats,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const int nuts_max_depth,
    const int iteration,
    NUTSAdaptationController& nuts_adapt,
    const bool learn_mass_matrix,
    const bool selection,
    SafeRNG& rng,
    const BaseParameterPrior& interaction_prior,
    const BaseParameterPrior& difference_prior,
    const BaseParameterPrior& threshold_prior
) {
  // Pre-convert observations to double once (avoids repeated conversion in gradient evaluations)
  const arma::mat obs_double = arma::conv_to<arma::mat>::from(observations);

  arma::vec current_state = vectorize_model_parameters_bgmcompare(
    main_effects, pairwise_effects, inclusion_indicator,
    main_effect_indices, pairwise_effect_indices, num_categories,
    is_ordinal_variable
  );

  arma::mat current_main = main_effects;
  arma::mat current_pair = pairwise_effects;

  auto index_maps = build_index_maps(
    main_effects, pairwise_effects,
    inclusion_indicator,
    main_effect_indices,
    pairwise_effect_indices, num_categories, is_ordinal_variable
  );
  auto& main_index = index_maps.first;
  auto& pair_index = index_maps.second;

  arma::vec grad_obs_act = gradient_observed_active(
    main_effect_indices, pairwise_effect_indices, projection,
    observations, group_indices, num_categories, inclusion_indicator,
    counts_per_category, blume_capel_stats, pairwise_stats,
    num_groups, is_ordinal_variable, baseline_category,
    main_index, pair_index
  );

  auto grad = [&](const arma::vec& theta_vec) {
    unvectorize_model_parameters_bgmcompare(
      theta_vec, current_main, current_pair, inclusion_indicator,
      main_effect_indices, pairwise_effect_indices, num_groups, num_categories,
      is_ordinal_variable
    );

    return gradient(
      current_main, current_pair, main_effect_indices, pairwise_effect_indices,
      projection, obs_double, group_indices, num_categories,
      counts_per_category, blume_capel_stats,
      pairwise_stats, num_groups, inclusion_indicator,
      is_ordinal_variable, baseline_category,
      main_index, pair_index,
      grad_obs_act,
      interaction_prior, difference_prior, threshold_prior
    );
  };

  auto joint = [&](const arma::vec& theta_vec) {
    unvectorize_model_parameters_bgmcompare(
      theta_vec, current_main, current_pair, inclusion_indicator,
      main_effect_indices, pairwise_effect_indices, num_groups, num_categories,
      is_ordinal_variable
    );

    return logp_and_gradient(
      current_main, current_pair, main_effect_indices, pairwise_effect_indices,
      projection, obs_double, group_indices, num_categories,
      counts_per_category, blume_capel_stats,
      pairwise_stats, num_groups, inclusion_indicator,
      is_ordinal_variable, baseline_category,
      main_index, pair_index, grad_obs_act,
      interaction_prior, difference_prior, threshold_prior
    );
  };

  //adapt
  arma::vec active_inv_mass = inv_mass_active(
    nuts_adapt.inv_mass_diag(), inclusion_indicator, num_groups, num_categories,
    is_ordinal_variable, main_index, pair_index, main_effect_indices,
    pairwise_effect_indices, selection
  );

  StepResult result = nuts_step(
    current_state, nuts_adapt.current_step_size(), joint,
    active_inv_mass, rng, nuts_max_depth
  );

  current_state = result.state;
  unvectorize_model_parameters_bgmcompare(
    current_state, main_effects, pairwise_effects, inclusion_indicator,
    main_effect_indices, pairwise_effect_indices, num_groups, num_categories,
    is_ordinal_variable
  );

  nuts_adapt.update(current_state, result.accept_prob, iteration);

  // If mass matrix was just updated, re-run the heuristic to find a good
  // step size for the new mass matrix. Use current step size as starting point.
  if (nuts_adapt.mass_matrix_just_updated()) {
    arma::vec new_inv_mass = inv_mass_active(
      nuts_adapt.inv_mass_diag(), inclusion_indicator, num_groups, num_categories,
      is_ordinal_variable, main_index, pair_index, main_effect_indices,
      pairwise_effect_indices, selection
    );
    double current_eps = nuts_adapt.current_step_size();
    double new_eps = heuristic_initial_step_size(
      current_state, grad, joint, new_inv_mass, rng,
      nuts_adapt.target_acceptance(),
      current_eps   // init_step: use current step size as starting point
    );
    nuts_adapt.reinit_stepsize(new_eps);
  }

  return result;
}



// Adapt proposal standard deviations (SDs) for main and pairwise effects
// during the warmup phase of the bgmCompare sampler.
//
// This function uses a Robbins–Monro stochastic approximation scheme to
// adjust proposal SDs toward a target acceptance rate. Adaptation occurs
// only when permitted by the current warmup schedule.
//
// Workflow:
//  1. For each main effect parameter (ordinal or Blume–Capel), run one
//     random-walk Metropolis (RWM) step, update the parameter, and adjust
//     the proposal SD.
//  2. For each pairwise effect parameter (overall and group differences),
//     do the same.
//  3. Proposal SDs are updated symmetrically across group columns if
//     differences are included.
//
// Inputs:
//  - proposal_sd_main_effects: Current SDs for main effects, updated in place.
//  - proposal_sd_pairwise_effects: Current SDs for pairwise effects, updated in place.
//  - main_effects, pairwise_effects: Parameter matrices, updated in place.
//  - main_effect_indices, pairwise_effect_indices: Index maps for main/pairwise parameters.
//  - inclusion_indicator: Marks which group differences are active.
//  - projection: Group projection matrix.
//  - num_categories: Categories per variable.
//  - observations: Data matrix [N × V].
//  - num_groups: Number of groups.
//  - group_indices: Row ranges per group in `observations`.
//  - counts_per_category, blume_capel_stats: Per-group sufficient statistics for main effects.
//  - pairwise_stats: Per-group sufficient statistics for pairwise effects.
//  - is_ordinal_variable: Marks ordinal vs. Blume–Capel variables.
//  - baseline_category: Reference category for Blume–Capel variables.
//  - pairwise_scale: Scale of Cauchy prior for pairwise effects.
//  - difference_scale: Scale of Cauchy prior for group differences.
//  - main_alpha, main_beta: Hyperparameters for Beta prior on main effects.
//  - iteration: Current iteration (to check schedule stage).
//  - rng: Random number generator.
//  - sched: Warmup schedule controlling when adaptation is active.
//  - target_accept: Desired acceptance probability (default 0.44).
//  - rm_decay: Robbins–Monro decay rate (default 0.75).
//
// Side effects:
//  - Updates `main_effects` and `pairwise_effects` with new parameter values.
//  - Updates `proposal_sd_main_effects` and `proposal_sd_pairwise_effects`.
//
// Notes:
//  - Adapts only when `sched.adapt_proposal_sd(iteration)` is true.
//  - Helps stabilize RWM acceptance rates before switching to sampling.
void tune_proposal_sd_bgmcompare(
    arma::mat& proposal_sd_main_effects,
    arma::mat& proposal_sd_pairwise_effects,
    arma::mat& main_effects,
    arma::mat& pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::imat& inclusion_indicator,
    const arma::mat& projection,
    const arma::ivec& num_categories,
    const arma::imat& observations,
    int num_groups,
    const arma::imat& group_indices,
    const std::vector<arma::imat>& counts_per_category,
    const std::vector<arma::imat>& blume_capel_stats,
    const std::vector<arma::mat>& pairwise_stats,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    int iteration,
    SafeRNG& rng,
    const WarmupSchedule& sched,
    const BaseParameterPrior& interaction_prior,
    const BaseParameterPrior& difference_prior,
    const BaseParameterPrior& threshold_prior,
    double target_accept = 0.44,
    double rm_decay = 0.75)
{
  if (!sched.adapt_proposal_sd(iteration)) return;

  // Robbins–Monro weight
  double t = iteration - sched.stage3b_start + 1;
  double rm_weight = std::pow(t, -rm_decay);

  const int V = observations.n_cols;

  // --- MAIN EFFECTS ---
  for (int var = 0; var < V; ++var) {
    int start = main_effect_indices(var, 0);
    bool group_differences = (inclusion_indicator(var, var) == 1);
    int hmax = group_differences ? num_groups : 1;

    if (is_ordinal_variable[var]) {
      // Ordinal variable: num_categories[var] free threshold parameters,
      // matching the Metropolis sweep and count_num_main_effects.
      int ncat = num_categories[var];
      for (int c = 0; c < ncat; ++c) {
        int row = start + c;
        for (int h = 0; h < hmax; ++h) {
          double& current = main_effects(row, h);
          double& prop_sd = proposal_sd_main_effects(row, h);

          auto log_post = [&](double theta) {
            main_effects(row, h) = theta;
            return log_pseudoposterior_main_component(
              main_effects, pairwise_effects,
              main_effect_indices, pairwise_effect_indices,
              projection, observations, group_indices,
              num_categories, counts_per_category, blume_capel_stats,
              num_groups, inclusion_indicator, is_ordinal_variable,
              baseline_category,
              var, c, -1, h,
              difference_prior, threshold_prior
            );
          };

          StepResult result = metropolis_step(current, prop_sd, log_post, rng);
          current = result.state[0];
          prop_sd = update_proposal_sd_with_robbins_monro(
            prop_sd, MY_LOG(result.accept_prob), rm_weight, target_accept
          );
        }
      }
    } else {
      // Non-ordinal variable: two parameters
      for (int par = 0; par < 2; ++par) {
        int row = start + par;
        for (int h = 0; h < hmax; ++h) {
          double& current = main_effects(row, h);
          double& prop_sd = proposal_sd_main_effects(row, h);

          auto log_post = [&](double theta) {
            main_effects(row, h) = theta;
            return log_pseudoposterior_main_component(
              main_effects, pairwise_effects,
              main_effect_indices, pairwise_effect_indices,
              projection, observations, group_indices,
              num_categories, counts_per_category, blume_capel_stats,
              num_groups, inclusion_indicator, is_ordinal_variable,
              baseline_category,
              var, -1, par, h,
              difference_prior, threshold_prior
            );
          };

          StepResult result = metropolis_step(current, prop_sd, log_post, rng);
          current = result.state[0];
          prop_sd = update_proposal_sd_with_robbins_monro(
            prop_sd, MY_LOG(result.accept_prob), rm_weight, target_accept
          );
        }
      }
    }
  }

  // --- PAIRWISE EFFECTS ---
  // Build group-specific pairwise matrices and residual matrices
  std::vector<arma::mat> pairwise_groups(num_groups);
  std::vector<arma::mat> residual_matrices(num_groups);
  std::vector<arma::mat> obs_double_groups(num_groups);

  for (int g = 0; g < num_groups; ++g) {
    const int r0 = group_indices(g, 0);
    const int r1 = group_indices(g, 1);
    const arma::vec proj_g = projection.row(g).t();

    arma::mat pairwise_group(V, V, arma::fill::zeros);
    for (int v = 0; v < V - 1; ++v) {
      for (int u = v + 1; u < V; ++u) {
        double w = compute_group_pairwise_effects(
          v, u, num_groups, pairwise_effects, pairwise_effect_indices,
          inclusion_indicator, proj_g
        );
        pairwise_group(v, u) = w;
        pairwise_group(u, v) = w;
      }
    }
    pairwise_groups[g] = pairwise_group;
    obs_double_groups[g] = arma::conv_to<arma::mat>::from(observations.rows(r0, r1));
    residual_matrices[g] = obs_double_groups[g] * pairwise_group;
  }

  for (int v1 = 0; v1 < V - 1; ++v1) {
    for (int v2 = v1 + 1; v2 < V; ++v2) {
      int idx = pairwise_effect_indices(v1, v2);
      bool group_differences = (inclusion_indicator(v1, v2) == 1);
      int hmax = group_differences ? num_groups : 1;

      for (int h = 0; h < hmax; ++h) {
        double& current = pairwise_effects(idx, h);
        double& prop_sd = proposal_sd_pairwise_effects(idx, h);

        auto log_post = [&](double theta) {
          double delta = theta - current;
          return log_pseudoposterior_pair_component(
            main_effects, pairwise_effects,
            main_effect_indices, pairwise_effect_indices,
            projection, observations, group_indices,
            num_categories, pairwise_stats, residual_matrices, num_groups,
            inclusion_indicator, is_ordinal_variable, baseline_category,
            v1, v2, h, delta,
            interaction_prior, difference_prior
          );
        };

        StepResult result = metropolis_step(current, prop_sd, log_post, rng);
        double value = result.state[0];

        if (current != value) {
          double delta = value - current;
          current = value;

          for (int g = 0; g < num_groups; ++g) {
            const arma::vec proj_g = projection.row(g).t();
            double delta_g = (h == 0) ? delta : delta * proj_g(h - 1);

            pairwise_groups[g](v1, v2) += delta_g;
            pairwise_groups[g](v2, v1) += delta_g;

            residual_matrices[g].col(v1) += obs_double_groups[g].col(v2) * delta_g;
            residual_matrices[g].col(v2) += obs_double_groups[g].col(v1) * delta_g;
          }
        }

        prop_sd = update_proposal_sd_with_robbins_monro(
          prop_sd, MY_LOG(result.accept_prob), rm_weight, target_accept
        );
      }
    }
  }
}



// Metropolis–Hastings updates for difference-inclusion indicators in bgmCompare.
//
// This function toggles whether group-level differences are included for
// main effects (diagonal entries of `inclusion_indicator`) and pairwise
// effects (off-diagonal entries). Each update proposes either:
//  - Turning a currently excluded difference “on” by drawing a new non-zero
//    value from a Gaussian proposal, or
//  - Turning an included difference “off” by setting its value(s) to zero.
//
// The acceptance probability combines:
//  - Pseudolikelihood ratio (data contribution),
//  - Prior ratio on inclusion indicators,
//  - Prior ratio on parameter values (Cauchy vs. point-mass-at-zero),
//  - Proposal density correction.
//
// Inputs:
//  - inclusion_probability_difference: Prior inclusion probabilities for
//    group differences [V × V].
//  - index: Matrix mapping pairwise interactions to variable indices.
//  - main_effects, pairwise_effects: Parameter matrices, updated in place.
//  - main_effect_indices, pairwise_effect_indices: Index maps for parameters.
//  - projection: Group projection matrix.
//  - observations: Data matrix [N × V].
//  - num_groups: Number of groups.
//  - group_indices: Row ranges per group in `observations`.
//  - num_categories: Categories per variable [V × G].
//  - inclusion_indicator: Indicator matrix for differences, updated in place.
//  - is_ordinal_variable: Marks ordinal vs. Blume–Capel variables [V].
//  - baseline_category: Reference category for Blume–Capel variables [V].
//  - proposal_sd_main, proposal_sd_pairwise: Proposal SD matrices for main
//    and pairwise effects.
//  - difference_scale: Scale of Cauchy prior for group differences.
//  - counts_per_category, blume_capel_stats: Per-group sufficient statistics
//    for main effects.
//  - pairwise_stats: Per-group sufficient statistics for pairwise effects.
//  - rng: Random number generator.
//
// Side effects:
//  - Updates `inclusion_indicator` entries for main/pairwise differences.
//  - Updates corresponding slices of `main_effects` and `pairwise_effects`.
//
// Notes:
//  - For main effects, differences correspond to columns 1..G-1 of the
//    parameter matrix.
//  - For pairwise effects, differences correspond to columns 1..G-1 of the
//    pairwise-effect matrix rows.
//  - Ensures symmetry of `inclusion_indicator` for pairwise updates.
void update_indicator_differences_metropolis_bgmcompare (
    const arma::mat& inclusion_probability_difference,
    const arma::imat& index,
    arma::mat& main_effects,
    arma::mat& pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::mat& projection,
    const arma::imat& observations,
    const int num_groups,
    const arma::imat& group_indices,
    const arma::imat& num_categories,
    arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const arma::mat& proposal_sd_main,
    const arma::mat& proposal_sd_pairwise,
    const std::vector<arma::imat>& counts_per_category,
    const std::vector<arma::imat>& blume_capel_stats,
    const std::vector<arma::mat>& pairwise_stats,
    const bool main_difference_selection,
    SafeRNG& rng,
    const BaseParameterPrior& difference_prior
) {
  const int num_variables = observations.n_cols;

  // --- main effects ---
  // Skip main effect indicator updates if main_difference_selection is disabled
  if (main_difference_selection) {
    for(int var = 0; var < num_variables; var++) {
      int start = main_effect_indices(var, 0);
    int stop = main_effect_indices(var, 1);

    arma::mat current_main_effects = main_effects;
    arma::mat proposed_main_effects = main_effects;
    int current_ind = inclusion_indicator(var, var);
    int proposed_ind = 1 - current_ind;

    for(int row = start; row <= stop; row++) {
      for(int h = 1; h < num_groups; h++) {
        if(proposed_ind == 0) {
          // Propose to set difference to zero value
          proposed_main_effects(row, h) = 0.0;
        } else {
          // Propose to set difference to non-zero value
          proposed_main_effects(row, h) = rnorm(
            rng, 0.0, proposal_sd_main(row, h)
          );
        }
      }
    }

    // Calculate log acceptance probability
    double log_accept = log_pseudolikelihood_ratio_main(
      current_main_effects, proposed_main_effects, pairwise_effects,
      main_effect_indices, pairwise_effect_indices, projection,
      observations, group_indices, num_categories, counts_per_category,
      blume_capel_stats, num_groups, inclusion_indicator,
      is_ordinal_variable, baseline_category, var
    );

    // Add prior inclusion probability contribution
    double inc_prob = inclusion_probability_difference(var, var);
    double logit_inc_prob = MY_LOG(inc_prob / (1 - inc_prob));
    if(proposed_ind == 1) {
      log_accept += logit_inc_prob;
    } else {
      log_accept -= logit_inc_prob;
    }

    // Add parameter prior contribution
    for(int row = start; row <= stop; row++) {
      if(proposed_ind == 1) {
        // Propose to set difference to non-zero
        for(int h = 1; h < num_groups; h++) {
          log_accept += difference_prior.logp(proposed_main_effects(row, h));
          log_accept -= R::dnorm(
            proposed_main_effects(row, h), current_main_effects(row, h),
            proposal_sd_main(row, h), true
          );
        }
      } else {
        // Propose to set difference to zero
        for(int h = 1; h < num_groups; h++) {
          log_accept -= difference_prior.logp(current_main_effects(row, h));
          log_accept += R::dnorm(
            current_main_effects(row, h), proposed_main_effects(row, h),
            proposal_sd_main(row, h), true
          );
        }
      }
    }

    // Perform Metropolis-Hastings step
    double U = runif(rng);
    if(MY_LOG(U) < log_accept) {
      inclusion_indicator(var, var) = proposed_ind;
      main_effects.rows(start, stop).cols(1, num_groups - 1) =
        proposed_main_effects.rows(start, stop).cols(1, num_groups - 1);
    }
  }
  } // end if (main_difference_selection)

  // --- pairwise effects ---
  const int num_pairwise = index.n_rows;
  for (int cntr = 0; cntr < num_pairwise; cntr++) {
    int var1 = index(cntr, 1);
    int var2 = index(cntr, 2);
    int int_index = pairwise_effect_indices(var1, var2);

    arma::mat current_pairwise_effects = pairwise_effects;
    arma::mat proposed_pairwise_effects = pairwise_effects;

    int current_ind = inclusion_indicator(var1, var2);
    int proposed_ind = 1 - current_ind;

    for(int h = 1; h < num_groups; h++) {
      if(proposed_ind == 0) {
        // Propose to set difference to zero value
        proposed_pairwise_effects(int_index, h) = 0.0;
      } else {
        // Propose to set difference to non-zero value
        proposed_pairwise_effects(int_index, h) = rnorm(
          rng, 0.0, proposal_sd_pairwise(int_index, h)
        );
      }
    }
    // Calculate log acceptance probability
    double log_accept = log_pseudolikelihood_ratio_pairwise(
      main_effects, current_pairwise_effects, proposed_pairwise_effects,
      main_effect_indices, pairwise_effect_indices, projection, observations,
      group_indices, num_categories, pairwise_stats, num_groups,
      inclusion_indicator, is_ordinal_variable, baseline_category, var1, var2
    );

    // Add prior inclusion probability contribution
    double inc_prob = inclusion_probability_difference(var1, var2);
    double logit_inc_prob = MY_LOG(inc_prob / (1 - inc_prob));
    if(proposed_ind == 1) {
      log_accept += logit_inc_prob;
    } else {
      log_accept -= logit_inc_prob;
    }

    // Add parameter prior contribution
    if(proposed_ind == 1) {
      // Propose to set difference to non-zero
      for(int h = 1; h < num_groups; h++) {
        log_accept += difference_prior.logp(
          proposed_pairwise_effects(int_index, h)
        );
        log_accept -= R::dnorm(
          proposed_pairwise_effects(int_index, h),
          current_pairwise_effects(int_index, h),
          proposal_sd_pairwise(int_index, h),
          true
        );
      }
    } else {
      // Propose to set difference to zero
      for(int h = 1; h < num_groups; h++) {
        log_accept -= difference_prior.logp(
          current_pairwise_effects(int_index, h)
        );
        log_accept += R::dnorm(
          current_pairwise_effects(int_index, h),
          proposed_pairwise_effects(int_index, h),
          proposal_sd_pairwise(int_index, h), true
        );
      }
    }

    // Metropolis-Hastings acceptance step
    double U = runif(rng);
    if (MY_LOG(U) < log_accept) {
      // Update inclusion inclusion_indicator
      inclusion_indicator(var1, var2) = proposed_ind;
      inclusion_indicator(var2, var1) = proposed_ind;

      // Update pairwise effects and rest matrix
      for (int h = 1; h < num_groups; h++) {
        pairwise_effects(int_index, h) = proposed_pairwise_effects(int_index, h);
      }
    }
  }
}



// Perform one Gibbs update step for the bgmCompare model.
//
// This function executes a single iteration of the Gibbs sampler, including:
//
//  Step 0: (optional) Initialize graph structure if difference selection
//          is enabled and the current iteration marks the start of Stage 3c.
//
//  Step 1: (optional) Update inclusion indicators for group differences
//          (main and pairwise effects) via Metropolis–Hastings proposals.
//
//  Step 2: Update model parameters according to the selected update method:
//    - "adaptive-metropolis": Update main and pairwise effects individually
//      with random-walk Metropolis and adaptive proposal SDs.
//    - "nuts": Update the full parameter vector using the No-U-Turn Sampler.
//      If past burn-in, store NUTS diagnostics (tree depth, divergences, energy).
//
//  Step 3: (Stage 3b only) Adapt proposal SDs for Metropolis updates using
//          Robbins–Monro tuning.
//
// Inputs:
//  - observations: Data matrix [N × V].
//  - num_categories: Number of categories per variable [V].
//  - pairwise_scale, difference_scale: Prior scale parameters.
//  - counts_per_category, blume_capel_stats: Sufficient statistics per group.
//  - main_alpha, main_beta: Hyperparameters for Beta prior on main effects.
//  - inclusion_indicator: Matrix of active group differences [V × V], updated in place.
//  - main_effects, pairwise_effects: Parameter matrices, updated in place.
//  - is_ordinal_variable: Marks ordinal vs. Blume–Capel variables.
//  - baseline_category: Reference categories for Blume–Capel variables.
//  - iteration: Current iteration index.
//  - pairwise_effect_indices, main_effect_indices: Index maps for parameters.
//  - pairwise_stats: Per-group pairwise sufficient statistics.
//  - nuts_max_depth: Maximum tree depth for NUTS.
//  - nuts_adapt: Adaptation controller for NUTS.
//  - metropolis_adapt_main, metropolis_adapt_pair: Adaptation controllers for RWM updates.
//  - learn_mass_matrix: Whether to adapt the mass matrix in NUTS.
//  - schedule: Warmup schedule, controls adaptation and selection phases.
//  - treedepth_samples, divergent_samples, energy_samples: Buffers for NUTS diagnostics.
//  - projection: Group projection matrix.
//  - num_groups: Number of groups.
//  - group_indices: Row ranges per group in `observations`.
//  - rng: Random number generator.
//  - inclusion_probability: Prior probabilities for including differences.
//  - update_method: Update strategy ("adaptive-metropolis", "nuts").
//  - proposal_sd_main, proposal_sd_pair: Proposal SD matrices for Metropolis updates.
//  - index: Index table for pairwise differences.
//
// Side effects:
//  - Updates parameters, inclusion indicators, and sufficient statistics.
//  - Updates adaptation controllers and (if NUTS) diagnostic buffers.
//
// Notes:
//  - This function encapsulates all update logic for bgmCompare.
//  - Choice of `update_method` governs whether updates are local (RWM) or
//    global (NUTS).
void gibbs_update_step_bgmcompare (
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const std::vector<arma::imat>& counts_per_category,
    const std::vector<arma::imat>& blume_capel_stats,
    arma::imat& inclusion_indicator,
    arma::mat& pairwise_effects,
    arma::mat& main_effects,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const int iteration,
    const arma::imat& pairwise_effect_indices,
    const std::vector<arma::mat>& pairwise_stats,
    const int nuts_max_depth,
    NUTSAdaptationController& nuts_adapt,
    MetropolisAdaptationController& metropolis_adapt_main,
    MetropolisAdaptationController& metropolis_adapt_pair,
    const bool learn_mass_matrix,
    WarmupSchedule const& schedule,
    arma::ivec& treedepth_samples,
    arma::ivec& divergent_samples,
    arma::vec& energy_samples,
    arma::vec& accept_prob_samples,
    const arma::imat& main_effect_indices,
    const arma::mat& projection,
    const int num_groups,
    const arma::imat group_indices,
    SafeRNG& rng,
    arma::mat& inclusion_probability,
    const UpdateMethod update_method,
    arma::mat& proposal_sd_main,
    arma::mat& proposal_sd_pair,
    const arma::imat& index,
    const bool main_difference_selection,
    const BaseParameterPrior& interaction_prior,
    const BaseParameterPrior& difference_prior,
    const BaseParameterPrior& threshold_prior
) {

  // Step 0: Initialise random graph structure when edge_selection = TRUE
  if (schedule.selection_enabled(iteration) && iteration == schedule.stage3c_start) {
    initialise_graph_bgmcompare(
      inclusion_indicator, main_effects, pairwise_effects, main_effect_indices,
      pairwise_effect_indices, inclusion_probability, main_difference_selection, rng
    );
  }

  // Step 1: Difference selection via MH indicator updates (if enabled)
  if (schedule.selection_enabled(iteration)) {
    update_indicator_differences_metropolis_bgmcompare (
        inclusion_probability, index, main_effects, pairwise_effects,
        main_effect_indices, pairwise_effect_indices, projection, observations,
        num_groups, group_indices, num_categories, inclusion_indicator,
        is_ordinal_variable, baseline_category, proposal_sd_main,
        proposal_sd_pair, counts_per_category,
        blume_capel_stats, pairwise_stats, main_difference_selection, rng,
        difference_prior
    );
  }

  // Step 2: Update parameters
  if(update_method == adaptive_metropolis) {
    update_main_effects_metropolis_bgmcompare (
        main_effects, pairwise_effects, main_effect_indices,
        pairwise_effect_indices, inclusion_indicator, projection,
        num_categories, observations, num_groups, group_indices,
        counts_per_category, blume_capel_stats, is_ordinal_variable,
        baseline_category, iteration,
        metropolis_adapt_main, rng, proposal_sd_main,
        difference_prior, threshold_prior
    );

    update_pairwise_effects_metropolis_bgmcompare (
        main_effects, pairwise_effects, main_effect_indices,
        pairwise_effect_indices, inclusion_indicator, projection,
        num_categories, observations, num_groups, group_indices,
        pairwise_stats, is_ordinal_variable, baseline_category,
        iteration, metropolis_adapt_pair, rng,
        proposal_sd_pair,
        interaction_prior, difference_prior
    );
  } else if (update_method == nuts) {
    StepResult result = update_nuts_bgmcompare(
      main_effects, pairwise_effects, main_effect_indices,
      pairwise_effect_indices, inclusion_indicator, projection, num_categories,
      observations, num_groups, group_indices, counts_per_category,
      blume_capel_stats, pairwise_stats, is_ordinal_variable,
      baseline_category,
      nuts_max_depth, iteration, nuts_adapt, learn_mass_matrix,
      schedule.selection_enabled(iteration), rng,
      interaction_prior, difference_prior, threshold_prior
    );

    if (iteration >= schedule.total_warmup) {
      int sample_index = iteration - schedule.total_warmup;
      if (auto diag = std::dynamic_pointer_cast<NUTSDiagnostics>(result.diagnostics)) {
        treedepth_samples(sample_index) = diag->tree_depth;
        divergent_samples(sample_index) = diag->divergent ? 1 : 0;
        energy_samples(sample_index) = diag->energy;
        accept_prob_samples(sample_index) = diag->accept_prob;
      }
    }
  }

  /* --- 2b.  proposal-sd tuning during Stage-3b ------------------------------ */
  tune_proposal_sd_bgmcompare(
    proposal_sd_main, proposal_sd_pair, main_effects,
    pairwise_effects, main_effect_indices, pairwise_effect_indices,
    inclusion_indicator, projection, num_categories, observations, num_groups,
    group_indices, counts_per_category, blume_capel_stats,
    pairwise_stats, is_ordinal_variable, baseline_category,
    iteration, rng, schedule,
    interaction_prior, difference_prior, threshold_prior
  );
}



// Run a full Gibbs sampler for the bgmCompare model.
//
// This function controls the full MCMC lifecycle for a single chain:
//  - Initializes parameter matrices, proposal SDs, and adaptation controllers.
//  - Optionally imputes missing data at each iteration.
//  - Executes Gibbs updates for main and pairwise effects, including
//    difference-selection if enabled.
//  - Adapts step size, mass matrix, and proposal SDs during warmup.
//  - Updates inclusion probabilities under the chosen prior
//    (e.g. Beta–Bernoulli).
//  - Collects posterior samples and diagnostics into a `SamplerOutput` struct.
//
// Inputs:
//  - chain_id: Identifier for this chain (1-based).
//  - observations: Data matrix [N × V].
//  - num_groups: Number of groups (G).
//  - counts_per_category: Per-group sufficient statistics (ordinal variables).
//  - blume_capel_stats: Per-group sufficient statistics (Blume–Capel variables).
//  - pairwise_stats: Per-group sufficient statistics for pairwise effects.
//  - num_categories: Number of categories per variable [V].
//  - main_alpha, main_beta: Hyperparameters for Beta prior on main effects.
//  - pairwise_scale: Scale parameter for overall pairwise priors.
//  - difference_scale: Scale parameter for group-difference priors.
//  - difference_selection_alpha, difference_selection_beta: Hyperparameters
//    for difference-selection prior.
//  - difference_prior: Prior type for difference-selection ("Beta-Bernoulli", ...).
//  - iter: Number of post–burn-in sampling iterations.
//  - warmup: Number of warmup iterations.
//  - na_impute: If true, impute missing observations at each iteration.
//  - missing_data_indices: Matrix of [person, variable] indices of missings.
//  - is_ordinal_variable: Marks ordinal vs. Blume–Capel variables.
//  - baseline_category: Reference categories for Blume–Capel variables.
//  - difference_selection: If true, include MH updates for group-difference indicators.
//  - main_effect_indices, pairwise_effect_indices: Index maps for parameter rows.
//  - target_accept: Target acceptance probability (NUTS).
//  - nuts_max_depth: Maximum tree depth for NUTS.
//  - learn_mass_matrix: Whether to adapt the mass matrix (NUTS).
//  - projection: Group projection matrix for contrasts.
//  - group_membership: Mapping of persons to groups.
//  - group_indices: Row ranges per group in `observations`.
//  - interaction_index_matrix: Index map for pairwise interactions.
//  - inclusion_probability: Matrix of prior inclusion probabilities, updated in place.
//  - rng: Random number generator.
//  - update_method: Update strategy ("adaptive-metropolis", "nuts").
//
// Returns:
//  - A `SamplerOutput` struct containing:
//      - main_samples: MCMC samples for main effects.
//      - pairwise_samples: MCMC samples for pairwise effects.
//      - indicator_samples: (optional) Inclusion indicator samples if
//        difference-selection is enabled.
//      - treedepth_samples, divergent_samples, energy_samples:
//        Diagnostics (for NUTS).
//      - chain_id: Identifier for this chain.
//
// Notes:
//  - Warmup is orchestrated via `WarmupSchedule`, which controls adaptation
//    phases and difference-selection activation.
//  - Proposal SDs are tuned via Robbins–Monro during Stage 3b.
//  - Difference-selection updates toggle inclusion indicators and adjust
//    associated parameters with MH proposals.
//  - This function runs entirely in C++ and is wrapped for parallel execution
//    via `GibbsCompareChainRunner`.
bgmCompareOutput run_gibbs_sampler_bgmCompare(
    int chain_id,
    arma::imat observations,
    const int num_groups,
    std::vector<arma::imat>& counts_per_category,
    std::vector<arma::imat>& blume_capel_stats,
    std::vector<arma::mat>& pairwise_stats,
    const arma::ivec& num_categories,
    const double difference_selection_alpha,
    const double difference_selection_beta,
    const std::string& difference_prior_type,
    const int iter,
    const int warmup,
    const bool na_impute,
    const arma::imat& missing_data_indices,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const bool difference_selection,
    const bool main_difference_selection,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const double target_accept,
    const int nuts_max_depth,
    const bool learn_mass_matrix,
    const arma::mat& projection,
    const arma::ivec& group_membership,
    const arma::imat& group_indices,
    const arma::imat& interaction_index_matrix,
    arma::mat inclusion_probability,
    SafeRNG& rng,
    const UpdateMethod update_method,
    ProgressManager& pm,
    const BaseParameterPrior& interaction_prior,
    const BaseParameterPrior& difference_prior,
    const BaseParameterPrior& threshold_prior,
    BaseEdgePrior& difference_edge_prior
) {
  // --- Setup: dimensions and storage structures
  const int num_variables = observations.n_cols;
  const int num_main = count_num_main_effects (
    num_categories, is_ordinal_variable
  );
  const int num_pair = num_variables * (num_variables - 1) / 2;

  // Initialize model parameter matrices
  arma::mat main_effects(num_main, num_groups, arma::fill::zeros);
  arma::mat pairwise_effects(num_pair, num_groups, arma::fill::zeros);
  arma::imat inclusion_indicator(num_variables, num_variables, arma::fill::ones);

  // Allocate storage for MCMC samples. Iterations that never run (user
  // interrupt) keep the fill values: NaN for floating samples, -1 for the
  // integer indicator and allocation samples.
  arma::mat main_effect_samples(iter, num_main * num_groups);
  main_effect_samples.fill(arma::datum::nan);
  arma::mat pairwise_effect_samples(iter, num_pair * num_groups);
  pairwise_effect_samples.fill(arma::datum::nan);
  arma::imat indicator_samples;

  if (difference_selection) {
    indicator_samples.set_size(iter, num_pair + num_variables);
    indicator_samples.fill(-1);
  }

  // SBM cluster allocation samples (only populated when difference prior is SBM).
  // Layout matches the bgm() side: (num_variables x iter), one column per draw.
  const bool is_sbm = difference_selection &&
                      (difference_prior_type == "Stochastic-Block");
  arma::imat allocation_samples;
  if (is_sbm) {
    allocation_samples.set_size(num_variables, iter);
    allocation_samples.fill(-1);
  }

  // For logging nuts performance
  arma::ivec treedepth_samples(iter, arma::fill::zeros);
  arma::ivec divergent_samples(iter, arma::fill::zeros);
  arma::vec energy_samples(iter, arma::fill::zeros);
  arma::vec accept_prob_samples(iter, arma::fill::zeros);

  // Edge update shuffling setup
  arma::uvec v = arma::regspace<arma::uvec>(0, num_pair - 1);
  arma::uvec order(num_pair);
  arma::imat index(num_pair, 3);

  // --- Initialize proposal SDs
  arma::mat proposal_sd_main(num_main, num_groups, arma::fill::ones);
  arma::mat proposal_sd_pair(num_pair, num_groups, arma::fill::ones);

  // --- Optional NUTS warmup stage
  double initial_step_size = 1.0;
  if (update_method == nuts) {
    initial_step_size = find_initial_stepsize_bgmcompare(
      main_effects, pairwise_effects, main_effect_indices,
      pairwise_effect_indices, inclusion_indicator, projection, num_categories,
      observations, num_groups, group_indices, counts_per_category,
      blume_capel_stats, pairwise_stats, is_ordinal_variable,
      baseline_category,
      target_accept, rng,
      interaction_prior, difference_prior, threshold_prior
    );
  }

  // --- Warmup scheduling + adaptation controller
  WarmupSchedule warmup_schedule(warmup, difference_selection, (update_method != adaptive_metropolis));

  NUTSAdaptationController nuts_adapt(
      (num_main + num_pair) * num_groups, initial_step_size, target_accept,
      warmup_schedule, learn_mass_matrix
  );

  MetropolisAdaptationController metropolis_adapt_main(
      proposal_sd_main, warmup_schedule, target_accept
  );
  MetropolisAdaptationController metropolis_adapt_pair(
      proposal_sd_pair, warmup_schedule, target_accept
  );

  const int total_iter = warmup_schedule.total_warmup + iter;

  // --- Main Gibbs sampling loop
  bool userInterrupt = false;
  for (int iteration = 0; iteration < total_iter; iteration++) {

    pm.update(chain_id - 1);
    if (pm.shouldExit()) {
      userInterrupt = true;
      break;
    }

    // Shuffle update order of edge indices
    order = arma_randperm(rng, num_pair);
    for (int i = 0; i < num_pair; i++) {
      index.row(i) = interaction_index_matrix.row(order(i));
    }

    // Optional imputation
    if (na_impute) {
      impute_missing_bgmcompare (
          main_effects, pairwise_effects, main_effect_indices,
          pairwise_effect_indices, inclusion_indicator, projection,
          observations, num_groups, group_membership, group_indices,
          counts_per_category, blume_capel_stats, pairwise_stats,
          num_categories, missing_data_indices, is_ordinal_variable,
          baseline_category, rng
      );
    }

    // Main Gibbs update step for parameters
    gibbs_update_step_bgmcompare (
        observations, num_categories, counts_per_category,
        blume_capel_stats, inclusion_indicator,
        pairwise_effects, main_effects, is_ordinal_variable, baseline_category,
        iteration, pairwise_effect_indices, pairwise_stats, nuts_max_depth,
        nuts_adapt, metropolis_adapt_main, metropolis_adapt_pair, learn_mass_matrix,
        warmup_schedule, treedepth_samples, divergent_samples, energy_samples,
        accept_prob_samples, main_effect_indices, projection, num_groups, group_indices,
        rng, inclusion_probability,
        update_method, proposal_sd_main, proposal_sd_pair, index,
        main_difference_selection,
        interaction_prior, difference_prior, threshold_prior
    );

    // --- Update difference probabilities under the prior (if difference selection is active)
    if (warmup_schedule.selection_enabled(iteration)) {
      if (difference_prior_type == "Beta-Bernoulli") {
        // Pool pairwise + (optional) main-effect indicators into a single
        // shared inclusion probability (legacy compare semantics).
        int sumG = 0;
        for (int i = 0; i < num_variables - 1; ++i) {
          for (int j = i + 1; j < num_variables; ++j) {
            sumG += inclusion_indicator(i, j);
          }
        }
        int num_main_selectable = 0;
        if (main_difference_selection) {
          for(int i = 0; i < num_variables; i++) {
            sumG += inclusion_indicator(i, i);
          }
          num_main_selectable = num_variables;
        }
        double prob = rbeta(rng,
                            difference_selection_alpha + sumG,
                            difference_selection_beta + num_pair + num_main_selectable - sumG);
        std::fill(inclusion_probability.begin(), inclusion_probability.end(), prob);
      } else if (difference_prior_type == "Stochastic-Block") {
        // Off-diagonal (pairwise differences): MFM-SBM block-allocation +
        // block-probability update. Mirrors the bgm() SBM path.
        difference_edge_prior.update(
          inclusion_indicator, inclusion_probability,
          num_variables, num_pair, rng
        );
        // Diagonal (main-effect differences): independent Beta-Bernoulli
        // conjugate update using the SBM's within-cluster (alpha, beta).
        if (main_difference_selection) {
          int sumD = 0;
          for (int i = 0; i < num_variables; ++i) {
            sumD += inclusion_indicator(i, i);
          }
          double prob_main = rbeta(rng,
                                   difference_selection_alpha + sumD,
                                   difference_selection_beta + num_variables - sumD);
          for (int i = 0; i < num_variables; ++i) {
            inclusion_probability(i, i) = prob_main;
          }
        }
      }
    }

    // --- Store states
    if (iteration >= warmup_schedule.total_warmup) {
      int sample_index = iteration - warmup_schedule.total_warmup;


      int cntr = 0;
      for (int col = 0; col < num_groups; ++col) {
        for (int row = 0; row < num_main; ++row) {
          main_effect_samples(sample_index, cntr) = main_effects(row, col);
          cntr++;
        }
      }

      cntr = 0;
      for (int col = 0; col < num_groups; ++col) {
        for (int row = 0; row < num_pair; ++row) {
          pairwise_effect_samples(sample_index, cntr) = pairwise_effects(row, col);
          cntr++;
        }
      }

      if (difference_selection) {
        int cntr = 0;
        for (int i = 0; i < num_variables; ++i) {
          for (int j = i; j < num_variables; ++j) {
            indicator_samples(sample_index, cntr) = inclusion_indicator(i, j);
            cntr++;
          }
        }
      }

      if (is_sbm && difference_edge_prior.has_allocations()) {
        allocation_samples.col(sample_index) =
          difference_edge_prior.get_allocations();
      }
    }
  }

  bgmCompareOutput out;
  out.chain_id = chain_id;
  out.main_samples = main_effect_samples;
  out.pairwise_samples = pairwise_effect_samples;
  out.treedepth_samples = treedepth_samples;
  out.divergent_samples = divergent_samples;
  out.energy_samples = energy_samples;
  out.accept_prob_samples = accept_prob_samples;
  out.has_indicator = difference_selection;
  if (difference_selection) {
    out.indicator_samples = indicator_samples;
  } else {
    out.indicator_samples = arma::imat();
  }
  out.has_allocations = is_sbm;
  if (is_sbm) {
    out.allocation_samples = allocation_samples;
  } else {
    out.allocation_samples = arma::imat();
  }
  out.userInterrupt = userInterrupt;

  return out;
}
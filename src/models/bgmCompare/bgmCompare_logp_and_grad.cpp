#include <RcppArmadillo.h>
#include "models/bgmCompare/bgmCompare_helper.h"
#include "models/bgmCompare/bgmCompare_logp_and_grad.h"
#include <cmath>
#include "math/explog_macros.h"
#include "utils/common_helpers.h"
#include "utils/variable_helpers.h"
#include "priors/parameter_prior.h"



// Compute the total length of the parameter vector in the bgmCompare model.
//
// The parameter vector consists of:
//  1. Main-effect overall parameters (column 0).
//  2. Pairwise-effect overall parameters (column 0).
//  3. Main-effect group-difference parameters (columns 1..G-1) for variables
//     with inclusion_indicator(v,v) == 1.
//  4. Pairwise-effect group-difference parameters (columns 1..G-1) for pairs
//     with inclusion_indicator(v1,v2) == 1.
//
// Inputs:
//  - num_variables: Number of observed variables.
//  - main_effect_indices: Row ranges [start,end] in main_effects for each variable.
//  - pairwise_effect_indices: Lookup table mapping (var1,var2) → row in pairwise_effects.
//  - inclusion_indicator: Symmetric binary matrix; diagonal entries control main-effect
//    differences, off-diagonal entries control pairwise-effect differences.
//  - num_categories: Vector of category counts per variable.
//  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel) for each variable.
//  - num_groups: Number of groups in the model.
//
// Returns:
//  - arma::uword: Total number of parameters in the vectorized model.
//
// Notes:
//  - This function must be consistent with vectorize_model_parameters_bgmcompare().
//  - Used to allocate gradient vectors, prior vectors, and mass matrices.
arma::uword total_length(
    const int num_variables,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::imat& inclusion_indicator,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable,
    const int num_groups
) {
  const int n_main_rows = count_num_main_effects(
    num_categories, is_ordinal_variable
  );
  const int n_pair_rows = num_variables * (num_variables - 1) / 2;

  arma::uword total_len = 0;
  total_len += n_main_rows;      // main col 0
  total_len += n_pair_rows;      // pair col 0
  for (int v = 0; v < num_variables; ++v) {
    if (inclusion_indicator(v, v) == 1) {
      const int r0 = main_effect_indices(v, 0);
      const int r1 = main_effect_indices(v, 1);
      total_len += static_cast<long long>(r1 - r0 + 1) * (num_groups - 1);
    }
  }
  for (int v2 = 0; v2 < num_variables - 1; ++v2) {
    for (int v1 = v2 + 1; v1 < num_variables; ++v1) {
      total_len += (inclusion_indicator(v1, v2) == 1) * (num_groups - 1);
    }
  }

  return total_len;
}



// Compute the observed-data contribution to the gradient vector
// in the bgmCompare model (active parameterization).
//
// This function accumulates observed sufficient statistics from the data
// and projects them into the parameter vector space. The output has the
// same length and ordering as `vectorize_model_parameters_bgmcompare()`,
// and includes:
//  1. Main-effect overall parameters (column 0).
//  2. Pairwise-effect overall parameters (column 0).
//  3. Main-effect group-difference parameters (columns 1..G-1) if
//     inclusion_indicator(v,v) == 1.
//  4. Pairwise-effect group-difference parameters (columns 1..G-1) if
//     inclusion_indicator(v1,v2) == 1.
//
// Inputs:
//  - main_effect_indices: Row ranges [start,end] in main_effects for each variable.
//  - pairwise_effect_indices: Lookup table mapping (var1,var2) → row in pairwise_effects.
//  - projection: Matrix of size (num_groups × (num_groups-1)) containing group projections.
//  - observations: Matrix of observed variable values (N × V).
//  - group_indices: Index ranges [start,end] defining which rows in `observations`
//    belong to each group.
//  - num_categories: Vector giving the number of categories per variable.
//  - inclusion_indicator: Symmetric binary matrix; diagonal entries control inclusion
//    of main-effect differences, off-diagonal entries control inclusion of pairwise
//    differences.
//  - counts_per_category_group: Per-group category count tables (list of matrices).
//  - blume_capel_stats_group: Per-group Blume–Capel sufficient statistics (list of matrices).
//  - pairwise_stats_group: Per-group pairwise sufficient statistics (list of matrices).
//  - num_groups: Number of groups.
//  - is_ordinal_variable: Indicator vector (1 = ordinal, 0 = Blume–Capel) per variable.
//  - baseline_category: Vector of baseline categories per variable (Blume–Capel).
//  - main_index: Index map for main effects (from build_index_maps()).
//  - pair_index: Index map for pairwise effects (from build_index_maps()).
//
// Returns:
//  - arma::vec: Observed-data contribution to the gradient (length = total_length()).
//
// Notes:
//  - This function computes the *data-dependent* part of the gradient only;
//    parameter-dependent expected statistics and priors must be added separately.
//  - The output ordering must remain consistent with `vectorize_model_parameters_bgmcompare()`.
arma::vec gradient_observed_active(
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::mat& projection,
    const arma::imat& observations,
    const arma::imat& group_indices,
    const arma::ivec& num_categories,
    const arma::imat& inclusion_indicator,
    const std::vector<arma::imat>& counts_per_category_group,
    const std::vector<arma::imat>& blume_capel_stats_group,
    const std::vector<arma::mat>&  pairwise_stats_group,
    const int num_groups,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const arma::imat main_index,
    const arma::imat pair_index
) {
  const int num_variables = observations.n_cols;
  arma::uword total_len = total_length(
    num_variables, main_effect_indices, pairwise_effect_indices,
    inclusion_indicator, num_categories, is_ordinal_variable, num_groups
  );

  arma::vec grad_obs(total_len, arma::fill::zeros);
  int off;

  // -------------------------------
  // Observed sufficient statistics
  // -------------------------------
  for (int g = 0; g < num_groups; g++) {
    // list access
    const arma::imat& counts_per_category = counts_per_category_group[g];
    const arma::imat& blume_capel_stats = blume_capel_stats_group[g];
    const arma::vec proj_g = projection.row(g).t(); // length = num_groups-1

    // Main effects
    for (int v = 0; v < num_variables; v++) {
      const int base = main_effect_indices(v, 0);
      const int num_cats = num_categories(v);

      if (is_ordinal_variable(v)) {
        for (int c = 0; c < num_cats; c++) {
          const int count = counts_per_category(c, v);
          // overall
          off = main_index(base + c, 0);
          grad_obs(off) += count;

          // diffs
          if(inclusion_indicator(v, v) != 0) {
            for (int k = 1; k < num_groups; k++) {
              off = main_index(base + c, k);
              grad_obs(off) += count * proj_g(k-1);
            }
          }
        }
      } else {
        const int bc_0 = blume_capel_stats(0, v);
        const int bc_1 = blume_capel_stats(1, v);

        // overall (2 stats)
        off = main_index(base, 0);
        grad_obs(off) += bc_0;

        off = main_index(base + 1, 0);
        grad_obs(off) += bc_1;

        // diffs
        if(inclusion_indicator(v, v) != 0) {
          for (int k = 1; k < num_groups; k++) {
            off = main_index(base, k);
            grad_obs(off) += bc_0 * proj_g(k-1);

            off = main_index(base + 1, k);
            grad_obs(off) += bc_1 * proj_g(k-1);
          }
        }
      }
    }

    // Pairwise (observed)
    const arma::mat& pairwise_stats = pairwise_stats_group[g];
    for (int v1 = 0; v1 < num_variables - 1; v1++) {
      for (int v2 = v1 + 1; v2 < num_variables; v2++) {
        const int row = pairwise_effect_indices(v1, v2);
        const double pw_stats = 2.0 * pairwise_stats(v1, v2);

        off = pair_index(row, 0);
        grad_obs(off) += pw_stats; // upper tri counted once

        if(inclusion_indicator(v1, v2) != 0){
          for (int k = 1; k < num_groups; k++) {
            off = pair_index(row, k);
            grad_obs(off) += pw_stats * proj_g(k-1);
          }
        }
      }
    }
  }

  return grad_obs;
}



// Computes the gradient of the log pseudoposterior for the bgmCompare model.
//
// The gradient combines three contributions:
//  1. Observed sufficient statistics (precomputed and supplied via `grad_obs`).
//  2. Expected sufficient statistics under the current parameter values
//     (computed using softmax probabilities for ordinal or Blume–Capel variables).
//  3. Prior contributions on main effects, pairwise effects, and group differences.
//
// Procedure:
//  - Initialize gradient with `grad_obs` (observed-data contribution).
//  - Loop over groups:
//    * Build group-specific main and pairwise effects using
//      `compute_group_main_effects()` and `compute_group_pairwise_effects()`.
//    * Compute expected sufficient statistics from residual scores and
//      subtract them from the gradient.
//  - Add prior contributions:
//    * Logistic–Beta prior gradient for main-effect baseline parameters.
//    * Parameter-prior gradient for group-difference parameters and pairwise effects.
//
// Inputs:
//  - main_effects: Matrix of main-effect parameters (rows = categories, cols = groups).
//  - pairwise_effects: Matrix of pairwise-effect parameters (rows = pairs, cols = groups).
//  - main_effect_indices: Index ranges [row_start,row_end] for each variable in main_effects.
//  - pairwise_effect_indices: Lookup table mapping (var1,var2) → row in pairwise_effects.
//  - projection: Group projection matrix (num_groups × (num_groups − 1)).
//  - observations_double: Observation matrix (N × V), pre-converted to double.
//  - group_indices: Row ranges [start,end] for each group in `observations_double`.
//  - num_categories: Number of categories per variable.
//  - counts_per_category_group: Per-group category counts (ordinal variables).
//  - blume_capel_stats_group: Per-group sufficient statistics (Blume–Capel variables).
//  - pairwise_stats_group: Per-group pairwise sufficient statistics.
//  - num_groups: Number of groups.
//  - inclusion_indicator: Symmetric binary matrix; diagonal entries control inclusion
//    of main-effect differences, off-diagonal entries control inclusion of pairwise
//    differences.
//  - is_ordinal_variable: Indicator vector (1 = ordinal, 0 = Blume–Capel).
//  - baseline_category: Reference categories for Blume–Capel variables.
//  - main_alpha, main_beta: Hyperparameters for Beta priors on main effects.
//  - interaction_scale: Scale parameter of the prior on baseline pairwise effects.
//  - difference_scale: Scale parameter of the difference prior.
//  - main_index: Index map for main-effect parameters (from build_index_maps()).
//  - pair_index: Index map for pairwise-effect parameters (from build_index_maps()).
//  - grad_obs: Precomputed observed-data contribution to the gradient
//    (output of `gradient_observed_active()`).
//
// Returns:
//  - arma::vec: Gradient of the log pseudoposterior with respect to all active
//    parameters, in the layout defined by `vectorize_model_parameters_bgmcompare()`.
//
// Notes:
//  - Must remain consistent with `vectorize_model_parameters_bgmcompare()` and
//    `unvectorize_model_parameters_bgmcompare()`.
//  - Expected sufficient statistics are computed on-the-fly, while observed
//    statistics are passed in via `grad_obs`.
//  - Priors are applied after observed and expected contributions.
arma::vec gradient(
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::mat& projection,
    const arma::mat& observations_double,
    const arma::imat& group_indices,
    const arma::ivec& num_categories,
    const std::vector<arma::imat>& counts_per_category_group,
    const std::vector<arma::imat>& blume_capel_stats_group,
    const std::vector<arma::mat>&  pairwise_stats_group,
    const int num_groups,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const arma::imat& main_index,
    const arma::imat& pair_index,
    const arma::vec& grad_obs,
    const BaseParameterPrior& interaction_prior,
    const BaseParameterPrior& difference_prior,
    const BaseParameterPrior& threshold_prior
) {
  const int num_variables  = observations_double.n_cols;
  const int max_num_categories = num_categories.max();

  arma::vec grad = grad_obs;

  int off;

  arma::mat main_group(num_variables, max_num_categories, arma::fill::none);
  arma::mat pairwise_group(num_variables, num_variables, arma::fill::none);

  for (int g = 0; g < num_groups; ++g) {
    const int r0 = group_indices(g, 0);
    const int r1 = group_indices(g, 1);

    const arma::vec proj_g = projection.row(g).t();
    main_group.zeros();
    pairwise_group.zeros();

    for (int v = 0; v < num_variables; ++v) {
      arma::vec me = compute_group_main_effects(
        v, num_groups, main_effects, main_effect_indices, proj_g
      );
      main_group(v, arma::span(0, me.n_elem - 1)) = me.t();

      for (int u = v + 1; u < num_variables; ++u) {
        double w = compute_group_pairwise_effects(
          v, u, num_groups, pairwise_effects, pairwise_effect_indices,
          inclusion_indicator, proj_g
        );
        pairwise_group(v, u) = w;
        pairwise_group(u, v) = w;
      }
    }

    // observations_double is already arma::mat, no conversion needed
    const arma::mat obs = observations_double.rows(r0, r1);
    const arma::mat obs_t = obs.t();  // Pre-transpose for BLAS vectorization
    const arma::mat residual_matrix = obs * pairwise_group;

    for (int v = 0; v < num_variables; ++v) {
      const int K = num_categories(v);
      const int ref = baseline_category(v);

      arma::vec rest_score = residual_matrix.col(v);
      arma::vec bound = K * rest_score;

      arma::mat probs;
      if (is_ordinal_variable(v)) {
        arma::vec main_param = main_group.row(v).cols(0, K - 1).t();
        probs = compute_probs_ordinal(
          main_param, rest_score, bound, K
        );
      } else {
        const double lin_effect = main_group(v, 0);
        const double quad_effect = main_group(v, 1);
        probs = compute_probs_blume_capel(
          rest_score, lin_effect, quad_effect, ref, K, bound
        );
      }

      // ---- MAIN expected ----
      const int base = main_effect_indices(v, 0);
      if (is_ordinal_variable(v)) {
        for (int s = 1; s <= K; s++) {
          const int j = s - 1;
          double sum_col_s = arma::accu(probs.col(s));

          off = main_index(base + j, 0);
          grad(off) -= sum_col_s;

          if (inclusion_indicator(v, v) == 0) continue;
          for (int k = 1; k < num_groups; k++) {
            off = main_index(base + j, k);
            grad(off) -= proj_g(k - 1) * sum_col_s;
          }
        }
      } else {
        arma::vec lin_score  = arma::regspace<arma::vec>(0 - ref, K - ref);
        arma::vec quad_score = arma::square(lin_score);

        double sum_lin  = arma::accu(probs * lin_score);
        double sum_quad = arma::accu(probs * quad_score);

        off = main_index(base, 0);
        grad(off) -= sum_lin;
        off = main_index(base + 1, 0);
        grad(off) -= sum_quad;

        if (inclusion_indicator(v, v) == 0) continue;
        for (int k = 1; k < num_groups; k++) {
          off = main_index(base, k);
          grad(off) -= proj_g(k - 1) * sum_lin;
          off = main_index(base + 1, k);
          grad(off) -= proj_g(k - 1) * sum_quad;
        }
      }

      // ---- PAIRWISE expected (BLAS vectorized) ----
      // Compute expected score E[i] = sum_{s} s * P(X=s|rest) for each observation i
      arma::vec E;
      if (is_ordinal_variable(v)) {
        arma::vec weights = arma::regspace<arma::vec>(1, K);
        E = probs.cols(1, K) * weights;
      } else {
        arma::vec score = arma::regspace<arma::vec>(0 - ref, K - ref);
        E = probs * score;
      }

      // Compute all pairwise contributions with BLAS matrix-vector multiplication
      arma::vec pw_grad = obs_t * E;

      for (int v2 = 0; v2 < num_variables; v2++) {
        if (v == v2) continue;

        double sum_expectation = pw_grad(v2);

        const int row = (v < v2) ? pairwise_effect_indices(v, v2)
          : pairwise_effect_indices(v2, v);

        off = pair_index(row, 0);
        grad(off) -= sum_expectation;

        if (inclusion_indicator(v, v2) == 0) continue;
        for (int k = 1; k < num_groups; k++) {
          off = pair_index(row, k);
          grad(off) -= proj_g(k - 1) * sum_expectation;

        }
      }
    }
  }

  // -------------------------------
  // Priors
  // -------------------------------
  // Main
  for (int v = 0; v < num_variables; ++v) {
    const int base     = main_effect_indices(v, 0);
    const int num_cats = num_categories(v);

    if (is_ordinal_variable(v)) {
      for (int c = 0; c < num_cats; ++c) {
        off = main_index(base + c, 0);
        double value = main_effects(base + c, 0);
        grad(off) += threshold_prior.grad(value);

        if (inclusion_indicator(v, v) == 0) continue;
        for (int k = 1; k < num_groups; ++k) {
          off = main_index(base + c, k);
          double value = main_effects(base + c, k);
          grad(off) += difference_prior.grad(value);
        }

      }
    } else {
      off = main_index(base, 0);
      double value = main_effects(base, 0);
      grad(off) += threshold_prior.grad(value);

      off = main_index(base + 1, 0);
      value = main_effects(base + 1, 0);
      grad(off) += threshold_prior.grad(value);


      if (inclusion_indicator(v, v) == 0) continue;
      for (int k = 1; k < num_groups; ++k) {
        off = main_index(base, k);
        double value = main_effects(base, k);
        grad(off) += difference_prior.grad(value);

        off = main_index(base + 1, k);
        value = main_effects(base + 1, k);
        grad(off) += difference_prior.grad(value);
      }
    }
  }

  // Pairwise
  for (int v1 = 0; v1 < num_variables - 1; ++v1) {
    for (int v2 = v1 + 1; v2 < num_variables; ++v2) {
      const int row = pairwise_effect_indices(v1, v2);

      off = pair_index(row, 0);
      double value = pairwise_effects(row, 0);
      grad(off) += interaction_prior.grad(value);


      if (inclusion_indicator(v1, v2) == 0) continue;
      for (int k = 1; k < num_groups; ++k) {
        off = pair_index(row, k);
        double value = pairwise_effects(row, k);
        grad(off) += difference_prior.grad(value);
      }
    }
  }

  return grad;
}


// Computes both log pseudoposterior and gradient in a single pass.
//
// Fuses the computations of `log_pseudoposterior()` and `gradient()`,
// sharing intermediate results (group-specific effects, residual matrices,
// and probability computations) to avoid redundant work during NUTS sampling.
//
// Returns:
//  - std::pair containing:
//    - first: log pseudoposterior value (scalar)
//    - second: gradient vector (same layout as gradient())
std::pair<double, arma::vec> logp_and_gradient(
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::mat& projection,
    const arma::mat& observations_double,
    const arma::imat& group_indices,
    const arma::ivec& num_categories,
    const std::vector<arma::imat>& counts_per_category_group,
    const std::vector<arma::imat>& blume_capel_stats_group,
    const std::vector<arma::mat>&  pairwise_stats_group,
    const int num_groups,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const arma::imat& main_index,
    const arma::imat& pair_index,
    const arma::vec& grad_obs,
    const BaseParameterPrior& interaction_prior,
    const BaseParameterPrior& difference_prior,
    const BaseParameterPrior& threshold_prior
) {
  const int num_variables  = observations_double.n_cols;
  const int max_num_categories = num_categories.max();

  double log_pp = 0.0;
  arma::vec grad = grad_obs;

  int off;

  arma::mat main_group(num_variables, max_num_categories, arma::fill::none);
  arma::mat pairwise_group(num_variables, num_variables, arma::fill::none);

  // --- per group ---
  for (int g = 0; g < num_groups; ++g) {
    const int r0 = group_indices(g, 0);
    const int r1 = group_indices(g, 1);

    const arma::imat& counts_per_category = counts_per_category_group[g];
    const arma::imat& blume_capel_stats = blume_capel_stats_group[g];
    const arma::vec proj_g = projection.row(g).t();

    main_group.zeros();
    pairwise_group.zeros();

    // ---- build group-specific main & pairwise effects ----
    for (int v = 0; v < num_variables; v++) {
      arma::vec me = compute_group_main_effects(
        v, num_groups, main_effects, main_effect_indices, proj_g
      );
      main_group(v, arma::span(0, me.n_elem - 1)) = me.t();

      for (int u = v + 1; u < num_variables; ++u) {
        double w = compute_group_pairwise_effects(
          v, u, num_groups, pairwise_effects, pairwise_effect_indices,
          inclusion_indicator, proj_g
        );
        pairwise_group(v, u) = w;
        pairwise_group(u, v) = w;
      }

      // ---- data contribution pseudolikelihood (linear terms) ----
      const int num_cats = num_categories(v);
      if (is_ordinal_variable(v)) {
        for (int c = 0; c < num_cats; c++) {
          const double val = main_group(v, c);
          log_pp += static_cast<double>(counts_per_category(c, v)) * val;
        }
      } else {
        log_pp += static_cast<double>(blume_capel_stats(0, v)) * main_group(v, 0);
        log_pp += static_cast<double>(blume_capel_stats(1, v)) * main_group(v, 1);
      }
    }

    // ---- data contribution pseudolikelihood (quadratic terms) ----
    const arma::mat obs = observations_double.rows(r0, r1);
    const arma::mat obs_t = obs.t();  // Pre-transpose for BLAS vectorization
    const arma::mat& pairwise_stats = pairwise_stats_group[g];

    log_pp += arma::accu(pairwise_group % pairwise_stats);

    // ---- pseudolikelihood normalizing constants & gradient (per variable) ----
    const arma::mat residual_matrix = obs * pairwise_group;

    for (int v = 0; v < num_variables; ++v) {
      const int K = num_categories(v);
      const int ref = baseline_category(v);
      const int base = main_effect_indices(v, 0);

      const arma::vec rest_score = residual_matrix.col(v);
      arma::vec bound = K * rest_score;

      // Persistent per-thread scratch — survives across NUTS leapfrog calls.
      // Each TBB chain worker has its own thread, so each chain has its own
      // scratch (no contention; isolated state).
      thread_local LogZAndProbs result;
      thread_local LogZScratch  scratch;
      if (is_ordinal_variable(v)) {
        arma::vec main_param = main_group.row(v).cols(0, K - 1).t();
        compute_logZ_and_probs_ordinal_into(
          main_param, rest_score, bound, K, result, scratch
        );
      } else {
        const double lin_effect = main_group(v, 0);
        const double quad_effect = main_group(v, 1);
        compute_logZ_and_probs_blume_capel_into(
          rest_score, lin_effect, quad_effect, ref, K, bound, result, scratch
        );
      }

      // log_pp contribution: subtract log normalizers
      log_pp -= arma::accu(result.log_Z);

      // ---- gradient: MAIN expected ----
      const arma::mat& probs = result.probs;

      if (is_ordinal_variable(v)) {
        for (int s = 1; s <= K; s++) {
          const int j = s - 1;
          double sum_col_s = arma::accu(probs.col(s));

          off = main_index(base + j, 0);
          grad(off) -= sum_col_s;

          if (inclusion_indicator(v, v) == 0) continue;
          for (int k = 1; k < num_groups; k++) {
            off = main_index(base + j, k);
            grad(off) -= proj_g(k - 1) * sum_col_s;
          }
        }
      } else {
        arma::vec lin_score  = arma::regspace<arma::vec>(0 - ref, K - ref);
        arma::vec quad_score = arma::square(lin_score);

        double sum_lin  = arma::accu(probs * lin_score);
        double sum_quad = arma::accu(probs * quad_score);

        off = main_index(base, 0);
        grad(off) -= sum_lin;
        off = main_index(base + 1, 0);
        grad(off) -= sum_quad;

        if (inclusion_indicator(v, v) != 0) {
          for (int k = 1; k < num_groups; k++) {
            off = main_index(base, k);
            grad(off) -= proj_g(k - 1) * sum_lin;
            off = main_index(base + 1, k);
            grad(off) -= proj_g(k - 1) * sum_quad;
          }
        }
      }

      // ---- gradient: PAIRWISE expected (BLAS vectorized) ----
      arma::vec E;
      if (is_ordinal_variable(v)) {
        arma::vec weights = arma::regspace<arma::vec>(1, K);
        E = probs.cols(1, K) * weights;
      } else {
        arma::vec score = arma::regspace<arma::vec>(0 - ref, K - ref);
        E = probs * score;
      }

      arma::vec pw_grad = obs_t * E;

      for (int v2 = 0; v2 < num_variables; v2++) {
        if (v == v2) continue;

        double sum_expectation = pw_grad(v2);

        const int row = (v < v2) ? pairwise_effect_indices(v, v2)
          : pairwise_effect_indices(v2, v);

        off = pair_index(row, 0);
        grad(off) -= sum_expectation;

        if (inclusion_indicator(v, v2) == 0) continue;
        for (int k = 1; k < num_groups; k++) {
          off = pair_index(row, k);
          grad(off) -= proj_g(k - 1) * sum_expectation;
        }
      }
    }
  }

  // -------------------------------
  // Priors (same as in log_pseudoposterior and gradient)
  // -------------------------------

  // Main effects prior
  for (int v = 0; v < num_variables; v++) {
    const int r0 = main_effect_indices(v, 0);
    const int num_cats = num_categories(v);
    const int base = main_effect_indices(v, 0);

    if (is_ordinal_variable(v)) {
      for (int c = 0; c < num_cats; ++c) {
        double value = main_effects(r0 + c, 0);
        log_pp += threshold_prior.logp(value);

        off = main_index(base + c, 0);
        grad(off) += threshold_prior.grad(value);

        if (inclusion_indicator(v, v) == 0) continue;
        for (int eff = 1; eff < num_groups; eff++) {
          double diff_val = main_effects(r0 + c, eff);
          log_pp += difference_prior.logp(diff_val);

          off = main_index(base + c, eff);
          grad(off) += difference_prior.grad(diff_val);
        }
      }
    } else {
      for (int par = 0; par < 2; ++par) {
        double value = main_effects(r0 + par, 0);
        log_pp += threshold_prior.logp(value);

        off = main_index(base + par, 0);
        grad(off) += threshold_prior.grad(value);

        if (inclusion_indicator(v, v) == 0) continue;
        for (int eff = 1; eff < num_groups; eff++) {
          double diff_val = main_effects(r0 + par, eff);
          log_pp += difference_prior.logp(diff_val);

          off = main_index(base + par, eff);
          grad(off) += difference_prior.grad(diff_val);
        }
      }
    }
  }

  // Pairwise effects prior
  for (int v1 = 0; v1 < num_variables - 1; v1++) {
    for (int v2 = v1 + 1; v2 < num_variables; v2++) {
      const int idx = pairwise_effect_indices(v1, v2);

      double value = pairwise_effects(idx, 0);
      log_pp += interaction_prior.logp(value);

      off = pair_index(idx, 0);
      grad(off) += interaction_prior.grad(value);

      if (inclusion_indicator(v1, v2) == 0) continue;
      for (int eff = 1; eff < num_groups; eff++) {
        double diff_val = pairwise_effects(idx, eff);
        log_pp += difference_prior.logp(diff_val);

        off = pair_index(idx, eff);
        grad(off) += difference_prior.grad(diff_val);
      }
    }
  }

  return {log_pp, grad};
}



// Computes the log pseudoposterior contribution of a single main-effect parameter (bgmCompare model).
//
// This function isolates the contribution of one main-effect parameter,
// either the overall (baseline) effect or one of its group-specific differences.
//
// Procedure:
//  - For each group:
//    * Construct group-specific main effects for the selected variable
//      with `compute_group_main_effects()`.
//    * Construct group-specific pairwise effects for the variable.
//    * Add linear contributions from sufficient statistics.
//    * Subtract log normalizing constants from the group-specific likelihood.
//  - Add prior contribution:
//    * Logistic–Beta prior for baseline (h == 0).
//    * Difference prior for group differences (h > 0), if included.
//
// Inputs:
//  - main_effects: Matrix of main-effect parameters (rows = categories, cols = groups).
//  - pairwise_effects: Matrix of pairwise-effect parameters (rows = pairs, cols = groups).
//  - main_effect_indices: Index ranges [row_start,row_end] for each variable.
//  - pairwise_effect_indices: Lookup table mapping (var1,var2) → row in pairwise_effects.
//  - projection: Group projection matrix (num_groups × (num_groups − 1)).
//  - observations: Observation matrix (persons × variables).
//  - group_indices: Row ranges [start,end] for each group in observations.
//  - num_categories: Number of categories per variable.
//  - counts_per_category_group: Per-group category counts (for ordinal variables).
//  - blume_capel_stats_group: Per-group sufficient statistics (for Blume–Capel variables).
//  - num_groups: Number of groups.
//  - inclusion_indicator: Symmetric binary matrix of active variables (diag) and pairs (off-diag).
//  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
//  - baseline_category: Reference categories for Blume–Capel variables.
//  - main_alpha, main_beta: Hyperparameters for Beta priors on main effects.
//  - difference_scale: Scale parameter of the difference prior.
//  - variable: Index of the variable of interest.
//  - category: Category index (only used if variable is ordinal).
//  - par: Parameter index (0 = linear, 1 = quadratic; used for Blume–Capel).
//  - h: Column index (0 = overall baseline, >0 = group difference).
//
// Returns:
//  - The scalar log pseudoposterior contribution of the selected parameter.
//
// Notes:
//  - If h > 0 but inclusion_indicator(variable, variable) == 0,
//    the function returns 0.0 (no contribution).
//  - This component function is used in parameter-wise Metropolis updates.
//  - Consistent with the full `log_pseudoposterior()` for bgmCompare.
double log_pseudoposterior_main_component(
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::mat& projection,
    const arma::imat& observations,
    const arma::imat& group_indices,
    const arma::ivec& num_categories,
    const std::vector<arma::imat>& counts_per_category_group,
    const std::vector<arma::imat>& blume_capel_stats_group,
    const int num_groups,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    int variable,
    int category, // for ordinal variables only
    int par, // for Blume-Capel variables only
    int h, // Overall = 0, differences are 1,2,...
    const BaseParameterPrior& difference_prior,
    const BaseParameterPrior& threshold_prior
) {
  if(h > 0 && inclusion_indicator(variable, variable) == 0) {
    return 0.0; // No contribution if differences not included
  }

  const int num_variables = observations.n_cols;
  double log_pp = 0.0;

  // group-specific pairwise weights of `variable` with all other variables
  arma::vec pairwise_col(num_variables, arma::fill::zeros);

  // --- per group ---
  for (int group = 0; group < num_groups; ++group) {
    const arma::imat& counts_per_category = counts_per_category_group[group];
    const arma::imat& blume_capel_stats = blume_capel_stats_group[group];

    const arma::vec proj_g = projection.row(group).t(); // length = num_groups-1

    // ---- build group-specific main & pairwise effects ----
    arma::vec me = compute_group_main_effects(
      variable, num_groups, main_effects, main_effect_indices, proj_g
    );

    // pairwise weights with the other variables; entry `variable` stays zero
    for (int u = 0; u < num_variables; u++) {
      if(u == variable) continue;
      pairwise_col(u) = compute_group_pairwise_effects(
        variable, u, num_groups, pairwise_effects, pairwise_effect_indices,
        inclusion_indicator, proj_g
      );
    }

    // ---- data contribution pseudolikelihood (linear terms) ----
    if (is_ordinal_variable(variable)) {
      log_pp += static_cast<double>(counts_per_category(category, variable)) *
        me(category);
    } else {
      log_pp += static_cast<double>(blume_capel_stats(par, variable)) * me(par);
    }

    // ---- data contribution pseudolikelihood (quadratic terms) ----
    const int r0 = group_indices(group, 0);
    const int r1 = group_indices(group, 1);
    const arma::mat obs = arma::conv_to<arma::mat>::from(observations.rows(r0, r1));

    // ---- pseudolikelihood normalizing constants (per variable) ----
    const arma::vec rest_score = obs * pairwise_col;
    const int num_cats = num_categories(variable);

    // bound to stabilize exp; clamp at 0 (the reference-category exponent) so
    // exp cannot overflow. compute_denom_blume_capel overwrites bound with its
    // own max, so the clamp only affects the ordinal path.
    arma::vec bound = arma::clamp(num_cats * rest_score, 0.0, arma::datum::inf);
    arma::vec denom(rest_score.n_elem, arma::fill::zeros);

    if (is_ordinal_variable(variable)) {
      denom = compute_denom_ordinal(
        rest_score, me, bound
      );
    } else {
      const double lin_effect  = me(0);
      const double quad_effect = me(1);
      const int ref = baseline_category(variable);

      denom = compute_denom_blume_capel(
        rest_score, lin_effect, quad_effect, ref, num_cats,
        /*updated in place:*/bound
      );
    }

    // - sum_i [ bound_i + log denom_i ]
    log_pp -= arma::accu(bound + ARMA_MY_LOG(denom));
  }

  // ---- priors ----
  if (h == 0) {
    // Main effects prior (baseline)
    if(is_ordinal_variable(variable)) {
      int r = main_effect_indices(variable, 0) + category;
      log_pp += threshold_prior.logp(main_effects(r, 0));
    } else {
      int r = main_effect_indices(variable, 0) + par;
      log_pp += threshold_prior.logp(main_effects(r, 0));
    }
  } else {
    // Group-difference prior
    if(is_ordinal_variable(variable)) {
      int r = main_effect_indices(variable, 0) + category;
      log_pp += difference_prior.logp(main_effects(r, h));
    } else {
      int r = main_effect_indices(variable, 0) + par;
      log_pp += difference_prior.logp(main_effects(r, h));
    }
  }

  return log_pp;
}


// Computes the log pseudoposterior contribution of a single pairwise-effect parameter (bgmCompare model).
//
// Isolates the contribution of one pairwise-effect parameter between two variables,
// either the baseline effect (h == 0) or a group-specific difference (h > 0).
//
// Procedure:
//  - For each group:
//    * Construct group-specific main effects for the two variables.
//    * Add linear contributions from the pairwise sufficient statistic.
//      - Baseline (h == 0): contribution = 2 * suff_pair * proposed_value.
//      - Difference (h > 0): scaled by projection value proj_g(h-1).
//    * Subtract log normalizing constants from both variables' likelihoods.
//  - Add prior contribution:
//    * Interaction prior for baseline (scale = interaction_scale).
//    * Difference prior for group differences (scale = difference_scale).
//
// Inputs:
//  - main_effects: Matrix of main-effect parameters (rows = categories, cols = groups).
//  - pairwise_effects: Matrix of pairwise-effect parameters (rows = pairs, cols = groups).
//  - main_effect_indices: Index ranges [row_start, row_end] for each variable.
//  - pairwise_effect_indices: Lookup table mapping (var1, var2) to row in pairwise_effects.
//  - projection: Group projection matrix (num_groups × (num_groups − 1)).
//  - obs_double_groups: Per-group observation matrices converted to double.
//  - num_categories: Number of categories per variable.
//  - pairwise_stats_group: Per-group pairwise sufficient statistics.
//  - residual_matrices: Per-group residual matrices (persons × variables).
//  - num_groups: Number of groups.
//  - inclusion_indicator: Symmetric binary matrix of active variables (diag) and pairs (off-diag).
//  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
//  - baseline_category: Reference categories for Blume–Capel variables.
//  - interaction_prior: Prior (BaseParameterPrior) on baseline pairwise effects.
//  - difference_prior: Prior (BaseParameterPrior) on group differences.
//  - variable1, variable2: Indices of the variable pair.
//  - h: Column index (0 = baseline, > 0 = group difference).
//  - delta: Parameter change (proposed - current).
//
// Returns:
//  - The log pseudoposterior value at the proposed state.
//
// Notes:
//  - If h > 0 but inclusion_indicator(variable1, variable2) == 0, returns 0.0.
//  - The proposed value is computed as pairwise_effects(idx, h) + delta.
//  - Residual scores are adjusted by delta without modifying residual_matrices.
double log_pseudoposterior_pair_component(
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::mat& projection,
    const std::vector<arma::mat>& obs_double_groups,
    const arma::ivec& num_categories,
    const std::vector<arma::mat>& pairwise_stats_group,
    const std::vector<arma::mat>& residual_matrices,
    const int num_groups,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    int variable1,
    int variable2,
    int h,
    double delta,
    const BaseParameterPrior& interaction_prior,
    const BaseParameterPrior& difference_prior
) {
  if(h > 0 && inclusion_indicator(variable1, variable2) == 0) {
    return 0.0;
  }

  const int num_variables = obs_double_groups[0].n_cols;
  const int max_num_categories = num_categories.max();
  double log_pp = 0.0;
  int idx = pairwise_effect_indices(variable1, variable2);

  // Compute the proposed value
  const double proposed_value = pairwise_effects(idx, h) + delta;

  // --- per group ---
  for (int group = 0; group < num_groups; ++group) {
    const arma::vec proj_g = projection.row(group).t();

    // Compute group-specific delta: how much pairwise_group(var1,var2) changes for this group
    double delta_g = (h == 0) ? delta : delta * proj_g(h - 1);

    // ---- build group-specific main effects (only for the two variables) ----
    arma::mat main_group(num_variables, max_num_categories, arma::fill::zeros);
    for (int v : {variable1, variable2}) {
      arma::vec me = compute_group_main_effects(
        v, num_groups, main_effects, main_effect_indices, proj_g
      );
      main_group(v, arma::span(0, me.n_elem - 1)) = me.t();
    }

    // ---- data contribution pseudolikelihood ----
    const arma::mat& pairwise_stats = pairwise_stats_group[group];
    const double suff_pair = pairwise_stats(variable1, variable2);

    if(h == 0) {
      log_pp += 2.0 * suff_pair * proposed_value;
    } else {
      log_pp += 2.0 * suff_pair * proj_g(h-1) * proposed_value;
    }

    // ---- pseudolikelihood normalizing constants (using residual matrix + delta) ----
    const arma::mat& obs_g = obs_double_groups[group];

    for (int v : {variable1, variable2}) {
      const int num_cats = num_categories(v);
      const int other = (v == variable1) ? variable2 : variable1;

      // Use residual_matrix with delta adjustment: O(n) instead of O(n*p)
      arma::vec rest_score = residual_matrices[group].col(v) + obs_g.col(other) * delta_g;

      // bound to stabilize exp; clamp at 0 so exp cannot overflow (the
      // Blume-Capel branch overwrites bound with its own max).
      arma::vec bound = arma::clamp(num_cats * rest_score, 0.0, arma::datum::inf);
      arma::vec denom(rest_score.n_elem, arma::fill::zeros);

      if (is_ordinal_variable(v)) {
        arma::vec main_eff = main_group.row(v).cols(0, num_cats - 1).t();
        denom = compute_denom_ordinal(rest_score, main_eff, bound);
      } else {
        const double lin_effect  = main_group(v, 0);
        const double quad_effect = main_group(v, 1);
        const int ref = baseline_category(v);
        denom = compute_denom_blume_capel(rest_score, lin_effect, quad_effect, ref, num_cats, bound);
      }

      log_pp -= arma::accu(bound + ARMA_MY_LOG(denom));
    }
  }

  // ---- priors ----
  if (h == 0) {
    log_pp += interaction_prior.logp(proposed_value);
  } else {
    log_pp += difference_prior.logp(proposed_value);
  }
  return log_pp;
}



// Computes the log-ratio of pseudolikelihood normalizing constants
// for a single variable under current vs. proposed parameters (bgmCompare model).
//
// This function is used in Metropolis–Hastings updates for main-effect parameters.
// It evaluates how the normalizing constant (denominator of the pseudolikelihood)
// changes when switching from the current to the proposed parameter values.
//
// Procedure:
//  - For each group:
//    * Construct group-specific main effects (current vs. proposed).
//    * Construct group-specific pairwise weights for the variable.
//    * Compute residual scores for observations under both models.
//    * Calculate denominators with stability bounds (ordinal vs. Blume–Capel cases).
//    * Accumulate the log-ratio contribution across all observations.
//
// Inputs:
//  - current_main_effects, proposed_main_effects: Matrices of main-effect parameters
//    (rows = categories, cols = groups).
//  - current_pairwise_effects, proposed_pairwise_effects: Matrices of pairwise-effect parameters
//    (rows = pairs, cols = groups).
//  - main_effect_indices: Index ranges [row_start,row_end] for each variable.
//  - pairwise_effect_indices: Lookup table mapping (var1,var2) → row in pairwise_effects.
//  - projection: Group projection matrix (num_groups × (num_groups − 1)).
//  - observations: Observation matrix (persons × variables).
//  - group_indices: Row ranges [start,end] for each group in observations.
//  - num_categories: Number of categories per variable.
//  - num_groups: Number of groups.
//  - inclusion_indicator: Symmetric binary matrix of active variables (diag) and pairs (off-diag).
//  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
//  - baseline_category: Reference categories for Blume–Capel variables.
//  - variable: Index of the variable being updated.
//
// Returns:
//  - The scalar log-ratio of pseudolikelihood constants
//    (current model vs. proposed model).
//
// Notes:
//  - For ordinal variables, denominators include exp(-bound) and category terms.
//  - For Blume–Capel variables, denominators use linear/quadratic scores
//    with baseline centering.
//  - Stability bounds (`bound_current`, `bound_proposed`) are applied to avoid overflow.
double log_ratio_pseudolikelihood_constant_variable(
    const arma::mat& current_main_effects,
    const arma::mat& current_pairwise_effects,
    const arma::mat& proposed_main_effects,
    const arma::mat& proposed_pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::mat& projection,
    const arma::imat& observations,
    const arma::imat& group_indices,
    const arma::ivec& num_categories,
    const int num_groups,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const int variable
) {
  const int num_cats = num_categories(variable);
  const int num_variables = observations.n_cols;

  double log_ratio = 0.0;

  // --- per group ---
  for (int group = 0; group < num_groups; ++group) {
    const arma::vec proj_g = projection.row(group).t();

    // --- group-specific main effects (current/proposed) ---
    const arma::vec main_current = compute_group_main_effects(
      variable, num_groups, current_main_effects, main_effect_indices, proj_g
    );
    const arma::vec main_proposed = compute_group_main_effects(
      variable, num_groups, proposed_main_effects, main_effect_indices, proj_g
    );

    // --- group-specific pairwise effects for this variable (column) ---
    arma::vec weights_current(num_variables, arma::fill::zeros);
    arma::vec weights_proposed(num_variables, arma::fill::zeros);
    for (int u = 0; u < num_variables; ++u) {
      if (u == variable) continue;
      weights_current(u) = compute_group_pairwise_effects(
        variable, u, num_groups, current_pairwise_effects,
        pairwise_effect_indices, inclusion_indicator, proj_g
      );
      weights_proposed(u) = compute_group_pairwise_effects(
        variable, u, num_groups, proposed_pairwise_effects,
        pairwise_effect_indices, inclusion_indicator, proj_g
      );
    }

    // --- group observations and rest scores ---
    const int r0 = group_indices(group, 0);
    const int r1 = group_indices(group, 1);
    const arma::mat obs = arma::conv_to<arma::mat>::from(observations.rows(r0, r1));

    const arma::vec rest_current = obs * weights_current;
    const arma::vec rest_proposed = obs * weights_proposed;

    // --- denominators with stability bounds ---
    arma::vec bound_current;
    arma::vec bound_proposed;
    arma::vec denom_current(rest_current.n_elem, arma::fill::zeros);
    arma::vec denom_proposed(rest_proposed.n_elem, arma::fill::zeros);

    if (is_ordinal_variable (variable)) {
      bound_current = arma::clamp(rest_current * num_cats, 0.0, arma::datum::inf);
      bound_proposed = arma::clamp(rest_proposed * num_cats, 0.0, arma::datum::inf);

      denom_current += compute_denom_ordinal(
        rest_current, main_current, bound_current
      );
      denom_proposed += compute_denom_ordinal(
        rest_proposed, main_proposed, bound_proposed
      );
    } else {
      // Binary or categorical variable: linear + quadratic score
      const int ref_cat = baseline_category (variable);
      bound_current = rest_current * num_cats;
      bound_proposed = rest_proposed * num_cats;

      denom_current = compute_denom_blume_capel(
        rest_current, main_current (0), main_current (1),
        ref_cat, num_cats, /*Updated in place:*/bound_current
      );

      denom_proposed = compute_denom_blume_capel(
        rest_proposed, main_proposed (0), main_proposed (1),
        ref_cat, num_cats, /*Updated in place:*/bound_proposed
      );
    }

    // --- accumulate contribution ---
    log_ratio += arma::accu((bound_current - bound_proposed) +
      ARMA_MY_LOG(denom_current) - ARMA_MY_LOG(denom_proposed));
  }

  return log_ratio;
}



// Computes the log pseudolikelihood ratio for updating a single main-effect parameter (bgmCompare model).
//
// This function is used in Metropolis–Hastings updates for main effects.
// It compares the likelihood of the data under the current vs. proposed
// value of a single variable’s main-effect parameter, while keeping
// all other parameters fixed.
//
// Procedure:
//  - For each group:
//    * Compute group-specific main effects for the variable (current vs. proposed).
//    * Add contributions from observed sufficient statistics
//      (category counts or Blume–Capel stats).
//  - Add the ratio of pseudolikelihood normalizing constants by calling
//    `log_ratio_pseudolikelihood_constant_variable()`.
//
// Inputs:
//  - current_main_effects: Matrix of main-effect parameters (current state).
//  - proposed_main_effects: Matrix of main-effect parameters (candidate state).
//  - current_pairwise_effects: Matrix of pairwise-effect parameters (fixed at current state).
//  - main_effect_indices: Index ranges [row_start,row_end] for each variable.
//  - pairwise_effect_indices: Lookup table mapping (var1,var2) → row in pairwise_effects.
//  - projection: Group projection matrix (num_groups × (num_groups − 1)).
//  - observations: Observation matrix (persons × variables).
//  - group_indices: Row ranges [start,end] for each group in observations.
//  - num_categories: Number of categories per variable.
//  - counts_per_category_group: Per-group category counts (for ordinal variables).
//  - blume_capel_stats_group: Per-group sufficient statistics (for Blume–Capel variables).
//  - num_groups: Number of groups.
//  - inclusion_indicator: Symmetric binary matrix of active variables (diag) and pairs (off-diag).
//  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
//  - baseline_category: Reference categories for Blume–Capel variables.
//  - variable: Index of the variable being updated.
//
// Returns:
//  - The scalar log pseudolikelihood ratio (proposed vs. current).
//
// Notes:
//  - A temporary copy of `inclusion_indicator` is made to ensure the
//    variable’s self-term (diagonal entry) is included.
//  - Only the variable under update changes between current and proposed states;
//    all other variables and pairwise effects remain fixed.
//  - This function does not add prior contributions — only pseudolikelihood terms.
double log_pseudolikelihood_ratio_main(
    const arma::mat& current_main_effects,
    const arma::mat& proposed_main_effects,
    const arma::mat& current_pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::mat&  projection,
    const arma::imat& observations,
    const arma::imat& group_indices,
    const arma::ivec& num_categories,
    const std::vector<arma::imat>& counts_per_category_group,
    const std::vector<arma::imat>& blume_capel_stats_group,
    const int num_groups,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const int variable
) {
  double lr = 0.0;
  arma::imat tmp_ind = inclusion_indicator;
  tmp_ind(variable, variable) = 1; // Ensure self-interaction is included

  // Add data contribution (group-specific parameters via projection)
  for (int g = 0; g < num_groups; ++g) {
    const arma::vec proj_g = projection.row(g).t();

    const arma::vec main_cur  = compute_group_main_effects(
      variable, num_groups, current_main_effects,  main_effect_indices, proj_g
    );
    const arma::vec main_prop = compute_group_main_effects(
      variable, num_groups, proposed_main_effects, main_effect_indices, proj_g
    );

    if (is_ordinal_variable(variable)) {
      const arma::imat& num_obs = counts_per_category_group[g];
      const int num_cats = num_categories(variable);
      for (int c = 0; c < num_cats; ++c) {
        lr += (main_prop(c) - main_cur(c)) * static_cast<double>(num_obs(c, variable));
      }
    } else {
      const arma::imat& suff = blume_capel_stats_group[g];
      lr += (main_prop(0) - main_cur(0)) * static_cast<double>(suff(0, variable));
      lr += (main_prop(1) - main_cur(1)) * static_cast<double>(suff(1, variable));
    }
  }

  // Add ratio of normalizing constants
  lr += log_ratio_pseudolikelihood_constant_variable(
    current_main_effects, current_pairwise_effects, proposed_main_effects,
    /* same */ current_pairwise_effects, main_effect_indices,
    pairwise_effect_indices, projection, observations, group_indices,
    num_categories, num_groups, tmp_ind, is_ordinal_variable,
    baseline_category, variable
  );

  return lr;
}


// Computes the log pseudolikelihood ratio for updating a single pairwise-effect parameter (bgmCompare model).
//
// This function is used in Metropolis–Hastings updates for pairwise effects.
// It compares the likelihood of the data under the current vs. proposed
// value of a single interaction (var1,var2), while keeping all other
// parameters fixed.
//
// Procedure:
//  - Ensure the interaction is included in a temporary copy of inclusion_indicator.
//  - For each group:
//    * Compute group-specific pairwise effect for (var1,var2), current vs. proposed.
//    * Add linear contribution from the pairwise sufficient statistic.
//  - Add the ratio of pseudolikelihood normalizing constants for both variables:
//    * Call `log_ratio_pseudolikelihood_constant_variable()` separately for var1 and var2,
//      comparing current vs. proposed pairwise weights.
//
// Inputs:
//  - main_effects: Matrix of main-effect parameters (fixed).
//  - current_pairwise_effects: Matrix of pairwise-effect parameters (current state).
//  - proposed_pairwise_effects: Matrix of pairwise-effect parameters (candidate state).
//  - main_effect_indices: Index ranges [row_start,row_end] for each variable.
//  - pairwise_effect_indices: Lookup table mapping (var1,var2) → row in pairwise_effects.
//  - projection: Group projection matrix (num_groups × (num_groups − 1)).
//  - observations: Observation matrix (persons × variables).
//  - group_indices: Row ranges [start,end] for each group in observations.
//  - num_categories: Number of categories per variable.
//  - pairwise_stats_group: Per-group pairwise sufficient statistics.
//  - num_groups: Number of groups.
//  - inclusion_indicator: Symmetric binary matrix of active variables (diag) and pairs (off-diag).
//  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
//  - baseline_category: Reference categories for Blume–Capel variables.
//  - var1, var2: Indices of the variable pair being updated.
//
// Returns:
//  - The scalar log pseudolikelihood ratio (proposed vs. current).
//
// Notes:
//  - A temporary copy of `inclusion_indicator` is used to force the edge (var1,var2) as active.
//  - Only the selected pair changes between current and proposed states;
//    all other effects remain fixed.
//  - This function does not add prior contributions — only pseudolikelihood terms.
double log_pseudolikelihood_ratio_pairwise(
    const arma::mat& main_effects,
    const arma::mat& current_pairwise_effects,
    const arma::mat& proposed_pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::mat& projection,
    const arma::imat& observations,
    const arma::imat& group_indices,
    const arma::ivec& num_categories,
    const std::vector<arma::mat>& pairwise_stats_group,
    const int num_groups,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const int var1,
    const int var2
) {
  double lr = 0.0;
  // Ensure interaction is included
  arma::imat tmp_ind = inclusion_indicator;
  tmp_ind(var1, var2) = 1;
  tmp_ind(var2, var1) = 1;

  // Add data contribution
  for (int g = 0; g < num_groups; ++g) {
    const arma::vec proj_g = projection.row(g).t();
    const arma::mat& suff  = pairwise_stats_group[g];

    const double w_cur = compute_group_pairwise_effects(
      var1, var2, num_groups, current_pairwise_effects,
      pairwise_effect_indices, tmp_ind, proj_g
    );
    const double w_prop = compute_group_pairwise_effects(
      var1, var2, num_groups, proposed_pairwise_effects,
      pairwise_effect_indices, tmp_ind, proj_g
    );

    lr += 2.0 * (w_prop - w_cur) * suff(var1, var2);
  }

  // Add ratio of normalizing constant for `var1`
  lr += log_ratio_pseudolikelihood_constant_variable(
    main_effects, current_pairwise_effects, /* same */ main_effects,
    proposed_pairwise_effects, main_effect_indices, pairwise_effect_indices,
    projection, observations, group_indices, num_categories, num_groups,
    tmp_ind, is_ordinal_variable, baseline_category, var1
  );

  // Add ratio of normalizing constant for `var2`
  lr += log_ratio_pseudolikelihood_constant_variable(
    main_effects, current_pairwise_effects, /* same */ main_effects,
    proposed_pairwise_effects, main_effect_indices, pairwise_effect_indices,
    projection, observations, group_indices, num_categories, num_groups,
    tmp_ind, is_ordinal_variable, baseline_category, var2
  );

  return lr;
}
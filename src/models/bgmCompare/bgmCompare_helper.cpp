#include <RcppArmadillo.h>
#include <cmath>
#include "models/bgmCompare/bgmCompare_helper.h"
#include "models/bgmCompare/bgmCompare_state.h"
#include "utils/common_helpers.h"



// Converts the integer observations to double, whole-matrix and per group.
// Weight and residual members are left untouched; call
// rebuild_sweep_state_weights() afterwards to make the state consistent.
void initialize_sweep_state_observations(
    CompareSweepState& state,
    const arma::imat& observations,
    const arma::imat& group_indices,
    const int num_groups
) {
  state.obs_double_all = arma::conv_to<arma::mat>::from(observations);
  state.obs_double.resize(num_groups);
  for (int g = 0; g < num_groups; g++) {
    const int r0 = group_indices(g, 0);
    const int r1 = group_indices(g, 1);
    state.obs_double[g] = state.obs_double_all.rows(r0, r1);
  }

  const int num_variables = observations.n_cols;
  state.log_normalizer.zeros(num_variables, num_groups);
  state.normalizer_valid.zeros(num_variables);
}



// Recomputes the per-group effective pairwise weights from the current
// pairwise effects and inclusion indicators, and the residual matrices as
// one matrix product per group. Pairwise effects are stored on the
// association scale, so a rest score carries a factor two. Invalidates the
// normalizer cache.
void rebuild_sweep_state_weights(
    CompareSweepState& state,
    const arma::mat& pairwise_effects,
    const arma::imat& pairwise_effect_indices,
    const arma::imat& inclusion_indicator,
    const arma::mat& projection,
    const int num_groups
) {
  const int num_variables = inclusion_indicator.n_rows;
  state.pairwise_group.resize(num_groups);
  state.residual.resize(num_groups);

  for (int g = 0; g < num_groups; g++) {
    const arma::vec proj_g = projection.row(g).t();

    arma::mat& pairwise_g = state.pairwise_group[g];
    pairwise_g.zeros(num_variables, num_variables);
    for (int v = 0; v < num_variables - 1; v++) {
      for (int u = v + 1; u < num_variables; u++) {
        double w = compute_group_pairwise_effects(
          v, u, num_groups, pairwise_effects, pairwise_effect_indices,
          inclusion_indicator, proj_g
        );
        pairwise_g(v, u) = w;
        pairwise_g(u, v) = w;
      }
    }

    state.residual[g] = 2.0 * state.obs_double[g] * pairwise_g;
  }

  state.normalizer_valid.zeros();
}



// Computes group-specific main effects for a given variable (bgmCompare model).
//
// For a variable, the main-effect parameters are stored with:
//  - One "baseline" column (shared across groups).
//  - Additional columns for group-specific deviations.
//
// This function extracts the rows corresponding to the variable’s categories
// and combines them with the group projection vector to yield the
// group-specific main effects.
//
// Inputs:
//  - variable: Index of the variable of interest.
//  - num_groups: Number of groups in the model.
//  - main_effects: Matrix of main-effect parameters for all variables.
//  - main_effect_indices: Index matrix giving [start_row, end_row] for each variable.
//  - proj_group: Projection vector selecting the group (length = num_groups − 1).
//
// Returns:
//  - A vector of group-specific main effects for the categories of the variable.
//
// Notes:
//  - The projection vector should match the encoding used for group effects
//    (e.g. dummy or contrast coding).
//  - This function is used in likelihood evaluations where group-specific
//    parameters are required.
arma::vec compute_group_main_effects(
    const int variable,
    const int num_groups,
    const arma::mat& main_effects,
    const arma::imat& main_effect_indices,
    const arma::vec& proj_group
) {
  // Base index for accessing main effects for this variable
  int base_category_index = main_effect_indices(variable, 0);
  int last_category_index = main_effect_indices(variable, 1);

  arma::vec group_main_effects =
    arma::conv_to<arma::vec>::from(
      main_effects.rows(base_category_index, last_category_index).col(0));
  group_main_effects += main_effects.rows(
    base_category_index, last_category_index
  ).cols(
      1, num_groups - 1
  ) *
    proj_group;
  return group_main_effects;
}



// Computes the group-specific pairwise effect for a variable pair (bgmCompare model).
//
// For each variable pair, the pairwise-effect parameters are stored with:
//  - One "baseline" column (shared across groups).
//  - Additional columns for group-specific deviations.
//
// This function extracts the baseline effect and, if the edge is active
// (per the inclusion indicator), adds the group-specific deviation obtained
// from the projection vector.
//
// Inputs:
//  - var1, var2: Indices of the two variables forming the pair.
//  - num_groups: Number of groups in the model.
//  - pairwise_effects: Matrix of pairwise-effect parameters (rows = pairs).
//  - pairwise_effect_indices: Lookup matrix mapping (var1, var2) → row index.
//  - inclusion_indicator: Symmetric binary matrix of active edges.
//  - proj_group: Projection vector selecting the group (length = num_groups − 1).
//
// Returns:
//  - The group-specific pairwise effect for (var1, var2).
//
// Notes:
//  - The index matrix must match the storage convention (typically var1 < var2).
//  - If `inclusion_indicator(var1, var2) == 0`, only the baseline effect is used.
//  - This function is used in likelihood evaluations and Gibbs updates.
double compute_group_pairwise_effects(
    const int var1,
    const int var2,
    const int num_groups,
    const arma::mat& pairwise_effects,
    const arma::imat& pairwise_effect_indices,
    const arma::imat& inclusion_indicator,
    const arma::vec& proj_group
) {
  // Look up the row index for this pair; assume indices are provided for var1<var2
  // If your index matrix is symmetric, either order works; otherwise enforce (min,max).
  const int i = pairwise_effect_indices(var1, var2);

  // Shared/base effect
  double w = pairwise_effects(i, 0);

  // Optional group contrast
  if (inclusion_indicator(var1, var2) != 0) {
    // cols(1..G-1) * projection[group]^T  → scalar
    w += arma::as_scalar(pairwise_effects.row(i).cols(1, num_groups - 1) *
      proj_group);
  }
  return w;
}



// Flattens main-effect and pairwise-effect parameters into a single vector (bgmCompare model).
//
// Layout of the output vector:
//  1. Main-effect overall parameters (column 0 of main_effects), stacked by variable.
//  2. Pairwise-effect overall parameters (column 0 of pairwise_effects), stacked by pair.
//  3. Main-effect group differences (columns 1..G-1), included only if
//     the variable is marked active in inclusion_indicator(v,v).
//  4. Pairwise-effect group differences (columns 1..G-1), included only if
//     the pair is marked active in inclusion_indicator(v1,v2).
//
// Inputs:
//  - main_effects: Matrix of main-effect parameters (rows = categories, cols = groups).
//  - pairwise_effects: Matrix of pairwise-effect parameters (rows = pairs, cols = groups).
//  - inclusion_indicator: Symmetric binary matrix of active variables (diag) and pairs (off-diag).
//  - main_effect_indices: Index ranges [row_start, row_end] for each variable in main_effects.
//  - pairwise_effect_indices: Lookup table mapping (var1,var2) → row in pairwise_effects.
//  - num_categories: Number of categories per variable.
//  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
//
// Returns:
//  - A flat vector of parameters containing overall effects and (if active) group differences.
//
// Notes:
//  - The order of pairs in `pairwise_effects` must match the upper-triangle order
//    of (var1,var2) pairs as constructed in R.
//  - The length of the output vector depends on both the number of groups
//    and the active entries in inclusion_indicator.
//  - This function is the inverse of `unvectorize_model_parameters_bgmcompare()`.
arma::vec vectorize_model_parameters_bgmcompare(
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& inclusion_indicator,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable
) {
  const int num_variables = inclusion_indicator.n_rows;
  const int num_groups = main_effects.n_cols;
  const int n_main_rows = count_num_main_effects(
    num_categories, is_ordinal_variable
  );
  const int n_pair_rows = num_variables * (num_variables - 1) / 2;

  // total length
  int total_len = 0;
  total_len += n_main_rows;      // MAIN overall (col 0)
  total_len += n_pair_rows;      // PAIR overall (col 0)

  // MAIN differences
  for (int v = 0; v < num_variables; ++v) {
    if (inclusion_indicator(v, v) == 1) {
      const int r0 = main_effect_indices(v, 0);
      const int r1 = main_effect_indices(v, 1);
      total_len += static_cast<int>(r1 - r0 + 1) * (num_groups - 1);
    }
  }
  // PAIRWISE differences
  for (int v1 = 0; v1 < num_variables - 1; ++v1) {
    for (int v2 = v1 + 1; v2 < num_variables; ++v2) {
      if (inclusion_indicator(v1, v2) == 1) total_len += (num_groups - 1);
    }
  }

  arma::vec param_vec(total_len, arma::fill::zeros);
  int off = 0;

  // 1) MAIN overall (col 0) — vectorized
  param_vec.subvec(off, off + n_main_rows - 1) = main_effects.col(0);
  off += n_main_rows;

  // 2) PAIRWISE overall (col 0) — vectorized
  // (Relies on rows being in the same upper-triangle order as constructed in R.)
  param_vec.subvec(off, off + n_pair_rows - 1) = pairwise_effects.col(0);
  off += n_pair_rows;

  // 3) MAIN differences (cols 1..G-1) for selected variables
  for (int v = 0; v < num_variables; ++v) {
    if (inclusion_indicator(v, v) != 1) continue;
    const int r0 = main_effect_indices(v, 0);
    const int r1 = main_effect_indices(v, 1);
    for (int r = r0; r <= r1; ++r) {
      for (int g = 1; g < num_groups; ++g) {
        param_vec(off++) = main_effects(r, g);
      }
    }
  }

  // 4) PAIRWISE differences (cols 1..G-1) for selected pairs
  for (int v1 = 0; v1 < num_variables - 1; ++v1) {
    for (int v2 = v1 + 1; v2 < num_variables; ++v2) {
      if (inclusion_indicator(v1, v2) != 1) continue;
      const int row = pairwise_effect_indices(v1, v2);
      for (int g = 1; g < num_groups; ++g) {
        param_vec(off++) = pairwise_effects(row, g);
      }
    }
  }

  return param_vec;
}



// Reconstructs main-effect and pairwise-effect matrices from a flat parameter vector (bgmCompare model).
//
// The input vector must follow the layout produced by `vectorize_model_parameters_bgmcompare()`:
//  1. Main-effect overall parameters (column 0 of main_effects), stacked by variable.
//  2. Pairwise-effect overall parameters (column 0 of pairwise_effects), stacked by pair.
//  3. Main-effect group differences (columns 1..G-1), included only if
//     the variable is active in inclusion_indicator(v,v).
//  4. Pairwise-effect group differences (columns 1..G-1), included only if
//     the pair is active in inclusion_indicator(v1,v2).
//
// Inputs:
//  - param_vec: Flattened parameter vector.
//  - inclusion_indicator: Symmetric binary matrix of active variables (diag) and pairs (off-diag).
//  - main_effect_indices: Index ranges [row_start, row_end] for each variable in main_effects.
//  - pairwise_effect_indices: Lookup table mapping (var1,var2) → row in pairwise_effects.
//  - num_groups: Number of groups (columns in main_effects / pairwise_effects).
//  - num_categories: Number of categories per variable.
//  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
//
// Outputs:
//  - main_effects_out: Matrix of main-effect parameters (rows = categories, cols = groups).
//  - pairwise_effects_out: Matrix of pairwise-effect parameters (rows = pairs, cols = groups).
//
// Notes:
//  - The vector must have exactly the length returned by `vectorize_model_parameters_bgmcompare()`.
//  - Diagonal entries in `inclusion_indicator` determine whether main-effect
//    group differences are included.
//  - Off-diagonal entries in `inclusion_indicator` determine whether
//    pairwise-effect group differences are included.
//  - This function is the inverse of `vectorize_model_parameters_bgmcompare()`.
void unvectorize_model_parameters_bgmcompare(
    const arma::vec& param_vec,
    arma::mat& main_effects_out,                 // [n_main_rows × G]
    arma::mat& pairwise_effects_out,             // [n_pair_rows × G]
    const arma::imat& inclusion_indicator,       // [V × V]
    const arma::imat& main_effect_indices,       // [V × 2], inclusive [start,end]
    const arma::imat& pairwise_effect_indices,   // [V × V]
    const int num_groups,                          // G
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable
) {
  const int num_variables = inclusion_indicator.n_rows;
  const int num_main = count_num_main_effects(
    num_categories, is_ordinal_variable
  );
  const int num_pair = num_variables * (num_variables - 1) / 2;

  // Only reallocate if sizes differ (optimization to avoid repeated allocation)
  if (main_effects_out.n_rows != static_cast<arma::uword>(num_main) ||
      main_effects_out.n_cols != static_cast<arma::uword>(num_groups)) {
    main_effects_out.set_size(num_main, num_groups);
  }
  if (pairwise_effects_out.n_rows != static_cast<arma::uword>(num_pair) ||
      pairwise_effects_out.n_cols != static_cast<arma::uword>(num_groups)) {
    pairwise_effects_out.set_size(num_pair, num_groups);
  }
  main_effects_out.zeros();
  pairwise_effects_out.zeros();

  int off = 0;

  // 1) MAIN overall (col 0) — vectorized
  main_effects_out.col(0) = param_vec.subvec(off, off + num_main - 1);
  off += num_main;

  // 2) PAIRWISE overall (col 0) — vectorized
  pairwise_effects_out.col(0) = param_vec.subvec(off, off + num_pair - 1);
  off += num_pair;

  // 3) MAIN differences (cols 1..G-1) for selected variables
  for (int v = 0; v < num_variables; v++) {
    if (inclusion_indicator(v, v) == 0) continue;
    const int r0 = main_effect_indices(v, 0);
    const int r1 = main_effect_indices(v, 1);
    for (int r = r0; r <= r1; ++r) {
      for (int g = 1; g < num_groups; ++g) {
        main_effects_out(r, g) = param_vec(off++);
      }
    }
  }

  // 4) PAIRWISE differences (cols 1..G-1) for selected pairs
  for (int v1 = 0; v1 < num_variables - 1; ++v1) {
    for (int v2 = v1 + 1; v2 < num_variables; ++v2) {
      if (inclusion_indicator(v1, v2) == 0) continue;
      const int row = pairwise_effect_indices(v1, v2);
      for (int g = 1; g < num_groups; ++g) {
        pairwise_effects_out(row, g) = param_vec(off++);
      }
    }
  }
}



// Builds index maps linking matrix entries to positions in the vectorized parameter vector (bgmCompare model).
//
// The index maps are used to quickly locate where each parameter (main-effect or pairwise-effect,
// across groups) sits inside the flattened parameter vector produced by
// `vectorize_model_parameters_bgmcompare()`.
//
// Layout:
//  1. Main-effect overall parameters (col 0).
//  2. Pairwise-effect overall parameters (col 0).
//  3. Main-effect group differences (cols 1..G-1), included only if
//     the variable is active in inclusion_indicator(v,v).
//  4. Pairwise-effect group differences (cols 1..G-1), included only if
//     the pair is active in inclusion_indicator(v1,v2).
//
// Inputs:
//  - main_effects: Matrix of main-effect parameters (used for dimension info).
//  - pairwise_effects: Matrix of pairwise-effect parameters (used for dimension info).
//  - inclusion_indicator: Symmetric binary matrix of active variables (diag) and pairs (off-diag).
//  - main_effect_indices: Index ranges [row_start, row_end] for each variable in main_effects.
//  - pairwise_effect_indices: Lookup table mapping (var1,var2) → row in pairwise_effects.
//  - num_categories: Number of categories per variable.
//  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
//
// Returns:
//  - A pair of integer matrices:
//      * main_index: [num_main × num_groups], with each entry giving the position
//        in the parameter vector for that main-effect parameter (or -1 if inactive).
//      * pair_index: [num_pair × num_groups], with each entry giving the position
//        in the parameter vector for that pairwise-effect parameter (or -1 if inactive).
//
// Notes:
//  - Entries are set to -1 when the corresponding parameter is inactive.
//  - The returned index maps must always be consistent with the ordering used
//    in vectorization/unvectorization.
//  - A final check (e.g. verifying that `off == param_vec.n_elem`) can help
//    catch mismatches between index maps and vectorizer logic.
std::pair<arma::imat, arma::imat> build_index_maps(
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& inclusion_indicator,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable
) {
  const int num_variables = inclusion_indicator.n_rows;
  const int num_groups    = main_effects.n_cols;
  const int num_main = count_num_main_effects(
    num_categories, is_ordinal_variable
  );
  const int num_pair = num_variables * (num_variables - 1) / 2;

  arma::imat main_index(num_main, num_groups, arma::fill::value(-1));
  arma::imat pair_index(num_pair, num_groups, arma::fill::value(-1));

  int off = 0;

  // 1) main overall (col 0)
  for (int r = 0; r < num_main; ++r) {
    main_index(r,0) = off++;
  }

  // 2) pair overall (col 0)
  for (int r = 0; r < num_pair; ++r) {
    pair_index(r,0) = off++;
  }

  // 3) main differences
  for (int v = 0; v < num_variables; ++v) {
    if (inclusion_indicator(v,v) == 0) continue;
    const int r0 = main_effect_indices(v,0);
    const int r1 = main_effect_indices(v,1);
    for (int r = r0; r <= r1; ++r) {
      for (int g = 1; g < num_groups; ++g) {
        main_index(r,g) = off++;
      }
    }
  }

  // 4) pairwise differences
  for (int v1 = 0; v1 < num_variables-1; ++v1) {
    for (int v2 = v1+1; v2 < num_variables; ++v2) {
      if (inclusion_indicator(v1,v2) == 0) continue;
      const int row = pairwise_effect_indices(v1,v2);
      for (int g = 1; g < num_groups; ++g) {
        pair_index(row,g) = off++;
      }
    }
  }

  return {main_index, pair_index};
}



// Extracts entries of the inverse mass matrix corresponding to active parameters (bgmCompare model).
//
// If `selection` is false, the full diagonal vector is returned unchanged.
// If `selection` is true, the output is restricted to:
//  1. Main-effect overall parameters (column 0).
//  2. Pairwise-effect overall parameters (column 0).
//  3. Main-effect group differences (columns 1..G-1) for variables with
//     inclusion_indicator(v,v) == 1.
//  4. Pairwise-effect group differences (columns 1..G-1) for pairs with
//     inclusion_indicator(v1,v2) == 1.
//
// Inputs:
//  - inv_diag: Full inverse mass diagonal (length = all parameters).
//  - inclusion_indicator: Symmetric binary matrix of active variables (diag) and pairs (off-diag).
//  - num_groups: Number of groups.
//  - num_categories: Number of categories per variable.
//  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
//  - main_index: Index map for main effects (from build_index_maps()).
//  - pair_index: Index map for pairwise effects (from build_index_maps()).
//  - main_effect_indices: Index ranges [row_start, row_end] for each variable in main_effects.
//  - pairwise_effect_indices: Lookup table mapping (var1,var2) → row in pairwise_effects.
//  - selection: If true, restrict to active parameters; if false, return full inv_diag.
//
// Returns:
//  - A vector containing inverse mass entries for active parameters only.
//
// Notes:
//  - Must be consistent with the layout in `vectorize_model_parameters_bgmcompare()`.
//  - Index maps (`main_index`, `pair_index`) are required to locate group-difference entries.
arma::vec inv_mass_active(
    const arma::vec& inv_diag,
    const arma::imat& inclusion_indicator,
    const int num_groups,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable,
    const arma::imat& main_index,
    const arma::imat& pair_index,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const bool& selection
) {
  if(selection == false)
    return inv_diag;

  const int num_variables = inclusion_indicator.n_rows;
  const int num_main = count_num_main_effects(num_categories, is_ordinal_variable);
  const int num_pair = num_variables * (num_variables - 1) / 2;

  // total length
  int total_len = 0;
  total_len += num_main;      // MAIN overall (col 0)
  total_len += num_pair;      // PAIR overall (col 0)

  // MAIN differences
  for (int v = 0; v < num_variables; ++v) {
    if (inclusion_indicator(v, v) == 0) continue;
    const int r0 = main_effect_indices(v, 0);
    const int r1 = main_effect_indices(v, 1);
    total_len += static_cast<int>(r1 - r0 + 1) * (num_groups - 1);
  }

  // PAIRWISE differences
  for (int v1 = 0; v1 < num_variables - 1; ++v1) {
    for (int v2 = v1 + 1; v2 < num_variables; ++v2) {
      if (inclusion_indicator(v1, v2) == 1) total_len += (num_groups - 1);
    }
  }

  arma::vec active_inv_diag(total_len, arma::fill::zeros);

  int off = 0;

  // 1) MAIN overall (col 0) — vectorized
  active_inv_diag.subvec(off, off + num_main - 1) = inv_diag.subvec(off, off + num_main - 1);
  off += num_main;

  // 2) PAIRWISE overall (col 0) — vectorized
  // (Relies on rows being in the same upper-triangle order as constructed in R.)
  active_inv_diag.subvec(off, off + num_pair - 1) = inv_diag.subvec(off, off + num_pair - 1);
  off += num_pair;

  // inv_diag was learned in stage 2 with all indicators on, i.e. in the full
  // parameter layout. Read it through full-layout positions (main_index /
  // pair_index reflect the selection-reduced layout and would misalign once
  // parameters are excluded). The full layout is: main overall, pair overall,
  // then per-row main differences and per-row pairwise differences, each row
  // repeated over the G-1 group contrasts.
  const int diff_base = num_main + num_pair;
  const int pair_diff_base = diff_base + num_main * (num_groups - 1);

  // 3) MAIN differences (cols 1..G-1) for selected variables
  for (int v = 0; v < num_variables; ++v) {
    if (inclusion_indicator(v, v) == 0) continue;
    const int r0 = main_effect_indices(v, 0);
    const int r1 = main_effect_indices(v, 1);
    for (int r = r0; r <= r1; ++r) {
      for (int g = 1; g < num_groups; ++g) {
        int idx = diff_base + r * (num_groups - 1) + (g - 1);
        active_inv_diag(off++) = inv_diag(idx);
      }
    }
  }

  // 4) PAIRWISE differences (cols 1..G-1) for selected pairs
  for (int v1 = 0; v1 < num_variables - 1; ++v1) {
    for (int v2 = v1 + 1; v2 < num_variables; ++v2) {
      if (inclusion_indicator(v1, v2) != 1) continue;
      const int row = pairwise_effect_indices(v1, v2);
      for (int g = 1; g < num_groups; ++g) {
        int idx = pair_diff_base + row * (num_groups - 1) + (g - 1);
        active_inv_diag(off++) = inv_diag(idx);
      }
    }
  }

  return active_inv_diag;
}
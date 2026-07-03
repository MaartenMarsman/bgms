# ==============================================================================
# Compute utilities
# ==============================================================================
#
# Pure helper functions for computing derived quantities from data
# (category counts, sufficient statistics). Each function is pure:
# input -> output (no side-effects).
# ==============================================================================


# ------------------------------------------------------------------------------
# compute_counts_per_category
# ------------------------------------------------------------------------------
#
# Compute per-group category counts for each variable. Used to build
# precomputed structures for the bgmCompare C++ backend.
#
# @param x  Numeric matrix: the recoded data.
# @param num_categories  Integer vector: max category per variable.
# @param group  Integer vector: group membership.
#
# Returns: list of matrices (one per group), each max_cat x num_variables.
# ------------------------------------------------------------------------------
compute_counts_per_category = function(x, num_categories, group = NULL) {
  counts_per_category = list()
  for(g in unique(group)) {
    counts_per_category_gr = matrix(0, nrow = max(num_categories), ncol = ncol(x))
    for(variable in seq_len(ncol(x))) {
      for(category in seq_len(num_categories[variable])) {
        counts_per_category_gr[category, variable] = sum(x[group == g, variable] == category)
      }
    }
    counts_per_category[[length(counts_per_category) + 1]] = counts_per_category_gr
  }
  return(counts_per_category)
}


# ------------------------------------------------------------------------------
# compute_blume_capel_stats
# ------------------------------------------------------------------------------
#
# Compute sufficient statistics for Blume-Capel variables (linear and
# quadratic deviations from baseline). Used to build precomputed
# structures for the bgmCompare C++ backend.
#
# @param x  Numeric matrix: the recoded data.
# @param baseline_category  Integer vector: baseline categories.
# @param ordinal_variable  Logical vector: TRUE = ordinal, FALSE = BC.
# @param group  Integer vector or NULL: group membership.
#
# Returns: matrix (one-group) or list of matrices (multi-group),
#   each 2 x num_variables (row 1 = linear, row 2 = quadratic).
# ------------------------------------------------------------------------------
compute_blume_capel_stats = function(x, baseline_category, ordinal_variable, group = NULL) {
  if(is.null(group)) { # One-group design
    sufficient_stats = matrix(0, nrow = 2, ncol = ncol(x))
    bc_vars = which(!ordinal_variable)
    for(i in bc_vars) {
      sufficient_stats[1, i] = sum(x[, i] - baseline_category[i])
      sufficient_stats[2, i] = sum((x[, i] - baseline_category[i])^2)
    }
    return(sufficient_stats)
  } else { # Multi-group design
    sufficient_stats = list()
    for(g in unique(group)) {
      sufficient_stats_gr = matrix(0, nrow = 2, ncol = ncol(x))
      bc_vars = which(!ordinal_variable)
      for(i in bc_vars) {
        sufficient_stats_gr[1, i] = sum(x[group == g, i] - baseline_category[i])
        sufficient_stats_gr[2, i] = sum((x[group == g, i] - baseline_category[i])^2)
      }
      sufficient_stats[[length(sufficient_stats) + 1]] = sufficient_stats_gr
    }
    return(sufficient_stats)
  }
}


# ------------------------------------------------------------------------------
# compute_pairwise_stats
# ------------------------------------------------------------------------------
#
# Compute sufficient statistics for pairwise interactions (cross-product
# of observations per group). Used to build precomputed structures for
# the bgmCompare C++ backend.
#
# @param x  Numeric matrix: the centered data.
# @param group  Integer vector: group membership.
#
# Returns: list of p x p cross-product matrices (one per group).
# ------------------------------------------------------------------------------
compute_pairwise_stats = function(x, group) {
  result = list()

  for(g in unique(group)) {
    obs = x[group == g, , drop = FALSE]
    # cross-product: gives number of co-occurrences of categories
    result[[length(result) + 1]] = t(obs) %*% obs
  }

  result
}

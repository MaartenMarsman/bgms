# ==============================================================================
# Data validation functions
# ==============================================================================
#
# Pure validation and reformatting functions for observational data.
# Each function takes input and returns validated output (or errors).
# ==============================================================================


# ------------------------------------------------------------------------------
# data_check
# ------------------------------------------------------------------------------
#
# Coerce user data to a numeric matrix and perform basic dimension checks.
# Used by bgm_spec() as the first validation step for all model types.
#
# @param data  Data frame or matrix: the raw user input.
# @param name  Character: name of the argument (for error messages).
#
# Returns: numeric matrix.
# ------------------------------------------------------------------------------
data_check = function(data, name) {
  if(!inherits(data, c("matrix", "data.frame"))) {
    stop(paste(name, "must be a matrix or data.frame."))
  }
  if(inherits(data, "data.frame")) {
    data = data.matrix(data)
  }
  if(nrow(data) < 2 || ncol(data) < 2) {
    stop(paste(name, "must have at least 2 rows and 2 columns."))
  }
  return(data)
}


# ------------------------------------------------------------------------------
# center_continuous_data
# ------------------------------------------------------------------------------
#
# Column-centers a numeric data matrix. This is required for GGM models
# because the likelihood is formulated in terms of the precision matrix
# of a zero-mean Gaussian:
#
#   log p(X | Omega) \propto (n/2) log|Omega| - (1/2) tr(Omega S)
#
# where S = X'X on centered data. Without centering, S conflates the
# mean structure with the precision matrix, biasing estimates when
# column means are non-zero.
#
# @param x  Numeric matrix: the data (after missing-data handling).
#
# Returns:
#   Column-centered matrix (same dimensions, same colnames).
# ------------------------------------------------------------------------------
center_continuous_data = function(x) {
  sweep(x, 2, colMeans(x))
}


# ------------------------------------------------------------------------------
# validate_missing_data
# ------------------------------------------------------------------------------
#
# Handles missing data for all model types (OMRF, GGM, bgmCompare).
# Either removes rows with missing values (listwise) or identifies
# missing entries for C++ imputation.
#
# @param x  Numeric matrix: the data.
# @param na_action  Character: "listwise" or "impute".
# @param is_continuous  Logical: TRUE for GGM (continuous) models.
#   Imputation is supported for all model types, including GGM. This
#   argument is currently unused by the function body.
# @param group  Optional integer vector: group indicators for bgmCompare.
#   If provided, listwise deletion also filters the group vector.
#   NULL for bgm() calls.
#
# Returns:
#   list(x, na_impute, missing_index, n_removed, [group])
#   - x: data matrix (rows removed for listwise, NAs imputed for impute)
#   - na_impute: logical
#   - missing_index: matrix(NA, 1, 1) if no imputation,
#       otherwise Nx2 matrix of 0-based (row, col) indices
#   - n_removed: integer count of rows removed (listwise only, 0 otherwise)
#   - group: filtered group vector (only present when group was non-NULL)
# ------------------------------------------------------------------------------
validate_missing_data = function(x,
                                 na_action,
                                 is_continuous = FALSE,
                                 group = NULL) {
  if(na_action == "listwise") {
    return(handle_listwise(x, group))
  }

  # --- impute path ---
  handle_impute(x, group)
}


# ------------------------------------------------------------------------------
# handle_listwise (internal helper)
# ------------------------------------------------------------------------------
handle_listwise = function(x, group = NULL) {
  missing_rows = apply(x, 1, anyNA)

  if(all(missing_rows)) {
    stop(paste0(
      "All rows in x contain at least one missing response.\n",
      "You could try option na_action = impute."
    ))
  }

  n_removed = sum(missing_rows)
  if(n_removed > 0 && isTRUE(getOption("bgms.verbose", TRUE))) {
    n_remaining = nrow(x) - n_removed
    message(
      n_removed, " row", if(n_removed > 1) "s" else "",
      " with missing values excluded (n = ", n_remaining, " remaining).\n",
      "To impute missing values instead, use na_action = \"impute\"."
    )
  }

  x = x[!missing_rows, , drop = FALSE]

  if(is.null(ncol(x)) || ncol(x) < 2) {
    stop(paste0(
      "After removing missing observations from the input matrix x,\n",
      "there were less than two columns left in x."
    ))
  }
  if(is.null(nrow(x)) || nrow(x) < 2) {
    stop(paste0(
      "After removing missing observations from the input matrix x,\n",
      "there were less than two rows left in x."
    ))
  }

  result = list(
    x             = x,
    na_impute     = FALSE,
    missing_index = matrix(NA, nrow = 1, ncol = 1),
    n_removed     = n_removed
  )

  if(!is.null(group)) {
    result$group = group[!missing_rows]
  }

  result
}


# ------------------------------------------------------------------------------
# handle_impute (internal helper)
# ------------------------------------------------------------------------------
handle_impute = function(x, group = NULL) {
  num_missings = sum(is.na(x))

  if(num_missings == 0) {
    result = list(
      x             = x,
      na_impute     = FALSE,
      missing_index = matrix(NA, nrow = 1, ncol = 1),
      n_removed     = 0L
    )
    if(!is.null(group)) result$group = group
    return(result)
  }

  # Guard: entire-column-missing
  for(v in seq_len(ncol(x))) {
    if(all(is.na(x[, v]))) {
      stop(
        "Variable '", colnames(x)[v], "' has no observed values. ",
        "Remove it before fitting."
      )
    }
  }

  num_variables = ncol(x)
  missing_index = matrix(0, nrow = num_missings, ncol = 2)
  cntr = 0
  for(node in seq_len(num_variables)) {
    mis = which(is.na(x[, node]))
    if(length(mis) > 0) {
      observed = x[-mis, node]
      for(i in seq_along(mis)) {
        cntr = cntr + 1
        missing_index[cntr, 1] = mis[i] - 1 # C++ 0-based index
        missing_index[cntr, 2] = node - 1 # C++ 0-based index
        # Index explicitly: sample(v, 1) treats a length-1 numeric v as 1:v.
        x[mis[i], node] = observed[sample.int(length(observed), 1)]
      }
    }
  }

  result = list(
    x             = x,
    na_impute     = TRUE,
    missing_index = missing_index,
    n_removed     = 0L
  )
  if(!is.null(group)) result$group = group

  result
}


# ------------------------------------------------------------------------------
# reformat_ordinal_data
# ------------------------------------------------------------------------------
#
# Per-variable recoding of ordinal and Blume-Capel data to 0-based
# contiguous categories. Single-group only (no group-conditional
# collapsing --- see collapse_categories_across_groups() for that).
#
# @param x  Numeric matrix: the data (after missing-data handling).
# @param is_ordinal  Logical vector of length ncol(x): TRUE = regular
#   ordinal variable, FALSE = Blume-Capel variable.
# @param baseline_category  Integer vector of length ncol(x): baseline
#   (reference) categories for Blume-Capel variables.
#
# Returns:
#   list(x, num_categories, baseline_category)
#   - x: matrix with recoded values
#   - num_categories: integer vector (max observed value per variable)
#   - baseline_category: possibly adjusted baseline categories
# ------------------------------------------------------------------------------
reformat_ordinal_data = function(x, is_ordinal, baseline_category) {
  check_fail_zero = FALSE
  num_variables = ncol(x)
  num_categories = vector(length = num_variables)
  # Per ordinal variable: the sorted unique ORIGINAL values, i.e. the recode
  # map. The recoded category of an original value v is match(v, levels) - 1.
  # predict() needs this to recode newdata the same way (NULL for Blume-Capel).
  category_levels = vector("list", num_variables)
  # Per Blume-Capel variable: the additive shift applied to reach the 0-based
  # internal scale (NA for regular ordinal and continuous variables). predict()
  # subtracts it from newdata and simulate() adds it back.
  blume_capel_shift = rep(NA_real_, num_variables)

  for(node in 1:num_variables) {
    unq_vls = sort(unique(x[, node]))
    mx_vl = max(unq_vls)

    # Recode data --------------------------------------------------------------
    if(is_ordinal[node]) { # Regular ordinal variable
      # Capture the recode map (original sorted values) before recoding x.
      category_levels[[node]] = unq_vls
      # A regular ordinal variable needs repeated values: its category
      # thresholds are unidentified if every response is distinct (the column
      # then looks continuous). Blume-Capel is parametric in the category score
      # and is exempt.
      if(length(unq_vls) == nrow(x)) {
        stop(paste0(
          "Only unique responses observed for variable ",
          node,
          ". We expect >= 1 observations per category."
        ))
      }
      if(length(unq_vls) != mx_vl + 1 || any(unq_vls != 0:mx_vl)) {
        y = x[, node]
        cntr = 0
        for(value in unq_vls) {
          x[y == value, node] = cntr
          cntr = cntr + 1
        }
      }
    } else { # Blume-Capel ordinal variable
      # Check if observations are integer or can be recoded --------------------
      if(any(abs(unq_vls - round(unq_vls)) > .Machine$double.eps)) {
        int_unq_vls = unique(as.integer(unq_vls))
        if(anyNA(int_unq_vls)) {
          stop(paste0(
            "The Blume-Capel model assumes that its observations are coded as integers, but \n",
            "the category scores for node ", node, " were not integer. An attempt to recode \n",
            "them to integer failed. Please inspect the documentation for the base R \n",
            "function as.integer(), which bgm uses for recoding category scores."
          ))
        }

        if(length(int_unq_vls) != length(unq_vls)) {
          stop(paste0(
            "The Blume-Capel model assumes that its observations are coded as integers. The \n",
            "category scores of the observations for node ", node, " were not integers. An \n",
            "attempt to recode these observations as integers failed because, after rounding, \n",
            "a single integer value was used for several observed score categories."
          ))
        }
        x[, node] = as.integer(x[, node])

        if(baseline_category[node] < 0 || baseline_category[node] > max(x[, node])) {
          stop(paste0(
            "The reference category for the Blume-Capel variable ", node, "is outside its \n",
            "range of observations."
          ))
        }
      }

      # Check if observations start at zero and recode otherwise ---------------
      blume_capel_shift[node] = min(x[, node])
      if(min(x[, node]) != 0) {
        baseline_category[node] = baseline_category[node] - min(x[, node])
        x[, node] = x[, node] - min(x[, node])

        if(check_fail_zero == FALSE) {
          check_fail_zero = TRUE
          failed_zeroes = c(node)
        } else {
          failed_zeroes = c(failed_zeroes, node)
        }
      }

      check_range = length(unique(x[, node]))
      if(check_range < 3) {
        stop(paste0(
          "The Blume-Capel is only available for variables with more than one category \n",
          "observed. There two or less categories observed for variable ",
          node,
          "."
        ))
      }
    }

    # Warn that maximum category value is large --------------------------------
    num_categories[node] = max(x[, node])
    if(!is_ordinal[node] && num_categories[node] > 10) {
      warning(
        "Blume-Capel variable ", node, " has ", num_categories[node], " categories. ",
        "This may slow computation. Empty categories are not collapsed.",
        call. = FALSE
      )
    }

    # Check to see if not all responses are in one category --------------------
    if(num_categories[node] == 0) {
      stop(paste0(
        "Only one value [",
        unq_vls,
        "] was observed for variable ",
        node,
        "."
      ))
    }
  }

  if(check_fail_zero == TRUE && isTRUE(getOption("bgms.verbose", TRUE))) {
    nodes_str = paste(failed_zeroes, collapse = ", ")
    message(
      "Variable", if(length(failed_zeroes) > 1) "s" else "", " ", nodes_str,
      " recoded to start at 0 (baseline categor",
      if(length(failed_zeroes) > 1) "ies" else "y", " adjusted)."
    )
  }

  list(
    x                 = x,
    num_categories    = num_categories,
    baseline_category = baseline_category,
    category_levels   = category_levels,
    blume_capel_shift = blume_capel_shift
  )
}


# ------------------------------------------------------------------------------
# collapse_categories_across_groups
# ------------------------------------------------------------------------------
#
# For bgmCompare: recodes each regular ordinal variable onto the contiguous
# *union* of the category values observed in any group. A value that no group
# observes -- a true gap in the scale -- carries no information anywhere, so
# it is dropped and the remaining categories are renumbered contiguously
# (0-based). A value that at least one group observes is always retained,
# even when some other group never observes it.
#
# Retaining such a category makes it a structural zero for the groups that do
# not observe it: those groups contribute no observations to its threshold, so
# the corresponding group difference is driven by the prior rather than by
# data. That is the honest representation --- the alternative, merging the
# category away, silently redefines the variable in the groups that *do*
# observe it. The function warns whenever this happens.
#
# Blume-Capel variables are exempt. Their thresholds are parametric functions
# of the category *score*, so renumbering categories would change the model
# rather than relabel it, and an unobserved score is still a meaningful point
# on the scale. See the note in man/bgmCompare.Rd.
#
# Called immediately after reformat_ordinal_data() in the compare path.
# reformat_ordinal_data() already maps the pooled data onto contiguous codes,
# so on that input the ordinal branch here is a relabel-free pass that
# computes the per-group support; it is written to be correct standalone.
#
# @param x  Numeric matrix: data already recoded by reformat_ordinal_data().
# @param group  Integer vector of length nrow(x): group membership (1:K).
# @param is_ordinal  Logical vector of length ncol(x).
# @param num_categories  Integer vector from reformat_ordinal_data().
# @param baseline_category  Integer vector from reformat_ordinal_data().
#
# Returns:
#   list(x, num_categories, baseline_category, category_support)
#   - category_support: list, one entry per variable (NULL for Blume-Capel),
#     each a (num_categories + 1) x num_groups integer matrix of observation
#     counts on the final 0-based category codes.
# ------------------------------------------------------------------------------
collapse_categories_across_groups = function(x,
                                             group,
                                             is_ordinal,
                                             num_categories,
                                             baseline_category) {
  num_variables = ncol(x)
  group_ids = sort(unique(group))
  num_groups = length(group_ids)

  variable_names = colnames(x)
  variable_label = function(node) {
    if(is.null(variable_names) || is.na(variable_names[node]) ||
      !nzchar(variable_names[node])) {
      paste0("variable ", node)
    } else {
      paste0("variable '", variable_names[node], "'")
    }
  }
  group_label = function(g) {
    paste0("group ", group_ids[g])
  }

  category_support = vector("list", num_variables)
  renumbered = integer(0) # variables where a true gap was closed
  gaps_dropped = integer(0) # how many gap values each of those lost
  ref_cells = character(0) # empty cells in the reference category
  zero_cells = character(0) # "<variable>, category <c>, <group>"

  for(node in seq_len(num_variables)) {
    if(!is_ordinal[node]) next # Blume-Capel variables: exempt, see above

    unq_vls = sort(unique(x[, node]))
    n_unique = length(unq_vls)

    # A variable with a single observed value carries no threshold information.
    if(n_unique < 2) {
      stop(paste0("Only one value was observed for variable ", node, "."))
    }

    # Per-group counts on the retained (union) categories.
    support = matrix(0L, nrow = n_unique, ncol = num_groups)
    for(g in seq_len(num_groups)) {
      x_g = x[group == group_ids[g], node]
      support[, g] = tabulate(match(x_g, unq_vls), nbins = n_unique)
    }

    # Recode onto the contiguous union. Only true gaps -- values no group
    # observes -- move a category's code; with no gaps this is the identity.
    target = seq_len(n_unique) - 1L
    if(!isTRUE(all.equal(as.numeric(unq_vls), as.numeric(target)))) {
      renumbered = c(renumbered, node)
      gaps_dropped = c(gaps_dropped, as.integer(max(unq_vls) + 1 - n_unique))
      original = x[, node]
      for(i in seq_along(unq_vls)) {
        x[original == unq_vls[i], node] = target[i]
      }
    }

    num_categories[node] = n_unique - 1L

    dimnames(support) = list(
      paste0("category ", target),
      paste0("group ", group_ids)
    )
    category_support[[node]] = support

    # Structural zeros: a retained category with no observations in some group.
    # An empty cell in the REFERENCE category is the consequential one. Every
    # threshold is identified relative to category 0, so a group that never
    # used it has nothing fixing the level of its whole threshold vector, not
    # just one threshold. Those cells say so, and are listed first so the
    # head() below cannot drop them behind ordinary ones.
    empty = which(support == 0L, arr.ind = TRUE)
    if(nrow(empty) > 0) {
      is_reference = empty[, 1] == 1L
      cells = paste0(
        variable_label(node), ", category ", target[empty[, 1]], ", ",
        vapply(empty[, 2], group_label, character(1))
      )
      cells[is_reference] = paste0(
        cells[is_reference],
        " -- the reference category; every threshold of this variable is",
        " affected for that group"
      )
      ref_cells = c(ref_cells, cells[is_reference])
      zero_cells = c(zero_cells, cells[!is_reference])
    }
  }

  # --- Reporting ---------------------------------------------------------------
  # Two conditions, two volumes: renumbering a gap is benign bookkeeping, an
  # empty cell in a retained category changes what the group difference means.
  if(length(renumbered) > 0 && isTRUE(getOption("bgms.verbose", TRUE))) {
    message(
      "Some category values were not used by any group. They were dropped ",
      "and the remaining categories renumbered, for ",
      paste0(
        vapply(renumbered, variable_label, character(1)),
        " (", gaps_dropped, " dropped)",
        collapse = ", "
      ),
      ". No observed category was merged."
    )
  }

  affected_cells = c(ref_cells, zero_cells)
  if(length(affected_cells) > 0) {
    shown = utils::head(affected_cells, 10L)
    extra = length(affected_cells) - length(shown)
    # Classed, so a caller can catch this one condition without muffling every
    # warning the fit might raise. The class is part of the user-facing API.
    #
    # Kept tight on purpose: R truncates a warning at getOption("warning.length")
    # = 1000 characters by default, and the tail of this one is the part that
    # tells the reader where to look next.
    warning(warningCondition(
      paste0(
        "Some categories were not used by every group:\n",
        paste0("  ", shown, collapse = "\n"),
        if(extra > 0) paste0("\n  ... and ", extra, " more") else "",
        "\nThese categories are kept, because the other groups do use them. ",
        "But a group with no observations in a category has nothing to say ",
        "about where its threshold for that category lies, so the reported ",
        "difference for that group and that category is set by the prior, not ",
        "by the data. An empty reference category is worse: every threshold is ",
        "measured relative to category 0, so all of that variable's threshold ",
        "differences for that group rest on the prior, not just one. Expect ",
        "large, very uncertain numbers, and do not read them as evidence of a ",
        "group difference. Only the category thresholds are affected, not the ",
        "pairwise (edge) differences. The printed summary marks the rows; see ",
        "?summary.bgmCompare."
      ),
      class = "bgms_group_support_warning"
    ))
  }

  list(
    x                 = x,
    num_categories    = num_categories,
    baseline_category = baseline_category,
    category_support  = category_support
  )
}

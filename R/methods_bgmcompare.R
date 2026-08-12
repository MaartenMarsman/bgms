# R/methods_bgmcompare.R

# ---- group labels on human-facing displays ----------------------------------
# bgmCompare numbers its groups 1..K by first appearance in the group indicator
# and every extractor keys on that number: the `group1` / `group2` column names
# are a downstream contract (easybgm, JASP) and must stay numeric. The original
# indicator values are stored alongside as `group_labels`, and the displays a
# person reads -- print/summary headers, plot titles, calibration panels,
# centrality labels -- name the group as well as number it. Fits made before
# the field existed carry no labels and degrade to the bare number.

compare_group_labels = function(arguments, num_groups = NULL) {
  labels = arguments$group_labels
  if(is.null(labels) || !length(labels)) {
    return(NULL)
  }
  labels = as.character(labels)
  if(anyNA(labels)) {
    return(NULL)
  }
  if(!is.null(num_groups) && length(labels) != num_groups) {
    return(NULL)
  }
  labels
}

# "group 3" or "group 3 (fr)", for a title or a label. A numeric group
# indicator gives labels that are the numbers themselves, and "group 1 (1)"
# tells a reader nothing they did not already have; the parenthetical is for
# the case where the original value carries information the number does not.
group_tag = function(labels, g) {
  if(is.null(labels) || g > length(labels)) {
    return(sprintf("group %d", g))
  }
  if(identical(labels[g], as.character(g))) {
    return(sprintf("group %d", g))
  }
  sprintf("group %d (%s)", g, labels[g])
}

# The one-line legend that ties the numbers to the labels, with group sizes
# when the fit kept its case-level group vector.
group_mapping_line = function(arguments, num_groups = NULL) {
  labels = compare_group_labels(arguments, num_groups)
  if(is.null(labels)) {
    return(NULL)
  }
  sizes = NULL
  if(!is.null(arguments$group)) {
    counts = tabulate(as.integer(arguments$group), nbins = length(labels))
    if(all(counts > 0L)) sizes = sprintf(" (n = %d)", counts)
  }
  if(is.null(sizes)) sizes = rep("", length(labels))
  paste0(
    "groups: ",
    paste0(seq_along(labels), " = ", labels, sizes, collapse = ", ")
  )
}

# ---- prior-only main-effect difference rows ---------------------------------
# bgmCompare keeps the union of the ordinal categories the groups observe, so a
# retained category can have no observations at all in some group. That group's
# data then say nothing about where its threshold for that category lies, and
# the reported threshold difference is whatever the prior says. The fit warns
# about this once (class `bgms_group_support_warning`) and stores the per-group
# counts in `category_support`; the reader of a saved fit never saw the warning,
# so the printed summary marks the affected rows as well.
#
# Rows are mapped to support cells by index, never by parsing the label. The
# main-effects matrix carries one row per free parameter in variable order --
# `num_categories[v]` thresholds for an ordinal variable, two (linear,
# quadratic) for a Blume-Capel one -- and threshold `k` of variable `v` is the
# threshold for category code `k`, i.e. row `k + 1` of `category_support[[v]]`
# (row 1 is the reference category 0). Blume-Capel variables are exempt from
# the union recode and carry no support matrix, so their rows are never marked.
#
# A row is marked on either of two conditions. Its OWN cell can be empty: that
# group observes no one in that category, so its threshold for the category has
# nothing behind it. Or the REFERENCE cell can be empty: every threshold is
# identified relative to category 0, so a group that never uses category 0 has
# no data fixing the level of its threshold vector at all, and the whole
# vector -- every threshold of that variable, not just one -- rests on the
# prior. The second condition marks all of the variable's rows.
#
# The difference rows are contrasts, not groups: a group's effect is
# `baseline + projection[g, ] %*% differences`, so an empty cell in any group
# reaches every contrast of that variable-by-category pair. A pair is therefore
# marked in all of its contrasts or in none of them.
#
# Returns a logical vector, one entry per row of
# `posterior_summary_main_differences`, or NULL when the fit carries no
# usable `category_support` (fits made before the field existed, and any
# layout the mapping cannot verify). NULL means "print exactly as before".
compare_prior_only_main_diff = function(arguments, num_rows) {
  support = arguments$category_support
  num_variables = arguments$num_variables
  num_groups = arguments$num_groups
  num_categories = arguments$num_categories
  is_ordinal = arguments$is_ordinal_variable

  if(is.null(support) || is.null(num_variables) || is.null(num_groups) ||
    is.null(num_categories) || is.null(is_ordinal)) {
    return(NULL)
  }
  num_variables = as.integer(num_variables)
  num_groups = as.integer(num_groups)
  if(num_groups < 2L) {
    return(NULL)
  }
  if(length(support) != num_variables) {
    return(NULL)
  }
  if(length(num_categories) != num_variables) {
    return(NULL)
  }
  if(length(is_ordinal) != num_variables) {
    return(NULL)
  }

  # one flag per row of the main-effects matrix, in variable order
  by_row = vector("list", num_variables)
  for(v in seq_len(num_variables)) {
    if(!isTRUE(is_ordinal[v])) {
      by_row[[v]] = c(FALSE, FALSE) # Blume-Capel: linear + quadratic
      next
    }
    num_thresholds = as.integer(num_categories[v])
    cells = support[[v]]
    if(!is.matrix(cells) || nrow(cells) != num_thresholds + 1L ||
      ncol(cells) != num_groups) {
      by_row[[v]] = rep(FALSE, num_thresholds)
      next
    }
    # drop the reference row (category 0); row k of the remainder is threshold k
    own = apply(cells[-1L, , drop = FALSE] == 0L, 1L, any)
    # an empty reference cell unfixes the level of the whole threshold vector
    by_row[[v]] = own | any(cells[1L, ] == 0L)
  }
  by_row = unname(unlist(by_row))

  num_contrasts = num_groups - 1L
  # Both summarizers now lay the difference rows out contrast-major: the whole
  # main-effects matrix once per contrast. summarize_main_diff_compare() used to
  # emit them variable-major while labelling them contrast-major, so this
  # branched on difference_selection to follow it.
  flags = rep(by_row, times = num_contrasts)

  if(length(flags) != num_rows) {
    return(NULL)
  }
  flags
}


#' @name print.bgmCompare
#' @title Print method for `bgmCompare` objects
#' @description Minimal console output for `bgmCompare` fit objects.
#' @param x An object of class `bgmCompare`.
#' @param ... Ignored.
#' @return Invisibly returns `x`.
#'
#' @examples
#' \donttest{
#' # See ?bgmCompare for a full example
#' }
#'
#' @seealso [bgmCompare()], [summary.bgmCompare()], [coef.bgmCompare()]
#' @family posterior-methods
#'
#' @export
print.bgmCompare = function(x, ...) {
  arguments = extract_arguments(x)

  # Model type
  if(isTRUE(arguments$difference_selection)) {
    prior_msg = switch(as.character(arguments$difference_prior),
      "Bernoulli" = "Bayesian Difference Selection (Bernoulli prior on inclusion)",
      "Beta-Bernoulli" = "Bayesian Difference Selection (Beta-Bernoulli prior on inclusion)",
      "Bayesian Difference Selection"
    )
    cat(prior_msg, "\n")
  } else {
    cat("Bayesian Estimation (multi-group)\n")
  }

  # Dataset info
  cat(paste0(" Number of variables: ", arguments$num_variables, "\n"))
  if(!is.null(arguments$num_groups)) {
    cat(paste0(" Number of groups: ", arguments$num_groups, "\n"))
    mapping = group_mapping_line(arguments, arguments$num_groups)
    if(!is.null(mapping)) cat(paste0(" ", mapping, "\n"))
  }
  if(!is.null(arguments$num_cases)) {
    # In our build_output_compare() we stored total cases in num_cases.
    if(isTRUE(arguments$na_impute)) {
      cat(paste0(" Number of cases: ", arguments$num_cases, " (missings imputed)\n"))
    } else {
      cat(paste0(" Number of cases: ", arguments$num_cases, "\n"))
    }
  }

  # Iterations and chains
  if(!is.null(arguments$num_chains)) {
    total_iter = arguments$iter * arguments$num_chains
    cat(paste0(" Number of post-burnin MCMC iterations: ", total_iter, "\n"))
    cat(paste0(" Number of MCMC chains: ", arguments$num_chains, "\n"))
  } else {
    cat(paste0(" Number of post-burnin MCMC iterations: ", arguments$iter, "\n"))
  }

  cat("Use the `summary()` function for posterior summaries and diagnostics.\n")
  cat("See the `easybgm` package for additional summaries and plotting.\n")
  invisible(x)
}


#' @name summary.bgmCompare
#' @title Summary method for `bgmCompare` objects
#'
#' @description Returns posterior summaries and diagnostics for a fitted `bgmCompare` model.
#'
#' @details In the printed "Group differences (main effects)" block, a
#'   threshold difference that rests on the prior rather than on the data is
#'   marked with a leading `*`, and the block gains the footnote *a group lacks
#'   observations in this category or in the reference category; the estimate
#'   reflects the prior, not the data*.
#'
#'   `bgmCompare()` keeps the union of the categories the groups observe, so a
#'   group can contribute no observations at all to a retained category. Two
#'   things follow. If the empty category is the row's own, that group's
#'   threshold for it has nothing behind it. If the empty category is the
#'   reference category, the group has no data fixing the level of its
#'   threshold vector at all, and *every* threshold of that variable is
#'   affected, not just one --- so the whole variable is marked. Either way the
#'   reported difference is large and very uncertain without being evidence of
#'   a group difference.
#'
#'   The per-group counts behind the mark are in
#'   `extract_arguments(fit)$category_support`, one matrix of category-by-group
#'   observation counts per variable (`NULL` for Blume-Capel variables, which
#'   are exempt from the union recode and are never marked). Fits made before
#'   that field existed print unmarked.
#'
#' @param object An object of class `bgmCompare`.
#' @param ... Currently ignored.
#'
#' @return An object of class `summary.bgmCompare` with posterior summaries.
#'
#' @examples
#' \donttest{
#' # See ?bgmCompare for a full example
#' }
#'
#' @seealso [bgmCompare()], [print.bgmCompare()], [coef.bgmCompare()]
#' @family posterior-methods
#'
#' @export
summary.bgmCompare = function(object, ...) {
  ensure_summaries(object)
  arguments = extract_arguments(object)

  if(!is.null(object$posterior_summary_main_baseline) &&
    !is.null(object$posterior_summary_pairwise_baseline)) {
    out = list(
      main      = object$posterior_summary_main_baseline,
      pairwise  = object$posterior_summary_pairwise_baseline
    )

    if(!is.null(object$posterior_summary_indicator)) {
      out$indicator = object$posterior_summary_indicator
    }

    if(!is.null(object$posterior_summary_main_differences)) {
      out$main_diff = object$posterior_summary_main_differences
    }

    if(!is.null(object$posterior_summary_pairwise_differences)) {
      out$pairwise_diff = object$posterior_summary_pairwise_differences
    }

    out$arguments = arguments
    class(out) = "summary.bgmCompare"
    return(out)
  }

  message(
    "No summary statistics available for this model object.\n",
    "Try fitting the model again using the latest bgms version,\n",
    "or use the `easybgm` package for diagnostic summaries and plotting."
  )
  invisible(NULL)
}


#' @export
print.summary.bgmCompare = function(x, digits = 3, ...) {
  cat("Posterior summaries from Bayesian grouped MRF estimation (bgmCompare):\n\n")

  mapping = group_mapping_line(x$arguments, x$arguments$num_groups)
  if(!is.null(mapping)) cat(paste0(mapping, "\n\n"))

  print_df = function(df, digits) {
    df2 = df
    if(ncol(df2) > 1) {
      df2[, -1] = lapply(df2[, -1, drop = FALSE], round, digits = digits)
    }
    # The labels live in the `parameter` column and, since the summary tables
    # gained the bgm row-name contract, in the row names as well. Suppress the
    # row names so the label prints once, as the blocks below already do.
    print(head(df2, 6), row.names = FALSE)
  }

  if(!is.null(x$main)) {
    cat("Category thresholds:\n")
    print_df(x$main, digits)
    if(nrow(x$main) > 6) {
      cat("... (use `summary(fit)$main` to see full output)\n")
    }
    cat("\n")
  }

  if(!is.null(x$pairwise)) {
    cat("Pairwise interactions:\n")
    print_df(x$pairwise, digits)
    if(nrow(x$pairwise) > 6) {
      cat("... (use `summary(fit)$pairwise` to see full output)\n")
    }
    cat("\n")
  }

  if(!is.null(x$indicator)) {
    cat("Inclusion probabilities:\n")
    # mean/mcse/sd/n_eff/Rhat are the Rao-Blackwellized inclusion estimate;
    # n0->1 / n1->0 are the raw directional flip counts, which record the
    # indicator's exploration beside them.
    ind = head(x$indicator, 6)

    ind_has_na = anyNA(ind)

    # round only numeric columns
    ind[] = lapply(ind, function(col) {
      if(is.numeric(col)) {
        round(col, digits)
      } else {
        col
      }
    })

    # replace NA with empty string for printing
    ind[] = lapply(ind, function(col) {
      ifelse(is.na(col), "", col)
    })

    print(ind, row.names = FALSE)
    if(nrow(x$indicator) > 6) {
      cat("... (use `summary(fit)$indicator` to see full output)\n")
    }
    if(ind_has_na) {
      cat("Note: blank mcse/n_eff/Rhat cells mark indicators whose inclusion draws never\n")
      cat("varied: the per-iteration evidence is so one-sided that the draws round to\n")
      cat("exactly 0 or 1, and no variation is left to estimate precision from. An\n")
      cat("indicator that is almost as certain still shows numbers; the difference is\n")
      cat("rounding, not evidence.\n")
      if(isFALSE(x$arguments$main_difference_selection)) {
        cat("Main-effect difference rows are blank when main_difference_selection = FALSE\n")
        cat("left their indicators unsampled.\n")
      }
      cat("All computed values remain in `summary(fit)$indicator`.\n")
    }
    cat("\n")
  }

  if(!is.null(x$main_diff)) {
    cat("Group differences (main effects):\n")

    maind = head(x$main_diff, 6)
    maind_has_na = anyNA(maind)

    # Only round numeric columns
    is_num = vapply(maind, is.numeric, logical(1))
    maind[is_num] = lapply(
      maind[is_num],
      function(col) ifelse(is.na(col), "", round(col, digits))
    )

    # Rows whose group-by-category cell is empty carry the prior, not the data.
    # The mark goes in front of the label rather than after the numbers: it
    # qualifies which parameter the row is, and a leading gutter keeps it from
    # reading as part of `mean`.
    prior_only = compare_prior_only_main_diff(x$arguments, nrow(x$main_diff))
    marked = if(is.null(prior_only) || is.null(maind$parameter)) {
      logical(nrow(maind))
    } else {
      prior_only[seq_len(nrow(maind))]
    }
    if(any(marked)) {
      maind$parameter = paste0(ifelse(marked, "* ", "  "), maind$parameter)
    }

    print(maind, row.names = FALSE)

    if(nrow(x$main_diff) > 6) {
      cat("... (use `summary(fit)$main_diff` to see full output)\n")
    }

    if(any(marked)) {
      cat(
        "* a group lacks observations in this category or in the reference",
        "category; the estimate reflects the prior, not the data\n"
      )
    }

    if(!is.null(x$indicator) && maind_has_na) {
      cat("Note: NA values are suppressed in the print table. They occur for differences\n")
      cat("that were never selected, so the composite ESS and share are undefined;\n")
      cat("`summary(fit)$main_diff` still contains the NA values.\n")
    }
    cat("\n")
  }

  if(!is.null(x$pairwise_diff)) {
    cat("Group differences (pairwise effects):\n")

    pairwised = head(x$pairwise_diff, 6)
    pairwised_has_na = anyNA(pairwised)

    # Only round numeric columns
    is_num = vapply(pairwised, is.numeric, logical(1))
    pairwised[is_num] = lapply(
      pairwised[is_num],
      function(col) ifelse(is.na(col), "", round(col, digits))
    )

    print(pairwised, row.names = FALSE)

    if(nrow(x$pairwise_diff) > 6) {
      cat("... (use `summary(fit)$pairwise_diff` to see full output)\n")
    }

    if(!is.null(x$indicator) && pairwised_has_na) {
      cat("Note: NA values are suppressed in the print table. They occur for differences\n")
      cat("that were never selected, so the composite ESS and share are undefined;\n")
      cat("`summary(fit)$pairwise_diff` still contains the NA values.\n")
    }
    cat("\n")
  }

  cat("Use `summary(fit)$<component>` to access full results.\n")
  cat("See the `easybgm` package for other summary and plotting tools.\n")
}


#' @title Extract Coefficients from a bgmCompare Object
#' @name coef.bgmCompare
#' @description Returns posterior means for raw parameters (baseline + differences)
#' and group-specific effects from a \code{bgmCompare} fit, as well as inclusion indicators.
#'
#' @param object An object of class \code{bgmCompare}.
#' @param ... Ignored.
#'
#' @return A list with components:
#' \describe{
#'   \item{main_effects_raw}{Posterior means of the raw main-effect parameters
#'   (variables x (baseline + differences)).}
#'   \item{pairwise_effects_raw}{Posterior means of the raw pairwise-effect parameters
#'   (pairs x (baseline + differences)).}
#'   \item{main_effects_groups}{Posterior means of group-specific main effects
#'   (variables x groups), computed as baseline plus projected differences.}
#'   \item{pairwise_effects_groups}{Posterior means of group-specific pairwise effects
#'   (pairs x groups), computed as baseline plus projected differences.}
#'   \item{indicators}{Posterior mean inclusion probabilities as a symmetric matrix,
#'   with diagonals corresponding to main effects and off-diagonals to pairwise effects.}
#' }
#'
#' @examples
#' \donttest{
#' # See ?bgmCompare for a full example
#' }
#'
#' @seealso [bgmCompare()], [print.bgmCompare()], [summary.bgmCompare()]
#' @family posterior-methods
#'
#' @export
coef.bgmCompare = function(object, ...) {
  args = extract_arguments(object)
  raw = get_raw_samples(object)

  var_names = args$data_columnnames
  num_variables = as.integer(args$num_variables)

  # ---- baseline + group-specific main and pairwise effects ----
  gp = .compute_group_param_matrices(args, raw)
  main_mat = gp$main_mat
  pairwise_mat = gp$pairwise_mat
  main_effects_groups = gp$main_effects_groups
  pairwise_effects_groups = gp$pairwise_effects_groups

  # ============================================================
  # ---- indicators (present only if selection was used) ----
  indicators = NULL
  array3d_ind = samples_to_array3d(raw$indicator)
  if(!is.null(array3d_ind)) {
    mean_ind = apply(array3d_ind, 3, mean)

    # reconstruct VxV matrix using the sampler's interleaved order:
    # (1,1),(1,2),...,(1,V),(2,2),...,(2,V),...,(V,V)
    V = num_variables
    stopifnot(length(mean_ind) == V * (V + 1L) / 2L)

    ind_mat = matrix(0,
      nrow = V, ncol = V,
      dimnames = list(var_names, var_names)
    )
    pos = 1L
    for(i in seq_len(V)) {
      # diagonal (main indicator)
      ind_mat[i, i] = mean_ind[pos]
      pos = pos + 1L
      if(i < V) {
        for(j in (i + 1L):V) {
          val = mean_ind[pos]
          pos = pos + 1L
          ind_mat[i, j] = val
          ind_mat[j, i] = val
        }
      }
    }
    indicators = ind_mat
  }

  # ============================================================
  # ---- return both raw + group-specific ----
  list(
    main_effects_raw        = main_mat,
    pairwise_effects_raw    = pairwise_mat,
    main_effects_groups     = main_effects_groups,
    pairwise_effects_groups = pairwise_effects_groups,
    indicators              = indicators
  )
}


#' Access elements of a bgmCompare object
#'
#' @description Provides \code{$} access to S7 properties. Lazy
#'   \code{posterior_summary_*} properties trigger computation on first
#'   access via S7 property getters. Also supports legacy S3 list-based
#'   fit objects.
#'
#' @param x A \code{bgmCompare} object.
#' @param name Name of the element to access.
#'
#' @return The requested element.
#'
#' @method $ bgmCompare
#' @export
`$.bgmCompare` = function(x, name) {
  if(inherits(x, "S7_object")) {
    S7::prop(x, name)
  } else {
    if(startsWith(name, "posterior_summary_")) {
      cache = .subset2(x, "cache")
      if(!is.null(cache)) {
        ensure_summaries(x)
        val = cache[[name]]
        if(!is.null(val)) {
          return(val)
        }
      }
    }
    .subset2(x, name)
  }
}


#' @rdname cash-.bgmCompare
#' @param ... Ignored.
#' @method [[ bgmCompare
#' @export
`[[.bgmCompare` = function(x, name, ...) {
  if(inherits(x, "S7_object")) {
    if(is.character(name)) {
      S7::prop(x, name)
    } else {
      stop("numeric indexing is not supported for bgmCompare objects")
    }
  } else {
    if(is.character(name) && startsWith(name, "posterior_summary_")) {
      cache = .subset2(x, "cache")
      if(!is.null(cache)) {
        ensure_summaries(x)
        val = cache[[name]]
        if(!is.null(val)) {
          return(val)
        }
      }
    }
    .subset2(x, name)
  }
}


#' @method names bgmCompare
#' @export
names.bgmCompare = function(x) {
  if(inherits(x, "S7_object")) {
    S7::prop(x, ".field_names")
  } else {
    NextMethod()
  }
}

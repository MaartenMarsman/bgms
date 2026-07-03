# R/methods_bgmcompare.R

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

  print_df = function(df, digits) {
    df2 = df
    if(ncol(df2) > 1) {
      df2[, -1] = lapply(df2[, -1, drop = FALSE], round, digits = digits)
    }
    print(head(df2, 6))
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
    ind = head(x$indicator, 6)

    # Suppress n_eff_mixt where fewer than 5 transitions observed
    if(all(c("n0->1", "n1->0", "n_eff_mixt") %in% names(ind))) {
      few = ind[["n0->1"]] + ind[["n1->0"]] < 5
      few[is.na(few)] = TRUE
      ind[["n_eff_mixt"]][few] = NA
    }

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
    cat("Note: NA values are suppressed in the print table. They occur when an indicator\n")
    cat("was constant or had fewer than 5 transitions, so n_eff_mixt is unreliable;\n")
    cat("`summary(fit)$indicator` still contains all computed values.\n\n")
  }

  if(!is.null(x$main_diff)) {
    cat("Group differences (main effects):\n")

    maind = head(x$main_diff, 6)

    # Only round numeric columns
    is_num = vapply(maind, is.numeric, logical(1))
    maind[is_num] = lapply(
      maind[is_num],
      function(col) ifelse(is.na(col), "", round(col, digits))
    )

    print(maind, row.names = FALSE)

    if(nrow(x$main_diff) > 6) {
      cat("... (use `summary(fit)$main_diff` to see full output)\n")
    }

    if(!is.null(x$indicator)) {
      cat("Note: NA values are suppressed in the print table. They occur here when an\n")
      cat("indicator was zero across all iterations, so mcse/n_eff/n_eff_mixt/Rhat are undefined;\n")
      cat("`summary(fit)$main_diff` still contains the NA values.\n")
    }
    cat("\n")
  }

  if(!is.null(x$pairwise_diff)) {
    cat("Group differences (pairwise effects):\n")

    pairwised = head(x$pairwise_diff, 6)

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

    if(!is.null(x$indicator)) {
      cat("Note: NA values are suppressed in the print table. They occur here when an\n")
      cat("indicator was zero across all iterations, so mcse/n_eff/n_eff_mixt/Rhat are undefined;\n")
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
#' @keywords internal
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
#' @keywords internal
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
#' @keywords internal
names.bgmCompare = function(x) {
  if(inherits(x, "S7_object")) {
    S7::prop(x, ".field_names")
  } else {
    NextMethod()
  }
}

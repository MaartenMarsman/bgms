# Number of rows shown in summary print previews before truncating with "...".
.summary_preview_rows = 6L


#' @name print.bgms
#' @title Print method for `bgms` objects
#'
#' @description Minimal console output for `bgms` fit objects.
#'
#' @param x An object of class `bgms`.
#' @param ... Ignored.
#'
#' @return Invisibly returns `x`.
#'
#' @examples
#' \donttest{
#' fit = bgm(x = Wenchuan[, 1:3])
#' print(fit)
#' }
#'
#' @seealso [bgm()], [summary.bgms()], [coef.bgms()]
#' @family posterior-methods
#'
#' @export
print.bgms = function(x, ...) {
  arguments = extract_arguments(x)

  # Estimation method
  if(isTRUE(arguments$edge_selection)) {
    prior_msg = switch(arguments$edge_prior,
      "Bernoulli" = "Bayesian Edge Selection using a Bernoulli prior on edge inclusion",
      "Beta-Bernoulli" = "Bayesian Edge Selection using a Beta-Bernoulli prior on edge inclusion",
      "Stochastic-Block" = "Bayesian Edge Selection using a Stochastic Block prior on edge inclusion",
      "Bayesian Edge Selection"
    )
    cat(prior_msg, "\n")
  } else {
    cat("Bayesian Estimation\n")
  }

  # Model type
  mt = arguments$model_type
  if(!is.null(mt)) {
    mt_label = switch(mt,
      ggm = "GGM (Gaussian Graphical Model)",
      omrf = "OMRF (Ordinal Markov Random Field)",
      mixed_mrf = sprintf(
        "Mixed MRF (%d discrete, %d continuous)",
        arguments$num_discrete, arguments$num_continuous
      ),
      mt
    )
    cat(paste0(" Model: ", mt_label, "\n"))
  }

  # Dataset info
  cat(paste0(" Number of variables: ", arguments$num_variables, "\n"))
  if(isTRUE(arguments$na_impute)) {
    cat(paste0(" Number of cases: ", arguments$num_cases, " (missings imputed)\n"))
  } else {
    cat(paste0(" Number of cases: ", arguments$num_cases, "\n"))
  }

  # Iterations and chains
  if(!is.null(arguments$num_chains)) {
    total_iter = arguments$iter * arguments$num_chains
    cat(paste0(" Number of post-burnin MCMC iterations: ", total_iter, "\n"))
    cat(paste0(" Number of MCMC chains: ", arguments$num_chains, "\n"))
  } else {
    cat(paste0(" Number of post-burnin MCMC iterations: ", arguments$iter, "\n"))
  }

  cat("Use the `summary()` function for posterior summaries and chain diagnostics.\n")
  cat("See the `easybgm` package for summary and plotting tools.\n")
  invisible(x)
}


#' @name summary.bgms
#' @title Summary method for `bgms` objects
#'
#' @description Returns posterior summaries and diagnostics for a fitted `bgms` model.
#'
#' @param object An object of class `bgms`.
#' @param ... Currently ignored.
#'
#' @return An object of class `summary.bgms` with posterior summaries.
#'
#' @examples
#' \donttest{
#' fit = bgm(x = Wenchuan[, 1:3])
#' summary(fit)
#' }
#'
#' @seealso [bgm()], [print.bgms()], [coef.bgms()]
#' @family posterior-methods
#'
#' @export
summary.bgms = function(object, ...) {
  ensure_summaries(object)
  arguments = extract_arguments(object)

  has_main = !is.null(object$posterior_summary_main)
  has_quad = !is.null(object$posterior_summary_quadratic)
  has_pair = !is.null(object$posterior_summary_pairwise)

  if((has_main || has_quad) && has_pair) {
    mt = arguments$model_type
    main_label = switch(mt,
      ggm       = NULL,
      omrf      = "Category thresholds:",
      mixed_mrf = "Main effects (discrete thresholds and continuous means):",
      "Main effects:"
    )
    out = list(
      main = object$posterior_summary_main,
      quadratic = object$posterior_summary_quadratic,
      pairwise = object$posterior_summary_pairwise
    )
    attr(out, "main_label") = main_label

    if(!is.null(object$posterior_summary_indicator)) {
      out$indicator = object$posterior_summary_indicator
    }

    if(!is.null(object$posterior_summary_pairwise_allocations)) {
      out$allocations = object$posterior_summary_pairwise_allocations
      out$mean_allocations = get_posterior_mean(object, "allocations")
      out$mode_allocations = object$posterior_mode_allocations
      out$num_blocks = object$posterior_num_blocks
    }

    class(out) = "summary.bgms"
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
print.summary.bgms = function(x, digits = 3, ...) {
  cat("Posterior summaries from Bayesian estimation:\n\n")

  if(!is.null(x$main)) {
    main_label = attr(x, "main_label") %||% "Main effects:"
    cat(main_label, "\n")
    print(round(head(x$main, .summary_preview_rows), digits = digits))
    if(nrow(x$main) > .summary_preview_rows) cat("... (use `summary(fit)$main` to see full output)\n")
    cat("\n")
  }

  if(!is.null(x$quadratic)) {
    cat("Residual variances:\n")
    print(round(head(x$quadratic, .summary_preview_rows), digits = digits))
    if(nrow(x$quadratic) > .summary_preview_rows) cat("... (use `summary(fit)$quadratic` to see full output)\n")
    cat("\n")
  }

  if(!is.null(x$pairwise)) {
    cat("Pairwise interactions:\n")
    pair = head(x$pairwise, .summary_preview_rows)
    pair_has_na = anyNA(pair)
    pair[] = lapply(pair, function(col) ifelse(is.na(col), "", round(col, digits)))
    print(pair)
    if(nrow(x$pairwise) > .summary_preview_rows) cat("... (use `summary(fit)$pairwise` to see full output)\n")
    if(!is.null(x$indicator) && pair_has_na) {
      cat("Note: NA values are suppressed in the print table. They occur when an indicator\n")
      cat("is constant across iterations, so mcse/n_eff_mixt are undefined;\n")
      cat("`summary(fit)$pairwise` still contains the NA values.\n")
    }
    cat("\n")
  }

  if(!is.null(x$indicator)) {
    cat("Inclusion probabilities:\n")
    ind = head(x$indicator, .summary_preview_rows)
    # Suppress n_eff_mixt where fewer than 5 transitions observed
    if(all(c("n0->1", "n1->0", "n_eff_mixt") %in% names(ind))) {
      few = ind[["n0->1"]] + ind[["n1->0"]] < 5
      few[is.na(few)] = TRUE
      ind[["n_eff_mixt"]][few] = NA
    }
    ind_has_na = anyNA(ind)
    ind[] = lapply(ind, function(col) ifelse(is.na(col), "", round(col, digits)))
    print(ind)
    if(nrow(x$indicator) > .summary_preview_rows) cat("... (use `summary(fit)$indicator` to see full output)\n")
    if(ind_has_na) {
      cat("Note: NA values are suppressed in the print table. They occur when an indicator\n")
      cat("was constant or had fewer than 5 transitions, so n_eff_mixt is unreliable;\n")
      cat("`summary(fit)$indicator` still contains all computed values.\n")
    }
    cat("\n")
  }

  if(!is.null(x$allocations)) {
    cat("Pairwise node co-clustering proportion:\n")
    print(round(head(x$allocations, .summary_preview_rows), digits = digits))
    if(nrow(x$allocations) > .summary_preview_rows) cat("... (use `summary(fit)$allocations` to see full output)\n")
    cat("\n")
  }

  if(!is.null(x$mean_allocations)) {
    cat("Mean posterior node allocation vector:\n")
    print(round(head(x$mean_allocations, .summary_preview_rows), digits = digits))
    cat("Mode posterior node allocation vector:\n")
    print(round(head(x$mode_allocations, .summary_preview_rows), digits = digits))
    cat("\n")
  }

  if(!is.null(x$num_blocks)) {
    cat("Number of blocks and their posterior probability :\n")
    print(round(head(x$num_blocks, .summary_preview_rows), digits = digits))
    if(nrow(x$num_blocks) > .summary_preview_rows) cat("... (use `summary(fit)$num_blocks` to see full output)\n")
    cat("\n")
  }

  cat("Use `summary(fit)$<component>` to access full results.\n")
  if(!is.null(x$quadratic) && !is.null(x$main)) {
    cat("Use `extract_partial_correlations(fit)` for partial correlations,\n")
    cat("`extract_precision(fit)` for the precision matrix,\n")
    cat("and `extract_log_odds(fit)` for log odds ratios.\n")
  } else if(!is.null(x$quadratic)) {
    cat("Use `extract_partial_correlations(fit)` for partial correlations,\n")
    cat("and `extract_precision(fit)` for the precision matrix.\n")
  } else if(!is.null(x$main)) {
    cat("Use `extract_log_odds(fit)` for log odds ratios.\n")
  }
  cat("See the `easybgm` package for other summary and plotting tools.\n")
}


#' @title Extract Coefficients from a bgms Object
#' @name coef.bgms
#' @description Returns the posterior mean main effects, pairwise effects, and edge inclusion indicators from a \code{bgms} model fit.
#'
#' @param object An object of class \code{bgms}.
#' @param ... Ignored.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{main}{Posterior mean of the main-effect parameters. \code{NULL} for
#'     GGM models (no main effects). For OMRF models this is a numeric matrix
#'     (p x max_categories) of category thresholds. For mixed MRF models this
#'     is a list with \code{$discrete} (p x max_categories matrix) and
#'     \code{$continuous} (q x 1 matrix of means).}
#'   \item{pairwise}{Posterior mean of the partial association matrix (zero
#'     diagonal). Use [extract_precision()] for the full precision matrix.}
#'   \item{indicator}{Posterior mean of the edge inclusion indicators (if available).}
#'   \item{mean_allocations}{Posterior mean block allocations. Present only for
#'     Stochastic-Block edge-prior fits.}
#'   \item{mode_allocations}{Posterior mode block allocations. Present only for
#'     Stochastic-Block edge-prior fits.}
#'   \item{num_blocks}{Data frame of the posterior distribution over the number
#'     of blocks. Present only for Stochastic-Block edge-prior fits.}
#' }
#'
#' @examples
#' \donttest{
#' fit = bgm(x = Wenchuan[, 1:3])
#' coef(fit)
#' }
#'
#' @seealso [bgm()], [print.bgms()], [summary.bgms()]
#' @family posterior-methods
#'
#' @export
coef.bgms = function(object, ...) {
  out = list(
    main = get_posterior_mean(object, "main"),
    pairwise = get_posterior_mean(object, "pairwise")
  )
  if(!is.null(get_posterior_mean(object, "indicator"))) {
    out$indicator = get_posterior_mean(object, "indicator")
  }

  if(!is.null(get_posterior_mean(object, "allocations"))) {
    out$mean_allocations = get_posterior_mean(object, "allocations")
    out$mode_allocations = object$posterior_mode_allocations
    out$num_blocks = object$posterior_num_blocks
  }

  return(out)
}


#' Access elements of a bgms object
#'
#' @description Provides \code{$} access to S7 properties. Lazy
#'   \code{posterior_summary_*} properties trigger computation on first
#'   access via S7 property getters. Also supports legacy S3 list-based
#'   fit objects.
#'
#' @param x A \code{bgms} object.
#' @param name Name of the element to access.
#'
#' @return The requested element.
#'
#' @method $ bgms
#' @export
`$.bgms` = function(x, name) {
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


#' @rdname cash-.bgms
#' @param ... Ignored.
#' @method [[ bgms
#' @export
`[[.bgms` = function(x, name, ...) {
  if(inherits(x, "S7_object")) {
    if(is.character(name)) {
      S7::prop(x, name)
    } else {
      stop("numeric indexing is not supported for bgms objects")
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


#' @method names bgms
#' @export
names.bgms = function(x) {
  if(inherits(x, "S7_object")) {
    S7::prop(x, ".field_names")
  } else {
    NextMethod()
  }
}


.warning_state = new.env(parent = emptyenv())
.warning_state$issued = FALSE

warning_once = function(msg) {
  if(!.warning_state$issued) {
    warning(msg, call. = FALSE)
    .warning_state$issued = TRUE
  }
}

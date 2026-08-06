# ==============================================================================
# Extractor Functions - S3 Generics and Methods
# ==============================================================================


# ------------------------------------------------------------------------------
# samples_to_array3d()
# ------------------------------------------------------------------------------
# Stack a per-chain list of [iter x param] sample matrices into an
# [iter x chain x param] array. Returns NULL for a NULL input (e.g. absent
# indicator samples) and coerces any non-matrix element to a one-column
# matrix. Several bgmCompare extractors defined this as a local closure; the
# variants had drifted (only some handled NULL / vector inputs), so this is
# the single hardened version.
#
# @param xlist  List of per-chain matrices, or NULL.
#
# Returns: [iter x chain x param] numeric array, or NULL.
# ------------------------------------------------------------------------------
samples_to_array3d = function(xlist) {
  if(is.null(xlist)) {
    return(NULL)
  }
  stopifnot(length(xlist) >= 1)
  mats = lapply(xlist, function(x) {
    m = as.matrix(x)
    if(is.null(dim(m))) m = matrix(m, ncol = 1L)
    m
  })
  niter = nrow(mats[[1]])
  nparam = ncol(mats[[1]])
  arr = array(NA_real_, dim = c(niter, length(mats), nparam))
  for(c in seq_along(mats)) arr[, c, ] = mats[[c]]
  arr
}


#' @title Extract Model Arguments
#'
#' @description
#' Retrieves the arguments used when fitting a model with [bgm()] or
#' [bgmCompare()].
#'
#' @param bgms_object A fitted model object of class `bgms` (from [bgm()])
#'   or `bgmCompare` (from [bgmCompare()]).
#'
#' @return A named list containing all arguments passed to the fitting
#'   function, including data dimensions, prior settings, and MCMC
#'   configuration.
#'
#'   For `bgmCompare` fits the list additionally carries
#'   `main_effect_indices`: an integer matrix with one row per variable and
#'   two columns giving the zero-based first and last column of that
#'   variable's block in the baseline main-effect parameters returned by
#'   [extract_main_effects()]. Variables do not occupy a fixed number of
#'   columns -- an ordinal variable contributes one per category and a
#'   Blume-Capel variable two -- so this layout is what maps parameter
#'   columns back to variables. It is not available for `bgms` fits, whose
#'   main effects are already returned as one row per variable.
#'
#' @seealso [bgm()], [bgmCompare()], [summary.bgms()], [summary.bgmCompare()]
#' @family extractors
#' @export
extract_arguments = function(bgms_object) {
  UseMethod("extract_arguments")
}

#' @inheritParams extract_arguments
#' @exportS3Method
#' @noRd
extract_arguments.bgms = function(bgms_object) {
  if(is.null(bgms_object$arguments)) {
    stop("Fit object predates bgms version 0.1.3. Upgrade the model output.")
  }
  return(bgms_object$arguments)
}

#' @inheritParams extract_arguments
#' @exportS3Method
#' @noRd
extract_arguments.bgmCompare = function(bgms_object) {
  arguments = bgms_object$arguments
  if(is.null(arguments)) {
    stop("Fit object predates bgms version 0.1.3. Upgrade the model output.")
  }

  # The main-effect row layout is built by the spec and kept in the internal
  # cache, never in $arguments, so callers that need to map main-effect
  # parameter rows back to variables had no way to get at it. Surface it here
  # rather than in the stored object: it is derived, and $arguments is the
  # documented place users look. Fits written before the cache carried it fall
  # through unchanged.
  if(is.null(arguments$main_effect_indices)) {
    cached = tryCatch(bgms_object$cache$main_effect_indices, error = function(e) NULL)
    if(!is.null(cached)) {
      arguments$main_effect_indices = cached
    }
  }

  return(arguments)
}

#' @title Extract Indicator Samples
#'
#' @description
#' Retrieves posterior samples of inclusion indicators from a model fitted
#' with [bgm()] (edge inclusion indicators) or [bgmCompare()] (difference
#' indicators).
#'
#' @param bgms_object A fitted model object of class `bgms` (from [bgm()])
#'   or `bgmCompare` (from [bgmCompare()]).
#'
#' @return A matrix with one row per post-warmup iteration and one column per
#'   indicator, containing binary (0/1) samples.
#'   \describe{
#'     \item{bgms}{One column per edge. Requires `edge_selection = TRUE`.}
#'     \item{bgmCompare}{Columns for main-effect and pairwise difference
#'       indicators. Requires `difference_selection = TRUE`.}
#'   }
#'
#' @seealso [bgm()], [bgmCompare()],
#'   [extract_posterior_inclusion_probabilities()]
#' @family extractors
#' @export
extract_indicators = function(bgms_object) {
  UseMethod("extract_indicators")
}

#' @inheritParams extract_indicators
#' @exportS3Method
#' @noRd
extract_indicators.bgms = function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  if(!isTRUE(arguments$edge_selection)) {
    stop("To access edge indicators, the model must be run with edge_selection = TRUE.")
  }

  # Current format (0.1.6.0+)
  raw = get_raw_samples(bgms_object)
  if(!is.null(raw$indicator)) {
    indicators_list = raw$indicator
    indicator_samples = do.call(rbind, indicators_list)
    param_names = raw$parameter_names$indicator
    stopifnot("parameter_names$indicator missing in fit object" = !is.null(param_names))
    colnames(indicator_samples) = param_names
    return(indicator_samples)
  }

  # Deprecated format (0.1.4--0.1.5): $indicator stored at top level
  if(!is.null(bgms_object$indicator)) {
    lifecycle::deprecate_warn(
      "0.1.6.0",
      I("The '$indicator' field is deprecated; please refit with bgms >= 0.1.6.0")
    )
    return(bgms_object$indicator)
  }

  # Defunct format (pre-0.1.4): $gamma field
  lifecycle::deprecate_stop(
    "0.1.4",
    I("The '$gamma' field is defunct; please refit with bgms >= 0.1.6.0")
  )
}

#' @inheritParams extract_indicators
#' @exportS3Method
#' @noRd
extract_indicators.bgmCompare = function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  if(!isTRUE(arguments$difference_selection)) {
    stop("To access difference indicators, the model must be run with difference_selection = TRUE.")
  }

  # Current format (0.1.6.0+)
  raw = get_raw_samples(bgms_object)
  if(!is.null(raw$indicator)) {
    indicator_samples = do.call(rbind, raw$indicator)
    param_names = raw$parameter_names$indicators
    if(!is.null(param_names)) {
      colnames(indicator_samples) = param_names
    }
    return(indicator_samples)
  }

  # Deprecated format (0.1.4--0.1.5): $pairwise_difference_indicator at top level
  if(!is.null(bgms_object$pairwise_difference_indicator)) {
    lifecycle::deprecate_warn(
      "0.1.6.0",
      I("The '$pairwise_difference_indicator' field is deprecated; please refit with bgms >= 0.1.6.0")
    )
    return(bgms_object$pairwise_difference_indicator)
  }

  stop("No indicator samples found in fit object.")
}

#' @title Extract Posterior Inclusion Probabilities
#'
#' @description
#' Computes posterior inclusion probabilities from a model fitted with
#' [bgm()] (edge inclusion) or [bgmCompare()] (difference inclusion).
#'
#' Two estimators of the same posterior inclusion probability are available
#' through `estimator`. The default `"rb"` (Rao-Blackwellized) averages the
#' one-step draw
#' \eqn{J_t = \gamma_t + (1 - 2 \gamma_t)\,\alpha_t}, where \eqn{\gamma_t} is
#' the indicator state before the move and \eqn{\alpha_t} is the acceptance
#' probability of the birth/death proposal; `"raw"` instead averages the
#' indicator draws. Averaging \eqn{J_t} is a
#' lower-variance estimator and in exact arithmetic lies strictly inside
#' \eqn{(0, 1)}, so even indicators whose raw average saturates at 0 or 1
#' receive an interior estimate. In double precision, however, the average of
#' the stored \eqn{J_t} draws still rounds to exactly 0 or 1 for edges with
#' overwhelming per-iteration evidence, because \eqn{1 - \alpha_t} underflows
#' once \eqn{\alpha_t} drops below about `1e-16`. For inclusion Bayes factors,
#' use [extract_inclusion_bf()], which accumulates the odds on the
#' acceptance-probability scale and stays finite far beyond that ceiling; its
#' `log = TRUE` return carries evidence beyond what the Bayes factor scale can
#' represent in double precision. The
#' `"rb"` estimator changes only the summary, not the sampler; it inherits the
#' chain's mixing, does not rescue a chain that has failed to explore the model
#' space, and requires a fit from bgms >= 0.2.0.0. Because the RB draw is
#' continuous, the standard MCSE/ESS/split-R-hat machinery applies to it, and
#' the fit summary's inclusion table reports `mcse`, `n_eff`, and `Rhat` on the
#' RB draws. They quantify precision *conditional on exploration*: a stuck chain
#' can show a beautifully converged `J` chain with a high `n_eff`, so read them
#' beside the per-direction flip counts (`n0->1`, `n1->0`), which record the
#' exploration itself and whose asymmetry no symmetric summary recovers. The RB
#' draws vary on almost every edge, including edges whose indicator never
#' flipped; the three columns are `NA` only where the RB draws are constant to
#' double precision, which places the inclusion probability at its numerical
#' bound and the verdict beyond any threshold. All of this is in the fit
#' summary's inclusion table, `summary(fit)$indicator`.
#'
#' @param bgms_object A fitted model object of class `bgms` (from [bgm()])
#'   or `bgmCompare` (from [bgmCompare()]).
#' @param estimator Character; which estimator of the posterior inclusion
#'   probability to return. `"rb"` (default) returns the lower-variance
#'   Rao-Blackwellized average; `"raw"` averages the indicator draws.
#'
#' @return A symmetric p x p matrix of posterior inclusion probabilities,
#'   with variable names as row and column names.
#'   \describe{
#'     \item{bgms}{Off-diagonal entries are edge inclusion probabilities.
#'       Requires `edge_selection = TRUE`.}
#'     \item{bgmCompare}{Diagonal entries are main-effect inclusion
#'       probabilities; off-diagonal entries are pairwise difference
#'       inclusion probabilities. Requires `difference_selection = TRUE`.
#'       With `estimator = "rb"`, indicators that were not selected (e.g.
#'       main-effect differences when `main_difference_selection = FALSE`)
#'       are returned as `NA`.}
#'   }
#'
#' @seealso [extract_inclusion_bf()], [bgm()], [bgmCompare()],
#'   [extract_indicators()]
#' @family extractors
#' @export
extract_posterior_inclusion_probabilities = function(bgms_object,
                                                     estimator = c("rb", "raw")) {
  UseMethod("extract_posterior_inclusion_probabilities")
}

#' @inheritParams extract_posterior_inclusion_probabilities
#' @exportS3Method
#' @noRd
extract_posterior_inclusion_probabilities.bgms = function(bgms_object,
                                                          estimator = c("rb", "raw")) {
  estimator_missing = missing(estimator)
  arguments = extract_arguments(bgms_object)

  if(!isTRUE(arguments$edge_selection)) {
    stop("To estimate posterior inclusion probabilities, run bgm() with edge_selection = TRUE.")
  }

  # Handle legacy field name (no_variables -> num_variables in 0.1.6.0)
  num_vars = arguments$num_variables %||% arguments$no_variables
  data_columnnames = arguments$data_columnnames

  raw = get_raw_samples(bgms_object)
  # Default to the Rao-Blackwellized estimate, but fall back to the raw
  # indicator average for fits without RB draws (bgms < 0.2.0.0), so reading a
  # legacy fit with the default call still works; an explicit estimator = "rb"
  # errors below instead.
  estimator = if(estimator_missing && is.null(raw$rb_inclusion)) {
    "raw"
  } else {
    match.arg(estimator)
  }
  if(estimator == "rb") {
    # Rao-Blackwellized average of the stored one-step draws J.
    if(is.null(raw$rb_inclusion)) {
      stop("No Rao-Blackwellized inclusion draws found in fit object; refit with bgms >= 0.2.0.0.")
    }
    edge_means = colMeans(do.call(rbind, raw$rb_inclusion), na.rm = TRUE)
  } else if(!is.null(raw$indicator)) {
    # Current format (0.1.6.0+): raw indicator average.
    edge_means = colMeans(extract_indicators(bgms_object))
  } else if(!is.null(bgms_object$indicator)) {
    # Deprecated format (0.1.4--0.1.5): $indicator at top level
    lifecycle::deprecate_warn(
      "0.1.6.0",
      I("The '$indicator' field is deprecated; please refit with bgms >= 0.1.6.0")
    )
    edge_means = colMeans(bgms_object$indicator)
  } else {
    # Defunct format (pre-0.1.4)
    lifecycle::deprecate_stop(
      "0.1.4.2",
      I("The '$gamma' field is defunct; please refit with bgms >= 0.1.6.0")
    )
  }

  # Mixed-MRF fits store indicators in block order (discrete-discrete,
  # continuous-continuous, cross) over internally reordered variables, not
  # in global pair order: map them through the same block filler as
  # posterior_mean_indicator.
  spec = get_fit_spec(bgms_object)
  if(!is.null(spec) && identical(spec$model_type, "mixed_mrf")) {
    d = spec$data
    return(fill_mixed_symmetric(
      edge_means, d$num_discrete, d$num_continuous,
      d$discrete_indices, d$continuous_indices,
      list(data_columnnames, data_columnnames)
    ))
  }

  pip_matrix = matrix(0, num_vars, num_vars)
  pip_matrix[lower.tri(pip_matrix)] = edge_means
  pip_matrix = pip_matrix + t(pip_matrix)

  colnames(pip_matrix) = data_columnnames
  rownames(pip_matrix) = data_columnnames

  return(pip_matrix)
}


# ------------------------------------------------------------------
# .rb_log_odds_from_counts (internal)
# ------------------------------------------------------------------
# Per-edge log posterior inclusion odds from the RB odds accumulators,
# pooled over chains. counts is a per-chain list of (n_edges x 4) matrices
# with columns [n01, n10, n0_visits, n1_visits] on the acceptance-probability
# scale. Uses the exact identity
#   mean(J) / (1 - mean(J))
#     = (n01 + n1_visits - n10) / (n0_visits - n01 + n10),
# which avoids forming 1 - alpha per draw and so stays finite down to
# log-acceptances of about -745. Returns log(odds) (natural log): NA where an
# edge was never updated, +Inf when the denominator is exactly zero, -Inf when
# the numerator is exactly zero.
# ------------------------------------------------------------------
rb_log_odds_from_counts = function(counts) {
  pooled = Reduce(`+`, counts)
  n01 = pooled[, 1]
  n10 = pooled[, 2]
  n0_visits = pooled[, 3]
  n1_visits = pooled[, 4]

  num = n01 + n1_visits - n10
  den = n0_visits - n01 + n10

  out = log(num) - log(den)
  out[num == 0] = -Inf
  out[den == 0] = Inf
  out[(n0_visits + n1_visits) == 0] = NA_real_
  out
}


# ------------------------------------------------------------------
# rb_bf_scale
# ------------------------------------------------------------------
# Puts log inclusion Bayes factors on the scale requested by
# extract_inclusion_bf()'s log argument.
#
# @param log_bf  Matrix of natural-log inclusion Bayes factors.
# @param log     Logical. TRUE keeps the log scale, FALSE exponentiates.
#
# Returns: log_bf unchanged, or exp(log_bf), where -Inf maps to 0, +Inf and NA
# are preserved, and values above about 709.78 nats overflow to +Inf.
# ------------------------------------------------------------------
rb_bf_scale = function(log_bf, log) {
  if(length(log) != 1L || !is.logical(log) || is.na(log)) {
    stop("The log argument must be a single logical value, but not NA.")
  }
  if(log) log_bf else exp(log_bf)
}


#' @title Extract Rao-Blackwellized Inclusion Bayes Factors
#'
#' @description
#' Computes inclusion Bayes factors from a model fitted with [bgm()] (edge
#' inclusion) or [bgmCompare()] (difference inclusion), using the
#' Rao-Blackwellized odds accumulators recorded during sampling. For each
#' indicator the sampler sums the birth/death acceptance probability on the
#' acceptance-probability scale, so the posterior inclusion odds follow from
#' the exact identity
#' \deqn{\frac{\bar{J}}{1 - \bar{J}}
#'   = \frac{n_{01} + n_{1} - n_{10}}{n_{0} - n_{01} + n_{10}},}
#' where \eqn{n_{01}} and \eqn{n_{10}} sum the acceptance probabilities of birth
#' and death proposals and \eqn{n_0}, \eqn{n_1} count them. Because \eqn{1 -
#' \alpha} is never formed per draw, the odds stay finite down to log
#' acceptances of about -745, so edges that saturate the naive average of the
#' RB draws (which rounds to 0 or 1 near the boundary) still receive a finite
#' Bayes factor here.
#'
#' The prior inclusion odds are removed edge by edge, so the returned value is
#' the inclusion Bayes factor rather than the posterior odds: the two coincide
#' only at a prior inclusion probability of \eqn{1/2} (the default). For `bgm()`
#' fits the prior odds come from [extract_prior_inclusion_probabilities()]; for
#' continuous or stochastic-block models that call may run and cache a short
#' prior-only chain. For `bgmCompare()` fits the exchangeable difference prior
#' supplies a single prior inclusion probability (Bernoulli or Beta-Bernoulli);
#' a stochastic-block difference prior has no single marginal, so the result is
#' posterior odds there. The `log` argument applies to that return unchanged.
#'
#' The accumulators are exact on the log scale everywhere, while the Bayes
#' factor scale saturates at double precision: an entry whose log exceeds about
#' 709.78 nats (a Bayes factor beyond about 1.8e308) is `+Inf` under
#' `log = FALSE` even though its log-scale value is finite. Use `log = TRUE` for
#' workflows that must separate such extreme evidence.
#'
#' @param bgms_object A fitted model object of class `bgms` (from [bgm()])
#'   or `bgmCompare` (from [bgmCompare()]).
#' @param log Logical. If `FALSE` (default), return inclusion Bayes factors; if
#'   `TRUE`, return their natural logarithm.
#'
#' @return A symmetric p x p matrix of inclusion Bayes factors, or of their
#'   natural logarithms when `log = TRUE`, with variable names as row and column
#'   names. Entries are `NA` for indicators that were never updated, `+Inf` when
#'   no exclusion evidence remains (denominator exactly zero), and `-Inf` when no
#'   inclusion evidence remains. On the Bayes factor scale the latter is `0`, and
#'   `+Inf` also covers entries that overflow double precision. For `bgms` the
#'   diagonal is `NA`; for `bgmCompare` the diagonal holds main-effect difference
#'   Bayes factors.
#'
#'   "Never updated" means never *proposed*, which is narrower than it sounds.
#'   An edge indicator in [bgm()] that stays included for the whole run is still
#'   proposed at every iteration, so the accumulators see its conditional
#'   inclusion odds and it gets a finite (possibly very large) Bayes factor. An
#'   unselected main-effect difference in [bgmCompare()]
#'   (`main_difference_selection = FALSE`) is never proposed at all, so no
#'   Rao-Blackwellized quantity exists for it and `NA` is the honest entry
#'   rather than a lost number.
#'
#' @examples
#' \donttest{
#' fit = bgm(x = Wenchuan[, 1:3])
#' extract_inclusion_bf(fit)
#'
#' # log = TRUE keeps evidence that saturates the Bayes factor scale readable.
#' extract_inclusion_bf(fit, log = TRUE)
#' }
#'
#' @seealso [extract_posterior_inclusion_probabilities()],
#'   [extract_prior_inclusion_probabilities()]
#' @family extractors
#' @export
extract_inclusion_bf = function(bgms_object, log = FALSE) {
  UseMethod("extract_inclusion_bf")
}

#' @inheritParams extract_inclusion_bf
#' @exportS3Method
#' @noRd
extract_inclusion_bf.bgms = function(bgms_object, log = FALSE) {
  arguments = extract_arguments(bgms_object)

  if(!isTRUE(arguments$edge_selection)) {
    stop("To estimate Rao-Blackwellized inclusion Bayes factors, run bgm() with edge_selection = TRUE.")
  }

  raw = get_raw_samples(bgms_object)
  if(is.null(raw$rb_counts)) {
    stop("No Rao-Blackwellized odds accumulators found in fit object; refit with bgms >= 0.2.0.0.")
  }

  num_vars = arguments$num_variables %||% arguments$no_variables
  data_columnnames = arguments$data_columnnames

  log_odds = rb_log_odds_from_counts(raw$rb_counts)

  spec = get_fit_spec(bgms_object)
  if(!is.null(spec) && identical(spec$model_type, "mixed_mrf")) {
    d = spec$data
    post_log_odds = fill_mixed_symmetric(
      log_odds, d$num_discrete, d$num_continuous,
      d$discrete_indices, d$continuous_indices,
      list(data_columnnames, data_columnnames)
    )
  } else {
    post_log_odds = matrix(NA_real_, num_vars, num_vars)
    post_log_odds[lower.tri(post_log_odds)] = log_odds
    post_log_odds[upper.tri(post_log_odds)] =
      t(post_log_odds)[upper.tri(post_log_odds)]
    colnames(post_log_odds) = data_columnnames
    rownames(post_log_odds) = data_columnnames
  }

  # Turn posterior inclusion odds into a Bayes factor by removing the prior
  # inclusion odds edge by edge. Under a Beta-Bernoulli, stochastic-block, or
  # user-set edge prior the prior odds are not 1, so posterior odds and the
  # Bayes factor differ; they coincide only at a prior inclusion probability of
  # 1/2. extract_prior_inclusion_probabilities() returns a matrix in the same
  # shape and orientation, so the subtraction is per-edge.
  prior_pip = extract_prior_inclusion_probabilities(bgms_object)
  prior_log_odds = log(prior_pip) - log1p(-prior_pip)
  rb_bf_scale(post_log_odds - prior_log_odds, log)
}

#' @inheritParams extract_inclusion_bf
#' @exportS3Method
#' @noRd
extract_inclusion_bf.bgmCompare = function(bgms_object, log = FALSE) {
  arguments = extract_arguments(bgms_object)

  if(!isTRUE(arguments$difference_selection)) {
    stop("To estimate Rao-Blackwellized inclusion Bayes factors, run bgmCompare() with difference_selection = TRUE.")
  }

  raw = get_raw_samples(bgms_object)
  if(is.null(raw$rb_counts)) {
    stop("No Rao-Blackwellized odds accumulators found in fit object; refit with bgms >= 0.2.0.0.")
  }

  var_names = arguments$data_columnnames
  num_variables = as.integer(arguments$num_variables %||% arguments$no_variables)

  log_odds = rb_log_odds_from_counts(raw$rb_counts)

  # Reconstruct the VxV matrix using the sampler's interleaved order:
  # (1,1),(1,2),...,(1,V),(2,2),...,(2,V),...,(V,V).
  V = num_variables
  stopifnot(length(log_odds) == V * (V + 1L) / 2L)

  bf_mat = matrix(NA_real_,
    nrow = V, ncol = V,
    dimnames = list(var_names, var_names)
  )
  pos = 1L
  for(i in seq_len(V)) {
    bf_mat[i, i] = log_odds[pos]
    pos = pos + 1L
    if(i < V) {
      for(j in (i + 1L):V) {
        val = log_odds[pos]
        pos = pos + 1L
        bf_mat[i, j] = val
        bf_mat[j, i] = val
      }
    }
  }

  # Remove the difference prior inclusion odds so the result is a Bayes factor
  # rather than posterior odds. The difference prior is exchangeable across
  # difference indicators, so one prior inclusion probability applies to every
  # entry (main-effect differences on the diagonal, pairwise differences off
  # it). A stochastic-block difference prior has no single marginal, so the
  # result is left as posterior odds there.
  prior_p = difference_prior_inclusion(bgms_object)
  if(length(prior_p) == 1L && !is.na(prior_p)) {
    bf_mat = bf_mat - (log(prior_p) - log1p(-prior_p))
  }

  return(rb_bf_scale(bf_mat, log))
}


# ------------------------------------------------------------------
# difference_prior_inclusion
# ------------------------------------------------------------------
# Marginal prior inclusion probability of a bgmCompare() difference indicator.
#
# The difference prior is exchangeable across difference indicators, so one
# probability applies to every entry, main-effect differences on the diagonal
# and pairwise differences off it. A stochastic-block difference prior has no
# single marginal and returns NA, which leaves its caller reporting posterior
# odds rather than a Bayes factor.
#
# @param bgms_object  A fitted bgmCompare object.
#
# Returns: a single probability, or NA_real_.
# ------------------------------------------------------------------
difference_prior_inclusion = function(bgms_object) {
  arguments = extract_arguments(bgms_object)
  switch(as.character(arguments$difference_prior),
    "Bernoulli" = {
      dp = arguments$inclusion_probability
      if(is.matrix(dp)) dp[upper.tri(dp)][1L] else dp[1L]
    },
    "Beta-Bernoulli" = arguments$difference_selection_alpha /
      (arguments$difference_selection_alpha +
        arguments$difference_selection_beta),
    NA_real_
  )
}


#' @title Extract Stochastic Block Model Summaries
#'
#' @description
#' Retrieves posterior summaries from a model fitted with the Stochastic
#' Block prior. Works on both `bgms` fits (where SBM governs edge inclusion)
#' and `bgmCompare` fits (where SBM governs the off-diagonal pairwise
#' difference inclusions).
#'
#' @param bgms_object A fitted model object of class `bgms` (from [bgm()])
#'   or `bgmCompare` (from [bgmCompare()]).
#'
#' @return A list with elements `posterior_num_blocks`,
#'   `posterior_mean_allocations`, `posterior_mode_allocations`, and
#'   `posterior_mean_coclustering_matrix`. For `bgms`, requires
#'   `edge_selection = TRUE` and `edge_prior = sbm_prior(...)`. For
#'   `bgmCompare`, requires `difference_selection = TRUE` and
#'   `difference_prior = sbm_prior(...)`.
#'
#' @seealso [bgm()], [bgmCompare()], [extract_indicators()],
#'   [extract_posterior_inclusion_probabilities()]
#' @family extractors
#' @export
extract_sbm = function(bgms_object) {
  UseMethod("extract_sbm")
}

#' @inheritParams extract_sbm
#' @exportS3Method
#' @noRd
extract_sbm.bgms = function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  if(!isTRUE(arguments$edge_selection)) {
    stop("To extract SBM summaries, run bgm() with edge_selection = TRUE.")
  }
  if(!identical(arguments$edge_prior, "Stochastic-Block")) {
    stop(paste0(
      "edge_prior must be 'Stochastic-Block' (got '",
      as.character(arguments$edge_prior), "')."
    ))
  }

  return(list(
    posterior_num_blocks               = bgms_object$posterior_num_blocks,
    posterior_mean_allocations         = bgms_object$posterior_mean_allocations,
    posterior_mode_allocations         = bgms_object$posterior_mode_allocations,
    posterior_mean_coclustering_matrix = bgms_object$posterior_mean_coclustering_matrix
  ))
}

#' @inheritParams extract_sbm
#' @exportS3Method
#' @noRd
extract_sbm.bgmCompare = function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  if(!isTRUE(arguments$difference_selection)) {
    stop("To extract SBM summaries, run bgmCompare() with difference_selection = TRUE.")
  }
  if(!identical(arguments$difference_prior, "Stochastic-Block")) {
    stop(paste0(
      "difference_prior must be 'Stochastic-Block' (got '",
      as.character(arguments$difference_prior), "')."
    ))
  }

  return(list(
    posterior_num_blocks               = bgms_object$posterior_num_blocks,
    posterior_mean_allocations         = bgms_object$posterior_mean_allocations,
    posterior_mode_allocations         = bgms_object$posterior_mode_allocations,
    posterior_mean_coclustering_matrix = bgms_object$posterior_mean_coclustering_matrix
  ))
}


#' @inheritParams extract_posterior_inclusion_probabilities
#' @exportS3Method
#' @noRd
extract_posterior_inclusion_probabilities.bgmCompare = function(bgms_object,
                                                                estimator = c("rb", "raw")) {
  estimator_missing = missing(estimator)
  arguments = extract_arguments(bgms_object)

  if(!isTRUE(arguments$difference_selection)) {
    stop("To estimate posterior inclusion probabilities, run bgmCompare() with difference_selection = TRUE.")
  }

  var_names = arguments$data_columnnames
  # Handle legacy field name (no_variables -> num_variables in 0.1.6.0)
  num_variables = as.integer(arguments$num_variables %||% arguments$no_variables)

  raw = get_raw_samples(bgms_object)
  # Default to the Rao-Blackwellized estimate, but fall back to the raw
  # indicator average for fits without RB draws (bgms < 0.2.0.0), so reading a
  # legacy fit with the default call still works; an explicit estimator = "rb"
  # errors below instead.
  estimator = if(estimator_missing && is.null(raw$rb_inclusion)) {
    "raw"
  } else {
    match.arg(estimator)
  }
  if(estimator == "rb") {
    # Rao-Blackwellized average of the stored one-step draws J. Unselected
    # difference indicators have no draws and average to NA.
    if(is.null(raw$rb_inclusion)) {
      stop("No Rao-Blackwellized inclusion draws found in fit object; refit with bgms >= 0.2.0.0.")
    }
    mean_vals = apply(samples_to_array3d(raw$rb_inclusion), 3, mean, na.rm = TRUE)
  } else if(!is.null(raw$indicator)) {
    # Current format (0.1.6.0+): raw indicator average.
    mean_vals = apply(samples_to_array3d(raw$indicator), 3, mean)
  } else if(!is.null(bgms_object$pairwise_difference_indicator)) {
    # Deprecated format (0.1.4--0.1.5): $pairwise_difference_indicator at top level
    lifecycle::deprecate_warn(
      "0.1.6.0",
      I("The '$pairwise_difference_indicator' field is deprecated; please refit with bgms >= 0.1.6.0")
    )
    edge_means = colMeans(bgms_object$pairwise_difference_indicator)
    V = num_variables
    ind_mat = matrix(0, nrow = V, ncol = V)
    ind_mat[lower.tri(ind_mat)] = edge_means
    ind_mat = ind_mat + t(ind_mat)
    dimnames(ind_mat) = list(var_names, var_names)
    return(ind_mat)
  } else {
    stop("No indicator samples found in fit object.")
  }

  # Reconstruct the VxV matrix using the sampler's interleaved order:
  # (1,1),(1,2),...,(1,V),(2,2),...,(2,V),...,(V,V). Every cell is assigned,
  # so the NA init only survives where a mean is NA (unselected, rb estimator).
  V = num_variables
  stopifnot(length(mean_vals) == V * (V + 1L) / 2L)

  ind_mat = matrix(NA_real_,
    nrow = V, ncol = V,
    dimnames = list(var_names, var_names)
  )
  pos = 1L
  for(i in seq_len(V)) {
    ind_mat[i, i] = mean_vals[pos]
    pos = pos + 1L
    if(i < V) {
      for(j in (i + 1L):V) {
        val = mean_vals[pos]
        pos = pos + 1L
        ind_mat[i, j] = val
        ind_mat[j, i] = val
      }
    }
  }

  return(ind_mat)
}

#' @title Extract Indicator Prior Structure
#'
#' @description
#' Retrieves the prior specification used for inclusion indicators in a
#' model fitted with [bgm()] (edge indicators) or [bgmCompare()]
#' (difference indicators).
#'
#' @param bgms_object A fitted model object of class `bgms` (from [bgm()])
#'   or `bgmCompare` (from [bgmCompare()]).
#'
#' @return A named list describing the prior structure, including the prior
#'   type and any hyperparameters.
#'   \describe{
#'     \item{bgms}{Requires `edge_selection = TRUE`. Returns a list with the
#'       prior type (`"Bernoulli"`, `"Beta-Bernoulli"`, or
#'       `"Stochastic-Block"`) and associated hyperparameters.}
#'     \item{bgmCompare}{Requires `difference_selection = TRUE`. Returns the
#'       difference prior specification.}
#'   }
#'
#' @seealso [bgm()], [bgmCompare()], [extract_indicators()]
#' @family extractors
#' @export
extract_indicator_priors = function(bgms_object) {
  UseMethod("extract_indicator_priors")
}

#' @inheritParams extract_indicator_priors
#' @exportS3Method
#' @noRd
extract_indicator_priors.bgms = function(bgms_object) {
  arguments = extract_arguments(bgms_object)
  if(!isTRUE(arguments$edge_selection)) stop("No edge selection performed.")

  switch(arguments$edge_prior,
    "Bernoulli" = list(type = "Bernoulli", prior_inclusion_probability = arguments$inclusion_probability),
    "Beta-Bernoulli" = list(type = "Beta-Bernoulli", alpha = arguments$beta_bernoulli_alpha, beta = arguments$beta_bernoulli_beta),
    "Stochastic-Block" = list(
      type = "Stochastic-Block",
      beta_bernoulli_alpha = arguments$beta_bernoulli_alpha,
      beta_bernoulli_beta = arguments$beta_bernoulli_beta,
      dirichlet_alpha = arguments$dirichlet_alpha
    )
  )
}


#' @inheritParams extract_indicator_priors
#' @exportS3Method
#' @noRd
extract_indicator_priors.bgmCompare = function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  if(!isTRUE(arguments$difference_selection)) {
    stop("The model ran without selection, so there are no indicator priors specified.")
  }

  return(arguments$difference_prior)
}


#' @title Extract Pairwise Interaction Samples
#'
#' @description
#' Retrieves posterior samples of pairwise interaction parameters from a
#' model fitted with [bgm()] or [bgmCompare()].
#'
#' @param bgms_object A fitted model object of class `bgms` (from [bgm()])
#'   or `bgmCompare` (from [bgmCompare()]).
#'
#' @return A matrix with one row per post-warmup iteration and one column per
#'   edge, containing posterior samples of interaction strengths.
#'   \describe{
#'     \item{bgms}{Columns correspond to all unique variable pairs.}
#'     \item{bgmCompare}{Columns correspond to the baseline pairwise
#'       interaction parameters.}
#'   }
#'
#' @seealso [bgm()], [bgmCompare()], [extract_main_effects()]
#' @family extractors
#' @export
extract_pairwise_interactions = function(bgms_object) {
  UseMethod("extract_pairwise_interactions")
}

#' @inheritParams extract_pairwise_interactions
#' @exportS3Method
#' @noRd
extract_pairwise_interactions.bgms = function(bgms_object) {
  arguments = extract_arguments(bgms_object)
  # Handle legacy field name (no_variables -> num_variables in 0.1.6.0)
  num_vars = arguments$num_variables %||% arguments$no_variables
  var_names = arguments$data_columnnames

  # Current format (0.1.6.0+): raw samples
  raw = get_raw_samples(bgms_object)
  if(!is.null(raw)) {
    mats = raw$pairwise
    mat = do.call(rbind, mats)

    # Use stored parameter names when available (correct for all model types
    # including mixed MRF where block order differs from upper-triangle order)
    stored_names = raw$parameter_names$pairwise
    if(!is.null(stored_names)) {
      edge_names = stored_names
    } else {
      edge_names = character()
      for(i in 1:(num_vars - 1)) {
        for(j in (i + 1):num_vars) {
          edge_names = c(edge_names, paste0(var_names[i], "-", var_names[j]))
        }
      }
    }

    dimnames(mat) = list(paste0("iter", seq_len(nrow(mat))), edge_names)

    # GGM raw samples are on precision scale; convert to association scale
    if(isTRUE(arguments$is_continuous)) {
      mat = -0.5 * mat
    }

    return(mat)
  }

  # Deprecated format (0.1.4--0.1.5): $interactions
  if(!is.null(bgms_object$interactions)) {
    lifecycle::deprecate_warn(
      "0.1.6.0",
      I("The '$interactions' field is deprecated; please refit with bgms >= 0.1.6.0")
    )
    edge_means = colMeans(bgms_object$interactions)
    mat = matrix(0, nrow = num_vars, ncol = num_vars)
    mat[lower.tri(mat)] = edge_means
    mat = mat + t(mat)
    dimnames(mat) = list(var_names, var_names)
    return(mat)
  }

  # Defunct format (pre-0.1.4)
  lifecycle::deprecate_stop(
    "0.1.4.2",
    I("The '$pairwise_effects' field is defunct; please refit with bgms >= 0.1.6.0")
  )
}


#' @inheritParams extract_pairwise_interactions
#' @exportS3Method
#' @noRd
extract_pairwise_interactions.bgmCompare = function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  # Current format (0.1.6.0+)
  raw = get_raw_samples(bgms_object)
  if(!is.null(raw$pairwise)) {
    pairwise_samples = do.call(rbind, raw$pairwise)

    num_vars = bgms_object$arguments$num_variables
    num_pairs = num_vars * (num_vars - 1) / 2

    pairwise_samples = pairwise_samples[, 1:num_pairs]
    colnames(pairwise_samples) = raw$parameter_names$pairwise_baseline

    return(pairwise_samples)
  }

  # Deprecated format (0.1.4--0.1.5): $interactions at top level
  if(!is.null(bgms_object$interactions)) {
    lifecycle::deprecate_warn(
      "0.1.6.0",
      I("The '$interactions' field is deprecated; please refit with bgms >= 0.1.6.0")
    )
    return(bgms_object$interactions)
  }

  stop("No pairwise interaction samples found in fit object.")
}

#' @title Extract Main Effect Estimates
#'
#' @description
#' Retrieves main-effect parameters from a model fitted with [bgm()]
#' (posterior means) or [bgmCompare()] (posterior samples of baseline main
#' effects). For OMRF models these are category thresholds; for mixed MRF
#' models these include discrete thresholds and continuous means. GGM models
#' have no main effects and return `NULL`.
#'
#' @param bgms_object A fitted model object of class `bgms` (from [bgm()])
#'   or `bgmCompare` (from [bgmCompare()]).
#'
#' @return The structure depends on the model type:
#'   \describe{
#'     \item{GGM (bgms)}{`NULL` (invisibly). GGM models have no main effects;
#'       use [extract_precision()] to obtain the precision matrix.}
#'     \item{OMRF (bgms)}{A numeric matrix with one row per variable and one
#'       column per category threshold, containing posterior means. Columns
#'       beyond the number of categories for a variable are `NA`.}
#'     \item{Mixed MRF (bgms)}{A list with two elements:
#'       \describe{
#'         \item{discrete}{A numeric matrix (p rows x max_categories columns)
#'           of posterior mean thresholds for discrete variables.}
#'         \item{continuous}{A numeric matrix (q rows x 1 column) of
#'           posterior mean continuous variable means.}
#'       }}
#'     \item{bgmCompare}{A matrix with one row per post-warmup iteration,
#'       containing posterior samples of baseline main-effect parameters.}
#'   }
#'
#' @examples
#' \donttest{
#' fit = bgm(x = Wenchuan[, 1:3])
#' extract_main_effects(fit)
#' }
#'
#' @seealso [bgm()], [bgmCompare()], [extract_pairwise_interactions()],
#'   [extract_category_thresholds()]
#' @family extractors
#' @export
extract_main_effects = function(bgms_object) {
  UseMethod("extract_main_effects")
}

#' @inheritParams extract_main_effects
#' @exportS3Method
#' @noRd
extract_main_effects.bgms = function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  # GGM: no main effects; use extract_precision() for the precision matrix
  if(isTRUE(arguments$is_continuous)) {
    return(invisible(NULL))
  }

  # Mixed MRF: return pre-built list from posterior_mean_main
  if(isTRUE(arguments$is_mixed)) {
    return(get_posterior_mean(bgms_object, "main"))
  }

  # OMRF: return pre-built threshold matrix
  pm_main = get_posterior_mean(bgms_object, "main")
  if(!is.null(pm_main)) {
    return(pm_main)
  }

  # Deprecated format (0.1.4--0.1.5): $thresholds
  if(!is.null(bgms_object$thresholds)) {
    lifecycle::deprecate_warn(
      "0.1.6.0",
      I("The '$thresholds' field is deprecated; please refit with bgms >= 0.1.6.0")
    )
    var_names = arguments$data_columnnames
    means = colMeans(bgms_object$thresholds)
    mat = matrix(means, nrow = length(means), ncol = 1)
    rownames(mat) = var_names
    return(mat)
  }

  # Defunct format (pre-0.1.4)
  lifecycle::deprecate_stop(
    "0.1.4.2",
    I("The '$main_effects' field is defunct; please refit with bgms >= 0.1.6.0")
  )
}

#' @inheritParams extract_main_effects
#' @exportS3Method
#' @noRd
extract_main_effects.bgmCompare = function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  # Current format (0.1.6.0+)
  raw = get_raw_samples(bgms_object)
  if(!is.null(raw$main)) {
    main_samples = do.call(rbind, raw$main)

    num_main = length(raw$parameter_names$main_baseline)

    main_samples = main_samples[, 1:num_main]
    colnames(main_samples) = raw$parameter_names$main_baseline

    return(main_samples)
  }

  # Deprecated format (0.1.4--0.1.5): $thresholds or $thresholds_gr1/$thresholds_gr2 at top level
  if(!is.null(bgms_object$thresholds)) {
    lifecycle::deprecate_warn(
      "0.1.6.0",
      I("The '$thresholds' field is deprecated; please refit with bgms >= 0.1.6.0")
    )
    return(bgms_object$thresholds)
  }

  # Alternative deprecated format (0.1.4.1+): $thresholds_gr1, $thresholds_gr2
  if(!is.null(bgms_object$thresholds_gr1)) {
    lifecycle::deprecate_warn(
      "0.1.6.0",
      I("The '$thresholds_gr*' fields are deprecated; please refit with bgms >= 0.1.6.0")
    )
    return(cbind(bgms_object$thresholds_gr1, bgms_object$thresholds_gr2))
  }

  stop("No main effect samples found in fit object.")
}


#' @title Extract Category Threshold Estimates
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `extract_category_thresholds()` was renamed to [extract_main_effects()] to
#' reflect that main effects include continuous means and precision diagonal
#' (mixed MRF), not only category thresholds.
#'
#' @param bgms_object A fitted model object of class `bgms` (from [bgm()])
#'   or `bgmCompare` (from [bgmCompare()]).
#'
#' @return See [extract_main_effects()] for details.
#'
#' @seealso [extract_main_effects()]
#' @family extractors
#' @export
extract_category_thresholds = function(bgms_object) {
  lifecycle::deprecate_warn(
    "0.2.0",
    "extract_category_thresholds()",
    "extract_main_effects()"
  )
  extract_main_effects(bgms_object)
}

#' @title Extract Group-Specific Parameters
#'
#' @description
#' Computes group-specific parameter estimates by combining baseline
#' parameters and group differences from a model fitted with [bgmCompare()].
#'
#' @param bgms_object A fitted model object of class `bgmCompare`
#'   (from [bgmCompare()]).
#'
#' @return A list with elements `main_effects_groups` (main effects per
#'   group) and `pairwise_effects_groups` (pairwise effects per group).
#'
#' @seealso [bgmCompare()], [extract_pairwise_interactions()],
#'   [extract_main_effects()]
#' @family extractors
#' @export
extract_group_params = function(bgms_object) {
  UseMethod("extract_group_params")
}

#' @inheritParams extract_group_params
#' @exportS3Method
#' @noRd
extract_group_params.bgmCompare = function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  # Current format (0.1.6.0+)
  if(!is.null(get_raw_samples(bgms_object)$main)) {
    return(.extract_group_params_current(bgms_object, arguments))
  }

  # Deprecated format (0.1.4--0.1.5): separate fields for baseline and differences
  if(!is.null(bgms_object$interactions) && !is.null(bgms_object$pairwise_difference)) {
    lifecycle::deprecate_warn(
      "0.1.6.0",
      I("The legacy bgmCompare format is deprecated; please refit with bgms >= 0.1.6.0")
    )
    return(.extract_group_params_legacy(bgms_object, arguments))
  }

  stop("No group parameter samples found in fit object.")
}

# Helper for current format (0.1.6+)
# ------------------------------------------------------------------------------
# .compute_group_param_matrices()
# ------------------------------------------------------------------------------
# Shared bgmCompare group-parameter reconstruction. From the posterior-mean
# main and pairwise samples, builds the [param x (baseline, diff...)] matrices
# and the projection-expanded [param x group] effect matrices. Used by both
# extract_group_params() (group effects only) and coef.bgmCompare() (which also
# returns the raw baseline+diff matrices alongside the indicators).
#
# @param arguments  Fit arguments (data_columnnames, num_categories,
#   is_ordinal_variable, num_groups, num_variables, projection).
# @param raw        get_raw_samples() output (uses $main, $pairwise).
#
# Returns: list(main_mat, pairwise_mat, main_effects_groups,
#   pairwise_effects_groups).
# ------------------------------------------------------------------------------
.compute_group_param_matrices = function(arguments, raw) {
  var_names = arguments$data_columnnames
  num_categories = as.integer(arguments$num_categories)
  is_ordinal = as.logical(arguments$is_ordinal_variable)
  num_groups = as.integer(arguments$num_groups)
  num_variables = as.integer(arguments$num_variables)
  projection = arguments$projection # [num_groups x (num_groups-1)]

  # ---- main effects ----
  array3d_main = samples_to_array3d(raw$main)
  stopifnot(!is.null(array3d_main))
  mean_main = apply(array3d_main, 3, mean)

  stopifnot(length(mean_main) %% num_groups == 0L)
  num_main = as.integer(length(mean_main) / num_groups)

  main_mat = matrix(mean_main, nrow = num_main, ncol = num_groups, byrow = FALSE)

  # row names in sampler row order
  rownames(main_mat) = unlist(lapply(seq_len(num_variables), function(v) {
    if(is_ordinal[v]) {
      paste0(var_names[v], "(c", seq_len(num_categories[v]), ")")
    } else {
      c(
        paste0(var_names[v], "(linear)"),
        paste0(var_names[v], "(quadratic)")
      )
    }
  }))
  colnames(main_mat) = c("baseline", paste0("diff", seq_len(num_groups - 1L)))

  # group-specific main effects: baseline + P %*% diffs
  main_effects_groups = matrix(NA_real_, nrow = num_main, ncol = num_groups)
  for(r in seq_len(num_main)) {
    baseline = main_mat[r, 1]
    diffs = main_mat[r, -1, drop = TRUE]
    main_effects_groups[r, ] = baseline + as.vector(projection %*% diffs)
  }
  rownames(main_effects_groups) = rownames(main_mat)
  colnames(main_effects_groups) = paste0("group", seq_len(num_groups))

  # ---- pairwise effects ----
  array3d_pair = samples_to_array3d(raw$pairwise)
  stopifnot(!is.null(array3d_pair))
  mean_pair = apply(array3d_pair, 3, mean)

  stopifnot(length(mean_pair) %% num_groups == 0L)
  num_pair = as.integer(length(mean_pair) / num_groups)

  pairwise_mat = matrix(mean_pair, nrow = num_pair, ncol = num_groups, byrow = FALSE)

  # row names in sampler row order (upper-tri i<j)
  pair_names = character()
  if(num_variables >= 2L) {
    for(i in 1L:(num_variables - 1L)) {
      for(j in (i + 1L):num_variables) {
        pair_names = c(pair_names, paste0(var_names[i], "-", var_names[j]))
      }
    }
  }
  rownames(pairwise_mat) = pair_names
  colnames(pairwise_mat) = c("baseline", paste0("diff", seq_len(num_groups - 1L)))

  # group-specific pairwise effects
  pairwise_effects_groups = matrix(NA_real_, nrow = num_pair, ncol = num_groups)
  for(r in seq_len(num_pair)) {
    baseline = pairwise_mat[r, 1]
    diffs = pairwise_mat[r, -1, drop = TRUE]
    pairwise_effects_groups[r, ] = baseline + as.vector(projection %*% diffs)
  }
  rownames(pairwise_effects_groups) = rownames(pairwise_mat)
  colnames(pairwise_effects_groups) = paste0("group", seq_len(num_groups))

  list(
    main_mat = main_mat,
    pairwise_mat = pairwise_mat,
    main_effects_groups = main_effects_groups,
    pairwise_effects_groups = pairwise_effects_groups
  )
}

.extract_group_params_current = function(bgms_object, arguments) {
  gp = .compute_group_param_matrices(arguments, get_raw_samples(bgms_object))
  list(
    main_effects_groups = gp$main_effects_groups,
    pairwise_effects_groups = gp$pairwise_effects_groups
  )
}

# Helper for legacy format (0.1.4--0.1.5)
# v0.1.4.x only supported 2 groups with parameterization:
#   group1 = baseline + diff, group2 = baseline - diff
.extract_group_params_legacy = function(bgms_object, arguments) {
  var_names = arguments$data_columnnames
  # Handle legacy field name (no_variables -> num_variables in 0.1.6.0)
  num_variables = as.integer(arguments$num_variables %||% arguments$no_variables)

  # v0.1.4 format: baseline interactions and differences are separate
  # $interactions: [iter x n_pairs] baseline pairwise effects
  # $pairwise_difference: [iter x n_pairs] pairwise differences
  # $thresholds or $thresholds_gr1/$thresholds_gr2: main effects
  # $main_difference: [iter x n_vars] main differences

  # Compute posterior means
  mean_interactions = colMeans(bgms_object$interactions)
  mean_pairwise_diff = colMeans(bgms_object$pairwise_difference)

  # Get thresholds (handles both v0.1.4 and v0.1.4.1+ formats)
  if(!is.null(bgms_object$thresholds)) {
    mean_thresholds = colMeans(bgms_object$thresholds)
  } else if(!is.null(bgms_object$thresholds_gr1)) {
    # v0.1.4.1+ stored group-specific thresholds directly
    mean_thresholds_gr1 = colMeans(bgms_object$thresholds_gr1)
    mean_thresholds_gr2 = colMeans(bgms_object$thresholds_gr2)
    # Return directly since we have group-specific values
    main_effects_groups = cbind(mean_thresholds_gr1, mean_thresholds_gr2)
    colnames(main_effects_groups) = c("group1", "group2")
    rownames(main_effects_groups) = var_names

    pairwise_effects_groups = cbind(
      mean_interactions + mean_pairwise_diff,
      mean_interactions - mean_pairwise_diff
    )
    colnames(pairwise_effects_groups) = c("group1", "group2")

    # Row names for pairs
    pair_names = character()
    if(num_variables >= 2L) {
      for(i in 1L:(num_variables - 1L)) {
        for(j in (i + 1L):num_variables) {
          pair_names = c(pair_names, paste0(var_names[i], "-", var_names[j]))
        }
      }
    }
    rownames(pairwise_effects_groups) = pair_names

    return(list(
      main_effects_groups = main_effects_groups,
      pairwise_effects_groups = pairwise_effects_groups
    ))
  } else {
    stop("No threshold samples found in legacy fit object.")
  }

  mean_main_diff = colMeans(bgms_object$main_difference)

  # v0.1.4 parameterization: group1 = baseline + diff, group2 = baseline - diff
  main_effects_groups = cbind(
    mean_thresholds + mean_main_diff,
    mean_thresholds - mean_main_diff
  )
  colnames(main_effects_groups) = c("group1", "group2")
  rownames(main_effects_groups) = var_names

  pairwise_effects_groups = cbind(
    mean_interactions + mean_pairwise_diff,
    mean_interactions - mean_pairwise_diff
  )
  colnames(pairwise_effects_groups) = c("group1", "group2")

  # Row names for pairs
  pair_names = character()
  if(num_variables >= 2L) {
    for(i in 1L:(num_variables - 1L)) {
      for(j in (i + 1L):num_variables) {
        pair_names = c(pair_names, paste0(var_names[i], "-", var_names[j]))
      }
    }
  }
  rownames(pairwise_effects_groups) = pair_names

  return(list(
    main_effects_groups = main_effects_groups,
    pairwise_effects_groups = pairwise_effects_groups
  ))
}

#' @title Deprecated: Use extract_indicators instead
#' @param bgms_object A bgms or bgmCompare object.
#' @keywords internal
#' @export
extract_edge_indicators = function(bgms_object) {
  lifecycle::deprecate_warn("0.1.4.2", "extract_edge_indicators()", "extract_indicators()")
  extract_indicators(bgms_object)
}

#' @title Deprecated: Use extract_main_effects instead
#' @param bgms_object A bgms or bgmCompare object.
#' @keywords internal
#' @export
extract_pairwise_thresholds = function(bgms_object) {
  lifecycle::deprecate_warn("0.1.4.2", "extract_pairwise_thresholds()", "extract_main_effects()")
  extract_main_effects(bgms_object)
}


# ------------------------------------------------------------------------------
# extract_rhat() - R-hat Convergence Diagnostics
# ------------------------------------------------------------------------------

#' @title Extract R-hat Convergence Diagnostics
#'
#' @description
#' Retrieves R-hat convergence diagnostics for all parameters from a
#' model fitted with [bgm()] or [bgmCompare()].
#'
#' @param bgms_object A fitted model object of class `bgms` (from [bgm()])
#'   or `bgmCompare` (from [bgmCompare()]).
#'
#' @return A named list with R-hat values for each parameter type present in
#'   the model (e.g., `main`, `pairwise`, `indicator`).
#'
#' @details
#' The `indicator` element is the split-R-hat of the Rao-Blackwellized
#' inclusion draws, matching the default of [extract_ess()]; the indicator
#' chain's transition-based effective sample size remains available as
#' `extract_ess(fit, estimator = "mixt")` (deprecated).
#'
#' One failure mode these diagnostics may miss: chains can agree that an edge
#' is included while disagreeing about the magnitude of its weight, for
#' example when the weight's posterior has more than one mode and different
#' chains settle in different ones. The indicator diagnostics cannot see this
#' by construction, and the pooled-draw `Rhat` responds only weakly when
#' inclusion is intermittent. When it matters, compare each edge's weight
#' across chains using only the draws in which the edge is included; a
#' dedicated diagnostic for this is planned for a future release.
#'
#' @seealso [bgm()], [bgmCompare()], [extract_ess()]
#' @family extractors
#' @export
extract_rhat = function(bgms_object) {
  UseMethod("extract_rhat")
}

#' @inheritParams extract_rhat
#' @exportS3Method
#' @noRd
extract_rhat.bgms = function(bgms_object) {
  ensure_summaries(bgms_object)
  result = list()

  # Main effect Rhat
  if(!is.null(bgms_object$posterior_summary_main)) {
    result$main = bgms_object$posterior_summary_main$Rhat
    names(result$main) = rownames(bgms_object$posterior_summary_main)
  }

  # Precision diagonal (quadratic) Rhat
  if(!is.null(bgms_object$posterior_summary_quadratic)) {
    result$quadratic = bgms_object$posterior_summary_quadratic$Rhat
    names(result$quadratic) = rownames(bgms_object$posterior_summary_quadratic)
  }

  # Pairwise interaction Rhat
  if(!is.null(bgms_object$posterior_summary_pairwise)) {
    result$pairwise = bgms_object$posterior_summary_pairwise$Rhat
    names(result$pairwise) = rownames(bgms_object$posterior_summary_pairwise)
  }

  # Indicator Rhat (if edge selection was used)
  if(!is.null(bgms_object$posterior_summary_indicator)) {
    result$indicator = bgms_object$posterior_summary_indicator$Rhat
    names(result$indicator) = rownames(bgms_object$posterior_summary_indicator)
  }

  if(length(result) == 0) {
    stop("No posterior summary information found in this object.")
  }

  return(result)
}

#' @inheritParams extract_rhat
#' @exportS3Method
#' @noRd
extract_rhat.bgmCompare = function(bgms_object) {
  ensure_summaries(bgms_object)
  result = list()

  # Main baseline Rhat
  if(!is.null(bgms_object$posterior_summary_main_baseline)) {
    result$main_baseline = bgms_object$posterior_summary_main_baseline$Rhat
    names(result$main_baseline) = rownames(bgms_object$posterior_summary_main_baseline)
  }

  # Main differences Rhat
  if(!is.null(bgms_object$posterior_summary_main_differences)) {
    result$main_differences = bgms_object$posterior_summary_main_differences$Rhat
    names(result$main_differences) = rownames(bgms_object$posterior_summary_main_differences)
  }

  # Pairwise baseline Rhat
  if(!is.null(bgms_object$posterior_summary_pairwise_baseline)) {
    result$pairwise_baseline = bgms_object$posterior_summary_pairwise_baseline$Rhat
    names(result$pairwise_baseline) = rownames(bgms_object$posterior_summary_pairwise_baseline)
  }

  # Pairwise differences Rhat
  if(!is.null(bgms_object$posterior_summary_pairwise_differences)) {
    result$pairwise_differences = bgms_object$posterior_summary_pairwise_differences$Rhat
    names(result$pairwise_differences) = rownames(bgms_object$posterior_summary_pairwise_differences)
  }

  # Indicator Rhat (if difference selection was used)
  if(!is.null(bgms_object$posterior_summary_indicator)) {
    result$indicator = bgms_object$posterior_summary_indicator$Rhat
    names(result$indicator) = rownames(bgms_object$posterior_summary_indicator)
  }

  if(length(result) == 0) {
    stop("No posterior summary information found in this object.")
  }

  return(result)
}


# ------------------------------------------------------------------------------
# extract_ess() - Effective Sample Size
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# indicator_ess_column
# ------------------------------------------------------------------------------
# Resolve the indicator ESS that extract_ess() returns. The default "rb" is the
# continuous ESS of the Rao-Blackwellized inclusion draws, read from the summary
# table, so the NA masking summarize_rb_inclusion() applies carries over
# unchanged. The deprecated "mixt" is the indicator chain's transition-based
# ESS, no longer a summary column; it is recomputed from the raw indicator draws
# for as long as it is accepted. The deprecation warning is raised in the
# generic, where the caller is the user rather than bgms itself.
#
# Fits without Rao-Blackwellized draws (bgms < 0.2.0.0) carry no n_eff column,
# so a default call falls back to the transition ESS while an explicit
# estimator = "rb" errors, as in
# extract_posterior_inclusion_probabilities().
#
# @param fit                The fit object (for the raw indicator draws).
# @param summary_indicator  The fit's posterior_summary_indicator table.
# @param estimator          Character: "rb" or "mixt" (or the unevaluated
#   default vector).
# @param estimator_missing  Logical: TRUE if the caller did not supply one.
#
# Returns: numeric vector of ESS values, one per indicator.
# ------------------------------------------------------------------------------
indicator_ess_column = function(fit, summary_indicator, estimator, estimator_missing) {
  # names(), not $, because $ partial-matches n_eff to a legacy table's
  # n_eff_mixt.
  has_rb = "n_eff" %in% names(summary_indicator)
  estimator = if(estimator_missing && !has_rb) {
    "mixt"
  } else {
    match.arg(estimator, choices = c("rb", "mixt"))
  }
  if(estimator == "rb") {
    if(!has_rb) {
      stop(
        "This fit carries no Rao-Blackwellized inclusion draws, so the ",
        "Rao-Blackwellized effective sample size is unavailable. Refit with ",
        "bgms >= 0.2.0.0."
      )
    }
    return(summary_indicator$n_eff)
  }

  if("n_eff_mixt" %in% names(summary_indicator)) {
    return(summary_indicator$n_eff_mixt)
  }
  indicator_transition_ess(fit)
}


# ------------------------------------------------------------------------------
# indicator_transition_ess
# ------------------------------------------------------------------------------
# Transition-based ESS of the binary indicator chains, recomputed from the raw
# draws. Retained only to keep the deprecated extract_ess(estimator = "mixt")
# working; the quantity is no longer a summary column.
#
# @param fit  A bgms or bgmCompare object.
#
# Returns: numeric vector of transition ESS values, one per indicator.
# ------------------------------------------------------------------------------
indicator_transition_ess = function(fit) {
  cache = get_fit_cache(fit)
  raw = if(is.null(cache)) NULL else cache$raw
  if(is.null(raw) || is.null(raw[[1]][["indicator_samples"]])) {
    stop(
      "This fit carries no raw indicator draws, so the transition-based ",
      "effective sample size cannot be recomputed."
    )
  }
  ind_stats = .compute_indicator_ess_cpp(combine_chains(raw, "indicator_samples"))
  unname(ind_stats[, "n_eff_mixt"])
}


#' @title Extract Effective Sample Size
#'
#' @description
#' Retrieves effective sample size estimates for all parameters from a
#' model fitted with [bgm()] or [bgmCompare()].
#'
#' @param bgms_object A fitted model object of class `bgms` (from [bgm()])
#'   or `bgmCompare` (from [bgmCompare()]).
#' @param estimator Character; which effective sample size to return for the
#'   edge (or difference) indicators. `"rb"` (default) returns the ESS of the
#'   Rao-Blackwellized inclusion draws. `"mixt"`
#'   `r lifecycle::badge("deprecated")` returns the indicator chain's
#'   transition-based ESS. Ignored for all other parameter types.
#'
#' @return A named list with ESS values for each parameter type present in
#'   the model (e.g., `main`, `pairwise`, `indicator`).
#'
#' @details
#' The `indicator` element is the effective sample size of the inclusion
#' inference, and also the `n_eff` column of the fit summary's inclusion table
#' (`summary(fit)$indicator`): the continuous ESS of the Rao-Blackwellized
#' inclusion draws, and therefore the ESS of the inclusion probability the fit
#' reports (see [extract_posterior_inclusion_probabilities()], which is
#' Rao-Blackwellized by default). It is `NA` for indicators whose
#' Rao-Blackwellized draws are constant to double precision, where the
#' inclusion probability is at its numerical bound.
#'
#' `estimator = "mixt"`, the transition-based ESS of the binary indicator chain,
#' is deprecated and no longer a summary column; a call recomputes it from the
#' raw indicator draws and warns. It converts flip counts to an effective sample
#' size through a two-state first-order Markov model that a substantial share of
#' edge chains violate, and the precision of the inclusion probability and its
#' Bayes factor is carried by the Rao-Blackwellized ESS and Monte Carlo standard
#' error. The directional flip counts themselves (`n0->1`, `n1->0`) remain in
#' the inclusion table.
#'
#' Fits made with bgms < 0.2.0.0 carry no Rao-Blackwellized draws; a default
#' call falls back to the transition ESS without warning, while an explicit
#' `estimator = "rb"` errors.
#'
#' @seealso [bgm()], [bgmCompare()], [extract_rhat()] for the
#'   Rao-Blackwellized indicator R-hat,
#'   [extract_posterior_inclusion_probabilities()]
#' @family extractors
#' @export
extract_ess = function(bgms_object, estimator = c("rb", "mixt")) {
  if(!missing(estimator) && identical(estimator, "mixt")) {
    lifecycle::deprecate_warn(
      "0.2.0.0", "extract_ess(estimator = 'no longer supports \"mixt\"')",
      details = paste(
        "The transition-based effective sample size no longer bears on any",
        "inclusion verdict: the two-state model behind it is rejected on a",
        "substantial share of edge chains, and the precision of the inclusion",
        "probability and its Bayes factor is carried by the Rao-Blackwellized",
        "effective sample size and Monte Carlo standard error. Use the default",
        "estimator = \"rb\"; the directional flip counts remain in",
        "summary(fit)$indicator."
      )
    )
  }
  UseMethod("extract_ess")
}

#' @inheritParams extract_ess
#' @exportS3Method
#' @noRd
extract_ess.bgms = function(bgms_object, estimator = c("rb", "mixt")) {
  estimator_missing = missing(estimator)
  ensure_summaries(bgms_object)
  result = list()

  # Main effect ESS
  if(!is.null(bgms_object$posterior_summary_main)) {
    result$main = bgms_object$posterior_summary_main$n_eff
    names(result$main) = rownames(bgms_object$posterior_summary_main)
  }

  # Precision diagonal (quadratic) ESS
  if(!is.null(bgms_object$posterior_summary_quadratic)) {
    result$quadratic = bgms_object$posterior_summary_quadratic$n_eff
    names(result$quadratic) = rownames(bgms_object$posterior_summary_quadratic)
  }

  # Pairwise interaction ESS
  if(!is.null(bgms_object$posterior_summary_pairwise)) {
    result$pairwise = bgms_object$posterior_summary_pairwise$n_eff
    names(result$pairwise) = rownames(bgms_object$posterior_summary_pairwise)
  }

  # Indicator ESS (if edge selection was used)
  if(!is.null(bgms_object$posterior_summary_indicator)) {
    result$indicator = indicator_ess_column(
      bgms_object, bgms_object$posterior_summary_indicator,
      estimator, estimator_missing
    )
    names(result$indicator) = rownames(bgms_object$posterior_summary_indicator)
  }

  if(length(result) == 0) {
    stop("No posterior summary information found in this object.")
  }

  return(result)
}

#' @inheritParams extract_ess
#' @exportS3Method
#' @noRd
extract_ess.bgmCompare = function(bgms_object, estimator = c("rb", "mixt")) {
  estimator_missing = missing(estimator)
  ensure_summaries(bgms_object)
  result = list()

  # Main baseline ESS
  if(!is.null(bgms_object$posterior_summary_main_baseline)) {
    result$main_baseline = bgms_object$posterior_summary_main_baseline$n_eff
    names(result$main_baseline) = rownames(bgms_object$posterior_summary_main_baseline)
  }

  # Main differences ESS
  if(!is.null(bgms_object$posterior_summary_main_differences)) {
    result$main_differences = bgms_object$posterior_summary_main_differences$n_eff
    names(result$main_differences) = rownames(bgms_object$posterior_summary_main_differences)
  }

  # Pairwise baseline ESS
  if(!is.null(bgms_object$posterior_summary_pairwise_baseline)) {
    result$pairwise_baseline = bgms_object$posterior_summary_pairwise_baseline$n_eff
    names(result$pairwise_baseline) = rownames(bgms_object$posterior_summary_pairwise_baseline)
  }

  # Pairwise differences ESS
  if(!is.null(bgms_object$posterior_summary_pairwise_differences)) {
    result$pairwise_differences = bgms_object$posterior_summary_pairwise_differences$n_eff
    names(result$pairwise_differences) = rownames(bgms_object$posterior_summary_pairwise_differences)
  }

  # Indicator ESS (if difference selection was used)
  if(!is.null(bgms_object$posterior_summary_indicator)) {
    result$indicator = indicator_ess_column(
      bgms_object, bgms_object$posterior_summary_indicator,
      estimator, estimator_missing
    )
    names(result$indicator) = rownames(bgms_object$posterior_summary_indicator)
  }

  if(length(result) == 0) {
    stop("No posterior summary information found in this object.")
  }

  return(result)
}


#' @title Extract Posterior Mean Precision Matrix
#'
#' @description
#' Retrieves the posterior mean precision matrix from a model fitted with
#' [bgm()]. For GGM models this is the full precision matrix. For
#' mixed MRF models this is the precision matrix of the continuous
#' (Gaussian) block. OMRF models have no precision matrix and return `NULL`.
#'
#' For mixed MRF models the precision matrix is reconstructed from the
#' internal association-scale parameterization.
#'
#' @param bgms_object A fitted model object of class `bgms` (from [bgm()]).
#'
#' @return A named numeric matrix containing the posterior mean precision
#'   matrix, or `NULL` for OMRF models.
#'   \describe{
#'     \item{GGM}{A symmetric matrix with one row and column per variable.}
#'     \item{Mixed MRF}{A symmetric matrix with one row and column per
#'       continuous variable.}
#'     \item{OMRF}{`NULL` (invisibly).}
#'   }
#'
#' @examples
#' \donttest{
#' fit = bgm(
#'   x = Wenchuan[, 1:3],
#'   variable_type = rep("continuous", 3)
#' )
#' extract_precision(fit)
#' }
#'
#' @seealso [bgm()], [coef.bgms()], [extract_partial_correlations()]
#' @family extractors
#' @export
extract_precision = function(bgms_object) {
  UseMethod("extract_precision")
}

#' @inheritParams extract_precision
#' @exportS3Method
#' @noRd
extract_precision.bgms = function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  # OMRF: no precision matrix
  if(!isTRUE(arguments$is_continuous) && !isTRUE(arguments$is_mixed)) {
    return(invisible(NULL))
  }

  rv = get_posterior_mean(bgms_object, "residual_variance")
  associations = get_posterior_mean(bgms_object, "pairwise")

  if(isTRUE(arguments$is_mixed)) {
    # Mixed MRF: extract the q x q continuous block, convert to precision
    cont_idx = arguments$continuous_indices
    cont_names = arguments$data_columnnames_continuous
    cont_block = associations[cont_idx, cont_idx]
    precision = -2 * cont_block
    diag(precision) = 1 / rv
    dimnames(precision) = list(cont_names, cont_names)
    return(precision)
  }

  # GGM: associations are stored at half precision scale; convert to precision
  precision = -2 * associations
  diag(precision) = 1 / rv
  return(precision)
}


#' @title Extract Posterior Mean Partial Correlations
#'
#' @description
#' Computes the posterior mean partial correlation matrix from a model fitted
#' with [bgm()]. For GGM models this is the full matrix. For mixed
#' MRF models this is the matrix for the continuous block. OMRF models
#' have no partial correlations and return `NULL`.
#'
#' Partial correlations are computed from the precision matrix as
#' \eqn{\rho_{ij} = -\Theta_{ij} / \sqrt{\Theta_{ii} \Theta_{jj}}}{rho_ij = -Theta_ij / sqrt(Theta_ii * Theta_jj)}.
#'
#' @param bgms_object A fitted model object of class `bgms` (from [bgm()]).
#'
#' @return A named numeric matrix containing posterior mean partial
#'   correlations, or `NULL` for OMRF models.
#'   \describe{
#'     \item{GGM}{A symmetric matrix with ones on the diagonal and one
#'       row and column per variable.}
#'     \item{Mixed MRF}{A symmetric matrix with ones on the diagonal and
#'       one row and column per continuous variable.}
#'     \item{OMRF}{`NULL` (invisibly).}
#'   }
#'
#' @examples
#' \donttest{
#' fit = bgm(
#'   x = Wenchuan[, 1:3],
#'   variable_type = rep("continuous", 3)
#' )
#' extract_partial_correlations(fit)
#' }
#'
#' @seealso [bgm()], [extract_precision()]
#' @family extractors
#' @export
extract_partial_correlations = function(bgms_object) {
  UseMethod("extract_partial_correlations")
}

#' @inheritParams extract_partial_correlations
#' @exportS3Method
#' @noRd
extract_partial_correlations.bgms = function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  # OMRF: no partial correlations
  if(!isTRUE(arguments$is_continuous) && !isTRUE(arguments$is_mixed)) {
    return(invisible(NULL))
  }

  # Derive from precision: rho_ij = -Theta_ij / sqrt(Theta_ii * Theta_jj)
  precision = extract_precision(bgms_object)
  d = sqrt(diag(precision))
  partial_corr = -precision / outer(d, d)
  diag(partial_corr) = 1
  return(partial_corr)
}


#' @title Extract Posterior Mean Log-Odds (Pairwise Interactions)
#'
#' @description
#' Retrieves the posterior mean pairwise interaction matrix for discrete
#' variables from a model fitted with [bgm()]. These are the log-odds
#' parameters of the discrete (Markov random field) block. GGM models have
#' no discrete variables and return `NULL`.
#'
#' @param bgms_object A fitted model object of class `bgms` (from [bgm()]).
#'
#' @return A named numeric matrix of posterior mean log-odds interactions, or
#'   `NULL` for GGM models.
#'   \describe{
#'     \item{OMRF}{A symmetric matrix with zero diagonal and one row and
#'       column per variable.}
#'     \item{Mixed MRF}{A symmetric matrix with zero diagonal and one row
#'       and column per discrete variable.}
#'     \item{GGM}{`NULL` (invisibly).}
#'   }
#'
#' @examples
#' \donttest{
#' fit = bgm(x = Wenchuan[, 1:3])
#' extract_log_odds(fit)
#' }
#'
#' @seealso [bgm()], [extract_pairwise_interactions()], [extract_precision()]
#' @family extractors
#' @export
extract_log_odds = function(bgms_object) {
  UseMethod("extract_log_odds")
}

#' @inheritParams extract_log_odds
#' @exportS3Method
#' @noRd
extract_log_odds.bgms = function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  # GGM: no discrete variables
  if(isTRUE(arguments$is_continuous)) {
    return(invisible(NULL))
  }

  associations = get_posterior_mean(bgms_object, "pairwise")

  if(isTRUE(arguments$is_mixed)) {
    # Mixed MRF: extract the p x p discrete block, convert to log-odds
    disc_idx = arguments$discrete_indices
    disc_names = arguments$data_columnnames_discrete
    log_odds = 2 * associations[disc_idx, disc_idx]
    dimnames(log_odds) = list(disc_names, disc_names)
    return(log_odds)
  }

  # OMRF: log adjacent-category odds ratio = 2 * association
  return(2 * associations)
}

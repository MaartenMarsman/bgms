# ==============================================================================
# Edge verdicts and the Monte Carlo fragility flag
# ==============================================================================
#
# The three-way reading of an inclusion Bayes factor (presence / undecided /
# absence at a threshold t and its reciprocal) plus the flag that says whether
# the reading would survive a rerun. Both standard errors behind the flag are
# on the log inclusion Bayes factor scale, which equals the logit inclusion
# probability scale shifted by the (constant) prior inclusion odds.
# ==============================================================================


# ------------------------------------------------------------------
# two_state_se_logit
# ------------------------------------------------------------------
# Jeffreys-smoothed two-state summary of binary indicator chains: the standard
# error of the logit inclusion probability implied by a first-order two-state
# Markov model of the indicator. Transitions are counted within chains, never
# across the join. The smoothing is what makes the standard error defined at
# zero flips, where the unsmoothed transition rates are 0/0.
#
# @param draws  List of chains, each an niter x nparam 0/1 matrix.
#
# Returns: numeric vector of logit-scale standard errors, one per indicator.
# ------------------------------------------------------------------
two_state_se_logit = function(draws) {
  nparam = ncol(draws[[1]])
  niter = nrow(draws[[1]])
  n_total = niter * length(draws)

  n0 = n1 = n01 = n10 = numeric(nparam)
  for(chain in draws) {
    from = chain[-niter, , drop = FALSE]
    to = chain[-1L, , drop = FALSE]
    n0 = n0 + colSums(from == 0)
    n1 = n1 + colSums(from == 1)
    n01 = n01 + colSums(from == 0 & to == 1)
    n10 = n10 + colSums(from == 1 & to == 0)
  }

  a = (n01 + 0.5) / (n0 + 1)
  b = (n10 + 0.5) / (n1 + 1)
  p = a / (a + b)
  ess = n_total * (a + b) / (2 - a - b)
  sqrt(1 / (ess * p * (1 - p)))
}


# ------------------------------------------------------------------
# rb_se_logit
# ------------------------------------------------------------------
# Standard error of the logit inclusion probability from the Rao-Blackwellized
# draws, by the delta method: d logit(p) / dp = 1 / (p (1 - p)). The input MCSE
# is the summary table's, so it is NA where the RB draws are constant to double
# precision, and the standard error is NA there too.
#
# @param mcse  MCSE of the Rao-Blackwellized inclusion probability.
# @param pip   The Rao-Blackwellized inclusion probability.
#
# Returns: numeric vector of logit-scale standard errors, one per indicator.
# ------------------------------------------------------------------
rb_se_logit = function(mcse, pip) {
  jacobian = 1 / (pip * (1 - pip))
  jacobian[!is.finite(jacobian)] = NA_real_
  mcse * jacobian
}


# ------------------------------------------------------------------
# boundary_distance
# ------------------------------------------------------------------
# Distance from a log10 Bayes factor to the nearer of the two verdict
# boundaries, measured in standard errors of that log10 Bayes factor. The
# standard errors arrive on the logit (natural log odds) scale, so they are
# divided by log(10) to match.
#
# @param log10_bf  Log10 inclusion Bayes factor per indicator.
# @param se_logit  Logit-scale standard error per indicator.
# @param lthr      log10 of the evidence threshold.
#
# Returns: numeric vector of distances in standard errors.
# ------------------------------------------------------------------
boundary_distance = function(log10_bf, se_logit, lthr) {
  gap = pmin(abs(log10_bf - lthr), abs(log10_bf + lthr))
  se10 = se_logit / log(10)
  out = gap / se10
  out[!is.finite(se10) | se10 <= 0] = NA_real_
  out
}


# ------------------------------------------------------------------
# compare_indicator_index
# ------------------------------------------------------------------
# Row/column positions of the bgmCompare difference indicators in the order
# generate_indicator_names() lays them out: for each variable i, the main
# difference (i, i) followed by the pairwise differences (i, j) for j > i.
#
# @param num_variables  Number of variables.
#
# Returns: a V(V+1)/2 x 2 integer matrix of (row, column) positions.
# ------------------------------------------------------------------
compare_indicator_index = function(num_variables) {
  out = vector("list", num_variables)
  for(i in seq_len(num_variables)) {
    j = if(i < num_variables) seq.int(i + 1L, num_variables) else integer(0)
    out[[i]] = cbind(i, c(i, j))
  }
  do.call(rbind, out)
}


# ------------------------------------------------------------------
# build_verdicts
# ------------------------------------------------------------------
# Shared body of the verdicts() methods: everything downstream of resolving the
# per-indicator name, Bayes factor, and inclusion probability.
#
# @param parameter  Indicator names, in the raw draws' order.
# @param log10_bf   Log10 inclusion Bayes factor per indicator.
# @param pip        Rao-Blackwellized inclusion probability per indicator.
# @param mcse       Its MCSE, from the fit summary's inclusion table.
# @param draws      List of raw indicator chains (niter x nparam 0/1 matrices).
# @param evidence_threshold  The Bayes-factor threshold.
#
# Returns: a bgms_verdicts data frame.
# ------------------------------------------------------------------
build_verdicts = function(parameter, log10_bf, pip, mcse, draws, evidence_threshold) {
  lthr = log10(evidence_threshold)

  se_two_state = two_state_se_logit(draws)
  se_rb = rb_se_logit(mcse, pip)

  d_two_state = boundary_distance(log10_bf, se_two_state, lthr)
  d_rb = boundary_distance(log10_bf, se_rb, lthr)

  # Union rule: an edge is fragile when EITHER standard error places a verdict
  # boundary within two of them. Neither alone catches every Monte Carlo
  # verdict error (two-state 0.74, Rao-Blackwellized 0.94 recall over 37,010
  # graded edge-fits); the union catches all 66 at a 3.0% false-alarm rate.
  # A missing distance abstains rather than votes.
  fragile = (!is.na(d_two_state) & d_two_state < 2) |
    (!is.na(d_rb) & d_rb < 2)

  out = data.frame(
    parameter = parameter,
    pip = pip,
    bf = 10^log10_bf,
    log10_bf = log10_bf,
    verdict = factor(
      verdict_from_lbf(log10_bf, lthr),
      levels = c("presence", "undecided", "absence")
    ),
    se_two_state = se_two_state,
    se_rb = se_rb,
    distance_two_state = d_two_state,
    distance_rb = d_rb,
    fragile = fragile,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  structure(
    out,
    class = c("bgms_verdicts", "data.frame"),
    evidence_threshold = evidence_threshold
  )
}


#' @title Edge Verdicts and Their Monte Carlo Fragility
#'
#' @description
#' Reads each edge (or difference) indicator's inclusion Bayes factor as a
#' three-way verdict -- evidence of presence, undecided, evidence of absence --
#' and flags the verdicts that a rerun of the sampler could change.
#'
#' @param bgms_object A fitted model object of class `bgms` (from [bgm()] with
#'   `edge_selection = TRUE`) or `bgmCompare` (from [bgmCompare()] with
#'   `difference_selection = TRUE`).
#' @param evidence_threshold Numeric > 1; the inclusion Bayes factor above
#'   which an edge is called present. Its reciprocal is the threshold below
#'   which an edge is called absent. Default `10`.
#' @param ... Passed to methods.
#'
#' @return A data frame of class `bgms_verdicts`, one row per indicator, in the
#'   order of the fit's raw indicator draws, with columns:
#'   \describe{
#'     \item{parameter}{Indicator name, as in `summary(fit)$indicator`.}
#'     \item{pip}{Rao-Blackwellized posterior inclusion probability.}
#'     \item{bf, log10_bf}{Inclusion Bayes factor and its base-10 logarithm,
#'       from [extract_inclusion_bf()]. `log10_bf` stays finite where `bf`
#'       saturates at `0` or `Inf`.}
#'     \item{verdict}{Factor with levels `presence`, `undecided`, `absence`,
#'       and `NA` for indicators that were never updated (main-effect
#'       differences under `main_difference_selection = FALSE`).}
#'     \item{se_two_state, se_rb}{Standard errors of the logit inclusion
#'       probability, from the Jeffreys-smoothed two-state model of the
#'       indicator chain and from the Rao-Blackwellized draws. `se_rb` is `NA`
#'       where the Rao-Blackwellized draws are constant to double precision.}
#'     \item{distance_two_state, distance_rb}{Distance from `log10_bf` to the
#'       nearer verdict boundary, in units of each standard error.}
#'     \item{fragile}{`TRUE` when either distance is below 2.}
#'   }
#'   The evidence threshold is attached as an attribute.
#'
#' @details
#' Monte Carlo verdict errors are a boundary phenomenon. In a known-truth
#' calibration study of 37,010 graded edge-fits across ordinal, binary, and
#' Gaussian graphical models, every one of the 66 verdict errors sat within
#' 0.25 of a threshold on the log10 Bayes factor scale, and no edge further out
#' was ever misclassified. The fragility flag turns that into a per-edge
#' statement: an edge is fragile when a verdict boundary lies within two
#' standard errors of the estimated evidence, which is the regime where the
#' reported verdict rests on Monte Carlo noise.
#'
#' Two standard errors are computed because neither catches every error alone.
#' The two-state standard error models the binary indicator chain as a
#' first-order two-state Markov chain with Jeffreys-smoothed transition rates,
#' which keeps it defined when the indicator never flips; on its own it caught
#' 74% of the study's verdict errors. The Rao-Blackwellized standard error is
#' the Monte Carlo standard error of the one-step inclusion draws, carried to
#' the logit scale by the delta method; on its own it caught 94%. Flagging when
#' either places a boundary within two standard errors caught all 66, at the
#' cost of also flagging 3.0% of correct verdicts. The union transfers across
#' model types, so the flag applies to Gaussian graphical models as it does to
#' ordinal and binary ones.
#'
#' A fragile verdict is not a wrong verdict; it is a verdict the run is too
#' short to settle. The remedy is more sampling iterations.
#'
#' @examples
#' \donttest{
#' fit = bgm(Wenchuan[, 1:5], display_progress = "none")
#' v = verdicts(fit)
#' v[v$fragile, ]
#'
#' # A stricter reading of the evidence moves the boundaries, and with them
#' # which edges sit close enough to one to be fragile.
#' verdicts(fit, evidence_threshold = 30)
#' }
#'
#' @seealso [extract_inclusion_bf()] for the Bayes factors,
#'   [extract_posterior_inclusion_probabilities()] for the inclusion
#'   probabilities, [prior_sensitivity_check()] for whether a verdict depends
#'   on the prior rather than on the run length
#' @family posterior-methods
#' @export
verdicts = function(bgms_object, evidence_threshold = 10, ...) {
  UseMethod("verdicts")
}

#' @inheritParams verdicts
#' @exportS3Method
#' @noRd
verdicts.bgms = function(bgms_object, evidence_threshold = 10, ...) {
  check_evidence_threshold(evidence_threshold)
  ensure_summaries(bgms_object)

  arguments = extract_arguments(bgms_object)
  if(!isTRUE(arguments$edge_selection)) {
    stop(
      "Edge verdicts require edge selection. Refit with bgm(edge_selection = ",
      "TRUE); without it every edge is in the model and there is no inclusion ",
      "Bayes factor to read."
    )
  }

  raw = get_raw_samples(bgms_object)
  summary_indicator = bgms_object$posterior_summary_indicator

  num_variables = nrow(bgms_object$posterior_mean_indicator)
  idx = which(upper.tri(matrix(0, num_variables, num_variables)), arr.ind = TRUE)
  idx = idx[order(idx[, "row"], idx[, "col"]), , drop = FALSE]

  log10_bf = extract_inclusion_bf(bgms_object, log = TRUE)[idx] / log(10)
  pip = bgms_object$posterior_mean_indicator[idx]

  build_verdicts(
    parameter = raw$parameter_names$indicator,
    log10_bf = log10_bf,
    pip = pip,
    mcse = summary_indicator[["mcse"]],
    draws = raw$indicator,
    evidence_threshold = evidence_threshold
  )
}

#' @inheritParams verdicts
#' @exportS3Method
#' @noRd
verdicts.bgmCompare = function(bgms_object, evidence_threshold = 10, ...) {
  check_evidence_threshold(evidence_threshold)
  ensure_summaries(bgms_object)

  arguments = extract_arguments(bgms_object)
  if(!isTRUE(arguments$difference_selection)) {
    stop(
      "Difference verdicts require difference selection. Refit with ",
      "bgmCompare(difference_selection = TRUE); without it every difference ",
      "is in the model and there is no inclusion Bayes factor to read."
    )
  }

  raw = get_raw_samples(bgms_object)
  summary_indicator = bgms_object$posterior_summary_indicator

  bf = extract_inclusion_bf(bgms_object, log = TRUE)
  idx = compare_indicator_index(nrow(bf))

  build_verdicts(
    parameter = raw$parameter_names$indicator,
    log10_bf = bf[idx] / log(10),
    pip = extract_posterior_inclusion_probabilities(bgms_object)[idx],
    mcse = summary_indicator[["mcse"]],
    draws = raw$indicator,
    evidence_threshold = evidence_threshold
  )
}


# ------------------------------------------------------------------
# check_evidence_threshold
# ------------------------------------------------------------------
# Validate the Bayes-factor threshold shared by verdicts() and plot().
#
# @param evidence_threshold  The value to check.
#
# Returns: invisible(TRUE), or stops.
# ------------------------------------------------------------------
check_evidence_threshold = function(evidence_threshold) {
  if(!is.numeric(evidence_threshold) || length(evidence_threshold) != 1L ||
    !is.finite(evidence_threshold) || evidence_threshold <= 1) {
    stop(
      "Argument 'evidence_threshold' must be a single number greater than 1. ",
      "It is the inclusion Bayes factor above which an edge counts as present; ",
      "its reciprocal is the threshold for absence."
    )
  }
  invisible(TRUE)
}


#' @title Print Edge Verdicts
#'
#' @description
#' Prints the verdict tally, the leading rows of the verdict table, and a
#' warning line when any verdict is Monte-Carlo fragile.
#'
#' @param x An object of class `bgms_verdicts`, from [verdicts()].
#' @param digits Number of digits for the printed numeric columns. Default `3`.
#' @param max_rows Number of rows to print. Default `10`.
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#'
#' @seealso [verdicts()]
#' @family posterior-methods
#' @export
print.bgms_verdicts = function(x, digits = 3, max_rows = 10L, ...) {
  # Subsetting columns keeps the class but not the table: print what is left as
  # the plain data frame it has become, rather than failing on a missing column.
  required = c("parameter", "pip", "log10_bf", "verdict", "fragile")
  if(!all(required %in% names(x))) {
    print(as.data.frame(x), ...)
    return(invisible(x))
  }

  threshold = attr(x, "evidence_threshold")
  cat(sprintf(
    "Edge verdicts at an inclusion Bayes factor of %g (and %g for absence):\n\n",
    threshold, 1 / threshold
  ))

  tally = table(x$verdict)
  cat(sprintf(
    "  presence %d | undecided %d | absence %d   (%d indicators)\n",
    tally[["presence"]], tally[["undecided"]], tally[["absence"]], nrow(x)
  ))
  n_unselected = sum(is.na(x$verdict))
  if(n_unselected > 0) {
    cat(sprintf(
      "  %d indicator(s) were never updated and carry no verdict.\n", n_unselected
    ))
  }
  cat("\n")

  shown = utils::head(x, max_rows)
  body = data.frame(
    parameter = shown$parameter,
    pip = round(shown$pip, digits),
    log10_bf = round(shown$log10_bf, digits),
    verdict = as.character(shown$verdict),
    fragile = shown$fragile,
    check.names = FALSE
  )
  print(body, row.names = FALSE)
  if(nrow(x) > max_rows) {
    cat(sprintf("... (%d more rows)\n", nrow(x) - max_rows))
  }

  n_fragile = sum(x$fragile)
  if(n_fragile > 0) {
    cat(sprintf(
      "\n%d %s Monte-Carlo fragile: a verdict boundary lies within two standard\nerrors of the evidence, so the verdict could change on a rerun. Run longer.\n",
      n_fragile, if(n_fragile == 1L) "verdict is" else "verdicts are"
    ))
  }
  invisible(x)
}

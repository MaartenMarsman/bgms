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
# Distance from a log Bayes factor to the nearer of the two verdict
# boundaries, measured in standard errors of that log Bayes factor. Both the
# Bayes factor and the standard errors are on the natural log (logit) scale,
# so the two are already commensurable.
#
# @param log_bf    Natural log inclusion Bayes factor per indicator.
# @param se_logit  Logit-scale standard error per indicator.
# @param lthr      Natural log of the evidence threshold.
#
# Returns: numeric vector of distances in standard errors.
# ------------------------------------------------------------------
boundary_distance = function(log_bf, se_logit, lthr) {
  gap = pmin(abs(log_bf - lthr), abs(log_bf + lthr))
  out = gap / se_logit
  out[!is.finite(se_logit) | se_logit <= 0] = NA_real_
  out
}


# ------------------------------------------------------------------
# format_log_bf
# ------------------------------------------------------------------
# Display form of a log Bayes factor: the rounded value inside the reporting
# cap, and an inequality outside it. Past the cap the digits carry no
# information a longer run would reproduce, and an infinite Bayes factor has
# no value to print at all.
#
# @param log_bf  Natural log inclusion Bayes factor (scalar).
# @param cap     Magnitude past which an inequality is printed.
#
# Returns: a length-one character string.
# ------------------------------------------------------------------
format_log_bf = function(log_bf, cap = 1e4) {
  if(is.na(log_bf)) {
    return("NA")
  }
  if(log_bf > cap) {
    return(sprintf("> %s", format(cap, big.mark = ",", scientific = FALSE)))
  }
  if(log_bf < -cap) {
    return(sprintf("< -%s", format(cap, big.mark = ",", scientific = FALSE)))
  }
  # A Bayes factor that rounds to nothing is nothing, not "-0.0": the sign of a
  # rounded-away quantity is not information the run established. Same rule
  # estimate_lines() already applies to a weight (R/plot_bgms.R).
  shown = round(log_bf, 1)
  if(shown == 0) shown = 0
  sprintf("= %.1f", shown)
}


# ------------------------------------------------------------------
# indicator_pair_index
# ------------------------------------------------------------------
# Row/column positions of a bgm fit's edge indicators, in the order its raw
# indicator draws lay them out.
#
# A GGM or ordinal fit lays its indicators out as the row-major upper triangle
# of the variable order. A mixed fit lays them out by block --
# discrete-discrete, then continuous-continuous, then cross -- which is a
# different permutation whenever the discrete and continuous columns interleave.
# Reading a symmetric matrix (Bayes factors, inclusion probabilities) at the
# wrong one attaches every number to the wrong edge.
#
# @param bgms_object   A bgms fit.
# @param num_variables Number of variables.
#
# Returns: an E x 2 integer matrix of (row, column) positions.
# ------------------------------------------------------------------
indicator_pair_index = function(bgms_object, num_variables) {
  spec = get_fit_spec(bgms_object)
  if(!is.null(spec) && identical(spec$model_type, "mixed_mrf")) {
    d = spec$data
    num_pairs = num_variables * (num_variables - 1L) / 2L
    blank = character(num_variables)
    # fill_mixed_symmetric() is the layout the raw draws follow, so filling it
    # with the draw positions reads that layout back off as an index.
    positions = fill_mixed_symmetric(
      seq_len(num_pairs), d$num_discrete, d$num_continuous,
      d$discrete_indices, d$continuous_indices, list(blank, blank)
    )
    upper = which(upper.tri(positions), arr.ind = TRUE)
    return(upper[order(positions[upper]), , drop = FALSE])
  }
  idx = which(upper.tri(matrix(0, num_variables, num_variables)), arr.ind = TRUE)
  idx[order(idx[, "row"], idx[, "col"]), , drop = FALSE]
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
# @param log_bf     Natural log inclusion Bayes factor per indicator.
# @param pip        Rao-Blackwellized inclusion probability per indicator.
# @param mcse       Its MCSE, from the fit summary's inclusion table.
# @param draws      List of raw indicator chains (niter x nparam 0/1 matrices).
# @param evidence_threshold  The Bayes-factor threshold.
# @param flag_validated  Whether the fragility flag's operating point was
#   established for this kind of indicator. TRUE for single-network edge
#   indicators; FALSE for bgmCompare difference indicators, which no arm of the
#   calibration study covered.
#
# Returns: a bgms_verdicts data frame.
# ------------------------------------------------------------------
build_verdicts = function(parameter, log_bf, pip, mcse, draws, evidence_threshold,
                          flag_validated) {
  lthr = log(evidence_threshold)

  se_two_state = two_state_se_logit(draws)
  se_rb = rb_se_logit(mcse, pip)

  d_two_state = boundary_distance(log_bf, se_two_state, lthr)
  d_rb = boundary_distance(log_bf, se_rb, lthr)

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
    bf = exp(log_bf),
    log_bf = log_bf,
    verdict = factor(
      verdict_from_lbf(log_bf, lthr),
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
    evidence_threshold = evidence_threshold,
    flag_validated = flag_validated
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
#'     \item{bf, log_bf}{Inclusion Bayes factor and its natural logarithm,
#'       from [extract_inclusion_bf()], which returns the same natural-log
#'       scale under `log = TRUE`. `log_bf` stays finite where `bf` saturates
#'       at `0` or `Inf`.}
#'     \item{verdict}{Factor with levels `presence`, `undecided`, `absence`,
#'       and `NA` for indicators that were never updated (main-effect
#'       differences under `main_difference_selection = FALSE`).}
#'     \item{se_two_state, se_rb}{Standard errors of the logit inclusion
#'       probability, from the Jeffreys-smoothed two-state model of the
#'       indicator chain and from the Rao-Blackwellized draws. `se_rb` is `NA`
#'       where the Rao-Blackwellized draws are constant to double precision.}
#'     \item{distance_two_state, distance_rb}{Distance from `log_bf` to the
#'       nearer verdict boundary, in units of each standard error.}
#'     \item{fragile}{`TRUE` when either distance is below 2.}
#'   }
#'   The evidence threshold is attached as the `evidence_threshold` attribute,
#'   and whether the fragility flag's operating point covers this kind of
#'   indicator as `flag_validated`. A `bgmCompare` fit additionally carries
#'   `difference_fit`, and, under `main_difference_selection = FALSE`,
#'   `unselected_main`, marking the main-effect rows the print leaves out of
#'   its counts because no indicator exists for them.
#'
#' @details
#' Monte Carlo verdict errors are a boundary phenomenon. In a known-truth
#' calibration study of 37,010 graded edge-fits across ordinal, binary, and
#' Gaussian graphical models, every one of the 66 verdict errors sat within
#' 0.58 of a threshold on the log Bayes factor scale, and no edge further out
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
#' Every arm of that study fitted a single network, so the operating point
#' applies to the edge indicators of [bgm()]. On the difference indicators of
#' [bgmCompare()] the flag still marks verdicts sitting near a boundary, but no
#' study has measured what share of difference-verdict errors it catches or how
#' many correct verdicts it rejects; the print method says so, and the returned
#' object carries a `flag_validated` attribute.
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
  idx = indicator_pair_index(bgms_object, num_variables)

  log_bf = extract_inclusion_bf(bgms_object, log = TRUE)[idx]
  pip = bgms_object$posterior_mean_indicator[idx]

  build_verdicts(
    parameter = raw$parameter_names$indicator,
    log_bf = log_bf,
    pip = pip,
    mcse = summary_indicator[["mcse"]],
    draws = raw$indicator,
    evidence_threshold = evidence_threshold,
    flag_validated = TRUE
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

  out = build_verdicts(
    # The compare producer names this element `indicators`
    # (generate_param_names_bgmCompare); spell it exactly rather than leaning on
    # `$` partial matching, which would break the day a sibling name shares the
    # prefix.
    parameter = raw$parameter_names$indicators,
    log_bf = bf[idx],
    pip = extract_posterior_inclusion_probabilities(bgms_object)[idx],
    mcse = summary_indicator[["mcse"]],
    draws = raw$indicator,
    evidence_threshold = evidence_threshold,
    flag_validated = FALSE
  )
  # The print method needs to know these are difference verdicts (for the
  # scale-contingency caveat) and, under main_difference_selection = FALSE,
  # which rows are main-effect differences with no indicator to read.
  attr(out, "difference_fit") = TRUE
  if(!isTRUE(arguments$main_difference_selection)) {
    attr(out, "unselected_main") = idx[, 1] == idx[, 2]
  }
  out
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
  # Subsetting keeps the class but not the table: print what is left as the
  # plain data frame it has become, rather than failing on what is missing.
  # Two things can go: `[.data.frame` drops a column when `j` is supplied, and
  # it drops the evidence_threshold attribute whenever `j` is supplied at all --
  # so subset(x, ...), which always supplies `j`, keeps every display column and
  # still leaves nothing to state the boundaries with. Both are checked.
  required = c("parameter", "pip", "log_bf", "verdict", "fragile")
  threshold = attr(x, "evidence_threshold")
  if(!all(required %in% names(x)) || is.null(threshold)) {
    print(as.data.frame(x), ...)
    return(invisible(x))
  }

  cat(sprintf(
    "Edge verdicts at an inclusion Bayes factor of %g (and %g for absence):\n",
    threshold, signif(1 / threshold, 3)
  ))
  # The table reports the evidence as a natural log Bayes factor, so the
  # boundaries are stated in that unit rather than left to be converted.
  cat(sprintf(
    "presence: log BF > %.2f; absence: log BF < -%.2f\n\n",
    log(threshold), log(threshold)
  ))

  # A main-effect difference outside selection has no indicator, so there is
  # no verdict to tabulate; a bgmCompare fit marks those rows and they are
  # left out of the counts and the table.
  hidden = attr(x, "unselected_main")
  rows = if(is.null(hidden)) x else x[!hidden, , drop = FALSE]

  tally = table(rows$verdict)
  cat(sprintf(
    "  presence %d | undecided %d | absence %d   (%d %s)\n",
    tally[["presence"]], tally[["undecided"]], tally[["absence"]], nrow(rows),
    if(nrow(rows) == 1L) "indicator" else "indicators"
  ))
  if(!is.null(hidden) && any(hidden)) {
    cat(
      "  Main-effect differences are not under selection",
      "(main_difference_selection = FALSE).\n"
    )
  }
  n_unselected = sum(is.na(rows$verdict))
  if(n_unselected > 0) {
    cat(sprintf(
      "  %d indicator(s) were never updated and carry no verdict.\n", n_unselected
    ))
  }
  cat("\n")

  shown = utils::head(rows, max_rows)
  body = data.frame(
    parameter = shown$parameter,
    pip = round(shown$pip, digits),
    log_bf = round(shown$log_bf, digits),
    verdict = as.character(shown$verdict),
    fragile = shown$fragile,
    check.names = FALSE
  )
  print(body, row.names = FALSE)
  if(nrow(rows) > max_rows) {
    cat(sprintf("... (%d more rows)\n", nrow(rows) - max_rows))
  }

  n_fragile = sum(rows$fragile)
  if(n_fragile > 0) {
    cat(sprintf(
      "\n%d %s Monte-Carlo fragile: a verdict boundary lies within two standard\nerrors of the evidence, so the verdict could change on a rerun. Consider a\nlonger run.\n",
      n_fragile, if(n_fragile == 1L) "verdict is" else "verdicts are"
    ))
  }
  # The flag's operating point was established on single-network edge
  # indicators. Borrowing those numbers for difference indicators would report a
  # calibration that no study arm has measured.
  if(isFALSE(attr(x, "flag_validated"))) {
    cat(
      "\nThe fragility flag is not validated for difference indicators: its\n",
      "operating point was established on single-network edge indicators only.\n",
      "Read it as an indication that a verdict sits near a boundary, not as a\n",
      "calibrated error rate.\n",
      sep = ""
    )
  }
  if(isTRUE(attr(x, "difference_fit"))) {
    cat(
      "\nDifference verdicts are scale-contingent: group differences are priced\n",
      "on the association scale through difference_scale, and the calibration\n",
      "of that default is under study, so a verdict close to a decision\n",
      "threshold can move with the scale.\n",
      sep = ""
    )
  }
  invisible(x)
}

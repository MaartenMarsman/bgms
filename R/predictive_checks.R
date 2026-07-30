# ==============================================================================
# Posterior predictive checks
# ==============================================================================
#
# Replicated datasets are drawn one per posterior draw, so parameter and
# structure uncertainty propagate into the replicates. Three statistics, from
# local to global: the item margins, the pairwise associations, and the
# sum-score distribution, which no single pair pins down and which is therefore
# sensitive to dependence the margins cannot show.
# ==============================================================================


# ------------------------------------------------------------------
# ppc_observed_data
# ------------------------------------------------------------------
# The data the model saw, on the scale simulate() returns: the fit stores it
# zero-based with the original levels alongside, while simulate() maps back.
#
# @param bgms_object  A bgms fit.
#
# Returns: an integer matrix of observed responses on the original scale.
# ------------------------------------------------------------------
ppc_observed_data = function(bgms_object) {
  spec = get_fit_spec(bgms_object)
  x = spec$data$x
  levels_list = spec$data$category_levels

  out = x
  for(j in seq_len(ncol(x))) {
    out[, j] = levels_list[[j]][x[, j] + 1L]
  }
  colnames(out) = spec$data$data_columnnames
  out
}


# ------------------------------------------------------------------
# ppc_statistic_fun
# ------------------------------------------------------------------
# The function computing one statistic from a dataset, together with the labels
# of its elements. Both the observed data and every replicate go through the
# same closure, so the element order is guaranteed to line up.
#
# @param statistic  One of "margins", "pairwise", "sumscore".
# @param observed   The observed data matrix (fixes the support and the labels).
# @param levels_list  Per-variable category levels.
#
# Returns: a list with `fun` (matrix -> numeric vector) and `labels` (a data
#   frame of element identifiers).
# ------------------------------------------------------------------
ppc_statistic_fun = function(statistic, observed, levels_list) {
  variables = colnames(observed)
  num_variables = ncol(observed)

  if(statistic == "margins") {
    labels = data.frame(
      variable = rep(variables, lengths(levels_list)),
      category = unlist(levels_list, use.names = FALSE),
      stringsAsFactors = FALSE
    )
    fun = function(x) {
      unlist(lapply(seq_len(num_variables), function(j) {
        vapply(levels_list[[j]], function(k) mean(x[, j] == k), numeric(1))
      }), use.names = FALSE)
    }
  } else if(statistic == "pairwise") {
    pairs = which(upper.tri(matrix(0, num_variables, num_variables)), arr.ind = TRUE)
    pairs = pairs[order(pairs[, "row"], pairs[, "col"]), , drop = FALSE]
    labels = data.frame(
      variable1 = variables[pairs[, 1]],
      variable2 = variables[pairs[, 2]],
      stringsAsFactors = FALSE
    )
    fun = function(x) stats::cor(x, method = "spearman")[cbind(pairs[, 1], pairs[, 2])]
  } else {
    support = seq.int(
      sum(vapply(levels_list, min, numeric(1))),
      sum(vapply(levels_list, max, numeric(1)))
    )
    labels = data.frame(score = support)
    fun = function(x) {
      tabulate(match(rowSums(x), support), nbins = length(support)) / nrow(x)
    }
  }
  list(fun = fun, labels = labels)
}


#' @title Posterior Predictive Check
#'
#' @description
#' Compares observed data summaries with their posterior predictive
#' distributions, drawing one replicated dataset per posterior draw so that
#' parameter and structure uncertainty propagate into the replicates.
#'
#' @param bgms_object A fitted model object of class `bgms` (from [bgm()]) with
#'   discrete variables.
#' @param statistic Character vector; which statistics to check. Any of
#'   `"margins"` (the category proportions of each variable), `"pairwise"` (the
#'   Spearman correlation of each variable pair), and `"sumscore"` (the
#'   distribution of the total score across variables). Defaults to all three.
#' @param ndraws Number of posterior draws, and therefore of replicated
#'   datasets. Default `500`.
#' @param probs Numeric of length two; the predictive-interval quantiles.
#'   Default `c(0.025, 0.975)`.
#' @param seed Optional integer seed for the replicate simulation.
#' @param cores Number of cores for the simulation. Defaults to
#'   `parallel::detectCores()`.
#' @param display_progress Progress display passed to [simulate()]; one of
#'   `"per-chain"`, `"total"`, or `"none"`.
#' @param ... Passed to methods.
#'
#' @return An object of class `bgms_ppc`, a list with:
#'   \describe{
#'     \item{statistics}{A named list, one data frame per requested statistic,
#'       with the element identifiers, the `observed` value, the predictive
#'       mean `predicted`, the interval bounds `lower` and `upper`, and
#'       `covered`.}
#'     \item{coverage}{One row per statistic: the number of elements, how many
#'       the intervals covered, and how many a calibrated model is expected to
#'       cover.}
#'     \item{ndraws, nsim, probs}{The settings the check ran under.}
#'   }
#'
#' @details
#' The check reads on two layers. Coverage is the count of elements whose
#' observed value falls outside its predictive interval; with a 95% interval,
#' about 5% of elements land outside one under a model that fits, so a handful
#' of misses is the expected picture rather than a failure. The second layer is
#' the pattern: a systematic tilt or offset of the observed values against
#' their predictive means is misfit even when coverage looks unremarkable, and
#' it is what the plot shows and a coverage count cannot.
#'
#' The statistics run from local to global. The margins involve one variable at
#' a time and are usually reproduced by any model that fits the data's
#' univariate shape. The pairwise correlations involve two. The sum-score
#' distribution depends on the joint dependence among all variables, so it is
#' the statistic that can fail while the other two pass.
#'
#' Discrete variables only in this release; a fit containing continuous
#' variables errors, because these three statistics have no direct analogue
#' there.
#'
#' @examples
#' \donttest{
#' fit = bgm(Wenchuan[, 1:5], display_progress = "none")
#' check = posterior_predictive_check(fit, ndraws = 100, display_progress = "none")
#' check
#' plot(check)
#' }
#'
#' @seealso [simulate.bgms()] for the replicate draws themselves,
#'   [calibration_check()] for the conditional predictions
#' @family diagnostics
#' @export
posterior_predictive_check = function(bgms_object,
                                      statistic = c("margins", "pairwise", "sumscore"),
                                      ndraws = 500,
                                      probs = c(0.025, 0.975),
                                      seed = NULL,
                                      cores = parallel::detectCores(),
                                      display_progress = c("per-chain", "total", "none"),
                                      ...) {
  UseMethod("posterior_predictive_check")
}

#' @inheritParams posterior_predictive_check
#' @exportS3Method
#' @noRd
posterior_predictive_check.bgms = function(bgms_object,
                                           statistic = c("margins", "pairwise", "sumscore"),
                                           ndraws = 500,
                                           probs = c(0.025, 0.975),
                                           seed = NULL,
                                           cores = parallel::detectCores(),
                                           display_progress = c("per-chain", "total", "none"),
                                           ...) {
  statistic = match.arg(statistic, several.ok = TRUE)
  check_interval_probs(probs)

  arguments = extract_arguments(bgms_object)
  if(any(arguments$variable_type == "continuous")) {
    stop(
      "Posterior predictive checks currently cover discrete variables only. ",
      "This fit contains continuous variables, for which the category ",
      "proportions and the sum-score distribution have no direct analogue. ",
      "Use simulate() with a statistic of your own in the meantime."
    )
  }

  observed = ppc_observed_data(bgms_object)
  if(anyNA(observed)) {
    stop(
      "The fitted data contain missing values, so the observed statistics and ",
      "the complete replicated datasets are not comparable. Refit with ",
      "na_action = \"listwise\" to run the check on complete cases."
    )
  }

  levels_list = get_fit_spec(bgms_object)$data$category_levels
  nsim = nrow(observed)

  replicates = stats::simulate(
    bgms_object,
    nsim = nsim, method = "posterior-sample", ndraws = ndraws,
    seed = seed, cores = cores, display_progress = display_progress
  )

  statistics = list()
  coverage = list()
  for(s in statistic) {
    spec = ppc_statistic_fun(s, observed, levels_list)
    obs = spec$fun(observed)
    rep_values = vapply(replicates, spec$fun, obs)
    if(is.null(dim(rep_values))) rep_values = matrix(rep_values, nrow = length(obs))

    bounds = apply(rep_values, 1, stats::quantile, probs = probs, na.rm = TRUE)
    covered = obs >= bounds[1, ] & obs <= bounds[2, ]

    statistics[[s]] = data.frame(
      spec$labels,
      observed = obs,
      predicted = rowMeans(rep_values, na.rm = TRUE),
      lower = bounds[1, ],
      upper = bounds[2, ],
      covered = covered,
      row.names = NULL,
      stringsAsFactors = FALSE
    )
    coverage[[s]] = data.frame(
      statistic = s,
      total = length(obs),
      covered = sum(covered),
      expected = diff(probs) * length(obs),
      stringsAsFactors = FALSE
    )
  }

  structure(
    list(
      statistics = statistics,
      coverage = do.call(rbind, c(coverage, list(make.row.names = FALSE))),
      ndraws = length(replicates),
      nsim = nsim,
      probs = probs
    ),
    class = "bgms_ppc"
  )
}


# ------------------------------------------------------------------
# check_interval_probs
# ------------------------------------------------------------------
# Validate an interval-quantile pair shared by the predictive checks.
#
# @param probs  The value to check.
#
# Returns: invisible(TRUE), or stops.
# ------------------------------------------------------------------
check_interval_probs = function(probs) {
  if(!is.numeric(probs) || length(probs) != 2L || anyNA(probs) ||
    any(probs < 0) || any(probs > 1) || probs[1] >= probs[2]) {
    stop(
      "Argument 'probs' must be two increasing probabilities, the lower and ",
      "upper quantiles of the interval, e.g. c(0.025, 0.975)."
    )
  }
  invisible(TRUE)
}


#' @title Print a Posterior Predictive Check
#'
#' @description
#' Prints the coverage of each checked statistic against the coverage a
#' calibrated model is expected to reach.
#'
#' @param x An object of class `bgms_ppc`, from [posterior_predictive_check()].
#' @param digits Number of digits for the printed percentages. Default `1`.
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#'
#' @seealso [posterior_predictive_check()]
#' @family diagnostics
#' @export
print.bgms_ppc = function(x, digits = 1, ...) {
  level = 100 * diff(x$probs)
  cat(sprintf(
    "Posterior predictive check: %d replicated datasets of %d cases, %g%% intervals\n\n",
    x$ndraws, x$nsim, level
  ))

  body = x$coverage
  body$percent = round(100 * body$covered / body$total, digits)
  body$expected = round(body$expected, digits)
  print(body[, c("statistic", "total", "covered", "expected", "percent")],
    row.names = FALSE
  )
  cat(sprintf(
    "\n'expected' is the coverage a calibrated model reaches: about %g%% of\nelements land outside their interval even when the model fits. Read the\nplot for the pattern, which a coverage count cannot show.\n",
    100 - level
  ))
  invisible(x)
}


# ------------------------------------------------------------------
# ppc_panel_scatter
# ------------------------------------------------------------------
# The interval-scatter panel, for statistics with many independent elements:
# observed against predictive mean, a segment for the predictive interval, and
# the identity line. Coverage is the count of segments missing the diagonal;
# misfit is a tilt or offset of the whole cloud away from it.
#
# @param df       One statistic's data frame.
# @param title    Panel title.
# @param level    Interval level, in percent.
# @param labels   Optional character vector of per-element labels, drawn on the
#   elements outside their interval only.
# @param max_labels  Cap on the number of labels drawn.
#
# Returns: invisible(NULL).
# ------------------------------------------------------------------
ppc_panel_scatter = function(df, title, level, labels = NULL, max_labels = 6L) {
  ink = "grey25"
  muted = "grey55"
  accent = mover_palette()[2]

  lim = range(c(df$observed, df$lower, df$upper), finite = TRUE)
  graphics::plot(NA, NA,
    xlim = lim, ylim = lim, axes = FALSE,
    xlab = sprintf("Posterior predictive value (%g%% interval)", level),
    ylab = "Observed value", main = title
  )
  graphics::axis(1, col = muted, col.ticks = muted)
  graphics::axis(2, col = muted, col.ticks = muted, las = 1)
  graphics::abline(0, 1, col = muted, lwd = 0.8)

  graphics::segments(df$lower, df$observed, df$upper, df$observed,
    col = grDevices::adjustcolor(muted, 0.5), lwd = 1
  )
  miss = !df$covered
  graphics::points(df$predicted[!miss], df$observed[!miss],
    pch = 16, cex = 0.7, col = ink
  )
  # Identity never rides on colour alone: a missed element takes its own
  # symbol as well as the accent.
  graphics::points(df$predicted[miss], df$observed[miss],
    pch = 17, cex = 0.9, col = accent
  )

  if(!is.null(labels) && any(miss) && max_labels > 0L) {
    idx = which(miss)
    idx = idx[order(abs(df$observed[idx] - df$predicted[idx]), decreasing = TRUE)]
    idx = utils::head(idx, max_labels)
    # Missed elements cluster, so their labels would land on top of one
    # another; nudge them apart vertically and connect each to its point.
    idx = idx[order(df$observed[idx])]
    usr = graphics::par("usr")
    label_y = spread_labels(df$observed[idx], gap = 0.05 * diff(usr[3:4]))
    graphics::segments(df$predicted[idx], df$observed[idx], df$predicted[idx], label_y,
      col = grDevices::adjustcolor(muted, 0.6), lwd = 0.7
    )
    graphics::text(df$predicted[idx], label_y, labels[idx],
      pos = 4, offset = 0.35, cex = 0.65, col = ink
    )
  }

  graphics::mtext(
    sprintf(
      "%d of %d inside the interval (about %d expected)",
      sum(df$covered), nrow(df), round(level / 100 * nrow(df))
    ),
    side = 3, line = 0.1, cex = 0.72, col = muted
  )
  invisible(NULL)
}


# ------------------------------------------------------------------
# ppc_panel_band
# ------------------------------------------------------------------
# The predictive-band panel, for a distribution over an ordered support: the
# interval as a band across the support with the observed proportions on top.
#
# @param df     The sumscore data frame.
# @param level  Interval level, in percent.
#
# Returns: invisible(NULL).
# ------------------------------------------------------------------
ppc_panel_band = function(df, level) {
  ink = "grey25"
  muted = "grey55"
  accent = mover_palette()[2]

  graphics::plot(NA, NA,
    xlim = range(df$score), ylim = c(0, max(c(df$upper, df$observed))),
    axes = FALSE, xlab = "Total score", ylab = "Proportion of cases",
    main = "Sum-score distribution"
  )
  graphics::axis(1, col = muted, col.ticks = muted)
  graphics::axis(2, col = muted, col.ticks = muted, las = 1)

  graphics::polygon(
    c(df$score, rev(df$score)), c(df$lower, rev(df$upper)),
    col = grDevices::adjustcolor(muted, 0.22), border = NA
  )
  miss = !df$covered
  graphics::points(df$score[!miss], df$observed[!miss], pch = 16, cex = 0.7, col = ink)
  graphics::points(df$score[miss], df$observed[miss], pch = 17, cex = 0.9, col = accent)

  graphics::mtext(
    sprintf("%d of %d scores inside the %g%% band", sum(df$covered), nrow(df), level),
    side = 3, line = 0.1, cex = 0.72, col = muted
  )
  invisible(NULL)
}


#' @title Plot a Posterior Predictive Check
#'
#' @description
#' One panel per checked statistic: an interval scatter of observed against
#' predicted for the element-wise statistics, and a predictive band for the
#' sum-score distribution.
#'
#' @param x An object of class `bgms_ppc`, from [posterior_predictive_check()].
#' @param max_labels Number of elements outside their interval to label by
#'   name, largest departure first. Default `6`.
#' @param ... Ignored.
#'
#' @return `x`, invisibly. Called for the side effect of drawing.
#'
#' @details
#' Elements whose observed value falls outside its predictive interval are
#' drawn in an accent colour and with their own plotting symbol, and each
#' panel's subtitle states how many a calibrated model is expected to miss, so
#' the expected handful of misses does not read as a failure. The reading that
#' matters is the departure of the cloud from the identity line.
#'
#' @examples
#' \donttest{
#' fit = bgm(Wenchuan[, 1:5], display_progress = "none")
#' plot(posterior_predictive_check(fit, ndraws = 100, display_progress = "none"))
#' }
#'
#' @seealso [posterior_predictive_check()]
#' @family diagnostics
#' @export
plot.bgms_ppc = function(x, max_labels = 6L, ...) {
  level = 100 * diff(x$probs)
  names_shown = names(x$statistics)

  old_par = graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(
    mfrow = c(1L, length(names_shown)),
    mar = c(4.3, 4.1, 3.2, 1.4), mgp = c(2.5, 0.6, 0),
    col.axis = "grey25", col.lab = "grey25", col.main = "grey25"
  )

  for(s in names_shown) {
    df = x$statistics[[s]]
    if(s == "margins") {
      ppc_panel_scatter(df, "Category proportions", level,
        labels = df$variable, max_labels = max_labels
      )
    } else if(s == "pairwise") {
      ppc_panel_scatter(df, "Pairwise associations", level,
        labels = paste(df$variable1, df$variable2, sep = "-"),
        max_labels = max_labels
      )
    } else {
      ppc_panel_band(df, level)
    }
  }
  invisible(x)
}

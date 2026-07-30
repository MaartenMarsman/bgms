# ==============================================================================
# Calibration of the conditional predictions
# ==============================================================================
#
# The model's person-level predictions are the conditional distributions
# P(x_ij | the rest of case i's responses), which is the model's own regression
# unit. Calibration is therefore assessed one variable at a time: pooling across
# variables lets miscalibrations in opposite directions cancel.
#
# For every case and every category threshold k the model issues P(x_ij <= k |
# rest) and the data record the outcome. Isotonic regression (pool-adjacent-
# violators) estimates the observed frequency as a monotone function of the
# predicted probability; a calibrated model tracks the diagonal.
#
# Origin: the CORP reliability diagram (Dimitriadis, Gneiting & Jordan) and the
# visual predictive checks of Sailynoja et al. (arXiv 2503.01509).
# ==============================================================================


# ------------------------------------------------------------------
# fitted_observed_data
# ------------------------------------------------------------------
# The data the model saw, on the scale simulate() and predict() work in: the
# fit stores it zero-based with the original levels alongside.
#
# @param bgms_object  A bgms fit.
#
# Returns: an integer matrix of observed responses on the original scale.
# ------------------------------------------------------------------
fitted_observed_data = function(bgms_object) {
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
# check_interval_probs
# ------------------------------------------------------------------
# Validate an interval-quantile pair.
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


# ------------------------------------------------------------------
# pav_curve
# ------------------------------------------------------------------
# Pool-adjacent-violators fit of y on p, read off on a common grid.
#
# @param p     Predicted probabilities.
# @param y     0/1 outcomes.
# @param grid  Grid to interpolate the fit onto.
#
# Returns: numeric vector of fitted frequencies, one per grid point.
# ------------------------------------------------------------------
pav_curve = function(p, y, grid) {
  fit = stats::isoreg(p, y)
  stats::approx(fit$x[fit$ord], fit$yf,
    xout = grid, method = "linear", rule = 2, ties = "ordered"
  )$y
}


#' @title Calibration Check
#'
#' @description
#' Reliability diagrams of the model's conditional predictions, one per
#' variable, with the consistency band a calibrated model would wander inside.
#'
#' @param bgms_object A fitted model object of class `bgms` (from [bgm()]) with
#'   discrete variables.
#' @param newdata Optional data to evaluate the predictions on, in the layout
#'   [bgm()] was given. Defaults to the data the model was fitted to.
#' @param nrep Number of resampled datasets behind the consistency band.
#'   Default `200`.
#' @param probs Numeric of length two; the band's quantiles. Default
#'   `c(0.025, 0.975)`.
#' @param grid_size Number of points the curves are read off on. Default `101`.
#' @param seed Optional integer seed for the band resampling.
#' @param ... Passed to methods.
#'
#' @return An object of class `bgms_calibration`, a list with:
#'   \describe{
#'     \item{curves}{A data frame with one row per variable and grid point:
#'       `variable`, `grid`, the isotonic fit `pav`, and the band bounds
#'       `lower` and `upper`.}
#'     \item{summary}{One row per variable, worst first: `mean_dev` and
#'       `max_dev`, the mean and maximum absolute distance of the fit from the
#'       diagonal, and `share_outside_band`.}
#'     \item{nrep, probs, grid}{The settings the check ran under.}
#'   }
#'
#' @details
#' The band resamples each case's *category* from its own predicted
#' distribution and refits the curve. It has to be built that way: the
#' cumulative threshold events of one case are nested, so resampling the events
#' independently would understate the band and make an ordinary wander look
#' like miscalibration.
#'
#' Evaluated on the fitted data the check is in-sample, and the band is the
#' reference a model that is calibrated by construction produces on the same
#' data. Supplying held-out `newdata` makes it out-of-sample and strictly
#' harder to pass.
#'
#' Calibrated conditional predictions do not imply that the model reproduces
#' the joint distribution: a model can predict each variable well from the
#' others and still understate how strongly they depend on one another. That
#' second question is answered by a display built on [simulate()], not by this
#' check.
#'
#' The check covers discrete variables, whose predictions are distributions
#' over categories that an observed category either falls in or does not. A
#' continuous variable's prediction is a density, which needs a different
#' construction.
#'
#' @examples
#' \donttest{
#' fit = bgm(Wenchuan[, 1:5], display_progress = "none")
#' check = calibration_check(fit, nrep = 50)
#' check
#' plot(check)
#' }
#'
#' @seealso [predict.bgms()] for the conditional predictions themselves,
#'   [simulate.bgms()] for replicated datasets to build a joint-level display on
#' @family diagnostics
#' @export
calibration_check = function(bgms_object,
                             newdata = NULL,
                             nrep = 200,
                             probs = c(0.025, 0.975),
                             grid_size = 101,
                             seed = NULL,
                             ...) {
  UseMethod("calibration_check")
}

#' @inheritParams calibration_check
#' @exportS3Method
#' @noRd
calibration_check.bgms = function(bgms_object,
                                  newdata = NULL,
                                  nrep = 200,
                                  probs = c(0.025, 0.975),
                                  grid_size = 101,
                                  seed = NULL,
                                  ...) {
  check_interval_probs(probs)
  check_positive_integer(nrep, "nrep")
  check_positive_integer(grid_size, "grid_size")

  arguments = extract_arguments(bgms_object)
  if(any(arguments$variable_type == "continuous")) {
    stop(
      "Calibration checks cover discrete variables. This fit contains ",
      "continuous variables, whose predictions are densities rather than ",
      "distributions over categories, so the isotonic curve and its ",
      "category-resampling band do not apply to them."
    )
  }

  if(is.null(newdata)) newdata = fitted_observed_data(bgms_object)
  if(anyNA(newdata)) {
    stop(
      "The data carry missing values, which have no observed category to ",
      "calibrate against. Supply complete data through 'newdata', or refit ",
      "with na_action = \"listwise\"."
    )
  }
  if(!is.null(seed)) set.seed(check_seed(seed))

  predicted = stats::predict(bgms_object,
    newdata = newdata, type = "probabilities", method = "posterior-mean"
  )
  levels_list = get_fit_spec(bgms_object)$data$category_levels
  names(levels_list) = get_fit_spec(bgms_object)$data$data_columnnames

  grid = seq(0, 1, length.out = grid_size)
  curves = list()
  summaries = list()

  for(v in names(predicted)) {
    probabilities = predicted[[v]]
    num_categories = ncol(probabilities)
    num_thresholds = num_categories - 1L
    if(num_thresholds < 1L) next

    cumulative = t(apply(probabilities, 1, cumsum))
    # The last cumulative probability is 1 for every case and carries no
    # information about calibration.
    p = as.vector(cumulative[, -num_categories, drop = FALSE])

    observed_index = match(newdata[, v], levels_list[[v]])
    y = as.integer(outer(observed_index, seq_len(num_thresholds), "<="))
    observed_curve = pav_curve(p, y, grid)

    band = vapply(seq_len(nrep), function(r) {
      # Inverse-CDF draw of one category per case, from that case's own
      # predicted distribution: nested threshold events by construction.
      u = stats::runif(nrow(probabilities))
      resampled = rowSums(u > cumulative) + 1L
      pav_curve(p, as.integer(outer(resampled, seq_len(num_thresholds), "<=")), grid)
    }, numeric(grid_size))
    bounds = apply(band, 1, stats::quantile, probs = probs)

    curves[[v]] = data.frame(
      variable = v, grid = grid, pav = observed_curve,
      lower = bounds[1, ], upper = bounds[2, ],
      row.names = NULL, stringsAsFactors = FALSE
    )
    summaries[[v]] = data.frame(
      variable = v,
      mean_dev = mean(abs(observed_curve - grid)),
      max_dev = max(abs(observed_curve - grid)),
      share_outside_band = mean(observed_curve < bounds[1, ] | observed_curve > bounds[2, ]),
      row.names = NULL, stringsAsFactors = FALSE
    )
  }

  summary_table = do.call(rbind, c(summaries, list(make.row.names = FALSE)))
  summary_table = summary_table[order(summary_table$max_dev, decreasing = TRUE), ,
    drop = FALSE
  ]

  structure(
    list(
      curves = do.call(rbind, c(curves, list(make.row.names = FALSE))),
      summary = summary_table,
      nrep = nrep,
      probs = probs,
      grid = grid
    ),
    class = "bgms_calibration"
  )
}


#' @title Print a Calibration Check
#'
#' @description
#' Prints the per-variable departure from the diagonal, worst variable first.
#'
#' @param x An object of class `bgms_calibration`, from [calibration_check()].
#' @param digits Number of digits for the printed columns. Default `3`.
#' @param max_rows Number of variables to print. Default `10`.
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#'
#' @seealso [calibration_check()]
#' @family diagnostics
#' @export
print.bgms_calibration = function(x, digits = 3, max_rows = 10L, ...) {
  cat(sprintf(
    "Calibration of the conditional predictions, %g%% consistency band from %d resamples:\n\n",
    100 * diff(x$probs), x$nrep
  ))

  body = utils::head(x$summary, max_rows)
  body[, -1] = lapply(body[, -1, drop = FALSE], round, digits)
  print(body, row.names = FALSE)
  if(nrow(x$summary) > max_rows) {
    cat(sprintf("... (%d more variables)\n", nrow(x$summary) - max_rows))
  }
  cat(
    "\nA calibrated variable tracks the diagonal and stays inside its band;",
    "\n'share_outside_band' is the share of the curve that does not.\n"
  )
  invisible(x)
}


#' @title Plot a Calibration Check
#'
#' @description
#' Reliability diagrams as small multiples, one panel per variable: the
#' isotonic fit of observed frequency on predicted probability, its consistency
#' band, and the diagonal a calibrated model tracks.
#'
#' @param x An object of class `bgms_calibration`, from [calibration_check()].
#' @param variables Optional character vector selecting which variables to
#'   draw. Defaults to all, worst departure first.
#' @param ... Ignored.
#'
#' @return `x`, invisibly. Called for the side effect of drawing.
#'
#' @examples
#' \donttest{
#' fit = bgm(Wenchuan[, 1:5], display_progress = "none")
#' plot(calibration_check(fit, nrep = 50))
#' }
#'
#' @seealso [calibration_check()]
#' @family diagnostics
#' @export
plot.bgms_calibration = function(x, variables = NULL, ...) {
  shown = if(is.null(variables)) x$summary$variable else variables
  unknown = setdiff(shown, x$summary$variable)
  if(length(unknown)) {
    stop(
      "These variables are not in the calibration check: ",
      paste(unknown, collapse = ", "), "."
    )
  }

  ink = "grey25"
  muted = "grey55"
  accent = mover_palette()[1]

  n = length(shown)
  ncol_panels = min(n, ceiling(sqrt(n)))
  nrow_panels = ceiling(n / ncol_panels)

  old_par = graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(
    mfrow = c(nrow_panels, ncol_panels),
    mar = c(2.6, 2.6, 2.0, 0.8), mgp = c(1.5, 0.5, 0), oma = c(2.2, 2.2, 0, 0),
    col.axis = ink, col.lab = ink, col.main = ink
  )

  for(v in shown) {
    df = x$curves[x$curves$variable == v, , drop = FALSE]
    graphics::plot(NA, NA,
      xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, asp = 1,
      xlab = "", ylab = "", main = v, cex.main = 0.9
    )
    graphics::axis(1,
      at = c(0, 0.5, 1), labels = c("0", ".5", "1"),
      col = muted, col.ticks = muted, cex.axis = 0.75
    )
    graphics::axis(2,
      at = c(0, 0.5, 1), labels = c("0", ".5", "1"),
      col = muted, col.ticks = muted, las = 1, cex.axis = 0.75
    )
    graphics::polygon(
      c(df$grid, rev(df$grid)), c(df$lower, rev(df$upper)),
      col = grDevices::adjustcolor(muted, 0.22), border = NA
    )
    graphics::abline(0, 1, col = muted, lwd = 0.7)
    graphics::lines(df$grid, df$pav, col = accent, lwd = 1.8)
  }

  graphics::mtext("Predicted probability", side = 1, outer = TRUE, line = 0.8, col = ink)
  graphics::mtext("Observed frequency (isotonic fit)",
    side = 2, outer = TRUE, line = 0.8, col = ink
  )
  invisible(x)
}

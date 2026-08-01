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
# A continuous variable's conditional prediction is a density, so there is no
# category to have fallen in. The same question is asked through the probability
# integral transform: u = F(y | rest) is Uniform(0, 1) exactly when the
# conditional predictive distribution is right, and the panel compares the
# empirical distribution of the u's with that uniform.
#
# Origin: the CORP reliability diagram (Dimitriadis, Gneiting & Jordan) and the
# visual predictive checks of Sailynoja et al. (arXiv 2503.01509).
# ==============================================================================


# ------------------------------------------------------------------
# fitted_observed_data
# ------------------------------------------------------------------
# The data the model saw, on the scale simulate() and predict() work in.
# Discrete columns are stored zero-based with their original levels alongside;
# a continuous block is stored centered for a GGM (with the training means
# kept) and uncentered for a mixed model.
#
# @param bgms_object  A bgms fit.
#
# Returns: a matrix of observed responses on the original scale.
# ------------------------------------------------------------------
fitted_observed_data = function(bgms_object) {
  spec = get_fit_spec(bgms_object)
  d = spec$data

  if(identical(spec$model_type, "ggm")) {
    out = sweep(d$x, 2, d$column_means, "+")
    colnames(out) = d$data_columnnames
    return(out)
  }

  # Regular ordinal columns carry a recode map, Blume-Capel columns an additive
  # shift, so the inverse of the training recode is the same one simulate()
  # uses to return data on the original scale.
  decode = function(x) {
    recode_simulated_to_original(x, d$category_levels, d$blume_capel_shift)
  }

  if(identical(spec$model_type, "mixed_mrf")) {
    out = matrix(NA_real_,
      nrow = d$num_cases, ncol = d$num_variables,
      dimnames = list(NULL, d$data_columnnames)
    )
    out[, d$discrete_indices] = decode(d$x_discrete)
    out[, d$continuous_indices] = d$x_continuous
    return(out)
  }

  out = decode(d$x)
  colnames(out) = d$data_columnnames
  out
}


# ------------------------------------------------------------------
# discrete_category_index
# ------------------------------------------------------------------
# 1-based index of each observed value among a discrete variable's categories,
# i.e. the column of the predicted-probability matrix the case fell in.
#
# @param values      Observed values on the original scale.
# @param levels_v    The variable's recode map, or NULL for Blume-Capel.
# @param shift_v     The variable's additive shift, or NA for regular ordinal.
# @param variable    Variable name, for the error message.
#
# Returns: integer vector of category indices.
# ------------------------------------------------------------------
discrete_category_index = function(values, levels_v, shift_v, variable) {
  index = if(!is.null(levels_v)) {
    match(values, levels_v)
  } else if(!is.na(shift_v)) {
    # Blume-Capel: the category score itself indexes the categories once the
    # training shift to the 0-based scale is removed.
    as.integer(round(values - shift_v)) + 1L
  } else {
    stop(
      "The fit carries neither a recode map nor a category shift for ",
      "variable '", variable, "', so its observed categories cannot be ",
      "located. Re-fit with the current bgms version."
    )
  }
  index
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


# ------------------------------------------------------------------
# conditional_pit
# ------------------------------------------------------------------
# Probability integral transform of the observed values of the continuous
# variables under their conditional predictive distributions.
#
# The predictive distribution is the mixture over posterior draws,
# u_ij = mean_d pnorm((y_ij - m_ij^(d)) / s_ij^(d)), so parameter uncertainty
# sits inside the distribution the observation is transformed by, which is what
# makes the uniform reference the right one. A single draw's Gaussian would be
# too narrow by exactly the amount the posterior is uncertain.
#
# @param bgms_object  A bgms fit with continuous variables.
# @param newdata      Data in the layout bgm() was given.
# @param cont_names   Names of the continuous variables.
# @param ndraws       Posterior draws in the mixture.
#
# Returns: named list of numeric vectors of PIT values, one per variable.
# ------------------------------------------------------------------
conditional_pit = function(bgms_object, newdata, cont_names, ndraws) {
  arguments = extract_arguments(bgms_object)
  data_columnnames = arguments$data_columnnames
  predict_vars = match(cont_names, data_columnnames)

  per_draw = if(isTRUE(arguments$is_mixed)) {
    predict_bgms_mixed(
      object = bgms_object, newdata = newdata, predict_vars = predict_vars,
      arguments = arguments, type = "probabilities",
      method = "posterior-sample", ndraws = ndraws, return_draws = TRUE
    )
  } else {
    predict_bgms_ggm(
      object = bgms_object, newdata = newdata, predict_vars = predict_vars,
      data_columnnames = data_columnnames,
      num_variables = arguments$num_variables, type = "probabilities",
      method = "posterior-sample", ndraws = ndraws, return_draws = TRUE
    )
  }

  out = vector("list", length(cont_names))
  names(out) = cont_names
  for(v in cont_names) {
    y = as.numeric(newdata[, v])
    acc = numeric(length(y))
    for(d in seq_along(per_draw)) {
      pars = per_draw[[d]][[v]]
      acc = acc + stats::pnorm(y, mean = pars[, "mean"], sd = pars[, "sd"])
    }
    out[[v]] = acc / length(per_draw)
  }
  out
}


# ------------------------------------------------------------------
# uniform_ecdf_band
# ------------------------------------------------------------------
# Pointwise quantiles of the empirical distribution function of n independent
# uniforms, read off on a grid.
#
# The probability integral transform removes the dependence on the predicted
# distribution: whatever the conditional density was, a value drawn from it
# transforms to an exact Uniform(0, 1). So the band does not resample from the
# model, as the discrete one must; it depends on n alone and is computed once
# for every continuous variable in a fit.
#
# @param n     Number of cases.
# @param nrep  Simulated uniform samples behind the band.
# @param probs The band's quantiles.
# @param grid  Grid to read the band off on.
#
# Returns: a 2-row matrix of lower and upper bounds, one column per grid point.
# ------------------------------------------------------------------
uniform_ecdf_band = function(n, nrep, probs, grid) {
  draws = vapply(seq_len(nrep), function(r) {
    u = stats::runif(n)
    vapply(grid, function(g) mean(u <= g), numeric(1))
  }, numeric(length(grid)))
  apply(draws, 1, stats::quantile, probs = probs)
}


#' @title Calibration Check
#'
#' @description
#' Reliability diagrams of the model's conditional predictions, one per
#' variable, with the consistency band a calibrated model would wander inside.
#'
#' @param bgms_object A fitted model object of class `bgms` (from [bgm()]).
#' @param newdata Optional data to evaluate the predictions on, in the layout
#'   [bgm()] was given. Defaults to the data the model was fitted to.
#' @param nrep Number of resampled datasets behind the consistency band.
#'   Default `200`.
#' @param probs Numeric of length two; the band's quantiles. Default
#'   `c(0.025, 0.975)`.
#' @param grid_size Number of points the curves are read off on. Default `101`.
#' @param seed Optional integer seed for the band resampling.
#' @param ndraws Posterior draws in the predictive mixture behind a continuous
#'   variable's panel. Ignored when every variable is discrete. Default `500`.
#' @param ... Passed to methods.
#'
#' @return An object of class `bgms_calibration`, a list with:
#'   \describe{
#'     \item{curves}{A data frame with one row per variable and grid point:
#'       `variable`, the panel `kind` (`"pav"` or `"pit"`), `grid`, the fitted
#'       `curve`, and the band bounds `lower` and `upper`.}
#'     \item{summary}{One row per variable, worst first: `kind`, `mean_dev` and
#'       `max_dev`, the mean and maximum absolute distance of the curve from the
#'       diagonal, and `share_outside_band`.}
#'     \item{nrep, probs, grid, ndraws}{The settings the check ran under.}
#'   }
#'
#' @details
#' \strong{Discrete variables} are checked against the categories they fall in.
#' For every case and every category threshold the model issues a cumulative
#' probability and the data record the outcome; isotonic regression
#' (pool-adjacent-violators) estimates the observed frequency as a monotone
#' function of the predicted one, and a calibrated variable tracks the diagonal.
#' The band resamples each case's *category* from its own predicted
#' distribution and refits the curve. It has to be built that way: the
#' cumulative threshold events of one case are nested, so resampling the events
#' independently would understate the band and make an ordinary wander look
#' like miscalibration. The predictions are the posterior-mean ones.
#'
#' \strong{Continuous variables} have a density rather than a distribution over
#' categories, and are checked through the probability integral transform
#' \eqn{u_{ij} = F_j(y_{ij} \mid y_{i,-j})}. The panel is the empirical
#' distribution function of the \eqn{u}'s against the uniform diagonal.
#' \eqn{F_j} is the predictive mixture over `ndraws` posterior draws, so
#' parameter uncertainty sits inside the distribution the observation is
#' transformed by; the two panel kinds therefore differ in what they condition
#' on, the isotonic one on posterior-mean probabilities and this one on the full
#' predictive distribution. The band is the simultaneous envelope of the
#' empirical distribution functions of \eqn{n} independent uniforms. It is not
#' resampled from the model, because it does not have to be: whatever the
#' conditional density was, a value drawn from it transforms to an exact
#' \eqn{\textrm{Uniform}(0, 1)}, so the null is known and one band serves every
#' continuous variable in the fit.
#'
#' Both panels live on the unit square with the diagonal as the calibrated
#' reference, so a mixed fit produces one figure and one summary table, with the
#' `kind` column recording which construction produced each row.
#'
#' Evaluated on the fitted data the check is in-sample, and the band is the
#' reference a model that is calibrated by construction produces on the same
#' data. In-sample the observation also entered the parameters it is judged
#' against, so the transform is mildly under-dispersed and the check is
#' conservative. Supplying held-out `newdata` removes that and makes the check
#' strictly harder to pass.
#'
#' Calibrated conditional predictions do not imply that the model reproduces
#' the joint distribution: a model can predict each variable well from the
#' others and still understate how strongly they depend on one another. That
#' second question is answered by a display built on [simulate()], not by this
#' check.
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
                             ndraws = 500,
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
                                  ndraws = 500,
                                  ...) {
  check_interval_probs(probs)
  check_positive_integer(nrep, "nrep")
  check_positive_integer(grid_size, "grid_size")
  check_positive_integer(ndraws, "ndraws")

  arguments = extract_arguments(bgms_object)
  variable_type = expand_variable_type(
    arguments$variable_type, arguments$num_variables
  )
  is_continuous = variable_type == "continuous"
  cont_names = arguments$data_columnnames[is_continuous]

  if(is.null(newdata)) newdata = fitted_observed_data(bgms_object)
  if(inherits(newdata, "data.frame")) newdata = data.matrix(newdata)
  if(ncol(newdata) != arguments$num_variables) {
    stop(
      "'newdata' must have ", arguments$num_variables,
      " columns, one per variable in the fit, but has ", ncol(newdata), "."
    )
  }
  # Columns are matched by position, as predict() matches them, so the fit's
  # own names label them here whatever the caller supplied.
  colnames(newdata) = arguments$data_columnnames
  if(anyNA(newdata)) {
    stop(
      "The data carry missing values, which have no observed value to ",
      "calibrate against. Supply complete data through 'newdata', or refit ",
      "with na_action = \"listwise\"."
    )
  }
  if(!is.null(seed)) set.seed(check_seed(seed))

  grid = seq(0, 1, length.out = grid_size)
  panels = list()

  # Both panel kinds are a curve on the unit square read against the diagonal,
  # so they share their table rows and their summary; only the construction of
  # the curve and its band differs.
  panel = function(v, kind, observed_curve, bounds) {
    list(
      curve = data.frame(
        variable = v, kind = kind, grid = grid, curve = observed_curve,
        lower = bounds[1, ], upper = bounds[2, ],
        row.names = NULL, stringsAsFactors = FALSE
      ),
      summary = data.frame(
        variable = v, kind = kind,
        mean_dev = mean(abs(observed_curve - grid)),
        max_dev = max(abs(observed_curve - grid)),
        share_outside_band = mean(
          observed_curve < bounds[1, ] | observed_curve > bounds[2, ]
        ),
        row.names = NULL, stringsAsFactors = FALSE
      )
    )
  }

  # --- Discrete variables: isotonic reliability curve -------------------------
  disc_names = arguments$data_columnnames[!is_continuous]
  if(length(disc_names)) {
    predicted = stats::predict(bgms_object,
      newdata = newdata, variables = disc_names,
      type = "probabilities", method = "posterior-mean"
    )
    # Both bookkeeping vectors are indexed over the discrete columns alone on
    # the mixed path and over all columns otherwise.
    discrete_names = if(isTRUE(arguments$is_mixed)) {
      arguments$data_columnnames_discrete
    } else {
      arguments$data_columnnames
    }
    levels_list = arguments$category_levels
    names(levels_list) = discrete_names
    shifts = arguments$blume_capel_shift
    if(is.null(shifts)) shifts = rep(NA_real_, length(discrete_names))
    names(shifts) = discrete_names

    for(v in names(predicted)) {
      probabilities = predicted[[v]]
      num_categories = ncol(probabilities)
      num_thresholds = num_categories - 1L
      if(num_thresholds < 1L) next

      cumulative = t(apply(probabilities, 1, cumsum))
      # The last cumulative probability is 1 for every case and carries no
      # information about calibration.
      p = as.vector(cumulative[, -num_categories, drop = FALSE])

      observed_index = discrete_category_index(
        newdata[, v], levels_list[[v]], shifts[[v]], v
      )
      y = as.integer(outer(observed_index, seq_len(num_thresholds), "<="))
      observed_curve = pav_curve(p, y, grid)

      band = vapply(seq_len(nrep), function(r) {
        # Inverse-CDF draw of one category per case, from that case's own
        # predicted distribution: nested threshold events by construction.
        u = stats::runif(nrow(probabilities))
        resampled = rowSums(u > cumulative) + 1L
        pav_curve(p, as.integer(outer(resampled, seq_len(num_thresholds), "<=")), grid)
      }, numeric(grid_size))
      panels[[v]] = panel(
        v, "pav", observed_curve, apply(band, 1, stats::quantile, probs = probs)
      )
    }
  }

  # --- Continuous variables: probability integral transform ------------------
  if(length(cont_names)) {
    pit = conditional_pit(bgms_object, newdata, cont_names, ndraws)
    # One band for all of them: under the transform the null is Uniform(0, 1)
    # whatever the conditional density was, so the band depends on n alone.
    bounds = uniform_ecdf_band(nrow(newdata), nrep, probs, grid)
    for(v in cont_names) {
      u = pit[[v]]
      ecdf_curve = vapply(grid, function(g) mean(u <= g), numeric(1))
      panels[[v]] = panel(v, "pit", ecdf_curve, bounds)
    }
  }

  curves = lapply(panels, `[[`, "curve")
  summaries = lapply(panels, `[[`, "summary")
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
      grid = grid,
      ndraws = if(length(cont_names)) ndraws else NA_integer_
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
  numeric_cols = vapply(body, is.numeric, logical(1))
  body[numeric_cols] = lapply(body[numeric_cols], round, digits)
  print(body, row.names = FALSE)
  if(nrow(x$summary) > max_rows) {
    cat(sprintf("... (%d more variables)\n", nrow(x$summary) - max_rows))
  }
  cat(
    "\nA calibrated variable tracks the diagonal and stays inside its band;",
    "\n'share_outside_band' is the share of the curve that does not.\n"
  )
  kinds = unique(x$summary$kind)
  if("pit" %in% kinds) {
    cat(sprintf(
      paste(
        "\n'pit' panels are continuous variables: the distribution of the",
        "probability\nintegral transform under the predictive mixture over %d",
        "posterior draws,\nagainst the uniform a calibrated density gives.\n"
      ),
      x$ndraws
    ))
  }
  if("pav" %in% kinds && "pit" %in% kinds) {
    cat(
      "'pav' panels are discrete variables: the isotonic fit of observed",
      "frequency\non posterior-mean predicted probability.\n"
    )
  }
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
    # The two panel kinds share the unit square and the diagonal but not their
    # construction, so each names its own.
    graphics::mtext(
      if(identical(df$kind[1], "pit")) "PIT" else "isotonic",
      side = 3, line = 0.05, cex = 0.6, col = muted
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
    graphics::lines(df$grid, df$curve, col = accent, lwd = 1.8)
  }

  kinds = unique(x$curves$kind[x$curves$variable %in% shown])
  xlab = if(setequal(kinds, "pit")) {
    "Probability integral transform"
  } else if(setequal(kinds, "pav")) {
    "Predicted probability"
  } else {
    "Predicted probability (isotonic) / probability integral transform (PIT)"
  }
  ylab = if(setequal(kinds, "pit")) {
    "Observed share at or below"
  } else if(setequal(kinds, "pav")) {
    "Observed frequency (isotonic fit)"
  } else {
    "Observed frequency / observed share"
  }
  graphics::mtext(xlab, side = 1, outer = TRUE, line = 0.8, col = ink, cex = 0.85)
  graphics::mtext(ylab, side = 2, outer = TRUE, line = 0.8, col = ink, cex = 0.85)
  invisible(x)
}

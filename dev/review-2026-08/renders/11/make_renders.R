# =============================================================================
# Brief 11 -- BEFORE/AFTER renders of every figure the package draws
# =============================================================================
#
# One script, run twice, once against each source tree:
#
#   Rscript make_renders.R before ~/bgms-review/wt-fix8-base
#   Rscript make_renders.R after  ~/bgms-review/wt-fix8
#
# Both runs draw the same figures from the same seeded fits at the same device
# size, so a pair of files differs only where the drawing code differs. The
# fits are cached per side under <outdir>/.cache-<side>, because an S7 fit
# object belongs to the namespace that created it and must not be carried
# across trees; deleting the cache re-fits from the same seeds.
#
# Usage: Rscript make_renders.R <side> <package path> [output directory]
# =============================================================================

args = commandArgs(trailingOnly = TRUE)
if(length(args) < 2L) {
  stop("usage: Rscript make_renders.R <before|after> <package path> [outdir]")
}
side = match.arg(args[1], c("before", "after"))
pkg = normalizePath(args[2], mustWork = TRUE)
outdir = if(length(args) >= 3L) {
  args[3]
} else {
  # Default to this script's own directory, so the renders land beside it.
  this = grep("^--file=", commandArgs(FALSE), value = TRUE)
  if(length(this)) dirname(sub("^--file=", "", this[1])) else "."
}
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

suppressMessages(pkgload::load_all(pkg, quiet = TRUE))
library(qgraph)

cache = file.path(outdir, paste0(".cache-", side))
dir.create(cache, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Fits: small, seeded, and reused across figures. Runtime is not the point;
# every figure below has to show real structure -- some edges present, some
# ruled out, some undecided, and for the comparison a real difference.
# -----------------------------------------------------------------------------
cached = function(name, expr) {
  file = file.path(cache, paste0(name, ".rds"))
  if(file.exists(file)) {
    return(readRDS(file))
  }
  value = force(expr)
  saveRDS(value, file)
  value
}

data("Wenchuan", package = "bgms")
data("Boredom", package = "bgms")

fit = cached("fit", bgm(
  Wenchuan[, 1:6],
  chains = 2, iter = 400, warmup = 400, cores = 2, seed = 1,
  display_progress = "none", verbose = FALSE
))

fit_nosel = cached("fit_nosel", bgm(
  Wenchuan[, 1:6], edge_selection = FALSE,
  chains = 2, iter = 400, warmup = 400, cores = 2, seed = 1,
  display_progress = "none", verbose = FALSE
))

# A mixed fit, so the calibration display can be shown with both panel kinds:
# PIT panels for the continuous variables, isotonic panels for the discrete.
mixed_data = local({
  set.seed(31)
  n = 250
  latent = rnorm(n)
  discrete = function(sd) round(pmin(pmax(latent + rnorm(n, sd = sd), -1.2), 1.2)) + 1
  cbind(
    d1 = discrete(0.6), c1 = latent + rnorm(n, sd = 0.5),
    d2 = discrete(0.6), c2 = latent + rnorm(n, sd = 0.5)
  )
})
fit_mixed = cached("fit_mixed", bgm(
  mixed_data,
  variable_type = c("ordinal", "continuous", "ordinal", "continuous"),
  chains = 2, iter = 300, warmup = 300, cores = 2, seed = 8,
  display_progress = "none", verbose = FALSE
))

compare = cached("compare", bgmCompare(
  x = Wenchuan[1:120, 1:5], group_indicator = rep(1:2, each = 60),
  iter = 300, warmup = 300, chains = 2, cores = 2, seed = 13,
  difference_selection = TRUE, display_progress = "none"
))

# main_difference_selection = TRUE is what gives the nodes a probability to
# carry, so the node-mark change has a figure to show it on.
compare_main = cached("compare_main", local({
  rows = c(1:60, 491:550)
  bgmCompare(
    x = Boredom[rows, 2:5], group_indicator = Boredom[rows, "language"],
    difference_selection = TRUE, main_difference_selection = TRUE,
    iter = 300, warmup = 300, chains = 2, cores = 2, seed = 44,
    display_progress = "none"
  )
}))

calibration = cached("calibration", calibration_check(fit, nrep = 60, seed = 5))
calibration_mixed = cached("calibration_mixed",
  calibration_check(fit_mixed, nrep = 60, seed = 5)
)
calibration_groups = cached("calibration_groups",
  calibration_check(compare, nrep = 60, seed = 5)
)

# Five refits each; kept small deliberately, and run one at a time.
sensitivity = cached("sensitivity", prior_sensitivity_check(
  cached("fit_small", bgm(
    Wenchuan[1:150, 1:5],
    chains = 2, iter = 300, warmup = 300, cores = 2, seed = 4,
    display_progress = "none", verbose = FALSE
  )),
  iter = 300, warmup = 300, cores = 2, seed = 4L
))
sensitivity_difference = cached("sensitivity_difference",
  prior_sensitivity_check(compare, iter = 200, warmup = 200, cores = 2,
    seed = 4L)
)

# -----------------------------------------------------------------------------
# Drawing
# -----------------------------------------------------------------------------
# Every figure is seeded immediately before it is drawn: qgraph's spring layout
# is randomised, and without this the same code would not produce the same file
# twice, which would make a BEFORE/AFTER comparison meaningless.
render = function(name, width = 7.5, height = 7, opts = NULL, expr) {
  # Options are set and restored here, in this frame. Setting them inside the
  # `expr` promise instead leaks them into every later figure, because the
  # promise is evaluated in the caller's environment and its on.exit() never
  # belongs to a frame that ends when the figure does.
  if(length(opts)) {
    old = options(opts)
    on.exit(options(old), add = TRUE)
  }
  file = file.path(outdir, sprintf("%s-%s.png", name, side))
  grDevices::png(file, width = width, height = height, units = "in", res = 140)
  on.exit({
    grDevices::dev.off()
    cat(sprintf("  %s\n", basename(file)))
  }, add = TRUE)
  set.seed(2026)
  result = tryCatch(force(expr), error = function(e) {
    # A figure this tree cannot draw is recorded as such rather than taking the
    # whole run down: that is itself a difference worth seeing in the pair.
    graphics::par(mar = c(1, 1, 1, 1))
    graphics::plot.new()
    graphics::text(0.5, 0.5, paste("not drawn:", conditionMessage(e)),
      cex = 0.8, col = "grey30"
    )
    invisible(NULL)
  })
  invisible(result)
}

cat(sprintf("rendering %s from %s\n", side, pkg))

# --- plot.bgms ---------------------------------------------------------------
render("bgms-network", expr = plot(fit))
# The variant legend keys the three lines to the Bayes factors that produce
# them instead of naming the verdicts. On the BEFORE tree the option does not
# exist and the figure is the ordinary one, which is the comparison.
render("bgms-network-legend-evidence",
  opts = list(bgms.network_legend = "evidence"), expr = plot(fit)
)
render("compare-difference-legend-evidence",
  opts = list(bgms.network_legend = "evidence"), expr = plot(compare)
)
render("bgms-network-sparse", expr = plot(fit, evidence_threshold = 1000))
render("bgms-centrality", expr = plot(fit, type = "centrality"))

# --- plot.bgmCompare ---------------------------------------------------------
render("compare-difference", expr = plot(compare))
render("compare-difference-main-selection", expr = plot(compare_main))
render("compare-difference-empty",
  expr = plot(compare, evidence_threshold = 1e6)
)
render("compare-groups", width = 15, height = 6.5,
  expr = plot(compare, type = "groups")
)
render("compare-centrality", expr = plot(compare, type = "centrality"))

# --- plot.bgms_centrality ----------------------------------------------------
render("centrality-panel", expr = plot(extract_centrality(fit)))

# --- plot.bgms_calibration ---------------------------------------------------
render("calibration-isotonic", expr = plot(calibration))
render("calibration-mixed", expr = plot(calibration_mixed))
render("calibration-groups", width = 9, height = 8,
  expr = plot(calibration_groups)
)
render("calibration-single", expr = plot(calibration, variables = "intrusion"))

# --- plot.bgms_prior_sensitivity ---------------------------------------------
render("sensitivity-edges", expr = plot(sensitivity))
render("sensitivity-differences", expr = plot(sensitivity_difference))

# --- plot_edge_posterior -----------------------------------------------------
# One edge of each kind, chosen by the evidence rather than by name, so the
# three panels really are the presence, absence and undecided cases.
table = verdicts(fit, evidence_threshold = 10)
pick = function(verdict) {
  rows = table[table$verdict == verdict, , drop = FALSE]
  if(!nrow(rows)) {
    return(NULL)
  }
  # The clearest instance of each case: the most extreme Bayes factor for
  # presence and absence, the one nearest zero for undecided.
  order_by = switch(verdict,
    presence = -rows$log_bf, absence = rows$log_bf, abs(rows$log_bf)
  )
  strsplit(rows$parameter[order(order_by)][1], "-", fixed = TRUE)[[1]]
}
for(verdict in c("presence", "undecided", "absence")) {
  ends = pick(verdict)
  render(paste0("edge-panel-", verdict), expr = {
    if(is.null(ends)) {
      stop("no edge in this fit has verdict '", verdict, "'")
    }
    plot_edge_posterior(fit, ends[1], ends[2])
  })
}
render("edge-panel-no-selection",
  expr = plot_edge_posterior(fit_nosel, "intrusion", "dreams")
)

cat("done\n")

test_that("verdict_edge_colors carries sign for present edges and greys the rest", {
  palette = mover_palette()
  colors = verdict_edge_colors(
    weight = c(0.2, -0.2, 0.2, -0.2, 0),
    verdict = c("presence", "presence", "undecided", "undecided", "presence")
  )
  expect_equal(colors[1], palette[1])
  expect_equal(colors[2], palette[2])
  expect_equal(colors[3:4], rep("grey65", 2))
  # A zero weight on a present edge is not negative.
  expect_equal(colors[5], palette[1])
})

test_that("resolve_variable accepts names and positions and rejects the rest", {
  variables = c("a", "b", "c")
  expect_equal(resolve_variable("b", variables, "variable1"), 2L)
  expect_equal(resolve_variable(3, variables, "variable2"), 3L)
  expect_error(resolve_variable("z", variables, "variable1"), "not one of the model's variables")
  expect_error(resolve_variable(0, variables, "variable1"), "column position")
  expect_error(resolve_variable(4, variables, "variable1"), "column position")
  expect_error(resolve_variable(1.5, variables, "variable1"), "column position")
})

test_that("plot.bgms draws the verdict-encoded network and the centrality panel", {
  skip_on_cran()
  skip_if_not_installed("qgraph")
  fit = get_bgms_fit_wenchuan6()

  path = withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_invisible(plot(fit))
  expect_invisible(plot(fit, legend = FALSE))
  expect_invisible(plot(fit, type = "centrality"))

  # The threshold is the one verdicts() uses, so the picture and the table move
  # together; an unreachable threshold leaves nothing to draw.
  expect_invisible(plot(fit, evidence_threshold = 3))
  expect_error(plot(fit, evidence_threshold = 1), "greater than 1")
  expect_error(plot(fit, type = "nonesuch"))
})

test_that("plot_edge_posterior draws one edge and validates its arguments", {
  skip_on_cran()
  fit = get_bgms_fit_wenchuan6()

  path = withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_invisible(plot_edge_posterior(fit, "intrusion", "upset"))
  # The edge is unordered: either order names the same edge.
  expect_invisible(plot_edge_posterior(fit, "upset", "intrusion"))
  expect_invisible(plot_edge_posterior(fit, 1, 4))

  expect_error(plot_edge_posterior(fit, "intrusion", "intrusion"), "two different variables")
  expect_error(plot_edge_posterior(fit, "intrusion", "nope"), "not one of the model's variables")
  expect_error(
    plot_edge_posterior(fit, "intrusion", "upset", evidence_threshold = 0.5),
    "greater than 1"
  )
})

test_that("the network needs at least one edge that is not ruled out", {
  skip_on_cran()
  skip_if_not_installed("qgraph")
  fit = get_bgms_fit_wenchuan6()
  edges = verdicts(fit, evidence_threshold = 10)
  skip_if(all(edges$verdict == "absence"))

  path = withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)

  # Every edge the picture draws is one verdicts() does not call absent, and
  # every edge it omits is one verdicts() does.
  drawn = edges$parameter[edges$verdict != "absence"]
  expect_gt(length(drawn), 0L)
  expect_invisible(plot(fit))
})

test_that("main_difference_nodes fills each node's wheel to its inclusion probability", {
  verdict = c("presence", "undecided", "absence")
  pip = c(0.97, 0.5, 0.02)
  nodes = main_difference_nodes(verdict, pip, main_selected = TRUE)

  # The filled fraction carries the number, so the encoding does not fail for
  # a reader who cannot tell the wheel colours apart.
  expect_equal(nodes$prob, pip)
  expect_equal(nodes$color, c(mover_palette()[1], "grey55", "grey80"))

  # A never-updated indicator has a NaN probability; its wheel stays empty
  # rather than poisoning the wedge arithmetic.
  partial = main_difference_nodes(c("presence", NA), c(0.9, NaN), main_selected = TRUE)
  expect_equal(partial$prob, c(0.9, 0))

  # Without main selection there is no indicator and no wheel at all.
  empty = main_difference_nodes(rep(NA_character_, 3), rep(NaN, 3), main_selected = FALSE)
  expect_null(empty$prob)
  expect_null(empty$color)
})

test_that("both network methods report evidence through the same band", {
  verdict = c("presence", "presence", "undecided", "absence", "absence")
  log_bf = c(4.2, Inf, 0.3, -3.1, -9.4)

  edges = evidence_band(verdict, log_bf, 10, network_unit("edge"))
  differences = evidence_band(verdict, log_bf, 10, network_unit("difference"))

  # The tally counts the same three classes on both sides; only the noun for
  # the thing being counted differs.
  expect_equal(edges$left, "2 present  |  1 undecided  |  2 ruled out")
  expect_equal(differences$left, "2 differing  |  1 undecided  |  2 ruled out")

  # The numbers are identical, in the same wording, on the same natural-log
  # scale, through the same formatters: this is what parity between the two
  # displays means in code.
  expect_equal(edges$right, differences$right)
  expect_equal(edges$right[1], paste("threshold: log BF", format_log_bf(log(10))))
  # A saturated Bayes factor prints as the reporting cap, exactly as it does
  # on an edge panel.
  expect_equal(edges$right[2], "strongest for: log BF > 10,000")
  expect_equal(edges$right[3], "strongest against: log BF = -9.4")

  # An all-NA column of Bayes factors has no extremes to print, and says so by
  # printing only the threshold rather than an invented range.
  quiet = evidence_band(verdict, rep(NA_real_, 5), 10, network_unit("edge"))
  expect_length(quiet$right, 1L)
})

test_that("the legend can key its lines to verdicts or to the Bayes factors", {
  unit = network_unit("difference")
  expect_equal(network_legend_keys(unit, 10), unit$keys)

  withr::local_options(bgms.network_legend = "evidence")
  keys = network_legend_keys(unit, 10)
  # The variant says nothing a threshold has not already been given for: the
  # same three lines, named by the evidence that produces them.
  expect_length(keys, 3L)
  expect_true(all(grepl("log BF", keys, fixed = TRUE)))
  expect_false(any(grepl("undecided", keys, fixed = TRUE)))
})

test_that("compare_difference_verdicts splits the two indicator families", {
  skip_on_cran()
  fit = get_bgmcompare_fit_wenchuan5()
  found = compare_difference_verdicts(fit, evidence_threshold = 10)

  expect_equal(length(found$pairwise), 10L)
  expect_equal(length(found$main), 5L)
  expect_equal(nrow(found$pairs), 10L)
  # main_difference_selection is FALSE by default, so those indicators were
  # never updated and the node channel has nothing to say.
  expect_false(found$main_selected)
  expect_true(all(is.na(found$main)))

  # The split follows the indicator names, not an assumed ordering.
  names_all = get_raw_samples(fit)$parameter_names$indicator
  table = verdicts(fit, evidence_threshold = 10)
  expect_equal(
    found$pairwise,
    as.character(table$verdict[grepl("(pairwise)", names_all, fixed = TRUE)])
  )
})

test_that("plot.bgmCompare draws the difference network and the group panels", {
  skip_on_cran()
  skip_if_not_installed("qgraph")
  fit = get_bgmcompare_fit_wenchuan5()

  path = withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_invisible(plot(fit))
  expect_invisible(plot(fit, type = "groups"))
  expect_invisible(plot(fit, type = "centrality"))
  # The difference-centrality display is refused (F-063); its numbers remain
  # available through summary(extract_centrality()).
  expect_error(plot(fit, type = "centrality", group = c(1, 2)), "not offered")

  # An empty difference network is a result, not an error: at an unreachable
  # threshold every difference is ruled out and the nodes are drawn alone.
  expect_invisible(plot(fit, evidence_threshold = 1e6))
})

test_that("verdict_network_input aligns the encoding with qgraph's read order", {
  V = 5L
  pairs = which(upper.tri(matrix(0, V, V)), arr.ind = TRUE)
  pairs = pairs[order(pairs[, "row"], pairs[, "col"]), ]

  verdict = rep("absence", 10L)
  verdict[c(1L, 4L, 7L)] = c("presence", "presence", "undecided")
  weight = numeric(10L)
  weight[c(1L, 4L, 7L)] = c(0.4, -0.3, 0.05)

  network = verdict_network_input(weight, verdict, pairs, V)

  # A square matrix states the node set in its dimensions, so a network that
  # leaves a node out still draws that node.
  expect_equal(dim(network$weights), c(V, V))
  expect_true(isSymmetric(network$weights))
  expect_equal(network$weights[1, 2], 0.4)
  expect_equal(network$weights[1, 5], -0.3)
  expect_equal(sum(network$weights != 0), 6L)

  # qgraph reads the non-zero upper triangle in column-major order, which for
  # these three edges is (1,2), (1,5), (2,5) -- not the row-major order the
  # verdict table is in. The colour and line type must follow that read order,
  # or the encoding lands on the wrong edges.
  palette = mover_palette()
  expect_equal(unname(network$edge.color), c(palette[1], palette[2], "grey65"))
  expect_equal(network$lty, c(1L, 1L, 3L))

  # Nothing drawable is a matrix of zeros, and the caller draws the nodes.
  empty = verdict_network_input(weight, rep("absence", 10L), pairs, V)
  expect_true(all(empty$weights == 0))
  expect_length(empty$edge.color, 0L)

  # An exactly zero weight on a drawn edge would be read as a non-edge and
  # shift every later colour onto the wrong one.
  zero = verdict_network_input(
    c(0, rep(0, 9)), c("presence", rep("absence", 9L)), pairs, V
  )
  expect_equal(sum(zero$weights != 0), 2L)
  expect_length(zero$edge.color, 1L)
})

test_that("a sparse network draws rather than failing on the node set", {
  skip_on_cran()
  skip_if_not_installed("qgraph")
  fit = get_bgms_fit_wenchuan6()
  path = withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)

  # A threshold that leaves only a few edges drawable is the case that used to
  # fail: qgraph derives the node set from a weighted edgelist, and a network
  # missing a node took the whole call down.
  for(threshold in c(10, 100, 1000)) {
    expect_invisible(plot(fit, evidence_threshold = threshold))
  }
})

test_that("naming one variable twice says which one and how to fix it", {
  skip_on_cran()
  fit = get_bgms_fit_wenchuan6()
  expect_error(
    plot_edge_posterior(fit, 1, 1),
    "both resolve to 'intrusion'"
  )
  expect_error(plot_edge_posterior(fit, 1, 1), "plot_edge_posterior\\(fit, ")
})

test_that("a mixed network draws its weights at the pairs they belong to", {
  skip_on_cran()
  skip_if_not_installed("qgraph")
  set.seed(31)
  n = 200
  latent = rnorm(n)
  x = cbind(
    round(pmin(pmax(latent + rnorm(n, sd = 0.6), -1.2), 1.2)) + 1,
    latent + rnorm(n, sd = 0.5),
    round(pmin(pmax(latent + rnorm(n, sd = 0.6), -1.2), 1.2)) + 1,
    latent + rnorm(n, sd = 0.5),
    latent + rnorm(n, sd = 0.5)
  )
  colnames(x) = c("d1", "c1", "d2", "c2", "c3")
  fit = bgm(x,
    variable_type = c(
      "ordinal", "continuous", "ordinal", "continuous",
      "continuous"
    ),
    chains = 2, iter = 300, warmup = 300, cores = 2, seed = 8,
    display_progress = "none", verbose = FALSE
  )

  # The drawn weights come from the pairwise draws, which are in the fit's own
  # indicator order; the pair index has to follow it or every weight lands on
  # the wrong edge.
  weight = colMeans(extract_pairwise_interactions(fit))
  pairs = indicator_pair_index(fit, 5L)
  names_out = colnames(x)
  labels = paste(names_out[pairs[, 1]], names_out[pairs[, 2]], sep = "-")
  flipped = paste(names_out[pairs[, 2]], names_out[pairs[, 1]], sep = "-")
  expect_true(all(labels == names(weight) | flipped == names(weight)))

  path = withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(fit, type = "network", evidence_threshold = 3))
})


# ==============================================================================
# The JASP-style edge posterior panel (F-066)
# ==============================================================================

# What a panel says, as text: the wording and the structure, with the numbers
# at a precision a rerun reproduces. Snapshotting this rather than a raster
# keeps the reviewable content of the figure under test without pinning pixels.
describe_panel = function(panel) {
  present = function(part) if(is.null(part)) "no" else "yes"
  cat("label     : ", panel$label, "\n", sep = "")
  cat("subtitle  : ", panel$subtitle %||% "(none)", "\n", sep = "")
  cat("evidence  : ", paste(panel$evidence, collapse = " | "), "\n", sep = "")
  cat("estimate  : ", paste(panel$estimate %||% "(none)", collapse = " | "),
    "\n", sep = "")
  cat("wheel     : ", sprintf("%.3f", panel$wheel_prob), "\n", sep = "")
  cat("wheel tags: ", paste(panel$wheel_labels %||% "(none)", collapse = " / "),
    "\n", sep = "")
  cat("posterior : ", present(panel$posterior), "\n", sep = "")
  cat("prior     : ", present(panel$prior), "\n", sep = "")
  cat("dots      : ", if(is.null(panel$dots)) {
    "none"
  } else {
    paste(sprintf("%.3f", panel$dots$y), collapse = ", ")
  }, "\n", sep = "")
  cat("window    : ", paste(sprintf("%.2f", edge_panel_window(panel)),
    collapse = " to "), "\n", sep = "")
  cat("caption   : ", panel$caption, "\n", sep = "")
  invisible(NULL)
}

# One prior, fixed, so the snapshots test the panel and not the fit.
test_slab_prior = function(scale = 1) {
  list(
    density = function(w) stats::dnorm(w, 0, scale),
    quantile = function(p) stats::qnorm(p, 0, scale),
    family = "normal", scale = scale
  )
}

test_that("the panel of a decisive edge carries the wheel, not a stem", {
  set.seed(101)
  draws = c(rnorm(1990, 0.32, 0.04), rep(0, 10))
  panel = edge_panel_selection(
    "intrusion-dreams", draws, test_slab_prior(),
    pip = 0.995, log_bf = 12.4, verdict = "presence"
  )
  # The evidence is the Rao-Blackwellized indicator Bayes factor, so there are
  # no Savage-Dickey ordinates to mark.
  expect_null(panel$dots)
  expect_null(panel$wheel_labels)
  expect_equal(panel$wheel_prob, 0.995)
  expect_snapshot(describe_panel(panel))
})

test_that("the panel of an undecided edge splits its wheel", {
  set.seed(102)
  draws = c(rnorm(1200, 0.09, 0.03), rep(0, 800))
  panel = edge_panel_selection(
    "intrusion-upset", draws, test_slab_prior(),
    pip = 0.6, log_bf = 0.41, verdict = "undecided"
  )
  # "undecided" is not evidence of anything, so the panel does not say it is.
  expect_equal(panel$subtitle, "undecided")
  expect_snapshot(describe_panel(panel))
})

test_that("a saturated edge prints the capped Bayes factor and a full wheel", {
  set.seed(103)
  draws = rnorm(2000, 0.41, 0.03)
  panel = edge_panel_selection(
    "upset-physior", draws, test_slab_prior(),
    pip = 1, log_bf = Inf, verdict = "presence"
  )
  expect_equal(panel$evidence[1], "PIP > .99")
  expect_equal(panel$evidence[2], "log BF > 10,000")
  # The estimate is of the conditional posterior and is unaffected by the cap.
  expect_snapshot(describe_panel(panel))
})

test_that("a decisive absence with no included draw is a figure, not an error", {
  panel = edge_panel_selection(
    "a-b", rep(0, 2000), test_slab_prior(),
    pip = 0.0004, log_bf = -7.2, verdict = "absence"
  )
  # No conditional posterior exists, so the prior and the wheel carry the
  # panel; nothing is claimed about a weight that was never sampled.
  expect_null(panel$posterior)
  expect_null(panel$estimate)
  expect_null(panel$interval)
  expect_false(is.null(panel$prior))
  expect_match(panel$caption, "No retained draw included this edge")
  expect_snapshot(describe_panel(panel))
})

test_that("without edge selection the panel is the Savage-Dickey figure", {
  set.seed(104)
  # A posterior that has left zero: the ordinate at zero is near nothing and
  # the Bayes factor for the edge is decisive.
  decisive = edge_panel_savage_dickey(
    "intrusion-dreams", rnorm(4000, 0.32, 0.04), test_slab_prior()
  )
  expect_equal(decisive$wheel_labels, c("data|H1", "data|H0"))
  expect_false(is.null(decisive$dots))
  # The prior ordinate is exact; the first dot is dnorm(0, 0, 1).
  expect_equal(decisive$dots$y[1], stats::dnorm(0), tolerance = 1e-12)
  expect_snapshot(describe_panel(decisive))

  # A posterior piled on zero: the posterior ordinate exceeds the prior one,
  # the log Bayes factor is negative, and the wheel is mostly pale.
  absent = edge_panel_savage_dickey(
    "intrusion-avoidth", rnorm(4000, 0, 0.03), test_slab_prior()
  )
  expect_lt(absent$wheel_prob, 0.05)
  expect_gt(absent$dots$y[2], absent$dots$y[1])
  expect_snapshot(describe_panel(absent))
})

test_that("the panel reads a Blume-Capel fit like any other", {
  skip_on_cran()
  data("Wenchuan", package = "bgms")
  fit = bgm(Wenchuan[, 1:4],
    variable_type = "blume-capel", baseline_category = 2,
    chains = 2, iter = 300, warmup = 300, cores = 2, seed = 3,
    display_progress = "none", verbose = FALSE
  )
  label = "intrusion-dreams"
  evidence = edge_selection_evidence(fit, label, 10)
  panel = edge_panel_selection(
    label, extract_pairwise_interactions(fit)[, label],
    edge_slab_prior(fit), evidence$pip, evidence$log_bf, evidence$verdict
  )
  # Weights are continuous whatever the variable type, so the panel is the
  # ordinary one; what is asserted here is that nothing about the parameter
  # layout of a Blume-Capel fit reaches the figure.
  expect_false(is.null(panel$posterior))
  expect_null(panel$dots)
  expect_equal(panel$prior$family, "normal")
  expect_snapshot({
    cat("subtitle  : ", panel$subtitle, "\n", sep = "")
    cat("wheel tags: ", panel$wheel_labels %||% "(none)", "\n", sep = "")
    cat("caption   : ", panel$caption, "\n", sep = "")
  })

  path = withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot_edge_posterior(fit, "intrusion", "dreams"))
})

test_that("the slab prior is the fit's own, not the package default", {
  skip_on_cran()
  data("Wenchuan", package = "bgms")
  cauchy = bgm(Wenchuan[1:150, 1:4],
    interaction_prior = cauchy_prior(scale = 2.5),
    chains = 2, iter = 200, warmup = 200, cores = 2, seed = 6,
    display_progress = "none", verbose = FALSE
  )
  prior = edge_slab_prior(cauchy)
  expect_equal(prior$family, "cauchy")
  expect_equal(prior$scale, 2.5)
  # The default changed in 0.2.0, so a hardcoded Normal(0, 1) would be wrong
  # here by a factor of three at the origin.
  expect_equal(prior$density(0), stats::dcauchy(0, 0, 2.5))
  expect_equal(prior$density(0.4), stats::dcauchy(0.4, 0, 2.5))

  # extract_pairwise_interactions() reports every model type in the frame the
  # slab applies to, which is what makes one closed form serve them all: it is
  # the same theta the anchored sensitivity curve reweights on.
  expect_equal(
    unname(extract_pairwise_interactions(cauchy)[, 1]),
    unname(do.call(rbind, anchor_draws(cauchy)$theta)[, 1])
  )
})

test_that("a beta-prime slab is drawn through its logistic Jacobian", {
  skip_on_cran()
  data("Wenchuan", package = "bgms")
  fit = bgm(Wenchuan[1:150, 1:4],
    interaction_prior = beta_prime_prior(alpha = 0.5, beta = 0.5),
    chains = 2, iter = 200, warmup = 200, cores = 2, seed = 6,
    display_progress = "none", verbose = FALSE
  )
  prior = edge_slab_prior(fit)
  expect_equal(prior$family, "beta-prime")
  at = c(-0.8, 0, 0.35)
  p = stats::plogis(at)
  expect_equal(prior$density(at), stats::dbeta(p, 0.5, 0.5) * p * (1 - p))
  # It is a density: it integrates to one. The tails are evaluated rather than
  # returning NaN where the Beta density overflows and the Jacobian underflows.
  expect_equal(stats::integrate(prior$density, -60, 60)$value, 1,
    tolerance = 1e-4
  )
  expect_equal(prior$density(c(-800, 800)), c(0, 0))
})

test_that("a fit without edge selection takes the Savage-Dickey branch", {
  skip_on_cran()
  fit = get_bgms_fit_wenchuan6_noselection()

  path = withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)

  # verdicts() refuses such a fit, so the panel must not be reaching for it.
  expect_error(verdicts(fit), "require edge selection")
  expect_invisible(plot_edge_posterior(fit, "intrusion", "dreams"))

  panel = edge_panel_savage_dickey(
    "intrusion-dreams",
    extract_pairwise_interactions(fit)[, "intrusion-dreams"],
    edge_slab_prior(fit)
  )
  expect_equal(panel$wheel_labels, c("data|H1", "data|H0"))
  expect_length(panel$dots$y, 2L)
})

test_that("binwidth is deprecated rather than silently dropped", {
  skip_on_cran()
  fit = get_bgms_fit_wenchuan6()

  path = withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)

  withr::local_options(lifecycle_verbosity = "warning")
  warning = tryCatch(
    plot_edge_posterior(fit, "intrusion", "dreams", binwidth = 0.02),
    warning = function(w) w
  )
  expect_s3_class(warning, "lifecycle_warning_deprecated")
  expect_match(conditionMessage(warning), "binwidth")
  # The message says what replaced it rather than only that it went away.
  expect_match(conditionMessage(warning), "wheel")
  expect_match(conditionMessage(warning), "ignored")

  # Without it there is no warning at all.
  expect_no_warning(plot_edge_posterior(fit, "intrusion", "dreams"))
})

test_that("a compare fit is refused rather than drawn as one network's edge", {
  skip_on_cran()
  fit = get_bgmcompare_fit_wenchuan5()
  expect_error(
    plot_edge_posterior(fit, 1, 2),
    "parameterizes differences between networks"
  )
})

test_that("an edge is found in whichever order its fit names it", {
  columns = c("d1-d2", "c1-c2", "d1-c1", "d2-c1")
  # A mixed fit lays its pairwise draws out by block and names a cross edge by
  # its discrete end, which need not be the earlier column.
  expect_equal(edge_column_label("c1", "d2", columns), "d2-c1")
  expect_equal(edge_column_label("d2", "c1", columns), "d2-c1")
  expect_equal(edge_column_label("d1", "d2", columns), "d1-d2")
  expect_null(edge_column_label("d1", "nope", columns))
})

test_that("a mixed fit's cross edge is drawn rather than reported missing", {
  skip_on_cran()
  set.seed(31)
  n = 250
  latent = rnorm(n)
  x = cbind(
    d1 = round(pmin(pmax(latent + rnorm(n, sd = 0.6), -1.2), 1.2)) + 1,
    c1 = latent + rnorm(n, sd = 0.5),
    d2 = round(pmin(pmax(latent + rnorm(n, sd = 0.6), -1.2), 1.2)) + 1,
    c2 = latent + rnorm(n, sd = 0.5)
  )
  fit = bgm(x,
    variable_type = c("ordinal", "continuous", "ordinal", "continuous"),
    chains = 2, iter = 300, warmup = 300, cores = 2, seed = 8,
    display_progress = "none", verbose = FALSE
  )
  path = withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)

  # The stored name is "d2-c1"; asking for it in variable order used to fail.
  expect_true("d2-c1" %in% colnames(extract_pairwise_interactions(fit)))
  expect_invisible(plot_edge_posterior(fit, "c1", "d2"))
  expect_invisible(plot_edge_posterior(fit, "c1", "c2"))
})

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
})

test_that("a network with nothing left to support is still a figure", {
  skip_on_cran()
  skip_if_not_installed("qgraph")
  fit = get_bgms_fit_wenchuan6()

  path = withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)

  # At an unreachable threshold every edge lands in the absence panel. That is
  # the honest picture of such a fit, not an error: the second panel fills and
  # the first is empty.
  expect_invisible(plot(fit, evidence_threshold = 1e6))
  # And at a threshold nothing is ruled out at, the absence panel is the empty
  # one instead.
  expect_invisible(plot(fit, evidence_threshold = 1000))
})

test_that("main_difference_nodes fills each node's ring to its inclusion probability", {
  style = bgms_style()
  verdict = c("presence", "undecided", "absence")
  pip = c(0.97, 0.5, 0.02)
  nodes = main_difference_nodes(verdict, pip, main_selected = TRUE)

  # The ring is qgraph's pie channel, drawn around the node circle where a
  # reader of a network looks for a node's own quantity. The filled fraction
  # carries the number, so the encoding does not fail for a reader who cannot
  # tell the ring colours apart.
  expect_equal(nodes$pie, pip)
  expect_equal(nodes$pie_color,
    c(mover_palette()[1], style$muted, style$pale)
  )

  # A never-updated indicator has a NaN probability; its ring stays empty
  # rather than poisoning qgraph's arc arithmetic.
  partial = main_difference_nodes(c("presence", NA), c(0.9, NaN), main_selected = TRUE)
  expect_equal(partial$pie, c(0.9, 0))

  # Without main selection there is no indicator and no ring at all.
  empty = main_difference_nodes(rep(NA_character_, 3), rep(NaN, 3), main_selected = FALSE)
  expect_null(empty$pie)
  expect_null(empty$pie_color)
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

test_that("a panel matrix states its node set and reads its colours in order", {
  V = 5L
  pairs = which(upper.tri(matrix(0, V, V)), arr.ind = TRUE)
  pairs = pairs[order(pairs[, "row"], pairs[, "col"]), ]

  keep = rep(FALSE, 10L)
  keep[c(1L, 4L, 7L)] = TRUE
  weight = numeric(10L)
  weight[c(1L, 4L, 7L)] = c(0.4, -0.3, 0.05)

  m = panel_edge_matrix(weight, keep, pairs, V)

  # A square matrix states the node set in its dimensions, so a panel that
  # leaves a node out still draws that node.
  expect_equal(dim(m), c(V, V))
  expect_true(isSymmetric(m))
  expect_equal(m[1, 2], 0.4)
  expect_equal(m[1, 5], -0.3)
  expect_equal(sum(m != 0), 6L)

  # qgraph reads the non-zero upper triangle in column-major order, which for
  # these three edges is (1,2), (1,5), (2,5) -- not the row-major order the
  # verdict table is in. Reading the signs off the matrix is what makes the
  # colours land on the right edges whatever order the table was in.
  palette = mover_palette()
  expect_equal(matrix_edge_colors(m), c(palette[1], palette[2], palette[1]))

  # An empty panel is a matrix of zeros, and qgraph draws the nodes alone --
  # which is what an all-decided fit's undecided panel has to look like.
  empty = panel_edge_matrix(weight, rep(FALSE, 10L), pairs, V)
  expect_true(all(empty == 0))
  expect_length(matrix_edge_colors(empty), 0L)

  # An exactly zero weight on a drawn edge would be read as a non-edge and
  # shift every later colour onto the wrong one.
  zero = panel_edge_matrix(numeric(10L), c(TRUE, rep(FALSE, 9L)), pairs, V)
  expect_equal(sum(zero != 0), 2L)
  expect_length(matrix_edge_colors(zero), 1L)

  # A panel drawn at uniform width is the same construction with every value
  # set to one.
  uniform = panel_edge_matrix(rep(1, 10L), keep, pairs, V)
  expect_equal(sort(unique(uniform[uniform != 0])), 1)
})

test_that("the threshold rules read in Bayes factors, not their logarithms", {
  rules = threshold_rules(10)
  expect_equal(unname(rules[["presence"]]), "BF > 10")
  expect_equal(unname(rules[["absence"]]), "BF < 0.1")
  expect_equal(unname(rules[["undecided"]]), "0.1 < BF < 10")

  # A threshold a user chose is printed as they chose it.
  expect_equal(unname(threshold_rules(100)[["presence"]]), "BF > 100")
  expect_equal(unname(threshold_rules(100)[["absence"]]), "BF < 0.01")
  expect_equal(unname(threshold_rules(3)[["absence"]]), "BF < 0.33")
})

test_that("both network methods split their pairs with the same wording rule", {
  edge = network_unit("edge")
  difference = network_unit("difference")
  # The two displays are one routine; only the nouns differ, and both name the
  # same three classes in the same order.
  expect_equal(names(edge$panels), names(difference$panels))
  expect_equal(names(edge$panels), c("presence", "absence", "undecided"))
  expect_equal(unname(edge$panels[["presence"]]), "evidence of presence")
  expect_equal(unname(difference$panels[["presence"]]), "difference supported")
  expect_equal(edge$weights_label, "Edge weights")
  expect_equal(difference$weights_label, "Difference weights")
})

test_that("a pair that differs in several contrasts is laid out by its largest", {
  pairs = cbind(row = c(1L, 1L, 2L), col = c(2L, 3L, 3L))
  first = matrix(0, 3, 3)
  first[1, 2] = 0.4; first[1, 3] = -0.9; first[2, 3] = 0.1
  first = first + t(first)
  second = matrix(0, 3, 3)
  second[1, 2] = -0.7; second[1, 3] = 0.2; second[2, 3] = 0.05
  second = second + t(second)

  # The layout summary is the largest absolute difference over the contrasts,
  # so a pair that differs sharply in one of them is not averaged away.
  expect_equal(contrast_magnitude(list(first, second), pairs), c(0.7, 0.9, 0.1))

  # One contrast is the two-group case, and there it is just the magnitude.
  expect_equal(contrast_magnitude(list(first), pairs), c(0.4, 0.9, 0.1))
})

test_that("more than two groups keeps the panels and drops the width channel", {
  skip_on_cran()
  skip_if_not_installed("qgraph")
  data("ADHD", package = "bgms")
  fit = bgmCompare(
    x = ADHD[, 2:5],
    group_indicator = rep(1:3, length.out = nrow(ADHD)),
    iter = 50, warmup = 100, chains = 2, seed = 903,
    display_progress = "none"
  )
  # The premise of the display: one indicator per pair, several magnitudes.
  expect_true(is.list(fit@posterior_mean_pairwise_differences))
  expect_equal(length(fit@posterior_mean_pairwise_differences), 2L)

  path = withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)

  # The three-way split is defined for any K, so the figure draws.
  expect_invisible(plot(fit))
  expect_invisible(plot(fit, type = "groups"))

  # Without selection there is no split, and beyond two groups no single
  # magnitude either -- so the weighted display of such a fit is the groups
  # themselves, and plot() draws them rather than refusing.
  nosel = bgmCompare(
    x = ADHD[, 2:5],
    group_indicator = rep(1:3, length.out = nrow(ADHD)),
    difference_selection = FALSE,
    iter = 50, warmup = 100, chains = 2, seed = 904,
    display_progress = "none"
  )
  expect_invisible(plot(nosel))

  # It is the type = "groups" display, reached without asking for it: the two
  # routes behave alike, including the paging arguments, which is what makes
  # the default an alias for that display rather than a lookalike of it.
  withr::local_options(bgms.verbose = TRUE)
  expect_message(plot(nosel, max_panels = 2L), "Showing page 1 of 2")
  expect_message(plot(nosel, type = "groups", max_panels = 2L),
    "Showing page 1 of 2")
  expect_message(plot(nosel, max_panels = 2L, page = 2L), "Showing page 2 of 2")
  expect_error(plot(nosel, max_panels = 2L, page = 3L), "make 2 pages")
  expect_error(plot(nosel, max_panels = 0L), "max_panels")

  # Two groups keep the difference-weights network: there the pair has one
  # magnitude, so there is a weighted network to draw.
  two = bgmCompare(
    x = ADHD[, 2:5],
    group_indicator = rep(1:2, length.out = nrow(ADHD)),
    difference_selection = FALSE,
    iter = 50, warmup = 150, chains = 2, seed = 905,
    display_progress = "none"
  )
  expect_false(is.list(two@posterior_mean_pairwise_differences))
  expect_invisible(plot(two))
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
    pip = 0.995, log_bf = 12.4
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
    pip = 0.6, log_bf = 0.41
  )
  # No verdict word is printed, here least of all: the split wheel and the log
  # Bayes factor near zero are the statement, and "undecided" would have added
  # a threshold the reader did not choose.
  expect_null(panel$subtitle)
  expect_snapshot(describe_panel(panel))
})

test_that("a saturated edge prints the capped Bayes factor and a full wheel", {
  set.seed(103)
  draws = rnorm(2000, 0.41, 0.03)
  panel = edge_panel_selection(
    "upset-physior", draws, test_slab_prior(),
    pip = 1, log_bf = Inf
  )
  expect_equal(panel$evidence[1], "P(included) > .999")
  expect_equal(panel$evidence[2], "log BF > 10,000")
  # The estimate is of the conditional posterior and is unaffected by the cap.
  expect_snapshot(describe_panel(panel))
})

test_that("a decisive absence with no included draw is a figure, not an error", {
  panel = edge_panel_selection(
    "a-b", rep(0, 2000), test_slab_prior(),
    pip = 0.0004, log_bf = -7.2
  )
  # No conditional posterior exists, so the prior and the wheel carry the
  # panel; nothing is claimed about a weight that was never sampled.
  expect_null(panel$posterior)
  expect_null(panel$estimate)
  expect_null(panel$interval)
  expect_false(is.null(panel$prior))
  # No panel carries a caption: what each mark means is stated in the Rd, not
  # reprinted under every figure.
  expect_null(panel$caption)
  expect_snapshot(describe_panel(panel))
})

test_that("without edge selection the panel is the Savage-Dickey figure", {
  set.seed(104)
  # A posterior that has left zero: the ordinate at zero is near nothing and
  # the Bayes factor for the edge is decisive.
  decisive = edge_panel_savage_dickey(
    "intrusion-dreams", rnorm(4000, 0.32, 0.04), test_slab_prior()
  )
  # The wheel carries no notation and the panel carries no caption: what its
  # two shares are is said in the Rd.
  expect_null(decisive$wheel_labels)
  expect_null(decisive$caption)
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
  evidence = edge_selection_evidence(fit, label)
  panel = edge_panel_selection(
    label, extract_pairwise_interactions(fit)[, label],
    edge_slab_prior(fit), evidence$pip, evidence$log_bf
  )
  # Weights are continuous whatever the variable type, so the panel is the
  # ordinary one; what is asserted here is that nothing about the parameter
  # layout of a Blume-Capel fit reaches the figure.
  expect_false(is.null(panel$posterior))
  expect_null(panel$dots)
  expect_equal(panel$prior$family, "normal")
  expect_snapshot({
    cat("subtitle  : ", panel$subtitle %||% "(none)", "\n", sep = "")
    cat("wheel tags: ", panel$wheel_labels %||% "(none)", "\n", sep = "")
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
  expect_null(panel$wheel_labels)
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

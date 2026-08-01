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

test_that("main_difference_nodes fills each node's ring to its inclusion probability", {
  verdict = c("presence", "undecided", "absence")
  pip = c(0.97, 0.5, 0.02)
  nodes = main_difference_nodes(verdict, pip, main_selected = TRUE)

  # The ring fraction carries the number, so the encoding does not fail for a
  # reader who cannot tell the ring colours apart.
  expect_equal(nodes$pie, pip)
  expect_equal(nodes$pie_color, c(mover_palette()[1], "grey55", "grey80"))

  # A never-updated indicator has a NaN probability; its ring stays empty
  # rather than poisoning qgraph's arc arithmetic.
  partial = main_difference_nodes(c("presence", NA), c(0.9, NaN), main_selected = TRUE)
  expect_equal(partial$pie, c(0.9, 0))

  # Without main selection there is no indicator and no ring channel at all.
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

test_that("the edge panel title reports a capped natural log Bayes factor", {
  expect_equal(
    edge_panel_title("a-b", "presence", 301.44),
    "a-b\npresence, log BF = 301.4"
  )
  expect_equal(
    edge_panel_title("a-b", "absence", -2.34),
    "a-b\nabsence, log BF = -2.3"
  )
  # A saturated edge: exp(1e5) has no printable value, and exp(301) would put
  # 131 digits in the title.
  expect_equal(
    edge_panel_title("a-b", "presence", Inf),
    "a-b\npresence, log BF > 10,000"
  )
  expect_equal(
    edge_panel_title("a-b", "absence", -Inf),
    "a-b\nabsence, log BF < -10,000"
  )
  expect_snapshot(cat(edge_panel_title("intrusion-dreams", "presence", Inf)))
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

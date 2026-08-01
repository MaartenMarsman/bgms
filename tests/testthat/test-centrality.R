test_that("extract_centrality returns draws x nodes strength", {
  skip_on_cran()
  fit = get_bgms_fit_wenchuan6()
  samples = extract_pairwise_interactions(fit)
  strength = extract_centrality(fit)

  expect_s3_class(strength, "bgms_centrality")
  expect_equal(dim(strength), c(nrow(samples), 6L))
  expect_identical(attr(strength, "measure"), "strength")

  # Node order follows the data's columns, and the incidence comes from the
  # edge order rather than from splitting the "A-B" column names.
  expect_equal(colnames(strength), colnames(Wenchuan)[1:6])

  pairs = do.call(rbind, strsplit(colnames(samples), "-", fixed = TRUE))
  for(node in colnames(strength)) {
    incident = pairs[, 1] == node | pairs[, 2] == node
    expect_equal(
      unname(strength[, node]),
      unname(rowSums(abs(samples[, incident, drop = FALSE]))),
      tolerance = 1e-12
    )
  }

  # Every edge weight enters two nodes' strengths, so the row totals are twice
  # the summed absolute weights of the network in that draw.
  expect_equal(unname(rowSums(strength)), unname(2 * rowSums(abs(samples))),
    tolerance = 1e-10
  )

  expect_error(extract_centrality(fit, measure = "betweenness"))
})

test_that("centrality is model-averaged, so excluded edges contribute zero", {
  skip_on_cran()
  fit = get_bgms_fit_wenchuan6()
  samples = extract_pairwise_interactions(fit)
  strength = extract_centrality(fit)

  # The draws carry exact zeros for excluded edges; conditioning on inclusion
  # instead would give a strictly larger centrality for any node with an
  # uncertain edge, so the two must not coincide here.
  expect_true(any(samples == 0))

  node = colnames(strength)[1]
  pairs = do.call(rbind, strsplit(colnames(samples), "-", fixed = TRUE))
  incident = which(pairs[, 1] == node | pairs[, 2] == node)
  conditional = sum(vapply(incident, function(k) {
    w = samples[, k]
    mean(abs(w[w != 0]))
  }, numeric(1)))
  expect_lt(mean(strength[, node]), conditional)
})

test_that("summary.bgms_centrality orders by mean and normalises p_most_central", {
  skip_on_cran()
  fit = get_bgms_fit_wenchuan6()
  strength = extract_centrality(fit)
  summ = summary(strength)

  expect_named(summ, c("node", "mean", "lower", "upper", "p_most_central"))
  expect_equal(nrow(summ), 6L)
  expect_false(is.unsorted(rev(summ$mean)))
  expect_equal(sum(summ$p_most_central), 1, tolerance = 1e-12)
  expect_true(all(summ$lower <= summ$mean & summ$mean <= summ$upper))

  # The most central node by mean is the one that most often wins a draw, in
  # this fit; the two orderings can differ, so this is a value check, not a law.
  expect_equal(summ$mean, colMeans(strength)[summ$node], ignore_attr = TRUE)

  wide = summary(strength, probs = c(0.005, 0.995))
  expect_true(all(wide$lower <= summ$lower))
  expect_true(all(wide$upper >= summ$upper))
})

test_that("plot.bgms_centrality draws and returns its argument invisibly", {
  skip_on_cran()
  fit = get_bgms_fit_wenchuan5()
  strength = extract_centrality(fit)

  path = withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(strength))
})

test_that("extract_centrality.bgmCompare rebuilds each group per draw", {
  skip_on_cran()
  fit = get_bgmcompare_fit_wenchuan5()
  arguments = extract_arguments(fit)

  g1 = extract_centrality(fit, group = 1)
  g2 = extract_centrality(fit, group = 2)
  difference = extract_centrality(fit, group = c(1, 2))

  expect_s3_class(g1, "bgms_centrality")
  expect_equal(colnames(g1), arguments$data_columnnames)
  expect_equal(nrow(g1), sum(vapply(get_raw_samples(fit)$pairwise, nrow, integer(1))))

  # The pair is exactly the difference of the two singles, per draw.
  expect_equal(unclass(difference), unclass(g1) - unclass(g2), ignore_attr = TRUE)

  # The per-draw reconstruction averages to what extract_group_params()
  # computes from the posterior means, which is the reference it must match.
  projection = arguments$projection
  draws = do.call(rbind, get_raw_samples(fit)$pairwise)
  num_pairs = ncol(draws) / arguments$num_groups
  baseline = draws[, seq_len(num_pairs), drop = FALSE]
  contrast = draws[, num_pairs + seq_len(num_pairs), drop = FALSE]
  mine = sapply(seq_len(arguments$num_groups), function(g) {
    colMeans(baseline + projection[g, 1] * contrast)
  })
  expect_equal(
    mine, extract_group_params(fit)$pairwise_effects_groups,
    tolerance = 1e-10, ignore_attr = TRUE
  )

  expect_error(extract_centrality(fit, group = 3), "one group index")
  expect_error(extract_centrality(fit, group = c(1, 2, 1)), "one group index")
  expect_error(extract_centrality(fit, group = c(2, 2)), "must be different")
})

test_that("a centrality difference is summarized and drawn against zero", {
  skip_on_cran()
  fit = get_bgmcompare_fit_wenchuan5()
  difference = extract_centrality(fit, group = c(1, 2))
  summ = summary(difference)

  # "Most central" is not the question a difference answers.
  expect_named(summ, c("node", "mean", "lower", "upper", "p_positive"))
  expect_equal(summ$p_positive, colMeans(unclass(difference) > 0)[summ$node],
    ignore_attr = TRUE
  )

  # Draws in which every one of a node's difference indicators is excluded give
  # the two groups the same network, so the difference is exactly zero there.
  # That point mass is why p_positive and its mirror need not sum to one.
  zero_share = colMeans(unclass(difference) == 0)
  expect_true(any(zero_share > 0))
  expect_true(all(summ$p_positive + colMeans(unclass(difference) < 0)[summ$node] <= 1))

  # A single group keeps the level summary.
  expect_named(
    summary(extract_centrality(fit, group = 1)),
    c("node", "mean", "lower", "upper", "p_most_central")
  )

  path = withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(difference))
})

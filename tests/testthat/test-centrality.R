test_that("extract_centrality returns draws x nodes strength", {
  skip_on_cran()
  fit = bgm(Wenchuan[, 1:6],
    chains = 2, iter = 400, warmup = 400, seed = 1,
    display_progress = "none", verbose = FALSE
  )
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
  fit = bgm(Wenchuan[, 1:6],
    chains = 2, iter = 400, warmup = 400, seed = 1,
    display_progress = "none", verbose = FALSE
  )
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
  fit = bgm(Wenchuan[, 1:6],
    chains = 2, iter = 400, warmup = 400, seed = 1,
    display_progress = "none", verbose = FALSE
  )
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
  fit = bgm(Wenchuan[, 1:5],
    chains = 2, iter = 200, warmup = 200, seed = 4,
    display_progress = "none", verbose = FALSE
  )
  strength = extract_centrality(fit)

  path = withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(strength))
})

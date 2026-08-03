# ==============================================================================
# Unit tests for bgm_spec()
# Phase B.5 of the R scaffolding refactor.
# ==============================================================================

# These are internal (non-exported) functions <U+2014> bind via ::: for testing.
bgm_spec = bgms:::bgm_spec
validate_bgm_spec = bgms:::validate_bgm_spec
new_bgm_spec = bgms:::new_bgm_spec
sampler_sublist = bgms:::sampler_sublist

# ==============================================================================
# Shared helpers / fixtures
# ==============================================================================

# Minimal continuous data (GGM)
make_continuous_data = function(n = 20, p = 3) {
  set.seed(42)
  x = matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(x) = paste0("V", seq_len(p))
  x
}

# Minimal ordinal data (OMRF / Compare) <U+2014> values 0,1,2
make_ordinal_data = function(n = 30, p = 3, max_cat = 2) {
  set.seed(99)
  x = matrix(sample(0:max_cat, n * p, replace = TRUE), nrow = n, ncol = p)
  colnames(x) = paste0("V", seq_len(p))
  x
}

# Helper that calls bgm_spec() with sensible defaults, easy to override
spec = function(...) {
  defaults = list(
    x = make_ordinal_data(),
    model_type = "omrf",
    variable_type = "ordinal",
    na_action = "listwise",
    update_method = "nuts",
    iter = 500L,
    warmup = 100L,
    chains = 2L,
    cores = 1L,
    seed = 1L,
    display_progress = "none",
    verbose = FALSE
  )
  args = modifyList(defaults, list(...))
  do.call(bgm_spec, args)
}


# ==============================================================================
# 1.  S3 class & structure
# ==============================================================================

test_that("bgm_spec returns an object of class 'bgm_spec'", {
  s = spec()
  expect_s3_class(s, "bgm_spec")
})

test_that("bgm_spec has top-level elements", {
  s = spec()
  nms = names(s)
  expected = c(
    "model_type", "data", "variables", "missing",
    "prior", "sampler", "precomputed"
  )
  expect_true(all(expected %in% nms))
})


# ==============================================================================
# 2.  GGM path
# ==============================================================================

test_that("GGM spec: model_type resolved correctly", {
  s = spec(
    x = make_continuous_data(), model_type = "ggm",
    variable_type = "continuous",
    update_method = "adaptive-metropolis"
  )
  expect_equal(s$model_type, "ggm")
  expect_true(s$variables$is_continuous)
})

test_that("GGM spec: sampler forced to adaptive-metropolis", {
  s = spec(
    x = make_continuous_data(), model_type = "ggm",
    variable_type = "continuous",
    update_method = "adaptive-metropolis"
  )
  expect_equal(s$sampler$update_method, "adaptive-metropolis")
})

test_that("GGM spec: data sub-list fields", {
  x = make_continuous_data(n = 15, p = 4)
  s = spec(
    x = x, model_type = "ggm", variable_type = "continuous",
    update_method = "adaptive-metropolis"
  )
  expect_true(is.matrix(s$data$x))
  expect_equal(s$data$num_variables, 4L)
  expect_equal(s$data$num_cases, 15L)
  expect_equal(s$data$data_columnnames, colnames(x))
})

test_that("GGM spec: edge_selection = FALSE sets edge_prior to 'Not Applicable'", {
  s = spec(
    x = make_continuous_data(), model_type = "ggm",
    variable_type = "continuous",
    update_method = "adaptive-metropolis",
    edge_selection = FALSE
  )
  expect_false(s$prior$edge_selection)
  expect_equal(s$prior$edge_prior, "Not Applicable")
})

test_that("GGM spec: edge_selection = TRUE with Bernoulli prior", {
  s = spec(
    x = make_continuous_data(), model_type = "ggm",
    variable_type = "continuous",
    update_method = "adaptive-metropolis",
    edge_selection = TRUE,
    edge_prior = "Bernoulli",
    inclusion_probability = 0.3
  )
  expect_true(s$prior$edge_selection)
  expect_equal(s$prior$edge_prior, "Bernoulli")
  expect_true(is.matrix(s$prior$inclusion_probability))
  # Off-diagonal should be 0.3
  nv = s$data$num_variables
  ip = s$prior$inclusion_probability
  expect_true(all(ip[upper.tri(ip)] == 0.3))
})

test_that("GGM spec: missing data fields", {
  s = spec(
    x = make_continuous_data(), model_type = "ggm",
    variable_type = "continuous",
    update_method = "adaptive-metropolis",
    na_action = "listwise"
  )
  expect_equal(s$missing$na_action, "listwise")
  expect_false(s$missing$na_impute)
})

test_that("auto-promotion: omrf + continuous => ggm", {
  s = spec(
    x = make_continuous_data(), model_type = "omrf",
    variable_type = "continuous",
    update_method = "adaptive-metropolis"
  )
  expect_equal(s$model_type, "ggm")
})


# ==============================================================================
# 3.  OMRF path
# ==============================================================================

test_that("OMRF spec: model_type = 'omrf' for ordinal data", {
  s = spec()
  expect_equal(s$model_type, "omrf")
})

test_that("OMRF spec: data sub-list has num_categories", {
  s = spec()
  expect_true(!is.null(s$data$num_categories))
  expect_equal(length(s$data$num_categories), s$data$num_variables)
})

test_that("OMRF spec: precomputed has num_thresholds", {
  s = spec()
  expect_true(!is.null(s$precomputed$num_thresholds))
  expect_true(is.integer(s$precomputed$num_thresholds))
})

test_that("OMRF spec: sampler fields are properly coerced", {
  s = spec(iter = 500L, warmup = 100L, chains = 2L)
  expect_true(is.integer(s$sampler$iter))
  expect_true(is.integer(s$sampler$warmup))
  expect_true(is.integer(s$sampler$chains))
  expect_true(is.integer(s$sampler$cores))
  expect_true(is.integer(s$sampler$seed))
})

test_that("OMRF spec: Blume-Capel variable types work", {
  x = make_ordinal_data(n = 30, p = 3, max_cat = 2)
  s = spec(
    x = x, variable_type = c("ordinal", "blume-capel", "ordinal"),
    baseline_category = c(0L, 1L, 0L)
  )
  expect_equal(s$variables$is_ordinal, c(TRUE, FALSE, TRUE))
})

test_that("OMRF spec: edge_prior = 'Beta-Bernoulli' preserves alpha/beta", {
  s = spec(
    edge_selection = TRUE, edge_prior = "Beta-Bernoulli",
    beta_bernoulli_alpha = 2, beta_bernoulli_beta = 3
  )
  expect_equal(s$prior$beta_bernoulli_alpha, 2)
  expect_equal(s$prior$beta_bernoulli_beta, 3)
})

test_that("OMRF spec: edge_prior = 'Stochastic-Block' preserves all params", {
  s = spec(
    edge_selection = TRUE, edge_prior = "Stochastic-Block",
    beta_bernoulli_alpha_between = 1.5,
    beta_bernoulli_beta_between = 2.5,
    dirichlet_alpha = 0.5, lambda = 2
  )
  expect_equal(s$prior$beta_bernoulli_alpha_between, 1.5)
  expect_equal(s$prior$beta_bernoulli_beta_between, 2.5)
  expect_equal(s$prior$dirichlet_alpha, 0.5)
  expect_equal(s$prior$lambda, 2)
})


# ==============================================================================
# 4.  Compare path (x + y)
# ==============================================================================

test_that("Compare spec: basic two-group x/y", {
  x1 = make_ordinal_data(n = 20, p = 3, max_cat = 2)
  x2 = make_ordinal_data(n = 25, p = 3, max_cat = 2)
  s = spec(
    x = x1, y = x2, model_type = "compare",
    update_method = "nuts"
  )
  expect_equal(s$model_type, "compare")
  expect_equal(s$data$num_groups, 2L)
  expect_true(is.matrix(s$data$group_indices))
  expect_true(is.matrix(s$data$projection))
})

test_that("Compare spec: group_indicator path", {
  x = make_ordinal_data(n = 40, p = 3, max_cat = 2)
  gi = rep(c(1L, 2L), each = 20)
  s = spec(
    x = x, group_indicator = gi, model_type = "compare",
    update_method = "nuts"
  )
  expect_equal(s$model_type, "compare")
  expect_equal(s$data$num_groups, 2L)
})

test_that("Compare spec: three groups via group_indicator", {
  x = make_ordinal_data(n = 60, p = 3, max_cat = 2)
  gi = rep(1:3, each = 20)
  s = spec(
    x = x, group_indicator = gi, model_type = "compare",
    update_method = "nuts"
  )
  expect_equal(s$data$num_groups, 3L)
  expect_equal(nrow(s$data$group_indices), 3L)
  expect_equal(ncol(s$data$projection), 2L) # num_groups - 1
})

test_that("Compare spec: prior sub-list has difference fields", {
  x1 = make_ordinal_data(n = 20, p = 3, max_cat = 2)
  x2 = make_ordinal_data(n = 25, p = 3, max_cat = 2)
  s = spec(
    x = x1, y = x2, model_type = "compare",
    update_method = "nuts", difference_selection = TRUE
  )
  expect_true(s$prior$difference_selection)
  expect_true(s$prior$difference_prior %in% c("Bernoulli", "Beta-Bernoulli"))
  expect_true(is.matrix(s$prior$inclusion_probability_difference))
})

test_that("Compare spec: difference_selection = FALSE", {
  x1 = make_ordinal_data(n = 20, p = 3, max_cat = 2)
  x2 = make_ordinal_data(n = 25, p = 3, max_cat = 2)
  s = spec(
    x = x1, y = x2, model_type = "compare",
    update_method = "nuts", difference_selection = FALSE
  )
  expect_false(s$prior$difference_selection)
})

test_that("Compare spec: precomputed structures present", {
  x1 = make_ordinal_data(n = 20, p = 3, max_cat = 2)
  x2 = make_ordinal_data(n = 25, p = 3, max_cat = 2)
  s = spec(
    x = x1, y = x2, model_type = "compare",
    update_method = "nuts"
  )
  pc = s$precomputed
  expect_true(!is.null(pc$counts_per_category))
  expect_true(!is.null(pc$blume_capel_stats))
  expect_true(!is.null(pc$pairwise_stats))
  expect_true(!is.null(pc$main_effect_indices))
  expect_true(!is.null(pc$pairwise_effect_indices))
})


# ==============================================================================
# 5.  Error paths (user-facing bgm_spec)
# ==============================================================================

test_that("error: non-matrix input", {
  expect_error(spec(x = list(1, 2, 3)), "matrix|data.frame")
})

test_that("error: single row", {
  x = matrix(0:2, nrow = 1)
  expect_error(spec(x = x), "2 rows")
})

test_that("error: compare without y or group_indicator", {
  expect_error(
    spec(model_type = "compare", update_method = "nuts"),
    "group_indicator|y"
  )
})

test_that("error: compare with y having different ncol", {
  x = make_ordinal_data(n = 20, p = 3)
  y = make_ordinal_data(n = 20, p = 4)
  expect_error(
    spec(
      x = x, y = y, model_type = "compare",
      update_method = "nuts"
    ),
    "same number of columns"
  )
})

test_that("error: compare group_indicator length mismatch", {
  x = make_ordinal_data(n = 20, p = 3)
  gi = rep(1:2, each = 5) # length 10, not 20
  expect_error(
    spec(
      x = x, group_indicator = gi, model_type = "compare",
      update_method = "nuts"
    ),
    "Length of group_indicator"
  )
})

test_that("error: compare group_indicator with one group", {
  x = make_ordinal_data(n = 20, p = 3)
  gi = rep(1L, 20)
  expect_error(
    spec(
      x = x, group_indicator = gi, model_type = "compare",
      update_method = "nuts"
    ),
    "only one group"
  )
})


# ==============================================================================
# 6.  validate_bgm_spec cross-field invariants
# ==============================================================================

test_that("validate: GGM requires is_continuous = TRUE", {
  # Build a minimal spec that violates the invariant
  expect_error(
    validate_bgm_spec(
      structure(list(
        model_type = "ggm",
        variables  = list(is_continuous = FALSE)
      ), class = "bgm_spec")
    ),
    "is_continuous"
  )
})

test_that("validate: GGM allows nuts", {
  expect_silent(
    validate_bgm_spec(
      structure(list(
        model_type = "ggm",
        variables  = list(is_continuous = TRUE),
        sampler    = list(update_method = "nuts"),
        prior      = list(edge_selection = TRUE, edge_prior = "Bernoulli")
      ), class = "bgm_spec")
    )
  )
})

test_that("validate: compare requires num_groups >= 2", {
  expect_error(
    validate_bgm_spec(
      structure(list(
        model_type = "compare",
        data = list(
          group = 1L, num_groups = 1L,
          num_variables = 3L, num_categories = c(3L, 3L, 3L)
        ),
        variables = list(is_continuous = FALSE),
        sampler = list(update_method = "nuts"),
        prior = list(
          difference_selection = TRUE,
          difference_prior = "Bernoulli"
        )
      ), class = "bgm_spec")
    ),
    "num_groups"
  )
})

test_that("validate: edge_selection + edge_prior inconsistency", {
  expect_error(
    validate_bgm_spec(
      structure(list(
        model_type = "omrf",
        data = list(num_variables = 3L, num_categories = c(3L, 3L, 3L)),
        variables = list(is_continuous = FALSE),
        sampler = list(update_method = "nuts"),
        prior = list(
          edge_selection = TRUE,
          edge_prior = "Not Applicable"
        )
      ), class = "bgm_spec")
    ),
    "Not Applicable"
  )
})

# ==============================================================================
# 7.  new_bgm_spec type assertions
# ==============================================================================

test_that("new_bgm_spec: rejects invalid model_type", {
  expect_error(
    new_bgm_spec(
      model_type = "unknown",
      data = list(
        x = matrix(0, 2, 2), data_columnnames = c("A", "B"),
        num_variables = 2L, num_cases = 2L
      ),
      variables = list(
        variable_type = "ordinal", is_ordinal = c(TRUE, TRUE),
        is_continuous = FALSE, baseline_category = c(0L, 0L)
      ),
      missing = list(
        na_action = "listwise", na_impute = FALSE,
        missing_index = NULL
      ),
      prior = list(),
      sampler = list(
        update_method = "nuts", target_accept = 0.8,
        iter = 100L, warmup = 50L, chains = 2L, cores = 1L,
        nuts_max_depth = 10L,
        learn_mass_matrix = TRUE, seed = 1L, progress_type = 0L
      )
    ),
    "model_type"
  )
})

test_that("new_bgm_spec: rejects non-matrix x", {
  expect_error(
    new_bgm_spec(
      model_type = "omrf",
      data = list(
        x = data.frame(a = 1, b = 2), data_columnnames = c("a", "b"),
        num_variables = 2L, num_cases = 1L,
        num_categories = c(2L, 2L)
      ),
      variables = list(
        variable_type = "ordinal", is_ordinal = c(TRUE, TRUE),
        is_continuous = FALSE, baseline_category = c(0L, 0L)
      ),
      missing = list(
        na_action = "listwise", na_impute = FALSE,
        missing_index = NULL
      ),
      prior = list(
        pairwise_scale = 1, main_alpha = 0.5, main_beta = 0.5,
        edge_selection = TRUE, edge_prior = "Bernoulli",
        inclusion_probability = matrix(0.5, 2, 2)
      ),
      sampler = list(
        update_method = "nuts", target_accept = 0.8,
        iter = 100L, warmup = 50L, chains = 2L, cores = 1L,
        nuts_max_depth = 10L,
        learn_mass_matrix = TRUE, seed = 1L, progress_type = 0L
      )
    ),
    "is.matrix"
  )
})


# ==============================================================================
# 8.  sampler_sublist
# ==============================================================================

test_that("sampler_sublist coerces to integer", {
  s = list(
    update_method     = "nuts",
    target_accept     = 0.8,
    iter              = 1000, # numeric, not integer
    warmup            = 250,
    chains            = 2,
    cores             = 1,
    nuts_max_depth    = 10,
    learn_mass_matrix = TRUE,
    seed              = 42,
    progress_type     = 0
  )
  res = sampler_sublist(s)
  expect_true(is.integer(res$iter))
  expect_true(is.integer(res$warmup))
  expect_true(is.integer(res$chains))
  expect_true(is.integer(res$seed))
  expect_true(is.integer(res$progress_type))
})


# ==============================================================================
# 9.  print.bgm_spec
# ==============================================================================

test_that("print.bgm_spec outputs key info", {
  s = spec()
  out = capture.output(print(s))
  expect_true(any(grepl("bgm_spec object", out)))
  expect_true(any(grepl("model_type.*omrf", out)))
  expect_true(any(grepl("sampler.*nuts", out)))
})

test_that("print.bgm_spec for compare shows group info", {
  x1 = make_ordinal_data(n = 20, p = 3, max_cat = 2)
  x2 = make_ordinal_data(n = 25, p = 3, max_cat = 2)
  s = spec(
    x = x1, y = x2, model_type = "compare",
    update_method = "nuts"
  )
  out = capture.output(print(s))
  expect_true(any(grepl("groups:", out)))
})

test_that("print.bgm_spec returns invisible", {
  s = spec()
  capture.output(expect_invisible(print(s)))
})


# ==============================================================================
# 10.  Missing data handling
# ==============================================================================

test_that("OMRF spec: listwise deletion with NAs", {
  x = make_ordinal_data(n = 30, p = 3)
  x[1, 1] = NA
  x[5, 2] = NA
  s = spec(x = x, na_action = "listwise")
  expect_equal(s$missing$na_action, "listwise")
  expect_false(s$missing$na_impute)
  expect_equal(s$data$num_cases, 28L)
})

test_that("OMRF spec: impute path", {
  x = make_ordinal_data(n = 30, p = 3)
  x[1, 1] = NA
  s = spec(x = x, na_action = "impute")
  expect_equal(s$missing$na_action, "impute")
  expect_true(s$missing$na_impute)
  expect_equal(s$data$num_cases, 30L)
  expect_true(is.matrix(s$missing$missing_index))
  expect_equal(nrow(s$missing$missing_index), 1L)
})


# ==============================================================================
# 11.  Mixed MRF: baseline_category subsetting
# ==============================================================================
# Regression tests for a bug where build_spec_mixed_mrf passed a full-length
# baseline_category vector (num_variables) to validate_baseline_category,
# which received x_disc (only discrete columns). The length mismatch caused
# a spurious "needs to be a single integer or a vector of integers of length p"
# error.

make_mixed_data = function(n = 30) {
  set.seed(42)
  x = cbind(
    sample(0:2, n, replace = TRUE),
    sample(0:2, n, replace = TRUE),
    rnorm(n),
    rnorm(n)
  )
  colnames(x) = c("d1", "d2", "c1", "c2")
  x
}

test_that("mixed MRF: full-length baseline_category with NAs for continuous", {
  s = spec(
    x = make_mixed_data(),
    variable_type = c("blume-capel", "blume-capel", "continuous", "continuous"),
    baseline_category = c(1, 1, NA, NA),
    model_type = "omrf"
  )
  expect_equal(s$model_type, "mixed_mrf")
  expect_equal(s$variables$baseline_category[1:2], c(1L, 1L))
})

test_that("mixed MRF: full-length baseline_category with dummy values for continuous", {
  s = spec(
    x = make_mixed_data(),
    variable_type = c("blume-capel", "blume-capel", "continuous", "continuous"),
    baseline_category = c(1, 1, 1, 1),
    model_type = "omrf"
  )
  expect_equal(s$model_type, "mixed_mrf")
  expect_equal(s$variables$baseline_category[1:2], c(1L, 1L))
})

test_that("mixed MRF: scalar baseline_category still works", {
  s = spec(
    x = make_mixed_data(),
    variable_type = c("blume-capel", "blume-capel", "continuous", "continuous"),
    baseline_category = 1,
    model_type = "omrf"
  )
  expect_equal(s$model_type, "mixed_mrf")
  expect_equal(s$variables$baseline_category[1:2], c(1L, 1L))
})

test_that("mixed MRF: disc-length baseline_category works", {
  s = spec(
    x = make_mixed_data(),
    variable_type = c("blume-capel", "blume-capel", "continuous", "continuous"),
    baseline_category = c(1, 1),
    model_type = "omrf"
  )
  expect_equal(s$model_type, "mixed_mrf")
  expect_equal(s$variables$baseline_category[1:2], c(1L, 1L))
})

# ==============================================================================
# 12.  Mixed MRF: determinant-tilt default
# ==============================================================================
# The delta = NULL default is 0.5 * log(p) with p the dimension of the
# continuous precision matrix, so only continuous variables count toward it.

test_that("mixed MRF: default delta counts only continuous variables", {
  s = spec(
    x = make_mixed_data(),
    variable_type = c("blume-capel", "blume-capel", "continuous", "continuous"),
    baseline_category = 1,
    model_type = "omrf"
  )
  expect_equal(s$model_type, "mixed_mrf")
  expect_equal(s$prior$delta, 0.5 * log(2))
})

# ==============================================================================
# bgm_spec() and bgmCompare() must share the same difference_scale default so
# that a direct bgm_spec(model_type = "compare") call matches the public API.
# ==============================================================================

test_that("bgm_spec and bgmCompare share the difference_scale default", {
  expect_equal(
    eval(formals(bgm_spec)$difference_scale),
    eval(formals(bgmCompare)$difference_scale)
  )
})


# ==============================================================================
# 13.  Mixed MRF: degenerate-block guard (F-123)
# ==============================================================================
# The mixed model is the two-block model: its parameter layout, its indicator
# layout and its gradients all assume a non-empty discrete block and a
# non-empty continuous block. bgm_spec() routes pure-type data to the OMRF or
# the GGM, so no exported path reaches the builder with an empty block -- the
# guard is defensive, and the test drives the internal directly.

test_that("mixed MRF: the degenerate-block guard rejects a one-type block", {
  build_spec_mixed_mrf = bgms:::build_spec_mixed_mrf
  # R evaluates arguments lazily and the guard fires before the block
  # machinery is touched, so only the arguments it reads have to be supplied.
  degenerate = function(variable_type) {
    build_spec_mixed_mrf(
      x = matrix(0, 5L, length(variable_type)),
      variable_type = variable_type,
      pairwise_scale = 1, scale_rate = 1, scale_eta = NA_real_
    )
  }

  expect_error(
    degenerate(rep("continuous", 3L)),
    "at least one discrete and at least one continuous"
  )
  expect_error(degenerate(rep("continuous", 3L)), "0 discrete and 3 continuous")

  expect_error(
    degenerate(rep("ordinal", 3L)),
    "at least one discrete and at least one continuous"
  )
  expect_error(degenerate(rep("ordinal", 3L)), "3 discrete and 0 continuous")
})

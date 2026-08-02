# --------------------------------------------------------------------------- #
# SBC for the edge-prior normalizing-constant correction (with data)
#
# Simulation-based calibration run through bgm() end to end. Each replicate
# draws the edge-prior hyperparameters from their hyperprior, draws
# (K, Gamma) from the tilted joint prior AT THOSE FIXED hyperparameters,
# simulates data, fits the corrected model, and ranks the true values among
# the posterior draws. Uniform ranks certify the deployed correction: the
# fixed-hyperparameter prior chains run the plain Bernoulli path and never
# touch the correction tables, so the generator and the corrected fit share
# no correction code, and a wrong logC curve or a mis-wired hyperparameter
# update skews the ranks.
#
# Rank statistics per replicate, covering the full parameter block:
#   - the hyperparameter (theta for beta-bernoulli; the number of occupied
#     blocks for the block model),
#   - the edge count m = sum(Gamma) (the graph-level indicator statistic),
#   - every K diagonal entry,
#   - every K off-diagonal entry (mixed cell: every pairwise effect).
# Discrete statistics and spike-and-slab off-diagonals carry point masses
# that tie plain ranks, so those use the randomized rank
# r = #(draws < true) + U * (#(draws == true) + 1), which is exactly
# Uniform(0, L + 1) under calibration for any posterior with atoms.
#
# Cells:
#   - GGM beta-bernoulli, eta in {1, 2, 3}: the correction table depends on
#     eta through the diagonal rate, so the grid covers the eta-keyed table
#     family.
#   - GGM stochastic block model at delta = 2, with the block structure
#     drawn from the MFM prior and the graph drawn at fixed per-pair
#     inclusion probabilities.
#   - Mixed MRF beta-bernoulli: a prior-only identity at the normal slab.
#     The with-data cells are restricted to the all-continuous models
#     because those fit the exact Gaussian likelihood; the mixed sampler
#     fits a pseudolikelihood, whose pseudo-posterior is not calibrated in
#     the exact-SBC sense (rank uniformity fails by construction, and a
#     sandwich/Godambe recalibration is out of scope for a package test).
#     With no data the pseudolikelihood factor drops out, so the prior-only
#     chain is exact and gates the mixed correction sharply.
#
# The prior-only identity tests (test-mixed-correction.R, test-sbm-
# correction.R) carry further uncorrected controls that show what a missing
# correction does to these statistics.
#
# Weekly certification tier (T2, BGMS_RUN_CERTIFICATION): SBC replicate
# suites. BGMS_SBC_R overrides the replicate
# counts (all cells) for quick wiring runs; BGMS_SBC_CORES caps the fork
# pool (default 5).
# --------------------------------------------------------------------------- #

sbc_replicates = function(default) {
  v = Sys.getenv("BGMS_SBC_R", "")
  if(nzchar(v)) as.integer(v) else default
}

sbc_cores = function() {
  v = as.integer(Sys.getenv("BGMS_SBC_CORES", "5"))
  if(.Platform$OS.type == "unix") max(1L, v) else 1L
}

sbc_map = function(R, f) {
  if(sbc_cores() > 1L) {
    parallel::mclapply(
      seq_len(R), f,
      mc.cores = sbc_cores(), mc.preschedule = FALSE
    )
  } else {
    lapply(seq_len(R), f)
  }
}

sbc_correction_cache = function() {
  options(
    bgms.correction_cache_dir = file.path(tempdir(), "bgms-ctable-sbc")
  )
}

# Randomized rank: exactly Uniform(0, L + 1) under calibration for any
# posterior, including point masses (indicator counts, spike-and-slab
# off-diagonals at zero).
sbc_rank_smooth = function(draws, truth) {
  sum(draws < truth) + runif(1) * (sum(draws == truth) + 1)
}

# Gates on a rank matrix (replicates x statistics): per-statistic KS at
# alpha = 0.01 with the test-sbc-ggm.R false-positive allowance, plus a
# pooled 20-bin chi-squared across all statistics.
sbc_expect_calibrated = function(cell, ranks, L, R) {
  expect_gte(nrow(ranks), ceiling(0.9 * R))

  u = ranks / (L + 1)
  ks_p = apply(u, 2, function(x) {
    suppressWarnings(ks.test(x, "punif")$p.value)
  })
  n_fail = sum(ks_p <= 0.01)
  max_fail = max(1, ceiling(ncol(u) * 0.01 * 2))

  bins = cut(
    as.vector(u),
    breaks = seq(0, 1, length.out = 21), include.lowest = TRUE
  )
  chisq_p = suppressWarnings(chisq.test(tabulate(bins, nbins = 20))$p.value)

  worst = which.min(ks_p)
  cat(sprintf(
    paste0(
      "[sbc-correction] %s: %d/%d reps | %s mean rank %.1f (exp %.1f) ",
      "KS p=%.3f | KS fails %d/%d (worst %s p=%.4f) | pooled chisq p=%.3f\n"
    ),
    cell, nrow(ranks), R,
    colnames(ranks)[1L], mean(ranks[, 1L]), L / 2, ks_p[1L],
    n_fail, ncol(u), colnames(ranks)[worst], ks_p[worst], chisq_p
  ))

  expect_lte(n_fail, max_fail)
  expect_gt(chisq_p, 0.001)
}

# Reconstruct a symmetric p x p K matrix from the (off-diag, diag) rows
# returned by sample_ggm_prior(spec = "joint"); off-diag columns follow the
# (i < j) upper-triangle visit order.
sbc_reconstruct_K = function(off_row, diag_row, p) {
  K = matrix(0, p, p)
  idx = 1L
  for(i in seq_len(p - 1)) {
    for(j in (i + 1):p) {
      K[i, j] = off_row[idx]
      K[j, i] = off_row[idx]
      idx = idx + 1L
    }
  }
  diag(K) = diag_row
  K
}

# Ranks shared by the GGM cells: K diagonal (plain ranks), K off-diagonal
# and edge count (randomized ranks over the point masses).
sbc_ggm_param_ranks = function(fit, thin_idx, k_diag_true, k_off_true,
                               m_true, p) {
  main = fit$raw_samples$main[[1L]][thin_idx, , drop = FALSE]
  pw = fit$raw_samples$pairwise[[1L]][thin_idx, , drop = FALSE]
  ind = fit$raw_samples$indicator[[1L]][thin_idx, , drop = FALSE]

  n_edges = length(k_off_true)
  kdiag = vapply(
    seq_len(p),
    function(i) sum(main[, i] < k_diag_true[i]),
    numeric(1)
  )
  koff = vapply(
    seq_len(n_edges),
    function(e) sbc_rank_smooth(pw[, e], k_off_true[e]),
    numeric(1)
  )
  c(
    m = sbc_rank_smooth(rowSums(ind), m_true),
    setNames(kdiag, paste0("kdiag", seq_len(p))),
    setNames(koff, paste0("koff", seq_len(n_edges)))
  )
}


# --------------------------------------------------------------------------- #
# GGM beta-bernoulli across the eta-keyed table family
# --------------------------------------------------------------------------- #

run_ggm_bb_theta_sbc = function(eta, R, seed_base) {
  p = 5
  n = 1000
  delta = 0.5 * log(p)
  iter = 2000L
  warmup = 1000L
  thin_idx = seq(10L, iter, by = 10L)
  L = length(thin_idx)

  # Build the correction table once in the parent so forked fits read the
  # cached file instead of racing to build it.
  invisible(ggm_correction_table(
    p = p, delta = delta,
    interaction_prior = cauchy_prior(scale = 2.5),
    precision_scale_prior = gamma_prior(shape = 1, eta = eta),
    update_method = "gibbs"
  ))

  one_rep = function(r) {
    set.seed(seed_base + r)
    theta_true = rbeta(1, 1, 1)

    draw = tryCatch(
      sample_ggm_prior(
        p = p, n_samples = 1L, n_warmup = 1000L,
        interaction_prior = cauchy_prior(scale = 2.5),
        precision_scale_prior = gamma_prior(shape = 1, eta = eta),
        delta = delta, spec = "joint",
        edge_inclusion_prob = theta_true,
        update_method = "gibbs",
        verbose = FALSE, seed = seed_base + 7L * r
      ),
      error = function(e) NULL
    )
    if(is.null(draw)) {
      return(NULL)
    }

    K_true = sbc_reconstruct_K(draw$K_offdiag[1L, ], draw$K_diag[1L, ], p)
    Sigma = tryCatch(solve(K_true), error = function(e) NULL)
    if(is.null(Sigma)) {
      return(NULL)
    }

    X = MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
    dat = as.data.frame(X)
    colnames(dat) = paste0("V", seq_len(p))

    fit = tryCatch(
      bgm(dat,
        variable_type = "continuous",
        iter = iter, warmup = warmup, chains = 1, cores = 1,
        edge_selection = TRUE, update_method = "gibbs",
        edge_prior = beta_bernoulli_prior(1, 1),
        interaction_prior = cauchy_prior(scale = 2.5),
        precision_scale_prior = gamma_prior(shape = 1, eta = eta),
        delta = delta,
        display_progress = "none", seed = seed_base + 13L * r
      ),
      error = function(e) NULL
    )
    if(is.null(fit)) {
      return(NULL)
    }

    th = fit$inclusion_parameter_samples[[1L]][thin_idx]
    c(
      theta = sum(th < theta_true),
      sbc_ggm_param_ranks(
        fit, thin_idx,
        k_diag_true = draw$K_diag[1L, ],
        k_off_true = draw$K_offdiag[1L, ],
        m_true = sum(draw$edge_indicators[1L, ]),
        p = p
      )
    )
  }

  ranks = do.call(rbind, Filter(Negate(is.null), sbc_map(R, one_rep)))
  list(ranks = ranks, L = L)
}

test_that("SBC: corrected GGM beta-bernoulli ranks are uniform (eta = 1)", {
  skip_unless_certification()
  old = sbc_correction_cache()
  on.exit(options(old), add = TRUE)

  R = sbc_replicates(200L)
  res = run_ggm_bb_theta_sbc(eta = 1, R = R, seed_base = 51000L)
  sbc_expect_calibrated("ggm-bb eta=1", res$ranks, res$L, R)
})

test_that("SBC: corrected GGM beta-bernoulli ranks are uniform (eta = 2)", {
  skip_unless_certification()
  old = sbc_correction_cache()
  on.exit(options(old), add = TRUE)

  R = sbc_replicates(200L)
  res = run_ggm_bb_theta_sbc(eta = 2, R = R, seed_base = 52000L)
  sbc_expect_calibrated("ggm-bb eta=2", res$ranks, res$L, R)
})

test_that("SBC: corrected GGM beta-bernoulli ranks are uniform (eta = 3)", {
  skip_unless_certification()
  old = sbc_correction_cache()
  on.exit(options(old), add = TRUE)

  R = sbc_replicates(200L)
  res = run_ggm_bb_theta_sbc(eta = 3, R = R, seed_base = 53000L)
  sbc_expect_calibrated("ggm-bb eta=3", res$ranks, res$L, R)
})


# --------------------------------------------------------------------------- #
# Mixed MRF beta-bernoulli: prior-only identity at the normal slab
# --------------------------------------------------------------------------- #

# Prior-only mixed chain under the beta-bernoulli prior: with no data the
# pseudolikelihood factor is absent, so the corrected chain must return the
# Beta(1, 1) hyperprior on theta exactly.
run_mixed_bb_prior_chain = function(apply_correction, n_samples, seed) {
  p_disc = 2L
  q_cont = 3L
  q_tot = p_disc + q_cont
  delta = 0.5 * log(q_cont)
  slab_scale = 1
  ep = unpack_indicator_prior(
    beta_bernoulli_prior(1, 1),
    num_variables = q_tot
  )
  correction = NULL
  if(apply_correction) {
    table = ggm_correction_table(
      p = q_cont, delta = delta,
      interaction_prior = normal_prior(scale = slab_scale),
      precision_scale_prior = gamma_prior(shape = 1, eta = 1),
      update_method = "gibbs"
    )
    correction = correction_list_from_table(
      table, ep$edge_prior,
      is_continuous = c(rep(0L, p_disc), rep(1L, q_cont))
    )
  }
  input = list(
    discrete_observations = matrix(0L, 0, p_disc),
    continuous_observations = matrix(0, 0, q_cont),
    num_categories = rep(1L, p_disc),
    is_ordinal_variable = rep(1L, p_disc),
    baseline_category = rep(0L, p_disc),
    pairwise_scale = slab_scale,
    interaction_prior_type = "normal",
    scale_prior_type = "gamma",
    scale_shape = 1.0,
    scale_rate = 1 / slab_scale
  )
  res = sample_mixed_mrf(
    inputFromR = input,
    prior_inclusion_prob = ep$inclusion_probability,
    initial_edge_indicators = matrix(1L, q_tot, q_tot),
    no_iter = as.integer(n_samples),
    no_warmup = 1000L,
    no_chains = 1L,
    edge_selection = TRUE,
    seed = as.integer(seed),
    no_threads = 1L,
    progress_type = 0L,
    edge_prior = ep$edge_prior,
    beta_bernoulli_alpha = ep$beta_bernoulli_alpha,
    beta_bernoulli_beta = ep$beta_bernoulli_beta,
    beta_bernoulli_alpha_between = bb_between_or_sentinel(ep$beta_bernoulli_alpha_between),
    beta_bernoulli_beta_between = bb_between_or_sentinel(ep$beta_bernoulli_beta_between),
    dirichlet_alpha = ep$dirichlet_alpha,
    lambda = ep$lambda,
    sampler_type = "adaptive-metropolis",
    delta = delta,
    edge_prior_correction = correction
  )
  stopifnot(length(res) == 1L, !isTRUE(res[[1L]]$error))
  as.numeric(res[[1L]]$inclusion_parameter_samples)
}

test_that("SBC: corrected mixed prior chain returns Beta(1, 1) at the normal slab", {
  skip_unless_certification()
  old = sbc_correction_cache()
  on.exit(options(old), add = TRUE)

  th = run_mixed_bb_prior_chain(
    apply_correction = TRUE, n_samples = 20000, seed = 42
  )
  ks_d = unname(suppressWarnings(ks.test(th, "punif")$statistic))
  cat(sprintf(
    "[sbc-correction] mixed-bb identity (normal slab): corrected mean %.4f (exp 0.5), KS D=%.4f\n",
    mean(th), ks_d
  ))
  expect_lt(abs(mean(th) - 0.5), 0.02)
  expect_lt(ks_d, 0.02)
})

test_that("SBC: uncorrected mixed prior chain misses the hyperprior at the normal slab", {
  skip_unless_certification()

  th = run_mixed_bb_prior_chain(
    apply_correction = FALSE, n_samples = 20000, seed = 42
  )
  ks_d = unname(suppressWarnings(ks.test(th, "punif")$statistic))
  cat(sprintf(
    "[sbc-correction] mixed-bb identity (normal slab): uncorrected mean %.4f, KS D=%.4f\n",
    mean(th), ks_d
  ))
  expect_lt(mean(th), 0.45)
  expect_gt(ks_d, 0.1)
})


# --------------------------------------------------------------------------- #
# GGM stochastic block model
# --------------------------------------------------------------------------- #

# Ancestral draw from the MFM-SBM hyperprior: component count, weights,
# allocations, and per-block-pair inclusion probabilities.
sbc_draw_sbm_partition = function(q, lambda, dirichlet_alpha) {
  K = rpois(1, lambda) + 1L
  w = rgamma(K, dirichlet_alpha)
  sample.int(K, q, replace = TRUE, prob = w)
}

sbc_sbm_theta_matrix = function(z, a_within, b_within, a_between, b_between) {
  labs = sort(unique(z))
  nl = length(labs)
  th_rs = matrix(0, nl, nl)
  for(r in seq_len(nl)) {
    for(s in r:nl) {
      v = if(r == s) {
        rbeta(1, a_within, b_within)
      } else {
        rbeta(1, a_between, b_between)
      }
      th_rs[r, s] = v
      th_rs[s, r] = v
    }
  }
  zi = match(z, labs)
  q = length(z)
  tm = matrix(0.5, q, q)
  for(i in seq_len(q - 1)) {
    for(j in (i + 1):q) {
      tm[i, j] = th_rs[zi[i], zi[j]]
      tm[j, i] = tm[i, j]
    }
  }
  tm
}

run_sbm_nc_sbc = function(R, seed_base) {
  p = 5
  n = 1000
  delta = 2
  iter = 3000L
  warmup = 1500L
  thin_idx = seq(15L, iter, by = 15L)
  L = length(thin_idx)

  invisible(ggm_correction_table(
    p = p, delta = delta,
    interaction_prior = cauchy_prior(scale = 2.5),
    precision_scale_prior = gamma_prior(shape = 1, eta = 1),
    update_method = "gibbs"
  ))

  one_rep = function(r) {
    set.seed(seed_base + r)
    z_true = sbc_draw_sbm_partition(p, lambda = 1, dirichlet_alpha = 1)
    nc_true = length(unique(z_true))
    theta_matrix = sbc_sbm_theta_matrix(z_true, 1, 1, 1, 1)

    draw = tryCatch(
      sample_ggm_prior(
        p = p, n_samples = 1L, n_warmup = 1000L,
        interaction_prior = cauchy_prior(scale = 2.5),
        precision_scale_prior = gamma_prior(shape = 1, eta = 1),
        delta = delta, spec = "joint",
        edge_prior = bernoulli_prior(theta_matrix),
        update_method = "gibbs",
        verbose = FALSE, seed = seed_base + 7L * r
      ),
      error = function(e) NULL
    )
    if(is.null(draw)) {
      return(NULL)
    }

    K_true = sbc_reconstruct_K(draw$K_offdiag[1L, ], draw$K_diag[1L, ], p)
    Sigma = tryCatch(solve(K_true), error = function(e) NULL)
    if(is.null(Sigma)) {
      return(NULL)
    }

    X = MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
    dat = as.data.frame(X)
    colnames(dat) = paste0("V", seq_len(p))

    fit = tryCatch(
      bgm(dat,
        variable_type = "continuous",
        iter = iter, warmup = warmup, chains = 1, cores = 1,
        edge_selection = TRUE, update_method = "gibbs",
        edge_prior = sbm_prior(),
        interaction_prior = cauchy_prior(scale = 2.5),
        precision_scale_prior = gamma_prior(shape = 1, eta = 1),
        delta = delta,
        display_progress = "none", seed = seed_base + 13L * r
      ),
      error = function(e) NULL
    )
    if(is.null(fit)) {
      return(NULL)
    }

    alloc = fit$raw_samples$allocations[[1L]][thin_idx, , drop = FALSE]
    nc_draws = apply(alloc, 1, function(z) length(unique(z)))
    c(
      nc = sbc_rank_smooth(nc_draws, nc_true),
      sbc_ggm_param_ranks(
        fit, thin_idx,
        k_diag_true = draw$K_diag[1L, ],
        k_off_true = draw$K_offdiag[1L, ],
        m_true = sum(draw$edge_indicators[1L, ]),
        p = p
      )
    )
  }

  ranks = do.call(rbind, Filter(Negate(is.null), sbc_map(R, one_rep)))
  list(ranks = ranks, L = L)
}

test_that("SBC: corrected GGM SBM ranks are uniform", {
  skip_unless_certification()
  old = sbc_correction_cache()
  on.exit(options(old), add = TRUE)

  R = sbc_replicates(150L)
  res = run_sbm_nc_sbc(R = R, seed_base = 55000L)
  sbc_expect_calibrated("ggm-sbm", res$ranks, res$L, R)
})

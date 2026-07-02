# ==============================================================================
# Correction tables for hierarchical graph priors on the GGM (joint spec).
#
# Under the determinant-tilted joint spike-and-slab prior the per-graph
# normalizer Z(Gamma) tilts the graph law, and the graph-prior mean
# C(theta) = E[Z(Gamma)] enters the full conditional of every graph-prior
# hyperparameter. The conjugate updates (beta-binomial theta draw, stochastic
# block updates) omit the 1/C(theta) factor; the tables built here supply the
# correction as ratios read off one curve per model cell
# (q, delta, eta, slab family).
#
# The estimator is score integration on the edge count: at fixed theta the
# tilted prior chain gives the mean edge density d(theta) = E[m | theta] / E,
# and
#
#   d log C / d theta = E * (d(theta) - theta) / (theta * (1 - theta)),
#
# so the per-pair curve f(theta) = log C(theta) / E is the cumulative
# trapezoid of the per-pair score, up to an additive constant that cancels in
# every ratio. The implied per-edge slope at local density d is
# f'(d) = log[(d / (1 - d)) * (1 - theta) / theta], and the per-edge read-off
# table is fed(theta, d) = log(1 - theta + theta * exp(f'(d))).
#
# The sweep runs bgms's own tilted prior sampler (sample_ggm_prior,
# spec = "joint"), so the tables are consistent with the deployed model by
# construction. At high theta the chain can hit the positive-definite cone
# boundary and return spuriously sparse graphs; the true d(theta) is monotone
# increasing, so isotonic regression repairs those dips before the logs.
# ==============================================================================


# ------------------------------------------------------------------
# ggm_correction_theta_grid (internal)
# ------------------------------------------------------------------
# Theta grid for the score sweep, dense at both ends: the score is
# singular at 0 and 1, and the corrected hyperparameter draw explores
# the tails, so about a quarter of the points go below 0.1 and a
# quarter above 0.9 (log-spaced toward the ends).
# ------------------------------------------------------------------
ggm_correction_theta_grid = function(n_grid = 120L, lower = 0.002,
                                     upper = 0.995) {
  n_tail = as.integer(round(n_grid / 4))
  n_mid = n_grid - 2L * n_tail
  lower_tail = exp(seq(log(lower), log(0.1), length.out = n_tail + 1L))
  upper_tail = 1 - exp(seq(log(0.1), log(1 - upper), length.out = n_tail + 1L))
  sort(unique(c(
    lower_tail[-(n_tail + 1L)],
    seq(0.1, 0.9, length.out = n_mid),
    upper_tail[-1L]
  )))
}


# ------------------------------------------------------------------
# correction_table_from_edens (internal)
# ------------------------------------------------------------------
# Pure table math: from a raw edge-density sweep to the correction
# curves. Kept free of MCMC so it is testable against closed forms.
#
#  - isotonic repair of PD-boundary crash dips (edens is monotone
#    increasing in theta), then clamp away from {0, 1} for the logs
#  - per-pair score -> trapezoid -> f(theta) (per-pair log C)
#  - implied per-edge slope fprime vs local density
#  - fed(theta, d) read-off grid
#
# The integrated f/logC curve uses the full theta grid (the score has
# no boundary singularity in the estimate). The slope curve is
# fprime = logit(d) - logit(theta), so near the boundaries the logit
# amplifies Monte-Carlo noise in d and the resolution floor 0.5/E can
# exceed the true density; fprime and the fed table are therefore
# built only from grid points with theta <= fprime_theta_cap whose
# density is resolvable (inside (0.5/E, 1 - 0.5/E)), and extended as
# constants beyond the covered density range.
#
# Returns the table as a list; logC is the whole-graph curve
# num_pairs * f used by the beta-binomial draw. fprime is tabulated
# against local density (fprime_density, fprime).
# ------------------------------------------------------------------
correction_table_from_edens = function(theta, edens_raw, num_pairs,
                                       fprime_theta_cap = 0.90) {
  stopifnot(length(theta) == length(edens_raw), !is.unsorted(theta))

  edens = stats::isoreg(theta, edens_raw)$yf
  num_repaired = sum(edens_raw < edens - 1e-6)

  score = (edens - theta) / (theta * (1 - theta))
  n = length(score)
  f = c(0, cumsum(diff(theta) * (score[-n] + score[-1]) / 2))
  stopifnot(all(is.finite(f)))

  keep = theta <= fprime_theta_cap &
    edens > 0.5 / num_pairs & edens < 1 - 0.5 / num_pairs
  stopifnot(any(keep))
  fprime_density = edens[keep]
  fprime = log(
    (fprime_density / (1 - fprime_density)) * (1 - theta[keep]) / theta[keep]
  )
  stopifnot(all(is.finite(fprime)))

  fed_theta = seq(0.005, 0.995, length.out = 120L)
  fed_density = seq(min(fprime_density), max(fprime_density),
    length.out = 60L
  )
  fprime_at = function(d) {
    stats::approx(fprime_density, fprime, d, rule = 2, ties = mean)$y
  }
  fed = outer(fed_theta, fed_density, function(tt, dd) {
    log(1 - tt + tt * exp(fprime_at(dd)))
  })

  list(
    theta = theta,
    edens = edens,
    edens_raw = edens_raw,
    num_repaired = num_repaired,
    f = f,
    logC = num_pairs * f,
    fprime_density = fprime_density,
    fprime = fprime,
    fprime_theta_cap = fprime_theta_cap,
    fed_theta = fed_theta,
    fed_density = fed_density,
    fed = fed,
    num_pairs = num_pairs
  )
}


# ------------------------------------------------------------------
# normalize_builder_cores (internal)
# ------------------------------------------------------------------
# The sweep forks full R sessions via mclapply: honor R CMD check's
# core limit, never fork on Windows, and never exceed the machine.
# ------------------------------------------------------------------
normalize_builder_cores = function(cores) {
  cores = max(1L, as.integer(cores))
  if(identical(.Platform$OS.type, "windows")) {
    return(1L)
  }
  check_limit = Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if(nzchar(check_limit) && !identical(tolower(check_limit), "false")) {
    cores = min(cores, 2L)
  }
  min(cores, parallel::detectCores())
}


# ------------------------------------------------------------------
# sweep_prior_edge_density (internal)
# ------------------------------------------------------------------
# Run the tilted prior chain at each theta on the grid and record the
# mean edge density, averaged over n_seeds independent chains per
# theta. cores > 1 forks over (theta, seed) cells.
# ------------------------------------------------------------------
sweep_prior_edge_density = function(p, theta, delta,
                                    interaction_prior,
                                    precision_scale_prior,
                                    n_samples = 2000L, n_warmup = 500L,
                                    n_seeds = 3L,
                                    update_method = "gibbs",
                                    cores = 1L, base_seed = 1L) {
  cores = normalize_builder_cores(cores)
  num_pairs = p * (p - 1) / 2
  cells = expand.grid(theta = theta, seed = seq_len(n_seeds))

  one_cell = function(k) {
    draws = sample_ggm_prior(
      p = p, n_samples = n_samples, n_warmup = n_warmup,
      interaction_prior = interaction_prior,
      precision_scale_prior = precision_scale_prior,
      seed = as.integer(base_seed + k),
      verbose = FALSE, delta = delta,
      spec = "joint",
      edge_inclusion_prob = cells$theta[k],
      update_method = update_method
    )
    mean(draws$edge_indicators)
  }

  edens_cells = if(cores > 1L) {
    unlist(parallel::mclapply(
      seq_len(nrow(cells)), one_cell,
      mc.cores = cores, mc.preschedule = FALSE
    ))
  } else {
    vapply(seq_len(nrow(cells)), one_cell, numeric(1))
  }
  if(anyNA(edens_cells)) {
    stop("Correction-table sweep: a prior chain failed.")
  }

  agg = stats::aggregate(
    edens ~ theta,
    data.frame(theta = cells$theta, edens = edens_cells),
    mean
  )
  agg = agg[order(agg$theta), ]
  list(theta = agg$theta, edens_raw = agg$edens, num_pairs = num_pairs)
}


# ------------------------------------------------------------------
# ggm_correction_cell (internal)
# ------------------------------------------------------------------
# The model cell that keys a correction table. The table depends on
# the tilted within-model prior only through (q, delta, eta, slab
# family, diagonal shape): the slab scale cancels in every
# between-graph ratio at fixed eta = pairwise_scale * scale_rate.
# ------------------------------------------------------------------
ggm_correction_cell = function(p, delta, interaction_prior,
                               precision_scale_prior) {
  ip = unpack_interaction_prior(interaction_prior)
  sp = unpack_scale_prior(precision_scale_prior)
  rate = resolve_scale_rate(sp$scale_rate, sp$scale_eta, ip$pairwise_scale)
  list(
    q = as.integer(p),
    delta = as.numeric(delta),
    eta = ip$pairwise_scale * rate,
    slab_family = ip$interaction_prior_type,
    scale_shape = sp$scale_shape
  )
}


# ------------------------------------------------------------------
# build_ggm_correction_table (internal)
# ------------------------------------------------------------------
# Sweep + table math + cell metadata. delta = NULL resolves to the
# bgm() default 0.5 * log(p).
# ------------------------------------------------------------------
build_ggm_correction_table = function(
  p, delta = NULL,
  interaction_prior = cauchy_prior(scale = 2.5),
  precision_scale_prior = gamma_prior(shape = 1, eta = 1),
  n_grid = 120L, n_samples = 2000L, n_warmup = 500L, n_seeds = 3L,
  update_method = c("gibbs", "adaptive-metropolis"),
  cores = 1L, base_seed = 1L
) {
  update_method = match.arg(update_method)
  if(is.null(delta)) {
    delta = 0.5 * log(p)
  }

  theta = ggm_correction_theta_grid(n_grid)
  sweep = sweep_prior_edge_density(
    p = p, theta = theta, delta = delta,
    interaction_prior = interaction_prior,
    precision_scale_prior = precision_scale_prior,
    n_samples = n_samples, n_warmup = n_warmup, n_seeds = n_seeds,
    update_method = update_method, cores = cores, base_seed = base_seed
  )

  table = correction_table_from_edens(
    sweep$theta, sweep$edens_raw, sweep$num_pairs
  )
  table$cell = ggm_correction_cell(
    p, delta, interaction_prior, precision_scale_prior
  )
  table$builder = list(
    n_grid = n_grid, n_samples = n_samples, n_warmup = n_warmup,
    n_seeds = n_seeds, update_method = update_method,
    base_seed = base_seed, version = 1L
  )
  table
}


# ------------------------------------------------------------------
# ggm_edge_prior_correction (internal)
# ------------------------------------------------------------------
# Resolve whether a GGM fit needs the normalizing-constant correction
# and get-or-build the table for its model cell. Applies to fits with
# edge selection and a Beta-Bernoulli edge prior (the stochastic block
# corrections are separate). The tilted prior sweep runs the same
# priors as the fit; the prior sampler does not support a beta-prime
# slab, so those fits keep the uncorrected update with a warning.
#
# Returns list(theta =, logC =), both NULL when no correction applies.
# ------------------------------------------------------------------
ggm_edge_prior_correction = function(prior, sampler, num_variables) {
  none = list(theta = NULL, logC = NULL)
  if(!isTRUE(prior$edge_selection)) {
    return(none)
  }
  if(!identical(prior$edge_prior, "Beta-Bernoulli")) {
    return(none)
  }
  if(!prior$interaction_prior_type %in% c("cauchy", "normal")) {
    warning(
      "The Beta-Bernoulli inclusion-probability update is run without the ",
      "normalizing-constant correction: the tilted prior sampler supports ",
      "only cauchy_prior() and normal_prior() interaction priors.",
      call. = FALSE
    )
    return(none)
  }

  interaction_prior = switch(prior$interaction_prior_type,
    cauchy = cauchy_prior(scale = prior$pairwise_scale),
    normal = normal_prior(scale = prior$pairwise_scale)
  )
  precision_scale_prior = if(identical(prior$scale_prior_type, "exponential")) {
    exponential_prior(rate = prior$scale_rate)
  } else {
    gamma_prior(shape = prior$scale_shape, rate = prior$scale_rate)
  }

  if(isTRUE(sampler$verbose)) {
    message(
      "Edge-prior correction: building or loading the normalizing-constant ",
      "table for this model (cached across fits)."
    )
  }
  table = ggm_correction_table(
    p = num_variables, delta = prior$delta,
    interaction_prior = interaction_prior,
    precision_scale_prior = precision_scale_prior,
    update_method = "gibbs",
    cores = sampler$cores
  )
  list(theta = table$theta, logC = table$logC)
}


# ------------------------------------------------------------------
# ggm_correction_table (internal)
# ------------------------------------------------------------------
# Cache wrapper: get-or-build the table for a model cell. Tables are
# cached on disk keyed by cell + builder settings, so repeat fits of
# the same cell skip the sweep. Controlled by
# options(bgms.correction_table_cache = FALSE) and
# options(bgms.correction_cache_dir = <path>).
# ------------------------------------------------------------------
ggm_correction_table = function(
  p, delta = NULL,
  interaction_prior = cauchy_prior(scale = 2.5),
  precision_scale_prior = gamma_prior(shape = 1, eta = 1),
  n_grid = 120L, n_samples = 2000L, n_warmup = 500L, n_seeds = 3L,
  update_method = c("gibbs", "adaptive-metropolis"),
  cores = 1L, base_seed = 1L, refresh = FALSE
) {
  update_method = match.arg(update_method)
  if(is.null(delta)) {
    delta = 0.5 * log(p)
  }

  use_cache = isTRUE(getOption("bgms.correction_table_cache", TRUE))
  cache_file = NULL
  if(use_cache) {
    cell = ggm_correction_cell(
      p, delta, interaction_prior, precision_scale_prior
    )
    key = sprintf(
      "ggm_ctable_v1_q%d_delta%.8g_eta%.8g_%s_shape%.8g_g%d_ns%d_nw%d_sd%d_%s.rds",
      cell$q, cell$delta, cell$eta, cell$slab_family, cell$scale_shape,
      n_grid, n_samples, n_warmup, n_seeds, update_method
    )
    cache_dir = getOption(
      "bgms.correction_cache_dir",
      tools::R_user_dir("bgms", which = "cache")
    )
    cache_file = file.path(cache_dir, key)
    if(!refresh && file.exists(cache_file)) {
      table = tryCatch(readRDS(cache_file), error = function(e) NULL)
      if(!is.null(table)) {
        return(table)
      }
    }
  }

  table = build_ggm_correction_table(
    p = p, delta = delta,
    interaction_prior = interaction_prior,
    precision_scale_prior = precision_scale_prior,
    n_grid = n_grid, n_samples = n_samples, n_warmup = n_warmup,
    n_seeds = n_seeds, update_method = update_method,
    cores = cores, base_seed = base_seed
  )

  if(use_cache) {
    dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(table, cache_file)
  }
  table
}

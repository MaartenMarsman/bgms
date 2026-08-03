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
#
# This table serves the joint spec only and is theta-dependent by nature:
# C(theta) = E[Z(Gamma) | theta] is a property of the graph law at a given
# inclusion probability. The hierarchical spec's Option-B absolute-moment
# surface corrects a different quantity -- the per-edge normalizer *ratio* of a
# single toggle, which is theta-free -- and so does not and cannot replace this
# table. Both corrections coexist: the surface on the hierarchical route, this
# table on the joint route.
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
# constants beyond the covered density range. With a single pair the
# resolvable window is empty, so the slope pieces are returned as
# NULL; the logC curve does not depend on them.
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
  fprime_density = NULL
  fprime = NULL
  fed_theta = NULL
  fed_density = NULL
  fed = NULL
  if(any(keep)) {
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
  }

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
  if(identical(.Platform$OS.type, "windows")) {
    return(1L)
  }
  normalize_parallel_cores(cores)
}


# ------------------------------------------------------------------
# new_correction_progress (internal)
# ------------------------------------------------------------------
# In-place progress bar matching the bgms MCMC bar (see the theme and
# bar-width logic in src/utils/progress_manager.cpp): tortoise-shell
# brackets, a heavy horizontal rule filled in blue and empty in gray,
# a sub-cell partial glyph, and a "prefix: <bar> cur/tot (xx.x%)"
# layout redrawn with a carriage return. Unicode + ANSI colour when the
# session is UTF-8, an ASCII "[=== ]" fallback otherwise. Returns
# update(current) and close(); the caller decides whether to draw it.
# ------------------------------------------------------------------
new_correction_progress = function(total, prefix = "Correction table") {
  unicode = isTRUE(l10n_info()[["UTF-8"]])
  is_rstudio = Sys.getenv("RSTUDIO") == "1"

  # Bar width, mirroring progress_manager.cpp so the bar does not wrap.
  console_width = if(is_rstudio) {
    max(0L, as.integer(getOption("width", 80L))) + 3L
  } else {
    80L
  }
  line_width = if(is_rstudio) {
    max(10L, min(console_width - 25L, 70L))
  } else {
    70L
  }
  bar_width = if(line_width <= 5L) {
    0L
  } else if(line_width < 20L) {
    line_width - 10L
  } else if(line_width < 40L) {
    line_width - 15L
  } else if(line_width > 70L) {
    40L
  } else {
    line_width - 30L
  }
  if(is_rstudio) {
    bar_width = if(bar_width > 30L) bar_width - 20L else 10L
  }

  if(unicode) {
    # Tortoise-shell brackets (U+2997/U+2998), heavy horizontal rule
    # (U+2501) filled/empty, heavy sub-cell (U+257A); ANSI 38;5;73 blue
    # and 37 gray, matching progress_manager.cpp.
    lhs = "\u2997"
    rhs = "\u2998"
    filled = "\u001b[38;5;73m\u2501\u001b[39m"
    partial_more = filled
    partial_less = "\u001b[37m\u257a\u001b[39m"
    empty = "\u001b[37m\u2501\u001b[39m"
  } else {
    lhs = "["
    rhs = "]"
    filled = "="
    partial_more = " "
    partial_less = " "
    empty = " "
  }

  draw = function(current) {
    frac = if(total > 0L) current / total else 1
    exact = frac * bar_width
    n_filled = min(as.integer(exact), bar_width)
    bar = strrep(filled, n_filled)
    if(n_filled < bar_width) {
      part = exact - n_filled
      if(part > 0) {
        bar = paste0(bar, if(part > 0.5) partial_more else partial_less)
        n_filled = n_filled + 1L
      }
    }
    if(n_filled < bar_width) {
      bar = paste0(bar, strrep(empty, bar_width - n_filled))
    }
    cat(sprintf(
      "\r%s: %s%s%s %d/%d (%.1f%%)", prefix, lhs, bar, rhs,
      current, total, 100 * frac
    ))
    utils::flush.console()
  }

  list(update = draw, close = function() cat("\n"))
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
                                    cores = 1L, base_seed = 1L,
                                    show_progress = FALSE) {
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

  # The bar follows display_progress (like the sampler's own bar), not the
  # advisory bgms.verbose flag. It redraws in place with a carriage return, so
  # it is drawn only in an interactive session; batch runs (scripts, R CMD
  # check, the test suite) rely on the one-line build announcement in
  # ggm_correction_table(). A parallel sweep forks the cells in batches of
  # `cores` and advances the bar as each batch of chains completes.
  n_cells = nrow(cells)
  draw_bar = show_progress && interactive()
  pb = NULL
  if(draw_bar) {
    pb = new_correction_progress(n_cells)
    on.exit(pb$close(), add = TRUE)
    pb$update(0L)
  }
  edens_cells = if(cores > 1L) {
    out = numeric(n_cells)
    done = 0L
    for(batch in split(seq_len(n_cells), ceiling(seq_len(n_cells) / cores))) {
      out[batch] = as.numeric(unlist(parallel::mclapply(
        batch, one_cell,
        mc.cores = cores, mc.preschedule = FALSE
      )))
      done = done + length(batch)
      if(draw_bar) {
        pb$update(done)
      }
    }
    out
  } else {
    vapply(seq_len(n_cells), function(k) {
      v = one_cell(k)
      if(draw_bar) {
        pb$update(k)
      }
      v
    }, numeric(1))
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
  precision_scale_prior = exponential_prior(eta = 1),
  n_grid = 120L, n_samples = 2000L, n_warmup = 500L, n_seeds = 3L,
  update_method = c("gibbs", "adaptive-metropolis"),
  cores = 1L, base_seed = 1L, show_progress = FALSE
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
    update_method = update_method, cores = cores, base_seed = base_seed,
    show_progress = show_progress
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
# correction_list_from_table (internal)
# ------------------------------------------------------------------
# Assemble the per-edge-prior correction curves handed to the C++
# chain. The beta-bernoulli update reads the whole-graph logC(theta)
# curve; the stochastic block updates read the slope vs local density
# plus the per-pair f(theta) curve on a uniform quadrature grid. On
# the mixed-MRF path is_continuous (0/1 per node, model order) marks
# the continuous block whose pairs carry the tilt; NULL means all
# nodes are continuous (the GGM path).
# ------------------------------------------------------------------
correction_list_from_table = function(table, edge_prior, is_continuous = NULL) {
  correction = list(theta = table$theta, logC = table$logC)
  if(identical(edge_prior, "Stochastic-Block")) {
    quad_theta = seq(0.0025, 0.9975, length.out = 200L)
    correction$fprime_density = table$fprime_density
    correction$fprime = table$fprime
    correction$quad_theta = quad_theta
    correction$quad_f = stats::approx(
      table$theta, table$f, quad_theta,
      rule = 2
    )$y
    if(!is.null(is_continuous)) {
      correction$is_continuous = as.integer(is_continuous)
    }
  }
  correction
}


# ------------------------------------------------------------------
# correction_interaction_prior / correction_scale_prior (internal)
# ------------------------------------------------------------------
# Rebuild the prior objects for the tilted prior sweep from a fit
# spec's flattened prior list. The interaction prior is NULL for slab
# families the prior sampler does not support (beta-prime).
# ------------------------------------------------------------------
correction_interaction_prior = function(prior) {
  switch(prior$interaction_prior_type,
    cauchy = cauchy_prior(scale = prior$pairwise_scale),
    normal = normal_prior(scale = prior$pairwise_scale),
    NULL
  )
}

correction_scale_prior = function(prior) {
  if(identical(prior$scale_prior_type, "exponential")) {
    exponential_prior(rate = prior$scale_rate)
  } else {
    gamma_prior(shape = prior$scale_shape, rate = prior$scale_rate)
  }
}


# ------------------------------------------------------------------
# ggm_edge_prior_correction (internal)
# ------------------------------------------------------------------
# Resolve whether a fit needs the normalizing-constant correction and
# get-or-build the table for its model cell. Applies to fits with
# edge selection and a hierarchical edge prior (Beta-Bernoulli or
# Stochastic-Block). The tilted prior sweep runs the same priors as
# the fit; the prior sampler does not support a beta-prime slab, so
# those fits keep the uncorrected updates with a warning.
#
# On the mixed-MRF path the determinant tilt acts on the continuous
# precision block alone, so the table is built for a GGM of dimension
# num_continuous and the correction list carries the continuous-block
# node mask (model order: discrete block first). With fewer than two
# continuous variables no pair is tilted and the plain conjugate
# updates are exact, so no correction applies.
#
# Returns the correction list for sample_ggm / sample_mixed_mrf, or
# NULL when no correction applies.
# ------------------------------------------------------------------
ggm_edge_prior_correction = function(prior, sampler, num_variables,
                                     num_continuous = num_variables) {
  if(!isTRUE(prior$edge_selection)) {
    return(NULL)
  }
  if(!prior$edge_prior %in% c("Beta-Bernoulli", "Stochastic-Block")) {
    return(NULL)
  }
  if(num_continuous < 2L) {
    return(NULL)
  }
  if(!prior$interaction_prior_type %in% c("cauchy", "normal")) {
    warning(
      "The ", prior$edge_prior, " updates are run without the ",
      "normalizing-constant correction: the tilted prior sampler supports ",
      "only cauchy_prior() and normal_prior() interaction priors.",
      call. = FALSE
    )
    return(NULL)
  }

  interaction_prior = correction_interaction_prior(prior)
  precision_scale_prior = correction_scale_prior(prior)

  # The build announcement follows the advisory bgms.verbose flag; the progress
  # bar follows display_progress (progress_type 0 is "none"), matching the
  # sampler's own bar.
  show_progress = is.null(sampler$progress_type) ||
    !identical(as.integer(sampler$progress_type), 0L)

  table = ggm_correction_table(
    p = num_continuous, delta = prior$delta,
    interaction_prior = interaction_prior,
    precision_scale_prior = precision_scale_prior,
    update_method = "gibbs",
    cores = sampler$cores,
    verbose = isTRUE(sampler$verbose),
    show_progress = show_progress
  )
  if(identical(prior$edge_prior, "Stochastic-Block") &&
    is.null(table$fprime)) {
    warning(
      "The Stochastic-Block updates are run without the ",
      "normalizing-constant correction: the slope curve is not resolvable ",
      "for this model cell (a single tilted pair).",
      call. = FALSE
    )
    return(NULL)
  }
  is_continuous = NULL
  if(num_continuous < num_variables) {
    is_continuous = c(
      rep(0L, num_variables - num_continuous),
      rep(1L, num_continuous)
    )
  }
  correction_list_from_table(table, prior$edge_prior, is_continuous)
}


# ------------------------------------------------------------------
# Shared cache-policy accessors (internal)
# ------------------------------------------------------------------
# Single source for the correction-cache option defaults so the
# correction-table build and the Z-ratio surface build agree.
# ------------------------------------------------------------------
correction_cache_enabled = function() {
  isTRUE(getOption("bgms.correction_table_cache", TRUE))
}

correction_cache_dir = function() {
  getOption("bgms.correction_cache_dir", tools::R_user_dir("bgms", which = "cache"))
}

# ------------------------------------------------------------------
# ggm_correction_table_key (internal)
# ------------------------------------------------------------------
# Disk-cache file name for one correction table: the model cell, the
# builder settings, and the package version.
#
# The version is part of the key on the same convention as the Z-ratio
# surface cache (zratio_surface_cache_key): the sweep that builds the
# table is code, so a release that changes it must not be served an
# earlier version's table out of the shared cache directory. Files
# under an older key are simply never read again.
#
# @param cell           Cell identity from ggm_correction_cell().
# @param n_grid,n_samples,n_warmup,n_seeds,update_method  Builder settings.
#
# Returns: the file name, including the .rds extension.
# ------------------------------------------------------------------
ggm_correction_table_key = function(cell, n_grid, n_samples, n_warmup,
                                    n_seeds, update_method) {
  sprintf(
    "ggm_ctable_v1_%s_q%d_delta%.8g_eta%.8g_%s_shape%.8g_g%d_ns%d_nw%d_sd%d_%s.rds",
    as.character(utils::packageVersion("bgms")),
    cell$q, cell$delta, cell$eta, cell$slab_family, cell$scale_shape,
    n_grid, n_samples, n_warmup, n_seeds, update_method
  )
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
  precision_scale_prior = exponential_prior(eta = 1),
  n_grid = 120L, n_samples = 2000L, n_warmup = 500L, n_seeds = 3L,
  update_method = c("gibbs", "adaptive-metropolis"),
  cores = 1L, base_seed = 1L, refresh = FALSE, verbose = FALSE,
  show_progress = FALSE
) {
  update_method = match.arg(update_method)
  if(is.null(delta)) {
    delta = 0.5 * log(p)
  }

  use_cache = correction_cache_enabled()
  cache_file = NULL
  if(use_cache) {
    cell = ggm_correction_cell(
      p, delta, interaction_prior, precision_scale_prior
    )
    key = ggm_correction_table_key(
      cell, n_grid, n_samples, n_warmup, n_seeds, update_method
    )
    cache_dir = correction_cache_dir()
    cache_file = file.path(cache_dir, key)
    if(!refresh && file.exists(cache_file)) {
      table = tryCatch(readRDS(cache_file), error = function(e) NULL)
      if(!is.null(table)) {
        return(table)
      }
    }
  }

  if(verbose) {
    message(
      "Building the edge-selection prior correction table (one-time for ",
      "this model size and prior; cached for later fits)."
    )
  }
  table = build_ggm_correction_table(
    p = p, delta = delta,
    interaction_prior = interaction_prior,
    precision_scale_prior = precision_scale_prior,
    n_grid = n_grid, n_samples = n_samples, n_warmup = n_warmup,
    n_seeds = n_seeds, update_method = update_method,
    cores = cores, base_seed = base_seed, show_progress = show_progress
  )

  if(use_cache) {
    dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(table, cache_file)
  }
  table
}

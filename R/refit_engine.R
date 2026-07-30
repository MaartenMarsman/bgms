# ==============================================================================
# Refit engine for prior_sensitivity_check(): warm-started refits at a grid of
# interaction slab scales.
# ==============================================================================
#
# A refit reuses the original fit's spec, changes the slab scale, shortens the
# warmup, and warm-starts each chain from the original fit's final per-chain
# state (parameters + edge indicators). Adaptation stays live during the short
# warmup so proposals re-tune to the new scale. See
# dev/audit/2026-07-28-refit-sensitivity-brief.md.
# ==============================================================================


# ------------------------------------------------------------------
# extract_warm_state
# ------------------------------------------------------------------
# Per-chain final state of a fit, in the exact vectorization the C++ sampler
# consumes (set_storage_vectorized_parameters + set_vectorized_indicator_
# parameters). For omrf the storage vector is c(main effects, all pairwise
# upper-triangle row-major); the indicator vector is the pairwise upper triangle.
#
# @param fit  A fitted bgms object (omrf) with edge selection.
#
# The returned object is sampler-agnostic: `parameters` is universal (every
# sampler consumes it), while the tuning fields are per method and optional --
# `step_sizes` and `inv_mass` are the NUTS metric (NULL/NA for a non-NUTS
# origin). A refit sampler that does not recognize a field simply ignores it.
# Warm starts use a dense all-edges-active graph, so no indicator state is
# carried (warm-starting a sparse graph would desynchronize the active-parameter
# dimension; see refit_at_scale).
#
# Returns: list(parameters = per-chain numeric vectors,
#               step_sizes  = per-chain NUTS step sizes (NA if none),
#               inv_mass    = per-chain NUTS diagonal metric (NULL if none)).
# ------------------------------------------------------------------
extract_warm_state = function(fit) {
  spec = get_fit_spec(fit)
  raw = get_raw_samples(fit)
  if(!identical(spec$model_type, "omrf")) {
    stop(
      "Warm-started refits are currently implemented for ordinal (omrf) ",
      "models only; got model type '", spec$model_type, "'."
    )
  }
  nchains = length(raw$pairwise)
  parameters = vector("list", nchains)
  for(c in seq_len(nchains)) {
    ni = nrow(raw$pairwise[[c]])
    # c(main, all pairwise) is the storage vector, row for the last iteration.
    parameters[[c]] = as.numeric(c(raw$main[[c]][ni, ], raw$pairwise[[c]][ni, ]))
  }
  # Per-chain final NUTS step size (NA for non-NUTS fits) and diagonal metric
  # (NULL for non-NUTS fits), for warm-starting a NUTS refit's tuning state. A
  # fit saved before these properties existed lacks them entirely; treat that as
  # absent tuning state (the refit falls back to a cold metric/step size).
  get_prop = function(name) {
    tryCatch(S7::prop(fit, name), error = function(e) {
      tryCatch(fit[[name]], error = function(e2) NULL)
    })
  }
  step_sizes = get_prop("refit_step_sizes")
  inv_mass = get_prop("refit_inv_mass")
  list(parameters = parameters, step_sizes = step_sizes, inv_mass = inv_mass)
}


# ------------------------------------------------------------------
# refit_at_scale
# ------------------------------------------------------------------
# One warm-started refit of `fit` at absolute slab scale `scale`, reusing the
# original spec with a fixed scale, a (short) warmup, a warm start, and a fresh
# seed. Returns a fitted bgms object.
#
# @param fit          Original fitted bgms object.
# @param scale        Absolute interaction slab scale for this refit.
# @param warm_state   Output of extract_warm_state(fit), or NULL for a cold refit.
# @param warmup       Warmup iterations for the refit.
# @param iter         Sampling iterations for the refit.
# @param seed         Seed for the refit.
# @param cores        Threads for the refit (chains run in parallel internally).
# @param sampler      Optional update-method override ("nuts",
#                     "adaptive-metropolis", "gibbs"); NULL keeps the fit's.
# @param show_progress  Show the sampler's native progress bar for this refit.
# @param scale_field  Name of the prior field the scale is written to:
#                     "pairwise_scale" for bgm(), "difference_scale" for
#                     bgmCompare().
# @param diagonal_rate  Raw rate for the precision diagonal at this scale, or
#                     NULL to leave the fit's own rate untouched. The caller
#                     owns this policy; see resolve_vary().
# ------------------------------------------------------------------
refit_at_scale = function(fit, scale, warm_state, warmup, iter, seed,
                          cores = 1L, sampler = NULL, show_progress = FALSE,
                          scale_field = "pairwise_scale",
                          diagonal_rate = NULL) {
  spec = get_fit_spec(fit)

  spec$prior[[scale_field]] = scale
  # The spec stores the diagonal rate already resolved to the raw frame, so a
  # new slab scale leaves it alone unless the caller asks otherwise.
  if(!is.null(diagonal_rate)) {
    spec$prior$scale_rate = diagonal_rate
  }

  if(!is.null(sampler)) {
    spec$sampler$update_method = sampler
  }
  spec$sampler$warmup = as.integer(warmup)
  spec$sampler$iter = as.integer(iter)
  spec$sampler$seed = as.integer(seed)
  spec$sampler$cores = as.integer(cores)
  # Render the sampler's own native display (a bar per chain, with the
  # warmup/sampling stage and ETA) for this refit when asked; otherwise stay
  # silent. The manager redraws its block in place, so a refit's chains do not
  # pile up as they run.
  spec$sampler$display_progress = if(show_progress) "per-chain" else "none"
  spec$sampler$progress_type = if(show_progress) 2L else 0L
  spec$sampler$progress_callback = NULL

  # Warm-start the continuous parameters only, on the default dense
  # (all-edges-active) start. Warm-starting a sparse indicator configuration
  # would desynchronize the sampler's active-parameter dimension from the full
  # storage dimension before stage-3c; a dense start matches the cold-start
  # dimension exactly, and the short warmup re-settles the graph (which shifts
  # with the scale anyway).
  #
  # A NUTS refit additionally carries the origin fit's per-chain tuning state
  # (step size + diagonal metric): dual averaging stays live for the step size
  # while the metric is held fixed, so the short warmup does not pay to
  # re-estimate the mass matrix. Non-NUTS refits ignore these fields.
  spec$initial_state = if(is.null(warm_state)) {
    NULL
  } else {
    state = list(parameters = warm_state$parameters)
    if(identical(spec$sampler$update_method, "nuts")) {
      if(!is.null(warm_state$step_sizes) && all(is.finite(warm_state$step_sizes))) {
        state$step_sizes = warm_state$step_sizes
      }
      if(!is.null(warm_state$inv_mass)) {
        state$inv_mass = warm_state$inv_mass
      }
    }
    state
  }

  raw = run_sampler(spec)
  build_output(spec, raw)
}


# ------------------------------------------------------------------
# resolve_vary
# ------------------------------------------------------------------
# Which prior a scale sweep moves on a model that has a prior on the precision
# diagonal. The slab scale and that diagonal are tied through the standardized
# frame (raw rate = eta / s), so a sweep of the slab either holds the raw rate
# fixed or holds eta fixed and lets the raw rate follow. The two answer
# different questions and neither is the other's approximation, so the mode is
# resolved once and named in the report rather than left to the fit's frame.
#
# @param spec  The original fit's spec.
# @param vary  "auto", "slab", or "slab-and-diagonal".
#
# Returns: list(
#   mode        The resolved mode, one of "slab", "slab-and-diagonal", or
#               "none" for a model with no precision diagonal.
#   eta         Standardized rate held fixed under "slab-and-diagonal", or NA.
#   standardized  Whether the fit itself specified the diagonal in the
#               standardized frame.
# )
# ------------------------------------------------------------------
resolve_vary = function(spec, vary) {
  vary = match.arg(vary, c("auto", "slab", "slab-and-diagonal"))
  prior = spec$prior
  has_diagonal = !is.null(prior$scale_rate) && !is.na(prior$scale_rate)
  if(!has_diagonal) {
    return(list(mode = "none", eta = NA_real_, standardized = FALSE))
  }
  standardized = !is.null(prior$scale_eta) && !is.na(prior$scale_eta)
  mode = if(identical(vary, "auto")) {
    if(standardized) "slab-and-diagonal" else "slab"
  } else {
    vary
  }
  # Under "slab-and-diagonal" the quantity held fixed is eta. A fit specified in
  # the raw frame carries none, so use the eta its own scale implies: the same
  # one-parameter family, anchored at the chosen scale.
  eta = if(identical(mode, "slab-and-diagonal")) {
    if(standardized) prior$scale_eta else prior$scale_rate * prior$pairwise_scale
  } else {
    NA_real_
  }
  list(mode = mode, eta = eta, standardized = standardized)
}


# ------------------------------------------------------------------
# vary_diagonal_rate
# ------------------------------------------------------------------
# Raw diagonal rate a refit at `scale` runs with, given a resolved `vary`.
# NULL leaves the fit's own rate in place, which is what "slab" and a model
# without a precision diagonal both want.
# ------------------------------------------------------------------
vary_diagonal_rate = function(vary, scale) {
  if(!identical(vary$mode, "slab-and-diagonal")) {
    return(NULL)
  }
  vary$eta / scale
}


# ------------------------------------------------------------------
# resolve_refit_sampler
# ------------------------------------------------------------------
# Resolve the requested refit sampler against the original fit's update method.
# "same-as-fit" inherits; otherwise the named method is used. Returns the method
# and whether to recommend NUTS for speed (an AM/Gibbs original inherited as-is).
# ------------------------------------------------------------------
resolve_refit_sampler = function(fit, refit_sampler) {
  fit_method = get_fit_spec(fit)$sampler$update_method
  method = if(identical(refit_sampler, "same-as-fit")) fit_method else refit_sampler
  valid = c("nuts", "adaptive-metropolis", "gibbs")
  if(!method %in% valid) {
    stop(
      "'refit_sampler' must be \"same-as-fit\" or one of ",
      paste(sQuote(valid), collapse = ", "), "; got ", sQuote(method), "."
    )
  }
  # Recommend NUTS only when the user inherited a slower same-as-fit method.
  recommend_nuts = identical(refit_sampler, "same-as-fit") && method != "nuts"
  list(method = method, recommend_nuts = recommend_nuts, fit_method = fit_method)
}


# ------------------------------------------------------------------
# refit_run_length
# ------------------------------------------------------------------
# Warmup/iter for a refit. A NUTS refit that carries the origin's metric needs
# only a short warmup (the validated w500/i1000 default); every other path
# inherits the original fit's run length (no warm-start speed-up, cost ~1x).
# Explicit `warmup`/`iter` always win.
# ------------------------------------------------------------------
refit_run_length = function(fit, method, warm_metric, warmup, iter) {
  s = get_fit_spec(fit)$sampler
  if(is.null(warmup)) {
    warmup = if(identical(method, "nuts") && warm_metric) 500L else s$warmup
  }
  if(is.null(iter)) {
    iter = if(identical(method, "nuts") && warm_metric) 1000L else s$iter
  }
  list(warmup = as.integer(warmup), iter = as.integer(iter))
}


# ------------------------------------------------------------------
# data_preferred_scale
# ------------------------------------------------------------------
# The scale the data prefer, read off the ORIGINAL fit's draws with no refit:
# s_hat is the root-mean-square of the included interactions, with an approximate
# 1/sqrt(2m) log-scale standard error (m = mean number of included edges). For
# the GGM the slab sits on Kyy = -0.5 * Omega, matching the reweighting path.
#
# Returns: list(s_hat, lo, hi, log_sd, m). NA-valued if no edge is ever included.
# ------------------------------------------------------------------
data_preferred_scale = function(fit) {
  spec = get_fit_spec(fit)
  raw = get_raw_samples(fit)
  slab_factor = if(identical(spec$model_type, "ggm")) -0.5 else 1
  theta = slab_factor * do.call(rbind, raw$pairwise)
  gamma = do.call(rbind, raw$indicator)
  incl = gamma == 1
  n_incl = sum(incl)
  if(n_incl == 0L) {
    return(list(
      s_hat = NA_real_, lo = NA_real_, hi = NA_real_,
      log_sd = NA_real_, m = 0
    ))
  }
  s_hat = sqrt(mean(theta[incl]^2))
  m = mean(rowSums(incl))
  log_sd = 1 / sqrt(2 * max(m, 1))
  list(
    s_hat = s_hat, lo = s_hat * exp(-1.96 * log_sd),
    hi = s_hat * exp(1.96 * log_sd), log_sd = log_sd, m = m
  )
}


# ------------------------------------------------------------------
# verdict_from_lbf
# ------------------------------------------------------------------
# Map a log10 inclusion Bayes factor to a verdict at threshold lthr = log10(t).
# ------------------------------------------------------------------
verdict_from_lbf = function(lbf, lthr) {
  out = rep("undecided", length(lbf))
  out[lbf >= lthr] = "presence"
  out[lbf <= -lthr] = "absence"
  out[is.na(lbf)] = NA_character_
  out
}


# ------------------------------------------------------------------
# refit_edge_stats
# ------------------------------------------------------------------
# Per-edge quantities from one refit: the Rao-Blackwellized inclusion
# probability (canonical), its log10 inclusion Bayes factor and verdict, the
# within-chain MCSE of the log10 BF (RB machinery, PR #182), the between-chain
# half-band (2 * SEM of the grand mean from chain spread, propagated to the log10
# BF), per-chain verdict unanimity, and the zero-flip mask. Edge order is the
# native row-major upper triangle.
# ------------------------------------------------------------------
refit_edge_stats = function(fit, evidence_threshold) {
  raw = get_raw_samples(fit)
  enm = raw$parameter_names$indicator %||% raw$parameter_names$pairwise
  M = length(raw$rb_inclusion)
  pc = sapply(raw$rb_inclusion, colMeans) # n_edges x M
  if(is.null(dim(pc))) pc = matrix(pc, ncol = M)
  pbar = rowMeans(pc)

  prm = extract_prior_inclusion_probabilities(fit)
  nv = nrow(prm)
  prior_q = prm[upper.tri(prm)][order_upper_tri_rowmajor(nv)]
  prior_odds = prior_q / (1 - prior_q)

  ind = fit@posterior_summary_indicator
  mcse_pip = ind[, "mcse"]
  lthr = log10(evidence_threshold)

  lbf = function(p) log10((p / (1 - p)) / prior_odds)
  dfac = 1 / (log(10) * pmax(pbar * (1 - pbar), 1e-6)) # d log10 BF / d p
  lbf_bar = lbf(pbar)
  mcse_lbf = mcse_pip * dfac
  se_between = apply(pc, 1, stats::sd) / sqrt(M)
  band_half = 2 * se_between * dfac

  verdict = verdict_from_lbf(lbf_bar, lthr)
  vc = apply(pc, 2, function(p) verdict_from_lbf(lbf(p), lthr)) # n_edges x M
  if(is.null(dim(vc))) vc = matrix(vc, ncol = M)
  unanimous = apply(vc, 1, function(v) length(unique(v[!is.na(v)])) <= 1L)
  zeroflip = apply(pc, 1, function(p) all(p <= 0) || all(p >= 1))

  list(
    edge = enm, pip = pbar, prior_odds = prior_odds, lbf = lbf_bar,
    verdict = verdict, mcse_lbf = mcse_lbf, band_half = band_half,
    unanimous = unanimous, zeroflip = zeroflip
  )
}


# ------------------------------------------------------------------
# refit_convergence_gate
# ------------------------------------------------------------------
# The refit-level gate: a global-failure detector on BULK quantities. It gates
# on the MEDIAN continuous split-Rhat, not the max: the max is routinely
# inflated by a few flat-likelihood directions (near-unidentified extreme-
# category thresholds, and interactions of near-saturated edges) whose Rhat is
# an identification artifact, not a verdict-corrupting failure. Global slow
# mixing (e.g. an undermixed adaptive-Metropolis original) lifts the median and
# the RB-inclusion median together and is caught; a couple of flat directions
# are not. The RB-inclusion tail lives at the edge level (per-chain verdict
# unanimity), because warm starts weaken split-Rhat's between-chain signal. A
# refit failing this gate is reported as unusable, never silently pooled.
#
#   median continuous split-Rhat < 1.01 (main effects + pairwise),
#   bulk RB-inclusion median Rhat < 1.01,
#   E-BFMI > 0.3 and first/second-half energy variance ratio < 2 (NUTS).
# The energy-slope warmup heuristic is deliberately NOT gated (it flags healthy
# cold fits; see the diagnostics doctrine). The max continuous Rhat and the
# smallest Rao-Blackwellized inclusion ESS are reported for transparency.
# ------------------------------------------------------------------
refit_convergence_gate = function(fit) {
  rhats = numeric(0)
  esss = numeric(0)
  mn = fit@posterior_summary_main
  if(!is.null(mn)) {
    rhats = c(rhats, mn[, "Rhat"])
    esss = c(esss, mn[, "n_eff"])
  }
  pw = fit@posterior_summary_pairwise
  if(!is.null(pw)) {
    rhats = c(rhats, pw[, "Rhat"])
    esss = c(esss, pw[, "n_eff"])
  }
  ind = fit@posterior_summary_indicator

  rhat_cont = stats::median(rhats, na.rm = TRUE)
  rhat_cont_max = max(rhats, na.rm = TRUE)
  ess_cont = min(esss, na.rm = TRUE)
  ess_incl = min(ind[, "n_eff"], na.rm = TRUE)
  rb_med_rhat = stats::median(ind[, "Rhat"], na.rm = TRUE)

  wc = tryCatch(fit@nuts_diag$warmup_check, error = function(e) NULL)
  ebfmi = if(is.null(wc)) {
    Inf
  } else {
    min(c(wc$ebfmi_first_half, wc$ebfmi_second_half), na.rm = TRUE)
  }
  var_ratio = if(is.null(wc)) 0 else max(wc$var_ratio, na.rm = TRUE)

  usable = is.finite(rhat_cont) && rhat_cont < 1.01 &&
    rb_med_rhat < 1.01 && ebfmi > 0.3 && var_ratio < 2
  list(
    usable = usable, rhat_cont = rhat_cont, rhat_cont_max = rhat_cont_max,
    ess_cont = ess_cont, ess_incl = ess_incl, rb_med_rhat = rb_med_rhat,
    min_ebfmi = ebfmi, max_var_ratio = var_ratio
  )
}


# ------------------------------------------------------------------
# gate_failure_reason
# ------------------------------------------------------------------
# Plain-language reason a convergence gate failed, from a gate result. The
# criteria are checked in the gate's own order; the first failing one is
# reported. Returns NA when the gate passed.
# ------------------------------------------------------------------
gate_failure_reason = function(g) {
  if(!is.finite(g$rhat_cont) || g$rhat_cont >= 1.01) {
    return("the parameter chains have not converged (split R-hat above 1.01)")
  }
  if(g$rb_med_rhat >= 1.01) {
    return("the edge-inclusion chains have not converged (R-hat above 1.01)")
  }
  if(g$min_ebfmi <= 0.3) {
    return("the sampler explored the posterior poorly (low energy efficiency)")
  }
  if(g$max_var_ratio >= 2) {
    return("the warmup did not settle (the sampling energy kept drifting)")
  }
  NA_character_
}

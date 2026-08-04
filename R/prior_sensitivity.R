# ==============================================================================
# prior_sensitivity_check: scale-robustness of edge inclusion verdicts
# ==============================================================================
#
# Refit the model at a small grid of interaction slab scales (multipliers of the
# fit's chosen scale) and classify every edge at every scale, with a Monte Carlo
# wobble guard calibrated from a replicate refit at the chosen scale. Warm starts
# and the refit engine live in refit_engine.R. See the package vignette
# "Checking prior sensitivity" and Bartos et al. (arXiv 2604.21596).
#
# The previous single-fit conditional-density reweighting curve was removed here:
# on real network data the scale posterior collapses and its importance ESS is ~1
# across the advertised band (phase-B review, 2026-07-28). Git history and
# dev/audit/2026-07-28-phase-b-review.md are the archive if it is ever revived.
# ==============================================================================


# ------------------------------------------------------------------
# verdict_from_bf
# ------------------------------------------------------------------
# Map an inclusion Bayes factor to an evidence verdict.
#
# @param bf         Numeric inclusion Bayes factor(s).
# @param threshold  Upper evidence threshold; lower is 1 / threshold.
#
# Returns: character vector of "presence" / "absence" / "undecided".
# ------------------------------------------------------------------
verdict_from_bf = function(bf, threshold) {
  out = rep("undecided", length(bf))
  out[bf >= threshold] = "presence"
  out[bf <= 1 / threshold] = "absence"
  out[is.na(bf)] = NA_character_
  out
}


#' @title Prior Sensitivity of Inclusion Verdicts
#'
#' @description
#' Recovers each indicator's continuous inclusion-Bayes-factor curve
#' \eqn{\mathrm{BF}_e(s)} across the scale \eqn{s} of the slab that gates it,
#' and classifies every verdict trajectory along it. The curve is anchored at a
#' handful of fixed-scale fits --- the chosen-scale anchor is the original fit
#' itself --- and filled in between anchors by importance reweighting. A
#' built-in Monte Carlo noise guard, calibrated from a repeated refit at one
#' anchor, keeps the check from crying wolf on boundary cases.
#'
#' On a \code{\link{bgm}()} fit with edge selection the unit is the edge
#' indicator and the swept prior is the interaction slab. On a
#' \code{\link{bgmCompare}()} fit with difference selection the unit is the
#' difference indicator and the swept prior is the difference slab
#' (\code{difference_scale}), which covers the pairwise and the main-effect
#' difference families alike, since \code{bgmCompare()} gives them one scale.
#'
#' @details
#' \strong{The anchored curve.} The model is refit at the non-unit
#' \code{anchors} (multipliers of the chosen scale \eqn{s_0}; the \code{1x}
#' anchor is the original fit and is never refit). A fit at fixed anchor
#' scale \eqn{s_a} is reweighted to a nearby scale \eqn{s} with per-draw
#' slab-density ratios over the currently included edges; the likelihood
#' cancels, so no refit is needed between anchors. The posterior-density-ratio
#' reweighting identity of Bartos et al. (2026) is used locally around each
#' refit anchor; the anchored construction itself is specific to \pkg{bgms}.
#' Each point of a dense
#' log-spaced display grid pools every anchor that clears \code{ess_floor}
#' there, weighting each anchor's inclusion-probability estimate by its
#' inverse variance (\eqn{\mathrm{ESS} / (p(1-p))}); the pooling is on the
#' inclusion-probability scale and is then transformed to the natural log
#' Bayes factor, which keeps the curve continuous across
#' anchor hand-offs and finite at capped edges. Points where no anchor
#' clears the floor are \code{NA} rather than extrapolated, and
#' non-overlapping anchor radii trigger a warning to add anchors; log-spaced
#' default anchors make the radii overlap. Exactness is kept off the pooled
#' curve and on the anchor fits themselves: the per-anchor verdict columns
#' and every chosen-scale quantity are read straight from each fit's own
#' Rao-Blackwellized statistics, so the \code{1x} column is exactly the
#' original fit's reported analysis. Between the anchors the curve is therefore
#' importance-reweighted rather than refit, and can deviate from a refit at that
#' scale by up to roughly \code{0.01} in inclusion probability at the
#' extrapolation ends; the anchors themselves -- including the \code{1x} anchor,
#' which is the user's own fit -- are exact.
#'
#' \strong{Warm starts.} For ordinal (omrf) fits each refit starts from the
#' original fit's per-chain final state, and a NUTS refit additionally carries
#' the adapted step size and diagonal mass matrix, so a short warmup suffices
#' and the whole check costs about one original fit. Because the warm starts sit
#' near the chosen-scale posterior, cross-chain dispersion is reduced by
#' construction, which weakens split-\eqn{\hat R} as a between-chain diagnostic;
#' the refit gate therefore reads the \emph{median} split-\eqn{\hat R} over the
#' continuous parameters and over the Rao-Blackwellized inclusion probabilities
#' together with the NUTS energy diagnostics, and reports the smallest
#' Rao-Blackwellized inclusion \code{n_eff} as \code{$grid$inclusion_ess_min};
#' the per-chain agreement that guards an individual verdict is the edge-level
#' sufficiency check below. Continuous (GGM) and
#' mixed-MRF fits refit cold (full warmup), costing about one fit per scale.
#'
#' \strong{The mover rule.} An edge is flagged scale-sensitive only if its
#' verdict differs somewhere along the curve \emph{and} its natural log
#' Bayes-factor change across scales exceeds
#' \code{max(tolerance, 2 * MCSE, noise)}, with the per-point MCSE from the
#' chain-level spread of the reweighted estimate (so it carries both the
#' between-chain and the importance-sampling uncertainty) and the noise band
#' the 95th percentile of the spread between one anchor refit and its
#' repeat, over threshold-relevant edges
#' (\eqn{|\log \mathrm{BF}| \le 3 \log 10 \approx 6.91}; near-saturated edges would
#' inflate it). An edge that is threshold-relevant in one of the two refits and
#' saturated in the other has a censored rather than an infinite spread, and is
#' left out of that percentile; the printed report counts them. When no
#' threshold-relevant edge has a measurable spread the noise band is \code{NA}
#' and the mover rule falls back to \code{max(tolerance, 2 * MCSE)}.
#' The \code{$edges$mover} column stores \code{stable},
#' \code{indistinguishable-from-wobble}, or \code{moved-beyond-wobble}; the
#' printed report shows the same categories in plain language ("robust",
#' "changed, within run-to-run noise", "changed, beyond run-to-run noise",
#' with edges failing the sufficiency check below printed as "not
#' certifiable"). A bare verdict flip inside the replicate noise is never
#' reported as a move.
#'
#' \strong{Edge-level sufficiency.} An edge whose per-chain verdicts disagree, or
#' whose between-chain-inflated log-BF band straddles a verdict
#' threshold, is marked \code{insufficient} at that scale: the refit cannot
#' certify its verdict. This errs toward caution --- disagreement widens the
#' band rather than vanishing into a pooled estimate.
#'
#' \strong{Two-level gate.} Each refit passes a refit-level gate (continuous
#' split-\eqn{\hat R} below 1.01, bulk RB-inclusion median \eqn{\hat R} below
#' 1.01, E-BFMI and energy variance-ratio for NUTS) before its verdicts are
#' used; a refit that fails is reported as unusable rather than silently pooled.
#'
#' @param bgms_object A fitted \code{bgms} object from \code{\link{bgm}()} run
#'   with \code{edge_selection = TRUE}, or a \code{bgmCompare} object from
#'   \code{\link{bgmCompare}()} run with \code{difference_selection = TRUE}.
#' @param anchors Numeric vector of positive anchor multipliers of the chosen
#'   scale. Default \code{c(0.4, 0.63, 1, 1.6, 2.5)} (log-spaced, so adjacent
#'   anchors' usable reweighting radii overlap). The multiplier \code{1} is
#'   always included and is the original fit; each other anchor is one warm
#'   refit.
#' @param evidence_threshold Positive numeric. Inclusion Bayes factor threshold
#'   for a presence verdict; \code{1 / evidence_threshold} is the absence
#'   threshold. Default: \code{10}.
#' @param vary One of \code{"auto"} (default), \code{"slab"}, or
#'   \code{"slab-and-diagonal"}. Only relevant for models with a prior on the
#'   precision diagonal (continuous and mixed); ignored for discrete ones, which
#'   have none. The slab scale \eqn{s} and the diagonal rate are tied through
#'   the standardized frame (raw rate \eqn{= \eta / s}), so a sweep of \eqn{s}
#'   has to hold one of the two fixed. \code{"slab"} holds the raw diagonal rate
#'   at the fitted value and moves the interaction prior alone, answering how
#'   much the verdicts depend on how wide an edge is allowed to be.
#'   \code{"slab-and-diagonal"} holds \eqn{\eta} fixed and lets the raw rate
#'   follow, answering how much they depend on the overall prior scale with its
#'   shape held fixed. \code{"auto"} follows the frame the fit itself used:
#'   \code{"slab-and-diagonal"} when the diagonal prior was given as \code{eta},
#'   \code{"slab"} when it was given as a raw \code{rate}. The resolved mode is
#'   named in the printed report.
#' @param refit_sampler One of \code{"same-as-fit"} (default; inherit the
#'   original fit's update method) or an explicit \code{"nuts"},
#'   \code{"adaptive-metropolis"}, or \code{"gibbs"}. NUTS refits of an ordinal
#'   fit carry the adapted metric and run fastest; when an inherited slower
#'   sampler makes the check cost more than about a minute, a message suggests
#'   the switch. Refitting with a different sampler than the original fit is
#'   sound because every comparison the check makes runs refit-against-refit
#'   under one identical pipeline; the original fit never enters a comparison.
#' @param iter,warmup Integer sampling and warmup iterations per refit, or
#'   \code{NULL} (default) to use the validated short schedule for warm NUTS
#'   refits and inherit the original fit's schedule otherwise.
#' @param tolerance Numeric. Minimum natural log BF change across scales for
#'   a verdict flip to count as a move, before the MCSE and noise floors.
#'   Default: \code{0.5 * log(10)} (about 1.15, the same decision boundary the
#'   former 0.5 carried in \eqn{\log_{10}} units).
#' @param ess_floor Positive numeric. Minimum pooled importance effective
#'   sample size for a curve point to be reported; points below it are
#'   \code{NA}. Default: \code{400}.
#' @param include_preferred_scale Logical. Add an extra anchor at the
#'   data-preferred scale \eqn{\hat s}. Default: \code{FALSE}.
#' @param cores Integer thread count for each refit's chains. Default: the
#'   original fit's core count.
#' @param seed Integer base seed for the refits. Default: \code{1}.
#' @param keep_fits Logical. Retain the full refit objects in the result (for
#'   power users); the default keeps only per-scale summaries. Default:
#'   \code{FALSE}.
#' @param verbose Logical. If \code{TRUE}, print each internal refit's raw
#'   sampler notes (energy, tree-depth, dropped chains) live as it runs. By
#'   default these are captured, not printed: the convergence gate adjudicates
#'   them, a failed anchor is reported once in plain language, and the raw
#'   text stays available in \code{$refit_diagnostics}. Default: \code{FALSE}.
#'
#' @return An object of class \code{"bgms_prior_sensitivity"}: a list with the
#'   per-edge \code{edges} table (chosen-scale verdict from the original fit,
#'   per-anchor verdict columns, dense-grid stability range, mover category,
#'   an \code{insufficient} flag with its two subcauses
#'   \code{insufficient_noisy} (Bayes factor within Monte Carlo error of a
#'   threshold) and \code{insufficient_disagree} (chains disagree on the
#'   verdict)), a \code{grid} data frame (one row per anchor fit and
#'   the replicate, with convergence gates), the \code{multipliers} display
#'   grid with the \code{log_bf}, \code{log_bf_mcse}, and \code{verdict}
#'   curve matrices (grid-by-edge, \code{NA} where masked), the \code{curve}
#'   bookkeeping (per-point importance ESS, anchor used, \code{ess_floor},
#'   per-point chain unanimity), the \code{wobble} noise yardstick,
#'   \code{refit_diagnostics} (the captured raw sampler notes per anchor and
#'   the replicate), the data-\code{preferred_scale} (\code{NA} where the swept
#'   prior is not the one on the pairwise interactions), the \code{unit} the
#'   check reported on, the resolved \code{vary} mode, and the settings used.
#'
#' @seealso \code{\link{bgm}()}, \code{\link{extract_posterior_inclusion_probabilities}()}
#' @family diagnostics
#' @references \insertRef{bartos2026}{bgms}
#' @export
#' @examples
#' \donttest{
#' fit = bgm(Wenchuan[, 1:6], chains = 2)
#' ps = prior_sensitivity_check(fit)
#' ps
#' plot(ps)
#' }
prior_sensitivity_check = function(bgms_object,
                                   anchors = c(0.4, 0.63, 1, 1.6, 2.5),
                                   evidence_threshold = 10,
                                   vary = c("auto", "slab", "slab-and-diagonal"),
                                   refit_sampler = "same-as-fit",
                                   iter = NULL,
                                   warmup = NULL,
                                   tolerance = 0.5 * log(10),
                                   ess_floor = 400,
                                   include_preferred_scale = FALSE,
                                   cores = NULL,
                                   seed = 1L,
                                   keep_fits = FALSE,
                                   verbose = FALSE) {
  UseMethod("prior_sensitivity_check")
}

#' @inheritParams prior_sensitivity_check
#' @exportS3Method
#' @noRd
prior_sensitivity_check.bgms = function(bgms_object,
                                        anchors = c(0.4, 0.63, 1, 1.6, 2.5),
                                        evidence_threshold = 10,
                                        vary = c("auto", "slab", "slab-and-diagonal"),
                                        refit_sampler = "same-as-fit",
                                        iter = NULL,
                                        warmup = NULL,
                                        tolerance = 0.5 * log(10),
                                        ess_floor = 400,
                                        include_preferred_scale = FALSE,
                                        cores = NULL,
                                        seed = 1L,
                                        keep_fits = FALSE,
                                        verbose = FALSE) {
  spec = get_fit_spec(bgms_object)
  if(is.null(spec) || !isTRUE(spec$prior$edge_selection)) {
    stop(
      "prior_sensitivity_check() needs edge selection. Refit with ",
      "edge_selection = TRUE."
    )
  }
  prior_sensitivity_engine(bgms_object,
    unit = single_network_unit(spec),
    anchors = anchors, evidence_threshold = evidence_threshold, vary = vary,
    refit_sampler = refit_sampler, iter = iter, warmup = warmup,
    tolerance = tolerance, ess_floor = ess_floor,
    include_preferred_scale = include_preferred_scale, cores = cores,
    seed = seed, keep_fits = keep_fits, verbose = verbose
  )
}

#' @inheritParams prior_sensitivity_check
#' @exportS3Method
#' @noRd
prior_sensitivity_check.bgmCompare = function(bgms_object,
                                              anchors = c(0.4, 0.63, 1, 1.6, 2.5),
                                              evidence_threshold = 10,
                                              vary = c("auto", "slab", "slab-and-diagonal"),
                                              refit_sampler = "same-as-fit",
                                              iter = NULL,
                                              warmup = NULL,
                                              tolerance = 0.5 * log(10),
                                              ess_floor = 400,
                                              include_preferred_scale = FALSE,
                                              cores = NULL,
                                              seed = 1L,
                                              keep_fits = FALSE,
                                              verbose = FALSE) {
  spec = get_fit_spec(bgms_object)
  if(is.null(spec) || !isTRUE(spec$prior$difference_selection)) {
    stop(
      "prior_sensitivity_check() needs difference selection. Refit with ",
      "bgmCompare(difference_selection = TRUE); without it every difference ",
      "is in the model and there is no inclusion Bayes factor to trace."
    )
  }
  if(is.na(difference_prior_inclusion(bgms_object))) {
    stop(
      "A stochastic-block difference prior has no single marginal inclusion ",
      "probability, so the curve would report posterior odds rather than ",
      "Bayes factors. Refit with bernoulli_prior() or beta_bernoulli_prior() ",
      "to trace the difference scale."
    )
  }
  prior_sensitivity_engine(bgms_object,
    unit = difference_unit(spec),
    anchors = anchors, evidence_threshold = evidence_threshold, vary = vary,
    refit_sampler = refit_sampler, iter = iter, warmup = warmup,
    tolerance = tolerance, ess_floor = ess_floor,
    include_preferred_scale = include_preferred_scale, cores = cores,
    seed = seed, keep_fits = keep_fits, verbose = verbose
  )
}


# ------------------------------------------------------------------
# single_network_unit / difference_unit
# ------------------------------------------------------------------
# What the check sweeps and what it reports on. bgm() traces edge indicators
# across the interaction slab scale; bgmCompare() traces difference indicators
# across the difference slab scale. Everything between -- the anchored curve,
# the convergence gate, the wobble yardstick, the reporting -- is shared.
#
# Returns: list(
#   scale_field  Prior field the refits rescale.
#   chosen_scale The fit's own value of it.
#   noun         Singular name of one reported unit.
#   nouns        Its plural.
#   headline     The question the printed report opens with.
#   preferred    Whether a data-preferred scale can be estimated.
# )
# ------------------------------------------------------------------
single_network_unit = function(spec) {
  list(
    scale_field = "pairwise_scale",
    chosen_scale = spec$prior$pairwise_scale,
    noun = "edge", nouns = "edges",
    headline = "are the edge verdicts robust to the slab scale?",
    preferred = TRUE
  )
}

difference_unit = function(spec) {
  list(
    scale_field = "difference_scale",
    chosen_scale = spec$prior$difference_scale,
    noun = "difference", nouns = "differences",
    headline = "are the difference verdicts robust to the difference scale?",
    preferred = FALSE
  )
}


# ------------------------------------------------------------------
# prior_sensitivity_engine
# ------------------------------------------------------------------
# The shared anchored-curve machinery behind both methods. `unit` names the
# prior field to sweep and the vocabulary of the report; see single_network_unit
# and difference_unit. Every other argument is the method's, unchanged.
# ------------------------------------------------------------------
prior_sensitivity_engine = function(bgms_object,
                                    unit,
                                    anchors,
                                    evidence_threshold,
                                    vary,
                                    refit_sampler,
                                    iter,
                                    warmup,
                                    tolerance,
                                    ess_floor,
                                    include_preferred_scale,
                                    cores,
                                    seed,
                                    keep_fits,
                                    verbose) {
  if(!is.numeric(evidence_threshold) || length(evidence_threshold) != 1L ||
    evidence_threshold <= 1) {
    stop("'evidence_threshold' must be a single number greater than 1.")
  }
  if(!is.numeric(anchors) || any(anchors <= 0)) {
    stop("'anchors' must be positive multipliers of the chosen scale.")
  }
  if(!is.numeric(ess_floor) || length(ess_floor) != 1L || ess_floor <= 0) {
    stop("'ess_floor' must be a single positive number.")
  }

  spec = get_fit_spec(bgms_object)
  chosen_scale = unit$chosen_scale
  lthr = log(evidence_threshold)
  if(is.null(cores)) cores = spec$sampler$chains

  # Which prior the sweep moves; named in the printed report.
  vary = resolve_vary(spec, vary)

  # --- Sampler resolution -----------------------------------------------------
  rs = resolve_refit_sampler(bgms_object, refit_sampler)

  # --- Warm state (omrf) or cold refits (other models) ------------------------
  is_omrf = identical(spec$model_type, "omrf")
  warm_state = if(is_omrf) extract_warm_state(bgms_object) else NULL
  warm_metric = is_omrf && identical(rs$method, "nuts") &&
    !is.null(warm_state$inv_mass)

  # --- Data-preferred scale (no refit) ----------------------------------------
  # Only defined where the swept prior is the one on the pairwise interactions;
  # the difference slab has no comparable plug-in estimate.
  preferred = if(isTRUE(unit$preferred)) {
    data_preferred_scale(bgms_object)
  } else {
    list(s_hat = NA_real_, lo = NA_real_, hi = NA_real_)
  }

  # --- Anchor set: 1x is the original fit, never refit ------------------------
  anchors = sort(unique(c(anchors, 1)))
  if(include_preferred_scale && is.finite(preferred$s_hat)) {
    anchors = sort(unique(c(anchors, preferred$s_hat / chosen_scale)))
  }
  n_anchor = length(anchors)
  s0_idx = which(anchors == 1)
  non_unit = setdiff(seq_len(n_anchor), s0_idx)
  if(length(non_unit) == 0L) {
    stop("'anchors' must include at least one multiplier other than 1.")
  }
  # The replicate refit that calibrates run-to-run noise sits at the non-unit
  # anchor nearest 1.6x, so its noise is measured at the refit settings.
  rep_anchor = non_unit[which.min(abs(log(anchors[non_unit] / 1.6)))]
  rl = refit_run_length(bgms_object, rs$method, warm_metric, warmup, iter)

  # One refit per non-unit anchor, plus the replicate.
  jobs = data.frame(
    multiplier = c(anchors[non_unit], anchors[rep_anchor]),
    replicate = c(rep(FALSE, length(non_unit)), TRUE)
  )
  refit_cores = min(as.integer(cores), spec$sampler$chains)

  refits = vector("list", nrow(jobs))
  walls = numeric(nrow(jobs))
  # Per-refit sampler chatter (energy/tree-depth notes, dropped-chain warnings)
  # is captured, not printed: a user watching the check cannot act on per-refit
  # alarms mid-run, and the convergence gate below re-tests everything they
  # gesture at. The captured text stays inspectable in $refit_diagnostics, and
  # verbose = TRUE re-enables live printing.
  refit_notes = vector("list", nrow(jobs))
  run_one_refit = function(i, show_progress = FALSE) {
    s = jobs$multiplier[i] * chosen_scale
    refit_at_scale(
      bgms_object,
      scale = s,
      warm_state = warm_state, warmup = rl$warmup, iter = rl$iter,
      seed = seed + i, cores = refit_cores, sampler = rs$method,
      show_progress = show_progress,
      scale_field = unit$scale_field,
      diagonal_rate = vary_diagonal_rate(vary, s)
    )
  }
  # Each anchor refit (the check's only real cost) renders the sampler's own
  # native per-chain display, with the whole refit grid announced up front so
  # the user knows what is coming. A refit's own cat() diagnostics (NUTS
  # settle-in notes and the like) are the check's to adjudicate, not the
  # sampler's to print: suppress them through bgms.verbose unless the user asked
  # to see them, which leaves only the progress bar (the option does not gate
  # it). Under verbose = TRUE the bar yields and the captured notes are reported.
  if(!verbose) {
    old_bgms_verbose = options(bgms.verbose = FALSE)
    on.exit(options(old_bgms_verbose), add = TRUE)
  }
  show_bar = interactive() && !verbose
  if(show_bar) {
    refit_scales = sprintf("%.2gx", jobs$multiplier[!jobs$replicate])
    rep_scale = sprintf("%.2gx", jobs$multiplier[jobs$replicate][1])
    n_chains = length(get_raw_samples(bgms_object)$pairwise)
    message(sprintf(
      paste0(
        "Refitting at %d scale%s (%s), plus a %s repeat for the noise band.\n",
        "Each refit runs %d chain%s; warmup and sampling are shown per chain."
      ),
      length(refit_scales), if(length(refit_scales) == 1L) "" else "s",
      paste(refit_scales, collapse = ", "), rep_scale,
      n_chains, if(n_chains == 1L) "" else "s"
    ))
  }
  for(i in seq_len(nrow(jobs))) {
    t0 = Sys.time()
    warns_i = character(0)
    fit_i = NULL
    absorb_warning = function(w) {
      warns_i <<- c(warns_i, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
    if(show_bar) {
      # The native display writes to the console, so it must not be captured;
      # announce which refit is starting, then let its per-chain bars animate.
      lbl = if(jobs$replicate[i]) {
        sprintf("%.2gx scale (noise repeat)", jobs$multiplier[i])
      } else {
        sprintf("%.2gx scale", jobs$multiplier[i])
      }
      message(sprintf("\n[%d/%d] %s", i, nrow(jobs), lbl))
      fit_i = withCallingHandlers(
        run_one_refit(i, show_progress = TRUE),
        warning = absorb_warning
      )
      notes = warns_i
    } else if(verbose) {
      # The check reports these notes itself, so capture them rather than let
      # the refit print them mid-run.
      cat_i = utils::capture.output({
        fit_i = withCallingHandlers(run_one_refit(i), warning = absorb_warning)
      })
      notes = c(cat_i[nzchar(cat_i)], warns_i)
    } else {
      # Diagnostics are suppressed via bgms.verbose above; nothing to capture.
      fit_i = withCallingHandlers(run_one_refit(i), warning = absorb_warning)
      notes = warns_i
    }
    refits[[i]] = fit_i
    refit_notes[[i]] = notes
    walls[i] = as.numeric(Sys.time() - t0, units = "secs")
    if(verbose && length(notes) > 0) {
      cat(sprintf("[%.2gx refit] sampler notes:\n", jobs$multiplier[i]))
      cat(paste0("  ", notes), sep = "\n")
    }
    # Suggest the faster sampler only when the cost is material (projected
    # total above ~60 s); see the man page for why the switch is sound.
    if(i == 1L && rs$recommend_nuts && walls[1] * nrow(jobs) > 60) {
      message(
        "Refits use the fit's ", rs$method, " sampler (~",
        round(walls[1]), " s each); refit_sampler = \"nuts\" is usually faster."
      )
    }
  }

  # --- Anchor fits: the 1x anchor is the original fit -------------------------
  anchor_fits = vector("list", n_anchor)
  anchor_fits[[s0_idx]] = bgms_object
  for(k in seq_along(non_unit)) anchor_fits[[non_unit[k]]] = refits[[k]]
  rep_fit = refits[[length(refits)]]

  gates = lapply(anchor_fits, refit_convergence_gate)
  rep_gate = refit_convergence_gate(rep_fit)
  usable = vapply(gates, `[[`, logical(1), "usable")

  # Adjudicate the captured refit chatter: the gate is the arbiter. A failed
  # anchor is reported once, in plain voice (which anchor, which criterion,
  # what the curve loses); a passed anchor's notes are absorbed (settle-in
  # transients the gate already cleared), inspectable in $refit_diagnostics.
  for(k in seq_along(non_unit)) {
    a = non_unit[k]
    if(!usable[a]) {
      message(sprintf(
        "The %.2gx-scale refit did not converge: %s. The curve does not cover the scales nearest %.2gx.",
        anchors[a], gate_failure_reason(gates[[a]]), anchors[a]
      ))
    }
  }
  if(!rep_gate$usable) {
    message(sprintf(
      "The repeated %.2gx-scale refit did not converge: %s. The run-to-run noise band could not be measured, so the mover rule falls back to the tolerance and Monte Carlo error floors.",
      anchors[rep_anchor], gate_failure_reason(rep_gate)
    ))
  }
  if(!usable[s0_idx]) {
    message(
      "The original fit does not pass the convergence gate; its chosen-scale ",
      "verdicts are still reported (they are the analysis under check), but ",
      "read the whole curve with caution."
    )
    usable[s0_idx] = TRUE
  }

  # Each anchor fit's own RB statistics carry the exact per-anchor verdicts.
  # Exactness lives here, not on the pooled curve: verdict columns and every
  # chosen-scale quantity are read straight from the fit that produced them.
  anchor_stats = lapply(anchor_fits, refit_edge_stats,
    evidence_threshold = evidence_threshold
  )
  s0 = anchor_stats[[s0_idx]]
  edge_names = s0$edge
  n_edges = length(edge_names)
  prior_odds = s0$prior_odds
  # Exact per-anchor verdicts (n_anchor x n_edge); a gate-failed anchor is NA.
  anchor_verdict = t(vapply(seq_len(n_anchor), function(a) {
    if(usable[a]) anchor_stats[[a]]$verdict else rep(NA_character_, n_edges)
  }, character(n_edges)))

  # --- Display grid and the stitched curve ------------------------------------
  grid_mult = exp(seq(log(min(anchors)), log(max(anchors)), length.out = 41L))
  grid_mult = sort(unique(c(grid_mult, anchors, 1)))
  chosen_idx = which(grid_mult == 1)
  anchor_index = match(anchors, grid_mult)
  s_grid = grid_mult * chosen_scale

  draws = lapply(anchor_fits, anchor_draws)
  reweights = lapply(seq_len(n_anchor), function(a) {
    anchor_reweight(draws[[a]], anchors[a] * chosen_scale, s_grid)
  })
  curve = assemble_curve(reweights, usable, ess_floor, anchor_index)

  # A masked point between two usable anchors means their usable radii do not
  # overlap; the curve has a gap there.
  for(k in seq_len(n_anchor - 1L)) {
    if(!usable[k] || !usable[k + 1L]) next
    between = grid_mult > anchors[k] & grid_mult < anchors[k + 1L]
    if(any(between & is.na(curve$anchor_used))) {
      warning(sprintf(
        paste0(
          "The usable radii of the %.2gx and %.2gx anchors do not overlap; ",
          "the curve has a gap between them. Add an anchor in that range."
        ),
        anchors[k], anchors[k + 1L]
      ), call. = FALSE)
    }
  }

  # --- Curve quantities: log BF, per-point MCSE, verdicts ---------------------
  # Cap the pooled PIP away from 0 and 1 (the same 1e-6 floor the RB machinery
  # uses) so an edge that saturates at some scale stays a large finite BF
  # rather than an infinity: the verdict is unchanged (still decisive), and the
  # curve, its range, and the plot stay finite.
  lbf_of = function(pip_mat) {
    t(apply(pip_mat, 1, function(p) {
      pc = pmin(pmax(p, 1e-6), 1 - 1e-6)
      log((pc / (1 - pc)) / prior_odds)
    }))
  }
  lbf_mat = lbf_of(curve$pip)
  chain_lbf = simplify2array(lapply(curve$chain_pip, lbf_of)) # P x E x C
  mcse_mat = apply(chain_lbf, c(1, 2), function(v) {
    v = v[is.finite(v)]
    if(length(v) >= 2L) stats::sd(v) / sqrt(length(v)) else NA_real_
  })
  unanimous_mat = apply(chain_lbf, c(1, 2), function(v) {
    vv = verdict_from_lbf(v[is.finite(v)], lthr)
    length(unique(vv)) <= 1L
  })
  verdict_mat = matrix(
    verdict_from_lbf(as.vector(lbf_mat), lthr),
    nrow = nrow(lbf_mat), ncol = n_edges
  )
  rownames(lbf_mat) = rownames(verdict_mat) = sprintf("m%.3g", grid_mult)

  # --- Run-to-run noise from the replicate anchor pair ------------------------
  # Same estimator as the curve (reweighting at the anchor's own scale is the
  # identity, so this is the plain indicator average through the prior odds).
  # The average runs over the INDICATOR draws, the unit the whole check reports
  # on: a bgmCompare fit's gamma draws are per gated parameter (one indicator
  # column repeated for every threshold difference it gates), which would both
  # misalign with prior_odds and overweight the main-effect differences in the
  # noise quantile.
  own_lbf = function(fit) {
    d = anchor_draws(fit)
    p = colMeans(do.call(rbind, d$indicator))
    log((p / (1 - p)) / prior_odds)
  }
  if(usable[rep_anchor] && rep_gate$usable) {
    wob = wobble_yardstick(own_lbf(anchor_fits[[rep_anchor]]), own_lbf(rep_fit))
    d_wobble = wob$per_edge
    wobble_q95 = wob$q95
    wobble_med = wob$median
    wobble_censored = wob$censored
  } else {
    # The failure was already reported once above, in plain voice.
    d_wobble = rep(NA_real_, n_edges)
    wobble_q95 = NA_real_
    wobble_med = NA_real_
    wobble_censored = 0L
  }

  # --- Per-edge chosen-scale sufficiency (from the original fit) --------------
  # A near-constant edge has NA MCSE; treat missing uncertainty as zero (it is
  # a saturated/zero-flip edge, masked below anyway).
  band = s0$band_half
  band[is.na(band)] = 0
  mcse0 = s0$mcse_lbf
  mcse0[is.na(mcse0)] = 0
  hw = pmax(band, 2 * mcse0)
  straddle = (!s0$zeroflip) &
    (((s0$lbf - hw < lthr) & (s0$lbf + hw > lthr)) |
      ((s0$lbf - hw < -lthr) & (s0$lbf + hw > -lthr)))
  disagree = (!s0$zeroflip) & !s0$unanimous
  insufficient = (!s0$zeroflip) & (disagree | straddle)
  insufficient[is.na(insufficient)] = FALSE
  # Record why each insufficient edge was flagged, for the report.
  insufficient_noisy = straddle & insufficient
  insufficient_noisy[is.na(insufficient_noisy)] = FALSE
  insufficient_disagree = disagree & insufficient
  insufficient_disagree[is.na(insufficient_disagree)] = FALSE

  # --- Mover category + stability interval on the dense curve -----------------
  mover = character(n_edges)
  stability_lower = rep(NA_real_, n_edges)
  stability_upper = rep(NA_real_, n_edges)
  for(e in seq_len(n_edges)) {
    v = verdict_mat[, e]
    moved = length(unique(v[!is.na(v)])) > 1L
    lbf_e = lbf_mat[, e]
    lbf_e = lbf_e[is.finite(lbf_e)]
    dlbf = if(length(lbf_e)) diff(range(lbf_e)) else 0
    mcse_e = mcse_mat[, e]
    mcse_e = mcse_e[is.finite(mcse_e)]
    move_thr = max(c(
      tolerance, if(length(mcse_e)) 2 * max(mcse_e) else 0,
      wobble_q95
    ), na.rm = TRUE)
    mover[e] = if(!moved) {
      "stable"
    } else if(dlbf > move_thr) {
      "moved-beyond-wobble"
    } else {
      "indistinguishable-from-wobble"
    }
    si = stability_interval(verdict_mat[, e], grid_mult, chosen_idx)
    stability_lower[e] = si[1]
    stability_upper[e] = si[2]
  }

  # --- Edges table ------------------------------------------------------------
  # Chosen-scale columns are the original fit's own reported analysis (RB
  # machinery), not a reweighted estimate: the 1x anchor is that fit.
  edges = data.frame(
    edge = edge_names,
    prior_inclusion_probability = prior_odds / (1 + prior_odds),
    chosen_scale_pip = s0$pip,
    chosen_scale_log_bf = s0$lbf,
    chosen_scale_mcse = s0$mcse_lbf,
    chosen_scale_verdict = s0$verdict,
    stability_lower = stability_lower,
    stability_upper = stability_upper,
    mover = mover,
    insufficient = insufficient,
    insufficient_noisy = insufficient_noisy,
    insufficient_disagree = insufficient_disagree,
    saturated = s0$zeroflip,
    stringsAsFactors = FALSE, row.names = NULL
  )
  # Per-anchor verdict columns are each fit's own exact RB verdict, not a
  # pooled-curve row.
  vcols = as.data.frame(t(anchor_verdict), stringsAsFactors = FALSE)
  names(vcols) = paste0("verdict_x", anchors)
  edges = cbind(edges, vcols)

  # --- Anchor table (per anchor fit + replicate, with convergence gates) ------
  all_gates = c(gates, list(rep_gate))
  secs = numeric(n_anchor)
  secs[non_unit] = walls[seq_along(non_unit)]

  # Captured sampler notes per grid row (original fit carries none of its own
  # here; its warnings surfaced when the user fit it).
  grid_notes = vector("list", n_anchor + 1L)
  grid_notes[[s0_idx]] = character(0)
  for(k in seq_along(non_unit)) grid_notes[[non_unit[k]]] = refit_notes[[k]]
  grid_notes[[n_anchor + 1L]] = refit_notes[[length(refit_notes)]]
  names(grid_notes) = c(sprintf("%gx", anchors), sprintf("%gx (replicate)", anchors[rep_anchor]))
  warmup_incomplete_of = function(f) {
    wi = tryCatch(f@nuts_diag$summary$warmup_incomplete, error = function(e) NA)
    if(is.null(wi)) NA else isTRUE(wi)
  }

  grid = data.frame(
    multiplier = c(anchors, anchors[rep_anchor]),
    scale = c(anchors, anchors[rep_anchor]) * chosen_scale,
    replicate = c(rep(FALSE, n_anchor), TRUE),
    original_fit = c(seq_len(n_anchor) == s0_idx, FALSE),
    usable = c(usable, rep_gate$usable),
    rhat_continuous = vapply(all_gates, `[[`, numeric(1), "rhat_cont"),
    rhat_continuous_max = vapply(all_gates, `[[`, numeric(1), "rhat_cont_max"),
    ess_continuous = vapply(all_gates, `[[`, numeric(1), "ess_cont"),
    inclusion_ess_min = vapply(all_gates, `[[`, numeric(1), "ess_incl"),
    rb_median_rhat = vapply(all_gates, `[[`, numeric(1), "rb_med_rhat"),
    warmup_incomplete = c(
      vapply(anchor_fits, warmup_incomplete_of, logical(1)),
      warmup_incomplete_of(rep_fit)
    ),
    sampler_notes = lengths(grid_notes) > 0,
    seconds = c(secs, walls[length(walls)]),
    row.names = NULL
  )

  structure(
    list(
      edges = edges,
      grid = grid,
      anchors = anchors,
      replicate_anchor = anchors[rep_anchor],
      multipliers = grid_mult,
      scales = s_grid,
      chosen_scale = chosen_scale,
      chosen_index = chosen_idx,
      anchor_index = anchor_index,
      log_bf = lbf_mat,
      log_bf_mcse = mcse_mat,
      verdict = verdict_mat,
      anchor_verdict = anchor_verdict,
      curve = list(
        ess = curve$ess, anchor_used = curve$anchor_used,
        ess_floor = ess_floor, unanimous = unanimous_mat
      ),
      wobble = list(
        q95 = wobble_q95, median = wobble_med, per_edge = d_wobble,
        censored = wobble_censored, anchor = anchors[rep_anchor]
      ),
      preferred_scale = preferred,
      unit = unit,
      vary = vary,
      evidence_threshold = evidence_threshold,
      tolerance = tolerance,
      refit_sampler = rs$method,
      warm = warm_metric,
      model_type = spec$model_type,
      runtime_seconds = sum(walls),
      edge_names = edge_names,
      refit_diagnostics = grid_notes,
      fits = if(keep_fits) refits else NULL
    ),
    class = "bgms_prior_sensitivity"
  )
}


# ------------------------------------------------------------------
# order_upper_tri_rowmajor
# ------------------------------------------------------------------
# Permutation that reorders the column-major upper-triangle extraction
# matrix[upper.tri(matrix)] into the row-major edge order used by the
# raw pairwise / indicator columns and edge names.
#
# @param p  Number of variables.
#
# Returns: integer permutation of length p(p-1)/2.
# ------------------------------------------------------------------
order_upper_tri_rowmajor = function(p) {
  # Column-major upper triangle visits (i, j) with j outer, i < j inner.
  cm = which(upper.tri(matrix(0, p, p)), arr.ind = TRUE)
  # Row-major target order sorts by row then column.
  order(cm[, "row"], cm[, "col"])
}


# ------------------------------------------------------------------
# wobble_yardstick
# ------------------------------------------------------------------
# Monte Carlo wobble from the chosen-scale refit and its replicate. The
# q95 pools threshold-relevant edges only (|log BF| <= 3 log(10) at s0,
# the former 3-log10 window in natural log units):
# near-saturated edges have an exploding log-odds derivative, and their
# replicate spread would inflate the yardstick past any genuine
# scale-driven move. The per-edge spread and its median cover all edges.
#
# @param lbf_s0   Per-edge natural log BF from the chosen-scale refit.
# @param lbf_rep  Per-edge natural log BF from the s0 replicate.
#
# Returns: list(q95, median, per_edge); q95 is NA when no edge is
# threshold-relevant.
# ------------------------------------------------------------------
wobble_yardstick = function(lbf_s0, lbf_rep) {
  per_edge = abs(lbf_s0 - lbf_rep)
  relevant = is.finite(lbf_s0) & abs(lbf_s0) <= 3 * log(10)
  # An edge that is threshold-relevant in one refit and saturated in the other
  # has a censored spread, not an infinite one. Pooling it would carry the
  # yardstick to Inf, which classes every verdict move as run-to-run noise.
  measurable = relevant & is.finite(per_edge)
  list(
    q95 = if(any(measurable)) {
      stats::quantile(per_edge[measurable], 0.95, names = FALSE)
    } else {
      NA_real_
    },
    median = stats::median(per_edge[is.finite(per_edge)]),
    per_edge = per_edge,
    censored = sum(relevant & !is.finite(per_edge))
  )
}


# ------------------------------------------------------------------
# stability_interval
# ------------------------------------------------------------------
# Maximal contiguous relative-scale range around the chosen scale over
# which the verdict is unchanged, using only reliable (non-NA) grid
# points. NA verdicts break the run.
#
# @param verdict     Per-grid-point verdict (may contain NA).
# @param relative    Relative scale at each grid point.
# @param chosen_idx  Grid index of the chosen scale.
#
# Returns: c(lower, upper) relative scale, or c(NA, NA).
# ------------------------------------------------------------------
stability_interval = function(verdict, relative, chosen_idx) {
  target = verdict[chosen_idx]
  if(is.na(target)) {
    return(c(NA_real_, NA_real_))
  }
  n = length(verdict)
  lo = chosen_idx
  while(lo - 1L >= 1L && !is.na(verdict[lo - 1L]) && verdict[lo - 1L] == target) {
    lo = lo - 1L
  }
  hi = chosen_idx
  while(hi + 1L <= n && !is.na(verdict[hi + 1L]) && verdict[hi + 1L] == target) {
    hi = hi + 1L
  }
  c(relative[lo], relative[hi])
}


# ==============================================================================
# Methods for bgms_prior_sensitivity
# ==============================================================================

#' @title Print a Prior Sensitivity Check
#'
#' @description
#' Prints the answer first: how many edge verdicts are robust to the slab
#' scale, which edges genuinely depend on it (by name), which cannot be
#' certified, the verdict counts per scale, and how the chosen scale compares
#' with the size of the estimated interactions. Machinery (sampler, refit
#' cost, the run-to-run noise band) is confined to a closing Details line.
#'
#' @param x A \code{bgms_prior_sensitivity} object.
#' @param max_rows Integer. Maximum edges to name in the scale-dependent
#'   table; the rest are counted and left to \code{$edges}. Default: \code{10}.
#' @param ... Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{prior_sensitivity_check}()}
#' @family diagnostics
#' @export
print.bgms_prior_sensitivity = function(x, max_rows = 10L, ...) {
  edges = x$edges
  n_edges = nrow(edges)
  mlab = sprintf("%.2gx", x$anchors)
  gr = range(x$multipliers)

  # Wrap prose to the console width (capped) so long interpolated content does
  # not overflow; tables and the labeled footer keep their own layout.
  w = min(getOption("width", 80L), 80L)
  wrap = function(...) writeLines(strwrap(paste0(...), width = w))

  # bgm() traces edge indicators across the interaction slab; bgmCompare()
  # traces difference indicators across the difference slab.
  unit = x$unit %||% list(
    noun = "edge", nouns = "edges",
    headline = "are the edge verdicts robust to the slab scale?"
  )

  cat(sprintf("Prior sensitivity check: %s\n", unit$headline))
  wrap(sprintf(
    "Bayes-factor curve from %.2gx to %.2gx the chosen scale (anchors at %s; the 1x anchor is the original fit); %d %s.",
    gr[1], gr[2], paste(mlab, collapse = ", "), n_edges, unit$nouns
  ))
  # What moved along the curve. On a model with a prior on the precision
  # diagonal the slab scale and that diagonal are tied through the standardized
  # frame, so the sweep holds one of the two fixed and the curve means
  # different things depending on which.
  vary_mode = if(is.null(x$vary)) "none" else x$vary$mode
  if(identical(vary_mode, "slab")) {
    wrap(
      "Varied: the interaction slab scale alone; the prior on the precision ",
      "diagonal is the fitted one at every scale (vary = \"slab\")."
    )
  } else if(identical(vary_mode, "slab-and-diagonal")) {
    wrap(sprintf(
      paste(
        "Varied: the interaction slab scale together with the precision diagonal,",
        "holding the standardized rate eta at %.3g so the prior's shape is fixed",
        "and its overall scale moves (vary = \"slab-and-diagonal\")."
      ),
      x$vary$eta
    ))
  }
  cat("\n")

  # An indicator the sampler never updated has no verdict to be robust or
  # sensitive; bgmCompare() leaves the main-effect difference indicators there
  # unless main_difference_selection = TRUE.
  no_verdict = is.na(edges$chosen_scale_verdict)
  if(any(no_verdict)) {
    wrap(sprintf(
      paste(
        "%d of the %d %s were never updated by the sampler and carry no verdict",
        "at any scale; they are excluded from the counts below."
      ),
      sum(no_verdict), n_edges, unit$nouns
    ))
    cat("\n")
    edges = edges[!no_verdict, , drop = FALSE]
    x$anchor_verdict = x$anchor_verdict[, !no_verdict, drop = FALSE]
    n_edges = nrow(edges)
  }

  # Stability headline: the scale range over which the verdicts hold.
  full_span = !is.na(edges$stability_lower) & !is.na(edges$stability_upper) &
    edges$stability_lower <= gr[1] * 1.001 & edges$stability_upper >= gr[2] * 0.999
  if(all(full_span)) {
    wrap(sprintf(
      "All %d verdicts hold from %.2gx to %.2gx the chosen scale.",
      n_edges, gr[1], gr[2]
    ))
  } else {
    wrap(sprintf(
      "%d of %d verdicts hold across the whole %.2gx-%.2gx range; the exceptions are named below.",
      sum(full_span), n_edges, gr[1], gr[2]
    ))
  }
  cat("\n")

  # One category per edge: not-certifiable wins, then the mover category.
  uncert = edges$insufficient %in% TRUE
  robust = !uncert & edges$mover == "stable"
  within = !uncert & edges$mover == "indistinguishable-from-wobble"
  beyond = !uncert & edges$mover == "moved-beyond-wobble"
  counts = c(sum(robust), sum(within), sum(beyond), sum(uncert))
  labels = c(
    "robust (same verdict at every scale)",
    "changed, within run-to-run noise",
    "changed, beyond run-to-run noise",
    "not certifiable (too noisy to assess)"
  )
  keep = counts > 0 | seq_along(counts) <= 3L
  cat(sprintf("  %-38s %4d\n", labels[keep], counts[keep]), sep = "")
  cat("\n")

  # Name the units whose verdict genuinely depends on the scale.
  if(sum(beyond) > 0) {
    cat(if(sum(beyond) == 1L) {
      sprintf("1 %s's verdict genuinely depends on the scale:\n", unit$noun)
    } else {
      sprintf(
        "%d %s' verdicts genuinely depend on the scale:\n",
        sum(beyond), unit$nouns
      )
    })
    idx = which(beyond)
    show = utils::head(idx, max_rows)
    vmat = t(x$anchor_verdict[, show, drop = FALSE])
    tab = cbind(edge = edges$edge[show], vmat)
    colnames(tab) = c(unit$noun, mlab)
    rownames(tab) = rep("", nrow(tab))
    print(tab, quote = FALSE, print.gap = 2)
    if(length(idx) > length(show)) {
      cat(sprintf("  ...and %d more; see $edges.\n", length(idx) - length(show)))
    }
    cat("\n")
  }

  # Name the edges too noisy to assess, say why, and say what to do.
  if(sum(uncert) > 0) {
    nm = edges$edge[uncert]
    shown = utils::head(nm, max_rows)
    tail_txt = if(length(nm) > length(shown)) {
      sprintf(", and %d more (see $edges)", length(nm) - length(shown))
    } else {
      ""
    }
    one = length(nm) == 1L
    wrap(sprintf(
      "%d %s %s too noisy to assess: %s%s.",
      length(nm), if(one) unit$noun else unit$nouns,
      if(one) "is" else "are",
      paste(shown, collapse = ", "), tail_txt
    ))
    # Name only the subcause(s) that actually occur among these edges.
    noisy = any(edges$insufficient_noisy[uncert])
    disagree = any(edges$insufficient_disagree[uncert])
    subj = if(one) "Its" else "Their"
    obj = if(one) "it" else "them"
    these = if(one) "this verdict" else "these verdicts"
    cause = if(noisy && disagree) {
      sprintf(
        "%s Bayes factor sits within Monte Carlo error of an evidence threshold, or %s chains disagree on the verdict,",
        subj, tolower(subj)
      )
    } else if(disagree) {
      sprintf("%s chains disagree on the verdict", subj)
    } else {
      sprintf("%s Bayes factor sits within Monte Carlo error of an evidence threshold", subj)
    }
    wrap(sprintf(
      paste(
        "%s at the chosen scale itself; a rerun with a fresh seed could flip",
        "%s without any prior change. Run more iterations to settle %s before",
        "reading %s sensitivity."
      ),
      cause, obj, these, tolower(subj)
    ))
    cat("\n")
  }

  # Verdict counts at the anchor scales.
  cat("Verdict counts by scale (at the anchors):\n")
  vt = apply(x$anchor_verdict, 1, function(v) {
    c(
      presence = sum(v == "presence", na.rm = TRUE),
      undecided = sum(v == "undecided", na.rm = TRUE),
      absence = sum(v == "absence", na.rm = TRUE)
    )
  })
  colnames(vt) = mlab
  out = utils::capture.output(print(vt))
  cat(paste0("  ", out, collapse = "\n"), "\n", sep = "")
  ab = vt["absence", ]
  if(ab[length(ab)] > ab[1]) {
    wrap(sprintf(
      paste(
        "More absence at wider scales is expected: a wider slab strengthens",
        "evidence against borderline %s."
      ),
      unit$nouns
    ))
  }
  cat("\n")

  # Refits that failed their convergence check.
  bad = x$grid[!x$grid$usable & !x$grid$replicate, ]
  if(nrow(bad) > 0) {
    bl = sprintf("%.2gx", bad$multiplier)
    bl = if(length(bl) > 1L) {
      paste0(paste(bl[-length(bl)], collapse = ", "), " and ", bl[length(bl)])
    } else {
      bl
    }
    wrap(sprintf(
      "The %s refit%s did not converge and %s excluded from the verdicts.",
      bl, if(nrow(bad) == 1L) "" else "s", if(nrow(bad) == 1L) "is" else "are"
    ))
    cat("\n")
  }

  # How the chosen scale compares with the estimated interactions.
  ps = x$preferred_scale
  if(is.finite(ps$s_hat)) {
    ratio = x$chosen_scale / ps$s_hat
    if(ratio > 2) {
      wrap(sprintf(
        "Note: the chosen scale (%.3g) is much wider than the estimated interactions (about %.3g [%.3g, %.3g]); absence verdicts in particular depend on this choice.",
        x$chosen_scale, ps$s_hat, ps$lo, ps$hi
      ))
    } else if(ratio < 0.5) {
      wrap(sprintf(
        "Note: the chosen scale (%.3g) is much narrower than the estimated interactions (about %.3g [%.3g, %.3g]); presence verdicts in particular depend on this choice.",
        x$chosen_scale, ps$s_hat, ps$lo, ps$hi
      ))
    } else {
      wrap(sprintf(
        "The chosen scale (%.3g) matches the size of the estimated interactions (about %.3g).",
        x$chosen_scale, ps$s_hat
      ))
    }
    cat("\n")
  }

  secs = if(x$runtime_seconds < 10) {
    sprintf("%.1f s", x$runtime_seconds)
  } else {
    sprintf("%.0f s", x$runtime_seconds)
  }
  gr = range(x$anchors)
  n_refit = nrow(x$grid) - 1L
  start_txt = if(x$warm) "warm-started from the original fit" else "cold starts with full warmup"
  cat(sprintf(
    "Method:  %d-point curve from %d anchor fits (%.2gx to %.2gx the chosen scale),\n         joined by importance reweighting; the 1x anchor is the original fit.\n         Points with reweighting effective sample size below %g are not shown.\n",
    length(x$multipliers), length(x$anchors), gr[1], gr[2], x$curve$ess_floor
  ))
  cat(sprintf(
    "Refits:  %d %s refits, %s, %s total.\n",
    n_refit, x$refit_sampler, start_txt, secs
  ))
  if(is.na(x$wobble$q95)) {
    cat(sprintf(
      "Noise:   no threshold-relevant %s had a measurable spread between two identical\n         refits at %.2gx, so there is no run-to-run yardstick; verdict moves are\n         judged against the tolerance and Monte Carlo error alone.\n",
      unit$noun, x$wobble$anchor
    ))
  } else {
    cat(sprintf(
      "Noise:   two identical refits at %.2gx differed by up to %.2g log BF across\n         threshold-relevant %s; verdict moves smaller than that are reported\n         as run-to-run noise, not prior sensitivity.\n",
      x$wobble$anchor, x$wobble$q95, unit$nouns
    ))
    if(isTRUE(x$wobble$censored > 0)) {
      cat(sprintf(
        "         %d edge%s saturated in one of the two and %s left out of that spread.\n",
        x$wobble$censored, if(x$wobble$censored == 1L) "" else "s",
        if(x$wobble$censored == 1L) "was" else "were"
      ))
    }
  }
  cat("See ?prior_sensitivity_check for the full construction.\n")
  invisible(x)
}


# ------------------------------------------------------------------
# mover_palette
# ------------------------------------------------------------------
# Fixed-order categorical colors for the named mover lines (Okabe-Ito,
# colorblind-safe; warm/cool alternating for adjacent separation). Every
# mover line is also name-labeled, so identity never rides on color alone.
# ------------------------------------------------------------------
mover_palette = function() {
  c(
    "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
    "#56B4E9", "#882255", "#44AA99", "#661100", "#999933"
  )
}


# ------------------------------------------------------------------
# spread_labels
# ------------------------------------------------------------------
# Nudge label y-positions apart until adjacent labels are at least `gap`
# apart, preserving order. Returns the adjusted positions.
# ------------------------------------------------------------------
spread_labels = function(y, gap) {
  ord = order(y)
  ys = y[ord]
  for(k in seq_along(ys)[-1]) {
    if(ys[k] - ys[k - 1] < gap) ys[k] = ys[k - 1] + gap
  }
  out = y
  out[ord] = ys
  out
}


#' @title Plot a Prior Sensitivity Check
#'
#' @description
#' One panel. Each edge's natural log inclusion-Bayes-factor curve
#' (the evidence for the edge) is drawn across the anchored scale range, with
#' dots at the anchor scales; curve points masked for low importance ESS
#' leave visible gaps, and an edge that saturates at some scale is capped at a
#' large finite Bayes factor rather than running off to infinity. The
#' \code{max_labels} widest-swinging edges whose verdict genuinely depends on
#' the scale are colored and labeled by name; every other edge is the muted
#' background. The figure carries no counts of its own: the muted background
#' holds the robust edges and any scale-dependent edge past \code{max_labels}
#' alike, and \code{print()} is where the counts are and where the
#' scale-dependent edges are named. The shaded band is the undecided zone
#' between the evidence thresholds; the zones are labeled at the left edge.
#'
#' @details
#' The drawn curves clamp the pooled inclusion probability at
#' \eqn{1 - 10^{-6}}, which caps a plotted log Bayes factor at
#' \eqn{\ln(10^6) \approx 13.8}. A curve running flat along 13.8 has reached
#' that display cap; it is not evidence levelling off. The uncapped value at
#' the chosen scale is \code{x$edges$chosen_scale_log_bf}, which
#' \code{print()} reports and which can be far larger. \code{x$edges$saturated} does not mark capped curves: it
#' records that the edge's inclusion indicator never flipped in the chain,
#' which is a different condition.
#'
#' @param x A \code{bgms_prior_sensitivity} object.
#' @param max_labels Integer. Maximum scale-dependent edges to color and
#'   label by name, taken in order of evidence swing; the rest stay in the
#'   muted background and are named by \code{print()}. Default: \code{10}.
#' @param ... Ignored.
#'
#' @return \code{x}, invisibly. Called for the side effect of drawing.
#'
#' @seealso \code{\link{prior_sensitivity_check}()}
#' @family diagnostics
#' @export
plot.bgms_prior_sensitivity = function(x, max_labels = 10L, ...) {
  rel = x$multipliers
  thr = log(x$evidence_threshold)
  edges = x$edges
  n_edges = nrow(edges)

  # Same partition as print(): an edge the refits cannot certify is not
  # reported as scale-dependent.
  movers = which(
    edges$mover == "moved-beyond-wobble" &
      !(edges$saturated %in% TRUE) & !(edges$insufficient %in% TRUE)
  )
  # Label the largest evidence swings first.
  swing = apply(x$log_bf, 2, function(v) {
    v = v[is.finite(v)]
    if(length(v)) diff(range(v)) else 0
  })
  movers = movers[order(swing[movers], decreasing = TRUE)]
  named = utils::head(movers, max_labels)

  style = bgms_style()
  faint = grDevices::adjustcolor(style$muted, 0.35)
  name_cex = style$cex_caption
  # The leader from the end of a curve to its name, in inches, so the same gap
  # opens on any device.
  leader = 0.16

  # The right margin holds the edge names, so it is measured for the names
  # this check produced rather than set to a width that happened to fit an
  # earlier example. Anything shorter clipped them at the device edge.
  right = leader / graphics::par("csi") +
    margin_lines_for(edges$edge[named], cex = name_cex, pad = 0.8)
  style = bgms_panel_par(mar = c(5.6, 5.0, 3.4, right), scale = 1)
  on.exit(graphics::par(style$old_par), add = TRUE)

  # The y-range covers the decision band and every named line; edges that
  # live entirely outside it never change verdict and clip silently.
  named_bf = x$log_bf[, named, drop = FALSE]
  named_bf = named_bf[is.finite(named_bf)]
  ylim = range(c(-1.5 * thr, 1.5 * thr, named_bf))
  ylim = ylim + c(-0.4, 0.4)
  y_axis = bgms_axis_range(ylim)
  # The scale axis is logarithmic and its ticks are the anchors, so the data
  # range is opened on the log scale to hold the axis off the corner.
  log_rel = log10(range(rel))
  xlim = 10^(log_rel + c(-1, 1) * style$eps * diff(log_rel))

  unit = x$unit %||% list(noun = "edge", nouns = "edges")
  scale_label = if(identical(unit$noun, "difference")) {
    "Difference scale, relative to your fit"
  } else {
    "Slab scale, relative to your fit"
  }

  graphics::plot(NA, NA,
    xlim = xlim, ylim = y_axis$lim, log = "x", axes = FALSE,
    xlab = "", ylab = "", main = ""
  )
  bgms_axis(1, x$anchors, labels = sprintf("%.2gx", x$anchors), style = style)
  bgms_axis(2, y_axis$at,
    labels = formatC(y_axis$at, format = "g"), style = style
  )
  graphics::mtext(scale_label,
    side = 1, line = 2.9, cex = style$cex_lab, col = style$ink
  )
  graphics::mtext(
    sprintf("Evidence for the %s (log BF)", unit$noun),
    side = 2, line = 3.3, las = 0, cex = style$cex_lab, col = style$ink
  )
  # Verdict zones: shaded undecided band, dashed thresholds, margin labels.
  usr = graphics::par("usr")
  pin = graphics::par("pin")
  # Inches per unit on each axis, so every offset below is a device length
  # rather than a fraction of a range that changes with the data.
  per_log_unit = pin[1] / diff(usr[1:2])
  per_y_unit = pin[2] / diff(usr[3:4])

  graphics::rect(10^usr[1], -thr, 10^usr[2], thr,
    col = grDevices::adjustcolor("grey60", 0.12), border = NA
  )
  graphics::abline(h = c(-thr, thr), col = style$muted, lty = 2,
    lwd = style$lwd_axis
  )
  graphics::abline(v = 1, col = style$muted, lty = 3, lwd = style$lwd_axis)
  # "chosen" named a decision the reader had to reconstruct; the line is the
  # scale their own fit was run at, and that is what it now says.
  graphics::mtext("your fit",
    side = 3, at = 1, line = 0.3, cex = style$cex_annotation, col = style$muted
  )
  zone_x = 10^(usr[1] + 0.06 / per_log_unit)
  pad = 0.05 / per_y_unit
  graphics::text(zone_x, thr + pad, "presence",
    adj = c(0, 0), cex = style$cex_caption, col = style$muted
  )
  graphics::text(zone_x, thr - pad, "undecided",
    adj = c(0, 1), cex = style$cex_caption, col = style$muted
  )
  graphics::text(zone_x, -thr - pad, "absence",
    adj = c(0, 1), cex = style$cex_caption, col = style$muted
  )

  # Background: every other edge in one muted color.
  for(e in setdiff(seq_len(n_edges), named)) {
    if(isTRUE(edges$saturated[e])) next
    trajectory(rel, x$log_bf[, e], col = faint, lwd = 1)
  }

  # Foreground: the named movers, colored in fixed order and name-labeled.
  if(length(named)) {
    pal = mover_palette()
    # Anchor each name at the line's last finite point (0 if none is finite).
    anchor = vapply(named, function(e) {
      v = x$log_bf[, e]
      v = v[is.finite(v)]
      if(length(v)) v[length(v)] else 0
    }, numeric(1))
    end_y = spread_labels(anchor, gap = 0.05 * diff(y_axis$lim))
    stub_x = 10^(usr[2] + 0.35 * leader / per_log_unit)
    text_x = 10^(usr[2] + leader / per_log_unit)
    for(k in seq_along(named)) {
      e = named[k]
      trajectory(rel, x$log_bf[, e], col = pal[k], lwd = style$lwd_curve)
      graphics::points(rel[x$anchor_index], x$log_bf[x$anchor_index, e],
        col = pal[k], pch = 16, cex = 0.9
      )
      graphics::segments(
        max(rel), anchor[k], stub_x, end_y[k],
        col = grDevices::adjustcolor(pal[k], 0.5), lwd = 0.8, xpd = NA
      )
      graphics::text(text_x, end_y[k], edges$edge[e],
        xpd = NA, adj = 0, cex = name_cex, col = style$ink
      )
    }
  }

  invisible(x)
}


# ------------------------------------------------------------------
# trajectory
# ------------------------------------------------------------------
# One edge's evidence, as one line.
#
# The reweighted curve is masked wherever the importance ESS falls under the
# floor, so the raw series has holes in it, and drawing it directly leaves
# stubs of curve with the anchor estimates floating unattached between them --
# which reads as a broken plot rather than as a masked one. The finite points,
# anchors included, are joined in order instead: one continuous trajectory per
# edge, straight between the certified stretches it passes through.
#
# @param at      The scale grid.
# @param values  The log Bayes factor on that grid, with holes.
# @param ...     Passed to graphics::lines().
# ------------------------------------------------------------------
trajectory = function(at, values, ...) {
  keep = is.finite(values)
  if(sum(keep) < 2L) {
    return(invisible(NULL))
  }
  graphics::lines(at[keep], values[keep], ...)
  invisible(NULL)
}

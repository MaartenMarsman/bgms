# build_output_bgm() — unified GGM + OMRF output builder
#
# Split out of build_output.R (cleanup S4). Shared helpers and the build_output()
# dispatcher remain in build_output.R.


# build_output_bgm()  --- unified GGM + OMRF
# ==============================================================================
#
# The two paths share ~80% of logic. Differences:
#   1. Parameter naming: GGM uses "Var (precision)", OMRF uses "Var (k)"
#   2. Main posterior mean shape: GGM = px1, OMRF = pxmax_categories
# ==============================================================================
build_output_bgm = function(spec, raw) {
  d = spec$data
  v = spec$variables
  p = spec$prior
  s = spec$sampler

  is_continuous = v$is_continuous
  num_variables = d$num_variables
  data_columnnames = d$data_columnnames
  edge_selection = p$edge_selection
  edge_prior = p$edge_prior

  # Keep the raw chains for the Z-ratio trust gauge: it needs the untouched
  # indicator layout and the per-chain zratio block, both dropped by the
  # normalization below.
  zratio_chains = if(identical(p$precision_graph_prior, "hierarchical")) {
    raw
  } else {
    NULL
  }

  # --- Normalize raw C++ output -----------------------------------------------
  # The C++ GGM/OMRF backends return a flat `samples` matrix (params x iters)
  # via convert_results_to_list(). Split into main and pairwise components and
  # transpose to (iters x params) --- the layout that MCMC summary functions
  # expect.
  if(is_continuous) {
    # GGM: samples contain the upper triangle of the precision matrix
    # (row-major). Diagonal entries are "main"; off-diagonal are "pairwise".
    diag_idx = integer(num_variables)
    offdiag_idx = integer(num_variables * (num_variables - 1L) / 2L)
    pos = 0L
    di = 0L
    oi = 0L
    for(i in seq_len(num_variables)) {
      for(j in i:num_variables) {
        pos = pos + 1L
        if(i == j) {
          di = di + 1L
          diag_idx[di] = pos
        } else {
          oi = oi + 1L
          offdiag_idx[oi] = pos
        }
      }
    }
    raw = lapply(raw, function(chain) {
      samples_t = t(chain$samples)
      res = list(
        main_samples     = samples_t[, diag_idx, drop = FALSE],
        pairwise_samples = samples_t[, offdiag_idx, drop = FALSE],
        userInterrupt    = isTRUE(chain$userInterrupt),
        chain_id         = chain$chain_id
      )
      if(!is.null(chain$indicator_samples)) {
        res$indicator_samples = t(chain$indicator_samples)[, offdiag_idx, drop = FALSE]
      }
      if(!is.null(chain$rb_inclusion_samples)) {
        res$rb_inclusion_samples = t(chain$rb_inclusion_samples)[, offdiag_idx, drop = FALSE]
      }
      if(!is.null(chain$rb_counts)) {
        res$rb_counts = chain$rb_counts[offdiag_idx, , drop = FALSE]
      }
      if(!is.null(chain$allocation_samples)) {
        res$allocations = t(chain$allocation_samples)
      }
      if(!is.null(chain$inclusion_parameter_samples)) {
        res$inclusion_parameter = as.numeric(chain$inclusion_parameter_samples)
      }
      # Final adaptation-averaged NUTS step size (NaN otherwise), for warm starts.
      res$step_size = if(is.null(chain$step_size)) NA_real_ else as.numeric(chain$step_size)
      # Final adapted diagonal metric (NULL otherwise), for warm-starting refits.
      res$inv_mass = if(is.null(chain$inv_mass)) NULL else as.numeric(chain$inv_mass)
      attach_diagnostic_traces(res, chain)
    })
  } else {
    # OMRF: the first num_thresholds params are main effects, the rest are
    # pairwise. NUTS diagnostics use bare names from C++; rename to the
    # trailing-__ convention expected by summarize_nuts_diagnostics().
    num_thresholds = spec$precomputed$num_thresholds
    raw = lapply(raw, function(chain) {
      samples_t = t(chain$samples)
      n_params = ncol(samples_t)
      res = list(
        main_samples     = samples_t[, seq_len(num_thresholds), drop = FALSE],
        pairwise_samples = samples_t[, seq(num_thresholds + 1L, n_params), drop = FALSE],
        userInterrupt    = isTRUE(chain$userInterrupt),
        chain_id         = chain$chain_id
      )
      if(!is.null(chain$indicator_samples)) {
        res$indicator_samples = t(chain$indicator_samples)
      }
      if(!is.null(chain$rb_inclusion_samples)) {
        res$rb_inclusion_samples = t(chain$rb_inclusion_samples)
      }
      if(!is.null(chain$rb_counts)) {
        res$rb_counts = chain$rb_counts
      }
      if(!is.null(chain$allocation_samples)) {
        res$allocations = t(chain$allocation_samples)
      }
      if(!is.null(chain$inclusion_parameter_samples)) {
        res$inclusion_parameter = as.numeric(chain$inclusion_parameter_samples)
      }
      # Final adaptation-averaged NUTS step size (NaN otherwise), for warm starts.
      res$step_size = if(is.null(chain$step_size)) NA_real_ else as.numeric(chain$step_size)
      # Final adapted diagonal metric (NULL otherwise), for warm-starting refits.
      res$inv_mass = if(is.null(chain$inv_mass)) NULL else as.numeric(chain$inv_mass)
      attach_diagnostic_traces(res, chain)
    })
  }

  # --- Parameter names --------------------------------------------------------
  if(is_continuous) {
    # GGM: raw samples are on precision scale; label accordingly
    names_main = paste0(data_columnnames, " (precision)")
    is_ordinal_variable = NULL
    num_categories = NULL
  } else {
    # OMRF: per-category thresholds or BC linear/quadratic
    is_ordinal_variable = v$is_ordinal
    num_categories = d$num_categories
    names_main = character()
    for(vi in seq_len(num_variables)) {
      if(is_ordinal_variable[vi]) {
        cats = ordinal_threshold_labels(
          num_categories[vi],
          if(!is.null(d$category_levels)) d$category_levels[[vi]] else NULL
        )
        names_main = c(
          names_main,
          paste0(data_columnnames[vi], " (", cats, ")")
        )
      } else {
        names_main = c(
          names_main,
          paste0(data_columnnames[vi], " (linear)"),
          paste0(data_columnnames[vi], " (quadratic)")
        )
      }
    }
  }

  edge_names = character()
  for(i in seq_len(num_variables - 1)) {
    for(j in seq(i + 1, num_variables)) {
      edge_names = c(
        edge_names,
        paste0(data_columnnames[i], "-", data_columnnames[j])
      )
    }
  }

  # --- Lazy MCMC diagnostics cache --------------------------------------------
  # Store normalized raw chains and metadata so that ensure_summaries() can
  # compute ESS/Rhat/MCSE on first demand. Posterior_summary_* fields start
  # as NULL and are populated lazily.
  cache = new.env(parent = emptyenv())
  cache$raw = raw
  cache$edge_selection = edge_selection
  cache$names_main = names_main
  cache$edge_names = edge_names
  cache$is_continuous = is_continuous
  cache$model_type = if(is_continuous) "ggm" else "omrf"
  cache$summaries_computed = FALSE

  # --- Compute posterior means from raw samples (cheap) -----------------------
  pooled_main = do.call(rbind, lapply(raw, function(ch) ch$main_samples))
  pooled_pair = do.call(rbind, lapply(raw, function(ch) ch$pairwise_samples))
  main_means = colMeans(pooled_main)
  pair_means = colMeans(pooled_pair)

  results = list()

  # --- Edge selection summaries (deferred) ------------------------------------
  has_sbm = FALSE
  if(edge_selection) {
    has_sbm = identical(edge_prior, "Stochastic-Block") &&
      "allocations" %in% names(raw[[1]])

    if(has_sbm) {
      sbm_convergence = summarize_alloc_pairs(
        allocations = lapply(raw, `[[`, "allocations"),
        node_names  = data_columnnames
      )
      results$posterior_summary_pairwise_allocations = sbm_convergence$sbm_summary
      co_occur_matrix = sbm_convergence$co_occur_matrix
    }

    if("inclusion_parameter" %in% names(raw[[1]])) {
      results$inclusion_parameter_samples =
        lapply(raw, `[[`, "inclusion_parameter")
    }
  }

  # Per-chain final NUTS step size, retained on the fit so a refit can
  # warm-start it (NA for non-NUTS runs). Not user-facing.
  results$refit_step_sizes = vapply(
    raw, function(ch) if(is.null(ch$step_size)) NA_real_ else ch$step_size,
    numeric(1)
  )

  # Per-chain final NUTS diagonal metric (inverse mass), retained so a refit can
  # warm-start it (NULL for non-NUTS runs). Not user-facing.
  if(!is.null(raw[[1]]$inv_mass)) {
    results$refit_inv_mass = lapply(raw, `[[`, "inv_mass")
  }

  # --- Posterior mean: main ---------------------------------------------------
  if(is_continuous) {
    # GGM has no main effects
    results$posterior_mean_main = NULL
  } else {
    # OMRF: p x max_categories matrix
    num_params = ifelse(is_ordinal_variable, num_categories, 2L)
    max_num_categories = max(num_params)

    pmm = matrix(NA, nrow = num_variables, ncol = max_num_categories)
    start = 0L
    stop = 0L
    for(vi in seq_len(num_variables)) {
      if(is_ordinal_variable[vi]) {
        start = stop + 1L
        stop = start + num_categories[vi] - 1L
        pmm[vi, seq_len(num_categories[vi])] = main_means[start:stop]
      } else {
        start = stop + 1L
        stop = start + 1L
        pmm[vi, 1:2] = main_means[start:stop]
      }
    }
    results$posterior_mean_main = pmm
    rownames(results$posterior_mean_main) = data_columnnames
    colnames(results$posterior_mean_main) = paste0("cat (", seq_len(ncol(pmm)), ")")
  }

  # --- Posterior mean: associations -------------------------------------------
  # For GGM: C++ stores precision; convert to association scale (* -0.5).
  # For OMRF: C++ already stores association-scale values.
  associations = matrix(0,
    nrow = num_variables, ncol = num_variables,
    dimnames = list(data_columnnames, data_columnnames)
  )
  associations[lower.tri(associations)] = pair_means
  associations = associations + t(associations)
  if(is_continuous) {
    associations = -0.5 * associations
  }
  results$posterior_mean_pairwise = associations

  # --- Residual variance (GGM only) -------------------------------------------
  # C++ stores precision diagonal; residual variance = 1 / precision.
  # Average per-sample inversions to avoid Jensen's inequality bias.
  if(is_continuous) {
    results$posterior_mean_residual_variance = colMeans(1 / pooled_main)
    names(results$posterior_mean_residual_variance) = data_columnnames
  }

  # --- Posterior mean: indicator + SBM ----------------------------------------
  if(edge_selection) {
    # Report the Rao-Blackwellized inclusion probability as the canonical
    # estimate, matching posterior_summary_indicator and the default of
    # extract_posterior_inclusion_probabilities(); fall back to the raw
    # indicator average for fits without RB draws.
    if(!is.null(raw[[1]][["rb_inclusion_samples"]])) {
      pooled_ind = do.call(rbind, lapply(raw, function(ch) ch$rb_inclusion_samples))
    } else {
      pooled_ind = do.call(rbind, lapply(raw, function(ch) ch$indicator_samples))
    }
    indicator_means = colMeans(pooled_ind, na.rm = TRUE)
    results$posterior_mean_indicator = matrix(0,
      nrow = num_variables, ncol = num_variables,
      dimnames = list(data_columnnames, data_columnnames)
    )
    results$posterior_mean_indicator[lower.tri(results$posterior_mean_indicator)] =
      indicator_means
    results$posterior_mean_indicator = results$posterior_mean_indicator +
      t(results$posterior_mean_indicator)

    if(has_sbm) {
      results$posterior_mean_coclustering_matrix = co_occur_matrix
      results = attach_sbm_posterior_summary(results, raw, build_arguments(spec))
    }
  }

  # --- arguments + class ------------------------------------------------------
  results$arguments = build_arguments(spec)
  # Report the number of chains actually kept; failed chains are dropped
  # upstream, so raw holds only the survivors.
  results$arguments$num_chains = length(raw)
  class(results) = "bgms"

  # --- raw_samples ------------------------------------------------------------
  # Allocation samples are per node (one column per variable), so they carry
  # node names, not edge names.
  alloc_names = if(identical(edge_prior, "Stochastic-Block")) {
    data_columnnames
  } else {
    NULL
  }
  results$raw_samples = build_raw_samples_list(
    raw, edge_selection, edge_prior, names_main, edge_names, alloc_names
  )

  # --- Attach lazy cache for deferred diagnostics -----------------------------
  # NULL placeholders ensure names(fit) lists these fields for easybgm compat.
  # The $.bgms method routes actual access through the cache.
  # Use list(NULL) because results$x = NULL removes the element in R.
  results["posterior_summary_main"] = list(NULL)
  results["posterior_summary_pairwise"] = list(NULL)
  if(edge_selection) {
    results["posterior_summary_indicator"] = list(NULL)
  }
  results$cache = cache

  # --- Sampler diagnostics ----------------------------------------------------
  results = attach_sampler_diagnostics(
    results, raw, s$update_method, s$nuts_max_depth,
    names_main = names_main, names_pairwise = edge_names,
    target_accept = s$target_accept
  )

  # --- Z-ratio trust gauge (hierarchical graph-prior spec) ----------------------
  # Only when the in-chain gauge actually ran; it is on by default but can be
  # switched off (options(bgms.zratio_gauge_sweeps = 0)), so guard against
  # empty output.
  if(!is.null(zratio_chains) && zratio_gauge_present(zratio_chains)) {
    # Harm channel inputs: per-chain edge-inclusion probabilities from the
    # normalized chains (indicator_samples is iters x off-diagonal edges) and
    # the edge-prior identity from the spec.
    harm_inputs = if(edge_selection) {
      pip = lapply(raw, function(ch) {
        if(is.null(ch$indicator_samples)) NULL else colMeans(ch$indicator_samples)
      })
      zratio_harm_inputs(
        pip, p$edge_prior,
        a = p$beta_bernoulli_alpha, b = p$beta_bernoulli_beta
      )
    } else {
      NULL
    }
    results$zratio_diag = summarize_zratio_gauge(
      zratio_chains,
      verbose = TRUE, harm_inputs = harm_inputs
    )
  }

  # Extrapolation notice: independent of the gauge (which can be switched off), so
  # a fit that deployed the surface beyond its validated hull is never silent.
  if(!is.null(zratio_chains) && isTRUE(getOption("bgms.verbose", TRUE))) {
    zratio_extrapolation_notice(zratio_chains)
  }

  # Single vignette pointer covering both diagnostic blocks: print once if
  # either the NUTS diagnostics or the trust gauge reported issues.
  if(isTRUE(getOption("bgms.verbose", TRUE)) &&
    (isTRUE(results$nuts_diag$has_issues) ||
      isTRUE(results$zratio_diag$flagged))) {
    cat("See vignette('diagnostics') for guidance.\n")
  }

  results$.bgm_spec = spec
  if(needs_easybgm_s3_compat()) {
    results
  } else {
    s3_list_to_bgms(results)
  }
}

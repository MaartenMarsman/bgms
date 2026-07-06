# build_output_mixed_mrf() — Mixed-MRF output builder
#
# Split out of build_output.R (cleanup S4).


# ==============================================================================
# build_output_mixed_mrf()  --- Mixed MRF builder
# ==============================================================================
#
# Handles the mixed discrete + continuous parameter layout:
#   C++ flat vector: [main_discrete | pairwise_discrete_ut | main_continuous | pairwise_cross | pairwise_continuous_ut]
#   C++ indicators:  [Gxx_ut | Gyy_ut | Gxy]
#
# Splits into main (discrete thresholds, continuous means),
# quadratic (continuous diagonal), and pairwise (discrete,
# continuous off-diagonal, cross). The continuous diagonal is
# stored separately in posterior_mean_residual_variance;
# the diagonal of the pairwise interaction matrix is zero.
# ==============================================================================
build_output_mixed_mrf = function(spec, raw) {
  d = spec$data
  v = spec$variables
  pr = spec$prior
  s = spec$sampler

  p = d$num_discrete
  q = d$num_continuous
  num_variables = d$num_variables
  data_columnnames = d$data_columnnames
  disc_names = d$data_columnnames_discrete
  cont_names = d$data_columnnames_continuous
  disc_idx = d$discrete_indices
  cont_idx = d$continuous_indices
  is_ordinal = v$is_ordinal
  num_categories = d$num_categories
  edge_selection = pr$edge_selection

  # Keep the raw chains for the Z-ratio alarm suite: it needs the untouched
  # indicator layout and the per-chain zratio block, both dropped by the
  # normalization below.
  zratio_chains = if(identical(pr$precision_graph_prior, "hierarchical")) {
    raw
  } else {
    NULL
  }

  # --- Compute index layout in flat parameter vector --------------------------
  layout = compute_mixed_parameter_indices(
    num_thresholds = spec$precomputed$num_thresholds,
    p = p,
    q = q
  )
  nt = layout$num_thresholds
  main_idx = layout$main_idx
  pairwise_idx = layout$pairwise_idx

  # --- Indicator index layout -------------------------------------------------
  # C++ indicator vector: [Gxx_ut | Gyy_ut | Gxy]
  # All are pairwise, so indicator_samples maps directly to pairwise order:
  # Discrete, continuous, cross edges --- same order as pairwise_idx above.

  # --- Normalize raw output per chain -----------------------------------------
  raw = lapply(raw, function(chain) {
    samples_t = t(chain$samples)
    res = list(
      main_samples     = samples_t[, main_idx, drop = FALSE],
      pairwise_samples = samples_t[, pairwise_idx, drop = FALSE],
      userInterrupt    = isTRUE(chain$userInterrupt),
      chain_id         = chain$chain_id
    )
    if(!is.null(chain$indicator_samples)) {
      res$indicator_samples = t(chain$indicator_samples)
    }
    if(!is.null(chain$allocation_samples)) {
      res$allocations = t(chain$allocation_samples)
    }
    if(!is.null(chain$inclusion_parameter_samples)) {
      res$inclusion_parameter = as.numeric(chain$inclusion_parameter_samples)
    }
    attach_diagnostic_traces(res, chain)
  })

  # --- Parameter names --------------------------------------------------------
  # Main effect names (in internal order: discrete first, continuous second)
  names_main = character()
  for(si in seq_len(p)) {
    if(is_ordinal[si]) {
      cats = ordinal_threshold_labels(
        num_categories[si],
        if(!is.null(d$category_levels)) d$category_levels[[si]] else NULL
      )
      names_main = c(names_main, paste0(disc_names[si], " (", cats, ")"))
    } else {
      names_main = c(
        names_main,
        paste0(disc_names[si], " (linear)"),
        paste0(disc_names[si], " (quadratic)")
      )
    }
  }
  for(ji in seq_len(q)) {
    names_main = c(names_main, paste0(cont_names[ji], " (mean)"))
  }
  for(ji in seq_len(q)) {
    names_main = c(names_main, paste0(cont_names[ji], " (precision diag)"))
  }

  # Pairwise edge names --- internal order, mapped to original column names
  # We need a mapping from internal index to original variable name
  # Internal variables: [disc_1, ..., disc_p, cont_1, ..., cont_q]
  # Their original names: c(disc_names, cont_names)
  all_internal_names = c(disc_names, cont_names)

  edge_names = character()
  # Discrete-discrete edges
  if(p > 1) {
    for(i in seq_len(p - 1)) {
      for(j in seq(i + 1, p)) {
        edge_names = c(
          edge_names,
          paste0(disc_names[i], "-", disc_names[j])
        )
      }
    }
  }
  # Continuous-continuous edges (off-diagonal)
  if(q > 1) {
    for(i in seq_len(q - 1)) {
      for(j in seq(i + 1, q)) {
        edge_names = c(
          edge_names,
          paste0(cont_names[i], "-", cont_names[j])
        )
      }
    }
  }
  # Cross edges (discrete-continuous)
  if(p > 0 && q > 0) {
    for(i in seq_len(p)) {
      for(j in seq_len(q)) {
        edge_names = c(
          edge_names,
          paste0(disc_names[i], "-", cont_names[j])
        )
      }
    }
  }

  # --- Lazy MCMC diagnostics cache --------------------------------------------
  cache = new.env(parent = emptyenv())
  cache$raw = raw
  cache$edge_selection = edge_selection
  cache$names_main = names_main
  cache$edge_names = edge_names
  cache$is_continuous = FALSE
  cache$model_type = "mixed_mrf"
  cache$summaries_computed = FALSE

  # --- Compute posterior means from raw samples (cheap) -----------------------
  pooled_main = do.call(rbind, lapply(raw, function(ch) ch$main_samples))
  pooled_pair = do.call(rbind, lapply(raw, function(ch) ch$pairwise_samples))
  main_means = colMeans(pooled_main)
  pair_means = colMeans(pooled_pair)

  # Split main_means into true main effects and quadratic (precision diagonal)
  n_main = nt + q # thresholds + continuous means
  n_quad = layout$num_quadratic # precision diagonal entries

  cache$n_main = n_main
  cache$n_quad = n_quad

  results = list()

  # --- Edge selection summaries (deferred) ------------------------------------
  edge_prior = pr$edge_prior
  has_sbm = FALSE
  if(edge_selection) {
    has_sbm = identical(edge_prior, "Stochastic-Block") &&
      "allocations" %in% names(raw[[1]])

    if(has_sbm) {
      sbm_convergence = summarize_alloc_pairs(
        allocations = lapply(raw, `[[`, "allocations"),
        node_names  = all_internal_names
      )
      results$posterior_summary_pairwise_allocations = sbm_convergence$sbm_summary
      co_occur_matrix = sbm_convergence$co_occur_matrix
    }

    if("inclusion_parameter" %in% names(raw[[1]])) {
      results$inclusion_parameter_samples =
        lapply(raw, `[[`, "inclusion_parameter")
    }
  }

  # --- Posterior mean: main ---------------------------------------------------
  # Discrete main effects: p x max_cats matrix (like OMRF)
  num_params_disc = ifelse(is_ordinal, num_categories, 2L)
  max_num_cats = max(num_params_disc)
  pmm_disc = matrix(NA, nrow = p, ncol = max_num_cats)
  start = 0L
  stop = 0L
  for(si in seq_len(p)) {
    if(is_ordinal[si]) {
      start = stop + 1L
      stop = start + num_categories[si] - 1L
      pmm_disc[si, seq_len(num_categories[si])] = main_means[start:stop]
    } else {
      start = stop + 1L
      stop = start + 1L
      pmm_disc[si, 1:2] = main_means[start:stop]
    }
  }
  rownames(pmm_disc) = disc_names
  colnames(pmm_disc) = paste0("cat (", seq_len(max_num_cats), ")")

  # Continuous main effects: q x 1 matrix (means only)
  pmm_cont = matrix(main_means[nt + seq_len(q)],
    nrow = q, ncol = 1,
    dimnames = list(cont_names, "mean")
  )

  results$posterior_mean_main = list(
    discrete = pmm_disc,
    continuous = pmm_cont
  )

  # --- Posterior mean: associations (all blocks) -----------------------------
  dn = list(data_columnnames, data_columnnames)
  results$posterior_mean_pairwise = fill_mixed_symmetric(
    pair_means, p, q, disc_idx, cont_idx, dn
  )

  # --- Residual variance (continuous diagonal) --------------------------------
  # C++ stores negative association diagonal; residual variance = -1/(2*diag).
  # Average per-sample inversions to avoid Jensen's inequality bias.
  pooled_quad = pooled_main[, nt + q + seq_len(q), drop = FALSE]
  rv = colMeans(-1 / (2 * pooled_quad))
  names(rv) = cont_names
  results$posterior_mean_residual_variance = rv

  # --- Posterior mean: indicator -----------------------------------------------
  if(edge_selection) {
    pooled_ind = do.call(rbind, lapply(raw, function(ch) ch$indicator_samples))
    indicator_means = colMeans(pooled_ind)
    results$posterior_mean_indicator = fill_mixed_symmetric(
      indicator_means, p, q, disc_idx, cont_idx, dn
    )

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
  # NULL placeholders ensure names(fit) lists these fields for easybgm compat.
  # Use list(NULL) because results$x = NULL removes the element in R.
  results["posterior_summary_main"] = list(NULL)
  results["posterior_summary_pairwise"] = list(NULL)
  if(edge_selection) {
    results["posterior_summary_indicator"] = list(NULL)
  }
  results$cache = cache
  class(results) = "bgms"

  # --- raw_samples ------------------------------------------------------------
  alloc_names = if(identical(edge_prior, "Stochastic-Block")) {
    all_internal_names
  } else {
    NULL
  }
  results$raw_samples = build_raw_samples_list(
    raw, edge_selection, edge_prior, names_main, edge_names, alloc_names
  )

  # --- Sampler diagnostics ----------------------------------------------------
  results = attach_sampler_diagnostics(
    results, raw, s$update_method, s$nuts_max_depth,
    names_main = names_main, names_pairwise = edge_names,
    target_accept = s$target_accept
  )

  # --- Z-ratio alarm suite (hierarchical spec on the continuous block) ---------
  # The indicator vector is [Gxx upper, Gyy upper, Gxy row-major], both
  # triangles without diagonals; the audit reads the Gyy segment.
  if(!is.null(zratio_chains)) {
    zc = zratio_cell_constants(
      pr$delta, pr$pairwise_scale, pr$scale_rate, pr$scale_eta
    )
    results$zratio_diag = summarize_zratio_diagnostics(
      zratio_chains, zc,
      num_nodes = q,
      seed = s$seed,
      verbose = TRUE,
      layout_offset = as.integer(p * (p - 1) / 2),
      layout_diag = FALSE
    )
  }

  results$.bgm_spec = spec
  if(needs_easybgm_s3_compat()) {
    results
  } else {
    s3_list_to_bgms(results)
  }
}

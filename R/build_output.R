# ==============================================================================
# build_output: assemble fit objects from bgm_spec + raw sampler output
# ==============================================================================
#
# This file holds build_output() (the thin dispatcher) and the shared helpers
# (raw-sample assembly, diagnostic/SBM attachment, parameter-index and label
# helpers). The per-family builders live in sibling files (cleanup S4):
#   build_output_bgm()       --- build_output_bgm.R       (unified GGM + OMRF)
#   build_output_mixed_mrf() --- build_output_mixed_mrf.R (Mixed MRF)
#   build_output_compare()   --- build_output_compare.R   (bgmCompare)
#
# All builders use build_arguments() from build_arguments.R for the $arguments list.
# ==============================================================================


# ------------------------------------------------------------------
# fill_mixed_symmetric
# ------------------------------------------------------------------
# Fills a symmetric (p+q)x(p+q) matrix from a flat vector of edge
# values stored in discrete-discrete / continuous-continuous / cross
# block order. Used for both pairwise means and indicator means in
# the mixed MRF output builder.
#
# @param values     Flat numeric vector of edge values.
# @param p          Number of discrete variables.
# @param q          Number of continuous variables.
# @param disc_idx   Integer vector mapping discrete 1:p to original columns.
# @param cont_idx   Integer vector mapping continuous 1:q to original columns.
# @param dimnames   List of row/colnames for the result matrix.
#
# Returns: Symmetric matrix with values placed in original-column order.
# ------------------------------------------------------------------
fill_mixed_symmetric = function(values, p, q, disc_idx, cont_idx, dimnames) {
  n = length(dimnames[[1]])
  mat = matrix(0, nrow = n, ncol = n, dimnames = dimnames)
  idx = 0L

  # Discrete-discrete block (upper triangle)
  if(p > 1) {
    for(i in seq_len(p - 1)) {
      for(j in seq(i + 1, p)) {
        idx = idx + 1L
        oi = disc_idx[i]
        oj = disc_idx[j]
        mat[oi, oj] = values[idx]
        mat[oj, oi] = values[idx]
      }
    }
  }

  # Continuous-continuous block (upper triangle)
  if(q > 1) {
    for(i in seq_len(q - 1)) {
      for(j in seq(i + 1, q)) {
        idx = idx + 1L
        oi = cont_idx[i]
        oj = cont_idx[j]
        mat[oi, oj] = values[idx]
        mat[oj, oi] = values[idx]
      }
    }
  }

  # Cross block (all p x q pairs)
  if(p > 0 && q > 0) {
    for(i in seq_len(p)) {
      for(j in seq_len(q)) {
        idx = idx + 1L
        oi = disc_idx[i]
        oj = cont_idx[j]
        mat[oi, oj] = values[idx]
        mat[oj, oi] = values[idx]
      }
    }
  }

  mat
}


# ------------------------------------------------------------------
# compute_mixed_parameter_indices
# ------------------------------------------------------------------
# Computes slice indices for the mixed MRF flat parameter vector.
# Groups main-effect indices (discrete thresholds, continuous means),
# quadratic-effect indices (continuous precision diagonal), and pairwise indices
# (discrete edges, continuous off-diagonal, cross edges).
#
# @param num_thresholds  Total number of discrete threshold parameters.
# @param p               Number of discrete variables.
# @param q               Number of continuous variables.
#
# Returns: List with components num_thresholds, num_quadratic, main_idx,
#   pairwise_idx.
# ------------------------------------------------------------------
compute_mixed_parameter_indices = function(num_thresholds, p, q) {
  nt = num_thresholds
  nxx = as.integer(p * (p - 1) / 2)
  nyy_total = as.integer(q * (q + 1) / 2)
  nyy_offdiag = as.integer(q * (q - 1) / 2)
  nxy = as.integer(p * q)

  # Offsets in the flat vector (1-based)
  main_discrete_start = 1L
  main_discrete_end = nt
  pairwise_discrete_start = nt + 1L
  pairwise_discrete_end = nt + nxx
  main_continuous_start = nt + nxx + 1L
  main_continuous_end = nt + nxx + q
  pairwise_cross_start = nt + nxx + q + 1L
  pairwise_cross_end = nt + nxx + q + nxy
  pairwise_continuous_start = nt + nxx + q + nxy + 1L

  # Continuous diagonal vs off-diagonal within the continuous block
  precision_diag_within = integer(q)
  precision_offdiag_within = integer(nyy_offdiag)
  k_diag = 0L
  k_off = 0L
  pos = 0L
  for(i in seq_len(q)) {
    for(j in i:q) {
      pos = pos + 1L
      if(i == j) {
        k_diag = k_diag + 1L
        precision_diag_within[k_diag] = pos
      } else {
        k_off = k_off + 1L
        precision_offdiag_within[k_off] = pos
      }
    }
  }
  precision_diag_abs = pairwise_continuous_start - 1L + precision_diag_within
  precision_offdiag_abs = pairwise_continuous_start - 1L + precision_offdiag_within

  # Main: discrete thresholds + continuous means
  # Quadratic: precision diagonal (not a main effect)
  main_idx = c(
    seq(main_discrete_start, main_discrete_end),
    seq(main_continuous_start, main_continuous_end),
    precision_diag_abs
  )

  # Pairwise: discrete + precision off-diagonal + cross
  pairwise_idx = c(
    if(nxx > 0) seq(pairwise_discrete_start, pairwise_discrete_end) else integer(0),
    precision_offdiag_abs,
    if(nxy > 0) seq(pairwise_cross_start, pairwise_cross_end) else integer(0)
  )

  list(
    num_thresholds = nt,
    num_quadratic = q,
    main_idx = main_idx,
    pairwise_idx = pairwise_idx
  )
}


# ------------------------------------------------------------------
# build_raw_samples_list
# ------------------------------------------------------------------
# Assembles the $raw_samples list shared by all output builders.
#
# @param raw              Per-chain list (normalized).
# @param edge_selection   Logical.
# @param edge_prior       Character string naming the edge prior.
# @param names_main       Character vector of main-effect parameter names.
# @param edge_names       Character vector of edge parameter names.
# @param allocation_names Optional character vector; when non-NULL, added
#                         to $parameter_names$allocations.
#
# Returns: List with main, pairwise, indicator, allocations, nchains,
#          niter, parameter_names.
# ------------------------------------------------------------------
build_raw_samples_list = function(raw, edge_selection, edge_prior,
                                  names_main, edge_names,
                                  allocation_names = NULL) {
  list(
    main = lapply(raw, function(chain) chain$main_samples),
    pairwise = lapply(raw, function(chain) chain$pairwise_samples),
    indicator = if(edge_selection) {
      lapply(raw, function(chain) chain$indicator_samples)
    } else {
      NULL
    },
    allocations = if(edge_selection &&
      identical(edge_prior, "Stochastic-Block") &&
      "allocations" %in% names(raw[[1]])) {
      lapply(raw, `[[`, "allocations")
    } else {
      NULL
    },
    nchains = length(raw),
    niter = nrow(raw[[1]]$main_samples),
    parameter_names = list(
      main = names_main,
      pairwise = edge_names,
      indicator = if(edge_selection) edge_names else NULL,
      allocations = allocation_names
    )
  )
}


# ------------------------------------------------------------------
# attach_diagnostic_traces
# ------------------------------------------------------------------
# Copy NUTS / adaptive-Metropolis diagnostic traces from a raw C++
# chain onto the normalized chain list `res`, applying the trailing
# "__" naming convention expected by summarize_nuts_diagnostics().
# Each trace is attached only when present
# on the raw chain. Returns `res` with the traces attached.
#
# @param res    Normalized chain list under construction.
# @param chain  Raw per-chain C++ output.
#
# Returns: `res`, with any available diagnostic traces attached.
# ------------------------------------------------------------------
attach_diagnostic_traces = function(res, chain) {
  if(!is.null(chain$treedepth)) res[["treedepth__"]] = chain$treedepth
  if(!is.null(chain$divergent)) res[["divergent__"]] = chain$divergent
  if(!is.null(chain$energy)) res[["energy__"]] = chain$energy
  if(!is.null(chain$accept_prob)) res[["accept_prob__"]] = chain$accept_prob
  if(!is.null(chain$am_accept_prob)) res[["am_accept_prob__"]] = chain$am_accept_prob
  res
}


# ------------------------------------------------------------------
# attach_sampler_diagnostics
# ------------------------------------------------------------------
# Attach the sampler-specific diagnostics block to a results list:
# NUTS diagnostics under "nuts", adaptive-Metropolis diagnostics under
# "adaptive-metropolis", nothing otherwise. The three output builders
# differ only in the parameter-name vectors they pass to the AM
# summary, so those are arguments.
#
# @param results        Results list under construction.
# @param raw            Normalized per-chain list.
# @param update_method  Sampler name from spec$sampler$update_method.
# @param nuts_max_depth  Max tree depth (for the NUTS summary).
# @param names_main     Main-effect parameter names (for the AM summary).
# @param names_pairwise Pairwise parameter names (for the AM summary).
# @param target_accept  Target acceptance (for the AM summary).
#
# Returns: `results`, with $nuts_diag or $am_diag attached as applicable.
# ------------------------------------------------------------------
attach_sampler_diagnostics = function(results, raw, update_method,
                                      nuts_max_depth, names_main,
                                      names_pairwise, target_accept) {
  if(update_method == "nuts") {
    results$nuts_diag = summarize_nuts_diagnostics(
      raw,
      nuts_max_depth = nuts_max_depth
    )
  } else if(update_method == "adaptive-metropolis") {
    results$am_diag = summarize_am_diagnostics(
      raw,
      names_main = names_main,
      names_pairwise = names_pairwise,
      target_accept = target_accept
    )
  }
  results
}


# ------------------------------------------------------------------
# attach_sbm_posterior_summary
# ------------------------------------------------------------------
# Attach the Stochastic-Block posterior allocation summary (mean,
# mode, and block count) to a results list. All three output builders
# call posterior_summary_SBM() identically; they differ only in the
# `arguments` they pass (bgm/mixed forward build_arguments(spec);
# compare hand-builds list(dirichlet_alpha, lambda)), so that is an
# argument here.
#
# @param results    Results list under construction.
# @param raw        Normalized per-chain list (each carrying $allocations).
# @param arguments  Argument list forwarded to posterior_summary_SBM().
#
# Returns: `results`, with the three SBM allocation fields attached.
# ------------------------------------------------------------------
attach_sbm_posterior_summary = function(results, raw, arguments) {
  sbm_summary = posterior_summary_SBM(
    allocations = lapply(raw, `[[`, "allocations"),
    arguments = arguments
  )
  results$posterior_mean_allocations = sbm_summary$allocations_mean
  results$posterior_mode_allocations = sbm_summary$allocations_mode
  results$posterior_num_blocks = sbm_summary$blocks
  results
}


# ------------------------------------------------------------------
# needs_easybgm_s3_compat
# ------------------------------------------------------------------
# Returns TRUE when easybgm is loaded at a version that overwrites
# class(fit) and uses .subset2 directly, both of which are
# incompatible with S7 objects. In that case the builder returns a
# plain S3 list instead of converting to S7.
#
# easybgm >= 0.5.0 uses extractor functions and no longer overwrites
# class(fit), so it is S7-compatible.
# ------------------------------------------------------------------
needs_easybgm_s3_compat = function() {
  if(!"easybgm" %in% loadedNamespaces()) {
    return(FALSE)
  }
  ebgm_version = utils::packageVersion("easybgm")
  if(ebgm_version < "0.5.0") {
    warning(
      "easybgm ", ebgm_version, " is not compatible with S7-based bgms objects. ",
      "Running in S3 compatibility mode. ",
      "Please update easybgm to version 0.5.0 or later.",
      call. = FALSE
    )
    return(TRUE)
  }
  FALSE
}


# ==============================================================================
# build_output()  --- dispatcher
# ==============================================================================
build_output = function(spec, raw) {
  stopifnot(inherits(spec, "bgm_spec"))

  switch(spec$model_type,
    ggm       = build_output_bgm(spec, raw),
    omrf      = build_output_bgm(spec, raw),
    mixed_mrf = build_output_mixed_mrf(spec, raw),
    compare   = build_output_compare(spec, raw),
    stop("Unknown model_type: ", spec$model_type)
  )
}


# ==============================================================================
# ordinal_threshold_labels()
# ------------------------------------------------------------------------------
# Labels for an ordinal variable's category thresholds (recoded categories
# 1..K), in the user's original category scale when the fit carries a recode
# map (category_levels_v). Two map forms are supported:
#
#   * OMRF / mixed-MRF: an UNNAMED vector of the sorted original training
#     values, length K+1 with position 1 the reference category. The originals
#     for categories 1..K are returned directly.
#   * bgmCompare: a NAMED lookup, names = original values, values = the final
#     0-based category. Cross-group collapse makes this many-to-one, so for
#     each final category 1..K the contributing original values are joined with
#     "/" (e.g. "4/5" when originals 4 and 5 merged into one category).
#
# Without a usable map (NULL, or an unnamed vector whose length != K+1) the
# rescored indices 1..K are returned.
ordinal_threshold_labels = function(num_categories_v, category_levels_v = NULL) {
  if(is.null(category_levels_v)) {
    return(seq_len(num_categories_v))
  }
  if(!is.null(names(category_levels_v))) {
    finals = as.integer(category_levels_v)
    originals = names(category_levels_v)
    return(vapply(
      seq_len(num_categories_v),
      function(k) paste(originals[finals == k], collapse = "/"),
      character(1)
    ))
  }
  if(length(category_levels_v) == num_categories_v + 1L) {
    return(category_levels_v[-1L]) # drop the reference category (recoded 0)
  }
  seq_len(num_categories_v)
}

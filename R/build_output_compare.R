# build_output_compare() — bgmCompare output builder (+ parameter-name generator)
#
# Split out of build_output.R (cleanup S4).


# ==============================================================================
# build_output_compare()
# ==============================================================================
build_output_compare = function(spec, raw) {
  d = spec$data
  v = spec$variables
  p = spec$prior
  s = spec$sampler
  pc = spec$precomputed

  num_variables = d$num_variables
  num_groups = d$num_groups
  data_columnnames = d$data_columnnames
  num_categories = d$num_categories
  is_ordinal_variable = v$is_ordinal
  difference_selection = p$difference_selection
  difference_prior = p$difference_prior

  # Normalize SBM allocation samples (variables x iter -> iter x variables)
  # and the sampled difference slab scale.
  raw = lapply(raw, function(chain) {
    if(!is.null(chain$allocation_samples)) {
      chain$allocations = t(chain$allocation_samples)
    }
    if(!is.null(chain$scale_samples)) {
      chain$difference_scale = as.numeric(chain$scale_samples)
    }
    chain
  })

  # --- Parameter names --------------------------------------------------------
  names_all = generate_param_names_bgmCompare(
    data_columnnames    = data_columnnames,
    num_categories      = num_categories,
    is_ordinal_variable = is_ordinal_variable,
    num_variables       = num_variables,
    num_groups          = num_groups,
    category_levels     = d$category_levels
  )

  # --- Lazy MCMC diagnostics cache --------------------------------------------
  cache = new.env(parent = emptyenv())
  cache$raw = raw
  cache$difference_selection = difference_selection
  cache$main_effect_indices = pc$main_effect_indices
  cache$pairwise_effect_indices = pc$pairwise_effect_indices
  cache$num_variables = num_variables
  cache$num_groups = num_groups
  cache$names_all = names_all
  cache$model_type = "compare"
  cache$summaries_computed = FALSE

  # --- Compute baseline means from raw samples (cheap) -------------------------
  num_main = pc$main_effect_indices[nrow(pc$main_effect_indices), 2] + 1
  num_pair = pc$pairwise_effect_indices[
    nrow(pc$pairwise_effect_indices),
    nrow(pc$pairwise_effect_indices) - 1
  ] + 1
  pooled_main_bl = do.call(rbind, lapply(raw, function(ch) ch$main_samples[, 1:num_main, drop = FALSE]))
  pooled_pair_bl = do.call(rbind, lapply(raw, function(ch) ch$pairwise_samples[, 1:num_pair, drop = FALSE]))
  main_bl_means = colMeans(pooled_main_bl)
  pair_bl_means = colMeans(pooled_pair_bl)

  results = list()

  # --- Posterior mean: main baseline ------------------------------------------
  num_params = ifelse(is_ordinal_variable, num_categories, 2L)
  max_num_categories = max(num_params, na.rm = TRUE)

  pmm = matrix(NA, nrow = num_variables, ncol = max_num_categories)
  start = 0L
  stop = 0L
  for(vi in seq_len(num_variables)) {
    if(is_ordinal_variable[vi]) {
      start = stop + 1L
      stop = start + num_categories[vi] - 1L
      pmm[vi, seq_len(num_categories[vi])] = main_bl_means[start:stop]
    } else {
      start = stop + 1L
      stop = start + 1L
      pmm[vi, 1:2] = main_bl_means[start:stop]
    }
  }
  results$posterior_mean_main_baseline = pmm
  rownames(results$posterior_mean_main_baseline) = data_columnnames
  colnames(results$posterior_mean_main_baseline) =
    paste0("cat (", seq_len(ncol(pmm)), ")")

  # --- Posterior mean: associations baseline ----------------------------------
  results$posterior_mean_pairwise_baseline = matrix(0,
    nrow = num_variables, ncol = num_variables,
    dimnames = list(data_columnnames, data_columnnames)
  )
  results$posterior_mean_pairwise_baseline[
    lower.tri(results$posterior_mean_pairwise_baseline)
  ] = pair_bl_means
  results$posterior_mean_pairwise_baseline =
    results$posterior_mean_pairwise_baseline +
    t(results$posterior_mean_pairwise_baseline)

  # --- Posterior mean: main differences ---------------------------------------
  n_main_total = ncol(raw[[1]]$main_samples)
  pooled_main_diff = do.call(rbind, lapply(raw, function(ch) {
    ch$main_samples[, (num_main + 1):n_main_total, drop = FALSE]
  }))
  main_diff_means = colMeans(pooled_main_diff)
  num_contrasts = num_groups - 1L

  main_diff_list = vector("list", num_contrasts)
  for(g in seq_len(num_contrasts)) {
    offset = (g - 1L) * num_main
    pmm_diff = matrix(NA, nrow = num_variables, ncol = max_num_categories)
    start = 0L
    stop = 0L
    for(vi in seq_len(num_variables)) {
      if(is_ordinal_variable[vi]) {
        start = stop + 1L
        stop = start + num_categories[vi] - 1L
        pmm_diff[vi, seq_len(num_categories[vi])] =
          main_diff_means[offset + start:stop]
      } else {
        start = stop + 1L
        stop = start + 1L
        pmm_diff[vi, 1:2] = main_diff_means[offset + start:stop]
      }
    }
    rownames(pmm_diff) = data_columnnames
    colnames(pmm_diff) = paste0("cat (", seq_len(ncol(pmm_diff)), ")")
    main_diff_list[[g]] = pmm_diff
  }
  names(main_diff_list) = paste0("diff", seq_len(num_contrasts))

  if(num_contrasts == 1L) {
    results$posterior_mean_main_differences = main_diff_list[[1L]]
  } else {
    results$posterior_mean_main_differences = main_diff_list
  }

  # --- Posterior mean: associations differences --------------------------------
  n_pair_total = ncol(raw[[1]]$pairwise_samples)
  pooled_pair_diff = do.call(rbind, lapply(raw, function(ch) {
    ch$pairwise_samples[, (num_pair + 1):n_pair_total, drop = FALSE]
  }))
  pair_diff_means = colMeans(pooled_pair_diff)

  pair_diff_list = vector("list", num_contrasts)
  for(g in seq_len(num_contrasts)) {
    offset = (g - 1L) * num_pair
    mat = matrix(0,
      nrow = num_variables, ncol = num_variables,
      dimnames = list(data_columnnames, data_columnnames)
    )
    mat[lower.tri(mat)] = pair_diff_means[offset + seq_len(num_pair)]
    mat = mat + t(mat)
    pair_diff_list[[g]] = mat
  }
  names(pair_diff_list) = paste0("diff", seq_len(num_contrasts))

  if(num_contrasts == 1L) {
    results$posterior_mean_pairwise_differences = pair_diff_list[[1L]]
  } else {
    results$posterior_mean_pairwise_differences = pair_diff_list
  }

  # --- SBM allocation summaries (Stochastic-Block difference prior) -----------
  has_sbm = difference_selection &&
    identical(difference_prior, "Stochastic-Block") &&
    "allocations" %in% names(raw[[1]])

  if(has_sbm) {
    sbm_convergence = summarize_alloc_pairs(
      allocations = lapply(raw, `[[`, "allocations"),
      node_names  = data_columnnames
    )
    results$posterior_summary_pairwise_allocations = sbm_convergence$sbm_summary
    results$posterior_mean_coclustering_matrix = sbm_convergence$co_occur_matrix

    results = attach_sbm_posterior_summary(
      results, raw,
      list(dirichlet_alpha = p$dirichlet_alpha, lambda = p$lambda)
    )
  }

  # --- raw_samples ------------------------------------------------------------
  results$raw_samples = list(
    main = lapply(raw, function(chain) chain$main_samples),
    pairwise = lapply(raw, function(chain) chain$pairwise_samples),
    indicator = if(difference_selection) {
      lapply(raw, function(chain) chain$indicator_samples)
    } else {
      NULL
    },
    rb_inclusion = if(difference_selection) {
      lapply(raw, function(chain) chain$rb_inclusion_samples)
    } else {
      NULL
    },
    rb_counts = if(difference_selection) {
      lapply(raw, function(chain) chain$rb_counts)
    } else {
      NULL
    },
    allocations = if(has_sbm) {
      lapply(raw, `[[`, "allocations")
    } else {
      NULL
    },
    difference_scale = if("difference_scale" %in% names(raw[[1]])) {
      lapply(raw, `[[`, "difference_scale")
    } else {
      NULL
    },
    nchains = length(raw),
    niter = nrow(raw[[1]]$main_samples),
    parameter_names = names_all
  )

  # Difference-scale draws are collected independently of difference selection:
  # the random-scale hyperprior can be active with a dense difference model.
  if("difference_scale" %in% names(raw[[1]])) {
    results$difference_scale_samples =
      lapply(raw, `[[`, "difference_scale")
  }

  # --- arguments + class ------------------------------------------------------
  results$arguments = build_arguments(spec)
  # Report the number of chains actually kept; failed chains are dropped
  # upstream, so raw holds only the survivors.
  results$arguments$num_chains = length(raw)
  # NULL placeholders ensure names(fit) lists these fields for easybgm compat.
  # Use list(NULL) because results$x = NULL removes the element in R.
  results["posterior_summary_main_baseline"] = list(NULL)
  results["posterior_summary_pairwise_baseline"] = list(NULL)
  results["posterior_summary_main_differences"] = list(NULL)
  results["posterior_summary_pairwise_differences"] = list(NULL)
  if(difference_selection) {
    results["posterior_summary_indicator"] = list(NULL)
  }
  results$cache = cache
  class(results) = "bgmCompare"

  # --- Sampler diagnostics ----------------------------------------------------
  results = attach_sampler_diagnostics(
    results, raw, s$update_method, s$nuts_max_depth,
    names_main = c(names_all$main_baseline, names_all$main_diff),
    names_pairwise = c(names_all$pairwise_baseline, names_all$pairwise_diff),
    target_accept = s$target_accept
  )

  results$.bgm_spec = spec
  if(needs_easybgm_s3_compat()) {
    results
  } else {
    s3_list_to_bgmCompare(results)
  }
}


# ==============================================================================
# generate_param_names_bgmCompare()
# ==============================================================================
#
# Build parameter names for bgmCompare models. Used by build_output_compare().
#
# @param data_columnnames  Character vector: variable names.
# @param num_categories  Integer vector: max category per variable.
# @param is_ordinal_variable  Logical vector: TRUE = ordinal, FALSE = BC.
# @param num_variables  Integer: number of variables.
# @param num_groups  Integer: number of groups.
# @param category_levels  Optional list of per-variable recode maps (named
#   lookups from original category value to final 0-based category). When
#   supplied, ordinal threshold names use the original category scale.
#
# Returns: named list with main_baseline, main_diff, pairwise_baseline,
#   pairwise_diff, and indicators character vectors.
# ==============================================================================
generate_param_names_bgmCompare = function(
  data_columnnames,
  num_categories,
  is_ordinal_variable,
  num_variables,
  num_groups,
  category_levels = NULL
) {
  # --- main baselines
  names_main_baseline = character()
  for(v in seq_len(num_variables)) {
    if(is_ordinal_variable[v]) {
      cats = ordinal_threshold_labels(
        num_categories[v],
        if(!is.null(category_levels)) category_levels[[v]] else NULL
      )
      names_main_baseline = c(
        names_main_baseline,
        paste0(data_columnnames[v], " (", cats, ")")
      )
    } else {
      names_main_baseline = c(
        names_main_baseline,
        paste0(data_columnnames[v], " (linear)"),
        paste0(data_columnnames[v], " (quadratic)")
      )
    }
  }

  # --- main differences
  names_main_diff = character()
  for(g in 2:num_groups) {
    for(v in seq_len(num_variables)) {
      if(is_ordinal_variable[v]) {
        cats = ordinal_threshold_labels(
          num_categories[v],
          if(!is.null(category_levels)) category_levels[[v]] else NULL
        )
        names_main_diff = c(
          names_main_diff,
          paste0(data_columnnames[v], " (diff", g - 1, "; ", cats, ")")
        )
      } else {
        names_main_diff = c(
          names_main_diff,
          paste0(data_columnnames[v], " (diff", g - 1, "; linear)"),
          paste0(data_columnnames[v], " (diff", g - 1, "; quadratic)")
        )
      }
    }
  }

  # --- pairwise baselines
  names_pairwise_baseline = character()
  for(i in 1:(num_variables - 1)) {
    for(j in (i + 1):num_variables) {
      names_pairwise_baseline = c(
        names_pairwise_baseline,
        paste0(data_columnnames[i], "-", data_columnnames[j])
      )
    }
  }

  # --- pairwise differences
  names_pairwise_diff = character()
  for(g in 2:num_groups) {
    for(i in 1:(num_variables - 1)) {
      for(j in (i + 1):num_variables) {
        names_pairwise_diff = c(
          names_pairwise_diff,
          paste0(data_columnnames[i], "-", data_columnnames[j], " (diff", g - 1, ")")
        )
      }
    }
  }

  # --- indicators
  generate_indicator_names = function(data_columnnames) {
    V = length(data_columnnames)
    out = character()
    for(i in seq_len(V)) {
      # main (diagonal)
      out = c(out, paste0(data_columnnames[i], " (main)"))
      # then all pairs with i as the first index
      if(i < V) {
        for(j in seq.int(i + 1L, V)) {
          out = c(out, paste0(data_columnnames[i], "-", data_columnnames[j], " (pairwise)"))
        }
      }
    }
    # optional sanity check: length must be V*(V+1)/2
    stopifnot(length(out) == V * (V + 1L) / 2L)
    out
  }
  names_indicators = generate_indicator_names(data_columnnames)

  list(
    main_baseline = names_main_baseline,
    main_diff = names_main_diff,
    pairwise_baseline = names_pairwise_baseline,
    pairwise_diff = names_pairwise_diff,
    indicators = names_indicators
  )
}

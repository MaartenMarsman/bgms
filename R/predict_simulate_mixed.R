# Mixed-MRF helpers for the simulate()/predict() S3 methods
#
# Split out of simulate_predict.R (cleanup S4).


# ==============================================================================
#   Mixed MRF Simulation Helper
# ==============================================================================

# ------------------------------------------------------------------
# simulate_bgms_mixed
# ------------------------------------------------------------------
# Simulation implementation for mixed MRF models (called from simulate.bgms).
#
# @param object     Fitted bgms object (mixed MRF).
# @param nsim       Number of observations to simulate.
# @param seed       Random seed.
# @param method     "posterior-mean" or "posterior-sample".
# @param ndraws     Number of posterior draws (for posterior-sample).
# @param arguments  Output of extract_arguments().
# @param iter       Gibbs burn-in iterations.
# @param cores      Number of threads.
# @param progress_type  Progress bar type.
#
# Returns: matrix (posterior-mean) or list of matrices (posterior-sample).
# ------------------------------------------------------------------
simulate_bgms_mixed = function(object, nsim, seed, method, ndraws,
                               arguments, iter, cores, progress_type) {
  p = arguments$num_discrete
  q = arguments$num_continuous
  data_columnnames = arguments$data_columnnames
  disc_idx = arguments$discrete_indices
  cont_idx = arguments$continuous_indices
  num_categories = arguments$num_categories
  is_ordinal = arguments$is_ordinal
  category_levels = arguments$category_levels
  blume_capel_shift = arguments$blume_capel_shift
  baseline_category_disc = arguments$baseline_category

  disc_variable_type = ifelse(is_ordinal, "ordinal", "blume-capel")

  bc = integer(p)
  for(s in seq_len(p)) {
    if(is_ordinal[s]) {
      bc[s] = 0L
    } else {
      bc[s] = as.integer(baseline_category_disc[s])
    }
  }

  if(method == "posterior-mean") {
    params = build_mixed_params_mean(object, arguments)

    seed = check_seed(seed)

    result = sample_mixed_mrf_gibbs(
      num_states = as.integer(nsim),
      pairwise_disc_r = params$pairwise_disc,
      pairwise_cross_r = params$pairwise_cross,
      pairwise_cont_r = params$pairwise_cont,
      mux_r = params$mux,
      muy_r = params$muy,
      num_categories_r = as.integer(num_categories),
      variable_type_r = disc_variable_type,
      baseline_category_r = as.integer(bc),
      iter = as.integer(iter),
      seed = seed
    )

    result$x = recode_simulated_to_original(
      result$x, category_levels, blume_capel_shift
    )
    out = combine_mixed_result(result, disc_idx, cont_idx, data_columnnames)
    return(out)
  } else {
    sample_info = split_mixed_raw_samples(object, arguments)

    total_draws = sample_info$total_draws
    if(is.null(ndraws)) ndraws = total_draws
    ndraws = min(ndraws, total_draws)

    if(!is.null(seed)) set.seed(seed)
    draw_indices = sample.int(total_draws, ndraws)

    results = run_mixed_simulation_parallel(
      mux_samples = sample_info$mux_samples,
      disc_samples = sample_info$disc_samples,
      muy_samples = sample_info$muy_samples,
      cont_samples = sample_info$cont_samples,
      cross_samples = sample_info$cross_samples,
      draw_indices = as.integer(draw_indices),
      num_states = as.integer(nsim),
      p = as.integer(p),
      q = as.integer(q),
      num_categories = as.integer(num_categories),
      variable_type_r = disc_variable_type,
      baseline_category = as.integer(bc),
      iter = as.integer(iter),
      nThreads = cores,
      seed = check_seed(seed),
      progress_type = progress_type
    )

    for(i in seq_along(results)) {
      results[[i]]$x = recode_simulated_to_original(
        results[[i]]$x, category_levels, blume_capel_shift
      )
      results[[i]] = combine_mixed_result(
        results[[i]], disc_idx, cont_idx, data_columnnames
      )
    }

    return(results)
  }
}


# ==============================================================================
#   Mixed MRF Prediction Helper
# ==============================================================================

# ------------------------------------------------------------------
# predict_bgms_mixed
# ------------------------------------------------------------------
# Prediction implementation for mixed MRF models (called from predict.bgms).
#
# @param object       Fitted bgms object (mixed MRF).
# @param newdata      n x (p+q) matrix of observed data.
# @param predict_vars 1-based indices into the combined variable list.
# @param arguments    Output of extract_arguments().
# @param type         "probabilities" or "response".
# @param method       "posterior-mean" or "posterior-sample".
# @param ndraws       Number of posterior draws (for posterior-sample).
#
# Returns: Named list of prediction matrices.
# ------------------------------------------------------------------
predict_bgms_mixed = function(object, newdata, predict_vars, arguments,
                              type, method, ndraws) {
  p = arguments$num_discrete
  q = arguments$num_continuous
  data_columnnames = arguments$data_columnnames
  disc_idx = arguments$discrete_indices
  cont_idx = arguments$continuous_indices
  num_categories = arguments$num_categories
  is_ordinal = arguments$is_ordinal
  category_levels = arguments$category_levels
  blume_capel_shift = arguments$blume_capel_shift
  baseline_category_disc = arguments$baseline_category

  disc_variable_type = ifelse(is_ordinal, "ordinal", "blume-capel")

  bc = integer(p)
  for(s in seq_len(p)) {
    if(is_ordinal[s]) {
      bc[s] = 0L
    } else {
      bc[s] = as.integer(baseline_category_disc[s])
    }
  }

  # Split newdata into discrete and continuous parts. The discrete columns are
  # recoded from their original category values to the internal 0-based codes
  # (the mapping bgm() applied to the training data), as in the OMRF predict
  # path; unobserved categories warn and yield NA predictions.
  x_data = as.matrix(newdata[, disc_idx, drop = FALSE])
  x_data = recode_data_for_prediction(
    x_data, num_categories, is_ordinal, category_levels, blume_capel_shift
  )
  storage.mode(x_data) = "integer"
  y_data = as.matrix(newdata[, cont_idx, drop = FALSE])
  storage.mode(y_data) = "double"

  # Map user predict_vars (1-based, original order) to internal 0-based indices
  # Internal layout: [discrete_0..p-1, continuous_p..p+q-1]
  internal_predict_vars = integer(length(predict_vars))
  for(k in seq_along(predict_vars)) {
    orig_col = predict_vars[k]
    disc_match = match(orig_col, disc_idx)
    cont_match = match(orig_col, cont_idx)
    if(!is.na(disc_match)) {
      internal_predict_vars[k] = disc_match - 1L
    } else if(!is.na(cont_match)) {
      internal_predict_vars[k] = p + cont_match - 1L
    } else {
      stop(
        "Variable index ", orig_col,
        " not found in discrete or continuous indices."
      )
    }
  }

  compute_one_draw = function(pairwise_disc, pairwise_cross, pairwise_cont, mux, muy) {
    compute_conditional_mixed(
      x_observations = x_data,
      y_observations = y_data,
      predict_vars = as.integer(internal_predict_vars),
      pairwise_disc = pairwise_disc,
      pairwise_cross = pairwise_cross,
      pairwise_cont = pairwise_cont,
      mux = mux,
      muy = muy,
      num_categories = as.integer(num_categories),
      variable_type = disc_variable_type,
      baseline_category = as.integer(bc)
    )
  }

  if(method == "posterior-mean") {
    params = build_mixed_params_mean(object, arguments)
    raw_result = compute_one_draw(
      params$pairwise_disc, params$pairwise_cross, params$pairwise_cont, params$mux, params$muy
    )

    probs = format_mixed_predictions(
      raw_result, predict_vars, internal_predict_vars,
      p, num_categories, data_columnnames
    )
  } else {
    sample_info = split_mixed_raw_samples(object, arguments)
    total_draws = sample_info$total_draws
    if(is.null(ndraws)) ndraws = total_draws
    ndraws = min(ndraws, total_draws)

    draw_indices = sample.int(total_draws, ndraws)

    all_results = vector("list", ndraws)
    for(i in seq_len(ndraws)) {
      params_i = build_mixed_params_row(
        sample_info, draw_indices[i], p, q, num_categories, is_ordinal
      )
      all_results[[i]] = compute_one_draw(
        params_i$pairwise_disc, params_i$pairwise_cross, params_i$pairwise_cont, params_i$mux, params_i$muy
      )
    }

    # Average predictions across draws
    num_pv = length(predict_vars)
    probs = vector("list", num_pv)
    probs_sd = vector("list", num_pv)
    names(probs) = data_columnnames[predict_vars]
    names(probs_sd) = data_columnnames[predict_vars]

    for(k in seq_len(num_pv)) {
      avg = average_draws(all_results, k)
      probs[[k]] = avg$mean
      probs_sd[[k]] = avg$sd
    }

    probs = format_mixed_predictions(
      probs, predict_vars, internal_predict_vars,
      p, num_categories, data_columnnames
    )
    names(probs_sd) = names(probs)
    attr(probs, "sd") = probs_sd
  }

  if(type == "response") {
    return(format_mixed_response(
      probs, predict_vars, internal_predict_vars,
      p, data_columnnames
    ))
  }

  return(probs)
}


# ==============================================================================
#   Mixed MRF Internal Helpers
# ==============================================================================

# ------------------------------------------------------------------
# build_mixed_params_mean
# ------------------------------------------------------------------
# Reconstruct discrete, cross, and continuous interaction matrices,
# plus mux and muy from posterior mean summaries.
#
# @param object     Fitted bgms object (mixed MRF).
# @param arguments  Output of extract_arguments().
#
# Returns: list with pairwise_disc, pairwise_cross, pairwise_cont, mux, muy.
# ------------------------------------------------------------------
build_mixed_params_mean = function(object, arguments) {
  p = arguments$num_discrete
  q = arguments$num_continuous
  disc_idx = arguments$discrete_indices
  cont_idx = arguments$continuous_indices

  pmat = get_posterior_mean(object, "pairwise")

  pairwise_disc = matrix(0, p, p)
  for(i in seq_len(p)) {
    for(j in seq_len(p)) {
      if(i != j) pairwise_disc[i, j] = pmat[disc_idx[i], disc_idx[j]]
    }
  }

  pairwise_cross = matrix(0, p, q)
  for(i in seq_len(p)) {
    for(j in seq_len(q)) {
      pairwise_cross[i, j] = pmat[disc_idx[i], cont_idx[j]]
    }
  }

  pairwise_cont = matrix(0, q, q)
  for(i in seq_len(q)) {
    for(j in seq_len(q)) {
      pairwise_cont[i, j] = pmat[cont_idx[i], cont_idx[j]]
    }
  }
  # Convert residual variance back to association-scale diagonal
  rv = get_posterior_mean(object, "residual_variance")
  for(j in seq_len(q)) {
    pairwise_cont[j, j] = -1 / (2 * rv[j])
  }

  pm_main = get_posterior_mean(object, "main")
  mux = pm_main$discrete
  mux[is.na(mux)] = 0

  muy = as.numeric(pm_main$continuous[, "mean"])

  list(pairwise_disc = pairwise_disc, pairwise_cross = pairwise_cross, pairwise_cont = pairwise_cont, mux = mux, muy = muy)
}



# ------------------------------------------------------------------
# split_mixed_raw_samples
# ------------------------------------------------------------------
# Split raw main and pairwise sample matrices into separate component
# matrices for the C++ parallel simulation worker.
#
# @param object     Fitted bgms object (mixed MRF).
# @param arguments  Output of extract_arguments().
#
# Returns: list with mux_samples, disc_samples, muy_samples,
#   cont_samples, cross_samples, total_draws.
# ------------------------------------------------------------------
split_mixed_raw_samples = function(object, arguments) {
  p = arguments$num_discrete
  q = arguments$num_continuous
  num_categories = arguments$num_categories
  is_ordinal = arguments$is_ordinal

  raw = get_raw_samples(object)
  main_all = do.call(rbind, raw$main)
  pairwise_all = do.call(rbind, raw$pairwise)
  total_draws = nrow(main_all)

  # Main layout: [mux_flat | muy | cont_diag]
  num_mux = sum(ifelse(is_ordinal, num_categories, 2L))
  mux_cols = seq_len(num_mux)
  muy_cols = num_mux + seq_len(q)
  cont_diag_cols = num_mux + q + seq_len(q)

  mux_samples = main_all[, mux_cols, drop = FALSE]
  muy_samples = main_all[, muy_cols, drop = FALSE]
  cont_diag_values = main_all[, cont_diag_cols, drop = FALSE]

  # Pairwise layout: [disc_ut | cont_offdiag | cross]
  nxx = as.integer(p * (p - 1) / 2)
  nyy_offdiag = as.integer(q * (q - 1) / 2)
  nxy = as.integer(p * q)

  cont_off_end = nxx + nyy_offdiag
  cross_end = cont_off_end + nxy

  disc_samples = if(nxx > 0) {
    pairwise_all[, seq_len(nxx), drop = FALSE]
  } else {
    matrix(0, nrow = total_draws, ncol = 0)
  }

  cont_offdiag_values = if(nyy_offdiag > 0) {
    pairwise_all[, (nxx + 1):cont_off_end, drop = FALSE]
  } else {
    matrix(0, nrow = total_draws, ncol = 0)
  }

  cross_samples = if(nxy > 0) {
    pairwise_all[, (cont_off_end + 1):cross_end, drop = FALSE]
  } else {
    matrix(0, nrow = total_draws, ncol = 0)
  }

  # Combine continuous diagonal and off-diagonal into upper-triangle format
  # C++ expects column-major upper-triangle including diagonal
  nyy_total = as.integer(q * (q + 1) / 2)
  cont_samples = matrix(0, nrow = total_draws, ncol = nyy_total)
  diag_pos = 0L
  offdiag_pos = 0L
  out_pos = 0L
  for(col in seq_len(q)) {
    for(row in col:q) {
      out_pos = out_pos + 1L
      if(row == col) {
        diag_pos = diag_pos + 1L
        cont_samples[, out_pos] = cont_diag_values[, diag_pos]
      } else {
        offdiag_pos = offdiag_pos + 1L
        cont_samples[, out_pos] = cont_offdiag_values[, offdiag_pos]
      }
    }
  }

  list(
    mux_samples = mux_samples,
    disc_samples = disc_samples,
    muy_samples = muy_samples,
    cont_samples = cont_samples,
    cross_samples = cross_samples,
    total_draws = total_draws
  )
}



# ------------------------------------------------------------------
# build_mixed_params_row
# ------------------------------------------------------------------
# Reconstruct discrete, cross, and continuous interaction matrices,
# plus mux and muy from a single row of split
# sample matrices (used by predict posterior-sample).
#
# @param sample_info  Output of split_mixed_raw_samples().
# @param row_idx      1-based row index.
# @param p            Number of discrete variables.
# @param q            Number of continuous variables.
# @param num_categories  Categories per discrete variable.
# @param is_ordinal   Logical vector.
#
# Returns: list with pairwise_disc, pairwise_cross, pairwise_cont, mux, muy.
# ------------------------------------------------------------------
build_mixed_params_row = function(sample_info, row_idx,
                                  p, q, num_categories,
                                  is_ordinal) {
  mux_vec = sample_info$mux_samples[row_idx, ]
  num_params_disc = ifelse(is_ordinal, num_categories, 2L)
  max_cats = max(num_params_disc)
  mux = matrix(0, nrow = p, ncol = max_cats)
  pos = 1L
  for(s in seq_len(p)) {
    nc = num_params_disc[s]
    mux[s, seq_len(nc)] = mux_vec[pos:(pos + nc - 1L)]
    pos = pos + nc
  }

  pairwise_disc = matrix(0, p, p)
  if(p > 1) {
    disc_vec = sample_info$disc_samples[row_idx, ]
    idx = 0L
    for(col in seq_len(p - 1)) {
      for(row in (col + 1):p) {
        idx = idx + 1L
        pairwise_disc[row, col] = disc_vec[idx]
        pairwise_disc[col, row] = disc_vec[idx]
      }
    }
  }

  muy = sample_info$muy_samples[row_idx, ]

  cont_vec = sample_info$cont_samples[row_idx, ]
  pairwise_cont = matrix(0, q, q)
  idx = 0L
  for(col in seq_len(q)) {
    for(row in col:q) {
      idx = idx + 1L
      pairwise_cont[row, col] = cont_vec[idx]
      if(row != col) pairwise_cont[col, row] = cont_vec[idx]
    }
  }

  pairwise_cross = matrix(0, p, q)
  if(p > 0 && q > 0) {
    cross_vec = sample_info$cross_samples[row_idx, ]
    idx = 0L
    for(s in seq_len(p)) {
      for(j in seq_len(q)) {
        idx = idx + 1L
        pairwise_cross[s, j] = cross_vec[idx]
      }
    }
  }

  list(pairwise_disc = pairwise_disc, pairwise_cross = pairwise_cross, pairwise_cont = pairwise_cont, mux = mux, muy = muy)
}



# ------------------------------------------------------------------
# combine_mixed_result
# ------------------------------------------------------------------
# Combine $x (n x p integer) and $y (n x q double) into a single
# n x (p+q) matrix in the original column order.
#
# @param result    List with $x and $y matrices.
# @param disc_idx  Original column indices for discrete variables.
# @param cont_idx  Original column indices for continuous variables.
# @param colnames  Original data column names.
#
# Returns: n x (p+q) numeric matrix.
# ------------------------------------------------------------------
combine_mixed_result = function(result, disc_idx, cont_idx, colnames) {
  n = nrow(result$x)
  num_vars = length(disc_idx) + length(cont_idx)
  out = matrix(NA_real_, nrow = n, ncol = num_vars)
  out[, disc_idx] = result$x
  out[, cont_idx] = result$y
  colnames(out) = colnames
  out
}



# ------------------------------------------------------------------
# format_mixed_predictions
# ------------------------------------------------------------------
# Add names and column labels to C++ prediction output.
#
# @param raw_result         List from C++ compute_conditional_mixed.
# @param predict_vars       1-based user-facing variable indices.
# @param internal_predict_vars  0-based internal indices.
# @param p                  Number of discrete variables.
# @param num_categories     Categories per discrete variable.
# @param data_columnnames   Original data column names.
#
# Returns: Named list of prediction matrices.
# ------------------------------------------------------------------
format_mixed_predictions = function(raw_result, predict_vars,
                                    internal_predict_vars, p,
                                    num_categories, data_columnnames) {
  probs = raw_result
  names(probs) = data_columnnames[predict_vars]

  for(k in seq_along(predict_vars)) {
    int_idx = internal_predict_vars[k]
    if(int_idx < p) {
      s = int_idx + 1L
      n_cats = num_categories[s] + 1
      colnames(probs[[k]]) = paste0("cat_", 0:(n_cats - 1))
    } else {
      colnames(probs[[k]]) = c("mean", "sd")
    }
  }

  probs
}



# ------------------------------------------------------------------
# format_mixed_response
# ------------------------------------------------------------------
# Convert probability predictions to point predictions for mixed models.
#
# @param probs              Named list of prediction matrices.
# @param predict_vars       1-based user-facing variable indices.
# @param internal_predict_vars  0-based internal indices.
# @param p                  Number of discrete variables.
# @param data_columnnames   Original data column names.
#
# Returns: n x length(predict_vars) matrix of predicted values.
# ------------------------------------------------------------------
format_mixed_response = function(probs, predict_vars,
                                 internal_predict_vars, p,
                                 data_columnnames) {
  n = nrow(probs[[1]])
  out = matrix(NA_real_, nrow = n, ncol = length(predict_vars))
  colnames(out) = data_columnnames[predict_vars]

  for(k in seq_along(predict_vars)) {
    int_idx = internal_predict_vars[k]
    if(int_idx < p) {
      out[, k] = apply(probs[[k]], 1, which.max) - 1L
    } else {
      out[, k] = probs[[k]][, 1]
    }
  }

  out
}

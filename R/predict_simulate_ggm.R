# GGM helpers for the simulate()/predict() S3 methods
#
# Split out of simulate_predict.R (cleanup S4).


# ==============================================================================
#   GGM Prediction Helpers
# ==============================================================================

# Reconstruct the full precision matrix from association-scale values and
# residual variances.
#
# @param associations p x p symmetric matrix of pairwise
#   associations (zero diagonal).
# @param residual_variance Named numeric vector of residual variances.
#
# @return p x p precision matrix.
reconstruct_precision = function(associations, residual_variance) {
  omega = -2 * associations
  # Excluded edges (NA) have zero precision
  omega[is.na(omega)] = 0
  diag(omega) = 1 / residual_variance
  return(omega)
}



# Reconstruct precision matrix from a single posterior draw.
# GGM raw samples are already on precision scale.
#
# @param pairwise_vec Vector of p*(p-1)/2 off-diagonal precision elements
#   (lower-triangle order).
# @param main_vec Vector of p diagonal precision elements.
# @param p Number of variables.
#
# @return p x p precision matrix (Omega).
build_precision_from_draw = function(pairwise_vec, main_vec, p) {
  omega = matrix(0, nrow = p, ncol = p)
  omega[lower.tri(omega)] = pairwise_vec
  omega = omega + t(omega)
  diag(omega) = main_vec
  return(omega)
}



# ==============================================================================
#   GGM Simulation Helpers
# ==============================================================================

# GGM simulation implementation (called from simulate.bgms).
#
# @param object Fitted bgms object (GGM).
# @param nsim Number of observations to simulate.
# @param seed Random seed.
# @param method "posterior-mean" or "posterior-sample".
# @param ndraws Number of posterior draws (NULL = all).
# @param num_variables Number of variables p.
# @param data_columnnames Character vector of variable names.
# @param cores Number of parallel threads.
# @param progress_type Integer progress type (0/1/2).
#
# @return See simulate.bgms() documentation.
simulate_bgms_ggm = function(object, nsim, seed, method, ndraws,
                             num_variables, data_columnnames,
                             cores, progress_type) {
  # The model is fit on training-mean-centered data; simulate on the original
  # scale by using the stored means. Older fits without stored means fall
  # back to the centered (zero-mean) scale.
  train_means = extract_arguments(object)$column_means
  if(is.null(train_means)) {
    train_means = rep(0, num_variables)
  }

  if(method == "posterior-mean") {
    # Reconstruct precision matrix from off-diagonal + separate diagonal
    precision = reconstruct_precision(
      get_posterior_mean(object, "pairwise"),
      get_posterior_mean(object, "residual_variance")
    )

    # Call simulate_mrf with variable_type = "continuous"
    result = simulate_mrf(
      num_states = nsim,
      num_variables = num_variables,
      pairwise = precision,
      main = train_means,
      variable_type = "continuous",
      seed = seed
    )

    colnames(result) = data_columnnames
    return(result)
  } else {
    # Use posterior samples with parallel processing
    raw = get_raw_samples(object)
    pairwise_samples = do.call(rbind, raw$pairwise)
    main_samples = do.call(rbind, raw$main)

    total_draws = nrow(pairwise_samples)
    if(is.null(ndraws)) {
      ndraws = total_draws
    }
    ndraws = min(ndraws, total_draws)

    # Sample which draws to use
    if(!is.null(seed)) set.seed(seed)
    draw_indices = sample.int(total_draws, ndraws)

    # Call parallel C++ function for GGM
    results = run_ggm_simulation_parallel(
      pairwise_samples = pairwise_samples,
      main_samples = main_samples,
      draw_indices = as.integer(draw_indices),
      num_states = as.integer(nsim),
      num_variables = as.integer(num_variables),
      means = train_means,
      nThreads = cores,
      seed = seed,
      progress_type = progress_type
    )

    # Add column names
    for(i in seq_along(results)) {
      colnames(results[[i]]) = data_columnnames
    }

    return(results)
  }
}


# GGM prediction implementation (called from predict.bgms).
#
# @param object Fitted bgms object (GGM).
# @param newdata n x p numeric matrix of observed continuous data.
# @param predict_vars Integer vector of 1-based variable indices to predict.
# @param data_columnnames Character vector of variable names.
# @param num_variables Number of variables p.
# @param type "probabilities" or "response".
# @param method "posterior-mean" or "posterior-sample".
# @param ndraws Number of posterior draws (NULL = all).
#
# @return See predict.bgms() documentation for GGM return format.
predict_bgms_ggm = function(object, newdata, predict_vars, data_columnnames,
                            num_variables,
                            type, method, ndraws) {
  # Center newdata on the training column means so predictions match the
  # scale the model was fit on. Older fits without stored means fall back to
  # centering newdata by its own means.
  train_means = extract_arguments(object)$column_means
  if(is.null(train_means)) {
    train_means = colMeans(newdata)
  }
  newdata_centered = sweep(newdata, 2, train_means)

  if(method == "posterior-mean") {
    # Reconstruct precision matrix from posterior means
    omega = reconstruct_precision(
      get_posterior_mean(object, "pairwise"),
      get_posterior_mean(object, "residual_variance")
    )

    result = compute_conditional_ggm(
      observations = newdata_centered,
      predict_vars = predict_vars - 1L,
      precision = omega
    )

    # Add names and shift conditional means back to original scale
    names(result) = data_columnnames[predict_vars]
    for(v in seq_along(result)) {
      colnames(result[[v]]) = c("mean", "sd")
      result[[v]][, "mean"] =
        result[[v]][, "mean"] +
        train_means[predict_vars[v]]
    }
  } else {
    # Use posterior samples
    raw = get_raw_samples(object)
    pairwise_samples = do.call(rbind, raw$pairwise)
    main_samples = do.call(rbind, raw$main)

    total_draws = nrow(pairwise_samples)
    if(is.null(ndraws)) {
      ndraws = total_draws
    }
    ndraws = min(ndraws, total_draws)

    draw_indices = sample.int(total_draws, ndraws)

    # Collect predictions from each draw
    all_preds = vector("list", ndraws)

    for(i in seq_len(ndraws)) {
      idx = draw_indices[i]

      omega = build_precision_from_draw(
        pairwise_vec = pairwise_samples[idx, ],
        main_vec = main_samples[idx, ],
        p = num_variables
      )

      preds = compute_conditional_ggm(
        observations = newdata_centered,
        predict_vars = predict_vars - 1L,
        precision = omega
      )

      # Shift conditional means back to original scale
      for(v in seq_along(predict_vars)) {
        preds[[v]][, 1] = preds[[v]][, 1] + train_means[predict_vars[v]]
      }

      all_preds[[i]] = preds
    }

    # Average over draws
    result = vector("list", length(predict_vars))
    result_sd = vector("list", length(predict_vars))
    names(result) = data_columnnames[predict_vars]
    names(result_sd) = data_columnnames[predict_vars]

    for(v in seq_along(predict_vars)) {
      avg = average_draws(all_preds, v)
      result[[v]] = avg$mean
      result_sd[[v]] = avg$sd

      colnames(result[[v]]) = c("mean", "sd")
      colnames(result_sd[[v]]) = c("mean", "sd")
    }

    attr(result, "sd") = result_sd
  }

  if(type == "response") {
    # Return conditional means
    pred_matrix = sapply(result, function(m) m[, "mean"])
    if(is.vector(pred_matrix)) {
      pred_matrix = matrix(pred_matrix, ncol = 1)
    }
    colnames(pred_matrix) = data_columnnames[predict_vars]
    return(pred_matrix)
  }

  return(result)
}

# ==============================================================================
#   MRF Simulation and Prediction — S3 methods + shared helpers
#
#   This file contains the user-facing S3 methods and the helpers shared across
#   model families:
#     - simulate.bgms() / simulate.bgmCompare(): simulate from fitted models
#     - predict.bgms() / predict.bgmCompare(): conditional-probability prediction
#     - shared helpers: expand_variable_type, average_draws, reconstruct_main,
#       recode_data_for_prediction
#
#   Split out (cleanup S4):
#     - simulate_mrf() / mrfSampler()  --> simulate_mrf.R
#     - GGM helpers                    --> predict_simulate_ggm.R
#     - Mixed-MRF helpers              --> predict_simulate_mixed.R
# ==============================================================================


# ------------------------------------------------------------------
# expand_variable_type
# ------------------------------------------------------------------
# Recycles a scalar variable_type to length num_variables.
#
# @param variable_type  Character vector (possibly length 1).
# @param num_variables  Target length.
#
# Returns: Character vector of length num_variables.
# ------------------------------------------------------------------
expand_variable_type = function(variable_type, num_variables) {
  if(length(variable_type) == 1) {
    rep(variable_type, num_variables)
  } else {
    variable_type
  }
}


# ==============================================================================
#   simulate.bgms() - S3 Method for Simulating from Fitted Models
# ==============================================================================

#' Simulate Data from a Fitted bgms Model
#'
#' @description
#' Generates new observations from the Markov Random Field model using the
#' estimated parameters from a fitted \code{bgms} object. Supports ordinal,
#' Blume-Capel, continuous (GGM), and mixed MRF models.
#'
#' @param object An object of class \code{bgms}.
#' @param nsim Number of observations to simulate. Default: \code{500}.
#' @param seed Optional random seed for reproducibility.
#' @param method Character string specifying which parameter estimates to use:
#'   \describe{
#'     \item{\code{"posterior-mean"}}{Use posterior mean parameters (faster,
#'       single simulation).}
#'     \item{\code{"posterior-sample"}}{Sample from posterior draws, producing
#'       one dataset per draw (accounts for parameter uncertainty). This method
#'       uses parallel processing when \code{cores > 1}.}
#'   }
#' @param ndraws Number of posterior draws to use when
#'   \code{method = "posterior-sample"}. If \code{NULL},
#'   uses all available draws.
#' @param iter Number of Gibbs iterations for equilibration before collecting
#'   samples. Default: \code{1000}.
#' @param cores Number of CPU cores for parallel execution when
#'   \code{method = "posterior-sample"}.
#'   Default: \code{parallel::detectCores()}.
#' @param display_progress Character string specifying the type of progress bar.
#'   Options: \code{"per-chain"}, \code{"total"}, \code{"none"}.
#'   Default: \code{"per-chain"}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' If \code{method = "posterior-mean"}: A matrix with \code{nsim} rows and
#' \code{p} columns containing simulated observations.
#'
#' If \code{method = "posterior-sample"}: A list of matrices, one per posterior
#' draw, each with \code{nsim} rows and \code{p} columns.
#'
#' Discrete columns are returned on the original category scale of the training
#' data (the values supplied to \code{bgm()}), so the output can be passed
#' straight to \code{predict()}. For mixed MRF models, discrete columns contain
#' non-negative integers and continuous columns contain real-valued
#' observations, ordered as in the original data.
#'
#' @details
#' This function uses the estimated interaction and threshold
#' parameters to generate new data via Gibbs sampling. When
#' \code{method = "posterior-sample"}, parameter uncertainty is
#' propagated to the simulated data by using different
#' posterior draws. Parallel processing is available for this method via the
#' \code{cores} argument.
#'
#' @seealso \code{\link{predict.bgms}} for computing conditional probabilities,
#'   \code{\link{simulate_mrf}} for simulation with user-specified parameters.
#' @family prediction
#'
#' @examples
#' \donttest{
#' # Fit a model
#' fit = bgm(x = Wenchuan[, 1:5], chains = 2)
#'
#' # Simulate 100 new observations using posterior means
#' new_data = simulate(fit, nsim = 100)
#'
#' # Simulate with parameter uncertainty (10 datasets)
#' new_data_list = simulate(
#'   fit,
#'   nsim = 100,
#'   method = "posterior-sample", ndraws = 10
#' )
#'
#' # Use parallel processing for faster simulation
#' new_data_list = simulate(fit,
#'   nsim = 100, method = "posterior-sample",
#'   ndraws = 100, cores = 2
#' )
#' }
#'
#' @importFrom stats simulate
#' @export
simulate.bgms = function(object,
                         nsim = 500,
                         seed = NULL,
                         method = c("posterior-mean", "posterior-sample"),
                         ndraws = NULL,
                         iter = 1000,
                         cores = parallel::detectCores(),
                         display_progress = c("per-chain", "total", "none"),
                         ...) {
  method = match.arg(method)
  progress_type = progress_type_from_display_progress(display_progress)

  # Validate cores
  check_positive_integer(cores, "cores")
  cores = normalize_parallel_cores(cores)

  # Setting the seed
  seed = check_seed(seed)

  # Extract model information

  arguments = extract_arguments(object)
  num_variables = arguments$num_variables
  num_categories = arguments$num_categories
  variable_type = arguments$variable_type
  data_columnnames = arguments$data_columnnames

  # Handle variable_type
  variable_type = expand_variable_type(variable_type, num_variables)

  # Get baseline_category (for Blume-Capel variables)
  baseline_category = arguments$baseline_category
  if(is.null(baseline_category)) {
    baseline_category = rep(0L, num_variables)
  }

  # ============================================================================
  #   GGM (continuous) path
  # ============================================================================
  if(isTRUE(arguments$is_continuous)) {
    return(simulate_bgms_ggm(
      object = object,
      nsim = nsim,
      seed = seed,
      method = method,
      ndraws = ndraws,
      num_variables = num_variables,
      data_columnnames = data_columnnames,
      cores = cores,
      progress_type = progress_type
    ))
  }

  # ============================================================================
  #   Mixed MRF (discrete + continuous) path
  # ============================================================================
  if(isTRUE(arguments$is_mixed)) {
    return(simulate_bgms_mixed(
      object = object,
      nsim = nsim,
      seed = seed,
      method = method,
      ndraws = ndraws,
      arguments = arguments,
      iter = iter,
      cores = cores,
      progress_type = progress_type
    ))
  }

  # ============================================================================
  #   OMRF (ordinal / Blume-Capel) path
  # ============================================================================

  if(method == "posterior-mean") {
    # Use posterior mean parameters
    pairwise = get_posterior_mean(object, "pairwise")
    main = get_posterior_mean(object, "main")

    # Set R's RNG for simulate_mrf
    if(!is.null(seed)) set.seed(seed)

    # Call simulate_mrf
    result = simulate_mrf(
      num_states = nsim,
      num_variables = num_variables,
      num_categories = num_categories,
      pairwise = pairwise,
      main = main,
      variable_type = variable_type,
      baseline_category = baseline_category,
      iter = iter
    )

    result = recode_simulated_to_original(
      result, arguments$category_levels, arguments$blume_capel_shift
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

    # Call parallel C++ function
    results = run_simulation_parallel(
      pairwise_samples = pairwise_samples,
      main_samples = main_samples,
      draw_indices = as.integer(draw_indices),
      num_states = as.integer(nsim),
      num_variables = as.integer(num_variables),
      num_categories = as.integer(num_categories),
      variable_type_r = variable_type,
      baseline_category = as.integer(baseline_category),
      iter = as.integer(iter),
      nThreads = cores,
      seed = seed,
      progress_type = progress_type
    )

    # Map codes back to the original scale and add column names.
    for(i in seq_along(results)) {
      results[[i]] = recode_simulated_to_original(
        results[[i]], arguments$category_levels, arguments$blume_capel_shift
      )
      colnames(results[[i]]) = data_columnnames
    }

    return(results)
  }
}


# ============================================================
#   simulate.bgmCompare() - S3 Method for Group-Comparison
# ============================================================

#' Simulate Data from a Fitted bgmCompare Model
#'
#' @description
#' Generates new observations from the Markov Random Field model for a
#' specified group using the estimated parameters from a fitted
#' \code{bgmCompare} object.
#'
#' @param object An object of class \code{bgmCompare}.
#' @param nsim Number of observations to simulate. Default: \code{500}.
#' @param seed Optional random seed for reproducibility.
#' @param group Integer specifying which group to simulate from (1 to
#'   number of groups). Required argument.
#' @param method Character string specifying which parameter estimates to use:
#'   \describe{
#'     \item{\code{"posterior-mean"}}{Use posterior mean parameters (faster,
#'       single simulation).}
#'   }
#' @param iter Number of Gibbs iterations for equilibration before collecting
#'   samples. Default: \code{1000}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A matrix with \code{nsim} rows and \code{p} columns containing
#'   simulated observations for the specified group.
#'
#' @details
#' Group-specific parameters are obtained by applying the projection matrix
#' to convert baseline parameters and differences into group-level estimates:
#' \code{group_param = baseline + projection[group, ] \%*\% differences}.
#'
#' The function then uses these group-specific interaction and threshold
#' parameters to generate new data via Gibbs sampling.
#'
#' @seealso \code{\link{simulate.bgms}} for simulating from single-group models,
#'   \code{\link{predict.bgmCompare}} for computing conditional probabilities.
#' @family prediction
#'
#' @examples
#' \donttest{
#' # Fit a comparison model
#' x = Boredom[Boredom$language == "fr", 2:6]
#' y = Boredom[Boredom$language != "fr", 2:6]
#' fit = bgmCompare(x, y, chains = 2)
#'
#' # Simulate 100 observations from group 1
#' new_data_g1 = simulate(fit, nsim = 100, group = 1)
#'
#' # Simulate 100 observations from group 2
#' new_data_g2 = simulate(fit, nsim = 100, group = 2)
#' }
#'
#' @export
simulate.bgmCompare = function(object,
                               nsim = 500,
                               seed = NULL,
                               group,
                               method = c("posterior-mean"),
                               iter = 1000,
                               ...) {
  method = match.arg(method)

  # Validate group argument
  if(missing(group)) {
    stop(
      "Argument 'group' is required. ",
      "Specify which group to simulate from ",
      "(1 to num_groups)."
    )
  }

  arguments = extract_arguments(object)
  num_groups = arguments$num_groups

  invalid_group = !is.numeric(group) || length(group) != 1 ||
    is.na(group) || group < 1 || group > num_groups
  if(invalid_group) {
    stop(sprintf(
      "Argument 'group' must be an integer between 1 and %d.",
      num_groups
    ))
  }
  group = as.integer(group)

  # Setting the seed
  seed = check_seed(seed)

  # Extract model information
  num_variables = arguments$num_variables
  num_categories = arguments$num_categories
  is_ordinal = arguments$is_ordinal_variable
  data_columnnames = arguments$data_columnnames

  # Determine variable_type from is_ordinal
  variable_type = ifelse(is_ordinal, "ordinal", "blume-capel")

  # Get baseline_category (for Blume-Capel variables)
  baseline_category = arguments$baseline_category
  if(is.null(baseline_category)) {
    baseline_category = rep(0L, num_variables)
  }

  if(method == "posterior-mean") {
    # Extract group-specific parameters using projection
    group_params = extract_group_params(object)

    main_group = group_params$main_effects_groups[, group]
    pairwise_group = group_params$pairwise_effects_groups[, group]

    # Reconstruct threshold matrix (variable_type is ifelse(is_ordinal,
    # "ordinal", "blume-capel"), so reconstruct_main's blume-capel branch
    # matches the per-variable parameter count exactly).
    main = reconstruct_main(
      main_group, num_variables, num_categories, variable_type
    )

    # Reconstruct interaction matrix
    pairwise = matrix(0, nrow = num_variables, ncol = num_variables)
    pairwise[lower.tri(pairwise)] = pairwise_group
    pairwise = pairwise + t(pairwise)

    # Set R's RNG for simulate_mrf
    set.seed(seed)

    # Call simulate_mrf
    result = simulate_mrf(
      num_states = nsim,
      num_variables = num_variables,
      num_categories = num_categories,
      pairwise = pairwise,
      main = main,
      variable_type = variable_type,
      baseline_category = baseline_category,
      iter = iter
    )

    result = recode_simulated_to_original(
      result, arguments$category_levels, arguments$blume_capel_shift
    )
    colnames(result) = data_columnnames
    return(result)
  }
}


# ------------------------------------------------------------------------------
# average_draws()
# ------------------------------------------------------------------------------
# Posterior-sample prediction averaging. Given a per-draw list (one element
# per posterior draw, each itself a list indexed by predicted variable), stack
# variable `v`'s per-draw matrices into an n x k x ndraws array and reduce to
# posterior mean and sd matrices. Each per-draw matrix is n_obs x k, so the
# array dims are read from the matrix itself (rather than re-deriving n_obs /
# k at each call site, which had drifted across the three predict methods).
#
# @param per_draw_list  List of length ndraws; element i is a per-variable
#   list of n_obs x k prediction matrices.
# @param v              Index of the predicted variable to average.
#
# Returns: list(mean = n_obs x k matrix, sd = n_obs x k matrix).
# ------------------------------------------------------------------------------
average_draws = function(per_draw_list, v) {
  var_mats = lapply(per_draw_list, `[[`, v)
  arr = array(
    unlist(var_mats),
    dim = c(nrow(var_mats[[1]]), ncol(var_mats[[1]]), length(var_mats))
  )
  list(
    mean = apply(arr, c(1, 2), mean),
    sd = apply(arr, c(1, 2), sd)
  )
}


# ==============================================================================
#   predict.bgms() - S3 Method for Conditional Probability Prediction
# ==============================================================================

#' Predict Conditional Probabilities from a Fitted bgms Model
#'
#' @description
#' Computes conditional probability distributions for one or more variables
#' given the observed values of other variables in the data. Supports ordinal,
#' Blume-Capel, continuous (GGM), and mixed MRF models.
#'
#' @param object An object of class \code{bgms}.
#' @param newdata A matrix or data frame with \code{n} rows and \code{p} columns
#'   containing the observed data. Must have the same variables (columns) as
#'   the original data used to fit the model.
#' @param variables Which variables to predict. Can be:
#'   \itemize{
#'     \item A character vector of variable names
#'     \item An integer vector of column indices
#'     \item \code{NULL} (default) to predict all variables
#'   }
#' @param type Character string specifying the type of prediction:
#'   \describe{
#'     \item{\code{"probabilities"}}{Return the full conditional probability
#'       distribution for each variable and observation.}
#'     \item{\code{"response"}}{Return the predicted category (mode of the
#'       conditional distribution).}
#'   }
#' @param method Character string specifying which parameter estimates to use:
#'   \describe{
#'     \item{\code{"posterior-mean"}}{Use posterior mean parameters.}
#'     \item{\code{"posterior-sample"}}{Average predictions
#'       over posterior draws.}
#'   }
#' @param ndraws Number of posterior draws to use when
#'   \code{method = "posterior-sample"}. If \code{NULL},
#'   uses all available draws.
#' @param seed Optional random seed for reproducibility when
#'   \code{method = "posterior-sample"}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' \strong{Ordinal models:}
#'
#' For \code{type = "probabilities"}: A named list with one element per
#' predicted variable. Each element is a matrix with \code{n} rows and
#' \code{num_categories + 1} columns containing
#' \eqn{P(X_j = c | X_{-j})}{P(X_j = c | X_-j)}
#' for each observation and category.
#'
#' For \code{type = "response"}: A matrix with \code{n} rows and
#' \code{length(variables)} columns containing predicted categories.
#'
#' When \code{method = "posterior-sample"}, probabilities are averaged over
#' posterior draws, and an attribute \code{"sd"} is included containing the
#' standard deviation across draws.
#'
#' \strong{GGM (continuous) models:}
#'
#' For \code{type = "probabilities"}: A named list with one element per
#' predicted variable. Each element is a matrix with \code{n} rows and
#' 2 columns (\code{"mean"} and \code{"sd"}) containing the conditional
#' Gaussian parameters \eqn{E(X_j | X_{-j})}{E(X_j | X_{-j})} and
#' \eqn{\text{SD}(X_j | X_{-j})}{SD(X_j | X_{-j})}.
#'
#' For \code{type = "response"}: A matrix with \code{n} rows and
#' \code{length(variables)} columns containing conditional means.
#'
#' When \code{method = "posterior-sample"}, conditional parameters are
#' averaged over posterior draws, and an attribute \code{"sd"} is included.
#'
#' \strong{Mixed MRF models:}
#'
#' For mixed models, the return list contains elements for both discrete and
#' continuous predicted variables. Discrete variables return probability
#' matrices (as in ordinal models); continuous variables return conditional
#' mean and SD matrices (as in GGM models).
#'
#' @details
#' For each observation, the function computes the conditional distribution
#' of the target variable(s) given the observed values of all other variables.
#' This is the same conditional distribution used internally by the Gibbs
#' sampler.
#'
#' For GGM (continuous) models, the conditional distribution of
#' \eqn{X_j | X_{-j}}{X_j | X_{-j}} is Gaussian with mean
#' \eqn{-\omega_{jj}^{-1} \sum_{k \neq j}
#' \omega_{jk} x_k}{-omega_jj^{-1} sum_{k != j} omega_jk x_k}
#' and variance \eqn{\omega_{jj}^{-1}}{omega_jj^{-1}}, where \eqn{\Omega}{Omega}
#' is the precision matrix.
#'
#' @seealso \code{\link{simulate.bgms}} for generating new data from the model.
#' @family prediction
#'
#' @examples
#' \donttest{
#' # Fit a model
#' fit = bgm(x = Wenchuan[, 1:5], chains = 2)
#'
#' # Compute conditional probabilities for all variables
#' probs = predict(fit, newdata = Wenchuan[1:10, 1:5])
#'
#' # Predict the first variable only
#' probs_v1 = predict(fit, newdata = Wenchuan[1:10, 1:5], variables = 1)
#'
#' # Get predicted categories
#' pred_class = predict(fit, newdata = Wenchuan[1:10, 1:5], type = "response")
#' }
#'
#' @importFrom stats predict
#' @export
predict.bgms = function(object,
                        newdata,
                        variables = NULL,
                        type = c("probabilities", "response"),
                        method = c("posterior-mean", "posterior-sample"),
                        ndraws = NULL,
                        seed = NULL,
                        ...) {
  type = match.arg(type)
  method = match.arg(method)

  # Setting the seed (for R's RNG used by sample.int for draw selection)
  if(!is.null(seed)) {
    seed = check_seed(seed)
    set.seed(seed)
  }

  # Validate newdata
  if(missing(newdata)) {
    stop(
      "Argument 'newdata' is required. ",
      "Provide the data for which to ",
      "compute predictions."
    )
  }

  if(!inherits(newdata, "matrix") && !inherits(newdata, "data.frame")) {
    stop("'newdata' must be a matrix or data frame.")
  }

  if(inherits(newdata, "data.frame")) {
    newdata = data.matrix(newdata)
  }

  # Extract model information
  arguments = extract_arguments(object)
  num_variables = arguments$num_variables
  num_categories = arguments$num_categories
  variable_type = arguments$variable_type
  data_columnnames = arguments$data_columnnames

  # Validate dimensions

  if(ncol(newdata) != num_variables) {
    stop(paste0(
      "'newdata' must have ", num_variables,
      " columns (same as fitted model), ",
      "but has ", ncol(newdata), "."
    ))
  }

  # Handle variable_type
  variable_type = expand_variable_type(variable_type, num_variables)

  # Get baseline_category
  baseline_category = arguments$baseline_category
  if(is.null(baseline_category)) {
    baseline_category = rep(0L, num_variables)
  }

  # Convert variable_type to is_ordinal logical vector
  is_ordinal = variable_type != "blume-capel"

  # Determine which variables to predict
  if(is.null(variables)) {
    predict_vars = seq_len(num_variables)
  } else if(is.character(variables)) {
    predict_vars = match(variables, data_columnnames)
    if(anyNA(predict_vars)) {
      stop(
        "Variable names not found: ",
        paste(
          variables[is.na(predict_vars)],
          collapse = ", "
        )
      )
    }
  } else {
    predict_vars = as.integer(variables)
    if(any(predict_vars < 1 | predict_vars > num_variables)) {
      stop("Variable indices must be between 1 and ", num_variables)
    }
  }

  # ============================================================================
  #   GGM (continuous) path
  # ============================================================================
  if(isTRUE(arguments$is_continuous)) {
    return(predict_bgms_ggm(
      object = object,
      newdata = newdata,
      predict_vars = predict_vars,
      data_columnnames = data_columnnames,
      num_variables = num_variables,
      type = type,
      method = method,
      ndraws = ndraws
    ))
  }

  # ============================================================================
  #   Mixed MRF (discrete + continuous) path
  # ============================================================================
  if(isTRUE(arguments$is_mixed)) {
    return(predict_bgms_mixed(
      object = object,
      newdata = newdata,
      predict_vars = predict_vars,
      arguments = arguments,
      type = type,
      method = method,
      ndraws = ndraws
    ))
  }

  # ============================================================================
  #   OMRF (ordinal) path
  # ============================================================================

  # Recode data to 0-based integers (matching what bgm() did to the training
  # data) using the stored recode map when available.
  newdata_recoded = recode_data_for_prediction(
    newdata, is_ordinal,
    category_levels = arguments$category_levels,
    blume_capel_shift = arguments$blume_capel_shift
  )

  if(method == "posterior-mean") {
    # Use posterior mean parameters
    pairwise = get_posterior_mean(object, "pairwise")
    main = get_posterior_mean(object, "main")

    probs = compute_conditional_probs(
      observations = newdata_recoded,
      predict_vars = predict_vars - 1L, # C++ uses 0-based indexing
      pairwise = pairwise,
      main = main,
      num_categories = num_categories,
      variable_type = variable_type,
      baseline_category = baseline_category
    )

    # Add names
    names(probs) = data_columnnames[predict_vars]
    for(v in seq_along(probs)) {
      var_idx = predict_vars[v]
      n_cats = num_categories[var_idx] + 1
      colnames(probs[[v]]) = paste0("cat_", 0:(n_cats - 1))
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

    # Collect probabilities from each draw
    all_probs = vector("list", ndraws)

    for(i in seq_len(ndraws)) {
      idx = draw_indices[i]

      # Reconstruct interaction matrix
      pairwise = matrix(0, nrow = num_variables, ncol = num_variables)
      pairwise[lower.tri(pairwise)] = pairwise_samples[idx, ]
      pairwise = pairwise + t(pairwise)

      # Reconstruct threshold matrix
      main = reconstruct_main(
        main_samples[idx, ],
        num_variables,
        num_categories,
        variable_type
      )

      all_probs[[i]] = compute_conditional_probs(
        observations = newdata_recoded,
        predict_vars = predict_vars - 1L,
        pairwise = pairwise,
        main = main,
        num_categories = num_categories,
        variable_type = variable_type,
        baseline_category = baseline_category
      )
    }

    # Average over draws
    probs = vector("list", length(predict_vars))
    probs_sd = vector("list", length(predict_vars))
    names(probs) = data_columnnames[predict_vars]
    names(probs_sd) = data_columnnames[predict_vars]

    for(v in seq_along(predict_vars)) {
      avg = average_draws(all_probs, v)
      probs[[v]] = avg$mean
      probs_sd[[v]] = avg$sd

      var_idx = predict_vars[v]
      n_cats = num_categories[var_idx] + 1
      colnames(probs[[v]]) = paste0("cat_", 0:(n_cats - 1))
      colnames(probs_sd[[v]]) = paste0("cat_", 0:(n_cats - 1))
    }

    attr(probs, "sd") = probs_sd
  }

  if(type == "response") {
    # Return predicted categories (mode)
    pred_matrix = sapply(probs, function(p) {
      apply(p, 1, which.max) - 1L # Convert to 0-based category
    })
    if(is.vector(pred_matrix)) {
      pred_matrix = matrix(pred_matrix, ncol = 1)
    }
    colnames(pred_matrix) = data_columnnames[predict_vars]
    return(pred_matrix)
  }

  return(probs)
}


# ==============================================================================
#   predict.bgmCompare() - S3 Method for Group-Comparison Models
# ==============================================================================

#' Predict Conditional Probabilities from a Fitted bgmCompare Model
#'
#' @description
#' Computes conditional probability distributions for one or more variables
#' given the observed values of other variables in the data, using
#' group-specific parameters from a \code{bgmCompare} model.
#'
#' @param object An object of class \code{bgmCompare}.
#' @param newdata A matrix or data frame with \code{n} rows and \code{p} columns
#'   containing the observed data. Must have the same variables (columns) as
#'   the original data used to fit the model.
#' @param group Integer specifying which group's parameters to use for
#'   prediction (1 to number of groups). Required argument.
#' @param variables Which variables to predict. Can be:
#'   \itemize{
#'     \item A character vector of variable names
#'     \item An integer vector of column indices
#'     \item \code{NULL} (default) to predict all variables
#'   }
#' @param type Character string specifying the type of prediction:
#'   \describe{
#'     \item{\code{"probabilities"}}{Return the full conditional probability
#'       distribution for each variable and observation.}
#'     \item{\code{"response"}}{Return the predicted category (mode of the
#'       conditional distribution).}
#'   }
#' @param method Character string specifying which parameter estimates to use:
#'   \describe{
#'     \item{\code{"posterior-mean"}}{Use posterior mean parameters.}
#'   }
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' For \code{type = "probabilities"}: A named list with one
#' element per predicted variable. Each element is a matrix with
#' \code{n} rows and \code{num_categories + 1} columns containing
#' \eqn{P(X_j = c | X_{-j})}{P(X_j = c | X_-j)}
#' for each observation and category.
#'
#' For \code{type = "response"}: A matrix with \code{n} rows and
#' \code{length(variables)} columns containing predicted categories.
#'
#' @details
#' Group-specific parameters are obtained by applying the projection matrix
#' to convert baseline parameters and differences into group-level estimates.
#' The function then computes the conditional distribution of target variables
#' given the observed values of all other variables.
#'
#' @seealso \code{\link{predict.bgms}} for predicting
#'   from single-group models,
#'   \code{\link{simulate.bgmCompare}} for simulating
#'   from group-comparison models.
#' @family prediction
#'
#' @examples
#' \donttest{
#' # Fit a comparison model
#' x = Boredom[Boredom$language == "fr", 2:6]
#' y = Boredom[Boredom$language != "fr", 2:6]
#' fit = bgmCompare(x, y, chains = 2)
#'
#' # Predict conditional probabilities using group 1 parameters
#' probs_g1 = predict(fit, newdata = x[1:10, ], group = 1)
#'
#' # Predict responses using group 2 parameters
#' pred_g2 = predict(fit, newdata = y[1:10, ], group = 2, type = "response")
#' }
#'
#' @export
predict.bgmCompare = function(object,
                              newdata,
                              group,
                              variables = NULL,
                              type = c("probabilities", "response"),
                              method = c("posterior-mean"),
                              ...) {
  type = match.arg(type)
  method = match.arg(method)

  # Validate group argument
  if(missing(group)) {
    stop(
      "Argument 'group' is required. ",
      "Specify which group's parameters ",
      "to use (1 to num_groups)."
    )
  }

  arguments = extract_arguments(object)
  num_groups = arguments$num_groups

  invalid_group = !is.numeric(group) || length(group) != 1 ||
    is.na(group) || group < 1 || group > num_groups
  if(invalid_group) {
    stop(sprintf(
      "Argument 'group' must be an integer between 1 and %d.",
      num_groups
    ))
  }
  group = as.integer(group)

  # Validate newdata
  if(missing(newdata)) {
    stop(
      "Argument 'newdata' is required. ",
      "Provide the data for which to ",
      "compute predictions."
    )
  }

  if(!inherits(newdata, "matrix") && !inherits(newdata, "data.frame")) {
    stop("'newdata' must be a matrix or data frame.")
  }

  if(inherits(newdata, "data.frame")) {
    newdata = data.matrix(newdata)
  }

  # Extract model information
  num_variables = arguments$num_variables
  num_categories = arguments$num_categories
  is_ordinal = arguments$is_ordinal_variable
  data_columnnames = arguments$data_columnnames

  # Validate dimensions
  if(ncol(newdata) != num_variables) {
    stop(paste0(
      "'newdata' must have ", num_variables,
      " columns (same as fitted model), ",
      "but has ", ncol(newdata), "."
    ))
  }

  # Determine variable_type from is_ordinal
  variable_type = ifelse(is_ordinal, "ordinal", "blume-capel")

  # Get baseline_category (for Blume-Capel variables)
  baseline_category = arguments$baseline_category
  if(is.null(baseline_category)) {
    baseline_category = rep(0L, num_variables)
  }

  # Determine which variables to predict
  if(is.null(variables)) {
    predict_vars = seq_len(num_variables)
  } else if(is.character(variables)) {
    predict_vars = match(variables, data_columnnames)
    if(anyNA(predict_vars)) {
      stop(
        "Variable names not found: ",
        paste(
          variables[is.na(predict_vars)],
          collapse = ", "
        )
      )
    }
  } else {
    predict_vars = as.integer(variables)
    if(any(predict_vars < 1 | predict_vars > num_variables)) {
      stop("Variable indices must be between 1 and ", num_variables)
    }
  }

  # Recode data to 0-based integers using the stored recode map when available.
  newdata_recoded = recode_data_for_prediction(
    newdata, is_ordinal,
    category_levels = arguments$category_levels,
    blume_capel_shift = arguments$blume_capel_shift
  )

  if(method == "posterior-mean") {
    # Extract group-specific parameters using projection
    group_params = extract_group_params(object)

    main_group = group_params$main_effects_groups[, group]
    pairwise_group = group_params$pairwise_effects_groups[, group]

    # Reconstruct threshold matrix (variable_type is ifelse(is_ordinal,
    # "ordinal", "blume-capel"), so reconstruct_main's blume-capel branch
    # matches the per-variable parameter count exactly).
    main = reconstruct_main(
      main_group, num_variables, num_categories, variable_type
    )

    # Reconstruct interaction matrix
    pairwise = matrix(0, nrow = num_variables, ncol = num_variables)
    pairwise[lower.tri(pairwise)] = pairwise_group
    pairwise = pairwise + t(pairwise)

    probs = compute_conditional_probs(
      observations = newdata_recoded,
      predict_vars = predict_vars - 1L, # C++ uses 0-based indexing
      pairwise = pairwise,
      main = main,
      num_categories = num_categories,
      variable_type = variable_type,
      baseline_category = baseline_category
    )

    # Add names
    names(probs) = data_columnnames[predict_vars]
    for(v in seq_along(probs)) {
      var_idx = predict_vars[v]
      n_cats = num_categories[var_idx] + 1
      colnames(probs[[v]]) = paste0("cat_", 0:(n_cats - 1))
    }
  }

  if(type == "response") {
    # Return predicted categories (mode)
    pred_matrix = sapply(probs, function(p) {
      apply(p, 1, which.max) - 1L # Convert to 0-based category
    })
    if(is.vector(pred_matrix)) {
      pred_matrix = matrix(pred_matrix, ncol = 1)
    }
    colnames(pred_matrix) = data_columnnames[predict_vars]
    return(pred_matrix)
  }

  return(probs)
}


# ==============================================================================
#   Helper Functions
# ==============================================================================

# Helper function to reconstruct threshold matrix from flat vector
reconstruct_main = function(main_vec, num_variables,
                            num_categories,
                            variable_type) {
  max_cats = max(num_categories)
  main = matrix(NA, nrow = num_variables, ncol = max_cats)

  pos = 1
  for(v in seq_len(num_variables)) {
    if(variable_type[v] != "blume-capel") {
      k = num_categories[v]
      main[v, 1:k] = main_vec[pos:(pos + k - 1)]
      pos = pos + k
    } else {
      main[v, 1:2] = main_vec[pos:(pos + 1)]
      pos = pos + 2
    }
  }

  return(main)
}


# Helper function to recode newdata to the 0-based categories the model was
# fitted on. When the fit carries a recode map (category_levels) per ordinal
# variable, each newdata value is mapped through it exactly as the training data
# was recoded -- so non-contiguous categories and newdata with a different range
# are handled correctly. The map is either:
#   - an unnamed sorted vector of original values (bgm/OMRF): recoded category =
#     position - 1; or
#   - a named vector lookup (bgmCompare): names are original values, values are
#     the final (collapsed) categories, which may be many-to-one.
# Fits without a map (older fits) fall back to the legacy subtract-minimum shift.
recode_data_for_prediction = function(x, is_ordinal,
                                      category_levels = NULL,
                                      blume_capel_shift = NULL) {
  x = as.matrix(x)
  num_variables = ncol(x)

  for(v in seq_len(num_variables)) {
    if(!is_ordinal[v]) {
      # Blume-Capel: shift newdata onto the 0-based scale the model was fit on.
      # Continuous variables carry no shift (NA) and are left untouched.
      if(!is.null(blume_capel_shift) && !is.na(blume_capel_shift[v])) {
        x[, v] = x[, v] - blume_capel_shift[v]
      }
      next
    }

    levels_v = if(!is.null(category_levels)) category_levels[[v]] else NULL

    if(!is.null(levels_v)) {
      if(!is.null(names(levels_v))) {
        # Named lookup (bgmCompare): map original value -> final category.
        recoded = unname(levels_v[match(x[, v], as.numeric(names(levels_v)))])
      } else {
        # Sorted original values (OMRF): recoded category = position - 1.
        recoded = match(x[, v], levels_v) - 1L
      }
      observed = !is.na(x[, v])
      if(any(observed & is.na(recoded))) {
        warning(
          "newdata for variable ", v, " contains category values not ",
          "observed in the training data; predictions for those cells are NA.",
          call. = FALSE
        )
      }
      x[, v] = recoded
    } else {
      # Legacy fallback (fit has no recode map): shift to 0-based by the
      # per-column minimum.
      x[, v] = as.integer(x[, v])
      if(min(x[, v], na.rm = TRUE) > 0) {
        x[, v] = x[, v] - min(x[, v], na.rm = TRUE)
      }
    }
  }

  return(x)
}


# ------------------------------------------------------------------------------
# recode_simulated_to_original()
# ------------------------------------------------------------------------------
# Inverse of recode_data_for_prediction(): map the internal 0-based category
# codes the MRF sampler produces back to the original category values the fit
# was trained on, so simulate() returns data on the same scale predict() expects
# for newdata. Regular ordinal variables carry a recode map in category_levels;
# Blume-Capel variables carry an additive shift in blume_capel_shift, which is
# added back here. Continuous variables have neither and are returned as is.
#
# Two ordinal map forms, matching recode_data_for_prediction():
#   - unnamed sorted original values (bgm/OMRF): code k is the (k + 1)-th value.
#   - named lookup (bgmCompare): names are original values, entries the final
#     (possibly collapsed) codes; invert to the smallest original value mapping
#     to each code, which predict() recodes back to that same code.
recode_simulated_to_original = function(x, category_levels,
                                        blume_capel_shift = NULL) {
  if(is.null(category_levels) && is.null(blume_capel_shift)) {
    return(x)
  }
  for(v in seq_len(ncol(x))) {
    if(!is.null(blume_capel_shift) && !is.na(blume_capel_shift[v])) {
      x[, v] = x[, v] + blume_capel_shift[v]
      next
    }

    levels_v = if(!is.null(category_levels)) category_levels[[v]] else NULL
    if(is.null(levels_v)) next

    if(!is.null(names(levels_v))) {
      inverse = tapply(
        as.numeric(names(levels_v)), as.integer(unname(levels_v)), min
      )
      x[, v] = inverse[as.character(x[, v])]
    } else {
      x[, v] = levels_v[x[, v] + 1L]
    }
  }
  x
}

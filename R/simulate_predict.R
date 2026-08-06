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

  # Validate cores. parallel::detectCores() -- the default -- is documented to
  # return NA when it cannot tell, and check_positive_integer() turns that NA
  # into "missing value where TRUE/FALSE needed" before
  # normalize_parallel_cores() gets the chance to fall back to 1. The check runs
  # on everything the user can actually supply; only that one NA is left to the
  # normalizer.
  if(!(length(cores) == 1L && is.na(cores))) {
    check_positive_integer(cores, "cores")
  }
  cores = normalize_parallel_cores(cores)

  # nsim and iter reach four different simulation paths (OMRF mean, OMRF
  # sample, GGM, mixed); validating them here covers all of them at once.
  check_positive_integer(nsim, "nsim")
  check_positive_integer(iter, "iter")

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
#' for each observation and category. Columns are labelled with the original
#' category values of the training data.
#'
#' For \code{type = "response"}: A matrix with \code{n} rows and
#' \code{length(variables)} columns containing predicted categories, on the
#' original category scale that \code{\link{simulate.bgms}} returns.
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
#' \code{newdata} is matched to the fitted model by position. When it carries
#' column names they must be the model's variables in the model's order, or
#' \code{predict()} stops; when it carries none it is read positionally, with a
#' warning saying so.
#'
#' A discrete cell of \code{newdata} that is \code{NA}, or that holds a category
#' value never observed in the training data, leaves the conditional
#' distribution of every other variable in that row undefined. Those
#' predictions are returned as \code{NA}, with one warning giving the number of
#' rows affected; the variable whose own value is missing is unaffected, since
#' its conditional distribution does not use it.
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

  check_newdata_columns(newdata, data_columnnames)

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
      colnames(probs[[v]]) = probability_column_labels(original_category_values(
        predict_vars[v], num_categories[predict_vars[v]],
        arguments$category_levels, arguments$blume_capel_shift
      ))
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

      labels = probability_column_labels(original_category_values(
        predict_vars[v], num_categories[predict_vars[v]],
        arguments$category_levels, arguments$blume_capel_shift
      ))
      colnames(probs[[v]]) = labels
      colnames(probs_sd[[v]]) = labels
    }

    attr(probs, "sd") = probs_sd
  }

  # Blank the rows the kernel could not condition on, before anything is read
  # off them.
  mask = na_conditioning_mask(is.na(newdata_recoded), predict_vars)
  probs_sd = attr(probs, "sd")
  probs = apply_na_conditioning_mask(probs, mask)
  if(!is.null(probs_sd)) {
    attr(probs, "sd") = apply_na_conditioning_mask(probs_sd, mask)
  }
  warn_na_conditioning(mask)

  if(type == "response") {
    return(format_discrete_response(
      probs, predict_vars, data_columnnames, num_categories,
      arguments$category_levels, arguments$blume_capel_shift
    ))
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
#' for each observation and category. Columns are labelled with the original
#' category values of the training data.
#'
#' For \code{type = "response"}: A matrix with \code{n} rows and
#' \code{length(variables)} columns containing predicted categories, on the
#' original category scale that \code{\link{simulate.bgmCompare}} returns.
#'
#' @details
#' Group-specific parameters are obtained by applying the projection matrix
#' to convert baseline parameters and differences into group-level estimates.
#' The function then computes the conditional distribution of target variables
#' given the observed values of all other variables.
#'
#' The \code{newdata} column-name and missing-value rules of
#' \code{\link{predict.bgms}} apply here unchanged.
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

  check_newdata_columns(newdata, data_columnnames)

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
      colnames(probs[[v]]) = probability_column_labels(original_category_values(
        predict_vars[v], num_categories[predict_vars[v]],
        arguments$category_levels, arguments$blume_capel_shift
      ))
    }
  }

  # Blank the rows the kernel could not condition on, before anything is read
  # off them.
  mask = na_conditioning_mask(is.na(newdata_recoded), predict_vars)
  probs = apply_na_conditioning_mask(probs, mask)
  warn_na_conditioning(mask)

  if(type == "response") {
    return(format_discrete_response(
      probs, predict_vars, data_columnnames, num_categories,
      arguments$category_levels, arguments$blume_capel_shift
    ))
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
          "observed in the training data; those cells are treated as missing.",
          call. = FALSE
        )
      }
      x[, v] = recoded
    } else {
      # No recode map for an ordinal variable. Every fit bgm() and bgmCompare()
      # produce carries one, so this is only reachable for an object built by a
      # bgms old enough to predate it. The shift-by-the-column-minimum fallback
      # that used to stand here read the offset off newdata rather than off the
      # training data, which is the wrong number whenever newdata does not
      # happen to span the training range, and silently so.
      stop(
        "The fitted object carries no category recode map for variable ", v,
        ", so 'newdata' cannot be put on the scale the model was fitted on. ",
        "It predates the recode map; refit with the current bgms.",
        call. = FALSE
      )
    }
  }

  return(x)
}


# ------------------------------------------------------------------------------
# check_newdata_columns()
# ------------------------------------------------------------------------------
# newdata is matched to the fit by position, and a column count is not enough to
# establish that the match is the intended one: a data frame whose columns were
# reordered, or one built from a different subset of the same width, passes the
# count check and then predicts every variable from the wrong neighbours,
# silently. Named newdata therefore has to agree with the fit exactly. Unnamed
# newdata is still accepted -- it carries nothing to check -- but says so.
#
# @param newdata           The (already coerced) newdata matrix.
# @param data_columnnames  The fit's variable names.
# ------------------------------------------------------------------------------
check_newdata_columns = function(newdata, data_columnnames) {
  if(is.null(data_columnnames)) {
    return(invisible(NULL))
  }

  observed = colnames(newdata)
  if(is.null(observed)) {
    warning(
      "'newdata' has no column names, so its columns are matched to the ",
      "fitted model by position: ",
      paste(data_columnnames, collapse = ", "), ".",
      call. = FALSE
    )
    return(invisible(NULL))
  }

  if(identical(observed, data_columnnames)) {
    return(invisible(NULL))
  }

  if(setequal(observed, data_columnnames)) {
    stop(
      "'newdata' holds the fitted model's variables in a different order. ",
      "Columns are matched by position, so reorder 'newdata' to: ",
      paste(data_columnnames, collapse = ", "), ".",
      call. = FALSE
    )
  }

  absent = setdiff(data_columnnames, observed)
  unexpected = setdiff(observed, data_columnnames)
  stop(
    "'newdata' column names do not match the fitted model.",
    if(length(absent) > 0) {
      paste0(" Missing: ", paste(absent, collapse = ", "), ".")
    },
    if(length(unexpected) > 0) {
      paste0(" Not in the model: ", paste(unexpected, collapse = ", "), ".")
    },
    call. = FALSE
  )
}


# ------------------------------------------------------------------------------
# original_category_values()
# ------------------------------------------------------------------------------
# The original category values behind the internal codes 0..num_categories of
# one variable, in code order: the per-variable inverse of
# recode_data_for_prediction(). simulate() returns data on this scale
# (recode_simulated_to_original()), so predict() has to report responses and
# label probability columns on it too, or the round trip the two promise each
# other does not close for 1-based, non-contiguous or Blume-Capel variables.
#
# @param v                  Variable index.
# @param num_categories     Number of categories on top of the baseline, for v.
# @param category_levels    The fit's recode map (list, one entry per variable).
# @param blume_capel_shift  The fit's Blume-Capel shifts (NA elsewhere).
#
# Returns: numeric vector of length num_categories + 1.
# ------------------------------------------------------------------------------
original_category_values = function(v, num_categories,
                                    category_levels = NULL,
                                    blume_capel_shift = NULL) {
  codes = seq.int(0L, num_categories)

  if(!is.null(blume_capel_shift) && !is.na(blume_capel_shift[v])) {
    return(codes + blume_capel_shift[v])
  }

  levels_v = if(!is.null(category_levels)) category_levels[[v]] else NULL
  if(is.null(levels_v)) {
    return(as.numeric(codes))
  }

  if(!is.null(names(levels_v))) {
    # Named lookup (bgmCompare), possibly many-to-one: invert to the smallest
    # original value carrying each code, as recode_simulated_to_original() does,
    # so the value reported recodes back to the code it came from.
    inverse = tapply(
      as.numeric(names(levels_v)), as.integer(unname(levels_v)), min
    )
    return(unname(inverse[as.character(codes)]))
  }

  as.numeric(levels_v[codes + 1L])
}


# ------------------------------------------------------------------------------
# probability_column_labels()
# ------------------------------------------------------------------------------
# Column labels for a discrete variable's probability matrix, naming the
# ORIGINAL category values rather than the internal codes.
# ------------------------------------------------------------------------------
probability_column_labels = function(values) {
  paste0("cat_", values)
}


# ------------------------------------------------------------------------------
# na_conditioning_mask()
# ------------------------------------------------------------------------------
# Which rows leave a target variable's conditional distribution undefined,
# because some OTHER variable it conditions on is missing.
#
# The C++ kernels take the recoded observation matrix at face value: an NA cell
# reaches them as R's NA_integer_ sentinel and enters the rest score of every
# other variable in that row as a huge negative number, which comes back as a
# confident one-hot distribution rather than as missingness. Both origins of the
# NA -- a plain NA in newdata and a category value never seen in training -- are
# caught here, so the masking is done once in R after the kernel returns.
#
# @param is_missing    n x p logical matrix, TRUE where the conditioning value
#   is missing (in the kernel's own column order).
# @param target_cols   Column of `is_missing` each prediction targets.
#
# Returns: list of logical vectors, one per target, TRUE for rows to blank.
# ------------------------------------------------------------------------------
na_conditioning_mask = function(is_missing, target_cols) {
  total = rowSums(is_missing)
  lapply(target_cols, function(v) total - is_missing[, v] > 0)
}


# ------------------------------------------------------------------------------
# apply_na_conditioning_mask()
# ------------------------------------------------------------------------------
# Blank the masked rows of each prediction matrix and say how many rows of
# newdata lost a prediction. Warning is the caller's, once per call.
# ------------------------------------------------------------------------------
apply_na_conditioning_mask = function(predictions, mask) {
  for(k in seq_along(predictions)) {
    rows = mask[[k]]
    if(any(rows)) {
      predictions[[k]][rows, ] = NA_real_
    }
  }
  predictions
}


warn_na_conditioning = function(mask) {
  affected = sum(Reduce(`|`, mask))
  if(affected > 0) {
    warning(
      "newdata has ", affected, " row(s) in which a conditioning variable is ",
      "missing or carries a category value not observed in the training data. ",
      "The conditional distribution of the other variables is undefined in ",
      "those rows, so their predictions are NA.",
      call. = FALSE
    )
  }
  invisible(affected)
}


# ------------------------------------------------------------------------------
# row_modes()
# ------------------------------------------------------------------------------
# Index of the largest entry in each row (ties to the first), NA for rows that
# carry any NA -- where which.max() would return integer(0) and collapse the
# result to a list.
# ------------------------------------------------------------------------------
row_modes = function(probabilities) {
  out = rep(NA_integer_, nrow(probabilities))
  usable = !rowSums(is.na(probabilities)) > 0
  if(any(usable)) {
    out[usable] = max.col(
      probabilities[usable, , drop = FALSE],
      ties.method = "first"
    )
  }
  out
}


# ------------------------------------------------------------------------------
# format_discrete_response()
# ------------------------------------------------------------------------------
# Point predictions for discrete variables: the mode of each conditional
# distribution, reported on the ORIGINAL category scale.
#
# The n x L result is preallocated rather than built with sapply(), which
# collapses to a vector when n is 1 and then mis-shapes into an n x 1 matrix
# that the column names no longer fit.
# ------------------------------------------------------------------------------
format_discrete_response = function(probabilities, predict_vars,
                                    data_columnnames, num_categories,
                                    category_levels = NULL,
                                    blume_capel_shift = NULL) {
  out = matrix(
    NA_real_,
    nrow = nrow(probabilities[[1]]), ncol = length(predict_vars)
  )
  colnames(out) = data_columnnames[predict_vars]

  for(v in seq_along(predict_vars)) {
    var_idx = predict_vars[v]
    values = original_category_values(
      var_idx, num_categories[var_idx], category_levels, blume_capel_shift
    )
    out[, v] = values[row_modes(probabilities[[v]])]
  }

  out
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

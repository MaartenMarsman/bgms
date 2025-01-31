#' Bayesian Variable Selection or Bayesian Estimation for Differences in Markov Random Fields
#'
#' @description
#' The `bgmCompare` function estimates the pseudoposterior distribution of
#' the parameters of a Markov Random Field (MRF) model for mixed binary and ordinal
#' variables, as well as differences in pairwise interactions and category thresholds
#' across groups. Groups are assumed to be `G` independent samples.
#'
#' @details
#' ### Pairwise Interactions
#' Pairwise interactions between variables `i` and `j` are modeled as:
#' \deqn{\boldsymbol{\theta}_{ij} = \phi_{ij} + \boldsymbol{\delta}_{ij},}{
#' \boldsymbol{\theta}_{ij} = \phi_{ij} + \boldsymbol{\delta}_{ij},}
#' where:
#' - \eqn{\boldsymbol{\theta}_{ij}}{\boldsymbol{\theta}_{ij}} is the vector of pairwise interaction parameters of length `G`.
#' - \eqn{\phi_{ij}}{\phi_{ij}} is the overall pairwise interaction (nuisance parameter).
#' - \eqn{\boldsymbol{\delta}_{ij}}{\boldsymbol{\delta}_{ij}} represents group-specific differences constrained to sum to zero for identification.
#'
#' ### Ordinal Variables
#' The function supports two types of ordinal variables:
#' 1. **Regular ordinal variables**: Introduce a threshold parameter for each category except the lowest, modeled as:
#'    \deqn{\boldsymbol{\mu}_{ic} = \tau_{ic} + \boldsymbol{\epsilon}_{ic},}{\boldsymbol{\mu}_{ic} = \tau_{ic} + \boldsymbol{\epsilon}_{ic},}
#'    where:
#'    - \eqn{\tau_{ic}}{\tau_{ic}} denotes an overall effect (nuisance parameter).
#'    - \eqn{\boldsymbol{\epsilon}_{ic}}{\boldsymbol{\epsilon}_{ic}} represents group-specific differences constrained to sum to zero.
#'
#' 2. **Blume-Capel ordinal variables**: Assume a specific reference category and score responses based on distance to it:
#'    \deqn{\boldsymbol{\mu}_{ic} = (\tau_{i1} + \boldsymbol{\epsilon}_{i1}) \cdot c + (\tau_{i2} + \boldsymbol{\epsilon}_{i2}) \cdot (c - r)^2,}{
#'    \boldsymbol{\mu}_{ic} = (\tau_{i1} + \boldsymbol{\epsilon}_{i1}) \cdot c + (\tau_{i2} + \boldsymbol{\epsilon}_{i2}) \cdot (c - r)^2,}
#'    where:
#'    - `r` is the reference category.
#'    - \eqn{\tau_{i1}}{\tau_{i1}} and \eqn{\tau_{i2}}{\tau_{i2}} are nuisance parameters.
#'    - \eqn{\boldsymbol{\epsilon}_{i1}}{\boldsymbol{\epsilon}_{i1}} and \eqn{\boldsymbol{\epsilon}_{i2}}{\boldsymbol{\epsilon}_{i2}} represent group-specific differences.
#'
#' ### Variable Selection
#' Bayesian variable selection enables testing of parameter differences or equivalence across groups. Independent spike-and-slab priors are applied to difference parameters:
#' - **Bernoulli Model**: Assigns a fixed probability to parameter inclusion.
#' - **Beta-Bernoulli Model**: Incorporates a beta prior to model inclusion probabilities.
#'
#' ### Saving Options
#' Users can store sampled states for parameters (`main_effects`, `pairwise_effects`, `indicator`) during Gibbs sampling. Enabling these options (`save_main`, `save_pairwise`, `save_indicator`) increases output size and memory usage, so use them judiciously.
#'
#' @param x Data frame or matrix with binary and ordinal responses. Regular ordinal variables should be coded as integers starting from 0. Missing categories are collapsed for regular ordinal variables but retained for Blume-Capel variables.
#' @param y A data frame or matrix similar to `x`, used for two-group designs. `x` contains Group 1 data, and `y` contains Group 2 data. Ignored for multi-group designs.
#' @param g Group membership vector for rows in `x`. Required for multi-group designs; ignored if `y` is provided.
#' @param difference_selection Logical. Enables modeling of inclusion/exclusion of parameter differences (`TRUE`) or estimation of all differences (`FALSE`). Default: `TRUE`.
#' @param save_main,save_pairwise,save_indicator Logical. Enable saving sampled states for `main_effects`, `pairwise_effects`, and `indicator`, respectively. Default: `FALSE`.
#' @param main_difference_model Character. Specifies how to handle threshold differences when categories are unmatched. Options: `"Collapse"`, `"Free"`. The option "Collapse" collapses categories unobserved in one or more groups. The option "Free" option estimates thresholds separately for each group and does not model their difference. Default: `"Free"`.
#' @param variable_type Character or vector. Specifies the type of variables in `x` (`"ordinal"` or `"blume-capel"`). Default: `"ordinal"`.
#' @param reference_category Integer or vector. Reference category for Blume-Capel variables. Required if there is at least one Blume-Capel variable.
#' @param pairwise_difference_scale Double. Scale parameter for the Cauchy prior on pairwise differences. Default: `1`.
#' @param main_difference_scale Double. Scale parameter for the Cauchy prior on threshold differences. Default: `1`.
#' @param pairwise_difference_prior,main_difference_prior Character. Specifies the inclusion probability model (`"Bernoulli"` or `"Beta-Bernoulli"`). Default: `"Bernoulli"`.
#' @param iter,burnin Integer. Number of Gibbs iterations (`iter`) and burn-in iterations (`burnin`). Defaults: `iter = 1e4`, `burnin = 5e2`.
#' @param na.action Character. Specifies handling of missing data. `"listwise"` deletes rows with missing values; `"impute"` imputes values during Gibbs sampling. Default: `"listwise"`.
#' @param display_progress Logical. Show progress bar during computation. Default: `TRUE`.
#' @param threshold_alpha,threshold_beta Double. Shape parameters for the beta-prime prior on nuisance threshold parameters.
#' @param interaction_scale Double. Scale of the Cauchy prior for nuisance pairwise interactions. Default: `2.5`.
#' @param main_beta_bernoulli_alpha,main_beta_bernoulli_beta Double. Shape parameters for the Beta-Bernoulli prior on threshold differences.
#' @param pairwise_beta_bernoulli_alpha,pairwise_beta_bernoulli_beta Double. Shape parameters for the Beta-Bernoulli prior on pairwise differences.
#'
#' @return A list containing the posterior means and, optionally, sampled states based on the `save_*` options. The returned components include:
#' - `posterior_mean_main`, `posterior_mean_pairwise`, and `posterior_mean_indicator` for posterior means.
#' - If saving options are enabled:
#'   - `raw_samples_main` for sampled states of `main_effects`.
#'   - `raw_samples_pairwise` for sampled states of `pairwise_effects`.
#'   - `raw_samples_indicator` for sampled states of the inclusion indicators.
#'
#' In addition to the results of the analysis, the output lists some of the
#' arguments of its call. This is useful for post-processing the results.
#'
#' @importFrom utils packageVersion
#'
#' @export
bgmCompare_dev = function(x,
                          y,
                          g,
                          t.test = TRUE, #this is for checking software
                          difference_selection = TRUE,
                          main_difference_model = c("Free", "Collapse", "Constrain"),
                          variable_type = "ordinal",
                          reference_category,
                          pairwise_difference_scale = 1,
                          main_difference_scale = 1,
                          pairwise_difference_prior = c("Bernoulli", "Beta-Bernoulli"),
                          main_difference_prior = c("Bernoulli", "Beta-Bernoulli"),
                          pairwise_difference_probability = 0.5,
                          main_difference_probability = 0.5,
                          pairwise_beta_bernoulli_alpha = 1,
                          pairwise_beta_bernoulli_beta = 1,
                          main_beta_bernoulli_alpha = 1,
                          main_beta_bernoulli_beta = 1,
                          interaction_scale = 2.5,
                          threshold_alpha = 0.5,
                          threshold_beta = 0.5,
                          iter = 1e4,
                          burnin = 5e2,
                          na.action = c("listwise", "impute"),
                          save = FALSE,
                          save_main = FALSE,
                          save_pairwise = FALSE,
                          save_indicator = FALSE,
                          display_progress = TRUE) {

  # Validate inputs
  t.test <- hasArg(y) || t.test

  # Deprecation warning for save parameter
  if(hasArg(save)) {
    warning("`save` is deprecated. Use `save_main`, `save_pairwise`, or `save_indicator` instead.")
    save_main <- save_main || save
    save_pairwise <- save_pairwise || save
    save_indicator <- save_indicator || save
  }

  # Check and preprocess data
  x <- data_check(x, "x")
  if (hasArg(y)) {
    y <- data_check(y, "y")
    if (ncol(x) != ncol(y)) stop("x and y must have the same number of columns.")
  }

  if(!t.test & !hasArg(g))
    stop(paste0("For multi-group designs, the bgmCompare function requires input for\n",
                "either y (group 2 data) or g (group indicator)."))

  # Validate group indicators
  if (!t.test && hasArg(g)) {
    g = as.vector(g)
    if (anyNA(g)) stop("g cannot contain missing values.")
    if (length(g) != nrow(x)) stop("Length of g must match number of rows in x.")

    unique_g = unique(g)
    if(length(unique_g) == 2 && t.test == TRUE) {
      y = x[g == unique_g[2],]
      x = x[g == unique_g[1],]
      t.test = TRUE
    }
  }

  # Model and preprocessing
  if(t.test == FALSE) #True if either hasArg(y) or length(unique(g)) = 2
    y = NULL
  if(!hasArg(g))
    g = NULL

  model = check_compare_model_dev(
    x = x, y = y, g = g, ttest = t.test,
    difference_selection = difference_selection,variable_type = variable_type,
    reference_category = reference_category,
    pairwise_difference_scale = pairwise_difference_scale,
    main_difference_scale = main_difference_scale,
    pairwise_difference_prior = pairwise_difference_prior,
    main_difference_prior = main_difference_prior,
    pairwise_difference_probability = pairwise_difference_probability,
    main_difference_probability = main_difference_probability,
    main_beta_bernoulli_alpha = main_beta_bernoulli_alpha,
    main_beta_bernoulli_beta = main_beta_bernoulli_beta,
    pairwise_beta_bernoulli_alpha = pairwise_beta_bernoulli_alpha,
    pairwise_beta_bernoulli_beta = pairwise_beta_bernoulli_beta,
    interaction_scale = interaction_scale, threshold_alpha = threshold_alpha,
    threshold_beta = threshold_beta,
    main_difference_model = main_difference_model)

  ordinal_variable <- model$variable_bool
  reference_category <- model$reference_category
  independent_thresholds <- (model$main_difference_model == "Free")

  # Check Gibbs input
  check_positive_integer(iter, "iter")
  check_non_negative_integer(burnin, "burnin")

  # Check na.action
  na.action_input <- na.action
  na.action <- try(match.arg(na.action), silent = TRUE)
  if (inherits(na.action, "try-error")) {
    stop(sprintf("Invalid value for `na.action`. Expected 'listwise' or 'impute', got: %s", na.action_input))
  }

  # Check save options
  save_main <- check_logical(save_main, "save_main")
  save_pairwise <- check_logical(save_pairwise, "save_pairwise")
  save_indicator <- check_logical(save_indicator, "save_indicator")

  # Check display_progress
  display_progress <- check_logical(display_progress, "display_progress")

  ## Format data
  data <- compare_reformat_data_dev(
    x = x, y = y, g = g, ttest = t.test,
    na.action = na.action,
    variable_bool = ordinal_variable,
    reference_category = reference_category,
    main_difference_model = model$main_difference_model
  )

  x = data$x


  y = data$y

  no_categories_gr1 = data$no_categories[, 1]
  no_categories_gr2 = data$no_categories[, 2]
  missing_index_gr1 = data$missing_index_gr1
  missing_index_gr2 = data$missing_index_gr2

  no_obs_groups = c(nrow(x), nrow(y))

  na_impute = data$na_impute
  reference_category = data$reference_category
  no_variables = ncol(x)
  no_interactions = no_variables * (no_variables - 1) / 2

  # Compute `n_cat_obs`
  n_cat_obs_gr1 <- compute_n_cat_obs_dev(x, no_categories_gr1)
  n_cat_obs_gr2 <- compute_n_cat_obs_dev(y, no_categories_gr2)

  # Compute sufficient statistics for Blume-Capel variables
  sufficient_blume_capel_gr1 <- compute_sufficient_blume_capel_dev(x, reference_category,ordinal_variable)
  sufficient_blume_capel_gr2 <- compute_sufficient_blume_capel_dev(y, reference_category,ordinal_variable)

  # Index vector used to sample interactions in a random order -----------------
  Index = matrix(0, nrow = no_interactions, ncol = 3)
  counter = 0
  for(variable1 in 1:(no_variables - 1)) {
    for(variable2 in (variable1 + 1):no_variables) {
      counter =  counter + 1
      Index[counter, ] = c(counter, variable1 - 1, variable2 - 1)
    }
  }

  # Gibbs sampling
    #the old rcpp function (for two groups only)
    out <- compare_ttest_gibbs_sampler(
      observations_gr1 = data$x,
      observations_gr2 = data$y,
      no_categories_gr1 = data$no_categories[, 1],
      no_categories_gr2 = data$no_categories[, 2],
      interaction_scale = interaction_scale,
      pairwise_difference_scale = pairwise_difference_scale,
      main_difference_scale = main_difference_scale,
      pairwise_difference_prior = model$pairwise_difference_prior,
      main_difference_prior = model$main_difference_prior,
      inclusion_probability_difference = model$inclusion_probability_difference,
      pairwise_beta_bernoulli_alpha = pairwise_beta_bernoulli_alpha,
      pairwise_beta_bernoulli_beta = pairwise_beta_bernoulli_beta,
      main_beta_bernoulli_alpha = main_beta_bernoulli_alpha,
      main_beta_bernoulli_beta = main_beta_bernoulli_beta,
      Index = Index, iter = iter, burnin = burnin,
      n_cat_obs_gr1 = n_cat_obs_gr1, n_cat_obs_gr2 = n_cat_obs_gr2,
      sufficient_blume_capel_gr1 = sufficient_blume_capel_gr1,
      sufficient_blume_capel_gr2 = sufficient_blume_capel_gr2,
      threshold_alpha = threshold_alpha, threshold_beta = threshold_beta,
      na_impute = na_impute, missing_index_gr1 = missing_index_gr1,
      missing_index_gr2 = missing_index_gr2, ordinal_variable = ordinal_variable,
      reference_category = reference_category,
      independent_thresholds = independent_thresholds,
      save = save,
      display_progress = display_progress,
      difference_selection = difference_selection
    )

    # Main output handler in the wrapper function
    output <- prepare_output_bgmCompare_old(
      out = out, x = x, t.test = t.test,
      independent_thresholds = independent_thresholds, no_variables = no_variables,
      no_categories = if (t.test) cbind(no_categories_gr1, no_categories_gr2) else no_categories,
      group = group, iter = iter,
      data_columnnames = if (is.null(colnames(x))) paste0("Variable ", seq_len(ncol(x))) else colnames(x),
      save_options = list(save_main = save_main, save_pairwise = save_pairwise,
                          save_indicator = save_indicator),
      difference_selection = difference_selection, na_action = na.action,
      na_impute = na_impute, variable_type = variable_type, burnin = burnin,
      interaction_scale = interaction_scale, threshold_alpha = threshold_alpha,
      threshold_beta = threshold_beta,
      main_difference_model = model$main_difference_model,
      pairwise_difference_prior = model$pairwise_difference_prior,
      main_difference_prior = model$main_difference_prior,
      inclusion_probability_difference = model$inclusion_probability_difference,
      pairwise_beta_bernoulli_alpha = pairwise_beta_bernoulli_alpha,
      pairwise_beta_bernoulli_beta = pairwise_beta_bernoulli_beta,
      main_beta_bernoulli_alpha = main_beta_bernoulli_alpha,
      main_beta_bernoulli_beta = main_beta_bernoulli_beta,
      main_difference_scale = main_difference_scale,
      pairwise_difference_scale = pairwise_difference_scale,
      projection,
      is_ordinal_variable = ordinal_variable
    )



  return(output)
}
#' Compare the parameters of a graphical model for binary and ordinal variables
#' between two groups using Bayesian edge selection.
#'
#' @description
#' The \code{bgmCompare} function estimates the pseudoposterior distribution of
#' the parameters of a Markov Random Field model for mixed binary and ordinal
#' variables, and the differences in pairwise interactions and category thresholds
#' between two groups. The groups can be two independent samples, or paired.
#'
#' @details
#' In the first group, the pairwise interactions between the variables \eqn{i}{i}
#' and \eqn{j}{j} are modeled as
#' \deqn{\sigma_{\text{ij}} = \theta_{\text{ij}} + \delta_{\text{ij}} / 2,}{\sigma_{\text{ij}} = \theta_{\text{ij}} + \delta_{\text{ij}} / 2,}
#' and in the second group as
#' \deqn{\sigma_{\text{ij}} = \theta_{\text{ij}} - \delta_{\text{ij}} / 2,}{\sigma_{\text{ij}} = \theta_{\text{ij}} - \delta_{\text{ij}} / 2.}
#' The pairwise interaction parameter \eqn{\theta_{\text{ij}}}{\theta_{\text{ij}}}
#' denotes an overall effect that is considered nuisance, and attention is focused
#' on the pairwise difference parameter \eqn{\delta_{\text{ij}}}{\delta_{\text{ij}}},
#' which reflects the difference in the pairwise interaction between the two groups.
#'
#' The \code{bgmCompare} supports two types of ordinal variables, which can be mixed.
#' The default ordinal variable introduces a threshold parameter for each category
#' except the lowest category. For this variable type, the threshold parameter for
#' variable \eqn{i}{i}, category \eqn{c}{c}, is modeled as
#' \deqn{\mu_{\text{ic}} = \tau_{\text{ic}} + \epsilon_{\text{ic}} / 2,}{\mu_{\text{ic}} = \tau_{\text{ic}} + \epsilon_{\text{ic}} / 2,}
#' in the first group and in the second group as
#' \deqn{\mu_{\text{ic}} = \tau_{\text{ic}} - \epsilon_{\text{ic}} / 2,}{\mu_{\text{ic}} = \tau_{\text{ic}} - \epsilon_{\text{ic}} / 2.}
#' The category threshold parameter \eqn{\tau_{\text{ic}}}{\tau_{\text{ic}}} denotes
#' an overall effect that is considered nuisance, and attention is focused on the
#' threshold difference parameter \eqn{\epsilon_{\text{ic}}}{\epsilon_{\text{ic}}},
#' which reflects the difference in threshold of for variable \eqn{i}{i}, category
#' \eqn{c}{c} between the two groups.
#'
#' The Blume-Capel ordinal variable assumes that there is a specific reference
#' category, such as ``neutral'' in a Likert scale, and responses are scored
#' according to their distance from this reference category. In the first group,
#' the threshold parameters are modelled as
#' \deqn{\mu_{\text{ic}} = (\tau_{\text{i1}} + \epsilon_{\text{i1}} / 2) \times \text{c} + (\tau_{\text{i2}} + \epsilon_{\text{i2}} / 2) \times (\text{c} - \text{r})^2,}{ {\mu_{\text{ic}} = (\tau_{\text{i1}} + \epsilon_{\text{i1}} / 2) \times \text{c} + (\tau_{\text{i2}} + \epsilon_{\text{i2}} / 2) \times (\text{c} - \text{r})^2,}}
#' and in the second groups as
#' \deqn{\mu_{\text{ic}} = (\tau_{\text{i1}} - \epsilon_{\text{i1}} / 2) \times \text{c} + (\tau_{\text{i2}} - \epsilon_{\text{i2}} / 2) \times (\text{c} - \text{r})^2.}{ {\mu_{\text{ic}} = (\tau_{\text{i1}} - \epsilon_{\text{i1}} / 2) \times \text{c} + (\tau_{\text{i2}} - \epsilon_{\text{i2}} / 2) \times (\text{c} - \text{r})^2.}}
#' The linear and quadratic category threshold parameters
#' \eqn{\tau_{\text{i1}}}{\tau_{\text{i1}}} and \eqn{\tau_{\text{i2}}}{\tau_{\text{i2}}}
#' denote overall effects that are considered nuisance, and attention is focused
#' on the two threshold difference parameters
#' \eqn{\epsilon_{\text{i1}}}{\epsilon_{\text{i1}}} and
#' \eqn{\epsilon_{\text{i2}}}{\epsilon_{\text{i2}}}, which reflect the differences
#' in the quadratic model for the variable \eqn{i}{i} between the two groups.
#'
#' In a paired samples design, the pairwise interactions between the samples must
#' be modeled to account for the dependence in the repeated measures. This
#' dependence can also be modeled by assigning a random effect per case to each
#' variable in the model. The matrix of between-sample pairwise interactions can
#' then be viewed as the covariance matrix of the random effects.
#'
#' Bayesian variable selection is used to model the presence or absence of the
#' difference parameters \eqn{\delta}{\delta} and \eqn{\epsilon}{\epsilon}, which
#' allow us to assess parameter differences between the two groups. Independent
#' spike and slab priors are specified for these difference parameters. The spike
#' and slab priors use binary indicator variables to select the difference parameters,
#' assigning them a diffuse Cauchy prior with an optional scaling parameter if
#' selected, or setting the difference parameter to zero if not selected.
#'
#' The function offers two models for the probabilistic inclusion of parameter
#' differences:
#' \itemize{
#'   \item \strong{Bernoulli Model}: This model assigns a fixed probability of
#'   selecting a parameter difference, treating them as independent events. A
#'   probability of 0.5 indicates no preference, giving equal prior weight to
#'   all configurations.
#'   \item \strong{Beta-Bernoulli Model}: Introduces a beta distribution prior
#'   for the inclusion probability that models the complexity of the configuration
#'   of the difference indicators. When the alpha and beta shape parameters of the
#'   beta distribution are set to 1, the model assigns the same prior weight to
#'   the number of differences present (i.e., a configuration with two differences
#'   or with four differences is a priori equally likely).
#'   }
#'   Inclusion probabilities can be specified for pairwise interactions with
#'   \code{pairwise_difference_probability} and for category thresholds with
#'   \code{threshold_difference_probability}.
#'
#' The pairwise interaction parameters \eqn{\theta}{\theta}, the category
#' threshold parameters \eqn{\tau}{\tau}, and, in paired-samples designs,
#' the between-sample interactions \eqn{\omega}{\omega} are considered
#' nuisance parameters that are common to all models. The pairwise interaction
#' parameters \eqn{\theta}{\theta} and the between-sample interactions
#' \eqn{\omega}{\omega} are assigned a diffuse Cauchy prior with an optional
#' scaling parameter. The exponent of the category threshold parameters
#' \eqn{\tau}{\tau} are assigned beta-prime distribution with optional scale
#' values.
#'
#' @param x A data frame or matrix with \eqn{n_1}{n_1} rows and \code{p} columns
#' containing binary and ordinal responses for the first group. Regular ordinal
#' variables are recoded as non-negative integers \code{(0, 1, ..., m)} if not
#' already done. Unobserved categories are collapsed into other categories after
#' recoding (i.e., if category 1 is unobserved, the data are recoded from (0, 2)
#' to (0, 1)). Blume-Capel ordinal variables are also coded as non-negative
#' integers if not already done. However, since ``distance'' from the reference
#' category plays an important role in this model, unobserved categories are not
#' collapsed after recoding.
#' @param y A data frame or matrix with \eqn{n_2}{n_2} rows and \code{p} columns
#' containing binary and ordinal responses for the second group. The variables
#' or columns in \code{y} must match the variables or columns in \code{x}. In
#' the paired samples design, the rows in \code{x} must match the rows in
#' \code{y}. Note that \code{y} and \code{y} are recoded independently, although
#' the function checks that the number of different responses observed matches
#' between \code{y} and \code{x}.
#' @param paired Logical, if \code{TRUE} models the case-specific dependence using
#' a paired-samples design; if \code{FALSE} treats the groups as independent.
#' Default is \code{FALSE}.
#' @param variable_type A string or vector specifying the type of variables
#' in \code{x} (and \code{y}). Supported types are "ordinal" and "blume-capel",
#' with binary variables treated as "ordinal". Default is "ordinal".
#' @param reference_category The reference category in the Blume-Capel model.
#' Should be an integer within the range of integer values observed for the "blume-capel"
#' variable. Can be a single number that sets the reference category for all
#' Blume-Capel variables at once, or a vector of length \code{p}, where the
#' \code{i}-th element is the reference category for the \code{i} variable if
#' it is a Blume-Capel variable, and elements for other variable types are
#' ignored. The value of the reference category is also recoded when
#' `bgmCompare` recodes the corresponding observations. Required only if there
#' is at least one variable of type "blume-capel".
#' @param pairwise_difference_scale The scale of the Cauchy distribution that is
#' used as the prior for the pairwise difference parameters. Defaults to
#' \code{0.1}.
#' @param threshold_difference_scale The scale of the Cauchy distribution that
#' is used as the prior for the threshold difference parameters. Defaults to
#' \code{0.1}.
#' @param pairwise_difference_prior A character string that specifies the model
#' to use for the  inclusion probability of pairwise differences. Options are
#' "Bernoulli" or "Beta-Bernoulli". Default is "Bernoulli".
#' @param threshold_difference_prior A character string that specifies the model
#' to use for the  inclusion probability of threshold differences. Options are
#' "Bernoulli" or "Beta-Bernoulli". Default is "Bernoulli".
#' @param pairwise_difference_probability The inclusion probability for a
#' pairwise difference in the Bernoulli model. Can be a single probability or a
#' matrix of \code{p} rows and \code{p} columns specifying the probability of a
#' difference for each edge pair. Defaults to 0.5.
#' @param threshold_difference_probability The inclusion probability for a
#' threshold difference in the Bernoulli model. Defaults to 0.5, implying no
#' prior preference. Can be a single probability or a matrix of \code{p} rows and
#' \code{max(m)} columns specifying the probability of a difference for each
#' category threshold. Defaults to 0.5.
#' @param main_beta_bernoulli_alpha The alpha parameter of the beta distribution
#' for the Beta-Bernoulli model for group differences in category thresholds.
#' Default is 1.
#' @param main_beta_bernoulli_beta The beta parameter of the beta distribution
#' for the Beta-Bernoulli model for group differences in category thresholds.
#' Default is 1.
#' @param pairwise_beta_bernoulli_alpha The alpha parameter of the beta distribution
#' for the Beta-Bernoulli model for group differences in pairwise interactions.
#' Default is 1.
#' @param pairwise_beta_bernoulli_beta The beta parameter of the beta distribution
#' for the Beta-Bernoulli model for group differences in pairwise interactions.
#' Default is 1.
#' @param interaction_scale The scale of the Cauchy distribution that is used as
#' a prior for the nuisance pairwise interaction parameters. Defaults to
#' \code{2.5}.
#' @param cross_lagged_scale The scale of the Cauchy distribution that is used as
#' a prior for the nuisance cross-laggedinteraction parameters in paired samples
#' designs. Defaults to \code{2.5}.
#' @param threshold_alpha,threshold_beta The shape parameters of the beta-prime
#' prior density for the nuisance threshold parameters. Must be positive values.
#' If the two values are equal, the prior density is symmetric about zero. If
#' \code{threshold_beta} is greater than \code{threshold_alpha}, the distribution
#' is left-skewed, and if \code{threshold_beta} is less than \code{threshold_alpha},
#' it is right-skewed. Smaller values tend to result in more diffuse prior
#' distributions.
#' @param iter The function uses a Gibbs sampler to sample from the posterior
#' distribution of the model parameters and indicator variables. How many
#' iterations should this Gibbs sampler run? The default of \code{1e4} is for
#' illustrative purposes. For stable estimates, it is recommended to run the
#' Gibbs sampler for at least \code{1e5} iterations.
#' @param burnin The number of iterations of the Gibbs sampler before saving its
#' output. Since it may take some time for the Gibbs sampler to converge to the
#' posterior distribution, it is recommended not to set this number too low.
#' @param na.action How do you want the function to handle missing data? If
#' \code{na.action = "listwise"}, listwise deletion is used. If
#' \code{na.action = "impute"}, missing data will be imputed iteratively during
#' Gibbs sampling. Since imputation of missing data can have a negative impact
#' on the convergence speed of the Gibbs sampling procedure, it is recommended
#' to run the procedure for more iterations.
#' @param save Should the function collect and return all samples from the Gibbs
#' sampler (\code{save = TRUE})? Or should it only return the (model-averaged)
#' posterior means (\code{save = FALSE})? Defaults to \code{FALSE}.
#' @param display_progress Should the function show a progress bar
#' (\code{display_progress = TRUE})? Or not (\code{display_progress = FALSE})?
#' The default is \code{TRUE}.
#'
#' @return If \code{save = FALSE} (the default), the result is a list of class
#' ``bgmCompare'' containing the following matrices:
#' \itemize{
#'    \item \code{indicator}: A matrix with \code{p} rows and \code{p}
#'    columns containing the posterior inclusion probabilities of the differences
#'    in pairwise interactions on the off-diagonal and the posterior inclusion
#'    probabilities of the differences in category thresholds on the diagonal.
#'    \item \code{difference_pairwise}: A matrix with \code{p} rows and \code{p}
#'    columns, containing model-averaged posterior means of the differences in
#'    pairwise interactions.
#'    \item \code{difference_threshold}: A matrix with \code{p} rows and
#'    \code{max(m)} columns, containing model-averaged posterior means of the
#'    differences in category thresholds.
#'    \item \code{interactions}: A matrix with \code{p} rows and \code{p} columns,
#'    containing posterior means of the nuisance pairwise interactions.
#'    \item \code{thresholds}: A matrix with \code{p} rows and \code{max(m)}
#'    columns containing the posterior means of the nuisance category thresholds.
#'    In the case of ``blume-capel'' variables, the first entry is the parameter
#'    for the linear effect and the second entry is the parameter for the
#'    quadratic effect, which models the offset to the reference category.
#'  }
#' If \code{paired = TRUE}, the list will also contain
#' \itemize{
#'    \item \code{cross_lagged}: A matrix with \code{p} rows and
#'    \code{p} columns containing the posterior means of the estimated
#'    cross-lagged relations.
#' }
#'
#' If \code{save = TRUE}, the result is a list of class ``bgmCompare''
#' containing the following matrices:
#'  \itemize{
#'    \item \code{indicator_pairwise}: A matrix with \code{iter} rows and
#'    \code{p * (p - 1) / 2} columns containing the inclusion indicators for the
#'    differences in pairwise interactions from each iteration of the Gibbs
#'    sampler.
#'    \item \code{difference_pairwise}: A matrix with \code{iter} rows and
#'    \code{p * (p - 1) / 2} columns, containing parameter states for the
#'    differences in pairwise interactions from each iteration of the Gibbs
#'    sampler.
#'    \item \code{indicator_threshold}: A matrix with \code{iter} rows and
#'    \code{sum(m)} columns, containing the inclusion indicators for the
#'    differences in category thresholds from each iteration of the Gibbs
#'    sampler.
#'    \item \code{difference_threshold}: A matrix with \code{iter} rows and
#'    \code{sum(m)} columns, containing the parameter states for the differences
#'    in category thresholds from each iteration of the Gibbs sampler.
#'    \item \code{interactions}: A matrix with \code{iter} rows and
#'    \code{p * (p - 1) / 2} columns, containing parameter states for the
#'    nuisance pairwise interactions in each iteration of the Gibbs sampler.
#'    \item \code{thresholds}: A matrix with \code{iter} rows and \code{sum(m)}
#'    columns, containing parameter states for the nuisance category thresholds
#'    in each iteration of the Gibbs sampler.
#'  }
#' If \code{paired = TRUE}, the list will also contain
#' \itemize{
#' \item \code{random_effect_covariance}: A matrix with \code{iter} rows and
#' \code{p * (p - 1) / 2} columns, containing parameter states for the
#' cross-lagged relations.
#' }
#' Column averages of these matrices provide the model-averaged posterior means.
#'
#' In addition to the results of the analysis, the output lists some of the
#' arguments of its call. This is useful for post-processing the results.
#'
#' @importFrom utils packageVersion
#' @export
bgmCompare = function(x,
                      y,
                      paired = FALSE,
                      variable_type = "ordinal",
                      reference_category,
                      pairwise_difference_scale = 0.1,
                      main_difference_scale = 0.1,
                      pairwise_difference_prior = c("Bernoulli", "Beta-Bernoulli"),
                      main_difference_prior = c("Bernoulli", "Beta-Bernoulli"),
                      pairwise_difference_probability = 0.5,
                      main_difference_probability = 0.5,
                      pairwise_beta_bernoulli_alpha = 1,
                      pairwise_beta_bernoulli_beta = 1,
                      main_beta_bernoulli_alpha = 1,
                      main_beta_bernoulli_beta = 1,
                      interaction_scale = 2.5,
                      cross_lagged_scale = 2.5,
                      threshold_alpha = 0.5,
                      threshold_beta = 0.5,
                      iter = 1e4,
                      burnin = 1e3,
                      na.action = c("listwise", "impute"),
                      save = FALSE,
                      display_progress = TRUE) {

  #Check data input ------------------------------------------------------------
  if(!inherits(x, what = "matrix") && !inherits(x, what = "data.frame"))
    stop("The input x needs to be a matrix or dataframe.")
  if(inherits(x, what = "data.frame"))
    x = data.matrix(x)
  if(ncol(x) < 2)
    stop("The matrix x should have more than one variable (columns).")
  if(nrow(x) < 2)
    stop("The matrix x should have more than one observation (rows).")

  if(!inherits(y, what = "matrix") && !inherits(y, what = "data.frame"))
    stop("The input y needs to be a matrix or dataframe.")
  if(inherits(y, what = "data.frame"))
    y = data.matrix(y)
  if(ncol(y) < 2)
    stop("The matrix y should have more than one variable (columns).")
  if(nrow(y) < 2)
    stop("The matrix y should have more than one observation (rows).")

  if(ncol(x) != ncol(y))
    stop("The matrix x should have as many variables (columns) as the matrix y.")


  #Check model input -----------------------------------------------------------
  model = check_compare_model(x = x,
                              y = y,
                              paired = paired,
                              variable_type = variable_type,
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
                              interaction_scale = interaction_scale,
                              threshold_alpha = threshold_alpha,
                              threshold_beta = threshold_beta)

  # ----------------------------------------------------------------------------
  # The vector variable_type is now coded as boolean.
  # Ordinal (variable_bool == TRUE) or Blume-Capel (variable_bool == FALSE)
  # ----------------------------------------------------------------------------
  variable_bool = model$variable_bool
  # ----------------------------------------------------------------------------

  reference_category = model$reference_category
  main_difference_prior = model$main_difference_prior
  pairwise_difference_prior = model$pairwise_difference_prior
  inclusion_probability_difference = model$inclusion_probability_difference

  #Check Gibbs input -----------------------------------------------------------
  if(abs(iter - round(iter)) > .Machine$double.eps)
    stop("Parameter ``iter'' needs to be a positive integer.")
  if(iter <= 0)
    stop("Parameter ``iter'' needs to be a positive integer.")
  if(abs(burnin - round(burnin)) > .Machine$double.eps || burnin < 0)
    stop("Parameter ``burnin'' needs to be a non-negative integer.")
  if(burnin <= 0)
    stop("Parameter ``burnin'' needs to be a positive integer.")

  #Check na.action -------------------------------------------------------------
  na.action_input = na.action
  na.action = try(match.arg(na.action), silent = TRUE)
  if(inherits(na.action, what = "try-error"))
    stop(paste0("The na.action argument should equal listwise or impute, not ",
                na.action_input,
                "."))
  #Check save ------------------------------------------------------------------
  save_input = save
  save = as.logical(save)
  if(is.na(save))
    stop(paste0("The save argument should equal TRUE or FALSE, not ",
                save_input,
                "."))

  #Check display_progress ------------------------------------------------------
  display_progress = as.logical(display_progress)
  if(is.na(display_progress))
    stop("The display_progress argument should equal TRUE or FALSE.")

  #Format the data input -------------------------------------------------------
  data = compare_reformat_data(x = x,
                               y = y,
                               na.action = na.action,
                               variable_bool = variable_bool,
                               reference_category = reference_category,
                               paired = paired)
  x = data$x
  y = data$y

  no_categories = data$no_categories
  missing_index_gr1 = data$missing_index_gr1
  missing_index_gr2 = data$missing_index_gr2
  na_impute = data$na_impute
  reference_category = data$reference_category

  no_variables = ncol(x)
  no_interactions = no_variables * (no_variables - 1) / 2
  no_thresholds = sum(no_categories)

  #Precompute the number of observations per category for each variable --------
  n_cat_obs_gr1 = n_cat_obs_gr2 = matrix(0,
                                         nrow = max(no_categories) + 1,
                                         ncol = no_variables)
  for(variable in 1:no_variables) {
    for(category in 0:no_categories[variable]) {
      n_cat_obs_gr1[category + 1, variable] = sum(x[, variable] == category)
      n_cat_obs_gr2[category + 1, variable] = sum(y[, variable] == category)
    }
  }

  #Precompute the sufficient statistics for the two Blume-Capel parameters -----
  sufficient_blume_capel_gr1 = matrix(0, nrow = 2, ncol = no_variables)
  sufficient_blume_capel_gr2 = matrix(0, nrow = 2, ncol = no_variables)
  if(any(!variable_bool)) {
    # Ordinal (variable_bool == TRUE) or Blume-Capel (variable_bool == FALSE)
    bc_vars = which(!variable_bool)
    for(i in bc_vars) {
      sufficient_blume_capel_gr1[1, i] = sum(x[, i])
      sufficient_blume_capel_gr1[2, i] = sum((x[, i] - reference_category[i]) ^ 2)
      sufficient_blume_capel_gr2[1, i] = sum(y[, i])
      sufficient_blume_capel_gr2[2, i] = sum((y[, i] - reference_category[i]) ^ 2)
    }
  }

  # Index vector used to sample interactions in a random order -----------------
  Index = matrix(0,
                 nrow = no_variables * (no_variables - 1) / 2,
                 ncol = 3)
  cntr = 0
  for(variable1 in 1:(no_variables - 1)) {
    for(variable2 in (variable1 + 1):no_variables) {
      cntr =  cntr + 1
      Index[cntr, 1] = cntr
      Index[cntr, 2] = variable1
      Index[cntr, 3] = variable2
    }
  }

  #The Metropolis within Gibbs sampler -----------------------------------------
  out = compare_gibbs_sampler(observations_gr1 = x,
                              observations_gr2 = y,
                              no_categories = no_categories,
                              interaction_scale = interaction_scale,
                              cross_lagged_scale = cross_lagged_scale,
                              pairwise_difference_scale = pairwise_difference_scale,
                              main_difference_scale = main_difference_scale,
                              pairwise_difference_prior = pairwise_difference_prior,
                              main_difference_prior = main_difference_prior,
                              inclusion_probability_difference = inclusion_probability_difference,
                              pairwise_beta_bernoulli_alpha = pairwise_beta_bernoulli_alpha,
                              pairwise_beta_bernoulli_beta = pairwise_beta_bernoulli_beta,
                              main_beta_bernoulli_alpha = main_beta_bernoulli_alpha,
                              main_beta_bernoulli_beta = main_beta_bernoulli_beta,
                              Index = Index,
                              iter = iter,
                              burnin = burnin,
                              n_cat_obs_gr1 = n_cat_obs_gr1,
                              n_cat_obs_gr2 = n_cat_obs_gr2,
                              sufficient_blume_capel_gr1 = sufficient_blume_capel_gr1,
                              sufficient_blume_capel_gr2 = sufficient_blume_capel_gr2,
                              threshold_alpha = threshold_alpha,
                              threshold_beta = threshold_beta,
                              na_impute = na_impute,
                              missing_index_gr1 = missing_index_gr1,
                              missing_index_gr2 = missing_index_gr2,
                              variable_bool = variable_bool,
                              reference_category = reference_category,
                              paired = paired,
                              save = save,
                              display_progress = display_progress)

  #Preparing the output --------------------------------------------------------
  arguments = list(
    no_variables = no_variables,
    no_cases_gr1 = nrow(x),
    no_cases_gr2 = nrow(y),
    na_impute = na_impute,
    variable_type = variable_type,
    iter = iter,
    burnin = burnin,
    interaction_scale = interaction_scale,
    threshold_alpha = threshold_alpha,
    threshold_beta = threshold_beta,
    pairwise_difference_prior = pairwise_difference_prior,
    main_difference_prior = main_difference_prior,
    inclusion_probability_difference = inclusion_probability_difference,
    pairwise_beta_bernoulli_alpha = pairwise_beta_bernoulli_alpha,
    pairwise_beta_bernoulli_beta = pairwise_beta_bernoulli_beta,
    main_beta_bernoulli_alpha = main_beta_bernoulli_alpha,
    main_beta_bernoulli_beta = main_beta_bernoulli_beta,
    main_difference_scale = main_difference_scale,
    pairwise_difference_scale = pairwise_difference_scale,
    paired = paired,
    na.action = na.action,
    save = save,
    version = packageVersion("bgms")
  )

  if(save == FALSE) {
    indicator = out$pairwise_difference_indicator
    main_difference_indicator = out$main_difference_indicator
    diag(indicator) = main_difference_indicator
    interactions = out$interactions
    pairwise_difference = out$pairwise_difference
    thresholds = out$thresholds
    main_difference = out$main_difference

    if(is.null(colnames(x))){
      data_columnnames = paste0("variable ", 1:no_variables)
    } else {
      data_columnnames <- colnames(x)
    }
    colnames(interactions) = data_columnnames
    rownames(interactions) = data_columnnames
    colnames(pairwise_difference) = data_columnnames
    rownames(pairwise_difference) = data_columnnames
    colnames(indicator) = data_columnnames
    rownames(indicator) = data_columnnames
    rownames(thresholds) = data_columnnames
    rownames(main_difference) = data_columnnames

    if(paired == TRUE) {
      cross_lagged = out$cross_lagged
      colnames(cross_lagged) = data_columnnames
      rownames(cross_lagged) = data_columnnames
    }

    colnames(thresholds) = paste0("category ", 1:max(no_categories))
    colnames(main_difference) = paste0("category ", 1:max(no_categories))

    arguments$data_columnnames = data_columnnames

    if(paired == TRUE) {
      output = list(indicator = indicator,
                    interactions = interactions,
                    pairwise_difference = pairwise_difference,
                    main_difference = main_difference,
                    thresholds = thresholds,
                    cross_lagged = cross_lagged,
                    arguments = arguments)
    } else {
      output = list(indicator = indicator,
                    interactions = interactions,
                    pairwise_difference = pairwise_difference,
                    main_difference = main_difference,
                    thresholds = thresholds,
                    arguments = arguments)
    }


    class(output) = c("bgmCompare")
    return(output)
  } else {
    pairwise_difference_indicator = out$pairwise_difference_indicator
    main_difference_indicator = out$main_difference_indicator
    pairwise_difference = out$pairwise_difference
    main_difference = out$main_difference
    interactions = out$interactions
    thresholds = out$thresholds
    if(paired == TRUE) {
      cross_lagged = out$cross_lagged
    }

    if(is.null(colnames(x))){
      data_columnnames <- 1:ncol(x)
    } else {
      data_columnnames <- colnames(x)
    }

    p <- ncol(x)
    names_bycol <- matrix(rep(data_columnnames, each = p), ncol = p)
    names_byrow <- matrix(rep(data_columnnames, each = p), ncol = p, byrow = T)
    names_comb <- matrix(paste0(names_byrow, "-", names_bycol), ncol = p)
    names_vec <- names_comb[lower.tri(names_comb)]

    colnames(pairwise_difference_indicator) = names_vec
    colnames(interactions) = names_vec
    colnames(pairwise_difference) = names_vec

    if(paired == TRUE) {
      names_vec = vector(length = ncol(x) (ncol(x) + 1) / 2)
      cntr = 0
      for(variable1 in 1:ncol(x)) {
        for(variable2 in variable1:ncol(x)) {
          cntr = cntr + 1
          names_vec[cntr] = paste0(colnames[variable1], "--", colnames[variable2])
        }
      }
      colnames(cross_lagged) = names_vec
    }

    names = character(length = sum(no_categories))
    cntr = 0
    for(variable in 1:no_variables) {
      for(category in 1:no_categories[variable]) {
        cntr = cntr + 1
        names[cntr] = paste0("threshold(",variable, ", ",category,")")
      }
    }
    colnames(main_difference_indicator) = data_columnnames
    colnames(thresholds) = names
    colnames(main_difference) = names

    dimnames(pairwise_difference_indicator) = list(Iter. = 1:iter, colnames(pairwise_difference_indicator))
    dimnames(main_difference_indicator) = list(Iter. = 1:iter, colnames(main_difference_indicator))
    dimnames(pairwise_difference) = list(Iter. = 1:iter, colnames(pairwise_difference))
    dimnames(main_difference) = list(Iter. = 1:iter, colnames(main_difference))
    dimnames(interactions) = list(Iter. = 1:iter, colnames(interactions))
    dimnames(thresholds) = list(Iter. = 1:iter, colnames(thresholds))

    if(paired == TRUE) {
      dimnames(cross_lagged) = list(Iter. = 1:iter, colnames(cross_lagged))
    }

    arguments$data_columnnames = data_columnnames

    if(paired == TRUE) {
      output = list(pairwise_difference_indicator = pairwise_difference_indicator,
                    main_difference_indicator = main_difference_indicator,
                    interactions = interactions,
                    pairwise_difference = pairwise_difference,
                    main_difference = main_difference,
                    thresholds = thresholds,
                    cross_lagged = cross_lagged,
                    arguments = arguments)
    } else {
      output = list(pairwise_difference_indicator = pairwise_difference_indicator,
                    main_difference_indicator = main_difference_indicator,
                    interactions = interactions,
                    pairwise_difference = pairwise_difference,
                    main_difference = main_difference,
                    thresholds = thresholds,
                    arguments = arguments)
    }

    class(output) = "bgmCompare"
    return(output)
  }
}
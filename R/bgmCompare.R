#' Bayesian variable selection or Bayesian estimation for differences in the
#' Markov random field model for binary and/or ordinal variables in independent
#' samples.
#'
#' @description
#' The \code{bgmCompare} function estimates the pseudoposterior distribution of
#' the parameters of a Markov Random Field model for mixed binary and ordinal
#' variables, and the differences in pairwise interactions and category thresholds
#' between groups. The groups are assumed to be \code{G} independent samples.
#'
#' @details
#' The pairwise interactions between the variables \eqn{i}{i} and \eqn{j}{j} are
#' modeled as
#' \deqn{\boldsymbol{\theta}_{\text{ij}} = \phi_{\text{ij}} + \boldsymbol{\delta}_{\text{ij}},}{
#' \boldsymbol{\theta}_{\text{ij}} = \phi_{\text{ij}} + \boldsymbol{\delta}_{\text{ij}},}
#' where \eqn{\boldsymbol{\theta}_{\text{ij}}}{\boldsymbol{\theta}_{\text{ij}}}
#' is the vector of pairwise interaction parameters of length \eqn{G}{G},
#' \eqn{G}{G} is the number of groups, \eqn{\phi_{\text{ij}}}{\phi_{\text{ij}}}
#' is the overall pairwise interaction that is considered nuisance, and attention
#' is focused on the group differences from the overall pairwise interactions in
#' the vector \eqn{\boldsymbol{\delta}_{\text{ij}}}{\boldsymbol{\delta}_{\text{ij}}}.
#' The group differences from the overall pairwise interaction are constrained
#' to sum to zero for identification.
#'
#' The \code{bgmCompare} function supports two types of ordinal variables, which
#' can be mixed. The default ordinal variable introduces a threshold parameter
#' for each category except the lowest category. For this variable type,
#' the vector of threshold parameters for variable \eqn{i}{i}, category
#' \eqn{c}{c}, is modeled as
#' \deqn{\boldsymbol{\mu}_{\text{ic}} = \tau_{\text{ic}} + \boldsymbol{\epsilon}_{\text{ic}},}{\boldsymbol{\mu}_{\text{ic}} = \tau_{\text{ic}} + \boldsymbol{\epsilon}_{\text{ic}},}
#' where category threshold parameter \eqn{\tau_{\text{ic}}}{\tau_{\text{ic}}}
#' denotes an overall effect that is considered nuisance, and attention is
#' focused on the vector of group differences from the overall category
#' threshold parameter \eqn{\boldsymbol{\epsilon}_{\text{ic}}}{\boldsymbol{\epsilon}_{\text{ic}}}.
#' The group differences from the overall category threshold are constrained to
#' sum to zero for identification.
#'
#' The Blume-Capel ordinal variable assumes that there is a specific reference
#' category, such as ``neutral'' in a Likert scale, and responses are scored
#' according to their distance from this reference category. The vector of
#' category threshold parameters are modeled as
#' \deqn{\boldsymbol{\mu}_{\text{ic}} = (\tau_{\text{i1}} + \boldsymbol{\epsilon}_{\text{i1}}) \times \text{c} + (\tau_{\text{i2}} + \boldsymbol{\epsilon}_{\text{i2}}) \times (\text{c} - \text{r})^2,}{ \boldsymbol{\mu}_{\text{ic}} = (\tau_{\text{i1}} + \boldsymbol{\epsilon}_{\text{i1}}) \times \text{c} + (\tau_{\text{i2}} + \boldsymbol{\epsilon}_{\text{i2}}) \times (\text{c} - \text{r})^2,}
#' where \eqn{\tau_{\text{i1}}}{\tau_{\text{i1}}} and \eqn{\tau_{\text{i2}}}{\tau_{\text{i2}}}
#' denote overall effects that are considered nuisance, and attention is focused
#' on the vectors of group differences from these overall effects;
#' \eqn{\boldsymbol{\epsilon}_{\text{i1}}}{\boldsymbol{\epsilon}_{\text{i1}}} and
#' \eqn{\boldsymbol{\epsilon}_{\text{i2}}}{\boldsymbol{\epsilon}_{\text{i2}}}.
#' The group differences from the overall Blume-Capel parameters are constrained
#' to sum to zero for identification.
#'
#' Bayesian variable selection is used to model the presence or absence of the
#' vectors of difference parameters \eqn{\boldsymbol{\delta}}{\boldsymbol{\delta}}
#' and \eqn{\boldsymbol{\epsilon}}{\boldsymbol{\epsilon}}, which allow us to
#' test hypotheses about parameter differences or equivalence across groups.
#' Independent spike and slab priors are specified for these (vectors of)
#' difference parameters. The spike and slab prior uses a binary indicator
#' variable to select the vector of difference parameters, assigning them a
#' diffuse (multivariate) Cauchy prior with an optional scaling parameter if
#' selected, or setting the vector of difference parameters to zero if not
#' selected.
#'
#' The function offers two models for the probabilistic inclusion of parameter
#' differences:
#' \itemize{
#'   \item \strong{Bernoulli Model}: This model assigns a fixed probability of
#'   selecting a vector of parameter differences, treating them as independent
#'   events. A probability of 0.5 indicates no preference, giving equal prior
#'   weight to all configurations.
#'   \item \strong{Beta-Bernoulli Model}: Introduces a beta distribution prior
#'   for the inclusion probability that models the complexity of the
#'   configuration of the vector of difference indicators. When the shape
#'   parameters of the beta distribution are set to 1, the model assigns the
#'   same prior weight to the number of difference vectors present (i.e., a
#'   configuration with two difference or with four differences is a priori
#'   equally likely).
#' }
#' Inclusion probabilities can be specified for pairwise interactions with
#' \code{pairwise_difference_probability} and for category thresholds with
#' \code{threshold_difference_probability}.
#'
#' The pairwise interaction parameters \eqn{\theta}{\theta}, the category
#' threshold parameters \eqn{\tau}{\tau} are considered nuisance parameters that
#' are common to all models. The pairwise interaction parameters \eqn{\theta}{\theta}
#' are assigned a diffuse Cauchy prior with an optional scaling parameter. The
#' exponent of the category threshold parameters \eqn{\tau}{\tau} are assigned
#' beta-prime distribution with optional scale values.
#'
#' @param x A data frame or matrix with \eqn{n}{n} rows and \code{p} columns
#' containing binary and ordinal responses. Regular ordinal variables are
#' recoded as non-negative integers \code{(0, 1, ..., m)} if not already done.
#' Unobserved categories are collapsed into other categories after recoding
#' (i.e., if category 1 is unobserved, the data are recoded from (0, 2) to
#' (0, 1)). Blume-Capel ordinal variables are also coded as non-negative
#' integers if not already done. However, since ``distance'' from the reference
#' category plays an important role in this model, unobserved categories are not
#' collapsed after recoding.
#' @param y A data frame or matrix with \eqn{n}{n} rows and \code{p} columns
#' containing binary and ordinal responses, similar to \code{x}. Optional
#' argument for two-group designs, where \code{x} contains the data for Group 1
#' and \code{y} contains the data for Group 2.
#' @param g A vector of length \eqn{n}{n} containing group membership indicators
#' for the rows of \code{x}. Optional argument for two-group designs (see the
#' documentation for \code{y}), but required for multi-group designs. This
#' argument is ignored if there is an input for \code{y}.
#' @param difference_selection Logical. If \code{TRUE}, \code{bgmCompare} will
#' model the inclusion or exclusion of the between samples parameter
#' differences; if \code{FALSE}, it will estimate all between-sample
#' parameter differences. Default is \code{TRUE}.
#' @param main_difference_model A string specifying options for how the
#' bgmCompare function should handle the comparison of threshold parameters when
#' the observed categories in the samples do not match. The "Collapse" option
#' tells bgmCompare to collapse the two categories into one (for the data set
#' where both categories were observed). The "Constrain" option sets the
#' difference between the category thresholds in the two data sets to zero if
#' the category is not observed in one of the two data sets. The "Free" option
#' tells bgmCompare to estimate a separate set of thresholds in the two samples
#' and to not model their differences.
#' @param variable_type A string or vector specifying the type of variables
#' in \code{x} (and \code{y}). Supported types are "ordinal" and "blume-capel",
#' with binary variables treated as "ordinal". Default is "ordinal".
#' @param reference_category The reference category in the Blume-Capel model.
#' Should be an integer within the range of integer values observed for the
#' "blume-capel" variable. Can be a single number that sets the reference
#' category for all Blume-Capel variables at once, or a vector of length
#' \code{p}, where the \code{i}-th element is the reference category for the
#' \code{i} variable if it is a Blume-Capel variable, and elements for other
#' variable types are ignored. The value of the reference category is also
#' recoded when `bgmCompare` recodes the corresponding observations. Required
#' only if there is at least one variable of type "blume-capel".
#' @param pairwise_difference_scale The scale of the Cauchy distribution that is
#' used as the prior for the pairwise difference parameters. Defaults to
#' \code{1}.
#' @param main_difference_scale The scale of the Cauchy distribution that
#' is used as the prior for the threshold difference parameters. Defaults to
#' \code{1}.
#' @param pairwise_difference_prior A character string that specifies the model
#' to use for the  inclusion probability of pairwise differences. Options are
#' "Bernoulli" or "Beta-Bernoulli". Default is "Bernoulli".
#' @param main_difference_prior A character string that specifies the model
#' to use for the  inclusion probability of threshold differences. Options are
#' "Bernoulli" or "Beta-Bernoulli". Default is "Bernoulli".
#' @param pairwise_difference_probability The inclusion probability for a
#' pairwise difference in the Bernoulli model. Can be a single probability or a
#' matrix of \code{p} rows and \code{p} columns specifying the probability of a
#' difference for each edge pair. Defaults to 0.5.
#' @param main_difference_probability The inclusion probability for a
#' threshold difference in the Bernoulli model. Defaults to 0.5, implying no
#' prior preference. Can be a single probability or a vector of length \code{p}
#' specifying the probability of a difference between the category thresholds
#' for each variable. Defaults to 0.5.
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
#' @param threshold_alpha,threshold_beta The shape parameters of the beta-prime
#' prior density for the nuisance threshold parameters. Must be positive values.
#' If the two values are equal, the prior density is symmetric about zero. If
#' \code{threshold_beta} is greater than \code{threshold_alpha}, the distribution
#' is left-skewed, and if \code{threshold_beta} is less than
#' \code{threshold_alpha}, it is right-skewed. Smaller values tend to result in
#' more diffuse prior distributions.
#' @param iter The function uses a Gibbs sampler to sample from the posterior
#' distribution of the model parameters and indicator variables. How many
#' iterations should this Gibbs sampler run? The default of \code{1e4} is for
#' illustrative purposes. For stable estimates, it is recommended to run the
#' Gibbs sampler for at least \code{1e5} iterations.
#' @param burnin The number of iterations of the Gibbs sampler before saving its
#' output. Since it may take some time for the Gibbs sampler to converge to the
#' posterior distribution, it is recommended not to set this number too low.
#' When \code{difference_selection = TRUE}, the bgm function will perform
#' \code{2 * burnin} iterations, first \code{burnin} iterations without
#' difference selection, then \code{burnin} iterations with difference
#' selection. This helps ensure that the Markov chain used for estimation starts
#' with good parameter values and that the adaptive MH proposals are properly
#' calibrated.
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
#'    \item \code{difference_pairwise}: A list of length \code{G}, the number of
#'    groups, each containing a matrix with \code{p} rows and \code{p} columns
#'    containing model-averaged posterior means of the
#'    differences in pairwise interactions for that group.
#'    \item \code{difference_threshold}: A list of length \code{G}, the number of
#'    groups, each containing a matrix with \code{p} rows and \code{max(m)}
#'    columns containing model-averaged posterior means of the
#'    differences in category thresholds for that group.
#'    \item \code{interactions}: A matrix with \code{p} rows and \code{p} columns,
#'    containing posterior means of the nuisance pairwise interactions.
#'    \item \code{thresholds}: A matrix with \code{p} rows and \code{max(m)}
#'    columns containing the posterior means of the nuisance category thresholds.
#'    In the case of ``blume-capel'' variables, the first entry is the parameter
#'    for the linear effect and the second entry is the parameter for the
#'    quadratic effect, which models the offset to the reference category.
#'  }
#'
#' If \code{save = TRUE}, the result is a list of class ``bgmCompare''
#' containing the following matrices:
#'  \itemize{
#'    \item \code{indicator_pairwise}: A matrix with \code{iter} rows and
#'    \code{p * (p - 1) / 2} columns containing the inclusion indicators for the
#'    differences in pairwise interactions from each iteration of the Gibbs
#'    sampler.
#'    \item \code{difference_pairwise}: A list of length \code{G}, the number of
#'    groups, each containing a matrix with \code{iter} rows and
#'    \code{p * (p - 1) / 2} columns, containing parameter states for the
#'    differences in pairwise interactions from each iteration of the Gibbs
#'    sampler for that group.
#'    \item \code{indicator_threshold}: A matrix with \code{iter} rows and
#'    \code{sum(m)} columns, containing the inclusion indicators for the
#'    differences in category thresholds from each iteration of the Gibbs
#'    sampler.
#'    \item \code{difference_threshold}: A list of length \code{G}, the number of
#'    groups, each containing a matrix with \code{iter} rows and
#'    \code{sum(m)} columns, containing the parameter states for the differences
#'    in category thresholds from each iteration of the Gibbs sampler for that
#'    group.
#'    \item \code{interactions}: A matrix with \code{iter} rows and
#'    \code{p * (p - 1) / 2} columns, containing parameter states for the
#'    nuisance pairwise interactions in each iteration of the Gibbs sampler.
#'    \item \code{thresholds}: A matrix with \code{iter} rows and \code{sum(m)}
#'    columns, containing parameter states for the nuisance category thresholds
#'    in each iteration of the Gibbs sampler.
#'  }
#' Column averages of these matrices provide the model-averaged posterior means.
#'
#' In addition to the results of the analysis, the output lists some of the
#' arguments of its call. This is useful for post-processing the results.
#'
#' @importFrom utils packageVersion
#' @export
bgmCompare = function(x,
                      y,
                      g,
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
                      display_progress = TRUE) {

  ttest = hasArg(y)

  #Check data input ------------------------------------------------------------
  if(!ttest & !hasArg(g))
    stop(paste0("For multi-group designs, the bgmCompare function requires input for\n",
                "either y (group 2 data) or g (group indicator)."))

  if(!inherits(x, what = "matrix") && !inherits(x, what = "data.frame"))
    stop("The input x needs to be a matrix or dataframe.")
  if(inherits(x, what = "data.frame"))
    x = data.matrix(x)
  if(ncol(x) < 2)
    stop("The matrix x should have more than one variable (columns).")
  if(nrow(x) < 2)
    stop("The matrix x should have more than one observation (rows).")

  if(ttest) {
    if(!inherits(y, what = "matrix") && !inherits(y, what = "data.frame"))
      stop("The input y needs to be a matrix or dataframe.")
    if(inherits(y, what = "data.frame"))
      y = data.matrix(y)

    if(ncol(x) != ncol(y))
      stop("The matrix x should have as many variables (columns) as the matrix y.")

    if(nrow(y) < 2)
      stop("The matrix y should have more than one observation (rows).")
  }

  if(!ttest & hasArg(g)) {
    g = as.vector(g)
    if(anyNA(g))
      stop("The input g has missing values.")
    if(length(g) != nrow(x))
      stop("The input g needs to be a vector of length nrow(x).")

    unique_g = unique(g)
    if(length(unique_g) == 2) {
      y = x[g == unique_g[2]]
      x = x[g == unique_g[1]]
      ttest = TRUE
    }
  }

  #Check model input -----------------------------------------------------------
  if(!ttest) #True if either hasArg(y) or length(unique(g)) = 2
    y = NULL
  if(!hasArg(g))
    g = NULL

  model = check_compare_model(x = x,
                              y = y,
                              g = g,
                              ttest = ttest,
                              difference_selection = difference_selection,
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
                              threshold_beta = threshold_beta,
                              main_difference_model = main_difference_model)

  # ----------------------------------------------------------------------------
  # The vector variable_type is now coded as boolean.
  # Ordinal (variable_bool == TRUE) or Blume-Capel (variable_bool == FALSE)
  # ----------------------------------------------------------------------------
  ordinal_variable = model$variable_bool
  # ----------------------------------------------------------------------------

  reference_category = model$reference_category
  main_difference_prior = model$main_difference_prior
  pairwise_difference_prior = model$pairwise_difference_prior
  inclusion_probability_difference = model$inclusion_probability_difference
  main_difference_model = model$main_difference_model
  independent_thresholds = (main_difference_model == "Free")

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
                               g = g,
                               ttest = ttest,
                               na.action = na.action,
                               variable_bool = ordinal_variable,
                               reference_category = reference_category,
                               main_difference_model = main_difference_model)
  x = data$x
  if(ttest == TRUE) {
    y = data$y

    no_categories_gr1 = data$no_categories[, 1]
    no_categories_gr2 = data$no_categories[, 2]
    missing_index_gr1 = data$missing_index_gr1
    missing_index_gr2 = data$missing_index_gr2

    no_obs_groups = c(nrow(x), nrow(y))
  } else {
    group = data$group
    no_obs_groups = tabulate(group)
    missing_index = data$missing_index
  }

  na_impute = data$na_impute
  reference_category = data$reference_category

  no_variables = ncol(x)
  no_interactions = no_variables * (no_variables - 1) / 2

  #Precompute the number of observations per category for each variable --------
  if(ttest == TRUE) {
    if(main_difference_model == "Free") {
      n_cat_obs_gr1 = n_cat_obs_gr2 = matrix(0,
                                             nrow = max(c(no_categories_gr1,
                                                          no_categories_gr2)) + 1,
                                             ncol = no_variables)
      for(variable in 1:no_variables) {
        for(category in 0:no_categories_gr1[variable]) {
          n_cat_obs_gr1[category + 1, variable] = sum(x[, variable] == category)
        }
        for(category in 0:no_categories_gr2[variable]) {
          n_cat_obs_gr2[category + 1, variable] = sum(y[, variable] == category)
        }
      }
    } else {
      n_cat_obs_gr1 = n_cat_obs_gr2 = matrix(0,
                                             nrow = max(no_categories_gr1) + 1,
                                             ncol = no_variables)
      for(variable in 1:no_variables) {
        for(category in 0:no_categories_gr1[variable]) {
          n_cat_obs_gr1[category + 1, variable] = sum(x[, variable] == category)
          n_cat_obs_gr2[category + 1, variable] = sum(y[, variable] == category)
        }
      }
    }

    #Precompute the sufficient statistics for the two Blume-Capel parameters -----
    sufficient_blume_capel_gr1 = matrix(0, nrow = 2, ncol = no_variables)
    sufficient_blume_capel_gr2 = matrix(0, nrow = 2, ncol = no_variables)
    if(any(!ordinal_variable)) {
      # Ordinal (variable_bool == TRUE) or Blume-Capel (variable_bool == FALSE)
      bc_vars = which(!ordinal_variable)
      for(i in bc_vars) {
        sufficient_blume_capel_gr1[1, i] = sum(x[, i])
        sufficient_blume_capel_gr1[2, i] = sum((x[, i] - reference_category[i]) ^ 2)
        sufficient_blume_capel_gr2[1, i] = sum(y[, i])
        sufficient_blume_capel_gr2[2, i] = sum((y[, i] - reference_category[i]) ^ 2)
      }
    }
  } else {
    n_cat_obs = list()
    for(g in group) {
      n_cat_obs_gr = matrix(0,
                            nrow = max(no_categories[, g]) + 1,
                            ncol = no_variables)

      for(variable in 1:no_variables) {
        for(category in 0:no_categories_gr[variable, g]) {
          n_cat_obs_gr[category + 1, variable] = sum(x[group == g, variable] == category)
        }
      }

      n_cat_obs[[g]] = n_cat_obs_gr
    }

    #Precompute the sufficient statistics for the two Blume-Capel parameters -----
    sufficient_blume_capel = list()
    for(g in group) {
      sufficient_blume_capel_gr = matrix(0, nrow = 2, ncol = no_variables)
      if(any(!ordinal_variable)) {
        # Ordinal (variable_bool == TRUE) or Blume-Capel (variable_bool == FALSE)
        bc_vars = which(!ordinal_variable)
        for(i in bc_vars) {
          sufficient_blume_capel_gr[1, i] = sum(x[group == g, i])
          sufficient_blume_capel_gr[2, i] = sum((x[group == g, i] - reference_category[i]) ^ 2)
        }
      }
      sufficient_blume_capel[[g]] = sufficient_blume_capel_gr
    }
  }

  # Index vector used to sample interactions in a random order -----------------
  Index = matrix(0,
                 nrow = no_interactions,
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
  if(ttest == TRUE) {
    out = compare_gibbs_sampler(observations_gr1 = x,
                                observations_gr2 = y,
                                no_categories_gr1 = no_categories_gr1,
                                no_categories_gr2 = no_categories_gr2,
                                interaction_scale = interaction_scale,
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
                                ordinal_variable = ordinal_variable,
                                reference_category = reference_category,
                                independent_thresholds = independent_thresholds,
                                save = save,
                                display_progress = display_progress,
                                difference_selection = difference_selection)
  } else {
    out = compare_gibbs_sampler_anova(observations = x,
                                      group = group,
                                      no_categories = no_categories,
                                      no_obs_groups = no_obs_groups,
                                      interaction_scale = interaction_scale,
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
                                      n_cat_obs = n_cat_obs,
                                      sufficient_blume_capel = sufficient_blume_capel,
                                      threshold_alpha = threshold_alpha,
                                      threshold_beta = threshold_beta,
                                      na_impute = na_impute,
                                      missing_index = missing_index,
                                      ordinal_variable = ordinal_variable,
                                      reference_category = reference_category,
                                      independent_thresholds = independent_thresholds,
                                      save = save,
                                      display_progress = display_progress,
                                      difference_selection = difference_selection)
  }


  #Preparing the output --------------------------------------------------------
  arguments = list(
    no_variables = no_variables,
    no_cases = no_obs_groups,
    na_impute = na_impute,
    variable_type = variable_type,
    iter = iter,
    burnin = burnin,
    difference_selection = difference_selection,
    interaction_scale = interaction_scale,
    threshold_alpha = threshold_alpha,
    threshold_beta = threshold_beta,
    main_difference_model = main_difference_model,
    pairwise_difference_prior = pairwise_difference_prior,
    main_difference_prior = main_difference_prior,
    inclusion_probability_difference = inclusion_probability_difference,
    pairwise_beta_bernoulli_alpha = pairwise_beta_bernoulli_alpha,
    pairwise_beta_bernoulli_beta = pairwise_beta_bernoulli_beta,
    main_beta_bernoulli_alpha = main_beta_bernoulli_alpha,
    main_beta_bernoulli_beta = main_beta_bernoulli_beta,
    main_difference_scale = main_difference_scale,
    pairwise_difference_scale = pairwise_difference_scale,
    na.action = na.action,
    save = save,
    version = packageVersion("bgms"),
    independent_thresholds = independent_thresholds
  )

  if(save == FALSE) {
    indicator = out$pairwise_difference_indicator
    if(independent_thresholds == TRUE) {
      if(ttest == TRUE) {
        thresholds_gr1 = out$thresholds_gr1
        thresholds_gr2 = out$thresholds_gr2
      } else {
        thresholds = out$thresholds
      }
    } else {
      main_difference_indicator = out$main_difference_indicator
      diag(indicator) = main_difference_indicator
      if(ttest == TRUE) {
        thresholds = out$thresholds
      } else {
        thresholds = out$thresholds[[1]]
      }
      main_difference = out$main_difference
    }

    interactions = out$interactions
    pairwise_difference = out$pairwise_difference

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

    if(independent_thresholds == TRUE) {
      if(ttest == TRUE) {
        rownames(thresholds_gr1) = data_columnnames
        rownames(thresholds_gr2) = data_columnnames
        colnames(thresholds_gr1) = paste0("category ", 1:max(no_categories_gr1))
        colnames(thresholds_gr2) = paste0("category ", 1:max(no_categories_gr2))
      } else {
        for(g in group) {
           rownames(thresholds[[g]]$thresholds) = data_columnnames
           colnames(thresholds[[g]]$thresholds) = paste0("category ", 1:max(no_categories[, g]))
        }
      }
    } else {
      rownames(thresholds) = data_columnnames
      rownames(main_difference) = data_columnnames
      if(ttest == TRUE) {
        colnames(thresholds) = paste0("category ", 1:max(no_categories_gr1))
        colnames(main_difference) = paste0("category ", 1:max(no_categories_gr1))
      } else {
        colnames(thresholds) = paste0("category ", 1:max(no_categories))
        colnames(main_difference) = paste0("category ", 1:max(no_categories))
      }
    }
    arguments$data_columnnames = data_columnnames

    if(independent_thresholds == TRUE) {
      if(ttest == TRUE) {
        output = list(indicator = indicator,
                      interactions = interactions,
                      pairwise_difference = pairwise_difference,
                      thresholds_gr1 = thresholds_gr1,
                      thresholds_gr2 = thresholds_gr2,
                      arguments = arguments)
      } else {
        output = list(indicator = indicator,
                      interactions = interactions,
                      pairwise_difference = pairwise_difference,
                      thresholds = thresholds,
                      arguments = arguments)
      }

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
    pairwise_difference = out$pairwise_difference
    interactions = out$interactions

    if(independent_thresholds == TRUE) {
      if(ttest == TRUE) {
        thresholds_gr1 = out$thresholds_gr1
        thresholds_gr2 = out$thresholds_gr2
      } else {
        thresholds = out$thresholds
      }
    } else {
      main_difference_indicator = out$main_difference_indicator
      main_difference = out$main_difference
      thresholds = out$thresholds
    }

    if(is.null(colnames(x))){
      data_columnnames <- 1:ncol(x)
    } else {
      data_columnnames <- colnames(x)
    }
    arguments$data_columnnames = data_columnnames

    p <- ncol(x)
    names_bycol <- matrix(rep(data_columnnames, each = p), ncol = p)
    names_byrow <- matrix(rep(data_columnnames, each = p), ncol = p, byrow = T)
    names_comb <- matrix(paste0(names_byrow, "-", names_bycol), ncol = p)
    names_vec <- names_comb[lower.tri(names_comb)]

    colnames(pairwise_difference_indicator) = names_vec
    colnames(interactions) = names_vec
    if(ttest == TRUE) {
      dimnames(pairwise_difference) = list(Iter. = 1:iter, names_vec)
    } else {
      for(g in group) {
        dimnames(pairwise_difference[[g]]$pairwise_difference) = list(Iter. = 1:iter, names_vec)
      }
    }

    dimnames(pairwise_difference_indicator) = list(Iter. = 1:iter, colnames(pairwise_difference_indicator))
    dimnames(interactions) = list(Iter. = 1:iter, colnames(interactions))

    if(independent_thresholds == TRUE) {
      if(ttest == TRUE) {
        names = character(length = sum(no_categories_gr1))
        cntr = 0
        for(variable in 1:no_variables) {
          for(category in 1:no_categories_gr1[variable]) {
            cntr = cntr + 1
            names[cntr] = paste0("threshold(",variable, ", ",category,")")
          }
        }
        dimnames(thresholds_gr1) = list(Iter. = 1:iter, names)

        names = character(length = sum(no_categories_gr2))
        cntr = 0
        for(variable in 1:no_variables) {
          for(category in 1:no_categories_gr2[variable]) {
            cntr = cntr + 1
            names[cntr] = paste0("threshold(",variable, ", ",category,")")
          }
        }
        dimnames(thresholds_gr2) = list(Iter. = 1:iter, names)

      } else {
        for(g in group) {
          names = character(length = sum(no_categories[, g]))
          cntr = 0
          for(variable in 1:no_variables) {
            for(category in 1:no_categories[variable, g]) {
              cntr = cntr + 1
              names[cntr] = paste0("threshold(",variable, ", ",category,")")
            }
          }
          dimnames(thresholds[[g]]$thresholds) = list(Iter. = 1:iter, names)
        }
      }
    } else {
      if(ttest == TRUE) {
        names = character(length = sum(no_categories_gr1))
        cntr = 0
        for(variable in 1:no_variables) {
          for(category in 1:no_categories_gr1[variable]) {
            cntr = cntr + 1
            names[cntr] = paste0("threshold(",variable, ", ",category,")")
          }
        }
      } else {
        names = character(length = sum(no_categories[1,]))
        cntr = 0
        for(variable in 1:no_variables) {
          for(category in 1:no_categories[variable, 1]) {
            cntr = cntr + 1
            names[cntr] = paste0("threshold(",variable, ", ",category,")")
          }
        }

      }
      colnames(main_difference_indicator) = data_columnnames
      dimnames(thresholds) = list(Iter. = 1:iter, names)

      dimnames(main_difference_indicator) = list(Iter. = 1:iter, colnames(main_difference_indicator))

      if(ttest == TRUE) {
        dimnames(main_difference) = list(Iter. = 1:iter, colnames(main_difference))
      } else {
        for(g in group) {
          dimnames(main_difference[[g]]$main_difference) = list(Iter. = 1:iter, names)
        }
      }
    }

    if(independent_thresholds == TRUE) {
      if(ttest == TRUE) {
        output = list(pairwise_difference_indicator = pairwise_difference_indicator,
                      interactions = interactions,
                      pairwise_difference = pairwise_difference,
                      thresholds_gr1 = thresholds_gr1,
                      thresholds_gr2 = thresholds_gr2,
                      arguments = arguments)
      } else {
        output = list(pairwise_difference_indicator = pairwise_difference_indicator,
                      interactions = interactions,
                      pairwise_difference = pairwise_difference,
                      thresholds = thresholds,
                      arguments = arguments)
      }
    } else {
      output = list(pairwise_difference_indicator = pairwise_difference_indicator,
                    main_difference_indicator = main_difference_indicator,
                    interactions = interactions,
                    pairwise_difference = pairwise_difference,
                    main_difference = main_difference,
                    thresholds = thresholds,
                    arguments = arguments)
    }

    class(output) = c("bgmCompare", "bgms")
    return(output)
  }
}
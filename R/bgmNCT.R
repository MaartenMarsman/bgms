#' A Bayesian comparison of the network structure for two groups for a graphical
#' model of mixed binary and ordinal variables using MCMC.
#'
#' The function \code{bgmNCT} explores the joint pseudoposterior distribution of
#' the parameters of a Markov Random Field model for mixed binary and ordinal
#' variables, group differences in pairwise interactions for a two-sample
#' design, and indicator variables for these group differences.
#'
#' The pairwise interactions between variables \eqn{i}{i} and \eqn{j}{j} are
#' modeled as follows
#' \deqn{\sigma_{\text{ij}} = \eta_{\text{ij}} - \delta_{\text{ij}} / 2,}{\sigma_{\text{ij}} = \eta_{\text{ij}} -  \delta_{\text{ij}} / 2,}
#' in one group and as
#' \deqn{\sigma_{\text{ij}} = \eta_{\text{ij}} +  \delta_{\text{ij}} / 2,}{\sigma_{\text{ij}} = \eta_{\text{ij}} + \delta_{\text{ij}} / 2,}
#' in the other group. The parameter \eqn{\eta_{\text{ij}}}{\eta_{\text{ij}}}
#' denotes an overall interaction parameter and is considered to be a nuisance
#' parameter, and \eqn{\delta_{\text{ij}}}{\delta_{\text{ij}}} denotes the a
#' pairwise difference parameter, the difference in the pairwise interaction
#' between the two groups.
#'
#' Bayesian variable selection is used to model the presence or absence of the
#' difference in pairwise interactions between the two groups. It imposes a
#' discrete spike and slab prior distribution on the differences \eqn{\delta_{\text{ij}}}{\delta_{\text{ij}}}.
#' By formulating it as a mixture of mutually singular distributions, the function
#' can use a combination of Metropolis-Hastings and Gibbs sampling to create a
#' Markov chain that has the joint posterior distribution as an invariant. The
#' current option for the slab distribution is a Cauchy with an optional scaling
#' parameter. Two prior distributions are implemented for the indicator variables
#' for the pairwise difference parameters (i.e., the prior probability that there
#' is a pairwise difference); the Bernoulli prior and the Beta-Bernoulli prior.
#'
#' Currently, \code{bgmNCT} supports two types of ordinal variables. The regular,
#' default, ordinal variable type has gets its own threshold parameter. The
#' Blume-Capel ordinal variable assumes that there is a specific reference
#' category, such as``neutral'' in a Likert scale, and responses are scored
#' according to their distance from that reference category. Specifically, the
#' Blume-Capel model specifies the following quadratic model for the threshold
#' parameters:
#' \deqn{\mu_{\text{c}} = \alpha \times \text{c} + \beta \times (\text{c} - \text{r})^2,}{{\mu_{\text{c}} = \alpha \times \text{c} + \beta \times (\text{c} - \text{r})^2,}}
#' where \eqn{\mu_{\text{c}}}{\mu_{\text{c}}} is the threshold for category c.
#' The parameter \eqn{\alpha}{\alpha} models a linear trend across categories,
#' such that \eqn{\alpha > 0}{\alpha > 0} results in an increasing number of
#' observations in higher response categories and \eqn{\alpha <0}{\alpha <0}
#' results in a decreasing number of observations in higher response categories.
#' The parameter \eqn{\beta}{\beta} models the response style in terms of an
#' offset with respect to the reference category \eqn{r}{r}; if \eqn{\beta<0}{\beta<0}
#' there is a preference to respond in the reference category (i.e., the model
#' introduces a penalty for responding in a category further away from the
#' reference_category category \code{r}), while if \eqn{\beta > 0}{\beta > 0}
#' there is preference to score in the extreme categories further away from the
#' reference_category category.
#'
#' The pairwise interaction parameters and the category parameters are considered
#' nuisance parameters that are common to all models. The prior distribution for
#' the pairwise interaction parameters is a Cauchy distribution with an optional
#' scaling value. The prior for the exponent of the category parameters is a
#' beta-prime distribution, which also has two optional scaling parameters.
#'
#' @param x A data frame or matrix with \code{n} rows and \code{p} columns
#' containing binary and ordinal variables for \code{n} independent observations
#' and \code{p} variables in the network. Regular binary and ordinal variables
#' are recoded as non-negative integers \code{(0, 1, ..., m)} if not already
#' done. Unobserved categories are collapsed into other categories after
#' recoding (i.e., if category 1 is unobserved, the data are recoded from
#' (0, 2) to (0, 1)). Blume-Capel ordinal variables are also coded as
#' non-negative integers if not already done. However, since ``distance'' to the
#' reference category plays an important role in this model, unobserved
#' categories are not collapsed after recoding.
#' @param variable_type What kind of variables are there in \code{x}? Can be a
#' single character string specifying the variable type of all \code{p}
#' variables at once or a vector of character strings of length \code{p}
#' specifying the type for each variable in \code{x} separately. Currently, bgm
#' supports ``ordinal'' and ``blume-capel''. Binary variables are automatically
#' treated as ``ordinal’’. Defaults to \code{variable_type = "ordinal"}.
#' @param group_indicator A vector of length \code{p} that indicates, for each
#' row in \code{x}, to which of two groups it belongs.
#' @param reference_category The reference category in the Blume-Capel model.
#' Should be an integer within the range of integer scores observed for the
#' ``blume-capel'' variable. Can be a single number specifying the reference
#' category for all Blume-Capel variables at once, or a vector of length
#' \code{p} where the \code{i}-th element contains the reference category for
#' variable \code{i} if it is Blume-Capel, and bgm ignores its elements for
#' other variable types. The value of the reference category is also recoded
#' when bgm recodes the corresponding observations. Only required if there is at
#' least one variable of type ``blume-capel''.
#' @param iter How many iterations should the Gibbs sampler run? The default of
#' \code{1e4} is for illustrative purposes. For stable estimates, it is
#' recommended to run the Gibbs sampler for at least \code{1e5} iterations.
#' @param burnin The number of iterations of the Gibbs sampler before saving its
#' output. Since it may take some time for the Gibbs sampler to converge to
#' the posterior distribution, it is recommended not to set this number too low.
#' @param interaction_scale The scale of the Cauchy distribution that is used as a
#' prior for the pairwise interaction parameters. Defaults to \code{2.5}.
#' @param difference_scale The scale of the Cauchy distribution that is used as the
#' prior for the pairwise difference parameters. Defaults to \code{0.1}.
#' @param threshold_alpha,threshold_beta The shape parameters of the beta-prime
#' prior density for the threshold parameters. Must be positive values. If the
#' two values are equal, the prior density is symmetric about zero. If
#' \code{threshold_beta} is greater than \code{threshold_alpha}, the
#' distribution is skewed to the left, and if \code{threshold_beta} is less than
#' \code{threshold_alpha}, it is skewed to the right. Smaller values tend to
#' lead to more diffuse prior distributions.
#' @param difference_prior The inclusion or exclusion of pairwise differences in
#' the pairwise interactions between the two groups are modeled with binary
#' indicator variables. The argument \code{difference_prior} is used to set a
#' prior distribution for these indicator variables, i.e., the structure of the
#' network. Currently, two options are implemented: The Bernoulli model
#' \code{difference_prior = "Bernoulli"} assumes that the probability that there
#' is a pairwise difference is equal to \code{difference_probability}
#' and independent of other pairwise differences. When
#' \code{difference_probability = 0.5}, this means that each possible
#' configuration of the presence or absence of pairwise differences is given the
#' same prior weight. The Beta-Bernoulli model
#' \code{difference_prior = "Beta-Bernoulli"} assumes a beta prior for the unknown
#' difference probability with shape parameters \code{beta_bernoulli_alpha} and
#' \code{beta_bernoulli_beta}. If \code{beta_bernoulli_alpha = 1} and
#' \code{beta_bernoulli_beta = 1}, this means that different configurations of
#' present or absent diffeences that have the same complexity (number of present
#' differences) receive the same prior weight. The default is
#' \code{difference_prior = "Bernoulli"}.
#' @param difference_probability The prior probability for a pairwise difference
#' in the Bernoulli model. Can be a single probability, or a matrix of \code{p} rows
#' and \code{p} columns specifying the probability of a difference for each edge
#' pair. The default is \code{difference_probability = 0.5}.
#' @param beta_bernoulli_alpha,beta_bernoulli_beta The two shape parameters of
#' the Beta prior density for the Bernoulli inclusion probability. Must be
#' positive numbers. Defaults to \code{beta_bernoulli_alpha = 1} and
#' \code{beta_bernoulli_beta = 1}.
#' @param na.action How do you want the function to handle missing data? If
#' \code{na.action = "listwise"}, listwise deletion is used. If
#' \code{na.action = "impute"}, missing data are imputed iteratively during the
#' MCMC procedure. Since imputation of missing data can have a negative impact
#' on the convergence speed of the MCMC procedure, it is recommended to run the
#' MCMC for more iterations. Also, since the numerical routines that search for
#' the mode of the posterior do not have an imputation option, the bgm function
#' will automatically switch to \code{interaction_prior = "Cauchy"} and
#' \code{adaptive = TRUE}.
#' @param save Should the function collect and return all samples from the Gibbs
#' sampler (\code{save = TRUE})? Or should it only return the (model-averaged)
#' posterior means (\code{save = FALSE})? Defaults to \code{FALSE}.
#' @param display_progress Should the function show a progress bar
#' (\code{display_progress = TRUE})? Or not (\code{display_progress = FALSE})?
#' The default is \code{TRUE}.
#'
#' @return If \code{save = FALSE} (the default), the result is a list of class
#' ``bgmNCT'' containing the following matrices:
#'  \itemize{
#'    \item \code{gamma}: A matrix with \code{p} rows and \code{p} columns,
#'    containing posterior inclusion probabilities of pairwise differences.
#'    \item \code{interactions}: A matrix with \code{p} rows and \code{p} columns,
#'    containing model-averaged posterior means of the pairwise associations.
#'    \item \code{interactions_difference}: A matrix with \code{p} rows and \code{p}
#'    columns,containing model-averaged posterior means of the pairwise differences.
#'    \item \code{thresholds}: A matrix with \code{p} rows and \code{max(m)}
#'    columns, containing model-averaged category thresholds. In the case of
#'    ``blume-capel'' variables, the first entry is the parameter for the linear
#'    effect and the second entry is the parameter for the quadratic effect, which
#'    models the offset to the reference category.
#'  }
#'
#' If \code{save = TRUE}, the result is a list of class ``bgmNCT'' containing:
#' \itemize{
#'  \item \code{gamma}: A matrix with \code{iter} rows and
#'   \code{p * (p - 1) / 2} columns, containing the inclusion indicators of the
#'    pairwise difference from every iteration of the Gibbs sampler.
#'  \item \code{interactions}: A matrix with \code{iter} rows and
#'    \code{p * (p - 1) / 2} columns, containing parameter states from every
#'    iteration of the Gibbs sampler for the pairwise associations.
#'  \item \code{interactions_difference}: A matrix with \code{iter} rows and
#'    \code{p * (p - 1) / 2} columns, containing parameter states from every
#'    iteration of the Gibbs sampler for the pairwise differences.
#'  \item \code{thresholds}: A matrix with \code{iter} rows and
#'    \code{sum(m)} columns, containing parameter states from every iteration of
#'    the Gibbs sampler for the category thresholds.
#'  }
#' Column averages of these matrices provide the model-averaged posterior means.
#'
#' In addition to the analysis results, the bgmNCT output lists some of the
#' arguments of its call. This is useful for post-processing the results.
#'
#' @importFrom utils packageVersion
#' @export
bgmNCT = function(x,
                  variable_type = "ordinal",
                  group_indicator,
                  reference_category,
                  iter = 1e4,
                  burnin = 1e3,
                  interaction_scale = 2.5,
                  difference_scale = 0.1,
                  threshold_alpha = 0.5,
                  threshold_beta = 0.5,
                  difference_prior = c("Bernoulli", "Beta-Bernoulli"),
                  difference_probability = 0.5,
                  beta_bernoulli_alpha = 1,
                  beta_bernoulli_beta = 1,
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

  #Check model input -----------------------------------------------------------
  model = nct_check_model(x = x,
                          variable_type = variable_type,
                          reference_category = reference_category,
                          interaction_scale = interaction_scale,
                          difference_scale = difference_scale,
                          threshold_alpha = threshold_alpha,
                          threshold_beta = threshold_beta,
                          difference_prior = difference_prior,
                          difference_probability = difference_probability,
                          beta_bernoulli_alpha = beta_bernoulli_alpha,
                          beta_bernoulli_beta = beta_bernoulli_beta)

  # ----------------------------------------------------------------------------
  # The vector variable_type is now coded as boolean.
  # Ordinal (variable_bool == TRUE) or Blume-Capel (variable_bool == FALSE)
  # ----------------------------------------------------------------------------
  variable_bool = model$variable_bool
  # ----------------------------------------------------------------------------

  reference_category = model$reference_category
  difference_prior = model$difference_prior
  theta = model$theta

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
  data = nct_reformat_data(x = x,
                           na.action = na.action,
                           variable_bool = variable_bool,
                           reference_category = reference_category,
                           group_indicator = group_indicator)
  x = data$x
  no_categories = data$no_categories
  missing_index = data$missing_index
  na.impute = data$na.impute
  reference_category = data$reference_category
  group_indicator = data$group_indicator

  no_variables = ncol(x)
  no_interactions = no_variables * (no_variables - 1) / 2
  no_thresholds = sum(no_categories)

  #Specify the variance of the (normal) proposal distribution ------------------
  proposal_sd = matrix(1,
                       nrow = no_variables,
                       ncol = no_variables)
  proposal_sd_blumecapel = matrix(1,
                                  nrow = no_variables,
                                  ncol = 2)
  proposal_sd_int_diff = matrix(1,
                                nrow = no_variables,
                                ncol = no_variables)

  # Starting value of model matrix ---------------------------------------------
  gamma = matrix(1,
                 nrow = no_variables,
                 ncol = no_variables)


  #Starting values of interactions and thresholds (posterior mode) -------------
  interactions = matrix(0, nrow = no_variables, ncol = no_variables)
  interactions_difference = matrix(0, nrow = no_variables, ncol = no_variables)
  thresholds = matrix(0, nrow = no_variables, ncol = max(no_categories))

  #Precompute the number of observations per category for each variable --------
  n_cat_obs = matrix(0,
                     nrow = max(no_categories) + 1,
                     ncol = no_variables)
  for(variable in 1:no_variables) {
    for(category in 0:no_categories[variable]) {
      n_cat_obs[category + 1, variable] = sum(x[, variable] == category)
    }
  }

  #Precompute the sufficient statistics for the two Blume-Capel parameters -----
  sufficient_blume_capel = matrix(0, nrow = 2, ncol = no_variables)
  if(any(!variable_bool)) {
    # Ordinal (variable_bool == TRUE) or Blume-Capel (variable_bool == FALSE)
    bc_vars = which(!variable_bool)
    for(i in bc_vars) {
      sufficient_blume_capel[1, i] = sum(x[, i])
      sufficient_blume_capel[2, i] = sum((x[, i] - reference_category[i]) ^ 2)
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
  out = nct_gibbs_sampler(observations = x,
                          gamma = gamma,
                          interactions = interactions,
                          thresholds = thresholds,
                          interactions_difference = interactions_difference,
                          group_indicator = group_indicator,
                          no_categories  = no_categories,
                          interaction_scale = interaction_scale,
                          difference_scale = difference_scale,
                          proposal_sd = proposal_sd,
                          proposal_sd_blumecapel = proposal_sd_blumecapel,
                          proposal_sd_int_diff = proposal_sd_int_diff,
                          difference_prior = difference_prior,
                          theta = theta,
                          beta_bernoulli_alpha = beta_bernoulli_alpha,
                          beta_bernoulli_beta = beta_bernoulli_beta,
                          Index = Index,
                          iter = iter,
                          burnin = burnin,
                          n_cat_obs = n_cat_obs,
                          sufficient_blume_capel = sufficient_blume_capel,
                          threshold_alpha = threshold_alpha,
                          threshold_beta = threshold_beta,
                          na_impute = na.impute,
                          missing_index = missing_index,
                          variable_bool = variable_bool,
                          reference_category = reference_category,
                          save = save,
                          display_progress = display_progress)

  #Preparing the output --------------------------------------------------------
  arguments = list(
    no_variables = no_variables,
    no_cases = nrow(x),
    na_impute = na.impute,
    variable_type = variable_type,
    iter = iter,
    burnin = burnin,
    interaction_scale = interaction_scale,
    threshold_alpha = threshold_alpha,
    threshold_beta = threshold_beta,
    difference_prior = difference_prior,
    difference_probability = theta,
    beta_bernoulli_alpha = beta_bernoulli_alpha ,
    beta_bernoulli_beta =  beta_bernoulli_beta,
    na.action = na.action,
    save = save,
    group_indicator = group_indicator,
    version = packageVersion("bgms")
  )

  if(save == FALSE) {
    gamma = out$gamma
    interactions = out$interactions
    interactions_difference = out$interactions_difference
    tresholds = out$thresholds

    if(is.null(colnames(x))){
      data_columnnames = paste0("variable ", 1:no_variables)
      colnames(interactions) = data_columnnames
      rownames(interactions) = data_columnnames
      colnames(interactions_difference) = data_columnnames
      rownames(interactions_difference) = data_columnnames

      colnames(gamma) = data_columnnames
      rownames(gamma) = data_columnnames
      rownames(thresholds) = data_columnnames
    } else {
      data_columnnames <- colnames(x)
      colnames(interactions) = data_columnnames
      rownames(interactions) = data_columnnames
      colnames(interactions_difference) = data_columnnames
      rownames(interactions_difference) = data_columnnames

      colnames(gamma) = data_columnnames
      rownames(gamma) = data_columnnames
      rownames(thresholds) = data_columnnames
    }

    colnames(tresholds) = paste0("category ", 1:max(no_categories))

    arguments$data_columnnames = data_columnnames

    output = list(gamma = gamma,
                  interactions = interactions,
                  interactions_difference = interactions_difference,
                  thresholds = thresholds,
                  arguments = arguments)

    class(output) = c("bgmNCT")
    return(output)
  } else {
    gamma = out$gamma
    interactions = out$interactions
    interactions_difference = out$interactions_difference
    thresholds = out$thresholds

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

    colnames(gamma) = colnames(interactions) = colnames(interactions_difference) = names_vec
    names = character(length = sum(no_categories))
    cntr = 0
    for(variable in 1:no_variables) {
      for(category in 1:no_categories[variable]) {
        cntr = cntr + 1
        names[cntr] = paste0("threshold(",variable, ", ",category,")")
      }
    }
    colnames(thresholds) = names

    dimnames(gamma) = list(Iter. = 1:iter, colnames(gamma))
    dimnames(interactions) = list(Iter. = 1:iter, colnames(interactions))
    dimnames(interactions_difference) = list(Iter. = 1:iter, colnames(interactions_difference))
    dimnames(thresholds) = list(Iter. = 1:iter, colnames(thresholds))

    arguments$data_columnnames = data_columnnames

    output = list(gamma = gamma,
                  interactions = interactions,
                  interactions_difference = interactions_difference,
                  thresholds = thresholds,
                  arguments = arguments)
    class(output) = "bgmsNCT"
    return(output)
  }
}
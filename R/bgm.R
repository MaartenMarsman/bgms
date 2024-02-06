#' Bayesian edge selection or Bayesian estimation for Markov Random Fields of
#' mixed binary and ordinal variables using MCMC.
#'
#' The function \code{bgm} explores the joint pseudoposterior distribution of
#' parameters and possibly edge indicators for a Markov Random Field model for
#' mixed binary and ordinal variables.
#'
#' Currently, bgm supports two types of ordinal variables. The regular, default,
#' ordinal variable type has no restrictions on its distribution. Every response
#' category except the first receives its own threshold parameter. The
#' Blume-Capel ordinal variable assumes that there is a specific reference
#' category, such as the ``neutral'' in a Likert scale, and responses are scored
#' in terms of their distance to this reference category. Specifically, the
#' Blume-Capel model specifies the following quadratic model for the threshold
#' parameters:
#' \deqn{\mu_{\text{c}} = \alpha \times \text{c} + \beta \times (\text{c} - \text{r})^2,}{{\mu_{\text{c}} = \alpha \times \text{c} + \beta \times (\text{c} - \text{r})^2,}}
#' where \eqn{\mu_{\text{c}}}{\mu_{\text{c}}} is the threshold for category c.
#' The parameter \eqn{\alpha}{\alpha} models a linear trend across categories,
#' such that \eqn{\alpha > 0}{\alpha > 0} leads to an increasing number of
#' observations in higher response categories and \eqn{\alpha <0}{\alpha <0}
#' leads to a decreasing number of observations in higher response categories.
#' The parameter \eqn{\beta}{\beta} models the response style in terms of an
#' offset with respect to the reference category \eqn{r}{r}; if \eqn{\beta<0}{\beta<0}
#' there is a preference to respond in the reference category (i.e., the model
#' introduces a penalty for responding in a category further away from the
#' reference_category category \code{r}), while if \eqn{\beta > 0}{\beta > 0}
#' there is preference to score in the extreme categories further away from the
#' reference_category category.
#'
#' The Bayesian estimation procedure (\code{edge_selection = FALSE}) simply
#' estimates the threshold and pairwise interaction parameters of the ordinal
#' MRF, while the Bayesian edge selection procedure
#' (\code{edge_selection = TRUE}) also models the probability that individual
#' edges should be included or excluded from the model. Bayesian edge selection
#' imposes a discrete spike and slab prior distribution on the pairwise
#' interactions. By formulating it as a mixture of mutually singular
#' distributions, the function can use a combination of Metropolis-Hastings and
#' Gibbs sampling to create a Markov chain that has the joint posterior
#' distribution as an invariant. The current option for the slab distribution is
#' a Cauchy with an optional scaling parameter. If there are no missing data and
#' no Blume-Capel variables, there is also an option to use a unit-information
#' prior for the slab distribution instead. The slab distribution is also used
#' as the prior for the interaction parameters for Bayesian estimation. A
#' beta-prime distribution is used for the exponent of the category parameters.
#' For Bayesian edge selection, two prior distributions are implemented for the
#' edge inclusion variables (i.e., the prior probability that an edge is
#' included); the Bernoulli prior and the Beta-Bernoulli prior.
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
#' @param interaction_prior The prior distribution to use for the pairwise
#' interaction parameters. If \code{edge_selection = TRUE} this prior is the
#' slab distribution. The current option for this prior is a Cauchy distribution
#' (\code{interaction_prior = "Cauchy"}). If there are no missing data and no
#' Blume-Capel variables, there is also an option to use a unit information
#' prior instead (\code{interaction_prior = "UnitInfo"}). The default is
#' (\code{interaction_prior = "Cauchy"}).
#' @param cauchy_scale The scale of the Cauchy distribution that is used as a
#' prior for the pairwise interaction parameters. Defaults to \code{2.5}.
#' @param threshold_alpha,threshold_beta The shape parameters of the beta-prime
#' prior density for the threshold parameters. Must be positive values. If the
#' two values are equal, the prior density is symmetric about zero. If
#' \code{threshold_beta} is greater than \code{threshold_alpha}, the
#' distribution is skewed to the left, and if \code{threshold_beta} is less than
#' \code{threshold_alpha}, it is skewed to the right. Smaller values tend to
#' lead to more diffuse prior distributions.
#' @param edge_selection Should the function perform Bayesian edge selection on
#' the edges of the MRF in addition to estimating its parameters
#' (\code{edge_selection = TRUE}), or should it just estimate the parameters
#' (\code{edge_selection = FALSE})? The default is \code{edge_selection = TRUE}.
#' @param edge_prior The inclusion or exclusion of individual edges in the
#' network is modeled with binary indicator variables that capture the structure
#' of the network. The argument \code{edge_prior} is used to set a prior
#' distribution for the edge indicator variables, i.e., the structure of the
#' network. Currently, two options are implemented: The Bernoulli model
#' \code{edge_prior = "Bernoulli"} assumes that the probability that an edge
#' between two variables is included is equal to \code{inclusion_probability}
#' and independent of other edges or variables. When
#' \code{inclusion_probability = 0.5}, this means that each possible network
#' structure is given the same prior weight. The Beta-Bernoulli model
#' \code{edge_prior = "Beta-Bernoulli"} assumes a beta prior for the unknown
#' inclusion probability with shape parameters \code{beta_bernoulli_alpha} and
#' \code{beta_bernoulli_beta}. If \code{beta_bernoulli_alpha = 1} and
#' \code{beta_bernoulli_beta = 1}, this means that networks with the same
#' complexity (number of edges) get the same prior weight. The default is
#' \code{edge_prior = "Bernoulli"}.
#' @param inclusion_probability The prior edge inclusion probability for the
#' Bernoulli model. Can be a single probability, or a matrix of \code{p} rows
#' and \code{p} columns specifying an inclusion probability for each edge pair.
#' The default is \code{inclusion_probability = 0.5}.
#' @param beta_bernoulli_alpha,beta_bernoulli_beta The two shape parameters of
#' the Beta prior density for the Bernoulli inclusion probability. Must be
#' positive numbers. Defaults to \code{beta_bernoulli_alpha = 1} and
#' \code{beta_bernoulli_beta = 1}.
#' @param adaptive A random walk Metropolis algorithm is used to sample from the
#' fully conditional posterior distributions of the pairwise interaction
#' parameters. This requires a variance to be specified for the prior
#' distribution. If there are no missing values and no ``blume-capel'' variables,
#' there is an option to set this variance equal to the curvature around the
#' posterior mode. This is the default of the bgm function, but requires the
#' second derivative of the pseudoposterior at its mode, which bgm cannot
#' determine in the case of missing data, and is not implemented for
#' ``blume-capel'' variables. In other cases, bgm switches to an adaptive
#' Metropolis algorithm, which adjusts the proposal variance to the acceptance
#' probability of the random walk Metropolis algorithm to be close to the
#' optimum of \code{.234} using a Robbins-Monro-type algorithm. The user can
#' also select the adaptive Metropolis algorithm by default
#' (\code{adaptive = TRUE}).
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
#' ``bgms'' containing the following matrices:
#' \itemize{
#' \item \code{gamma}: A matrix with \code{p} rows and \code{p} columns,
#' containing posterior inclusion probabilities of individual edges.
#' \item \code{interactions}: A matrix with \code{p} rows and \code{p} columns,
#' containing model-averaged posterior means of the pairwise associations.
#' \item \code{thresholds}: A matrix with \code{p} rows and \code{max(m)}
#' columns, containing model-averaged category thresholds. In the case of
#' ``blume-capel'' variables, the first entry is the parameter for the linear
#' effect and the second entry is the parameter for the quadratic effect, which
#' models the offset to the reference category.
#' }
#'
#' If \code{save = TRUE}, the result is a list of class ``bgms'' containing:
#' \itemize{
#' \item \code{gamma}: A matrix with \code{iter} rows and
#' \code{p * (p - 1) / 2} columns, containing the edge inclusion indicators from
#' every iteration of the Gibbs sampler.
#' \item \code{interactions}: A matrix with \code{iter} rows and
#' \code{p * (p - 1) / 2} columns, containing parameter states from every
#' iteration of the Gibbs sampler for the pairwise associations.
#' \item \code{thresholds}: A matrix with \code{iter} rows and
#' \code{sum(m)} columns, containing parameter states from every iteration of
#' the Gibbs sampler for the category thresholds.
#' }
#' Column averages of these matrices provide the model-averaged posterior means.
#'
#' In addition to the analysis results, the bgm output lists some of the
#' arguments of its call. This is useful for post-processing the results.
#'
#' @examples
#' \donttest{
#'  #Store user par() settings
#'  op <- par(no.readonly = TRUE)
#'
#'  ##Analyse the Wenchuan dataset
#'
#'  # Here, we use 1e4 iterations, for an actual analysis please use at least
#'  # 1e5 iterations.
#'  fit = bgm(x = Wenchuan)
#'
#'
#'  #------------------------------------------------------------------------------|
#'  # INCLUSION - EDGE WEIGHT PLOT
#'  #------------------------------------------------------------------------------|
#'
#'  par(mar = c(6, 5, 1, 1))
#'  plot(x = fit$interactions[lower.tri(fit$interactions)],
#'       y = fit$gamma[lower.tri(fit$gamma)], ylim = c(0, 1),
#'       xlab = "", ylab = "", axes = FALSE, pch = 21, bg = "gray", cex = 1.3)
#'  abline(h = 0, lty = 2, col = "gray")
#'  abline(h = 1, lty = 2, col = "gray")
#'  abline(h = .5, lty = 2, col = "gray")
#'  mtext("Posterior Mode Edge Weight", side = 1, line = 3, cex = 1.7)
#'  mtext("Posterior Inclusion Probability", side = 2, line = 3, cex = 1.7)
#'  axis(1)
#'  axis(2, las = 1)
#'
#'
#'  #------------------------------------------------------------------------------|
#'  # EVIDENCE - EDGE WEIGHT PLOT
#'  #------------------------------------------------------------------------------|
#'
#'  #For the default choice of the structure prior, the prior odds equal one:
#'  prior.odds = 1
#'  posterior.inclusion = fit$gamma[lower.tri(fit$gamma)]
#'  posterior.odds = posterior.inclusion / (1 - posterior.inclusion)
#'  log.bayesfactor = log(posterior.odds / prior.odds)
#'  log.bayesfactor[log.bayesfactor > 5] = 5
#'
#'  par(mar = c(5, 5, 1, 1) + 0.1)
#'  plot(fit$interactions[lower.tri(fit$interactions)], log.bayesfactor, pch = 21, bg = "#bfbfbf",
#'       cex = 1.3, axes = FALSE, xlab = "", ylab = "", ylim = c(-5, 5.5),
#'       xlim = c(-0.5, 1.5))
#'  axis(1)
#'  axis(2, las = 1)
#'  abline(h = log(1/10), lwd = 2, col = "#bfbfbf")
#'  abline(h = log(10), lwd = 2, col = "#bfbfbf")
#'
#'  text(x = 1, y = log(1 / 10), labels = "Evidence for Exclusion", pos = 1,
#'       cex = 1.7)
#'  text(x = 1, y = log(10), labels = "Evidence for Inclusion", pos = 3, cex = 1.7)
#'  text(x = 1, y = 0, labels = "Absence of Evidence", cex = 1.7)
#'  mtext("Log-Inclusion Bayes Factor", side = 2, line = 3, cex = 1.5, las = 0)
#'  mtext("Posterior Mean Interactions ", side = 1, line = 3.7, cex = 1.5, las = 0)
#'
#'
#'  #------------------------------------------------------------------------------|
#'  # THE LOCAL MEDIAN PROBABILITY NETWORK
#'  #------------------------------------------------------------------------------|
#'
#'  tmp = fit$interactions[lower.tri(fit$interactions)]
#'  tmp[posterior.inclusion < 0.5] = 0
#'
#'  median.prob.model = matrix(0, nrow = ncol(Wenchuan), ncol = ncol(Wenchuan))
#'  median.prob.model[lower.tri(median.prob.model)] = tmp
#'  median.prob.model = median.prob.model + t(median.prob.model)
#'
#'  rownames(median.prob.model) = colnames(Wenchuan)
#'  colnames(median.prob.model) = colnames(Wenchuan)
#'
#'  library(qgraph)
#'  qgraph(median.prob.model,
#'         theme = "TeamFortress",
#'         maximum = .5,
#'         fade = FALSE,
#'         color = c("#f0ae0e"), vsize = 10, repulsion = .9,
#'         label.cex = 1.1, label.scale = "FALSE",
#'         labels = colnames(Wenchuan))
#'
#'  #Restore user par() settings
#'  par(op)
#' }
#' @export
bgm = function(x,
               variable_type = "ordinal",
               reference_category,
               iter = 1e4,
               burnin = 1e3,
               interaction_prior = c("Cauchy", "UnitInfo"),
               cauchy_scale = 2.5,
               threshold_alpha = 0.5,
               threshold_beta = 0.5,
               edge_selection = TRUE,
               edge_prior = c("Bernoulli", "Beta-Bernoulli"),
               inclusion_probability = 0.5,
               beta_bernoulli_alpha = 1,
               beta_bernoulli_beta = 1,
               adaptive = FALSE,
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
  model = check_bgm_model(x = x,
                          variable_type = variable_type,
                          reference_category = reference_category,
                          interaction_prior = interaction_prior,
                          cauchy_scale = cauchy_scale,
                          threshold_alpha = threshold_alpha,
                          threshold_beta = threshold_beta,
                          edge_selection = edge_selection,
                          edge_prior = edge_prior,
                          inclusion_probability = inclusion_probability,
                          beta_bernoulli_alpha = beta_bernoulli_alpha,
                          beta_bernoulli_beta = beta_bernoulli_beta,
                          adaptive = adaptive)

  variable_type = model$variable_type
  reference_category = model$reference_category
  interaction_prior = model$interaction_prior
  edge_prior = model$edge_prior
  adaptive = model$adaptive
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
  na.action = match.arg(na.action)

  #Format the data input -------------------------------------------------------
  data = reformat_data_bgm(x = x,
                           na.action = na.action,
                           variable_type = variable_type,
                           reference_category = reference_category)
  x = data$x
  no_categories = data$no_categories
  missing_index = data$missing_index
  na.impute = data$na.impute
  reference_category = data$reference_category

  if(na.impute == TRUE) {
    if(interaction_prior != "Cauchy")
      warning(paste0("There were missing responses and na.action was set to ``impute''. The bgm \n",
                     "function must switch the interaction_prior to ``Cauchy''."))
    adaptive = TRUE
    interaction_prior = "Cauchy"
    if(cauchy_scale <= 0 || is.na(cauchy_scale) || is.infinite(cauchy_scale))
      stop("The scale of the Cauchy prior needs to be positive.")
  }

  no_variables = ncol(x)
  no_interactions = no_variables * (no_variables - 1) / 2
  no_thresholds = sum(no_categories)

  #Proposal set-up for the interaction parameters ------------------------------
  if(interaction_prior == "UnitInfo") {
    pps = try(mppe(x = x,
                   interaction_prior = interaction_prior),
              silent = TRUE)
    if(inherits(pps, what = "try-error"))
      stop(paste0(
        "For the Unit Information prior we need to estimate the posterior mode. \n",
        "Unfortunately, bgm could not find this mode for your data. Please try the \n",
        "Cauchy prior option."))
    unit_info = sqrt(pps$unit_info)
  } else {
    if(!na.impute && !any(variable_type == "blume-capel")) {
      pps = try(mppe(x = x,
                     interaction_prior = interaction_prior,
                     cauchy_scale = cauchy_scale),
                silent = TRUE)
      if(inherits(pps, what = "try-error") & adaptive == FALSE) {
        stop(paste0("By default, the MCMC procedure underlying the bgm function uses a Metropolis \n",
                    "algorithm with a fixed proposal distribution. The bgm function attempts to fit \n",
                    "this proposal distribution to the target posterior distribution by locating the \n",
                    "posterior mode and using information about the curvature around that model to \n",
                    "set the variance of the proposal distributions. Unfortunately, bgm was unable \n",
                    "to locate the posterior mode for your data. Please try again with ``adaptive = \n",
                    "TRUE''."))
      }
    }
    unit_info = matrix(data = NA, nrow = 1, ncol = 1)
  }

  #Specify the variance of the (normal) proposal distribution ------------------
  proposal_sd = matrix(1,
                       nrow = no_variables,
                       ncol = no_variables)
  proposal_sd_blumecapel = matrix(1,
                                 nrow = no_variables,
                                 ncol = 2)
  if(adaptive == FALSE && !na.impute) {
    hessian = pps$hessian[-c(1:no_thresholds), -c(1:no_thresholds)]
    cntr = 0
    if(no_variables == 2) {
      proposal_sd[1, 2] = sqrt(-1 / hessian[cntr])
      proposal_sd[2, 1] = proposal_sd[1, 2]
    } else {
      for(variable1 in 1:(no_variables - 1)) {
        for(variable2 in (variable1 + 1):no_variables) {
          cntr = cntr + 1
          proposal_sd[variable1, variable2] = sqrt(-1 / hessian[cntr, cntr])
          proposal_sd[variable2, variable1] = proposal_sd[variable1, variable2]
        }
      }
    }
  }

  # Starting value of model matrix ---------------------------------------------
  gamma = matrix(1,
                 nrow = no_variables,
                 ncol = no_variables)


  #Starting values of interactions and thresholds (posterior mode) -------------
  if(!na.impute && !any(variable_type == "blume-capel") && !inherits(pps, what = "try-error")) {
    interactions = pps$interactions
    thresholds = pps$thresholds
  } else {
    interactions = matrix(0, nrow = no_variables, ncol = no_variables)
    thresholds = matrix(0, nrow = no_variables, ncol = max(no_categories))
  }

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
  sufficient_statistics_blume_capel = matrix(0, nrow = 2, ncol = no_variables)
  if(any(variable_type == "blume-capel")) {
    bc_vars = which(variable_type == "blume-capel")
    for(i in bc_vars) {
      sufficient_statistics_blume_capel[1, i] = sum(x[, i])
      sufficient_statistics_blume_capel[2, i] = sum((x[, i] - reference_category[i]) ^ 2)
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
  out = gibbs_sampler(observations = x,
                      gamma = gamma,
                      interactions = interactions,
                      thresholds = thresholds,
                      no_categories  = no_categories,
                      interaction_prior = interaction_prior,
                      cauchy_scale = cauchy_scale,
                      unit_info = unit_info,
                      proposal_sd = proposal_sd,
                      proposal_sd_blumecapel = proposal_sd_blumecapel,
                      edge_prior = edge_prior,
                      theta = theta,
                      beta_bernoulli_alpha = beta_bernoulli_alpha,
                      beta_bernoulli_beta = beta_bernoulli_beta,
                      Index = Index,
                      iter = iter,
                      burnin = burnin,
                      n_cat_obs = n_cat_obs,
                      sufficient_statistics_blume_capel = sufficient_statistics_blume_capel,
                      threshold_alpha = threshold_alpha,
                      threshold_beta = threshold_beta,
                      na_impute = na.impute,
                      missing_index = missing_index,
                      variable_type = variable_type,
                      reference_category = reference_category,
                      adaptive = adaptive,
                      save = save,
                      display_progress = display_progress,
                      edge_selection = edge_selection)


  #Preparing the output --------------------------------------------------------
  bgm_arguments = list(
    iter = iter,
    burnin = burnin,
    interaction_prior = interaction_prior,
    cauchy_scale = cauchy_scale,
    threshold_alpha = threshold_alpha,
    threshold_beta = threshold_beta,
    edge_selection = edge_selection,
    edge_prior = edge_prior,
    inclusion_probability = inclusion_probability,
    beta_bernoulli_alpha = beta_bernoulli_alpha ,
    beta_bernoulli_beta =  beta_bernoulli_beta,
    adaptive = adaptive,
    na.action = na.action,
    save = save
  )

  if(save == FALSE) {
    if(edge_selection == TRUE) {
      gamma = out$gamma
    }
    interactions = out$interactions
    tresholds = out$thresholds

    if(is.null(colnames(x))){
      data_columnnames = paste0("variable ", 1:no_variables)
      colnames(interactions) = data_columnnames
      rownames(interactions) = data_columnnames
      if(edge_selection == TRUE) {
        colnames(gamma) = data_columnnames
        rownames(gamma) = data_columnnames
      }
      rownames(thresholds) = data_columnnames
    } else {
      data_columnnames <- colnames(x)
      colnames(interactions) = data_columnnames
      rownames(interactions) = data_columnnames
      if(edge_selection == TRUE) {
        colnames(gamma) = data_columnnames
        rownames(gamma) = data_columnnames
      }
      rownames(thresholds) = data_columnnames
    }

    colnames(tresholds) = paste0("category ", 1:max(no_categories))

    if(edge_selection == TRUE) {
      output = list(gamma = gamma,
                    interactions = interactions,
                    thresholds = thresholds,
                    bgm_arguments = bgm_arguments,
                    colnames = data_columnnames)
    } else {
      output = list(interactions = interactions,
                    thresholds = thresholds,
                    bgm_arguments = bgm_arguments,
                    colnames = data_columnnames)
    }

    class(output) = "bgms"
    return(output)
  } else {
    if(edge_selection == TRUE) {
      gamma = out$gamma
    }
    interactions = out$interactions
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

    if(edge_selection == TRUE) {
      colnames(gamma) = colnames(interactions) = names_vec
    }
    names = character(length = sum(no_categories))
    cntr = 0
    for(variable in 1:no_variables) {
      for(category in 1:no_categories[variable]) {
        cntr = cntr + 1
        names[cntr] = paste0("threshold(",variable, ", ",category,")")
      }
    }
    colnames(thresholds) = names

    if(edge_selection == TRUE) {
      dimnames(gamma) = list(Iter. = 1:iter, colnames(gamma))
    }
    dimnames(interactions) = list(Iter. = 1:iter, colnames(interactions))
    dimnames(thresholds) = list(Iter. = 1:iter, colnames(thresholds))

    if(edge_selection == TRUE) {
      output = list(gamma = gamma,
                    interactions = interactions,
                    thresholds = thresholds,
                    bgm_arguments = bgm_arguments,
                    colnames = data_columnnames)
    } else {
      output = list(interactions = interactions,
                    thresholds = thresholds,
                    bgm_arguments = bgm_arguments,
                    colnames = data_columnnames)
    }
    class(output) = "bgms"
    return(output)
  }
}


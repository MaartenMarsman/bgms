#' Bayesian edge selection or Bayesian estimation for a Markov random field
#' model for binary and/or ordinal variables.
#'
#' The function \code{bgm} explores the joint pseudoposterior distribution of
#' parameters and possibly edge indicators for a Markov Random Field model for
#' mixed binary and ordinal variables.
#'
#' Currently, \code{bgm} supports two types of ordinal variables. The regular, default,
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
#' a Cauchy with an optional scaling parameter. The slab distribution is also used
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
#' output. Since it may take some time for the Gibbs sampler to converge to the
#' posterior distribution, it is recommended not to set this number too low.
#' When \code{edge_selection = TRUE}, the bgm function will perform
#' \code{2 * burnin} iterations, first \code{burnin} iterations without edge
#' selection, then \code{burnin} iterations with edge selection. This helps
#' ensure that the Markov chain used for estimation starts with good parameter
#' values and that the adaptive MH proposals are properly calibrated.
#' @param interaction_scale The scale of the Cauchy distribution that is used as
#' a prior for the pairwise interaction parameters. Defaults to \code{2.5}.
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
#' network. Currently, three options are implemented: The Bernoulli model
#' \code{edge_prior = "Bernoulli"} assumes that the probability that an edge
#' between two variables is included is equal to \code{inclusion_probability}
#' and independent of other edges or variables. When
#' \code{inclusion_probability = 0.5}, this means that each possible network
#' structure is given the same prior weight. The Beta-Bernoulli model
#' \code{edge_prior = "Beta-Bernoulli"} assumes a beta prior for the unknown
#' inclusion probability with shape parameters \code{beta_bernoulli_alpha} and
#' \code{beta_bernoulli_beta}. If \code{beta_bernoulli_alpha = 1} and
#' \code{beta_bernoulli_beta = 1}, this means that networks with the same
#' complexity (number of edges) get the same prior weight. The Stochastic Block
#' model \code{edge_prior = "Stochastic-Block"} assumes that nodes can be
#' organized into blocks or clusters. In principle, the assignment of nodes to
#' such clusters is unknown, and the model as implemented here considers all
#' possible options \insertCite{@i.e., specifies a Dirichlet prior on the probability
#' of allocations as described by @GengEtAl_2019}{bgms}. This model is advantageous
#' when nodes are expected to fall into distinct clusters. The inclusion
#' probabilities for the edges are defined at the level of the clusters, with a
#' beta prior for the unknown inclusion probability with shape parameters
#' \code{beta_bernoulli_alpha} and \code{beta_bernoulli_beta}, and a Dirichlet
#' prior on the cluster assignment probabilities with a common concentration
#' parameter \code{dirichlet_alpha}. The default is
#' \code{edge_prior = "Bernoulli"}.
#' @param inclusion_probability The prior edge inclusion probability for the
#' Bernoulli model. Can be a single probability, or a matrix of \code{p} rows
#' and \code{p} columns specifying an inclusion probability for each edge pair.
#' The default is \code{inclusion_probability = 0.5}.
#' @param beta_bernoulli_alpha,beta_bernoulli_beta The two shape parameters of
#' the Beta prior density for the Bernoulli inclusion probability. Must be
#' positive numbers. Defaults to \code{beta_bernoulli_alpha = 1} and
#' \code{beta_bernoulli_beta = 1}.
#' @param dirichlet_alpha The shape of the Dirichlet prior on the node-to-block
#' allocation probabilities for the Stochastic Block model.
#' @param na_action How do you want the function to handle missing data? If
#' \code{na_action = "listwise"}, listwise deletion is used. If
#' \code{na_action = "impute"}, missing data are imputed iteratively during the
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
#' ``bgms'' containing the following matrices with model-averaged quantities:
#' \itemize{
#' \item \code{indicator}: A matrix with \code{p} rows and \code{p} columns,
#' containing the posterior inclusion probabilities of individual edges.
#' \item \code{interactions}: A matrix with \code{p} rows and \code{p} columns,
#' containing model-averaged posterior means of the pairwise associations.
#' \item \code{thresholds}: A matrix with \code{p} rows and \code{max(m)}
#' columns, containing model-averaged category thresholds. In the case of
#' ``blume-capel'' variables, the first entry is the parameter for the linear
#' effect and the second entry is the parameter for the quadratic effect, which
#' models the offset to the reference category.
#' \item In the case of
#' \code{edge_prior = "Stochastic-Block"} two additional elements are returned:
#' a vector \code{allocations} with the estimated cluster assignments of the nodes and an
#' table \code{clusters} with the estimated posterior probability of the number
#' of clusters in the network. The vector of node allocations is calculated using
#' a method proposed by \insertCite{Dahl2009}{bgms} and also used by
#' \insertCite{GengEtAl_2019}{bgms}.
#' }
#' If \code{save = TRUE}, the result is a list of class ``bgms'' containing:
#' \itemize{
#' \item \code{indicator}: A matrix with \code{iter} rows and
#' \code{p * (p - 1) / 2} columns, containing the edge inclusion indicators from
#' every iteration of the Gibbs sampler.
#' \item \code{interactions}: A matrix with \code{iter} rows and
#' \code{p * (p - 1) / 2} columns, containing parameter states from every
#' iteration of the Gibbs sampler for the pairwise associations.
#' \item \code{thresholds}: A matrix with \code{iter} rows and
#' \code{sum(m)} columns, containing parameter states from every iteration of
#' the Gibbs sampler for the category thresholds.
#' \item In the case of
#' \code{edge_prior = "Stochastic-Block"} a matrix \code{allocations} with the
#' cluster assignments of the nodes from each iteration is returned. This matrix
#' can be used to calculate the posterior probability of the number of clusters
#' by utilizing the \code{summary_SBM(bgm_output[["allocations"]])} function.
#' }
#' Column averages of these matrices provide the model-averaged posterior means.
#' Except for the \code{allocations} matrix, for which the \code{summary_SBM}
#' needs to be utilized.
#'
#' In addition to the analysis results, the bgm output lists some of the
#' arguments of its call. This is useful for post-processing the results.
#'
#' @references
#'   \insertAllCited{}
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
#'       y = fit$indicator[lower.tri(fit$indicator)], ylim = c(0, 1),
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
#'  posterior.inclusion = fit$indicator[lower.tri(fit$indicator)]
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
#'  # THE MEDIAN PROBABILITY NETWORK
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
#'
#' @importFrom utils packageVersion
#' @importFrom Rcpp evalCpp
#' @importFrom Rdpack reprompt
#'
#' @export
bgm = function(x,
               variable_type = "ordinal",
               reference_category,
               iter = 1e4,
               burnin = 5e2,
               interaction_scale = 2.5,
               threshold_alpha = 0.5,
               threshold_beta = 0.5,
               edge_selection = TRUE,
               edge_prior = c("Bernoulli", "Beta-Bernoulli", "Stochastic-Block"),
               inclusion_probability = 0.5,
               beta_bernoulli_alpha = 1,
               beta_bernoulli_beta = 1,
               dirichlet_alpha = 1,
               na_action = c("listwise", "impute"),
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
  model = check_model(x = x,
                      variable_type = variable_type,
                      reference_category = reference_category,
                      interaction_scale = interaction_scale,
                      threshold_alpha = threshold_alpha,
                      threshold_beta = threshold_beta,
                      edge_selection = edge_selection,
                      edge_prior = edge_prior,
                      inclusion_probability = inclusion_probability,
                      beta_bernoulli_alpha = beta_bernoulli_alpha,
                      beta_bernoulli_beta = beta_bernoulli_beta,
                      dirichlet_alpha = dirichlet_alpha)

  # ----------------------------------------------------------------------------
  # The vector variable_type is now coded as boolean.
  # Ordinal (variable_bool == TRUE) or Blume-Capel (variable_bool == FALSE)
  # ----------------------------------------------------------------------------
  variable_bool = model$variable_bool
  # ----------------------------------------------------------------------------

  reference_category = model$reference_category
  edge_selection = model$edge_selection
  edge_prior = model$edge_prior
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

  #Check na_action -------------------------------------------------------------
  na_action_input = na_action
  na_action = try(match.arg(na_action), silent = TRUE)
  if(inherits(na_action, what = "try-error"))
    stop(paste0("The na_action argument should equal listwise or impute, not ",
                na_action_input,
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
  data = reformat_data(x = x,
                       na_action = na_action,
                       variable_bool = variable_bool,
                       reference_category = reference_category)
  x = data$x
  no_categories = data$no_categories
  missing_index = data$missing_index
  na_impute = data$na_impute
  reference_category = data$reference_category

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

  # Starting value of model matrix ---------------------------------------------
  indicator = matrix(1,
                 nrow = no_variables,
                 ncol = no_variables)


  #Starting values of interactions and thresholds (posterior mode) -------------
  interactions = matrix(0, nrow = no_variables, ncol = no_variables)
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
      Index[cntr, 2] = variable1 - 1
      Index[cntr, 3] = variable2 - 1
    }
  }

  #The Metropolis within Gibbs sampler -----------------------------------------
  out = gibbs_sampler(observations = x,
                      indicator = indicator,
                      interactions = interactions,
                      thresholds = thresholds,
                      no_categories  = no_categories,
                      interaction_scale = interaction_scale,
                      proposal_sd = proposal_sd,
                      proposal_sd_blumecapel = proposal_sd_blumecapel,
                      edge_prior = edge_prior,
                      theta = theta,
                      beta_bernoulli_alpha = beta_bernoulli_alpha,
                      beta_bernoulli_beta = beta_bernoulli_beta,
                      dirichlet_alpha = dirichlet_alpha,
                      Index = Index,
                      iter = iter,
                      burnin = burnin,
                      n_cat_obs = n_cat_obs,
                      sufficient_blume_capel = sufficient_blume_capel,
                      threshold_alpha = threshold_alpha,
                      threshold_beta = threshold_beta,
                      na_impute = na_impute,
                      missing_index = missing_index,
                      variable_bool = variable_bool,
                      reference_category = reference_category,
                      save = save,
                      display_progress = display_progress,
                      edge_selection = edge_selection)


  #Preparing the output --------------------------------------------------------
  arguments = list(
    no_variables = no_variables,
    no_cases = nrow(x),
    na_impute = na_impute,
    variable_type = variable_type,
    iter = iter,
    burnin = burnin,
    interaction_scale = interaction_scale,
    threshold_alpha = threshold_alpha,
    threshold_beta = threshold_beta,
    edge_selection = edge_selection,
    edge_prior = edge_prior,
    inclusion_probability = theta,
    beta_bernoulli_alpha = beta_bernoulli_alpha ,
    beta_bernoulli_beta =  beta_bernoulli_beta,
    dirichlet_alpha = dirichlet_alpha,
    na_action = na_action,
    save = save,
    version = packageVersion("bgms")
  )

  if(save == FALSE) {
    if(edge_selection == TRUE) {
      indicator = out$indicator
    }
    interactions = out$interactions
    tresholds = out$thresholds

    if(is.null(colnames(x))){
      data_columnnames = paste0("variable ", 1:no_variables)
      colnames(interactions) = data_columnnames
      rownames(interactions) = data_columnnames
      if(edge_selection == TRUE) {
        colnames(indicator) = data_columnnames
        rownames(indicator) = data_columnnames
      }
      rownames(thresholds) = data_columnnames
    } else {
      data_columnnames <- colnames(x)
      colnames(interactions) = data_columnnames
      rownames(interactions) = data_columnnames
      if(edge_selection == TRUE) {
        colnames(indicator) = data_columnnames
        rownames(indicator) = data_columnnames
      }
      rownames(thresholds) = data_columnnames
    }

    if(any(variable_bool)) {
      colnames(thresholds) = paste0("category ", 1:max(no_categories))
    } else {
      thresholds = thresholds[, 1:2]
      colnames(thresholds) = c("linear", "quadratic")
    }

    arguments$data_columnnames = data_columnnames

    if(edge_selection == TRUE) {
      if(edge_prior == "Stochastic-Block"){

        summarySbm <- summary_SBM(cluster_allocations = out$allocations)

        output = list(indicator = indicator,
                      interactions = interactions,
                      thresholds = thresholds,
                      allocations = summarySbm$allocations,
                      clusters = summarySbm$no_clusters,
                      arguments = arguments)
      } else {

        output = list(indicator = indicator,
                           interactions = interactions,
                           thresholds = thresholds,
                           arguments = arguments)
        }

    } else {
      output = list(interactions = interactions,
                    thresholds = thresholds,
                    arguments = arguments)
    }

    class(output) = "bgms"
    return(output)
  } else {
    if(edge_selection == TRUE) {
      indicator = out$indicator
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
      colnames(indicator) = names_vec
    }
    colnames(interactions) = names_vec
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
      dimnames(indicator) = list(Iter. = 1:iter, colnames(indicator))
    }
    dimnames(interactions) = list(Iter. = 1:iter, colnames(interactions))
    dimnames(thresholds) = list(Iter. = 1:iter, colnames(thresholds))

    arguments$data_columnnames = data_columnnames

    if(edge_selection == TRUE) {
      if(edge_prior == "Stochastic-Block"){
        output = list(indicator = indicator,
                      interactions = interactions,
                      thresholds = thresholds,
                      allocations = out$allocations,
                      arguments = arguments)
      } else {

        output = list(indicator = indicator,
                      interactions = interactions,
                      thresholds = thresholds,
                      arguments = arguments)
      }
    } else {
      output = list(interactions = interactions,
                    thresholds = thresholds,
                    arguments = arguments)
    }
    class(output) = "bgms"
    return(output)
  }
}
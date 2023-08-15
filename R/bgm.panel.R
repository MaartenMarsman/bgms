#' Bayesian structure learning in a cross-lagged Markov Random Field model for
#' mixed binary and ordinal variables using MCMC.
#'
#' The function \code{bgm.panel} explores the joint pseudoposterior distribution
#' of structures and parameters in a cross-lagged Markov Random Field model for
#' mixed binary and ordinal variables, which can be used to analyze mixed binary
#' and ordinal panel data.
#'
#' A discrete spike and slab prior distribution is stipulated on the pairwise
#' cross-sectional and cross-lagged interactions. By formulating it as a mixture
#' of mutually singular distributions, the function can use a combination of
#' Metropolis-Hastings and Gibbs sampling to create a Markov chain that has the
#' joint posterior distribution as invariant.
#'
#' @param x A data frame or matrix with \code{n} rows and \code{p} times \code{t}
#' columns containing binary and ordinal variables for \code{n} independent
#' observations and \code{p} variables in \code{t} measurement occasions.
#' Variables are recoded as non-negative integers \code{(0, 1, ..., m)} if not
#' already done. Unobserved categories are collapsed into other categories after
#' recoding (i.e., if category 1 is unobserved, the data will be recoded from
#' (0, 2) to (0, 1)).
#' @param no_nodes The number of nodes per timepoint.
#' @param no_timepoints The number of timepoints. Must be a positive number.
#' Does not include time t = 0.
#' @param iter The number of iterations of the Gibbs sampler. The default of
#' \code{1e4} is for illustrative purposes. For stable estimates, it is
#' recommended to run the Gibbs sampler for at least \code{1e5} iterations.
#' @param burnin The number of iterations of the Gibbs sampler before its output
#' is saved. Since it may take some time for the Gibbs sampler to converge to
#' the posterior distribution, it is recommended not to set this number too low.
#' @param cauchy_scale The scale of the Cauchy prior for interactions. Defaults
#' to \code{2.5}.
#' @param threshold_alpha,threshold_beta The shape parameters of the beta-prime
#' prior density for the threshold parameters. Must be positive values. If the
#' two values are equal, the prior density is symmetric about zero. If
#' \code{threshold_beta} is greater than \code{threshold_alpha}, the
#' distribution is skewed to the left, and if \code{threshold_beta} is less than
#' \code{threshold_alpha}, it is skewed to the right. Smaller values tend to
#' lead to more diffuse prior distributions.
#' @param cross_sectional_edge_prior,cross_lagged_edge_prior The prior
#' distribution for the cross_lagged edges of the panel network. Currently, two
#' prior distributions are implemented: The Bernoulli model
#' \code{edge_prior = "Bernoulli"} assumes that the probability that an edge
#' between two variables is included is equal to \code{inclusion_probability}
#' and independent of other edges or variables. When
#' \code{inclusion_probability = 0.5}, this implies that each network structure
#' receives the same prior weight. The Beta-Bernoulli model
#' \code{edge_prior = "Beta-Bernoulli"} assumes a beta prior for the unknown
#' inclusion probability with shape parameters \code{beta_bernoulli_alpha} and
#' \code{beta_bernoulli_beta}. If \code{beta_bernoulli_alpha = 1} and
#' \code{beta_bernoulli_beta = 1}, this means that networks with the same
#' complexity (number of edges) get the same prior weight. Defaults to
#' \code{cross_sectional_edge_prior = "Bernoulli"} and
#' \code{cross_lagged_edge_prior = "Bernoulli"}.
#' @param cross_sectional_inclusion_probability,cross_lagged_inclusion_probability
#' The prior edge inclusion probability for the Bernoulli model. Can be a single
#' probability, or a matrix of \code{p} rows and \code{p} columns specifying an
#' inclusion probability for each edge pair. Defaults to
#' \code{cross_sectional_inclusion_probability = 0.5} and
#' \code{cross_lagged_inclusion_probability = 0.5}.
#' @param cross_sectional_beta_bernoulli_alpha,cross_sectional_beta_bernoulli_beta,cross_lagged_beta_bernoulli_alpha,cross_lagged_beta_bernoulli_beta
#' The two shape parameters of the Beta prior density for the Bernoulli
#' inclusion probability. Must be positive numbers. Defaults to
#' \code{cross_sectional_beta_bernoulli_alpha = 1},
#' \code{cross_sectional_beta_bernoulli_beta = 1},
#' \code{cross_lagged_beta_bernoulli_alpha = 1}, and
#' \code{cross_lagged_beta_bernoulli_beta = 1}.
#' @param save Should the function collect and return all samples from the Gibbs
#' sampler (\code{save = TRUE})? Or should it only return the (model-averaged)
#' posterior means (\code{save = FALSE})? Defaults to \code{FALSE}.
#' @param display_progress Should the function show a progress bar
#' (\code{display_progress = TRUE})? Or not (\code{display_progress = FALSE})?
#' Defaults to \code{TRUE}.
#'
#' @return If \code{save = FALSE} (the default), the result is a list containing
#' the following matrices:
#' \itemize{
#' \item \code{gamma}: A symmetric matrix with \code{p} rows and \code{p} columns,
#' containing posterior inclusion probabilities of cross-sectional edges.
#' \item \code{delta}: A matrix with \code{p} rows and \code{p} columns,
#' containing posterior inclusion probabilities of cross-lagged edges.
#' \item \code{cross_sectional_interactions}: A matrix with \code{p} rows and
#' \code{p} columns, containing model-averaged posterior means of the pairwise
#' cross-sectional associations.
#' \item \code{cross_lagged_interactions}: A matrix with \code{p} rows and
#' \code{p} columns, containing model-averaged posterior means of the pairwise
#' cross-lagged associations.
#' \item \code{thresholds}: A matrix with \code{p * t} rows and \code{max(m)}
#' columns, containing model-averaged category thresholds.
#' }
#'
#' If \code{save = TRUE}, the result is a list containing:
#' \itemize{
#' \item \code{samples.gamma}: A matrix with \code{iter} rows and
#' \code{p * (p - 1) / 2} columns, containing the cross-sectional edge inclusion
#' indicators from every iteration of the Gibbs sampler.
#' \item \code{samples.delta}: A matrix with \code{iter} rows and \code{p * p}
#' columns, containing the cross-lagged edge inclusion indicators from every
#' iteration of the Gibbs sampler.
#' \item \code{samples.cross.sectional.interactions}: A matrix with \code{iter}
#' rows and \code{p * (p - 1) / 2} columns, containing parameter states from
#' every iteration of the Gibbs sampler for the pairwise cross-sectional
#' associations.
#' \item \code{samples.cross.lagged.interactions}: A matrix with \code{iter}
#' rows and \code{p * p} columns, containing parameter states from every
#' iteration of the Gibbs sampler for the pairwise cross-lagged associations.
#' \item \code{samples.thresholds}: A matrix with \code{iter} rows and
#' \code{sum(m) * t} columns, containing parameter states from every iteration of
#' the Gibbs sampler for the category thresholds.
#' }
#' Column averages of these matrices provide the model-averaged posterior means.
#'
#' @export
bgm.panel = function(x,
                     no_nodes,
                     no_timepoints,
                     iter = 1e4,
                     burnin = 1e3,
                     cauchy_scale = 2.5,
                     cross_sectional_edge_prior = c("Bernoulli", "Beta-Bernoulli"),
                     cross_lagged_edge_prior = c("Bernoulli", "Beta-Bernoulli"),
                     cross_sectional_inclusion_probability = 0.5,
                     cross_lagged_inclusion_probability = 0.5,
                     cross_sectional_beta_bernoulli_alpha = 1,
                     cross_sectional_beta_bernoulli_beta = 1,
                     cross_lagged_beta_bernoulli_alpha = 1,
                     cross_lagged_beta_bernoulli_beta = 1,
                     threshold_alpha = 1,
                     threshold_beta = 1,
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

  #Check design input
  if(no_nodes < 2 || is.null(no_nodes))
    stop("The design needs at least two variables.")
  if(no_timepoints <= 0 || is.null(no_timepoints))
    stop("The design needs at least one timepoint.")
  if(ncol(x) != no_nodes * (no_timepoints + 1))
    stop(paste0("The number of variables in x does not match the design input. \n",
                "It is expected that x has no_nodes * (1 + no_timepoints) = ",
                no_nodes * (no_timepoints + 1), " columns. Instead it has ",
                ncol(x), " columns."))

  #Check Gibbs input -----------------------------------------------------------
  if(abs(iter - round(iter)) > sqrt(.Machine$double.eps))
    stop("Parameter ``iter'' needs to be a positive integer.")
  if(abs(burnin - round(burnin)) > sqrt(.Machine$double.eps) || burnin < 0)
    stop("Parameter ``burnin'' needs to be a non-negative integer.")

  #Check prior set-up for the interaction parameters ---------------------------
  if(cauchy_scale <= 0 || is.na(cauchy_scale) || is.infinite(cauchy_scale))
    stop("The scale of the Cauchy prior needs to be positive.")

  #Check prior set-up for the cross-sectional edge indicators ------------------
  cross_sectional_edge_prior = match.arg(cross_sectional_edge_prior)
  if(cross_sectional_edge_prior == "Bernoulli") {
    if(length(cross_sectional_inclusion_probability) == 1) {
      crsec_theta = cross_sectional_inclusion_probability[1]
      if(is.na(crsec_theta) || is.null(crsec_theta))
        stop("There is no value specified for the cross-sectional inclusion probability.")
      if(crsec_theta <= 0)
        stop("The cross-sectional inclusion probability needs to be positive.")
      if(crsec_theta >= 1)
        stop("The cross-sectional inclusion probability cannot exceed the value one.")
      crsec_theta = matrix(crsec_theta, nrow = no_nodes, ncol = no_nodes)
    } else {
      if(!inherits(cross_sectional_inclusion_probability, what = "matrix") ||
         !inherits(cross_sectional_inclusion_probability, what = "data.frame"))
        stop("The input for the cross-sectional inclusion probability argument needs to be a single number, matrix, or dataframe.")

      if(inherits(cross_sectional_inclusion_probability, what = "data.frame")) {
        crsec_theta = data.matrix(cross_sectional_inclusion_probability)
      } else {
        crsec_theta = cross_sectional_inclusion_probability
      }
      if(!isSymmetric(crsec_theta))
        stop("The cross-sectional inclusion probability matrix needs to be symmetric.")
      if(ncol(crsec_theta) != ncol(x))
        stop("The cross-sectional inclusion probability matrix needs to have as many rows (columns) as there are variables in the data.")

      if(any(is.na(crsec_theta[lower.tri(crsec_theta)])) ||
         any(is.null(crsec_theta[lower.tri(crsec_theta)])))
        stop("One or more elements of the elements in cross-sectional inclusion probability matrix are not specified.")
      if(any(crsec_theta <= 0))
        stop("The cross-sectional inclusion probability matrix contains negative values, the values need to be positive.")
      if(theta >= 1)
        stop("The cross-sectional inclusion probability matrix contains values greater than one; inclusion probabilities cannot exceed the value one.")
    }
  }
  if(cross_sectional_edge_prior == "Beta-Bernoulli") {
    crsec_theta = matrix(0.5, nrow = no_nodes, ncol = no_nodes)
    if(cross_sectional_beta_bernoulli_alpha <= 0 ||
       cross_sectional_beta_bernoulli_beta <= 0)
      stop("The scale parameters of the beta distribution need to be positive.")
    if(!is.finite(cross_sectional_beta_bernoulli_alpha) ||
       !is.finite(cross_sectional_beta_bernoulli_beta))
      stop("The scale parameters of the beta distribution need to be finite.")
    if(is.na(cross_sectional_beta_bernoulli_alpha) ||
       is.na(cross_sectional_beta_bernoulli_beta) ||
       is.null(cross_sectional_beta_bernoulli_alpha) ||
       is.null(cross_sectional_beta_bernoulli_beta))
      stop("Values for both scale parameters of the beta distribution need to be specified.")
  }

  #Check prior set-up for the cross-lagged edge indicators ---------------------
  cross_lagged_edge_prior = match.arg(cross_lagged_edge_prior)
  if(cross_lagged_edge_prior == "Bernoulli") {
    if(length(cross_lagged_inclusion_probability) == 1) {
      crlag_theta = cross_lagged_inclusion_probability[1]
      if(is.na(crlag_theta) || is.null(crlag_theta))
        stop("There is no value specified for the cross-lagged inclusion probability.")
      if(crlag_theta <= 0)
        stop("The cross-lagged inclusion probability needs to be positive.")
      if(crlag_theta >= 1)
        stop("The cross-lagged inclusion probability cannot exceed the value one.")
      crlag_theta = matrix(crlag_theta, nrow = no_nodes, ncol = no_nodes)
    } else {
      if(!inherits(cross_lagged_inclusion_probability, what = "matrix") ||
         !inherits(cross_lagged_inclusion_probability, what = "data.frame"))
        stop("The input for the cross-lagged inclusion probability argument needs to be a single number, matrix, or dataframe.")

      if(inherits(cross_lagged_inclusion_probability, what = "data.frame")) {
        crlag_theta = data.matrix(cross_lagged_inclusion_probability)
      } else {
        crlag_theta = cross_lagged_inclusion_probability
      }
      if(!isSymmetric(crlag_theta))
        stop("The cross-lagged inclusion probability matrix needs to be symmetric.")
      if(ncol(crlag_theta) != ncol(x))
        stop("The cross-lagged inclusion probability matrix needs to have as many rows (columns) as there are variables in the data.")

      if(any(is.na(crlag_theta[lower.tri(crlag_theta)])) ||
         any(is.null(crlag_theta[lower.tri(crlag_theta)])))
        stop("One or more elements of the elements in cross-lagged inclusion probability matrix are not specified.")
      if(any(crlag_theta <= 0))
        stop("The cross-lagged inclusion probability matrix contains negative values, the values need to be positive.")
      if(theta >= 1)
        stop("The cross-lagged inclusion probability matrix contains values greater than one; inclusion probabilities cannot exceed the value one.")
    }
  }
  if(cross_lagged_edge_prior == "Beta-Bernoulli") {
    crlag_theta = matrix(0.5, nrow = no_nodes, ncol = no_nodes)
    if(cross_lagged_beta_bernoulli_alpha <= 0 ||
       cross_lagged_beta_bernoulli_beta <= 0)
      stop("The scale parameters of the beta distribution need to be positive.")
    if(!is.finite(cross_lagged_beta_bernoulli_alpha) ||
       !is.finite(cross_lagged_beta_bernoulli_beta))
      stop("The scale parameters of the beta distribution need to be finite.")
    if(is.na(cross_lagged_beta_bernoulli_alpha) ||
       is.na(cross_lagged_beta_bernoulli_beta) ||
       is.null(cross_lagged_beta_bernoulli_alpha) ||
       is.null(cross_lagged_beta_bernoulli_beta))
      stop("Values for both scale parameters of the beta distribution need to be specified.")
  }

  #Check prior set-up for the threshold parameters -----------------------------
  if(threshold_alpha <= 0  | !is.finite(threshold_alpha))
    stop("Parameter ``threshold_alpha'' needs to be positive.")
  if(threshold_beta <= 0  | !is.finite(threshold_beta))
    stop("Parameter ``threshold_beta'' needs to be positive.")

  #Format the data input -------------------------------------------------------
  ### FOR NOW "na.action = LISWISE"
  data = reformat_data_bgm_panel(x = x,
                                 no_nodes,
                                 no_timepoints,
                                 na.action = "listwise")
  x = data$x
  no_persons = nrow(x)
  no_categories = data$no_categories
  missing_index = data$missing_index
  na.impute = data$na.impute

  no_cross_sectional_interactions = no_nodes * (no_nodes - 1) / 2
  no_cross_lagged_interactions = no_nodes * no_nodes
  no_thresholds = sum(no_categories) * no_timepoints

  #Set up the variance of the (normal) proposal distribution
  crsec_proposal_sd = matrix(1,
                             nrow = no_nodes,
                             ncol = no_nodes)
  crlag_proposal_sd = matrix(1,
                             nrow = no_nodes,
                             ncol = no_nodes)

  #Starting value of model matrix:
  gamma = matrix(1,
                 nrow = no_nodes,
                 ncol = no_nodes)
  delta = matrix(1,
                 nrow = no_nodes,
                 ncol = no_nodes)

  crsec_interactions = matrix(0, nrow = no_nodes, ncol = no_nodes)
  crlag_interactions = matrix(0, nrow = no_nodes, ncol = no_nodes)
  thresholds = matrix(0,
                      nrow = no_nodes * no_timepoints,
                      ncol = max(no_categories))

  #Precomputing number of observations per category for each node.
  n_cat_obs = matrix(0,
                     nrow = max(no_categories) + 1,
                     ncol = no_nodes * no_timepoints)
  start = (0:no_timepoints) * no_nodes
  for(t in 1:no_timepoints) {
    for(node in 1:no_nodes) {
      for(category in 0:no_categories[node]) {
        n_cat_obs[category + 1, node + start[t]] =
          sum(x[, node + start[t + 1]] == category)
      }
    }
  }

  # Index vector used to sample interactions in a random order.
  crsec_Index = matrix(0,
                       nrow = no_nodes * (no_nodes - 1) / 2,
                       ncol = 3)
  cntr = 0
  for(node1 in 1:(no_nodes - 1)) {
    for(node2 in (node1 + 1):no_nodes) {
      cntr =  cntr + 1
      crsec_Index[cntr, 1] = cntr
      crsec_Index[cntr, 2] = node1
      crsec_Index[cntr, 3] = node2
    }
  }

  crlag_Index = matrix(0,
                       nrow = no_nodes * no_nodes,
                       ncol = 3)
  cntr = 0
  for(node1 in 1:no_nodes) {
    for(node2 in 1:no_nodes) {
      cntr =  cntr + 1
      crlag_Index[cntr, 1] = cntr
      crlag_Index[cntr, 2] = node1
      crlag_Index[cntr, 3] = node2
    }
  }

  #The Metropolis within Gibbs sampler -----------------------------------------
  out = gibbs_sampler_cross_lagged_mrf(observations = x,
                                       no_persons = no_persons,
                                       no_nodes = no_nodes,
                                       no_timepoints = no_timepoints,
                                       gamma = gamma,
                                       delta = delta,
                                       crsec_interactions = crsec_interactions,
                                       crlag_interactions = crlag_interactions,
                                       thresholds = thresholds,
                                       no_categories = no_categories,
                                       start = start,
                                       cauchy_scale = cauchy_scale,
                                       crsec_proposal_sd = crsec_proposal_sd,
                                       crlag_proposal_sd = crlag_proposal_sd,
                                       crsec_edge_prior = cross_sectional_edge_prior,
                                       crlag_edge_prior = cross_lagged_edge_prior,
                                       crsec_theta = crsec_theta,
                                       crlag_theta = crlag_theta,
                                       crsec_beta_bernoulli_alpha = cross_sectional_beta_bernoulli_alpha,
                                       crsec_beta_bernoulli_beta = cross_sectional_beta_bernoulli_beta,
                                       crlag_beta_bernoulli_alpha = cross_lagged_beta_bernoulli_alpha,
                                       crlag_beta_bernoulli_beta = cross_lagged_beta_bernoulli_beta,
                                       crsec_Index = crsec_Index,
                                       crlag_Index = crlag_Index,
                                       iter = iter,
                                       burnin = burnin,
                                       n_cat_obs = n_cat_obs,
                                       threshold_alpha = threshold_alpha,
                                       threshold_beta = threshold_beta,
                                       save = save,
                                       display_progress = display_progress)

  #Preparing the output --------------------------------------------------------
  ### FOR NOW NO NAMES
  if(save == FALSE) {
    gamma = out$gamma
    delta = out$delta
    crsec_interactions = out$crsec_interactions
    crlag_interactions = out$crlag_interactions
    tresholds = out$thresholds

    return(list(gamma = gamma,
                delta = delta,
                cross_sectional_interactions = crsec_interactions,
                cross_lagged_interactions = crlag_interactions,
                thresholds = tresholds))
  } else {
    gamma = out$gamma
    delta = out$delta
    crsec_interactions = out$crsec_interactions
    crlag_interactions = out$crlag_interactions
    thresholds = out$thresholds

    return(list(gamma = gamma,
                delta = delta,
                cross_sectional_interactions = crsec_interactions,
                cross_lagged_interactions = crlag_interactions,
                thresholds = thresholds))
  }
}
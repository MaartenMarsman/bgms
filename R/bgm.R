#' Bayesian structure learning in Markov Random Fields of mixed binary and
#' ordinal variables using MCMC.
#'
#' The function \code{bgm} explores the joint pseudoposterior distribution of
#' structures and parameters in a Markov Random Field for mixed binary and
#' ordinal variables.
#'
#' A discrete spike and slab prior distribution is stipulated on the pairwise
#' interactions. By formulating it as a mixture of mutually singular
#' distributions, the function can use a combination of Metropolis-Hastings and
#' Gibbs sampling to create a Markov chain that has the joint posterior
#' distribution as invariant. Current options for the slab distribution are the
#' unit-information prior or a Cauchy with an optional scaling parameter. A
#' Beta-prime distribution is used for the exponent of the category parameters.
#' A uniform prior is used for edge inclusion variables (i.e., the prior
#' probability that the edge is included is 0.5).
#'
#' @param x A data frame or matrix with \code{n} rows and \code{p} columns
#' containing binary and ordinal variables for \code{n} independent observations
#' and \code{p} variables in the network. Variables are recoded as non-negative
#' integers \code{(0, 1, ..., m)} if not already done. Unobserved categories are
#' collapsed into other categories after recoding (i.e., if category 1 is
#' unobserved, the data will be recoded from (0, 2) to (0, 1)).
#' @param iter The number of iterations of the Gibbs sampler. The default of
#' \code{1e4} is for illustrative purposes. For stable estimates, it is
#' recommended to run the Gibbs sampler for at least \code{1e5} iterations.
#' @param burnin The number of iterations of the Gibbs sampler before its output
#' is saved. Since it may take some time for the Gibbs sampler to converge to
#' the posterior distribution, it is recommended not to set this number too low.
#' @param interaction_prior The type of prior distribution that is used for the
#' interaction effects. Currently, two prior densities are implemented: The
#' Cauchy prior (\code{interaction_prior = "Cauchy"}) and the Unit Information
#' prior (\code{interaction_prior = "UnitInfo"}).
#' @param cauchy_scale The scale of the Cauchy prior for interactions. Defaults
#' to \code{2.5}.
#' @param edge_prior The prior distribution for the edges or structure of the
#' network. Two prior distributions are currently implemented: The Bernoulli
#' model \code{edge_prior = "Bernoulli"} assumes that the probability that an
#' edge between two variables is included is equal to
#' \code{inclusion_probability} and independent of other edges or variables.
#' When \code{inclusion_probability = 0.5}, this implies that each network
#' structure receives the same prior weight. The Beta-Bernoulli model
#' \code{edge_prior = "Beta-Bernoulli"} assumes a beta prior for the unknown
#' inclusion probability with shape parameters \code{beta_bernoulli_alpha} and
#' \code{beta_bernoulli_beta}. If \code{beta_bernoulli_alpha = 1} and
#' \code{beta_bernoulli_beta = 1}, this means that networks with the same
#' complexity (number of edges) receive the same prior weight. Defaults to
#' \code{edge_prior = "Bernoulli"}.
#' @param inclusion_probability The prior edge inclusion probability for the
#' Bernoulli model. Can be a single probability, or a matrix of \code{p} rows
#' and \code{p} columns specifying an inclusion probability for each edge pair.
#' Defaults to \code{inclusion_probability = 0.5}.
#' @param beta_bernoulli_alpha,beta_bernoulli_beta The two shape parameters of
#' the Beta prior density for the Bernoulli inclusion probability. Must be
#' positive numbers. Defaults to \code{beta_bernoulli_alpha = 1} and
#' \code{beta_bernoulli_beta = 1}.
#' @param threshold_alpha,threshold_beta The shape parameters of the beta-prime
#' prior density for the threshold parameters. Must be positive values. If the
#' two values are equal, the prior density is symmetric about zero. If
#' \code{threshold_beta} is greater than \code{threshold_alpha}, the
#' distribution is skewed to the left, and if \code{threshold_beta} is less than
#' \code{threshold_alpha}, it is skewed to the right. Smaller values tend to
#' lead to more diffuse prior distributions.
#' @param adaptive Should the function use an adaptive Metropolis algorithm to
#' update interaction parameters within the model? If \code{adaptive = TRUE},
#' bgm adjusts the proposal variance to match the acceptance probability of the
#' random walk Metropolis algorithm to be close to the optimum of \code{.234}
#' using a Robbins-Monro type algorithm. If \code{adaptive = FALSE}, it sets the
#' proposal variance to the inverse of the observed Fisher information matrix
#' (the second derivative at the posterior mode). Defaults to \code{FALSE}.
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
#' Defaults to \code{TRUE}.
#'
#' @return If \code{save = FALSE} (the default), the result is a list of class
#' ``bgms'' containing the following matrices:
#' \itemize{
#' \item \code{gamma}: A matrix with \code{p} rows and \code{p} columns,
#' containing posterior inclusion probabilities of individual edges.
#' \item \code{interactions}: A matrix with \code{p} rows and \code{p} columns,
#' containing model-averaged posterior means of the pairwise associations.
#' \item \code{thresholds}: A matrix with \code{p} rows and \code{max(m)}
#' columns, containing model-averaged category thresholds.
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
#'  #The bgms package currently assumes that the prior odds are 1:
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
               iter = 1e4,
               burnin = 1e3,
               interaction_prior = c("Cauchy", "UnitInfo"),
               cauchy_scale = 2.5,
               edge_prior = c("Bernoulli", "Beta-Bernoulli"),
               inclusion_probability = 0.5,
               beta_bernoulli_alpha = 1,
               beta_bernoulli_beta = 1,
               threshold_alpha = 1,
               threshold_beta = 1,
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

  #Check Gibbs input -----------------------------------------------------------
  if(abs(iter - round(iter)) > sqrt(.Machine$double.eps))
    stop("Parameter ``iter'' needs to be a positive integer.")
  if(abs(burnin - round(burnin)) > sqrt(.Machine$double.eps) || burnin < 0)
    stop("Parameter ``burnin'' needs to be a non-negative integer.")

  #Check prior set-up for the interaction parameters ---------------------------
  interaction_prior = match.arg(interaction_prior)
  if(interaction_prior == "Cauchy") {
    if(cauchy_scale <= 0 || is.na(cauchy_scale) || is.infinite(cauchy_scale))
      stop("The scale of the Cauchy prior needs to be positive.")
  }

  #Check prior set-up for the edge indicators ----------------------------------
  edge_prior = match.arg(edge_prior)
  if(edge_prior == "Bernoulli") {
    if(length(inclusion_probability) == 1) {
      theta = inclusion_probability[1]
      if(is.na(theta) || is.null(theta))
        stop("There is no value specified for the inclusion probability.")
      if(theta <= 0)
        stop("The inclusion probability needs to be positive.")
      if(theta >= 1)
        stop("The inclusion probability cannot exceed the value one.")
      theta = matrix(theta, nrow = ncol(x), ncol = ncol(x))
    } else {
      if(!inherits(inclusion_probability, what = "matrix") ||
         !inherits(inclusion_probability, what = "data.frame"))
        stop("The input for the inclusion probability argument needs to be a single number, matrix, or dataframe.")

      if(inherits(inclusion_probability, what = "data.frame")) {
        theta = data.matrix(inclusion_probability)
      } else {
        theta = inclusion_probability
      }
      if(!isSymmetric(theta))
        stop("The inclusion probability matrix needs to be symmetric.")
      if(ncol(theta) != ncol(x))
        stop("The inclusion probability matrix needs to have as many rows (columns) as there are variables in the data.")

      if(any(is.na(theta[lower.tri(theta)])) ||
         any(is.null(theta[lower.tri(theta)])))
        stop("One or more elements of the elements in inclusion probability matrix are not specified.")
      if(any(theta <= 0))
        stop("The inclusion probability matrix contains negative values, the values need to be positive.")
      if(theta >= 1)
        stop("The inclusion probability matrix contains values greater than one; inclusion probabilities cannot exceed the value one.")
    }
  }
  if(edge_prior == "Beta-Bernoulli") {
    theta = matrix(0.5, nrow = ncol(x), ncol = ncol(x))
    if(beta_bernoulli_alpha <= 0 || beta_bernoulli_beta <= 0)
      stop("The scale parameters of the beta distribution need to be positive.")
    if(!is.finite(beta_bernoulli_alpha) || !is.finite(beta_bernoulli_beta))
      stop("The scale parameters of the beta distribution need to be finite.")
    if(is.na(beta_bernoulli_alpha) || is.na(beta_bernoulli_beta) ||
       is.null(beta_bernoulli_alpha) || is.null(beta_bernoulli_beta))
      stop("Values for both scale parameters of the beta distribution need to be specified.")
  }

  #Check prior set-up for the threshold parameters -----------------------------
  if(threshold_alpha <= 0  | !is.finite(threshold_alpha))
    stop("Parameter ``threshold_alpha'' needs to be positive.")
  if(threshold_beta <= 0  | !is.finite(threshold_beta))
    stop("Parameter ``threshold_beta'' needs to be positive.")

  #Check na.action -------------------------------------------------------------
  na.action = match.arg(na.action)

  #Format the data input -------------------------------------------------------
  data = reformat_data_bgm(x = x, na.action)
  x = data$x
  no_categories = data$no_categories
  missing_index = data$missing_index
  na.impute = data$na.impute

  if(na.impute == TRUE) {
    if(interaction_prior != "Cauchy")
      warning(paste0(
        "There were missing responses and na.action was set to ``impute''. The \n",
        "bgm function must switch the interaction_prior to ``Cauchy''."))
    adaptive = TRUE
    interaction_prior = "Cauchy"
    if(cauchy_scale <= 0 || is.na(cauchy_scale) || is.infinite(cauchy_scale))
      stop("The scale of the Cauchy prior needs to be positive.")
  }

  no_nodes = ncol(x)
  no_interactions = no_nodes * (no_nodes - 1) / 2
  no_thresholds = sum(no_categories)

  #Proposal set-up for the interaction parameters ------------------------------
  if(interaction_prior == "UnitInfo") {
    pps = try(mppe(x = x,
                   interaction_prior = interaction_prior),
              silent = TRUE)
    if(inherits(pps, what = "try-error"))
      stop(paste0(
        "For the Unit Information prior we need to estimate the posterior mode.\n",
        "Unfortunately, we could not find this mode for your data. Please try the\n",
        "Cauchy prior option."))
    unit_info = sqrt(pps$unit_info)
  } else {
    if(!na.impute) {
      pps = try(mppe(x = x,
                     interaction_prior = interaction_prior,
                     cauchy_scale = cauchy_scale),
                silent = TRUE)
      if(inherits(pps, what = "try-error") & adaptive == FALSE) {
        stop(paste0(
          "By default, the MCMC procedure underlying the bgm function uses a \n",
          "Metropolis algorithm with a fixed proposal distribution. We attempt to \n",
          "fit this proposal distribution to the target posterior distribution by \n",
          "locating the posterior mode and using information about the curvature \n",
          "around that mode to set the variance of the proposal distributions. \n",
          "Unfortunately, we were unable to locate the posterior mode for your data.\n",
          "Please try again with ``adaptive = TRUE''."))
      }
    }
    unit_info = matrix(data = NA, nrow = 1, ncol = 1)
  }

  #Set up the variance of the (normal) proposal distribution
  proposal_sd = matrix(1,
                       nrow = no_nodes,
                       ncol = no_nodes)
  if(adaptive == FALSE && !na.impute) {
    hessian = pps$hessian[-c(1:no_thresholds), -c(1:no_thresholds)]
    cntr = 0
    for(node1 in 1:(no_nodes - 1)) {
      for(node2 in (node1 + 1):no_nodes) {
        cntr = cntr + 1
        proposal_sd[node1, node2] = sqrt(-1 / hessian[cntr, cntr])
        proposal_sd[node2, node1] = proposal_sd[node1, node2]
      }
    }
  }

  # # Starting value of model matrix:
  gamma = matrix(1,
                 nrow = no_nodes,
                 ncol = no_nodes)

  #Starting values of interactions and thresholds (posterior mode)
  if(!na.impute && !inherits(pps, what = "try-error")) {
    interactions = pps$interactions
    thresholds = pps$thresholds
  } else {
    interactions = matrix(0, nrow = no_nodes, ncol = no_nodes)
    thresholds = matrix(0, nrow = no_nodes, ncol = max(no_categories))
  }

  #Precomputing number of observations per category for each node.
  n_cat_obs = matrix(0,
                     nrow = max(no_categories) + 1,
                     ncol = no_nodes)
  for(node in 1:no_nodes) {
    for(category in 0:no_categories[node]) {
      n_cat_obs[category + 1, node] = sum(x[, node] == category)
    }
  }

  # Index vector used to sample interactions in a random order.
  Index = matrix(0,
                 nrow = no_nodes * (no_nodes - 1) / 2,
                 ncol = 3)
  cntr = 0
  for(node1 in 1:(no_nodes - 1)) {
    for(node2 in (node1 + 1):no_nodes) {
      cntr =  cntr + 1
      Index[cntr, 1] = cntr
      Index[cntr, 2] = node1
      Index[cntr, 3] = node2
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
                      edge_prior = edge_prior,
                      theta = theta,
                      beta_bernoulli_alpha = beta_bernoulli_alpha,
                      beta_bernoulli_beta = beta_bernoulli_beta,
                      Index = Index,
                      iter = iter,
                      burnin = burnin,
                      n_cat_obs = n_cat_obs,
                      threshold_alpha = threshold_alpha,
                      threshold_beta = threshold_beta,
                      na.impute,
                      missing_index,
                      adaptive = adaptive,
                      save = save,
                      display_progress = display_progress)

  #Preparing the output --------------------------------------------------------
  if(save == FALSE) {
    gamma = out$gamma
    interactions = out$interactions
    tresholds = out$thresholds

    colnames(interactions) = paste0("node ", 1:no_nodes)
    rownames(interactions) = paste0("node ", 1:no_nodes)
    colnames(gamma) = paste0("node ", 1:no_nodes)
    rownames(gamma) = paste0("node ", 1:no_nodes)
    colnames(tresholds) = paste0("category ", 1:max(no_categories))
    rownames(tresholds) = paste0("node ", 1:no_nodes)

    output = list(gamma = gamma,
                  interactions = interactions,
                  thresholds = thresholds,
                  edge_prior = edge_prior,
                  inclusion_probability = inclusion_probability,
                  beta_bernoulli_alpha = beta_bernoulli_alpha,
                  beta_bernoulli_beta = beta_bernoulli_beta,
                  save = save)
    class(output) = "bgms"
    return(output)
  } else {
    gamma = out$gamma
    interactions = out$interactions
    thresholds = out$thresholds

    names1 = names2 = character(length = no_nodes * (no_nodes - 1) / 2)
    cntr = 0
    for(node in 1:(no_nodes - 1)) {
      for(node_2 in (node + 1):no_nodes) {
        cntr = cntr + 1
        names1[cntr] = paste0("sigma(",node, ", ",node_2,")")
        names2[cntr] = paste0("gamma(",node, ", ",node_2,")")
      }
    }
    colnames(gamma) = names2
    colnames(interactions) = names1

    names = character(length = sum(no_categories))
    cntr = 0
    for(node in 1:no_nodes) {
      for(category in 1:no_categories[node]) {
        cntr = cntr + 1
        names[cntr] = paste0("threshold(",node, ", ",category,")")
      }
    }
    colnames(thresholds) = names

    rownames(gamma) = paste0("Iter. ", 1:iter)
    rownames(interactions) = paste0("Iter. ", 1:iter)
    rownames(thresholds) = paste0("Iter. ", 1:iter)

    output = list(gamma = gamma,
                  interactions = interactions,
                  thresholds = thresholds,
                  edge_prior = edge_prior,
                  inclusion_probability = inclusion_probability,
                  beta_bernoulli_alpha = beta_bernoulli_alpha,
                  beta_bernoulli_beta = beta_bernoulli_beta,
                  save = save)
    class(output) = "bgms"
    return(output)
  }
}


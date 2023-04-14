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
#' @param x A matrix with \code{n} rows and \code{p} columns, containing binary
#' and ordinal variables for \code{n} independent observations and \code{p}
#' variables in the network. Variables are recoded as non-negative integers
#' \code{(0, 1, ..., m)} if not done already. Unobserved categories are
#' collapsed into other categories after recoding. See
#' \code{reformat_data} for details.
#'
#' @param iter The number of iterations of the Gibbs sampler. Defaults to
#' \code{1e4}. For better estimates, it is recommended to run the procedure for
#' at least \code{1e5} iterations.
#' @param burnin The number of burnin iterations. The output of the Gibbs
#' sampler is stored after burnin iterations.
#' @param interaction_prior The prior distribution for the interaction effects.
#' Currently, two prior densities are implemented: The Unit Information prior
#' (\code{interaction_prior = "UnitInfo"}) and the Cauchy prior
#' (\code{interaction_prior = "Cauchy"}). Defaults to \code{"UnitInfo"}.
#' @param cauchy_scale The scale of the Cauchy prior for interactions. Defaults
#' to \code{2.5}.
#' @param threshold_alpha,threshold_beta The shape parameters of the Beta-prime
#' prior for the thresholds. Defaults to \code{1}.
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
#' \item \code{gamma}: A matrix with \code{p} rows and \code{p} columns,
#' containing posterior inclusion probabilities of individual edges.
#' \item \code{interactions}: A matrix with \code{p} rows and \code{p} columns,
#' containing model-averaged posterior means of the pairwise associations.
#' \item \code{thresholds}: A matrix with \code{p} rows and \code{max(m)}
#' columns, containing model-averaged category thresholds.
#' }
#'
#' If \code{save = TRUE}, the result is a list containing:
#' \itemize{
#' \item \code{samples.gamma}: A matrix with \code{iter} rows and
#' \code{p * (p - 1) / 2} columns, containing the edge inclusion indicators from
#' every iteration of the Gibbs sampler.
#' \item \code{samples.interactions}: A matrix with \code{iter} rows and
#' \code{p * (p - 1) / 2} columns, containing parameter states from every
#' iteration of the Gibbs sampler for the pairwise associations.
#' \item \code{samples.thresholds}: A matrix with \code{iter} rows and
#' \code{sum(m)} columns, containing parameter states from every iteration of
#' the Gibbs sampler for the category thresholds.
#' }
#' Column averages of these matrices provide the model-averaged posterior means.
#'
#' @examples
#' \dontrun{
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
#'  }
#' @export
bgm = function(x,
               iter = 1e4,
               burnin = 1e3,
               interaction_prior = c("UnitInfo", "Cauchy"),
               cauchy_scale = 2.5,
               threshold_alpha = 1,
               threshold_beta = 1,
               save = FALSE,
               display_progress = TRUE) {

  #Check Gibbs input -----------------------------------------------------------
  if(abs(iter - round(iter)) > sqrt(.Machine$double.eps))
    stop("Parameter ``iter'' needs to be a positive integer.")
  if(abs(burnin - round(burnin)) > sqrt(.Machine$double.eps) || burnin < 0)
    stop("Parameter ``burnin'' needs to be a non-negative integer.")

  #Check prior set-up for the interaction parameters ---------------------------
  interaction_prior = match.arg(interaction_prior)
  if(interaction_prior == "Cauchy") {
    if(cauchy_scale <= 0)
      stop("The scale of the Cauchy prior needs to be positive.")
  }

  #Check prior set-up for the threshold parameters -----------------------------
  if(threshold_alpha <= 0  | !is.finite(threshold_alpha))
    stop("Parameter ``threshold_alpha'' needs to be positive.")
  if(threshold_beta <= 0  | !is.finite(threshold_beta))
    stop("Parameter ``threshold_beta'' needs to be positive.")

  #Check data input ------------------------------------------------------------
  if(!inherits(x, what = "matrix"))
    stop("The input x is supposed to be a matrix.")

  if(ncol(x) < 2)
    stop("The matrix x should have more than one variable (columns).")
  if(nrow(x) < 2)
    stop("The matrix x should have more than one observation (rows).")

  #Format the data input -------------------------------------------------------
  data = reformat_data(x = x)
  x = data$x
  no_categories = data$no_categories

  no_nodes = ncol(x)
  no_interactions = no_nodes * (no_nodes - 1) / 2
  no_thresholds = sum(no_categories)

  #Proposal set-up for the interaction parameters ------------------------------
  if(interaction_prior == "Cauchy") {
    pps = try(mppe(x = x,
                   no_categories = no_categories,
                   interaction_prior = interaction_prior,
                   cauchy_scale = cauchy_scale),
              silent = TRUE)
  } else {
    pps = try(mppe(x = x,
                   no_categories = no_categories,
                   interaction_prior = interaction_prior),
              silent = TRUE)
  }
  if(inherits(pps, what = "try-error"))
    stop(paste0(
      "We use a Metropolis within Gibbs algorithm to learn the structure of the MRF.\n",
      "The proposal distribution used in the Metropolis algorithm is a normal distribution\n",
      "that takes as its variance a function of the second derivatives of the MRF evaluated\n",
      "at the posterior mode. Unfortunately, we could not find this mode for your data.\n",
      "Perhaps there were low category counts?"))

  if(interaction_prior == "UnitInfo") {
    unit_info = sqrt(pps$unit_info)
  } else {
    unit_info = matrix(data = NA, nrow = 1, ncol = 1)
  }

  #Set up the variance of the (normal) proposal distribution
  hessian = pps$hessian[-c(1:no_thresholds), -c(1:no_thresholds)]

  proposal_sd = matrix(0,
                       nrow = no_nodes,
                       ncol = no_nodes)
  cntr = 0
  for(node1 in 1:(no_nodes - 1)) {
    for(node2 in (node1 + 1):no_nodes) {
      cntr = cntr + 1
      proposal_sd[node1, node2] = sqrt(-1 / hessian[cntr, cntr])
      proposal_sd[node2, node1] = proposal_sd[node1, node2]
    }
  }

  # # Starting value of model matrix:
  gamma = matrix(1,
                 nrow = no_nodes,
                 ncol = no_nodes)

  #Starting values of interactions and thresholds (posterior mode)
  interactions = pps$interactions
  thresholds = pps$thresholds

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
                      Index = Index,
                      iter = iter,
                      burnin = burnin,
                      n_cat_obs = n_cat_obs,
                      threshold_alpha = threshold_alpha,
                      threshold_beta = threshold_beta,
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

    return(list(gamma = gamma,
                interactions = interactions,
                thresholds = tresholds))
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

    return(list(gamma = gamma,
                interactions = interactions,
                thresholds = thresholds))
  }
}


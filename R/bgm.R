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
#' @param x An \code{n} by \code{p} matrix containing the binary and ordinal 
#' variables for \code{n} independent observations on \code{p} variables in the 
#' network or graph. If not done already, \code{bgms} recodes the variables as 
#' non-negative integers (i.e., \code{0, 1, ..., m}). Unobserved categories are 
#' collapsed into other categories after recoding. See \code{reformat_data} for
#' details.
#' 
#' @param iter The number of iterations of the Gibbs sampler. Defaults to 
#' \code{1e4}. This usually gives a good idea. But for good estimates it is 
#' recommended to run the procedure for \code{1e5} iterations. 
#' 
#' @param interaction_prior The prior distribution for the interaction effects. 
#' Currently, two prior densities are implemented: The Unit Information prior
#' (\code{interaction_prior = "UnitInfo"}) and the Cauchy prior 
#' (\code{interaction_prior = "Cauchy"}). Defaults to \code{"UnitInfo"}.
#' 
#' @param cauchy_scale The scale of the Cauchy prior for interactions. Defaults 
#' to \code{2.5}. 
#' 
#' @param threshold_alpha,threshold_beta The shape parameters of the Beta-prime 
#' prior for the thresholds. Defaults to \code{1}.
#' 
#' @param save Should the function collect and return all samples from the 
#' Gibbs sampler (\code{save = TRUE})? Or should it only return the 
#' (model-averaged) posterior means (\code{save = FALSE})? Defaults to 
#' \code{FALSE}.
#' 
#' @param caching Some terms in the pseudolikelihood need to be evaluated often,
#' yet their values do not necessarily change, or in a predictable way. If 
#' \code{caching = TRUE} these terms are cached in an \code{n} by \code{p} 
#' matrix of real numbers. This matrix is computed once and updated only when 
#' its values change. This decreases the run time of the algorithm, but requires 
#' more memory. If \code{caching = FALSE}, the procedure does not cache terms. 
#' Defaults to \code{TRUE}. 
#'
#' @param display_progress Should the function show a progress bar 
#' (\code{display_progress = TRUE})? Or not (\code{display_progress = FALSE})?
#' Defauls to \code{TRUE}.
#' 
#' @param burnin The number of burnin iterations. The output of the Gibbs 
#' sampler is stored after \code{burnin} iterations.   
#' 
#' @return If \code{save = FALSE} (the default), a list containing the 
#' \code{p} by \code{p} matrices \code{gamma} and \code{interactions}, and the 
#' \code{p} by \code{max(no_categories)} matrix \code{thresholds}. The matrix 
#' \code{gamma} is a numeric matrix that contains the inclusion probabilities 
#' for individual edges. The matrices \code{interactions} and \code{thresholds} 
#' are numeric matrices that contain the (model or structure-averaged) posterior 
#' means (EAP estimates) of the pairwise associations, and category thresholds, 
#' respectively. If \code{save = TRUE}, a list containing the 
#' \code{iter} by \code{p *  (p - 1) / 2} matrices \code{samples.gamma} 
#' and \code{samples.interactions}, and the \code{iter} by 
#' \code{sum(no_categories)} matrix \code{samples.thresholds}. These contain the 
#' parameter states at every iteration of the Gibbs sampler. Column averages 
#' offer the EAP estimates.
bgm = function(x,
               iter = 1e4,
               burnin = 1e3,
               interaction_prior = c("UnitInfo", "Cauchy"),
               cauchy_scale = 2.5,
               threshold_alpha = 1,
               threshold_beta = 1,
               save = FALSE,
               caching = TRUE,
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
    stop("We use normal approximations to the posterior as proposal distribution 
  in a Metropolis within Gibbs approach. The normal approximation is based on 
  the curvature around the mode, estimated from optimizing the full
  pseudoposterior. For your data the pseudoposterior could not be optimized. 
  Please check your data for missing categories, or low category counts. Please 
  contact the package author if the data checks out and the data do not explain 
  why optimization failed.")
  
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
                      caching = caching,
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
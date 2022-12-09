#' Bayesian structure learning in Markov Random Fields of mixed binary and 
#' ordinal variables using MCMC. 
#'
#' The function \code{bsl} explores the joint pseudoposterior distribution of 
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
#' network or graph. If not done already, \code{bsl} recodes the variables as 
#' non-negative integers (i.e., \code{0, 1, ..., m}). Unobserved categories are 
#' collapsed into other categories after recoding.
#' 
#' @param no_iterations The number of iterations of the Gibbs sampler.
#' 
#' @param interaction_prior The prior distribution for the interaction effects. 
#' Currently, two prior densities are implemented: The Unit Information prior
#' (\code{interaction_prior = "UnitInfo"}) and the Cauchy prior 
#' (\code{interaction_prior = "Cauchy"}). Defaults to \code{"Cauchy"}.
#' 
#' @param cauchy_scale The scale of the Cauchy prior for interactions. Defaults 
#' to \code{2.5}. 
#' 
#' @param threshold_alpha,threshold_beta The shape parameters of the Beta-prime prior for the 
#' thresholds. Defaults to \code{1}.
#' 
#' @param samples Should the function collect and return all samples from the 
#' Gibbs sampler (\code{samples = TRUE})? Or should it only return the 
#' (model-averaged) posterior means (\code{samples = FALSE})? Defauls to 
#' \code{FALSE}.
#' 
#' @param display_progress Should the function show a progress bar 
#' (\code{display_progress = TRUE})? Or not (\code{display_progress = FALSE})?
#' Defauls to \code{FALSE}.
#' 
#' @param moms_method In what way should we update the edge indicator-
#'  interaction pairs with the metropolis algorithm. The add-delete scheme 
#'  tries to update one pair, followed by an update of all interactions that are 
#'  included at the time \code{"AddDelete"}. The Gibbs scheme tries to update 
#'  every pair in turn \code{"Gibbs"}. Finally, we can follow up the Gibbs 
#'  scheme with an update of all interactions that are included at the time 
#'  \code{"GibbsFull"}. Defaults to \code{"GibbsFull"}. 
#' 
#' @return If \code{samples = FALSE} (the default), a list containing the 
#' \code{p} by \code{p} matrices \code{gamma} and \code{interactions}, and the 
#' \code{p} by \code{max(no_categories)} matrix \code{thresholds}. The matrix 
#' \code{gamma} is a numeric matrix that contains the inclusion probabilities 
#' for individual edges. The matrices \code{interactions} and \code{thresholds} 
#' are numeric matrices that contain the (model or structure-averaged) posterior 
#' means (EAP estimates) of the pairwise associations, and category thresholds, 
#' respectively. If \code{samples = TRUE}, a list containing the 
#' \code{no_iterations} by \code{p *  (p - 1) / 2} matrices \code{samples.gamma} 
#' and \code{samples.interactions}, and the \code{no_iterations} by 
#' \code{sum(no_categories)} matrix \code{samples.thresholds}. These contain the 
#' parameter states at every iteration of the Gibbs sampler. Column averages 
#' offer the EAP estimates.
bsl = function(x,
               no_iterations = 1e5,
               interaction_prior = "Cauchy",
               cauchy_scale = 2.5,
               threshold_alpha = 1,
               threshold_beta = 1,
               samples = FALSE,
               display_progress = FALSE,
               moms_method = "GibbsFull") {
  
  #Check prior set-up for the interaction parameters ---------------------------
  if(interaction_prior != "Cauchy" & interaction_prior != "UnitInfo")
    stop("For the interaction effects we currently only have implemented the 
     Cauchy prior and the Unit Information prior. Please select one.")
  if(interaction_prior == "Cauchy") {
    if(cauchy_scale <= 0)
      stop("The scale of the Cauchy prior needs to be positive.")
  }  
  
  #Check prior set-up for the threshold parameters -----------------------------
  if(threshold_alpha <= 0  | !is.finite(threshold_alpha))
    stop("Parameter ``threshold_alpha'' needs to be positive.")
  if(threshold_beta <= 0  | !is.finite(threshold_beta))
    stop("Parameter ``threshold_beta'' needs to be positive.")
  
  #Check set-up for the mixtures of mutually singular distributions metropolis -
  if(!(moms_method %in% c("Gibbs", "GibbsFull", "AddDelete")))
    stop("The metropolis method, i.e., ``moms_method'', should be set to either
    ``Gibbs'', ``GibbsFull'', or ``AddDelete''.")
  
  #Check data input ------------------------------------------------------------
  if(class(x)[1] != "matrix") {
    stop("The input x is supposed to be a matrix.")
  }
  if(ncol(x) < 2)
    stop("The matrix x should have more than one variable (columns).")
  if(nrow(x) < 2)
    stop("The matrix x should have more than one observation (rows).")
  
  no_nodes = ncol(x)
  no_interactions = no_nodes * (no_nodes - 1) / 2
  
  #Format the data input -------------------------------------------------------
  no_categories = vector(length = no_nodes)
  for(node in 1:no_nodes) {
    unq_vls = sort(unique(x[,  node]))
    mx_vl = max(unq_vls)
    # Check
    if(mx_vl == nrow(x))
      stop(paste0("Only unique category responses observed for variable ", 
                  node, 
                  ". Expect discrete data with >= 1 observations per category."))
    if(length(unq_vls) != mx_vl + 1 || any(unq_vls != 0:mx_vl)) {
      y = x[, node]
      cntr = 0
      for(value in unq_vls) {
        x[y == value, node] = cntr
        cntr = cntr + 1
      }
    }
    no_categories[node] = max(x[,node])
    if(no_categories[node] == 0)
      stop(paste0("Only one value [", 
                  unq_vls,  
                  "] was observed for variable ", 
                  node, 
                  "."))
  }
  no_thresholds = sum(no_categories)
  
  #Proposal set-up for the interaction parameters ------------------------------
  pps = try(mppe(x = x, 
                 no_categories = no_categories, 
                 interaction_prior =  interaction_prior), 
            silent = TRUE)
  
  if(class(pps)[1] == "try-error")
    stop("We use normal approximations to the posterior as proposal distribution 
  in a Metropolis within Gibbs approach. The normal approximation is based on 
  the mode and curvature around the mode, estimated from optimizing the full
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
  
  # Starting value of model matrix:
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
  if(samples == FALSE) {
    out = gibbs_eap(observations = x,
                    no_categories  = no_categories,
                    interaction_prior = interaction_prior,
                    cauchy_scale = cauchy_scale,
                    unit_info = unit_info,
                    proposal_sd = proposal_sd,
                    Index = Index,
                    no_iterations = no_iterations,
                    n_cat_obs = n_cat_obs, 
                    threshold_alpha = threshold_alpha,
                    threshold_beta = threshold_beta,
                    display_progress = display_progress,
                    moms_method = moms_method)
    
    eap_edges = out$eap.edges
    eap_tresholds = out$eap.thresholds
    
    eap_interactions = matrix(0, 
                              nrow = no_nodes, 
                              ncol = no_nodes)
    eap_interactions[lower.tri(eap_interactions)] = 
      eap_edges[lower.tri(eap_edges)]
    eap_interactions = eap_interactions + t(eap_interactions) 
    
    eap_gamma = matrix(0, 
                       nrow = no_nodes, 
                       ncol = no_nodes)
    eap_gamma[upper.tri(eap_gamma)] = 
      eap_edges[upper.tri(eap_edges)]
    eap_gamma = eap_gamma + t(eap_gamma) 
    
    colnames(eap_interactions) = paste0("node ", 1:no_nodes)
    rownames(eap_interactions) = paste0("node ", 1:no_nodes)
    colnames(eap_gamma) = paste0("node ", 1:no_nodes)
    rownames(eap_gamma) = paste0("node ", 1:no_nodes)
    colnames(eap_tresholds) = paste0("category ", 1:max(no_categories))
    rownames(eap_tresholds) = paste0("node ", 1:no_nodes)
    
    return(list(gamma = eap_gamma, 
                interactions = eap_interactions,
                thresholds = eap_tresholds))
  } else {
    out = gibbs_samples(observations = x,
                        no_categories  = no_categories,
                        interaction_prior = interaction_prior,
                        cauchy_scale = cauchy_scale,
                        unit_info = unit_info,
                        proposal_sd = proposal_sd,
                        Index = Index,
                        no_iterations = no_iterations,
                        n_cat_obs = n_cat_obs, 
                        threshold_alpha = threshold_alpha,
                        threshold_beta = threshold_beta,
                        display_progress = display_progress,
                        moms_method = moms_method)
    
    samples_int = out$samples.interactions
    samples_gamma = out$samples.gamma
    samples_tre = out$samples.thresholds
    
    names1 = names2 = character(length = no_nodes * (no_nodes - 1) / 2)
    cntr = 0
    for(node in 1:(no_nodes - 1)) {
      for(node_2 in (node + 1):no_nodes) {
        cntr = cntr + 1
        names1[cntr] = paste0("sigma(",node, ", ",node_2,")")
        names2[cntr] = paste0("gamma(",node, ", ",node_2,")")
      }
    }
    colnames(samples_gamma) = names2
    colnames(samples_int) = names1
    
    names = character(length = sum(no_categories))
    cntr = 0
    for(node in 1:no_nodes) {
      for(category in 1:no_categories[node]) {
        cntr = cntr + 1
        names[cntr] = paste0("threshold(",node, ", ",category,")")
      }
    }
    colnames(samples_tre) = names
    
    rownames(samples_tre) = paste0("Iter. ", 1:no_iterations)
    rownames(samples_int) = paste0("Iter. ", 1:no_iterations)
    rownames(samples_gamma) = paste0("Iter. ", 1:no_iterations)
    
    return(list(samples.gamma = samples_gamma, 
                samples.interactions = samples_int,
                samples.thresholds = samples_tre))
  }
}
#' EM edge selection for a Markov Random Field model for ordinal variables. 
#'
#' The function \code{emes} selects promising edges for the ordinal 
#' MRF using the joint pseudolikelihood and a continuous spike and slab prior 
#' distribution stipulated on the MRF's interaction or association parameters.
#'
#' @param x An \code{no_persons} by \code{no_nodes} matrix containing the 
#'   categories coded as non-negative integers (i.e., coded 
#'   \code{0, 1, ..., no_categories}) for \code{no_persons} independent 
#'   observations on \code{no_nodes} variables in the network or graph.
#' 
#' @param precision A number between zero and one. The prior precision that is
#' desired for edge selection. Equal to one minus the desired type-1 error.
#' Defaults to \code{.975}.
#' 
#' @param convergence_criterion The convergence criterion for the 
#' pseudoposterior values in the EM algorithm. Defaults to 
#' \code{sqrt(.Machine$double.eps)}.
#'
#' @param theta The prior inclusion probability. The value \code{theta = 0.5},
#'   in combination with \code{hierarchical = FALSE} stipulates a uniform prior
#'   on the space of network structures.
#'
#' @param hierarchical Logical. If TRUE, a beta prior distribution is
#'  imposed on the prior inclusion probability \code{theta} with
#'  hyperparameters \code{alpha} and \code{beta}. A uniform prior on the
#'  inclusion probability, a beta with \code{alpha = beta = 1}, stipulates a
#'  uniform prior on network structure complexity.
#'
#' @param indicator_alpha,indicator_beta The hyperparameters of the beta prior 
#'  distribution stipulated on the prior inclusion probability \code{theta} if
#'  \code{hierarchical = TRUE}. Default to \code{1}.
#'   
#' @param maximum_iterations The maximum number of EM iterations used. Defaults 
#'   to \code{1e3}. A warning is issued if procedure has not converged in 
#'   \code{maximum_iterations} iterations.
#'
#' @param threshold_alpha,threshold_beta The shape parameters of the Beta-prime 
#'  prior for the thresholds. Defaults to \code{1}.
#' 
#' @return A list containing the \code{no_nodes} by \code{no_nodes} matrices 
#'  \code{interactions} and \code{gamma}, the \code{no_nodes} by 
#'  \code{no_categories} matrix \code{thresholds}, and, if 
#'  \code{hierarchical == TRUE}, a numeric valued \code{theta}. The matrix 
#'  \code{interactions} is a numeric matrix with pairwise
#'   association estimates on the off-diagonal elements. The matrix \code{gamma} 
#'   contains the expected values of edge inclusion variables (i.e., the local 
#'   posterior probability of edge inclusion. The matrix \code{thresholds} 
#'   contains the category thresholds per node. If \code{hierarchical = TRUE}, 
#'   the modal estimate of the prior inclusion probability \code{theta} is also 
#'   provided.
emes = function(x, 
                precision = 0.975,
                convergence_criterion = sqrt(.Machine$double.eps), 
                theta = 0.5, 
                hierarchical = FALSE, 
                indicator_alpha = 1, 
                indicator_beta = 1, 
                maximum_iterations = 1e3,
                threshold_alpha = 1,
                threshold_beta = 1) {
  
  #Check prior set-up for the interaction parameters ---------------------------
  if(precision < 0 || precision > 1)
    stop("The precision parameter needs to be between 0 and 1.")
  
  #Check prior set-up for the indicator variables ------------------------------
  if(theta < 0 || theta > 1) 
    stop("Parameter ``theta''is a probability and needs to be between 0 and 1.")
  if(indicator_alpha <= 0 | !is.finite(indicator_alpha))
    stop("Parameter ``indicator_alpha'' needs to be positive.")
  if(indicator_beta <= 0 | !is.finite(indicator_beta))
    stop("Parameter ``indicator_beta'' needs to be positive.")
  
  #Check prior set-up for the threshold parameters -----------------------------
  if(threshold_alpha <= 0 | !is.finite(threshold_alpha))
    stop("Parameter ``threshold_alpha'' needs to be positive.")
  if(threshold_beta <= 0 | !is.finite(threshold_beta))
    stop("Parameter ``threshold_beta'' needs to be positive.")
  
  #Check EM input --------------------------------------------------------------
  if(convergence_criterion <= 0) 
    stop("Parameter ``convergence_criterion'' needs to be positive.")
  if(maximum_iterations <= 0 || 
     abs(maximum_iterations - round(maximum_iterations)) > sqrt(.Machine$double.eps)) 
    stop("Parameter ``maximum_iterations'' needs to be a positive integer.")
  
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
  
  # Data reformatting
  no_categories = vector(length = no_nodes)
  for(node in 1:no_nodes) {
    unq_vls = sort(unique(x[,  node]))
    mx_vl = max(unq_vls)
    # Check
    if(mx_vl == nrow(x))
      stop(paste0("Only unique values observed for variable ", 
                  node, 
                  ". Expect discrete data with >= 1 observations per unique value."))
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
  no_parameters = no_thresholds + no_interactions
  no_persons = nrow(x)
  
  # Set spike and slab prior variances -----------------------------------------
  fit <- try(mple(x = x, no_categories = no_categories), 
             silent = TRUE)
  if(class(fit) != "try-error") {
    thresholds <- fit$thresholds
    interactions <- fit$interactions
  } else {
    stop("Pseudolikelihood optimization failed. Please check the data. If the 
         data checks out, please try different starting values.")
  }
  
  xi <- uniroot (f = xi_delta_matching,
                 interval = c(.Machine$double.eps,
                              no_persons - sqrt(.Machine$double.eps)),
                 delta = qnorm(precision, lower.tail = TRUE),
                 n = no_persons)$root
  
  slab_var <- set_slab(x = x, 
                       no_categories = no_categories, 
                       thresholds = thresholds, 
                       interactions = interactions)
  
  # EM ------------------------------------------------------------------------
  hessian <- matrix(data = NA, 
                    nrow = no_parameters,
                    ncol = no_parameters)
  gradient <- matrix(data = NA,
                     nrow = 1,
                     ncol = no_parameters)

  log_pseudoposterior <- 
    emvs_log_unnormalized_pseudoposterior(interactions = interactions,
                                          thresholds = thresholds,
                                          observations = x,
                                          no_categories = no_categories,
                                          xi = xi,
                                          slab_var = slab_var,
                                          theta = theta,
                                          hierarchical = hierarchical, 
                                          indicator_alpha = indicator_alpha, 
                                          indicator_beta = indicator_beta,
                                          threshold_alpha = threshold_alpha, 
                                          threshold_beta = threshold_beta)
  
  #starting values
  thresholds = matrix(0, 
                      nrow = no_nodes,
                      ncol = max(no_categories))
  interactions = matrix(0, 
                        nrow = no_nodes,
                        ncol = no_nodes)
  
  for(iteration in 1:maximum_iterations) {  
    old_log_pseudoposterior <- log_pseudoposterior
    
    # E-step - update selection variables -------------------------------------
    gamma <- em_gamma (interactions = interactions, 
                       slab_var = slab_var,
                       theta = theta,
                       xi = xi,
                       no_persons = no_persons)
    
    # M-step - update prior inclusion probability -----------------------------
    if(hierarchical == TRUE) {
      tmp <- sum(gamma[lower.tri(gamma)])
      theta <- (tmp + indicator_alpha - 1) / 
        (indicator_alpha + indicator_beta - 2 + no_interactions)
    }
    
    # M-step - update model parameters ----------------------------------------
    
    #Compute gradient vector --------------------------------------------------
    interaction_var <- em_interaction_var(gamma = gamma, 
                                          slab_var = slab_var,
                                          theta = theta,
                                          xi = xi,
                                          no_persons = no_persons)
    
    gradient[1:no_thresholds] <- 
      gradient_thresholds_pseudoposterior(interactions = interactions, 
                                          thresholds = thresholds, 
                                          observations = x, 
                                          no_categories = no_categories,
                                          threshold_alpha,
                                          threshold_beta)
    
    gradient[-c(1:no_thresholds)] <-
      gradient_interactions_pseudoposterior_normal(interactions = interactions, 
                                                   thresholds = thresholds, 
                                                   observations = x, 
                                                   no_categories = no_categories,
                                                   interaction_var = interaction_var)
    
    # Compute Hessian matrix (second order partial derivatives) ---------------
    hessian[1:no_thresholds, 1:no_thresholds] <- 
      hessian_thresholds_pseudoposterior(interactions = interactions, 
                                         thresholds = thresholds, 
                                         observations = x,
                                         no_categories = no_categories,
                                         threshold_alpha,
                                         threshold_beta)
    
    hessian[-(1:no_thresholds), -(1:no_thresholds)] <- 
      hessian_interactions_pseudoposterior_normal(interactions = interactions, 
                           thresholds = thresholds, 
                           observations = x, 
                           no_categories = no_categories,
                           interaction_var = interaction_var)

    hessian[-(1:no_thresholds), 1:no_thresholds] <- 
      hessian_crossparameters(interactions = interactions, 
                              thresholds = thresholds, 
                              observations = x,
                              no_categories = no_categories)
    
    hessian[1:no_thresholds, -(1:no_thresholds)] <- 
      t(hessian[-(1:no_thresholds), 1:no_thresholds])
    
    # Update parameter values (Newton-Raphson step) ---------------------------
    Delta <- gradient %*% solve(hessian)
    if(any(is.nan(Delta)) || any(is.infinite(Delta)))
      stop("Pseudoposterior optimization failed. Please check the data. If the 
           data checks out, please try different starting values.")
    
    cntr = 0
    for(node in 1:no_nodes) {
      for(category in 1:no_categories[node]) {
        cntr = cntr + 1
        thresholds[node, category] = thresholds[node, category] - Delta[cntr]
      }
    }
    for(node in 1:(no_nodes - 1)) {
      for(node_2 in (node + 1):no_nodes) {
        cntr = cntr + 1
        interactions[node, node_2] = interactions[node, node_2] - Delta[cntr]
        interactions[node_2, node] = interactions[node, node_2]
      }
    }
    
    # recompute log-pseudoposterior -------------------------------------------
    log_pseudoposterior <- 
      emvs_log_unnormalized_pseudoposterior(interactions = interactions,
                                            thresholds = thresholds,
                                            observations = x,
                                            no_categories = no_categories,
                                            xi = xi,
                                            slab_var = slab_var,
                                            theta = theta,
                                            hierarchical = hierarchical, 
                                            indicator_alpha = indicator_alpha, 
                                            indicator_beta = indicator_beta,
                                            threshold_alpha = threshold_alpha, 
                                            threshold_beta = threshold_beta)
    
    if(abs(log_pseudoposterior - old_log_pseudoposterior) < 
       convergence_criterion)
      break
  }
  
  if(abs(log_pseudoposterior - old_log_pseudoposterior) >= 
     convergence_criterion && 
     iteration == maximum_iterations)
    warning(paste("The optimization procedure did not convergence in", 
                  maximum_iterations, "iterations.",
                  sep = " "), 
            call. = FALSE)
  
  colnames(interactions) = paste0("node ", 1:no_nodes)
  rownames(interactions) = paste0("node ", 1:no_nodes)
  colnames(gamma) = paste0("node ", 1:no_nodes)
  rownames(gamma) = paste0("node ", 1:no_nodes)
  colnames(thresholds) = paste0("category ", 1:max(no_categories))
  rownames(thresholds) = paste0("node ", 1:no_nodes)
  
  if(!hierarchical)
    return(list(interactions = interactions, 
                thresholds = thresholds, 
                gamma = gamma))
  return(list(interactions = interactions, 
              thresholds = thresholds, 
              gamma = gamma, 
              theta = theta))
}
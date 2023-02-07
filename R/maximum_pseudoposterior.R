#' Optimization of the Pseudoposterior for a Markov Random Field model for 
#' ordinal variables. 
#'
#' The function \code{mppe} estimates the parameters for the ordinal MRF 
#' by optimizing the pseudoposterior with Newton-Raphson.
#'
#' @param x An \code{no_persons} by \code{no_nodes} matrix containing the 
#' categories coded as non-negative integers (i.e., coded 
#' \code{0, 1, ..., no_categories}) for \code{no_persons} independent 
#' observations on \code{no_nodes} variables in the network or graph.
#'
#' @param no_categories The maximum category.
#' 
#' @param interaction_prior The prior distribution for the interaction effects. 
#' Currently, two prior densities are implemented: The Unit Information prior
#' (\code{interaction_prior = "UnitInfo"}) and the Cauchy prior 
#' (\code{interaction_prior = "Cauchy"}). Defaults to \code{"Cauchy"}.
#' 
#' @param cauchy_scale The scale of the Cauchy prior for interactions. Defaults 
#' to \code{2.5}. 
#' 
#' @param threshold_alpha,threshold_beta The shape parameters of the Beta-prime 
#' prior for the thresholds. Default to \code{1}.
#' 
#' @param precision A number between zero and one. The prior precision that is
#' desired for edge selection. Equal to one minus the desired type-1 error.
#' Defaults to \code{.975}.
#' 
#' @param convergence_criterion The convergence criterion for the 
#' pseudoposterior values in the EM algorithm. Defaults to 
#' \code{sqrt(.Machine$double.eps)}.
#'   
#' @param maximum_iterations The maximum number of EM iterations used. Defaults 
#' to \code{1e3}. A warning is issued if procedure has not converged in 
#' \code{maximum_iterations} iterations.
#'
#' @param thresholds A \code{no_nodes} by \code{no_categories} matrix 
#' \code{thresholds}. Used as starting values in the Newton-Raphson procedure. 
#' Optional.
#'
#' @param interactions A \code{no_nodes} by \code{no_nodes} matrices 
#' \code{interactions}. Used as starting values in the Newton-Raphson procedure. 
#' Optional. 
#' 
#' @return A list containing the \code{no_nodes} by \code{no_nodes} matrices 
#' \code{interactions}, the \code{no_nodes} by \code{no_categories} matrix 
#' \code{thresholds}, the \code{hessian} matrix as computed in the final 
#' iteration of the optimization procedure, and, in case 
#' \code{interaction_prior = "UnitInfo"}, the \code{no_nodes} by \code{no_nodes} 
#' matrix \code{unit_info}. The matrix \code{interactions} is a numeric matrix 
#' that contains the maximum pseudoposterior estimates of the pairwise 
#' associations. The matrix \code{thresholds} contains the maximum 
#' pseudoposterior estimates of the category thresholds parameters. The 
#' \code{hessian} matrix has dimensions equal to the number of thresholds + 
#' associations. The topleft square contains the thresholds, the bottomright
#' square the associations (of the form \code{(1,2), (1, 3), ..., (2, 1), ...}).
#' The \code{unit_information} matrix contains the variances of the unit 
#' information prior for the association effects. 
mppe = function(x, 
                no_categories, 
                interaction_prior = "Cauchy",
                cauchy_scale = 2.5,
                threshold_alpha = 1,
                threshold_beta = 1,
                convergence_criterion = sqrt(.Machine$double.eps),
                maximum_iterations = 1e3, 
                thresholds, 
                interactions) {

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
  no_persons = nrow(x)
  no_thresholds = sum(no_categories)
  no_interactions = no_nodes * (no_nodes - 1) / 2
  no_parameters = no_thresholds + no_interactions
  
  #Check NR input --------------------------------------------------------------
  if(convergence_criterion <= 0) 
    stop("Parameter ``convergence_criterion'' needs to be positive.")
  if(maximum_iterations <= 0 || 
     abs(maximum_iterations - round(maximum_iterations)) > 
     sqrt(.Machine$double.eps)) 
    stop("Parameter ``maximum_iterations'' needs to be a positive integer.")

  # Starting values -----------------------------------------------------------
  if(!hasArg("thresholds")) {
    thresholds = matrix(0, 
                        nrow = no_nodes, 
                        ncol = max(no_categories))
  }
  if(!hasArg("interactions")) {
    interactions = matrix(0, 
                          nrow = no_nodes, 
                          ncol = no_nodes)
  }
  
  # Newton-Raphson -------------------------------------------------------------
  if(interaction_prior == "UnitInfo") {
    # Maximum pseudolikelihood -------------------------------------------------
    mpl = try(mple(x = x, no_categories = no_categories), 
              silent = TRUE)
    if(inherits(mpl, what = "try-error"))
      stop("You have chosen the unit information prior. To set-up this prior, 
      the log-pseudolikelihood needs to be optimized. This failed for your data. 
      Please check your data for missing categories, or low category counts.
      Please contact the package author if the data checks out and the data do 
      not explain why optimization failed.")
    
    # Asymptotic covariance ----------------------------------------------------
    hessian =  hessian_interactions_pseudolikelihood(interactions = mpl$interactions, 
                                                     thresholds = mpl$thresholds, 
                                                     observations = x,
                                                     no_categories)
    pr_var = no_persons * diag(-solve(hessian))
    unit_info = matrix(0, 
                       nrow = no_nodes, 
                       ncol = no_nodes)
    cntr = 0
    for(node1 in 1:(no_nodes - 1)) {
      for(node2 in (node1 + 1):no_nodes) {
        cntr = cntr + 1
        unit_info[node1, node2] = pr_var[cntr]
        unit_info[node2, node1] = unit_info[node1, node2]
      }
    }
    
    log_pseudoposterior = 
      log_unnormalized_pseudoposterior_normal(interactions, 
                                              thresholds, 
                                              observations = x,
                                              no_categories, 
                                              interaction_var = unit_info,
                                              threshold_alpha,
                                              threshold_beta)
    
    hessian = matrix(data = NA, 
                     nrow = no_parameters,
                     ncol = no_parameters)
    gradient = matrix(data = NA,
                      nrow = 1,
                      ncol = no_parameters)
    
    for(iteration in 1:maximum_iterations) {  
      old_log_pseudoposterior = log_pseudoposterior
      
      #Compute gradient vector (first order derivatives) ------------------------
      gradient[1:no_thresholds] = 
        gradient_thresholds_pseudoposterior(interactions = interactions, 
                                            thresholds = thresholds, 
                                            observations = x,
                                            no_categories,
                                            threshold_alpha,
                                            threshold_beta)
      
      gradient[-c(1:no_thresholds)] =
        gradient_interactions_pseudoposterior_normal(interactions = interactions, 
                                                     thresholds = thresholds, 
                                                     observations = x,
                                                     no_categories, 
                                                     interaction_var = unit_info)
      
      # Compute Hessian matrix (second order partial derivatives) ---------------
      hessian[1:no_thresholds, 1:no_thresholds] = 
        hessian_thresholds_pseudoposterior(interactions = interactions, 
                                           thresholds = thresholds, 
                                           observations = x,
                                           no_categories,
                                           threshold_alpha,
                                           threshold_beta)
      
      hessian[-(1:no_thresholds), -(1:no_thresholds)] = 
        hessian_interactions_pseudoposterior_normal(interactions = interactions, 
                                                    thresholds = thresholds, 
                                                    observations = x,
                                                    no_categories, 
                                                    interaction_var = unit_info)

      hessian[-(1:no_thresholds), 1:no_thresholds] = 
        hessian_crossparameters(interactions = interactions, 
                                thresholds = thresholds, 
                                observations = x,
                                no_categories)
      
      hessian[1:no_thresholds, -(1:no_thresholds)] = 
        t(hessian[-(1:no_thresholds), 1:no_thresholds])
      
      # Update parameter values (Newton-Raphson step) ---------------------------
      Delta = gradient %*% solve(hessian)
      if(any(is.nan(Delta)) || any(is.infinite(Delta)))
        stop("log_pseudoposterior optimization failed. Please check the data. 
             If the data checks out, please try different starting values.")
      
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
      
      # Check for convergence ---------------------------------------------------
      log_pseudoposterior = 
        log_unnormalized_pseudoposterior_normal(interactions, 
                                                thresholds, 
                                                observations = x,
                                                no_categories, 
                                                interaction_var = unit_info,
                                                threshold_alpha,
                                                threshold_beta)
      
      if(abs(log_pseudoposterior- old_log_pseudoposterior) < 
         convergence_criterion)
        break
    }  
  } else if (interaction_prior == "Cauchy") {
    log_pseudoposterior = 
      log_unnormalized_pseudoposterior_cauchy(interactions, 
                                              thresholds, 
                                              observations = x,
                                              no_categories, 
                                              cauchy_scale = cauchy_scale,
                                              threshold_alpha,
                                              threshold_beta)
    
    hessian = matrix(data = NA, 
                     nrow = no_parameters,
                     ncol = no_parameters)
    gradient = matrix(data = NA,
                      nrow = 1,
                      ncol = no_parameters)
    
    for(iteration in 1:maximum_iterations) {  
      old_log_pseudoposterior = log_pseudoposterior
      
      #Compute gradient vector (first order derivatives) -----------------------
      gradient[1:no_thresholds] = 
        gradient_thresholds_pseudoposterior(interactions = interactions, 
                                            thresholds = thresholds, 
                                            observations = x,
                                            no_categories,
                                            threshold_alpha,
                                            threshold_beta)
      
      gradient[-c(1:no_thresholds)] =
        gradient_interactions_pseudoposterior_cauchy(interactions = interactions, 
                                                     thresholds = thresholds, 
                                                     observations = x,
                                                     no_categories, 
                                                     cauchy_scale = cauchy_scale)
      
      # Compute Hessian matrix (second order partial derivatives) --------------
      hessian[1:no_thresholds, 1:no_thresholds] = 
        hessian_thresholds_pseudoposterior(interactions = interactions, 
                                           thresholds = thresholds, 
                                           observations = x,
                                           no_categories,
                                           threshold_alpha,
                                           threshold_beta)
      
      hessian[-(1:no_thresholds), -(1:no_thresholds)] = 
        hessian_interactions_pseudoposterior_cauchy(interactions = interactions, 
                                                    thresholds = thresholds, 
                                                    observations = x,
                                                    no_categories, 
                                                    cauchy_scale = cauchy_scale)
      
      hessian[-(1:no_thresholds), 1:no_thresholds] = 
        hessian_crossparameters(interactions = interactions, 
                                thresholds = thresholds, 
                                observations = x,
                                no_categories)
      
      hessian[1:no_thresholds, -(1:no_thresholds)] = 
        t(hessian[-(1:no_thresholds), 1:no_thresholds])
      
      # Update parameter values (Newton-Raphson step) --------------------------
      Delta = gradient %*% solve(hessian)
      if(any(is.nan(Delta)) || any(is.infinite(Delta)))
        stop("log_pseudoposterior optimization failed. Please check the data. 
             If the data checks out, please try different starting values.")
      
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
      
      # Check for convergence ---------------------------------------------------
      log_pseudoposterior = 
        log_unnormalized_pseudoposterior_cauchy(interactions, 
                                                thresholds, 
                                                observations = x,
                                                no_categories, 
                                                cauchy_scale = cauchy_scale,
                                                threshold_alpha,
                                                threshold_beta)
      
      if(abs(log_pseudoposterior - old_log_pseudoposterior) < 
         convergence_criterion)
        break
    }
  }
  
  if(abs(log_pseudoposterior - old_log_pseudoposterior) >= 
     convergence_criterion && iteration == maximum_iterations)
    warning(paste("The optimization procedure did not convergence in", 
                  maximum_iterations, "iterations.",
                  sep = " "), 
            call. = FALSE)
  
  colnames(interactions) = paste0("node ", 1:no_nodes)
  rownames(interactions) = paste0("node ", 1:no_nodes)
  colnames(thresholds) = paste0("category ", 1:max(no_categories))
  rownames(thresholds) = paste0("node ", 1:no_nodes)
  
  names = character(length = no_nodes * (no_nodes - 1) / 2 + 
                      sum(no_categories))
  cntr = 0
  for(node in 1:no_nodes) {
    for(category in 1:no_categories[node]) {
      cntr = cntr + 1
      names[cntr] = paste0("threshold(", node, ", ",category,")")
    }
  }
  for(node in 1:(no_nodes - 1)) {
    for(node_2 in (node + 1):no_nodes) {
      cntr = cntr + 1
      names[cntr] = paste0("sigma(", node, ", ",node_2,")")
    }
  }
  rownames(hessian) = names
  colnames(hessian) = names
  
  if(interaction_prior == "UnitInfo"){
    colnames(unit_info) = paste0("node ", 1:no_nodes)
    rownames(unit_info) = paste0("node ", 1:no_nodes)
    return(list(interactions = interactions, 
                thresholds = thresholds,
                hessian = hessian, 
                unit_info = unit_info))
  } 
  return(list(interactions = interactions, 
              thresholds = thresholds,
              hessian = hessian))
}
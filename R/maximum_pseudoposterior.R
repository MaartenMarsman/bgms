#' Optimize Pseudoposterior for an Ordinal Markov Random Field Model
#'
#' The function \code{mppe} estimates the parameters for the ordinal MRF
#' by optimizing the pseudoposterior with the Newton-Raphson method.
#'
#' @param x A dataframe or matrix with \code{n} rows and \code{p} columns,
#' containing binary and ordinal variables for \code{n} independent observations
#' and \code{p} variables in the network. Variables are recoded as non-negative
#' integers \code{(0, 1, ..., m)} if not done already. Unobserved categories are
#' collapsed into other categories after recoding. See \code{reformat_data} for
#' details.
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
#' @param convergence_criterion The convergence criterion for the
#' pseudoposterior values in the EM algorithm. Defaults to
#' \code{sqrt(.Machine$double.eps)}.
#'
#' @param maximum_iterations The maximum number of EM iterations used. Defaults
#' to \code{1e3}. A warning is issued if the procedure has not converged in
#' \code{maximum_iterations} iterations.
#'
#' @param thresholds A matrix with \code{p} rows and \code{max(m)} columns,
#' containing the category thresholds for each node. Used as starting values in
#' the Newton-Raphson procedure. Optional.
#'
#' @param interactions A matrix with \code{p} rows and \code{p} columns,
#' containing the pairwise association estimates in the off-diagonal elements.
#' Used as starting values in the Newton-Raphson procedure. Optional.
#'
#' @return A list containing:
#' \itemize{
#' \item \code{interactions}: A matrix with \code{p} rows and \code{p} columns,
#' containing the maximum pseudoposterior estimates of the pairwise
#' associations in the off-diagonal elements.
#' \item \code{thresholds}: A matrix with \code{p} rows and \code{max(m)}
#' columns, containing the maximum pseudoposterior estimates of the category
#' thresholds for each node.
#' \item \code{hessian}: A square matrix with \code{sum(m) + p(p-1)/2} rows and
#' columns, evaluated at the maximum pseudoposterior estimates. The top-left
#' square contains the thresholds, the bottom-right square the associations (of
#' the form \code{(1,2), (1, 3), ..., (2, 1), ...}).
#' }
#' In the case that \code{interaction_prior = "UnitInfo"}, the list also
#' contains the \code{p} by \code{p} matrix \code{unit_info}, which contains the
#' asymptotic variances that are used to set the unit information prior for the
#' association effects in the \code{bgms} function.
#'
#' @export
mppe = function(x,
                interaction_prior = c("Cauchy", "UnitInfo"),
                cauchy_scale = 2.5,
                threshold_alpha = 1,
                threshold_beta = 1,
                convergence_criterion = sqrt(.Machine$double.eps),
                maximum_iterations = 1e3,
                thresholds,
                interactions) {

  #Check data input ------------------------------------------------------------
  if(!inherits(x, what = "matrix") && !inherits(x, what = "data.frame"))
    stop("The input x needs to be a matrix or dataframe.")
  if(inherits(x, what = "data.frame"))
    x = data.matrix(x)
  if(ncol(x) < 2)
    stop("The matrix x should have more than one variable (columns).")
  if(nrow(x) < 2)
    stop("The matrix x should have more than one observation (rows).")

  #Check prior set-up for the interaction parameters ---------------------------
  interaction_prior = match.arg(interaction_prior)
  if(interaction_prior == "Cauchy") {
    if(cauchy_scale <= 0 || is.na(cauchy_scale) || is.infinite(cauchy_scale))
      stop("The scale of the Cauchy prior needs to be positive.")
  }
  if(interaction_prior == "Cauchy") {
    if(cauchy_scale <= 0)
      stop("The scale of the Cauchy prior needs to be positive.")
  }

  #Check prior set-up for the threshold parameters -----------------------------
  if(threshold_alpha <= 0  | !is.finite(threshold_alpha))
    stop("Parameter ``threshold_alpha'' needs to be positive.")
  if(threshold_beta <= 0  | !is.finite(threshold_beta))
    stop("Parameter ``threshold_beta'' needs to be positive.")

  #Format the data input -------------------------------------------------------
  data = reformat_data(x = x, fn.name = "mppe")
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
    mpl = try(mple(x = x), silent = TRUE)
    if(inherits(mpl, what = "try-error"))
      stop(paste0("You have chosen the unit information prior. To set-up this prior,\n",
                  "the log-pseudolikelihood needs to be optimized. This failed for \n",
                  "your data. Please switch to a different prior distribution."))

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

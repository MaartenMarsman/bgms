#' Maximum Pseudolikelihood Estimation for an Ordinal Markov Random Field Model
#'
#' The function \code{mple} estimates the parameters for the ordinal MRF
#' by optimizing the joint pseudolikelihood with the Newton-Raphson method.
#'
#' @param x A dataframe or matrix with \code{n} rows and \code{p} columns,
#' containing binary and ordinal variables for \code{n} independent observations
#' and \code{p} variables in the network. Variables are recoded as non-negative
#' integers \code{(0, 1, ..., m)} if not done already. Unobserved categories are
#' collapsed into other categories after recoding. See \code{reformat_data} for
#' details.
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
#' containing the maximum pseudolikelihood estimates of the pairwise
#' associations in the off-diagonal elements.
#' \item \code{thresholds}: A matrix with \code{p} rows and \code{max(m)}
#' columns, containing the maximum pseudolikelihood estimates of the category
#' thresholds for each node.
#' }
#' @export
mple = function(x,
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

  #Format the data input -------------------------------------------------------
  data = reformat_data(x = x, fn.name = "mple")
  x = data$x
  no_categories = data$no_categories
  no_nodes = ncol(x)
  no_interactions = no_nodes * (no_nodes - 1) / 2
  no_thresholds = sum(no_categories)
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

  # Newton-Raphson ------------------------------------------------------------
  log_pl = log_pseudolikelihood(interactions,
                                thresholds,
                                observations = x,
                                no_categories)

  hessian = matrix(data = NA,
                   nrow = no_parameters,
                   ncol = no_parameters)
  gradient = matrix(data = NA,
                    nrow = 1,
                    ncol = no_parameters)

  for(iteration in 1:maximum_iterations) {
    old_log_pl = log_pl

    #Compute gradient vector (first order derivatives) ------------------------
    gradient[1:no_thresholds] =
      gradient_thresholds_pseudolikelihood(interactions = interactions,
                                           thresholds = thresholds,
                                           observations = x,
                                           no_categories)
    gradient[-c(1:no_thresholds)] =
      gradient_interactions_pseudolikelihood(interactions = interactions,
                                             thresholds = thresholds,
                                             observations = x,
                                             no_categories)

    # Compute Hessian matrix (second order partial derivatives) ---------------
    hessian[1:no_thresholds, 1:no_thresholds] =
      hessian_thresholds_pseudolikelihood(interactions = interactions,
                                          thresholds = thresholds,
                                          observations = x,
                                          no_categories)

    hessian[-(1:no_thresholds), -(1:no_thresholds)] =
      hessian_interactions_pseudolikelihood(interactions = interactions,
                                            thresholds = thresholds,
                                            observations = x,
                                            no_categories)

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
      stop("log_pseudolikelihood optimization failed. Please check the data. If the data checks out, please try different starting values.")

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
    log_pl = log_pseudolikelihood(interactions,
                                  thresholds,
                                  observations = x,
                                  no_categories)

    if(abs(log_pl - old_log_pl) < convergence_criterion)
      break
  }

  if(abs(log_pl - old_log_pl) >= convergence_criterion &&
     iteration == maximum_iterations)
    warning(paste("The optimization procedure did not convergence in",
                  maximum_iterations, "iterations.",
                  sep = " "),
            call. = FALSE)

  #Preparing the output --------------------------------------------------------
  if(is.null(colnames(x))){
    data_columnnames = paste0("node ", 1:no_nodes)
    colnames(interactions) = data_columnnames
    rownames(interactions) = data_columnnames
    rownames(thresholds) = data_columnnames
  } else {
    data_columnnames <- colnames(x)
    colnames(interactions) = data_columnnames
    rownames(interactions) = data_columnnames
    rownames(thresholds) = data_columnnames
  }
  colnames(thresholds) = paste0("category ", 1:max(no_categories))

  return(list(interactions = interactions, thresholds = thresholds))
}

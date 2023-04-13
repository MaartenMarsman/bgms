#' EM variable selection for a Markov Random Field model for ordinal variables.
#'
#' The function \code{bgm.em} selects promising edges for the ordinal
#' MRF using the joint pseudolikelihood and a continuous spike and slab prior
#' distribution stipulated on the MRF's interaction or association parameters.
#'
#' @param x A matrix with \code{n} rows and \code{p} columns, containing binary
#' and ordinal variables for \code{n} independent observations and \code{p}
#' variables in the network. Variables are recoded as non-negative integers
#' \code{(0, 1, ..., m)} if not done already. Unobserved categories are
#' collapsed into other categories after recoding. See \code{reformat_data} for
#' details.
#' @param precision A value between 0 and 1 representing the desired precision
#' for edge selection, equal to one minus the desired type-1 error rate. Default
#' is 0.975.
#' @param convergence_criterion The criterion for the pseudoposterior values'
#' convergence in the EM algorithm. Default is \code{sqrt(.Machine$double.eps)}.
#' @param theta The prior inclusion probability. A value of \code{0.5}, combined
#' with \code{hierarchical = FALSE}, specifies a uniform prior on the network
#' structures' space.
#' @param hierarchical If TRUE, a beta prior distribution with hyperparameters
#' \code{indicator_alpha, indicator_beta} is imposed on the prior inclusion
#' probability \code{theta}. A uniform prior on inclusion probability, using a
#' beta with \code{indicator_alpha = indicator_beta = 1}, specifies a uniform
#' prior on network structure complexity.
#' @param indicator_alpha,indicator_beta Hyperparameters for the beta prior
#' distribution on the prior inclusion probability theta when
#' \code{hierarchical = TRUE}. Default is 1.
#' @param maximum_iterations Maximum number of EM iterations used. Default is
#' 1e3. A warning appears if the procedure hasn't converged within the maximum
#' number of iterations.
#' @param threshold_alpha,threshold_beta Shape parameters for the Beta-prime
#' prior on thresholds. Default is 1.
#'
# @return A list containing:
#' \itemize{
#' \item \code{interactions}: A matrix with \code{p} rows and \code{p} columns,
#' containing the pairwise association estimates in the off-diagonal elements.
#' \item \code{gamma}: A matrix with \code{p} rows and \code{p} columns,
#' containing the expected values of edge inclusion variables (local posterior
#' probabilities of edge inclusion).
#' \item \code{thresholds}: A matrix with \code{p} rows and \code{max(m)}
#' columns, containing the category thresholds for each node.
#' \item \code{theta} (if \code{hierarchical == TRUE}): A numeric value
#' representing the modal estimate of the prior inclusion probability.
#' }
#'
#' @examples
#' \dontrun{
#'  ##Analyse the Wenchuan dataset
#'  fit = bgm.em(x = Wenchuan)
#'
#'
#'  #------------------------------------------------------------------------------|
#'  # INCLUSION - EDGE WEIGHT PLOT
#'  #------------------------------------------------------------------------------|
#'
#'  par(mar = c(6, 5, 1, 1))
#'  plot(x = fit$interactions[lower.tri(fit$interactions)],
#'       y = fit$gamma[lower.tri(fit$gamma)], ylim = c(0, 1),
#'       xlab = "", ylab = "", axes = FALSE, pch = 21, bg = "#bfbfbf", cex = 1.3)
#'  abline(h = 0, lty = 2, col = "#bfbfbf")
#'  abline(h = 1, lty = 2, col = "#bfbfbf")
#'  abline(h = .5, lty = 2, col = "#bfbfbf")
#'  mtext("Posterior Inclusion Probability", side = 1, line = 3, cex = 1.7)
#'  mtext("Posterior Mode Edge Weight", side = 2, line = 3, cex = 1.7)
#'  axis(1)
#'  axis(2, las = 1)
#'
#'
#'  #------------------------------------------------------------------------------|
#'  # THE LOCAL MEDIAN PROBABILITY NETWORK
#'  #------------------------------------------------------------------------------|
#'
#'  library(qgraph) #For plotting the estimated network
#'
#'  posterior.inclusion = fit$gamma[lower.tri(fit$gamma)]
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
#'  qgraph(median.prob.model,
#'         theme = "TeamFortress",
#'         maximum = .5,
#'         fade = FALSE,
#'         color = c("#f0ae0e"), vsize = 10, repulsion = .9,
#'         label.cex = 1.1, label.scale = "FALSE",
#'         labels = colnames(Wenchuan))
#' }
#' @export
bgm.em = function(x,
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
  no_parameters = no_interactions + no_thresholds
  no_persons = nrow(x)

  # Set spike and slab prior variances -----------------------------------------
  fit <- try(mple(x = x, no_categories = no_categories),
             silent = TRUE)
  if(inherits(fit, "try-error")) {
    stop(paste0(
      "We use a continuous spike and slab prior for the pairwise interactions.\n",
      "And use a function of the second derivatives of the MRF evaluated at the maximum\n",
      "pseudolikelihood estimates to set the variance of this prior distribution. We\n",
      "could not find this maximum for your data. Perhaps there were low category\n",
      "counts?"))
  } else {
    thresholds <- fit$thresholds
    interactions <- fit$interactions
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

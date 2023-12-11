#' Bayesian estimation of a Markov Random Field of mixed binary and ordinal
#' variables using MCMC.
#'
#' The function \code{est_bgm} explores the joint posterior distribution
#' of the parameters of a Markov Random Field for mixed binary and ordinal
#' variables. It uses a combination of Metropolis-Hastings and Gibbs sampling to
#' create a Markov chain that has the posterior distribution as its invariant
#' distribution. The function uses a pseudolikelihood approach to circumvent
#' having to compute the intractable normalizing constant of the Markov Random
#' Field.
#'
#' The prior distribution for the interactions is a Cauchy or Unit-Info type
#' prior. A Beta-prime distribution is used for the exponent of the category
#' parameters.
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
#' \item \code{interactions}: A matrix with \code{p} rows and \code{p} columns,
#' containing model-averaged posterior means of the pairwise associations.
#' \item \code{thresholds}: A matrix with \code{p} rows and \code{max(m)}
#' columns, containing model-averaged category thresholds.
#' }
#'
#' If \code{save = TRUE}, the result is a list of class ``bgms'' containing:
#' \itemize{
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
#' @export
est_bgm = function(x,
               iter = 1e4,
               burnin = 1e3,
               interaction_prior = c("Cauchy", "UnitInfo"),
               cauchy_scale = 2.5,
               threshold_alpha = 0.5,
               threshold_beta = 0.5,
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
  if(iter <= 0)
    stop("Parameter ``iter'' needs to be a positive integer.")
  if(abs(burnin - round(burnin)) > sqrt(.Machine$double.eps) || burnin < 0)
    stop("Parameter ``burnin'' needs to be a non-negative integer.")
  if(burnin <= 0)
    stop("Parameter ``burnin'' needs to be a positive integer.")

  #Check prior set-up for the interaction parameters ---------------------------
  interaction_prior = match.arg(interaction_prior)
  if(interaction_prior == "Cauchy") {
    if(cauchy_scale <= 0 || is.na(cauchy_scale) || is.infinite(cauchy_scale))
      stop("The scale of the Cauchy prior needs to be positive.")
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
    if(no_nodes == 2) {
      proposal_sd[1, 2] = sqrt(-1 / hessian[cntr])
      proposal_sd[2, 1] = proposal_sd[1, 2]
    } else {
      for(node1 in 1:(no_nodes - 1)) {
        for(node2 in (node1 + 1):no_nodes) {
          cntr = cntr + 1
          proposal_sd[node1, node2] = sqrt(-1 / hessian[cntr, cntr])
          proposal_sd[node2, node1] = proposal_sd[node1, node2]
        }
      }
    }
  }

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
  out = est_gibbs_sampler(observations = x,
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
                          na.impute,
                          missing_index,
                          adaptive = adaptive,
                          save = save,
                          display_progress = display_progress)

  #Preparing the output --------------------------------------------------------
  if(save == FALSE) {
    interactions = out$interactions
    tresholds = out$thresholds

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

    colnames(tresholds) = paste0("category ", 1:max(no_categories))

    output = list(interactions = interactions,
                  thresholds = thresholds,
                  save = save,
                  colnames = data_columnnames)
    class(output) = "bgms"
    return(output)
  } else {
    interactions = out$interactions
    thresholds = out$thresholds

    if(is.null(colnames(x))){
      data_columnnames <- 1:ncol(x)
    } else {
      data_columnnames <- colnames(x)
    }
    p <- ncol(x)
    names_bycol <- matrix(rep(data_columnnames, each = p), ncol = p)
    names_byrow <- matrix(rep(data_columnnames, each = p), ncol = p, byrow = T)
    names_comb <- matrix(paste0(names_byrow, "-", names_bycol), ncol = p)
    names_vec <- names_comb[lower.tri(names_comb)]

    colnames(interactions) = names_vec

    names = character(length = sum(no_categories))
    cntr = 0
    for(node in 1:no_nodes) {
      for(category in 1:no_categories[node]) {
        cntr = cntr + 1
        names[cntr] = paste0("threshold(",node, ", ",category,")")
      }
    }
    colnames(thresholds) = names

    rownames(interactions) = paste0("Iter. ", 1:iter)
    rownames(thresholds) = paste0("Iter. ", 1:iter)

    output = list(interactions = interactions,
                  thresholds = thresholds,
                  save = save,
                  colnames = data_columnnames)
    class(output) = "bgms"
    return(output)
  }
}
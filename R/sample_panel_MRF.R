#' Sample states of an ordinal cross-lagged panel MRF
#'
#' This function samples states from the ordinal cross-lagged panel MRF using a
#' Gibbs sampler. The Gibbs sampler is initiated with random values from the
#' response options, after which it proceeds by simulating states for each node
#' from a logistic model using the other node states as predictor variables.
#'
#' @param no_states The number of states of the ordinal cross-lagged panel MRF
#' to be generated. Must be a positive integer.
#'
#' @param no_nodes The number of nodes in the ordinal cross-lagged panel MRF.
#' Must be an integer greater than or equal to two.
#'
#' @param no_timepoints The number of timepoints in the ordinal cross-lagged
#' panel MRF. Must be a positive integer. The sampler also generates a vector of
#' responses for time \code{t = 0}, which are conditioned on at time \code{t = 1}.
#'
#' @param no_categories Either a positive integer or a vector of positive
#' integers of length \code{no_nodes}. The number of response categories on top
#' of the base category: \code{no_categories = 1} generates binary states.
#'
#' @param cross_sectional_interactions A symmetric \code{no_nodes} by
#' \code{no_nodes} matrix of pairwise cross-sectional interactions. Only its
#' off-diagonal elements are used.
#'
#' @param cross_lagged_interactions A \code{no_nodes} by \code{no_nodes} matrix
#' of pairwise cross-lagged interactions.
#'
#' @param thresholds A \code{no_nodes * no_timepoints} by
#' \code{max(no_categories)} matrix of category thresholds. The first
#' \code{no_nodes} rows correspond to the thresholds at time \code{t = 1}, the
#' second \code{no_nodes} rows correspond to the thresholds at time \code{t = 2},
#' etc. If \code{no_categories} is a vector, only the first
#' \code{no_categories[r]} elements are used in row \code{r}.
#'
#' @param null_interactions A \code{no_nodes} by \code{no_nodes} matrix
#' of pairwise interactions used to generate observations at time \code{t = 0}.
#'
#' @param null_thresholds A \code{no_nodes} by \code{max(no_categories)} matrix
#' of category thresholds used to generate observations at time \code{t = 0}. If
#' \code{no_categories} is a vector, only the first \code{no_categories[r]}
#' elements are used in row \code{r}.
#'
#' @param iter The number of iterations used by the Gibbs sampler.
#' The function provides the last state of the Gibbs sampler as output. By
#' default set to \code{1e3}.
#'
#' @return A \code{no_states} by \code{no_nodes * (no_timepoints + 1)} matrix of
#' simulated states of the ordinal cross-lagged panel MRF. The first
#' \code{no_nodes} columns correspond to the observations at time \code{t = 0},
#' the second \code{no_nodes} columns to the observations at time \code{t = 1},
#' etc.
#'
#' @examples
#' # Generate responses from a network of five ordinal variables on five timepoints.
#' no_states = 1000
#' no_nodes = 5
#' no_timepoints = 5
#' no_categories = 2
#'
#' u = runif(no_nodes, min = -.5, max = 1)
#' cross_sectional_interactions = .1 * u %*% t(u)
#' u = runif(no_nodes, min = -.5, max = 1)
#' cross_lagged_interactions = u %*% t(u)
#' thresholds = matrix(-abs(rnorm(no_timepoints * no_nodes * no_categories, 2)),
#'   nrow = no_timepoints * no_nodes, ncol = no_categories)
#'
#' null_interactions = 2 * cross_sectional_interactions
#' null_thresholds = thresholds[1:no_nodes, ]
#' x = panelmrfSampler (no_states = no_states,
#'   no_nodes = no_nodes,
#'   no_timepoints = no_timepoints,
#'   no_categories,
#'   cross_sectional_interactions,
#'   cross_lagged_interactions,
#'   thresholds,
#'   null_interactions,
#'   null_thresholds,
#'   iter = 1e4)

#' @export
panelmrfSampler = function(no_states,
                           no_nodes,
                           no_timepoints,
                           no_categories,
                           cross_sectional_interactions,
                           cross_lagged_interactions,
                           thresholds,
                           null_interactions,
                           null_thresholds,
                           iter = 1e3) {
  #check no_states, no_nodes, iter
  if(no_states <= 0 ||
     abs(no_states - round(no_states)) > .Machine$double.eps^.5)
    stop("``no_states'' must be a positive integer.")
  if(no_nodes <= 0 ||
     abs(no_nodes - round(no_nodes)) > .Machine$double.eps^.5)
    stop("``no_nodes'' must be a positive integer.")
  if(iter <= 0 ||
     abs(iter - round(iter)) > .Machine$double.eps^.5)
    stop("``iter'' must be a positive integer.")
  if(no_timepoints <= 0 ||
     abs(no_timepoints - round(no_timepoints)) > .Machine$double.eps^.5)
    stop("``no_timepoints'' must be a positive integer.")

  #check no_categories
  if(length(no_categories) == 1) {
    if(no_categories <= 0 ||
       abs(no_categories - round(no_categories)) > .Machine$double.eps^.5)
      stop("``no_categories'' must be a (vector of) positive integer(s).")
    no_categories = rep(no_categories, no_nodes)
  } else {
    for(node in 1:no_nodes) {
      if(no_categories[node] <= 0 ||
         abs(no_categories[node] - round(no_categories[node])) >
         .Machine$double.eps^.5)
        stop(paste("For node", node, "``no_categories'' must be a positive integer."))
    }
  }

  #check interactions
  if(!isSymmetric(cross_sectional_interactions))
    stop("The matrix ``cross_sectional_interactions'' must be symmetric.")
  diag(cross_sectional_interactions) = 0

  if(!isSymmetric(null_interactions))
    stop("The matrix ``null_interactions'' must be symmetric.")
  diag(null_interactions) = 0

  if(nrow(cross_sectional_interactions) != no_nodes)
    stop("The matrix ``cross_sectional_interactions'' must have ``no_nodes'' rows and columns.")

  if(nrow(cross_lagged_interactions) != no_nodes ||
     ncol(cross_lagged_interactions) != no_nodes)
    stop("The matrix ``cross_lagged_interactions'' must have ``no_nodes'' rows and columns.")

  if(nrow(null_interactions) != no_nodes)
    stop("The matrix ``null_interactions'' must have ``no_nodes'' rows and columns.")

  #check thresholds
  if(!inherits(thresholds, what = "matrix")) {
    if(max(no_categories) == 1) {
      if(length(thresholds) == no_nodes * no_timepoints) {
        thresholds = matrix(thresholds, ncol = 1)
      } else {
        stop(paste0("The matrix ``thresholds'' has ",
                    length(thresholds),
                    " elements, but requires",
                    no_nodes * no_timepoints,
                    "."))
      }
    } else {
      stop("Argument ``thresholds'' must be a matrix.")
    }
  }

  if(!inherits(null_thresholds, what = "matrix")) {
    if(max(no_categories) == 1) {
      if(length(null_thresholds) == no_nodes) {
        null_thresholds = matrix(null_thresholds, ncol = 1)
      } else {
        stop(paste0("The matrix ``null_thresholds'' has ",
                    length(null_thresholds),
                    " elements, but requires",
                    no_nodes,
                    "."))
      }
    } else {
      stop("Argument ``null_thresholds'' must be a matrix.")
    }
  }

  if(nrow(thresholds) != no_nodes * no_timepoints)
    stop("The matrix ``thresholds'' must have ``no_nodes * no_timepoints'' rows.")

  if(nrow(null_thresholds) != no_nodes)
    stop("The matrix ``null_thresholds'' must have ``no_nodes'' rows.")

  for(t in 1:no_timepoints) {
    start = (t - 1) * no_nodes + 1
    for(node in 1:no_nodes) {
      if(any(is.na(thresholds[node + start - 1, 1:no_categories[node]]))) {
        stop(paste0("The matrix ``thresholds'' contains NA(s) for node ",
                    node,
                    "at time",
                    t,
                    " in categorie(s)",
                    which(is.na(thresholds[node, 1:no_categories[node]])),
                    ", where a numeric value is needed."))
      }
      if(ncol(thresholds) > no_categories[node]) {
        if(any(!is.na(thresholds[node + start - 1, (no_categories[node]+1):ncol(thresholds)]))){
          warning(paste0("The matrix ``thresholds'' contains numeric values for node ",
                         node,
                         "at time",
                         t,
                         " for categories(s) (i.e., columns) exceding the maximum of ",
                         no_categories[node],
                         ". These values will be ignored."))
        }
      }
    }

    for(node in 1:no_nodes) {
      for(category in 1:no_categories[node]) {
        if(!is.finite(thresholds[node + start - 1, category]))
          stop(paste("The threshold parameter for node", node, "and category",
                     category, "at t =", t, "is NA or not finite."))
      }
    }
  }

  for(node in 1:no_nodes) {
    if(any(is.na(null_thresholds[node, 1:no_categories[node]]))) {
      stop(paste0("The matrix ``null_thresholds'' contains NA(s) for node ",
                  node,
                  " in categorie(s)",
                  which(is.na(null_thresholds[node, 1:no_categories[node]])),
                  ", where a numeric value is needed."))
    }
    if(ncol(null_thresholds) > no_categories[node]) {
      if(any(!is.na(null_thresholds[node, (no_categories[node]+1):ncol(null_thresholds)]))){
        warning(paste0("The matrix ``null_thresholds'' contains numeric values for node ",
                       node,
                       " for categories(s) (i.e., columns) exceding the maximum of ",
                       no_categories[node],
                       ". These values will be ignored."))
      }
    }
  }

  for(node in 1:no_nodes) {
    for(category in 1:no_categories[node]) {
      if(!is.finite(thresholds[node, category]))
        stop(paste("The threshold parameter for node", node, "and category",
                   category, " at t = 0 is NA or not finite."))
    }
  }


  x = sample_panel_mrf_gibbs(no_states = no_states,
                             no_nodes = no_nodes,
                             no_timepoints = no_timepoints,
                             no_categories = no_categories,
                             cross_sectional_interactions = cross_sectional_interactions,
                             cross_lagged_interactions = cross_lagged_interactions,
                             thresholds = thresholds,
                             null_interactions = null_interactions,
                             null_thresholds = null_thresholds,
                             iter)

  return(x)
}

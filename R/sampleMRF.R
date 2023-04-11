#' Sample states of the ordinal MRF
#'
#' This function samples states from the ordinal MRF using a Gibbs sampler. The
#' Gibbs sampler is initiated with random values from the response options,
#' after which it proceeds by simulating states for each node from a logistic
#' model using the other node states as predictor variables.
#'
#' @param no_states The number of states of the ordinal MRF to be generated.
#'
#' @param no_nodes The number of nodes in the ordinal MRF.
#'
#' @param no_categories Either a positive integer or a vector of positive
#' integers of length \code{no_nodes}. The number of response categories on top
#' of the base category: \code{no_categories = 1} generates binary states.
#'
#' @param interactions A symmetric \code{no_nodes} by \code{no_nodes} matrix of
#' pairwise interactions. Only its off-diagonal elements are used.
#'
#' @param thresholds A \code{no_nodes} by \code{max(no_categories)} matrix of
#' category thresholds. The elements in row \code{r} indicate the thresholds of
#' node \code{r}. If \code{no_categories} is a vector, only the first
#' \code{no_categories[r]} elements are used in row \code{r}.
#'
#' @param iter The number of iterations used by the Gibbs sampler.
#' The function provides the last state of the Gibbs sampler as output. By
#' default set to \code{1e3}.
#'
#' @return A \code{no_states} by \code{no_nodes} matrix of simulated states of
#' the ordinal MRF.
#'
#' @examples
#' # Generate responses from a network of five binary and ordinal variables.
#' no_nodes = 5
#' no_categories = sample(1:5, size = no_nodes, replace = TRUE)
#'
#' Interactions = matrix(0, nrow = no_nodes, ncol = no_nodes)
#' Interactions[2, 1] = Interactions[4, 1] = Interactions[3, 2] =
#'   Interactions[5, 2] = Interactions[5, 4] = .25
#' Interactions = Interactions + t(Interactions)
#' Thresholds = matrix(0, nrow = no_nodes, ncol = max(no_categories))
#' x = mrfSampler(no_states = 1e3,
#'                no_nodes = no_nodes,
#'                no_categories = no_categories,
#'                interactions = Interactions,
#'                thresholds = Thresholds)
#' @export
mrfSampler = function(no_states,
                      no_nodes,
                      no_categories,
                      interactions,
                      thresholds,
                      iter = 1e3) {
  #check no_states, no_nodes, iter
  if(no_states <= 0 ||
     abs(no_states - round(no_states)) > .Machine$double.eps^.5)
    stop("``no_states'' needs be a positive integer.")
  if(no_nodes <= 0 ||
     abs(no_nodes - round(no_nodes)) > .Machine$double.eps^.5)
    stop("``no_nodes'' needs be a positive integer.")
  if(iter <= 0 ||
     abs(iter - round(iter)) > .Machine$double.eps^.5)
    stop("``iter'' needs be a positive integer.")

  #check no_categories
  if(length(no_categories) == 1) {
    if(no_categories <= 0 ||
       abs(no_categories - round(no_categories)) > .Machine$double.eps^.5)
      stop("``no_categories'' needs be a (vector of) positive integer(s).")
    no_categories = rep(no_categories, no_nodes)
  } else {
    for(node in 1:no_nodes) {
      if(no_categories[node] <= 0 ||
         abs(no_categories[node] - round(no_categories[node])) >
         .Machine$double.eps^.5)
        stop(paste("For node", node, "``no_categories'' was not a positive
                   integer."))
    }
  }

  #check interactions
  if(!isSymmetric(interactions))
    stop("The matrix ``interactions'' needs to be symmetric.")
  if(nrow(interactions) != no_nodes)
    stop("The matrix ``interactions'' needs to be have ``no_nodes'' rows and
         columns.")

  #check thresholds
  if(nrow(thresholds) != no_nodes)
    stop("The matrix ``thresholds'' needs to be have ``no_nodes'' rows.")

  for(node in 1:no_nodes) {
    for(category in 1:no_categories[node]) {
      if(!is.finite(thresholds[node, category]))
        stop(paste("The threshold parameter for node", node, "and category",
                   category, "is NA or not finite."))
    }
  }

  x <- sample_omrf_gibbs(no_states = no_states,
                         no_nodes = no_nodes,
                         no_categories = no_categories,
                         interactions = interactions,
                         thresholds = thresholds,
                         iter = iter)

  return(x)
}

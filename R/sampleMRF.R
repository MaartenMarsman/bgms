#' Sample states of the ordinal MRF
#'
#' This function samples states from the ordinal MRF using a Gibbs sampler. The
#' Gibbs sampler is initiated with random values from the response options,
#' after which it proceeds by simulating states for each node from a logistic
#' model using the other node states as predictor variables.
#'
#' There are two modeling options for the category thresholds. The default
#' option assumes that the category thresholds are free, except that the first
#' threshold is set to zero for identification. The user then only needs to
#' specify the thresholds for the remaining response categories. This option is
#' useful for any type of ordinal variable and gives the user the most freedom
#' in specifying their model.
#'
#' The Blume-Capel option is specifically designed for ordinal variables that
#' have a special type of reference_category category, such as the neutral category in a
#' Likert scale. The Blume-Capel model specifies the following quadratic model
#' for the threshold parameters:
#' \deqn{\mu_{\text{c}} = \alpha \times \text{c} + \beta \times (\text{c} - \text{r})^2,}{{\mu_{\text{c}} = \alpha \times \text{c} + \beta \times (\text{c} - \text{r})^2,}}
#' where \eqn{\mu_{\text{c}}}{\mu_{\text{c}}} is the threshold for category c
#' (which now includes zero), \eqn{\alpha}{\alpha} offers a linear trend
#' across categories (increasing threshold values if
#' \eqn{\alpha > 0}{\alpha > 0} and decreasing threshold values if
#' \eqn{\alpha <0}{\alpha <0}), if \eqn{\beta < 0}{\beta < 0}, it offers an
#' increasing penalty for responding in a category further away from the
#' reference_category category r, while \eqn{\beta > 0}{\beta > 0} suggests a preference_category
#' for responding in the reference_category category.
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
#' category thresholds. The elements in row \code{i} indicate the thresholds of
#' node \code{i}. If \code{no_categories} is a vector, only the first
#' \code{no_categories[i]} elements are used in row \code{i}. If the Blume-Capel
#' model is used for the category thresholds for node \code{i}, then row
#' \code{i} requires two values (details below); the first is
#' \eqn{\alpha}{\alpha}, the linear contribution of the Blume-Capel model and
#' the second is \eqn{\beta}{\beta}, the quadratic contribution.
#'
#' @param blume_capel Either a single logical value or a logical vector of
#' length \code{no_nodes} that indicates for node \code{i} whether it is modeled
#' with Blume-Capel (\code{blume_capel[i] = TRUE}) or if it is left free
#' (\code{blume_capel[i] = FALSE}). The Blume-Capel model is only available for
#' nodes with more than two categories (\code{no_categories > 1}). Defaults to
#' \code{blume_capel = rep(FALSE, no_nodes)}.
#'
#' @param reference_category An integer vector of length \code{no_nodes} specifying the
#' reference_category category that is used for the Blume-Capel model (details below).
#' Can be any integer value between \code{0} and \code{no_categories} (or
#' \code{no_categories[i]}).
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
                      blume_capel = rep(FALSE, no_nodes),
                      reference_category,
                      iter = 1e3) {
  # Check no_states, no_nodes, iter --------------------------------------------
  if(no_states <= 0 ||
     abs(no_states - round(no_states)) > .Machine$double.eps^.5)
    stop("``no_states'' needs be a positive integer.")
  if(no_nodes <= 0 ||
     abs(no_nodes - round(no_nodes)) > .Machine$double.eps^.5)
    stop("``no_nodes'' needs be a positive integer.")
  if(iter <= 0 ||
     abs(iter - round(iter)) > .Machine$double.eps^.5)
    stop("``iter'' needs be a positive integer.")

  # Check no_categories --------------------------------------------------------
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

  # Check Blume-Capel ----------------------------------------------------------
  if(length(blume_capel) == 1) {
    if(blume_capel == TRUE && any(no_categories < 2)) {
      stop(paste0("The Blume-Capel model only works for ordinal variables with more than two response options.\\
                  But nodes ", which(no_categories < 2), " are binary variables."))
    }

    blume_capel = rep(blume_capel, no_nodes)

  } else {
    if(length(blume_capel) != no_nodes) {
      stop("The argument ``blume_capel'' should be either a single logical value or a vector of length ``no_nodes''.")
    } else {

      if(any (blume_capel == TRUE & no_categories < 2)) {
        stop(paste0("The Blume-Capel model only works for ordinal variables with more than two response options.\\
                  But nodes ", which(blume_capel == TRUE & no_categories < 2), " are binary variables."))
      }
    }
  }

  # Check interactions ---------------------------------------------------------
  if(!inherits(interactions, what = "matrix"))
    interactions = as.matrix(interactions)
  if(!isSymmetric(interactions))
    stop("The matrix ``interactions'' needs to be symmetric.")
  if(nrow(interactions) != no_nodes)
    stop("The matrix ``interactions'' needs to have ``no_nodes'' rows and columns.")

  # Check thresholds, and alpha and beta for Blume-Capel -----------------------
  if(!inherits(thresholds, what = "matrix")) {
    if(max(no_categories) == 1) {
      if(length(thresholds) == no_nodes) {
        thresholds = matrix(thresholds, ncol = 1)
      } else {
        stop(paste0("The matrix ``thresholds'' has ",
                    length(thresholds),
                    " elements, but requires",
                    no_nodes,
                    "."))
      }
    } else {
      stop("``Thresholds'' needs to be a matrix.")
    }
  }

  if(nrow(thresholds) != no_nodes)
    stop("The matrix ``thresholds'' needs to be have ``no_nodes'' rows.")

  for(node in 1:no_nodes) {
    if(blume_capel[node] == FALSE) {
      if(any(is.na(thresholds[node, 1:no_categories[node]]))) {
        tmp = which(is.na(thresholds[node, 1:no_categories[node]]))

        string = paste(tmp, sep = ",")

        stop(paste0("The matrix ``thresholds'' contains NA(s) for node ",
                    node,
                    " in categorie(s) ",
                    paste(which(is.na(thresholds[node, 1:no_categories[node]])), collapse = ", "),
                    ", where a numeric value is needed."))
      }
      if(ncol(thresholds) > no_categories[node]) {
        if(any(!is.na(thresholds[node, (no_categories[node]+1):ncol(thresholds)]))){
          warning(paste0("The matrix ``thresholds'' contains numeric values for node ",
                         node,
                         " for categories(s) (i.e., columns) exceding the maximum of ",
                         no_categories[node],
                         ". These values will be ignored."))
        }
      }
    } else {
      if(any(is.na(thresholds[node, 1:2]))) {
        stop(paste0("The Blume-Capel model is chosen for the category thresholds of node ",
                    node,
                    ". This model has two parameters that need to be placed in columns 1 and 2, row ",
                    node,
                    ", of the ``thresholds'' input matrix. Currently, there are NA(s) in these entries, where a numeric value is needed."))
      }
      if(ncol(thresholds) > 2) {
        if(any(!is.na(thresholds[node, 3:ncol(thresholds)]))){
          warning(paste0("The Blume-Capel model is chosen for the category thresholds of node ",
                         node,
                         ". This model has two parameters that need to be placed in columns 1 and 2, row ",
                         node,
                         ", of the ``thresholds'' input matrix. However, there are numeric values in higher categories. These values will be ignored."))
        }
      }
    }

  }

  for(node in 1:no_nodes) {
    if(blume_capel[node] == FALSE) {
      for(category in 1:no_categories[node]) {
        if(!is.finite(thresholds[node, category]))
          stop(paste("The threshold parameter for node", node, "and category",
                     category, "is NA or not finite."))
      }
    } else {
      if(!is.finite(thresholds[node, 1]))
        stop(paste("The alpha parameter for the Blume-Capel model for node", node, "is NA or not finite."))
      if(!is.finite(thresholds[node, 2]))
        stop(paste("The beta parameter for the Blume-Capel model for node", node, "is NA or not finite."))
    }
  }

  # Check on reference_category category on Blume-Capel specification -------------------
  if(any(blume_capel == TRUE)) {
    if(length(reference_category) == 1) {
      reference_category = rep(reference_category, no_nodes)
    }
    if(any(reference_category < 0) || any(abs(reference_category - round(reference_category)) > .Machine$double.eps^.5)) {
      stop(paste0("For nodes ",
                  which(reference_category < 0),
                  " ``reference_category'' was not a non-negative integer."))
    }
    if(any(reference_category - no_categories > 0)) {
      stop(paste0("For nodes ",
                  which(reference_category - no_categories > 0),
                  " the ``reference_category'' category was larger than the maximum category value."))
    }
  }

  if(!any(blume_capel == TRUE)) {
    x <- sample_omrf_gibbs(no_states = no_states,
                           no_nodes = no_nodes,
                           no_categories = no_categories,
                           interactions = interactions,
                           thresholds = thresholds,
                           iter = iter)
  } else {
    x <- sample_bcomrf_gibbs(no_states = no_states,
                             no_nodes = no_nodes,
                             no_categories = no_categories,
                             interactions = interactions,
                             thresholds = thresholds,
                             blume_capel = blume_capel,
                             reference = reference_category,
                             iter = iter)
  }

  return(x)
}

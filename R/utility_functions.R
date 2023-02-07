#' Data checks and reformatting for Bayesian graphical models
#'
#' The function \code{reformat_data} does some data checks to see if variables 
#' do not have exclusively unique observations (e.g., continuous responses) or 
#' all responses in one category. If one of these situations occur the program
#' halts immediately and the issue is reported back. 
#' 
#' The model and associated software takes in category scores. The levels or 
#' responses are therefore coded to be consistent with this. The lowest level or
#' response is therefore coded \code{0}, the next option is coded \code{1}, etc.
#' Note that response categories that are missing in the data are automatically
#' collapsed. If the second category is not observed, the third category will be
#' coded as \code{2}, etc.   
#'
#' @param x An \code{n} by \code{p} matrix containing the binary and ordinal 
#' variables for \code{n} independent observations on \code{p} variables in the 
#' network or graph. If not done already, \code{bsl} recodes the variables as 
#' non-negative integers (i.e., \code{0, 1, ..., m}).
#' 
#' @return A list containing the recoded \code{n} by \code{p} matrix \code{x}, 
#' the vector \code{no_categories} of length \code{p} that specifies for each
#' variable how many categories it has.
reformat_data = function(x) {
  no_nodes = ncol(x)
  no_categories = vector(length = no_nodes)
  for(node in 1:no_nodes) {
    unq_vls = sort(unique(x[,  node]))
    mx_vl = max(unq_vls)
    
    # Check if observed responses are not all unique ---------------------------
    if(mx_vl == nrow(x))
      stop(paste0("Only unique responses observed for variable ", 
                  node, 
                  ". We expect >= 1 observations per category."))
    if(length(unq_vls) != mx_vl + 1 || any(unq_vls != 0:mx_vl)) {
      y = x[, node]
      cntr = 0
      for(value in unq_vls) {
        x[y == value, node] = cntr
        cntr = cntr + 1
      }
    }
    no_categories[node] = max(x[,node])
    
    # Check to see if not all responses are in one category --------------------
    if(no_categories[node] == 0)
      stop(paste0("Only one value [", 
                  unq_vls,  
                  "] was observed for variable ", 
                  node, 
                  "."))
  }
  return(list(x = x, no_categories = no_categories))
}

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
#' @param no_iterations The number of iterations used by the Gibbs sampler.
#' The function provides the last state of the Gibbs sampler as output. By
#' default set to \code{1e3}.
#'
#' @return A \code{no_states} by \code{no_nodes} matrix of simulated states of
#' the ordinal MRF.
mrfSampler = function(no_states,
                      no_nodes,
                      no_categories,
                      interactions,
                      thresholds,
                      no_iterations = 1e3) {
  #check no_states, no_nodes, no_iterations
  if(no_states <= 0 || 
     abs(no_states - round(no_states)) > .Machine$double.eps^.5)
    stop("``no_states'' needs be a positive integer.")
  if(no_nodes <= 0 || 
     abs(no_nodes - round(no_nodes)) > .Machine$double.eps^.5)
    stop("``no_nodes'' needs be a positive integer.")
  if(no_iterations <= 0 || 
     abs(no_iterations - round(no_iterations)) > .Machine$double.eps^.5)
    stop("``no_iterations'' needs be a positive integer.")
  
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
                         no_iterations = no_iterations)

  return(x)
}

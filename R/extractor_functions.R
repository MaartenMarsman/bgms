#' Extractor Functions.
#'
#' @rdname extractor_functions
#' @param bgms_object A fit object created by the bgms package or specifically
#' by the bgm function.
#' @details Extract results from bgm objects in a safe way. Mainly intended for
#' developers of packages that build on top of the bgms package.
#' @keywords internal
#' @export
extract_arguments <- function(bgms_object) {
  if(!inherits(bgms_object, what = "bgms"))
    stop(paste0("Expected an object with class bgms and not one with class ",
      class(bgms_object)))

  if(is.null(bgms_object$arguments)) {
    stop(paste0("Extractor functions have been defined for bgms versions 0.1.3 and up but not \n",
                "for older versions. The current fit object predates version 0.1.3."))
  } else {
    return(bgms_object$arguments)
  }
}

#' @rdname extractor_functions
#' @export
extract_edge_indicators <- function(bgms_object) {
  arguments = extract_arguments(bgms_object)
  if(arguments$save) {
    edge_indicators = bgms_object$gamma
    return(edge_indicators)
  } else {
    stop(paste0("To access the sampled edge indicators the bgms package needs to be run using \n",
                "save = TRUE."))
  }
}

#' @rdname extractor_functions
#' @export
extract_posterior_inclusion_probabilities <- function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  if(arguments$save) {
    edge_means = colMeans(bgms_object$gamma)
    no_variables = arguments$no_variables

    posterior_inclusion_probabilities = matrix(0, no_variables, no_variables)
    posterior_inclusion_probabilities[lower.tri(posterior_inclusion_probabilities)] = edge_means
    posterior_inclusion_probabilities = posterior_inclusion_probabilities +
      t(posterior_inclusion_probabilities)

    data_columnnames = arguments$data_columnnames
    colnames(posterior_inclusion_probabilities) = data_columnnames
    rownames(posterior_inclusion_probabilities) = data_columnnames

  } else {

    posterior_inclusion_probabilities = bgms_object$gamma

  }
  return(posterior_inclusion_probabilities)
}


#' @rdname extractor_functions
#' @export
extract_edge_priors <- function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  if(!arguments$edge_selection) {
    stop(paste0("The bgm function did not perform edge selection, so there are no edge priors \n",
                "specified."))
  } else {
    if(arguments$edge_prior == "Bernoulli") {
      edge_prior = list(type = "Bernoulli",
                        prior_inclusion_probability = arguments$inclusion_probability)
    } else {
      edge_prior = list(type = "Beta-Bernoulli",
                        alpha = arguments$beta_bernoulli_alpha,
                        beta = arguments$beta_bernoulli_beta)
    }
  }
  return(edge_prior)
}

extract_pairwise_interactions <- function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  return(bgms_object$interactions)
}

#' @rdname extractor_functions
#' @export
extract_pairwise_thresholds <- function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  return(bgms_object$thresholds)
}
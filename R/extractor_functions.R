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
  UseMethod("extract_arguments")
}

#' @rdname extractor_functions
#' @export
extract_arguments.bgms <- function(bgms_object) {
  if(is.null(bgms_object$arguments)) {
    stop(paste0("Extractor functions have been defined for bgms versions 0.1.3 and up but not \n",
                "for older versions. The current fit object predates version 0.1.3."))
  } else {
    return(bgms_object$arguments)
  }
}

#' @rdname extractor_functions
#' @export
extract_arguments.bgmCompare <- function(bgms_object) {
  if(is.null(bgms_object$arguments)) {
    stop(paste0("Extractor functions have been defined for bgms versions 0.1.3 and up but not \n",
                "for older versions. The current fit object predates version 0.1.3."))
  } else {
    return(bgms_object$arguments)
  }
}

#' @rdname extractor_functions
#' @export
extract_indicators <- function(bgms_object) {
  UseMethod("extract_indicators")
}

#' @rdname extractor_functions
#' @export
extract_indicators.bgms <- function(bgms_object) {
  arguments = extract_arguments(bgms_object)
  if(arguments$edge_selection & arguments$save) {
    if(arguments$version < "0.1.4") {
      edge_indicators = bgms_object$gamma
    } else {
      edge_indicators = bgms_object$indicator
    }
    return(edge_indicators)
  } else {
    stop(paste0("To access the sampled edge indicators the bgms package needs to be run using \n",
                "edge_selection = TRUE and save = TRUE."))
  }
}

#' @rdname extractor_functions
#' @export
extract_indicators.bgmCompare <- function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  if(arguments$difference_selection & arguments$save) {
    pairwise_difference_indicator = bgms_object$pairwise_difference_indicator
    if(arguments$independent_thresholds == FALSE) {
      main_difference_indicator = bgms_object$main_difference_indicator
    } else {
      main_difference_indicator = NULL
    }
    return(list(main_difference_indicator = main_difference_indicator,
                pairwise_difference_indicator = pairwise_difference_indicator))
  } else {
    stop(paste0("To access the sampled difference indicators the bgmCompare function needs to be run using \n",
                "difference_selection = TRUE and save = TRUE."))
  }
}

#' @rdname extractor_functions
#' @export
extract_posterior_inclusion_probabilities <- function(bgms_object) {
  UseMethod("extract_posterior_inclusion_probabilities")
}

#' @rdname extractor_functions
#' @export
extract_posterior_inclusion_probabilities.bgms <- function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  if(!arguments$edge_selection) {
    stop(paste0("To estimate the posterior edge inclusion probabilities, please run the bgm \n",
                "function with edge_selection = TRUE."))
  }

  if(arguments$save) {
    edge_means = colMeans(bgms_object$indicator)
    no_variables = arguments$no_variables

    posterior_inclusion_probabilities = matrix(0, no_variables, no_variables)
    posterior_inclusion_probabilities[lower.tri(posterior_inclusion_probabilities)] = edge_means
    posterior_inclusion_probabilities = posterior_inclusion_probabilities +
      t(posterior_inclusion_probabilities)

    data_columnnames = arguments$data_columnnames
    colnames(posterior_inclusion_probabilities) = data_columnnames
    rownames(posterior_inclusion_probabilities) = data_columnnames

  } else {
    posterior_inclusion_probabilities = bgms_object$indicator
  }
  return(posterior_inclusion_probabilities)
}

#' @rdname extractor_functions
#' @export
extract_posterior_inclusion_probabilities.bgmCompare <- function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  if(!arguments$difference_selection) {
    stop(paste0("To estimate the posterior inclusion probabilities for the between-group \n",
                "parameter differences , please run the bgmCompare function with \n",
                "difference_selection = TRUE."))
  }

  if(arguments$save) {
    pairwise_difference_means = colMeans(bgms_object$pairwise_difference_indicator)
    no_variables = arguments$no_variables

    posterior_inclusion_probabilities = matrix(0, no_variables, no_variables)
    posterior_inclusion_probabilities[lower.tri(posterior_inclusion_probabilities)] = pairwise_difference_means
    posterior_inclusion_probabilities = posterior_inclusion_probabilities +
      t(posterior_inclusion_probabilities)
    if(!arguments$independent_thresholds) {
      main_difference_means = colMeans(bgms_object$main_difference_indicator)
      diag(posterior_inclusion_probabilities) = main_difference_means
    }

    data_columnnames = arguments$data_columnnames
    colnames(posterior_inclusion_probabilities) = data_columnnames
    rownames(posterior_inclusion_probabilities) = data_columnnames

  } else {
    posterior_inclusion_probabilities = bgms_object$indicator
  }
  return(posterior_inclusion_probabilities)
}

#' @rdname extractor_functions
#' @export
extract_indicator_priors <- function(bgms_object) {
  UseMethod("extract_indicator_priors")
}

#' @rdname extractor_functions
#' @export
extract_indicator_priors.bgms <- function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  if(!arguments$edge_selection) {
    stop(paste0("The bgm function did not perform edge selection, so there are no indicator\n",
                "priors specified."))
  } else {
    if(arguments$edge_prior == "Bernoulli") {
      indicator_prior = list(type = "Bernoulli",
                             prior_inclusion_probability = arguments$inclusion_probability)
    } else if (arguments$edge_prior == "Beta-Bernoulli") {
      indicator_prior = list(type = "Beta-Bernoulli",
                             alpha = arguments$beta_bernoulli_alpha,
                             beta = arguments$beta_bernoulli_beta)
    } else if (arguments$edge_prior == "Stochastic-Block") {
      indicator_prior = list(type = "Stochastic-Block",
                             beta_bernoulli_alpha = arguments$beta_bernoulli_alpha,
                             beta_bernoulli_beta = arguments$beta_bernoulli_beta,
                             dirichlet_alpha = arguments$dirichlet_alpha)
    }
  }
  return(indicator_prior)
}


#' @rdname extractor_functions
#' @export
extract_indicator_priors.bgmCompare <- function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  if(!arguments$difference_selection) {
    stop(paste0("The bgmCompare function did not perform selection on the between-group\n",
                "differences, so there are no indicator priors specified."))
  } else {
    if(arguments$pairwise_difference_prior == "Bernoulli") {
      difference_prior = list(pairwise_type = "Bernoulli",
                              prior_inclusion_probability = arguments$inclusion_probability_difference)
    } else {
      difference_prior = list(pairwise_type = "Beta-Bernoulli",
                              pairwise_alpha = arguments$pairwise_beta_bernoulli_alpha,
                              pairwise_beta = arguments$pairwise_beta_bernoulli_beta)
    }
    if(!arguments$independent_thresholds) {
      if(arguments$main_difference_prior == "Bernoulli") {
        difference_prior$main_type = "Bernoulli"
      } else {
        difference_prior$main_type = "Beta-Bernoulli"
        difference_prior$main_alpha = arguments$beta_bernoulli_alpha
        difference_prior$main_beta = arguments$beta_bernoulli_beta
      }
    }
  }
  return(difference_prior)
}


#' @rdname extractor_functions
#' @export
extract_pairwise_interactions <- function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  return(bgms_object$interactions)
}

#' @rdname extractor_functions
#' @export
extract_category_thresholds <- function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  return(bgms_object$thresholds)
}

#' @rdname extractor_functions
#' @export
extract_pairwise_difference.bgmCompare <- function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  return(bgms_object$pairwise_difference)
}

#' @rdname extractor_functions
#' @export
extract_main_difference.bgmCompare <- function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  return(bgms_object$main_difference)
}

#' @rdname extractor_functions
#' @export
extract_edge_indicators <- function(bgms_object) {
  warning(paste0("The ``extract_edge_indicators'' function is deprecated and will be removed in a \n",
                 "future release of bgms. Please use the ``extract_indicators'' function instead."))
  return(extract_indicators(bgms_object))
}

#' @rdname extractor_functions
#' @export
extract_pairwise_thresholds <- function(bgms_object) {
  warning(paste0("The ``extract_pairwise_thresholds'' function is deprecated and will be removed in a \n",
                 "future release of bgms. Please use the ``extract_category_thresholds'' function instead."))
  return(extract_category_thresholds(bgms_object))
}
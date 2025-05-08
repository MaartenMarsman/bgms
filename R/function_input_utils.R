#' @importFrom methods hasArg

check_positive_integer = function(value, name) {
  if (!is.numeric(value) || abs(value - round(value)) > .Machine$double.eps || value <= 0) {
    stop(sprintf("Parameter `%s` must be a positive integer. Got: %s", name, value))
  }
}

# Helper function for validating non-negative integers
check_non_negative_integer = function(value, name) {
  if (!is.numeric(value) || abs(value - round(value)) > .Machine$double.eps || value < 0) {
    stop(sprintf("Parameter `%s` must be a non-negative integer. Got: %s", name, value))
  }
}

# Helper function for validating logical inputs
check_logical = function(value, name) {
  value = as.logical(value)
  if (is.na(value)) {
    stop(sprintf("Parameter `%s` must be TRUE or FALSE. Got: %s", name, value))
  }
  return(value)
}

check_model = function(x,
                       variable_type,
                       reference_category,
                       interaction_scale = 2.5,
                       threshold_alpha = 0.5,
                       threshold_beta = 0.5,
                       edge_selection = TRUE,
                       edge_prior = c("Bernoulli", "Beta-Bernoulli", "Stochastic-Block"),
                       inclusion_probability = 0.5,
                       beta_bernoulli_alpha = 1,
                       beta_bernoulli_beta = 1,
                       dirichlet_alpha = dirichlet_alpha,
                       lambda = lambda) {

  #Check variable type input ---------------------------------------------------
  if(length(variable_type) == 1) {
    variable_input = variable_type
    variable_type = try(match.arg(arg = variable_type,
                                  choices = c("ordinal", "blume-capel")),
                        silent = TRUE)
    if(inherits(variable_type, what = "try-error"))
      stop(paste0("The bgm function supports variables of type ordinal and blume-capel, \n",
                  "but not of type ",
                  variable_input, "."))
    variable_bool = (variable_type == "ordinal")
    variable_bool = rep(variable_bool, ncol(x))
  } else {
    if(length(variable_type) != ncol(x))
      stop(paste0("The variable type vector variable_type should be either a single character\n",
                  "string or a vector of character strings of length p."))

    variable_input = unique(variable_type)
    variable_type = try(match.arg(arg = variable_type,
                                  choices = c("ordinal", "blume-capel"),
                                  several.ok = TRUE), silent = TRUE)

    if(inherits(variable_type, what = "try-error"))
      stop(paste0("The bgm function supports variables of type ordinal and blume-capel, \n",
                  "but not of type ",
                  paste0(variable_input, collapse = ", "), "."))

    num_types = sapply(variable_input, function(type) {
      tmp = try(match.arg(arg = type,
                          choices = c("ordinal", "blume-capel")),
                silent = TRUE)
      inherits(tmp, what = "try-error")
    })

    if(length(variable_type) != ncol(x))
      stop(paste0("The bgm function supports variables of type ordinal and blume-capel, \n",
                  "but not of type ",
                  paste0(variable_input[num_types], collapse = ", "), "."))

    variable_bool = (variable_type == "ordinal")
  }

  #Check Blume-Capel variable input --------------------------------------------
  if(any(!variable_bool)) {
    # Ordinal (variable_bool == TRUE) or Blume-Capel (variable_bool == FALSE)

    if(!hasArg("reference_category"))
      stop("The argument reference_category is required for Blume-Capel variables.")

    if(length(reference_category) != ncol(x) && length(reference_category) != 1)
      stop(paste0("The argument reference_category for the Blume-Capel model needs to be a \n",
                  "single integer or a vector of integers of length p."))

    if(length(reference_category) == 1) {
      #Check if the input is integer -------------------------------------------
      integer_check = try(as.integer(reference_category), silent = TRUE)
      if(is.na(integer_check))
        stop(paste0("The reference_category argument for the Blume-Capel model contains either \n",
                    "a missing value or a value that could not be forced into an integer value."))
      integer_check = reference_category - round(reference_category)
      if(integer_check > .Machine$double.eps)
        stop("Reference category needs to an integer value or a vector of integers of length p.")
      reference_category = rep.int(reference_category, times = ncol(x))
    }

    #Check if the input is integer -------------------------------------------
    blume_capel_variables = which(!variable_bool)
    # Ordinal (variable_bool == TRUE) or Blume-Capel (variable_bool == FALSE)

    integer_check = try(as.integer(reference_category[blume_capel_variables]),
                        silent = TRUE)
    if(anyNA(integer_check))
      stop(paste0("The reference_category argument for the Blume-Capel model contains either \n",
                  "missing values or values that could not be forced into an integer value."))

    integer_check = reference_category[blume_capel_variables] -
      round(reference_category[blume_capel_variables])

    if(any(integer_check > .Machine$double.eps)) {
      non_integers = blume_capel_variables[integer_check > .Machine$double.eps]
      if(length(non_integers) > 1) {
        stop(paste0("The entries in reference_category for variables ",
                    paste0(non_integers, collapse = ", "), " need to be integer."))
      } else {
        stop(paste0("The entry in reference_category for variable ",
                    non_integers, " needs to be an integer."))
      }
    }

    variable_lower = apply(x, 2, min, na.rm = TRUE)
    variable_upper = apply(x, 2, max, na.rm = TRUE)

    if(any(reference_category < variable_lower) | any(reference_category > variable_upper)) {
      out_of_range = which(reference_category < variable_lower | reference_category > variable_upper)
      stop(paste0("The Blume-Capel model assumes that the reference category is within the range \n",
                  "of the observed category scores. This was not the case for variable(s) \n",
                  paste0(out_of_range, collapse =", "),
                  "."))
    }

  } else {
    reference_category = rep.int(0, times = ncol(x))
  }

  #Check prior set-up for the interaction parameters ---------------------------
  if(interaction_scale <= 0 || is.na(interaction_scale) || is.infinite(interaction_scale))
    stop("The scale of the Cauchy prior needs to be positive.")

  #Check prior set-up for the threshold parameters -----------------------------
  if(threshold_alpha <= 0 | !is.finite(threshold_alpha))
    stop("Parameter threshold_alpha needs to be positive.")
  if(threshold_beta <= 0 | !is.finite(threshold_beta))
    stop("Parameter threshold_beta needs to be positive.")

  #Check set-up for the Bayesian edge selection model --------------------------
  edge_selection = as.logical(edge_selection)
  if(is.na(edge_selection))
    stop("The parameter edge_selection needs to be TRUE or FALSE.")
  if(edge_selection == TRUE) {
    #Check prior set-up for the edge indicators --------------------------------
    edge_prior = match.arg(edge_prior)
    if(edge_prior == "Bernoulli") {
      if(length(inclusion_probability) == 1) {
        theta = inclusion_probability[1]
        if(is.na(theta) || is.null(theta))
          stop("There is no value specified for the inclusion probability.")
        if(theta <= 0)
          stop("The inclusion probability needs to be positive.")
        if(theta > 1)
          stop("The inclusion probability cannot exceed the value one.")
        if(theta == 1)
          stop("The inclusion probability cannot equal one.")

        theta = matrix(theta, nrow = ncol(x), ncol = ncol(x))
      } else {
        if(!inherits(inclusion_probability, what = "matrix") &&
           !inherits(inclusion_probability, what = "data.frame"))
          stop("The input for the inclusion probability argument needs to be a single number, matrix, or dataframe.")

        if(inherits(inclusion_probability, what = "data.frame")) {
          theta = data.matrix(inclusion_probability)
        } else {
          theta = inclusion_probability
        }
        if(!isSymmetric(theta))
          stop("The inclusion probability matrix needs to be symmetric.")
        if(ncol(theta) != ncol(x))
          stop("The inclusion probability matrix needs to have as many rows (columns) as there are variables in the data.")

        if(anyNA(theta[lower.tri(theta)]) ||
           any(is.null(theta[lower.tri(theta)])))
          stop("One or more elements of the elements in inclusion probability matrix are not specified.")
        if(any(theta[lower.tri(theta)] <= 0))
          stop(paste0("The inclusion probability matrix contains negative or zero values;\n",
                      "inclusion probabilities need to be positive."))
        if(any(theta[lower.tri(theta)] >= 1))
          stop(paste0("The inclusion probability matrix contains values greater than or equal to one;\n",
                      "inclusion probabilities cannot exceed or equal the value one."))
      }
    }
    if(edge_prior == "Beta-Bernoulli") {
      theta = matrix(0.5, nrow = ncol(x), ncol = ncol(x))
      if(beta_bernoulli_alpha <= 0 || beta_bernoulli_beta <= 0)
        stop("The scale parameters of the beta distribution need to be positive.")
      if(!is.finite(beta_bernoulli_alpha) || !is.finite(beta_bernoulli_beta))
        stop("The scale parameters of the beta distribution need to be finite.")
      if(is.na(beta_bernoulli_alpha) || is.na(beta_bernoulli_beta) ||
         is.null(beta_bernoulli_alpha) || is.null(beta_bernoulli_beta))
        stop("Values for both scale parameters of the beta distribution need to be specified.")
    }
    if(edge_prior == "Stochastic-Block") {
      theta = matrix(0.5, nrow = ncol(x), ncol = ncol(x))
      if(beta_bernoulli_alpha <= 0 || beta_bernoulli_beta <= 0 || dirichlet_alpha <= 0 || lambda <= 0)
        stop("The scale parameters of the beta and Dirichlet distribution need to be positive.")
      if(!is.finite(beta_bernoulli_alpha) || !is.finite(beta_bernoulli_beta) || !is.finite(dirichlet_alpha) || !is.finite(lambda))
        stop("The scale parameters of the beta distribution, the concentration parameter of the Dirichlet distribution, and the rate parameter of the Poisson distribution need to be finite.")
      if(is.na(beta_bernoulli_alpha) || is.na(beta_bernoulli_beta) ||
         is.null(beta_bernoulli_alpha) || is.null(beta_bernoulli_beta) ||
         is.null(dirichlet_alpha) || is.null(dirichlet_alpha) || is.null(lambda) || is.null(lambda))
        stop("Values for both scale parameters of the beta distribution, the concentration parameter of the Dirichlet distribution, and the rate parameter of the Poisson distribution need to be specified.")
    }
  } else {
    theta = matrix(0.5, nrow = 1, ncol = 1)
    edge_prior = "Not Applicable"
  }

  return(list(variable_bool = variable_bool,
              reference_category = reference_category,
              edge_selection = edge_selection,
              edge_prior = edge_prior,
              inclusion_probability = theta))
}

check_compare_model = function(x,
                               y,
                               g,
                               difference_selection,
                               variable_type,
                               reference_category,
                               pairwise_difference_scale = 2.5,
                               main_difference_scale = 2.5,
                               pairwise_difference_prior = c("Bernoulli", "Beta-Bernoulli"),
                               main_difference_prior = c("Bernoulli", "Beta-Bernoulli"),
                               pairwise_difference_probability = 0.5,
                               main_difference_probability = 0.5,
                               main_beta_bernoulli_alpha = 1,
                               main_beta_bernoulli_beta = 1,
                               pairwise_beta_bernoulli_alpha = 1,
                               pairwise_beta_bernoulli_beta = 1,
                               interaction_scale = 2.5,
                               threshold_alpha = 0.5,
                               threshold_beta = 0.5,
                               main_difference_model = c("Free", "Collapse", "Constrain")) {

  if(!is.null(g)) {
    unique_g = unique(g)
    if(length(unique_g) == 0)
      stop(paste0("The bgmCompare function expects at least two groups, but the input g contains\n",
                  "no group value."))
    if(length(unique_g) == 1)
      stop(paste0("The bgmCompare function expects at least two groups, but the input g contains\n",
                  "only one group value."))
    if(length(unique_g) == length(g))
      stop(paste0("The bgmCompare function expects at least two groups, but the input g contains\n",
                  "only unique group values."))

    group = g
    for(u in unique_g) {
      group[g == u] = which(unique_g == u)
    }
    tab = tabulate(group)

    if(any(tab < 2))
      stop("One or more groups only had one member in the input g.")
  } else {
    group = c(rep.int(1, times = nrow(x)), rep.int(2, times = nrow(y)))
    x = rbind(x, y)
  }

  #Check variable type input ---------------------------------------------------
  if(length(variable_type) == 1) {
    variable_input = variable_type
    variable_type = try(match.arg(arg = variable_type,
                                  choices = c("ordinal", "blume-capel")),
                        silent = TRUE)
    if(inherits(variable_type, what = "try-error"))
      stop(paste0("The bgm function supports variables of type ordinal and blume-capel, \n",
                  "but not of type ",
                  variable_input, "."))
    variable_bool = (variable_type == "ordinal")
    variable_bool = rep(variable_bool, ncol(x))
  } else {
    if(length(variable_type) != ncol(x))
      stop(paste0("The variable type vector variable_type should be either a single character\n",
                  "string or a vector of character strings of length p."))

    variable_input = unique(variable_type)
    variable_type = try(match.arg(arg = variable_type,
                                  choices = c("ordinal", "blume-capel"),
                                  several.ok = TRUE), silent = TRUE)

    if(inherits(variable_type, what = "try-error"))
      stop(paste0("The bgm function supports variables of type ordinal and blume-capel, \n",
                  "but not of type ",
                  paste0(variable_input, collapse = ", "), "."))

    num_types = sapply(variable_input, function(type) {
      tmp = try(match.arg(arg = type,
                          choices = c("ordinal", "blume-capel")),
                silent = TRUE)
      inherits(tmp, what = "try-error")
    })

    if(length(variable_type) != ncol(x))
      stop(paste0("The bgm function supports variables of type ordinal and blume-capel, \n",
                  "but not of type ",
                  paste0(variable_input[num_types], collapse = ", "), "."))

    variable_bool = (variable_type == "ordinal")
  }

  #Check Blume-Capel variable input --------------------------------------------
  if(any(!variable_bool)) {
    # Ordinal (variable_bool == TRUE) or Blume-Capel (variable_bool == FALSE)

    if(!hasArg("reference_category"))
      stop("The argument reference_category is required for Blume-Capel variables.")

    if(length(reference_category) != ncol(x) && length(reference_category) != 1)
      stop(paste0("The argument reference_category for the Blume-Capel model needs to be a \n",
                  "single integer or a vector of integers of length p."))

    if(length(reference_category) == 1) {
      #Check if the input is integer -------------------------------------------
      integer_check = try(as.integer(reference_category), silent = TRUE)
      if(is.na(integer_check))
        stop(paste0("The reference_category argument for the Blume-Capel model contains either \n",
                    "a missing value or a value that could not be forced into an integer value."))
      integer_check = reference_category - round(reference_category)
      if(integer_check > .Machine$double.eps)
        stop("Reference category needs to an integer value or a vector of integers of length p.")
      reference_category = rep.int(reference_category, times = ncol(x))
    }

    #Check if the input is integer -------------------------------------------
    blume_capel_variables = which(!variable_bool)
    # Ordinal (variable_bool == TRUE) or Blume-Capel (variable_bool == FALSE)

    integer_check = try(as.integer(reference_category[blume_capel_variables]),
                        silent = TRUE)
    if(anyNA(integer_check))
      stop(paste0("The reference_category argument for the Blume-Capel model contains either \n",
                  "missing values or values that could not be forced into an integer value."))

    integer_check = reference_category[blume_capel_variables] -
      round(reference_category[blume_capel_variables])

    if(any(integer_check > .Machine$double.eps)) {
      non_integers = blume_capel_variables[integer_check > .Machine$double.eps]
      if(length(non_integers) > 1) {
        stop(paste0("The entries in reference_category for variables ",
                    paste0(non_integers, collapse = ", "), " need to be integer."))
      } else {
        stop(paste0("The entry in reference_category for variable ",
                    non_integers, " needs to be an integer."))
      }
    }

    variable_lower = apply(x, 2, min, na.rm = TRUE)
    variable_upper = apply(x, 2, max, na.rm = TRUE)

    if(any(reference_category < variable_lower) | any(reference_category > variable_upper)) {
      out_of_range = which(reference_category < variable_lower | reference_category > variable_upper)
      stop(paste0("The Blume-Capel model assumes that the reference category is within the range \n",
                  "of the observed category scores. This was not the case for variable(s) \n",
                  paste0(out_of_range, collapse =", "),
                  "."))
    }
  } else {
    reference_category = rep.int(0, times = ncol(x))
  }

  #Check prior set-up for the interaction parameters ---------------------------
  if(interaction_scale <= 0 || is.na(interaction_scale) || is.infinite(interaction_scale))
    stop("The scale of the Cauchy prior for the interactions needs to be positive.")

  #Check prior set-up for the interaction differences --------------------------
  if(pairwise_difference_scale <= 0 || is.na(pairwise_difference_scale) || is.infinite(pairwise_difference_scale))
    stop("The scale of the Cauchy prior for the difference in pairwise interactions needs to be positive.")

  #Check prior set-up for the interaction differences --------------------------
  if(main_difference_scale <= 0 || is.na(main_difference_scale) || is.infinite(main_difference_scale))
    stop("The scale of the Cauchy prior for the difference in category thresholds needs to be positive.")

  #Check prior set-up for the threshold parameters -----------------------------
  if(threshold_alpha <= 0 | !is.finite(threshold_alpha))
    stop("Parameter threshold_alpha needs to be positive.")
  if(threshold_beta <= 0 | !is.finite(threshold_beta))
    stop("Parameter threshold_beta needs to be positive.")

  #Check set-up for the Bayesian difference selection model --------------------
  difference_selection = as.logical(difference_selection)
  if(is.na(difference_selection))
    stop("The parameter difference_selection needs to be TRUE or FALSE.")
  if(difference_selection == TRUE) {
    inclusion_probability_difference = matrix(0,
                                              nrow = ncol(x),
                                              ncol = ncol(x))

    main_difference_prior = match.arg(main_difference_prior)
    if(main_difference_prior == "Bernoulli") {
      if(length(main_difference_probability) == 1) {
        main_difference_inclusion_probability = main_difference_probability[1]
        if(is.na(main_difference_inclusion_probability) || is.null(main_difference_inclusion_probability))
          stop("There is no value specified for the inclusion probability for category threshold differences.")
        if(main_difference_inclusion_probability <= 0)
          stop("The inclusion probability for category threshold differences needs to be positive.")
        if(main_difference_inclusion_probability > 1)
          stop("The inclusion probability for category threshold differences cannot exceed the value one.")
        if(main_difference_inclusion_probability == 1)
          stop("The inclusion probability for category threshold differences cannot equal one.")

        diag(inclusion_probability_difference) = main_difference_inclusion_probability

      } else {
        if(!inherits(main_difference_probability, what = "matrix") &&
           !inherits(main_difference_probability, what = "data.frame"))
          stop(paste0("The input for the inclusion probability argument for category threshold \n",
                      "differences needs to be a single number, matrix, or dataframe."))

        if(length(as.vector(main_difference_probability)) != ncol(x)) {
          stop(paste0("The inclusion probability argument for the category threshold \n",
                      "parameters needs to have either one or as many elements as\n",
                      "there are variables in the data."))
        } else {
          diag(inclusion_probability_difference) = as.vector(main_difference_probability)
        }

        if(anyNA(diag(inclusion_probability_difference)) ||
           any(is.null(diag(inclusion_probability_difference))))
          stop(paste0("One or more of the inclusion probabilities for the category threshold\n",
                      "parameters are not specified."))
        if(any(diag(inclusion_probability_difference) <= 0))
          stop(paste0("One or more of the inclusion probabilities for the category threshold\n",
                      "parameters are non-positive."))
        if(any(diag(inclusion_probability_difference) >= 1))
          stop(paste0("One or more of the inclusion probabilities for the category threshold\n",
                      "parameters are larger than one."))
      }

    } else {
      diag(inclusion_probability_difference) = 0.5
      if(main_beta_bernoulli_alpha <= 0 || main_beta_bernoulli_beta <= 0)
        stop(paste0("The scale parameters of the beta distribution for the differences in category\n",
                    "thresholds need to be positive."))
      if(!is.finite(main_beta_bernoulli_alpha) || !is.finite(main_beta_bernoulli_beta))
        stop(paste0("The scale parameters of the beta distribution for the differences in category\n",
                    "thresholds need to be finite."))
      if(is.na(main_beta_bernoulli_alpha) || is.na(main_beta_bernoulli_beta) ||
         is.null(main_beta_bernoulli_alpha) || is.null(main_beta_bernoulli_beta))
        stop(paste0("The scale parameters of the beta distribution for the differences in category\n",
                    "thresholds need to be specified."))
    }

    pairwise_difference_prior = match.arg(pairwise_difference_prior)
    if(pairwise_difference_prior == "Bernoulli") {
      if(length(pairwise_difference_probability) == 1) {
        pairwise_difference_inclusion_probability = pairwise_difference_probability[1]
        if(is.na(pairwise_difference_inclusion_probability) || is.null(pairwise_difference_inclusion_probability))
          stop(paste0("There is no value specified for the inclusion probability for the differences \n",
                      "in the pairwise interactions."))
        if(pairwise_difference_inclusion_probability <= 0)
          stop(paste0("The inclusion probability for the differences in the pairwise interactions \n",
                      "needs to be positive."))
        if(pairwise_difference_inclusion_probability > 1)
          stop(paste0("The inclusion probability for the differences in the pairwise interactions\n","
                      cannot exceed the value one."))
        if(pairwise_difference_inclusion_probability == 1)
          stop(paste0("The inclusion probability for the differences in the pairwise interactions\n", "
                      cannot equal one."))

        inclusion_probability_difference[lower.tri(inclusion_probability_difference)] =
          pairwise_difference_inclusion_probability
        inclusion_probability_difference[upper.tri(inclusion_probability_difference)] =
          pairwise_difference_inclusion_probability

      } else {
        if(!inherits(pairwise_difference_probability, what = "matrix") &&
           !inherits(pairwise_difference_probability, what = "data.frame"))
          stop(paste0("The input for the inclusion probability argument for differences in the \n",
                      "pairwise interactions needs to be a single number, matrix, or dataframe."))

        if(inherits(pairwise_difference_probability, what = "data.frame")) {
          tmp = data.matrix(pairwise_difference_probability)
          diag(tmp) = diag(inclusion_probability_difference)
          inclusion_probability_difference = tmp
        } else {
          tmp = pairwise_difference_probability
          diag(tmp) = diag(inclusion_probability_difference)
          inclusion_probability_difference = tmp
        }

        if(!isSymmetric(inclusion_probability_difference))
          stop("The inclusion probability matrix needs to be symmetric.")
        if(ncol(inclusion_probability_difference) != ncol(x))
          stop(paste0("The inclusion probability matrix needs to have as many rows (columns) as there\n",
                      " are variables in the data."))

        if(anyNA(inclusion_probability_difference[lower.tri(inclusion_probability_difference)]) ||
           any(is.null(inclusion_probability_difference[lower.tri(inclusion_probability_difference)])))
          stop(paste0("One or more inclusion probabilities for differences in the pairwise \n",
                      "interactions are not specified."))
        if(any(inclusion_probability_difference[lower.tri(inclusion_probability_difference)] <= 0))
          stop(paste0("One or more inclusion probabilities for differences in the pairwise \n",
                      "interactions are negative or equal to zero."))
        if(any(inclusion_probability_difference[lower.tri(inclusion_probability_difference)] >= 1))
          stop(paste0("One or more inclusion probabilities for differences in the pairwise \n",
                      "interactions are larger than one."))
      }
    } else {
      inclusion_probability_difference[lower.tri(inclusion_probability_difference)] =
        0.5
      inclusion_probability_difference[upper.tri(inclusion_probability_difference)] =
        0.5

      if(pairwise_beta_bernoulli_alpha <= 0 || pairwise_beta_bernoulli_beta <= 0)
        stop(paste0("The scale parameters of the beta distribution for the differences in pairwise\n",
                    "interactions need to be positive."))
      if(!is.finite(pairwise_beta_bernoulli_alpha) || !is.finite(pairwise_beta_bernoulli_beta))
        stop(paste0("The scale parameters of the beta distribution for the differences in pairwise\n",
                    "interactions need to be finite."))
      if(is.na(pairwise_beta_bernoulli_alpha) || is.na(pairwise_beta_bernoulli_beta) ||
         is.null(pairwise_beta_bernoulli_alpha) || is.null(pairwise_beta_bernoulli_beta))
        stop(paste0("The scale parameters of the beta distribution for the differences in pairwise\n",
                    "differences need to be specified."))
    }
  } else {
    main_difference_prior = "Not applicable"
    pairwise_difference_prior = "Not applicable"
    inclusion_probability_difference = matrix(0.5, 1, 1)
  }

  main_difference_model = match.arg(main_difference_model)

  u_g = unique(g)
  if(length(u_g) > 2 & main_difference_model == "Constrain")
    stop(paste0("The option Constrain for the main difference model is currently only supported \n",
                "for designs with two groups. Please try again with Free or Collapse."))

  if(main_difference_model == "Constrain")
    warning(paste0("The option Constrain for the main difference model is deprecated and will be\n",
                   "removed in a future version. Please use the Collapse option."))

  return(list(x = x,
              group = group,
              variable_bool = variable_bool,
              reference_category = reference_category,
              main_difference_prior = main_difference_prior,
              pairwise_difference_prior = pairwise_difference_prior,
              inclusion_probability_difference = inclusion_probability_difference,
              main_difference_model = main_difference_model))
}
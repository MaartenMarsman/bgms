# ==============================================================================
# Model validation functions
# ==============================================================================
#
# Pure validation functions for model-level arguments (variable types,
# baseline categories, edge priors, difference priors). Each function
# takes input and returns validated output (or errors).
# ==============================================================================

# ------------------------------------------------------------------------------
# validate_variable_types
# ------------------------------------------------------------------------------
#
# Parses and validates the variable_type argument for both bgm() and
# bgmCompare(). Returns a canonical list with:
#   $variable_type  - character vector of length num_variables
#   $variable_bool  - logical vector (TRUE = ordinal, FALSE = blume-capel)
#   $is_continuous  - scalar logical
#
# @param variable_type Character: single string or vector of variable types.
# @param num_variables Integer: number of columns in the data matrix.
# @param allow_continuous Logical: whether "continuous" is a valid type.
#   TRUE for bgm(), FALSE for bgmCompare() (until GGM-Compare is added).
# @param caller Character: name of the calling function, used in error messages.
# ------------------------------------------------------------------------------
validate_variable_types = function(variable_type,
                                   num_variables,
                                   allow_continuous = TRUE,
                                   allow_mixed = FALSE,
                                   caller = "bgm") {
  valid_choices = if(allow_continuous) {
    c("ordinal", "blume-capel", "continuous")
  } else {
    c("ordinal", "blume-capel")
  }

  supported_str = paste(valid_choices, collapse = ", ")

  is_continuous = FALSE
  is_mixed = FALSE

  if(length(variable_type) == 1) {
    # --- Single string: replicate to all variables ---
    variable_input = variable_type
    variable_type = try(
      match.arg(arg = variable_type, choices = valid_choices),
      silent = TRUE
    )
    if(inherits(variable_type, what = "try-error")) {
      stop(paste0(
        "The ", caller, " function supports variables of type ", supported_str,
        ", but not of type ", variable_input, "."
      ))
    }

    if(variable_type == "continuous") {
      is_continuous = TRUE
      variable_bool = rep(TRUE, num_variables)
      variable_type = rep("continuous", num_variables)
    } else {
      variable_bool = (variable_type == "ordinal")
      variable_bool = rep(variable_bool, num_variables)
      variable_type = rep(variable_type, num_variables)
    }
  } else {
    # --- Vector of types: validate each element ---
    if(length(variable_type) != num_variables) {
      stop(paste0(
        "The variable type vector variable_type should be either a single character\n",
        "string or a vector of character strings of length p."
      ))
    }

    has_continuous = any(variable_type == "continuous")
    has_discrete = any(variable_type %in% c("ordinal", "blume-capel"))
    is_mixed = has_continuous && has_discrete

    if(has_continuous && !has_discrete) {
      invalid_if_cont = setdiff(unique(variable_type), "continuous")
      if(length(invalid_if_cont) > 0) {
        stop(paste0(
          "When using continuous variables, all variables must be of type ",
          "'continuous' or mixed with ordinal/blume-capel variables."
        ))
      }
    }

    if(has_continuous && !allow_continuous) {
      stop(paste0(
        "The ", caller, " function supports variables of type ", supported_str,
        ", but not of type continuous."
      ))
    }

    if(is_mixed && !allow_mixed) {
      stop(paste0(
        "When using continuous variables, all variables must be of type ",
        "'continuous'. Mixtures of continuous and ordinal/blume-capel ",
        "variables are not supported by ", caller, "()."
      ))
    }

    if(is_mixed) {
      # Mixed: validate each entry individually
      variable_type_checked = try(
        match.arg(
          arg = variable_type,
          choices = valid_choices,
          several.ok = TRUE
        ),
        silent = TRUE
      )

      if(inherits(variable_type_checked, what = "try-error")) {
        invalid = setdiff(unique(variable_type), valid_choices)
        stop(paste0(
          "The ", caller, " function supports variables of type ", supported_str,
          ", but not of type ", paste0(invalid, collapse = ", "), "."
        ))
      }

      variable_type = variable_type_checked
      # variable_bool: TRUE = ordinal; FALSE = blume-capel/continuous
      variable_bool = (variable_type == "ordinal")
      is_continuous = FALSE
    } else if(has_continuous) {
      if(!allow_continuous) {
        stop(paste0(
          "The ", caller, " function supports variables of type ", supported_str,
          ", but not of type continuous."
        ))
      }
      is_continuous = TRUE
      variable_bool = rep(TRUE, num_variables)
    } else {
      variable_input = unique(variable_type)
      non_continuous_choices = c("ordinal", "blume-capel")

      variable_type_checked = try(
        match.arg(
          arg = variable_type,
          choices = non_continuous_choices,
          several.ok = TRUE
        ),
        silent = TRUE
      )

      if(inherits(variable_type_checked, what = "try-error")) {
        # Identify which types are invalid
        num_types = sapply(variable_input, function(type) {
          tmp = try(
            match.arg(arg = type, choices = non_continuous_choices),
            silent = TRUE
          )
          inherits(tmp, what = "try-error")
        })

        stop(paste0(
          "The ", caller, " function supports variables of type ", supported_str,
          ", but not of type ",
          paste0(variable_input[num_types], collapse = ", "), "."
        ))
      }

      # The match.arg with several.ok may have normalized partial matches
      variable_type = variable_type_checked

      # Re-check length after match.arg (defensive)
      if(length(variable_type) != num_variables) {
        num_types = sapply(variable_input, function(type) {
          tmp = try(
            match.arg(arg = type, choices = non_continuous_choices),
            silent = TRUE
          )
          inherits(tmp, what = "try-error")
        })

        stop(paste0(
          "The ", caller, " function supports variables of type ", supported_str,
          ", but not of type ",
          paste0(variable_input[num_types], collapse = ", "), "."
        ))
      }

      variable_bool = (variable_type == "ordinal")
    }
  }

  list(
    variable_type  = variable_type,
    variable_bool  = variable_bool,
    is_continuous  = is_continuous,
    is_mixed       = is_mixed
  )
}


# ------------------------------------------------------------------------------
# validate_baseline_category
# ------------------------------------------------------------------------------
#
# Validates and normalizes the baseline_category argument for Blume-Capel
# variables. Shared by both bgm() and bgmCompare().
#
# @param baseline_category  The user-supplied baseline_category value.
# @param baseline_category_provided  Logical: whether the user actually
#   supplied the argument (i.e., `hasArg("baseline_category")` from the
#   calling scope). Needed because baseline_category has no default.
# @param x  Numeric matrix: the (validated) data.
# @param variable_bool  Logical vector: TRUE = ordinal, FALSE = blume-capel.
#
# Returns:
#   Integer vector of length ncol(x). For ordinal-only models, all zeros.
# ------------------------------------------------------------------------------
validate_baseline_category = function(baseline_category,
                                      baseline_category_provided,
                                      x,
                                      variable_bool) {
  num_variables = ncol(x)

  # If all ordinal (no Blume-Capel), return zeros
  if(!any(!variable_bool)) {
    return(rep.int(0, times = num_variables))
  }

  # --- Blume-Capel variables present ---

  if(!baseline_category_provided) {
    stop("The argument baseline_category is required for Blume-Capel variables.")
  }

  if(length(baseline_category) != num_variables && length(baseline_category) != 1) {
    stop(paste0(
      "The argument baseline_category for the Blume-Capel model needs to be a \n",
      "single integer or a vector of integers of length p."
    ))
  }

  # Scalar: validate then replicate
  if(length(baseline_category) == 1) {
    integer_check = try(as.integer(baseline_category), silent = TRUE)
    if(is.na(integer_check)) {
      stop(paste0(
        "The baseline_category argument for the Blume-Capel model contains either \n",
        "a missing value or a value that could not be forced into an integer value."
      ))
    }
    integer_check = abs(baseline_category - round(baseline_category))
    if(integer_check > .Machine$double.eps) {
      stop("Reference category needs to an integer value or a vector of integers of length p.")
    }
    baseline_category = rep.int(baseline_category, times = num_variables)
  }

  # Validate integer-ness for Blume-Capel variables
  blume_capel_variables = which(!variable_bool)

  integer_check = try(as.integer(baseline_category[blume_capel_variables]),
    silent = TRUE
  )
  if(anyNA(integer_check)) {
    stop(paste0(
      "The baseline_category argument for the Blume-Capel model contains either \n",
      "missing values or values that could not be forced into an integer value."
    ))
  }

  integer_check = abs(baseline_category[blume_capel_variables] -
    round(baseline_category[blume_capel_variables]))

  if(any(integer_check > .Machine$double.eps)) {
    non_integers = blume_capel_variables[integer_check > .Machine$double.eps]
    if(length(non_integers) > 1) {
      stop(paste0(
        "The entries in baseline_category for variables ",
        paste0(non_integers, collapse = ", "), " need to be integer."
      ))
    } else {
      stop(paste0(
        "The entry in baseline_category for variable ",
        non_integers, " needs to be an integer."
      ))
    }
  }

  # Validate within observed data range
  variable_lower = apply(x, 2, min, na.rm = TRUE)
  variable_upper = apply(x, 2, max, na.rm = TRUE)

  if(any(baseline_category < variable_lower) || any(baseline_category > variable_upper)) {
    out_of_range = which(baseline_category < variable_lower | baseline_category > variable_upper)
    stop(paste0(
      "The Blume-Capel model assumes that the reference category is within the range \n",
      "of the observed category scores. This was not the case for variable(s) \n",
      paste0(out_of_range, collapse = ", "),
      "."
    ))
  }

  baseline_category
}


# ------------------------------------------------------------------------------
# validate_edge_prior
# ------------------------------------------------------------------------------
#
# Validates and normalizes the edge selection prior setup for bgm().
# Handles Bernoulli, Beta-Bernoulli, and Stochastic-Block priors.
#
# @param edge_selection  Logical: whether to perform Bayesian edge selection.
# @param edge_prior  Character: one of "Bernoulli", "Beta-Bernoulli",
#   "Stochastic-Block".
# @param inclusion_probability  Scalar or matrix of inclusion probabilities.
# @param num_variables  Integer: number of variables (ncol of data).
# @param beta_bernoulli_alpha  Numeric: alpha shape parameter for Beta prior.
# @param beta_bernoulli_beta  Numeric: beta shape parameter for Beta prior.
# @param beta_bernoulli_alpha_between  Numeric: alpha for between-cluster
#   (Stochastic-Block only).
# @param beta_bernoulli_beta_between  Numeric: beta for between-cluster
#   (Stochastic-Block only).
# @param dirichlet_alpha  Numeric: concentration parameter for Dirichlet
#   (Stochastic-Block only).
# @param lambda  Numeric: rate parameter for Poisson
#   (Stochastic-Block only).
#
# Returns:
#   list(edge_selection, edge_prior, inclusion_probability)
#   where inclusion_probability is a num_variables x num_variables matrix.
# ------------------------------------------------------------------------------
validate_edge_prior = function(edge_selection,
                               edge_prior = c(
                                 "Bernoulli", "Beta-Bernoulli",
                                 "Stochastic-Block"
                               ),
                               inclusion_probability = 0.5,
                               num_variables,
                               beta_bernoulli_alpha = 1,
                               beta_bernoulli_beta = 1,
                               beta_bernoulli_alpha_between = 1,
                               beta_bernoulli_beta_between = 1,
                               dirichlet_alpha = 1,
                               lambda = 1) {
  edge_selection = as.logical(edge_selection)
  if(is.na(edge_selection)) {
    stop("The parameter edge_selection needs to be TRUE or FALSE.")
  }

  if(!edge_selection) {
    return(list(
      edge_selection = FALSE,
      edge_prior = "Not Applicable",
      inclusion_probability = matrix(0.5, nrow = 1, ncol = 1)
    ))
  }

  # --- edge_selection == TRUE ---
  edge_prior = match.arg(edge_prior)

  if(edge_prior == "Bernoulli") {
    theta = validate_bernoulli_inclusion(
      probability      = inclusion_probability,
      num_variables    = num_variables,
      include_diagonal = FALSE,
      context          = ""
    )
  }

  if(edge_prior == "Beta-Bernoulli") {
    theta = matrix(0.5, nrow = num_variables, ncol = num_variables)
    if(is.null(beta_bernoulli_alpha) || is.null(beta_bernoulli_beta) ||
      is.na(beta_bernoulli_alpha) || is.na(beta_bernoulli_beta)) {
      stop("Values for both scale parameters of the beta distribution need to be specified.")
    }
    if(beta_bernoulli_alpha <= 0 || beta_bernoulli_beta <= 0) {
      stop("The scale parameters of the beta distribution need to be positive.")
    }
    if(!is.finite(beta_bernoulli_alpha) || !is.finite(beta_bernoulli_beta)) {
      stop("The scale parameters of the beta distribution need to be finite.")
    }
  }

  if(edge_prior == "Stochastic-Block") {
    theta = matrix(0.5, nrow = num_variables, ncol = num_variables)

    # Check that all beta parameters are provided
    if(is.null(beta_bernoulli_alpha) || is.null(beta_bernoulli_beta) ||
      is.null(beta_bernoulli_alpha_between) || is.null(beta_bernoulli_beta_between)) {
      stop(
        "The Stochastic-Block prior requires all four beta parameters: ",
        "beta_bernoulli_alpha, beta_bernoulli_beta, ",
        "beta_bernoulli_alpha_between, and beta_bernoulli_beta_between."
      )
    }

    # Check for NAs first (NA in comparisons would crash if())
    if(is.na(beta_bernoulli_alpha) || is.na(beta_bernoulli_beta) ||
      is.na(beta_bernoulli_alpha_between) || is.na(beta_bernoulli_beta_between) ||
      is.na(dirichlet_alpha) || is.na(lambda)) {
      stop(
        "Values for all shape parameters of the beta distribution, the concentration parameter of the Dirichlet distribution, ",
        "and the rate parameter of the Poisson distribution cannot be NA."
      )
    }

    # Check that all parameters are positive
    if(beta_bernoulli_alpha <= 0 || beta_bernoulli_beta <= 0 ||
      beta_bernoulli_alpha_between <= 0 || beta_bernoulli_beta_between <= 0 ||
      dirichlet_alpha <= 0 || lambda <= 0) {
      stop("The parameters of the beta and Dirichlet distributions need to be positive.")
    }

    # Check that all parameters are finite
    if(!is.finite(beta_bernoulli_alpha) || !is.finite(beta_bernoulli_beta) ||
      !is.finite(beta_bernoulli_alpha_between) || !is.finite(beta_bernoulli_beta_between) ||
      !is.finite(dirichlet_alpha) || !is.finite(lambda)) {
      stop(
        "The shape parameters of the beta distribution, the concentration parameter of the Dirichlet distribution, ",
        "and the rate parameter of the Poisson distribution need to be finite."
      )
    }
  }

  list(
    edge_selection        = TRUE,
    edge_prior            = edge_prior,
    inclusion_probability = theta
  )
}


# ------------------------------------------------------------------------------
# validate_bernoulli_inclusion (internal helper)
# ------------------------------------------------------------------------------
#
# Shared validator for Bernoulli inclusion probabilities. Used by both
# validate_edge_prior() (edges) and validate_difference_prior() (differences).
#
# @param probability  Scalar, matrix, or data.frame of inclusion probabilities.
# @param num_variables  Integer: number of variables.
# @param include_diagonal  Logical: whether to check the diagonal in matrix
#   validation. TRUE for differences (main-effect differences), FALSE for
#   edges (pairwise only).
# @param context  Character: label for error messages, e.g. "" (edge) or
#   " for differences" (difference).
#
# Returns: num_variables x num_variables matrix of inclusion probabilities.
# ------------------------------------------------------------------------------
validate_bernoulli_inclusion = function(probability,
                                        num_variables,
                                        include_diagonal = FALSE,
                                        context = "") {
  if(length(probability) == 1) {
    theta = probability[1]
    if(is.na(theta) || is.null(theta)) {
      stop(paste0(
        "There is no value specified for the inclusion probability",
        context, "."
      ))
    }
    if(theta <= 0) {
      stop(paste0(
        "The inclusion probability", context, " needs to be positive."
      ))
    }
    if(theta >= 1) {
      stop(paste0(
        "The inclusion probability", context, " cannot equal or exceed the value one."
      ))
    }
    return(matrix(theta, nrow = num_variables, ncol = num_variables))
  }

  # --- Matrix / data.frame path ---
  if(!inherits(probability, what = "matrix") &&
    !inherits(probability, what = "data.frame")) {
    stop(paste0(
      "The input for the inclusion probability argument", context,
      " needs to be a single number, matrix, or dataframe."
    ))
  }

  if(inherits(probability, what = "data.frame")) {
    theta = data.matrix(probability)
  } else {
    theta = probability
  }
  if(!isSymmetric(theta)) {
    stop("The inclusion probability matrix needs to be symmetric.")
  }
  if(ncol(theta) != num_variables) {
    stop(paste0(
      "The inclusion probability matrix needs to have as many rows (columns) as there",
      if(nzchar(context)) "\n " else " ",
      "are variables in the data."
    ))
  }

  tri = lower.tri(theta, diag = include_diagonal)
  if(anyNA(theta[tri])) {
    if(nzchar(context)) {
      stop(paste0(
        "One or more inclusion probabilities", context, " are not specified."
      ))
    } else {
      stop("One or more elements of the elements in inclusion probability matrix are not specified.")
    }
  }
  if(any(theta[tri] <= 0)) {
    if(nzchar(context)) {
      stop(paste0(
        "One or more inclusion probabilities", context, " are negative or zero."
      ))
    } else {
      stop(paste0(
        "The inclusion probability matrix contains negative or zero values;\n",
        "inclusion probabilities need to be positive."
      ))
    }
  }
  if(any(theta[tri] >= 1)) {
    if(nzchar(context)) {
      stop(paste0(
        "One or more inclusion probabilities", context, " are one or larger."
      ))
    } else {
      stop(paste0(
        "The inclusion probability matrix contains values greater than or equal to one;\n",
        "inclusion probabilities cannot exceed or equal the value one."
      ))
    }
  }

  theta
}


# ------------------------------------------------------------------------------
# validate_difference_prior
# ------------------------------------------------------------------------------
#
# Validates and normalizes the difference selection prior setup for
# bgmCompare(). Handles Bernoulli, Beta-Bernoulli, and Stochastic-Block priors.
#
# Under Stochastic-Block, the prior governs the off-diagonal (pairwise)
# difference indicators only. The diagonal main-effect difference indicators
# (used when main_difference_selection = TRUE) are sampled separately via a
# Beta-Bernoulli conjugate update with the same (alpha, beta) hyperparameters.
#
# @param difference_selection  Logical: whether to perform Bayesian
#   difference selection.
# @param difference_prior  Character: one of "Bernoulli", "Beta-Bernoulli",
#   "Stochastic-Block".
# @param difference_probability  Scalar or matrix of inclusion probabilities
#   for differences.
# @param num_variables  Integer: number of variables (ncol of data).
# @param beta_bernoulli_alpha  Numeric: alpha for within-cluster Beta prior.
# @param beta_bernoulli_beta  Numeric: beta for within-cluster Beta prior.
# @param beta_bernoulli_alpha_between  Numeric: alpha for between-cluster
#   Beta prior (Stochastic-Block only).
# @param beta_bernoulli_beta_between  Numeric: beta for between-cluster
#   Beta prior (Stochastic-Block only).
# @param dirichlet_alpha  Numeric: Dirichlet concentration on cluster sizes
#   (Stochastic-Block only).
# @param lambda  Numeric: Poisson rate on number of clusters
#   (Stochastic-Block only).
#
# Returns:
#   list(difference_selection, difference_prior,
#        inclusion_probability_difference,
#        beta_bernoulli_alpha, beta_bernoulli_beta,
#        beta_bernoulli_alpha_between, beta_bernoulli_beta_between,
#        dirichlet_alpha, lambda)
#   where inclusion_probability_difference is a
#   num_variables x num_variables matrix.
# ------------------------------------------------------------------------------
validate_difference_prior = function(difference_selection,
                                     difference_prior = c(
                                       "Bernoulli",
                                       "Beta-Bernoulli",
                                       "Stochastic-Block"
                                     ),
                                     difference_probability = 0.5,
                                     num_variables,
                                     beta_bernoulli_alpha = 1,
                                     beta_bernoulli_beta = 1,
                                     beta_bernoulli_alpha_between = 1,
                                     beta_bernoulli_beta_between = 1,
                                     dirichlet_alpha = 1,
                                     lambda = 1) {
  difference_selection = as.logical(difference_selection)
  if(is.na(difference_selection)) {
    stop("The parameter difference_selection needs to be TRUE or FALSE.")
  }

  if(!difference_selection) {
    return(list(
      difference_selection = FALSE,
      difference_prior = "Not applicable",
      inclusion_probability_difference = matrix(0.5, 1, 1),
      beta_bernoulli_alpha = beta_bernoulli_alpha,
      beta_bernoulli_beta = beta_bernoulli_beta,
      beta_bernoulli_alpha_between = beta_bernoulli_alpha_between,
      beta_bernoulli_beta_between = beta_bernoulli_beta_between,
      dirichlet_alpha = dirichlet_alpha,
      lambda = lambda
    ))
  }

  # --- difference_selection == TRUE ---
  difference_prior = match.arg(difference_prior)

  if(difference_prior == "Bernoulli") {
    theta = validate_bernoulli_inclusion(
      probability      = difference_probability,
      num_variables    = num_variables,
      include_diagonal = TRUE,
      context          = " for differences"
    )
  } else if(difference_prior == "Beta-Bernoulli") {
    theta = matrix(0.5, nrow = num_variables, ncol = num_variables)
    if(is.null(beta_bernoulli_alpha) || is.null(beta_bernoulli_beta) ||
      is.na(beta_bernoulli_alpha) || is.na(beta_bernoulli_beta)) {
      stop("The scale parameters of the beta distribution for the differences need to be specified.")
    }
    if(beta_bernoulli_alpha <= 0 || beta_bernoulli_beta <= 0) {
      stop("The scale parameters of the beta distribution for the differences need to be positive.")
    }
    if(!is.finite(beta_bernoulli_alpha) || !is.finite(beta_bernoulli_beta)) {
      stop("The scale parameters of the beta distribution for the differences need to be finite.")
    }
  } else {
    # Stochastic-Block
    theta = matrix(0.5, nrow = num_variables, ncol = num_variables)

    if(is.null(beta_bernoulli_alpha) || is.null(beta_bernoulli_beta) ||
      is.null(beta_bernoulli_alpha_between) || is.null(beta_bernoulli_beta_between)) {
      stop(
        "The Stochastic-Block prior for differences requires all four beta parameters: ",
        "beta_bernoulli_alpha, beta_bernoulli_beta, ",
        "beta_bernoulli_alpha_between, and beta_bernoulli_beta_between."
      )
    }
    if(is.na(beta_bernoulli_alpha) || is.na(beta_bernoulli_beta) ||
      is.na(beta_bernoulli_alpha_between) || is.na(beta_bernoulli_beta_between) ||
      is.na(dirichlet_alpha) || is.na(lambda)) {
      stop(
        "Hyperparameters of the Stochastic-Block prior for differences cannot be NA."
      )
    }
    if(beta_bernoulli_alpha <= 0 || beta_bernoulli_beta <= 0 ||
      beta_bernoulli_alpha_between <= 0 || beta_bernoulli_beta_between <= 0 ||
      dirichlet_alpha <= 0 || lambda <= 0) {
      stop("Hyperparameters of the Stochastic-Block prior for differences need to be positive.")
    }
    if(!is.finite(beta_bernoulli_alpha) || !is.finite(beta_bernoulli_beta) ||
      !is.finite(beta_bernoulli_alpha_between) || !is.finite(beta_bernoulli_beta_between) ||
      !is.finite(dirichlet_alpha) || !is.finite(lambda)) {
      stop("Hyperparameters of the Stochastic-Block prior for differences need to be finite.")
    }
  }

  list(
    difference_selection = TRUE,
    difference_prior = difference_prior,
    inclusion_probability_difference = theta,
    beta_bernoulli_alpha = beta_bernoulli_alpha,
    beta_bernoulli_beta = beta_bernoulli_beta,
    beta_bernoulli_alpha_between = beta_bernoulli_alpha_between,
    beta_bernoulli_beta_between = beta_bernoulli_beta_between,
    dirichlet_alpha = dirichlet_alpha,
    lambda = lambda
  )
}

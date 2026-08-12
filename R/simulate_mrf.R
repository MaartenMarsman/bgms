# simulate_mrf() / mrfSampler() — standalone MRF simulation from parameters
#
# Split out of simulate_predict.R (cleanup S4). These generate data directly from
# a user-specified MRF (not from a fitted object); the simulate()/predict() S3
# methods and shared helpers stay in simulate_predict.R.


#' @title Simulate Observations from a Markov Random Field
#'
#' @description
#' `simulate_mrf()` generates observations from a Markov Random
#' Field using user-specified parameters. For ordinal and
#' Blume-Capel variables, observations are generated via Gibbs
#' sampling. For continuous variables (Gaussian graphical model),
#' observations are drawn directly from the multivariate normal
#' distribution implied by the precision matrix.
#'
#' @details
#' \strong{Ordinal / Blume-Capel variables:}
#' The Gibbs sampler is initiated with random values from the response options,
#' after which it proceeds by simulating states for each variable from its full
#' conditional distribution given the other variable states.
#'
#' \strong{Continuous variables (GGM):}
#' Observations are drawn from \eqn{N(\mu, \Omega^{-1})}{N(mu, Omega^{-1})}
#' where \eqn{\Omega}{Omega} is the precision matrix specified via
#' `pairwise` and \eqn{\mu}{mu} is the means vector specified via `main`.
#' No Gibbs sampling is needed; `iter` is ignored.
#'
#' There are two modeling options for the category thresholds. The default
#' option assumes that the category thresholds are free, except that the first
#' threshold is set to zero for identification. The user then only needs to
#' specify the thresholds for the remaining response categories. This option is
#' useful for any type of ordinal variable and gives the user the most freedom
#' in specifying their model.
#'
#' The Blume-Capel option is specifically designed for ordinal variables that
#' have a special type of baseline_category category, such as the neutral
#' category in a Likert scale. The Blume-Capel model specifies the following
#' quadratic model for the threshold parameters:
#' \deqn{\mu_{\text{c}} = \alpha (\text{c} - \text{r})
#'   + \beta (\text{c} - \text{r})^2}
#' where \eqn{\mu_{\text{c}}}{\mu_{\text{c}}} is the threshold for category c
#' (which now includes zero), \eqn{\alpha}{\alpha} offers a linear trend
#' across categories (increasing threshold values if
#' \eqn{\alpha > 0}{\alpha > 0} and decreasing threshold values if
#' \eqn{\alpha <0}{\alpha <0}), if \eqn{\beta < 0}{\beta < 0}, it offers an
#' increasing penalty for responding in a category further away from the
#' baseline_category category r, while \eqn{\beta > 0}{\beta > 0} suggests a
#' preference for responding in the baseline_category category.
#'
#' @param num_states The number of observations to be generated.
#'
#' @param num_variables The number of variables in the MRF.
#'
#' @param num_categories Either a positive integer or a vector
#' of positive integers of length \code{num_variables}. The
#' number of response categories on top of the base category:
#' \code{num_categories = 1} generates binary states.
#' Only used for ordinal and Blume-Capel variables; ignored when
#' \code{variable_type = "continuous"}.
#'
#' @param pairwise A symmetric \code{num_variables} by
#' \code{num_variables} matrix. For ordinal and Blume-Capel
#' variables, this contains the pairwise interaction parameters;
#' only the off-diagonal elements are used. For continuous
#' variables,
#' this is the precision matrix \eqn{\Omega}{Omega} (including diagonal) and
#' must be positive definite.
#'
#' @param main For ordinal and Blume-Capel variables: a
#' \code{num_variables} by \code{max(num_categories)} matrix of category
#' thresholds. The elements in row \code{i} indicate the thresholds of
#' variable \code{i}. If \code{num_categories} is a vector, only the first
#' \code{num_categories[i]} elements are used in row \code{i}.
#' If the Blume-Capel
#' model is used for the category thresholds for variable \code{i}, then row
#' \code{i} requires two values (details below); the first is
#' \eqn{\alpha}{\alpha}, the linear contribution of the Blume-Capel model and
#' the second is \eqn{\beta}{\beta}, the quadratic contribution.
#' For continuous variables: a numeric vector of length \code{num_variables}
#' containing the means \eqn{\mu}{mu} for each variable. Defaults to zeros
#' if not supplied (missing(main)).
#'
#' @param variable_type What kind of variables are simulated? Can be a single
#' character string specifying the variable type of all \code{p} variables at
#' once or a vector of character strings of length \code{p} specifying the type
#' for each variable separately. Currently, bgm supports \code{"ordinal"},
#' \code{"blume-capel"}, and \code{"continuous"}. Binary variables are automatically
#' treated as \code{"ordinal"}. Ordinal and Blume-Capel variables can be mixed
#' freely, but continuous variables cannot be mixed with ordinal or Blume-Capel
#' variables. When \code{variable_type = "continuous"}, the function simulates
#' from a Gaussian graphical model.
#' Defaults to \code{variable_type = "ordinal"}.
#'
#' @param baseline_category An integer vector of length
#' \code{num_variables} specifying the baseline_category
#' category that is used for the Blume-Capel model
#' (details below).
#' Can be any integer value between \code{0} and \code{num_categories} (or
#' \code{num_categories[i]}).
#'
#' @param iter The number of iterations used by the Gibbs sampler
#' (ordinal/Blume-Capel variables only). The function provides the last state
#' of the Gibbs sampler as output. Ignored for continuous variables.
#' By default set to \code{1e3}.
#'
#' @param seed Optional integer seed for reproducibility. If \code{NULL},
#' a seed is generated from R's random number generator (so \code{set.seed()}
#' can be used before calling this function).
#'
#' @return A \code{num_states} by \code{num_variables} matrix of simulated
#' observations. For ordinal/Blume-Capel variables, entries are non-negative
#' integers. For continuous variables, entries are real-valued.
#'
#' @examples
#' # Generate responses from a network of five binary and ordinal variables.
#' num_variables = 5
#' num_categories = sample(1:5, size = num_variables, replace = TRUE)
#'
#' Pairwise = matrix(0, nrow = num_variables, ncol = num_variables)
#' Pairwise[2, 1] = Pairwise[4, 1] = Pairwise[3, 2] =
#'   Pairwise[5, 2] = Pairwise[5, 4] = .25
#' Pairwise = Pairwise + t(Pairwise)
#' Main = matrix(0, nrow = num_variables, ncol = max(num_categories))
#'
#' x = simulate_mrf(
#'   num_states = 1e3,
#'   num_variables = num_variables,
#'   num_categories = num_categories,
#'   pairwise = Pairwise,
#'   main = Main
#' )
#'
#' # Generate responses from a network of 2 ordinal and 3 Blume-Capel variables.
#' num_variables = 5
#' num_categories = 4
#'
#' Pairwise = matrix(0, nrow = num_variables, ncol = num_variables)
#' Pairwise[2, 1] = Pairwise[4, 1] = Pairwise[3, 2] =
#'   Pairwise[5, 2] = Pairwise[5, 4] = .25
#' Pairwise = Pairwise + t(Pairwise)
#'
#' Main = matrix(NA, num_variables, num_categories)
#' Main[, 1] = -1
#' Main[, 2] = -1
#' Main[3, ] = sort(-abs(rnorm(4)), decreasing = TRUE)
#' Main[5, ] = sort(-abs(rnorm(4)), decreasing = TRUE)
#'
#' x = simulate_mrf(
#'   num_states = 1e3,
#'   num_variables = num_variables,
#'   num_categories = num_categories,
#'   pairwise = Pairwise,
#'   main = Main,
#'   variable_type = c("b", "b", "o", "b", "o"),
#'   baseline_category = 2
#' )
#'
#' # Generate responses from a Gaussian graphical model (GGM) with 4 variables.
#' num_variables = 4
#'
#' # Precision matrix (symmetric, positive definite)
#' Omega = diag(c(1, 1.2, 0.8, 1.5))
#' Omega[2, 1] = Omega[1, 2] = 0.3
#' Omega[3, 1] = Omega[1, 3] = 0.3
#' Omega[4, 2] = Omega[2, 4] = -0.2
#'
#' x = simulate_mrf(
#'   num_states = 500,
#'   num_variables = num_variables,
#'   pairwise = Omega,
#'   variable_type = "continuous"
#' )
#'
#' @seealso \code{\link{simulate.bgms}} for simulating from a fitted model.
#' @family prediction
#'
#' @export
simulate_mrf = function(num_states,
                        num_variables,
                        num_categories,
                        pairwise,
                        main,
                        variable_type = "ordinal",
                        baseline_category,
                        iter = 1e3,
                        seed = NULL) {
  # Check num_states, num_variables ------
  check_positive_integer(num_states, "num_states")
  check_positive_integer(num_variables, "num_variables")

  # Check variable specification -----------------------------------------------
  vt = validate_variable_types(
    variable_type = variable_type,
    num_variables = num_variables,
    allow_continuous = TRUE,
    caller = "simulate_mrf"
  )
  variable_type = vt$variable_type
  is_continuous = vt$is_continuous

  # Blume-Capel binary guard (simulation-specific: uses num_categories)
  if(any(variable_type == "blume-capel")) {
    bc_binary = variable_type == "blume-capel" & num_categories < 2
    if(any(bc_binary)) {
      stop(paste0(
        "The Blume-Capel model only works for ordinal ",
        "variables with more than two \n",
        "response options. But variables ",
        paste(which(bc_binary), collapse = ", "),
        " are binary variables."
      ))
    }
  }

  # ===========================================================================
  #   Continuous (GGM) path --- direct multivariate normal sampling
  # ===========================================================================
  if(is_continuous) {
    # Check pairwise (full precision matrix, including diagonal)
    if(!inherits(pairwise, what = "matrix")) {
      pairwise = as.matrix(pairwise)
    }
    # NAs indicate excluded edges (zero precision)
    pairwise[is.na(pairwise)] = 0
    if(!isSymmetric(pairwise)) {
      stop("The matrix 'pairwise' needs to be symmetric.")
    }
    if(nrow(pairwise) != num_variables) {
      stop(
        "The matrix 'pairwise' needs to have ",
        "'num_variables' rows and columns."
      )
    }
    if(any(diag(pairwise) <= 0)) {
      stop("The diagonal of the precision matrix 'pairwise' must be positive.")
    }

    precision = pairwise

    # Handle means (from 'main', default to zero)
    if(missing(main)) {
      means = rep(0, num_variables)
    } else {
      means = as.numeric(main)
      if(length(means) != num_variables) {
        stop(paste0(
          "'main' must have ", num_variables,
          " elements (one mean per variable), but has ",
          length(means), "."
        ))
      }
      if(any(!is.finite(means))) {
        stop("All elements of 'main' must be finite.")
      }
    }

    # Handle seed
    seed = check_seed(seed)

    x = sample_ggm_direct(
      num_states = num_states,
      precision = precision,
      means = means,
      seed = seed
    )

    return(x)
  }

  # ===========================================================================
  #   Ordinal / Blume-Capel path --- Gibbs sampling
  # ===========================================================================
  check_positive_integer(iter, "iter")

  # Check num_categories ------
  if(length(num_categories) == 1) {
    not_pos_int = num_categories <= 0 ||
      abs(num_categories - round(num_categories)) > .Machine$double.eps
    if(not_pos_int) {
      stop("``num_categories'' needs be a (vector of) positive integer(s).")
    }
    num_categories = rep(num_categories, num_variables)
  } else {
    for(variable in 1:num_variables) {
      nc = num_categories[variable]
      not_pos_int = nc <= 0 || abs(nc - round(nc)) > .Machine$double.eps
      if(not_pos_int) {
        stop(paste(
          "For variable", variable,
          "``num_categories'' was not a",
          "positive integer."
        ))
      }
    }
  }

  # Check the baseline_category for Blume-Capel variables ---------------------
  if(any(variable_type == "blume-capel")) {
    if(length(baseline_category) == 1) {
      baseline_category = rep(baseline_category, num_variables)
    }
    bc_diff = abs(baseline_category - round(baseline_category))
    not_valid = any(baseline_category < 0) ||
      any(bc_diff > .Machine$double.eps)
    if(not_valid) {
      stop(paste0(
        "For variables ",
        which(baseline_category < 0),
        " ``baseline_category'' was either negative or not integer."
      ))
    }
    if(any(baseline_category - num_categories > 0)) {
      stop(paste0(
        "For variables ",
        which(baseline_category - num_categories > 0),
        " the ``baseline_category'' category was larger",
        " than the maximum category value."
      ))
    }
  }

  # Check pairwise ---------------------------------------------------------
  if(!inherits(pairwise, what = "matrix")) {
    pairwise = as.matrix(pairwise)
  }
  if(!isSymmetric(pairwise)) {
    stop("The matrix ``pairwise'' needs to be symmetric.")
  }
  if(nrow(pairwise) != num_variables) {
    stop(
      "The matrix ``pairwise'' needs to have",
      " ``num_variables'' rows and columns."
    )
  }

  # Check the threshold values -------------------------------------------------
  if(!inherits(main, what = "matrix")) {
    if(max(num_categories) == 1) {
      if(length(main) == num_variables) {
        main = matrix(main, ncol = 1)
      } else {
        stop(paste0(
          "The matrix ``main'' has ",
          length(main),
          " elements, but requires",
          num_variables,
          "."
        ))
      }
    } else {
      stop("``main'' needs to be a matrix.")
    }
  }

  if(nrow(main) != num_variables) {
    stop("The matrix ``main'' needs to be have ``num_variables'' rows.")
  }

  for(variable in 1:num_variables) {
    if(variable_type[variable] != "blume-capel") {
      if(anyNA(main[variable, 1:num_categories[variable]])) {
        na_cats = which(is.na(main[variable, 1:num_categories[variable]]))
        stop(paste0(
          "The matrix ``main'' contains NA(s) for variable ",
          variable,
          " in category \n",
          "(categories) ",
          paste(na_cats, collapse = ", "),
          ", where a numeric value is needed."
        ))
      }
      if(ncol(main) > num_categories[variable]) {
        if(!anyNA(main[variable, (num_categories[variable] + 1):ncol(main)])) {
          warning(paste0(
            "The matrix ``main'' contains numeric values for variable ",
            variable,
            " for category \n",
            "(categories, i.e., columns) exceding the maximum of ",
            num_categories[variable],
            ". These values will \n",
            "be ignored."
          ))
        }
      }
    } else {
      if(anyNA(main[variable, 1:2])) {
        stop(paste0(
          "The Blume-Capel model is chosen for the ",
          "category thresholds of variable ",
          variable,
          ". \n",
          "This model has two parameters that need ",
          "to be placed in columns 1 and 2, row \n",
          variable,
          ", of the ``main'' input matrix. ",
          "Currently, there are NA(s) in these \n",
          "entries, where a numeric value is needed."
        ))
      }
      if(ncol(main) > 2) {
        if(!anyNA(main[variable, 3:ncol(main)])) {
          warning(paste0(
            "The Blume-Capel model is chosen for ",
            "the category thresholds of variable ",
            variable,
            ". \n",
            "This model has two parameters that ",
            "need to be placed in columns 1 and ",
            "2, row \n",
            variable,
            ", of the ``main'' input matrix. ",
            "However, there are numeric values \n",
            "in higher categories. These values will be ignored."
          ))
        }
      }
    }
  }

  for(variable in 1:num_variables) {
    if(variable_type[variable] != "blume-capel") {
      for(category in 1:num_categories[variable]) {
        if(!is.finite(main[variable, category])) {
          stop(paste(
            "The threshold parameter for variable", variable, "and category",
            category, "is NA or not finite."
          ))
        }
      }
    } else {
      if(!is.finite(main[variable, 1])) {
        stop(paste0(
          "The alpha parameter for the Blume-Capel model for variable ",
          variable,
          " is NA \n",
          " or not finite."
        ))
      }
      if(!is.finite(main[variable, 2])) {
        stop(paste0(
          "The beta parameter for the Blume-Capel model for variable",
          variable,
          "is NA \n",
          " or not finite."
        ))
      }
    }
  }

  # Handle seed ----------------------------------------------------------------
  seed = check_seed(seed)

  # The Gibbs sampler ----------------------------------------------------------
  if(!any(variable_type == "blume-capel")) {
    x = sample_omrf_gibbs(
      num_states = num_states,
      num_variables = num_variables,
      num_categories = num_categories,
      pairwise = pairwise,
      main = main,
      iter = iter,
      seed = seed
    )
  } else {
    x = sample_bcomrf_gibbs(
      num_states = num_states,
      num_variables = num_variables,
      num_categories = num_categories,
      pairwise = pairwise,
      main = main,
      variable_type_r = variable_type,
      baseline_category = baseline_category,
      iter = iter,
      seed = seed
    )
  }

  return(x)
}


# ==============================================================================
#   mrfSampler() - Deprecated Wrapper for simulate_mrf()
# ==============================================================================

#' @title Sample observations from the ordinal MRF
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `mrfSampler()` was renamed to [simulate_mrf()] as of bgms 0.1.6.3 to
#' follow the package's naming conventions.
#'
#' @inheritParams simulate_mrf
#'
#' @return A matrix of simulated observations (see [simulate_mrf()]).
#'
#' @seealso [simulate_mrf()] for the current function.
#'
#' @keywords internal
#' @export
mrfSampler = function(num_states,
                      num_variables,
                      num_categories,
                      pairwise,
                      main,
                      variable_type = "ordinal",
                      baseline_category,
                      iter = 1e3,
                      seed = NULL) {
  lifecycle::deprecate_warn("0.1.6.3", "mrfSampler()", "simulate_mrf()")

  simulate_mrf(
    num_states = num_states,
    num_variables = num_variables,
    num_categories = num_categories,
    pairwise = pairwise,
    main = main,
    variable_type = variable_type,
    baseline_category = baseline_category,
    iter = iter,
    seed = seed
  )
}

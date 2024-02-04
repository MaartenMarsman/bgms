#' Sample states of the ordinal MRF
#'
#' This function samples states from the ordinal MRF using a Gibbs sampler. The
#' Gibbs sampler is initiated with random values from the response options,
#' after which it proceeds by simulating states for each variable from a logistic
#' model using the other variable states as predictor variables.
#'
#' There are two modeling options for the category thresholds. The default
#' option assumes that the category thresholds are free, except that the first
#' threshold is set to zero for identification. The user then only needs to
#' specify the thresholds for the remaining response categories. This option is
#' useful for any type of ordinal variable and gives the user the most freedom
#' in specifying their model.
#'
#' The Blume-Capel option is specifically designed for ordinal variables that
#' have a special type of reference_category category, such as the neutral
#' category in a Likert scale. The Blume-Capel model specifies the following
#' quadratic model for the threshold parameters:
#' \deqn{\mu_{\text{c}} = \alpha \times \text{c} + \beta \times (\text{c} - \text{r})^2,}{{\mu_{\text{c}} = \alpha \times \text{c} + \beta \times (\text{c} - \text{r})^2,}}
#' where \eqn{\mu_{\text{c}}}{\mu_{\text{c}}} is the threshold for category c
#' (which now includes zero), \eqn{\alpha}{\alpha} offers a linear trend
#' across categories (increasing threshold values if
#' \eqn{\alpha > 0}{\alpha > 0} and decreasing threshold values if
#' \eqn{\alpha <0}{\alpha <0}), if \eqn{\beta < 0}{\beta < 0}, it offers an
#' increasing penalty for responding in a category further away from the
#' reference_category category r, while \eqn{\beta > 0}{\beta > 0} suggests a
#' preference for responding in the reference_category category.
#'
#' @param no_states The number of states of the ordinal MRF to be generated.
#'
#' @param no_variables The number of variables in the ordinal MRF.
#'
#' @param no_categories Either a positive integer or a vector of positive
#' integers of length \code{no_variables}. The number of response categories on top
#' of the base category: \code{no_categories = 1} generates binary states.
#'
#' @param interactions A symmetric \code{no_variables} by \code{no_variables} matrix of
#' pairwise interactions. Only its off-diagonal elements are used.
#'
#' @param thresholds A \code{no_variables} by \code{max(no_categories)} matrix of
#' category thresholds. The elements in row \code{i} indicate the thresholds of
#' variable \code{i}. If \code{no_categories} is a vector, only the first
#' \code{no_categories[i]} elements are used in row \code{i}. If the Blume-Capel
#' model is used for the category thresholds for variable \code{i}, then row
#' \code{i} requires two values (details below); the first is
#' \eqn{\alpha}{\alpha}, the linear contribution of the Blume-Capel model and
#' the second is \eqn{\beta}{\beta}, the quadratic contribution.
#'
#' @param variable_type What kind of variables are simulated? Can be a single
#' character string specifying the variable type of all \code{p} variables at
#' once or a vector of character strings of length \code{p} specifying the type
#' for each variable separately. Currently, bgm supports ``ordinal'' and
#' ``blume-capel''. Binary variables are automatically treated as ``ordinal’’.
#' Defaults to \code{variable_type = "ordinal"}.
#'
#' @param reference_category An integer vector of length \code{no_variables} specifying the
#' reference_category category that is used for the Blume-Capel model (details below).
#' Can be any integer value between \code{0} and \code{no_categories} (or
#' \code{no_categories[i]}).
#'
#' @param iter The number of iterations used by the Gibbs sampler.
#' The function provides the last state of the Gibbs sampler as output. By
#' default set to \code{1e3}.
#'
#' @return A \code{no_states} by \code{no_variables} matrix of simulated states of
#' the ordinal MRF.
#'
#' @examples
#' # Generate responses from a network of five binary and ordinal variables.
#' no_variables = 5
#' no_categories = sample(1:5, size = no_variables, replace = TRUE)
#'
#' Interactions = matrix(0, nrow = no_variables, ncol = no_variables)
#' Interactions[2, 1] = Interactions[4, 1] = Interactions[3, 2] =
#'   Interactions[5, 2] = Interactions[5, 4] = .25
#' Interactions = Interactions + t(Interactions)
#' Thresholds = matrix(0, nrow = no_variables, ncol = max(no_categories))
#' x = mrfSampler(no_states = 1e3,
#'                no_variables = no_variables,
#'                no_categories = no_categories,
#'                interactions = Interactions,
#'                thresholds = Thresholds)
#' @export
mrfSampler = function(no_states,
                      no_variables,
                      no_categories,
                      interactions,
                      thresholds,
                      variable_type = "ordinal",
                      reference_category,
                      iter = 1e3) {
  # Check no_states, no_variables, iter --------------------------------------------
  if(no_states <= 0 ||
     abs(no_states - round(no_states)) > .Machine$double.eps)
    stop("``no_states'' needs be a positive integer.")
  if(no_variables <= 0 ||
     abs(no_variables - round(no_variables)) > .Machine$double.eps)
    stop("``no_variables'' needs be a positive integer.")
  if(iter <= 0 ||
     abs(iter - round(iter)) > .Machine$double.eps)
    stop("``iter'' needs be a positive integer.")

  # Check no_categories --------------------------------------------------------
  if(length(no_categories) == 1) {
    if(no_categories <= 0 ||
       abs(no_categories - round(no_categories)) > .Machine$double.eps)
      stop("``no_categories'' needs be a (vector of) positive integer(s).")
    no_categories = rep(no_categories, no_variables)
  } else {
    for(variable in 1:no_variables) {
      if(no_categories[variable] <= 0 ||
         abs(no_categories[variable] - round(no_categories[variable])) >
         .Machine$double.eps)
        stop(paste("For variable", variable, "``no_categories'' was not a positive integer."))
    }
  }

  # Check variable specification -----------------------------------------------
  if(length(variable_type) == 1) {
    variable_type = match.arg(arg = variable_type,
                              choices = c("ordinal", "blume-capel"))

    if(variable_type == "blume-capel" && any(no_categories < 2)) {
      stop(paste0("The Blume-Capel model only works for ordinal variables with more than two \n",
                  "response options. But variables ", which(no_categories < 2), " are binary variables."))
    }
    variable_type = rep(variable_type, no_variables)

  } else {
    if(length(variable_type) != no_variables) {
      stop(paste0("The argument ``variable_type'' should be either a single character string or a \n",
                  "vector of character strings of length ``no_variables''."))
    } else {
      for(variable in 1:no_variables)
        variable_type[variable] = match.arg(arg = variable_type[variable],
                                            choices = c("ordinal", "blume-capel"))
      if(any (variable_type == "blume-capel" & no_categories < 2)) {
        stop(paste0("The Blume-Capel model only works for ordinal variables with more than two \n",
                    "response options. But variables ",
                    which(variable_type == "blume-capel" & no_categories < 2),
                    " are binary variables."))
      }
    }
  }

  # Check the reference_category for Blume-Capel variables ---------------------
  if(any(variable_type == "blume-capel")) {
    if(length(reference_category) == 1) {
      reference_category = rep(reference_category, no_variables)
    }
    if(any(reference_category < 0) || any(abs(reference_category - round(reference_category)) > .Machine$double.eps)) {
      stop(paste0("For variables ",
                  which(reference_category < 0),
                  " ``reference_category'' was either negative or not integer."))
    }
    if(any(reference_category - no_categories > 0)) {
      stop(paste0("For variables ",
                  which(reference_category - no_categories > 0),
                  " the ``reference_category'' category was larger than the maximum category value."))
    }
  }

  # Check interactions ---------------------------------------------------------
  if(!inherits(interactions, what = "matrix"))
    interactions = as.matrix(interactions)
  if(!isSymmetric(interactions))
    stop("The matrix ``interactions'' needs to be symmetric.")
  if(nrow(interactions) != no_variables)
    stop("The matrix ``interactions'' needs to have ``no_variables'' rows and columns.")

  # Check the threshold values -------------------------------------------------
  if(!inherits(thresholds, what = "matrix")) {
    if(max(no_categories) == 1) {
      if(length(thresholds) == no_variables) {
        thresholds = matrix(thresholds, ncol = 1)
      } else {
        stop(paste0("The matrix ``thresholds'' has ",
                    length(thresholds),
                    " elements, but requires",
                    no_variables,
                    "."))
      }
    } else {
      stop("``Thresholds'' needs to be a matrix.")
    }
  }

  if(nrow(thresholds) != no_variables)
    stop("The matrix ``thresholds'' needs to be have ``no_variables'' rows.")

  for(variable in 1:no_variables) {
    if(variable_type[variable] != "blume-capel") {
      if(any(is.na(thresholds[variable, 1:no_categories[variable]]))) {
        tmp = which(is.na(thresholds[variable, 1:no_categories[variable]]))

        string = paste(tmp, sep = ",")

        stop(paste0("The matrix ``thresholds'' contains NA(s) for variable ",
                    variable,
                    " in category \n",
                    "(categories) ",
                    paste(which(is.na(thresholds[variable, 1:no_categories[variable]])), collapse = ", "),
                    ", where a numeric value is needed."))
      }
      if(ncol(thresholds) > no_categories[variable]) {
        if(any(!is.na(thresholds[variable, (no_categories[variable]+1):ncol(thresholds)]))){
          warning(paste0("The matrix ``thresholds'' contains numeric values for variable ",
                         variable,
                         " for category \n",
                         "(categories, i.e., columns) exceding the maximum of ",
                         no_categories[variable],
                         ". These values will \n",
                         "be ignored."))
        }
      }
    } else {
      if(any(is.na(thresholds[variable, 1:2]))) {
        stop(paste0("The Blume-Capel model is chosen for the category thresholds of variable ",
                    variable,
                    ". \n",
                    "This model has two parameters that need to be placed in columns 1 and 2, row \n",
                    variable,
                    ", of the ``thresholds'' input matrix. Currently, there are NA(s) in these \n",
                    "entries, where a numeric value is needed."))
      }
      if(ncol(thresholds) > 2) {
        if(any(!is.na(thresholds[variable, 3:ncol(thresholds)]))){
          warning(paste0("The Blume-Capel model is chosen for the category thresholds of variable ",
                         variable,
                         ". \n",
                         "This model has two parameters that need to be placed in columns 1 and 2, row \n",
                         variable,
                         ", of the ``thresholds'' input matrix. However, there are numeric values \n",
                         "in higher categories. These values will be ignored."))
        }
      }
    }

  }

  for(variable in 1:no_variables) {
    if(variable_type[variable] != "blume-capel") {
      for(category in 1:no_categories[variable]) {
        if(!is.finite(thresholds[variable, category]))
          stop(paste("The threshold parameter for variable", variable, "and category",
                     category, "is NA or not finite."))
      }
    } else {
      if(!is.finite(thresholds[variable, 1]))
        stop(paste0(
          "The alpha parameter for the Blume-Capel model for variable ",
          variable,
          " is NA \n",
          " or not finite."))
      if(!is.finite(thresholds[variable, 2]))
        stop(paste0("The beta parameter for the Blume-Capel model for variable",
                    variable,
                    "is NA \n",
                    " or not finite."))
    }
  }

  # The Gibbs sampler ----------------------------------------------------------
  if(!any(variable_type == "blume-capel")) {
    x <- sample_omrf_gibbs(no_states = no_states,
                           no_variables = no_variables,
                           no_categories = no_categories,
                           interactions = interactions,
                           thresholds = thresholds,
                           iter = iter)
  } else {
    x <- sample_bcomrf_gibbs(no_states = no_states,
                             no_variables = no_variables,
                             no_categories = no_categories,
                             interactions = interactions,
                             thresholds = thresholds,
                             variable_type = variable_type,
                             reference_category = reference_category,
                             iter = iter)
  }

  return(x)
}

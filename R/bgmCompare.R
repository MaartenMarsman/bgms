#' @title Bayesian Estimation and Variable Selection for Group Differences in Markov Random Fields
#'
#' @description
#' The \code{bgmCompare} function estimates group differences in category
#' threshold parameters (main effects) and pairwise interactions (pairwise
#' effects) of a Markov Random Field (MRF) for binary and ordinal variables.
#' Groups can be defined either by supplying two separate datasets (\code{x} and
#' \code{y}) or by a group membership vector. Optionally, Bayesian variable
#' selection can be applied to identify differences across groups.
#'
#' @details
#' Group-specific parameters are decomposed into a shared baseline plus
#' group differences that sum to zero. Difference selection uses
#' spike-and-slab priors (Bernoulli or Beta-Bernoulli). Parameters are
#' sampled with NUTS (default) or adaptive Metropolis--Hastings, using the
#' same multi-stage warmup schedule as \code{\link{bgm}}.
#'
#' Groups are numbered \code{1, 2, ...} in the order they first appear in
#' \code{group_indicator}, whatever that vector's storage type: the first row's
#' group is group 1, the first row belonging to some other group is group 2, and
#' so on. (Note that this is first appearance, not sorted order -- an indicator
#' reading \code{c("fr", "fr", "en", ...)} makes \code{"fr"} group 1.) With
#' \code{x} and \code{y} instead, \code{x} is group 1 and \code{y} is group 2.
#' Every output keys on these numbers, and the extractor column names
#' (\code{group1}, \code{group2}) stay numeric. The original labels are carried
#' into the fit and shown on the displays a person reads -- the \code{print} and
#' \code{summary} headers, plot panel titles, calibration panel titles, and
#' centrality labels -- so that, for example, group 2 prints as
#' \code{"group 2 (en)"}. Fits made with earlier versions of \pkg{bgms} carry no
#' labels and display the bare numbers.
#'
#' For full details on model specification, prior choices, and output
#' interpretation, see the package website at
#' \url{https://bayesian-graphical-modelling-lab.github.io/bgms-docs/}.
#'
#' @section Categories across groups:
#'
#' Groups being compared often differ in which categories of an ordinal
#' variable they actually use --- a clinical group may never give the lowest
#' answer, a control group never the highest. \code{bgmCompare()} models the
#' \strong{union} of the categories observed across the groups: a category that
#' any group uses is kept for all of them. Only a category value that
#' \emph{no} group uses is dropped, after which the remaining categories are
#' renumbered contiguously from 0. That renumbering is reported with a
#' \code{message()}; no category anyone observed is ever merged into another.
#'
#' Keeping a category that some group never uses has a consequence worth
#' knowing. That group contributes no observations to the category, so its
#' data say nothing about where its threshold for that category lies, and the
#' estimated difference for that group-by-category combination is determined
#' by the prior rather than by the data: it will be large and very uncertain,
#' and it is not evidence of a group difference. \code{bgmCompare()} raises a
#' \code{warning()} naming every variable, category, and group this affects,
#' and records the per-group category counts in the fitted object
#' (\code{extract_arguments(fit)$category_support}) so they can be checked
#' afterwards. Only the category thresholds are affected this way; the
#' pairwise (edge) parameters and their differences are estimated from all the
#' data and are not.
#'
#' Blume--Capel variables are exempt from all of this. Their two parameters
#' are functions of the numeric category \emph{score}, so renumbering the
#' categories would change the model rather than relabel it, and a score no
#' one happened to observe is still a meaningful point on the scale. Their
#' categories are therefore always retained exactly as supplied.
#'
#' @seealso \code{vignette("comparison", package = "bgms")} for a worked example.
#' @family model-fitting
#'
#' @param x A data frame or matrix of binary and ordinal responses for
#'   Group 1. Variables should be coded as nonnegative integers starting at
#'   0. See the \emph{Categories across groups} section for how the category
#'   codes of an ordinal variable are treated when the groups do not observe
#'   the same ones.
#' @param y Optional data frame or matrix for Group 2 (two-group designs).
#'   Must have the same variables (columns) as \code{x}.
#' @param group_indicator Optional integer vector of group memberships for
#'   rows of \code{x} (multi-group designs). Ignored if \code{y} is supplied.
#' @param difference_selection Logical. If \code{TRUE}, spike-and-slab priors
#'   are applied to difference parameters. Default: \code{TRUE}.
#' @param main_difference_selection Logical. If \code{TRUE}, apply spike-and-slab
#'   selection to main effect (threshold) differences. If \code{FALSE}, main
#'   effect differences are always included (no selection). Since main effects
#'   are often nuisance parameters and their selection can interfere with
#'   pairwise selection under the Beta-Bernoulli prior, the default is
#'   \code{FALSE}. Only used when \code{difference_selection = TRUE}.
#' @param variable_type Character vector specifying type of each variable:
#'   \code{"ordinal"} (default) or \code{"blume-capel"}.
#' @param baseline_category Integer or vector giving the baseline category
#'   for Blume--Capel variables.
#' @param difference_scale Double. Scale of the prior for difference
#'   parameters. Default: \code{1}.
#' @param difference_family Character. Distributional family of the prior on
#'   difference parameters, one of \code{"Normal"} (default) or \code{"Cauchy"}.
#'   Governs both the pairwise-interaction differences and the main-effect
#'   (threshold) differences; under \code{difference_selection = TRUE} it is
#'   the slab of the spike-and-slab. Independent of \code{interaction_prior},
#'   which governs the baseline interactions.
#' @param difference_prior An indicator prior specification object for
#'   difference selection, created by one of:
#'   \itemize{
#'     \item \code{\link{bernoulli_prior}()}: Fixed inclusion probability (default).
#'     \item \code{\link{beta_bernoulli_prior}()}: Beta-distributed inclusion.
#'     \item \code{\link{sbm_prior}()}: Stochastic Block Model.
#'   }
#'   Legacy character strings \code{"Bernoulli"} and \code{"Beta-Bernoulli"}
#'   are still accepted but deprecated.
#'   Default: \code{bernoulli_prior(0.5)}.
#' @param difference_probability `r lifecycle::badge("deprecated")` Numeric.
#'   Use \code{difference_prior = bernoulli_prior(probability)} instead.
#'   Default: \code{0.5}.
#' @param interaction_prior A prior specification object for baseline pairwise
#'   interaction parameters, created by one of the prior constructor functions:
#'   \itemize{
#'     \item \code{\link{normal_prior}()}: Normal(0, scale) prior (default).
#'     \item \code{\link{cauchy_prior}()}: Cauchy(0, scale) prior.
#'   }
#'   When supplied, overrides \code{pairwise_scale}.
#'   Default: \code{normal_prior(scale = 1)}, matching \code{\link{bgm}}.
#'   Governs the baseline pairwise interactions only; the group differences
#'   are governed by \code{difference_family} and \code{difference_scale}.
#' @param threshold_prior A prior specification object for threshold (main
#'   effect) parameters, created by one of the prior constructor functions:
#'   \itemize{
#'     \item \code{\link{beta_prime_prior}()}: Beta-prime prior (default).
#'     \item \code{\link{normal_prior}()}: Normal(0, scale) prior.
#'   }
#'   When supplied, overrides \code{main_alpha} and \code{main_beta}.
#'   Default: \code{beta_prime_prior(alpha = 0.5, beta = 0.5)}.
#' @param beta_bernoulli_alpha,beta_bernoulli_beta Doubles. Shape parameters
#'   of the Beta prior for inclusion probabilities in the Beta--Bernoulli
#'   model. Defaults: \code{1}.
#' @param pairwise_scale `r lifecycle::badge("deprecated")` Double. Scale of the
#'   baseline pairwise interaction prior. Retained for backward compatibility,
#'   it sets a Cauchy prior at that scale, which is not the current default.
#'   Use \code{interaction_prior} instead.
#' @param main_alpha,main_beta `r lifecycle::badge("deprecated")` Doubles. Shape
#'   parameters of the beta-prime prior for baseline threshold parameters.
#'   Use \code{threshold_prior = beta_prime_prior(alpha, beta)} instead.
#' @param iter Integer. Number of post--warmup iterations per chain.
#'   Default: \code{2e3}.
#' @param warmup Integer. Number of warmup iterations before sampling.
#'   Default: \code{2e3}.
#' @param na_action Character. How to handle missing data:
#'   \code{"listwise"} (drop rows) or \code{"impute"} (impute within Gibbs).
#'   Default: \code{"listwise"}.
#' @param display_progress Character. Controls progress reporting:
#'   \code{"per-chain"}, \code{"total"}, or \code{"none"}.
#'   Default: \code{"per-chain"}.
#' @param progress_callback An optional R function with signature
#'   \code{function(completed, total)} that is called at regular intervals
#'   during sampling, where \code{completed} is the number of iterations
#'   completed across all chains and \code{total} is the total number of
#'   iterations. Useful for external front-ends (e.g., JASP) that supply
#'   their own progress reporting.
#'   When \code{NULL} (the default), no callback is invoked.
#' @param verbose Logical. If \code{TRUE}, prints informational messages
#'   during data processing (e.g., missing data handling, variable recoding).
#'   Defaults to \code{getOption("bgms.verbose", TRUE)}. Set
#'   \code{options(bgms.verbose = FALSE)} to suppress messages globally.
#' @param update_method Character. Sampling algorithm:
#'   \code{"adaptive-metropolis"} or \code{"nuts"}. Default: \code{"nuts"}.
#' @param target_accept Numeric between 0 and 1. Target acceptance rate.
#'   Defaults: 0.44 (Metropolis), 0.80 (NUTS).
#' @param nuts_max_depth Integer. Maximum tree depth for NUTS. Default: \code{10}.
#' @param learn_mass_matrix Logical. If \code{TRUE}, adapts a diagonal mass
#' matrix during warmup (NUTS only). Default: \code{TRUE}.
#' @param chains Integer. Number of parallel chains. Default: \code{4}.
#' @param cores Integer. Number of CPU cores. Default:
#'   \code{parallel::detectCores()}.
#' @param seed Optional integer. Random seed for reproducibility.
#' @param main_difference_model,reference_category,pairwise_difference_scale,main_difference_scale,pairwise_difference_prior,main_difference_prior,pairwise_difference_probability,main_difference_probability,pairwise_beta_bernoulli_alpha,pairwise_beta_bernoulli_beta,main_beta_bernoulli_alpha,main_beta_bernoulli_beta,interaction_scale,threshold_alpha,threshold_beta,burnin,save
#'   `r lifecycle::badge("deprecated")`
#'   Deprecated arguments as of \strong{bgms 0.1.6.0}.
#'   Use `difference_scale`, `difference_prior`, `difference_probability`,
#'   `beta_bernoulli_alpha`, `beta_bernoulli_beta`, `baseline_category`,
#'   `interaction_prior`, `threshold_prior`, and `warmup` instead.
#' @param standardize `r lifecycle::badge("deprecated")` Logical. Deprecated
#'   as of \strong{bgms 0.2.0.0}. Through 0.1.6.3, \code{TRUE} adjusted each
#'   pair's baseline and difference prior scale by the product of the two
#'   variables' maximum scores. Pairwise interactions are now on the
#'   association scale and share one prior scale, so the per-pair adjustment is
#'   gone: \code{standardize = FALSE} (the old default) warns and proceeds,
#'   while \code{standardize = TRUE} errors and points to setting the scales
#'   directly through \code{interaction_prior} and \code{difference_scale}.
#' @return
#' An S7 object of class \code{bgmCompare} supporting list-style \code{$} /
#' \code{[[} access for backward compatibility, containing posterior summaries,
#' posterior mean matrices, and raw MCMC samples. Some
#' \code{posterior_summary_*} fields are computed lazily on first access:
#' \itemize{
#'   \item \code{posterior_summary_main_baseline},
#'     \code{posterior_summary_pairwise_baseline}: summaries of baseline
#'     thresholds and pairwise interactions.
#'   \item \code{posterior_summary_main_differences},
#'     \code{posterior_summary_pairwise_differences}: summaries of group
#'     differences in thresholds and pairwise interactions.
#'   \item \code{posterior_summary_indicator}: summaries of inclusion
#'     indicators (if \code{difference_selection = TRUE}).
#'   \item \code{posterior_mean_main_baseline},
#'     \code{posterior_mean_pairwise_baseline}: posterior mean matrices
#'     (legacy style).
#'   \item \code{raw_samples}: list of raw draws per chain for main,
#'     pairwise, and indicator parameters.
#'   \item \code{arguments}: list of function call arguments and metadata.
#' }
#'
#' The \code{summary()} method prints formatted summaries, and
#' \code{coef()} extracts posterior means.
#'
#' NUTS diagnostics (tree depth, divergences, energy, E-BFMI) are included
#' in \code{fit$nuts_diag} if \code{update_method = "nuts"}.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' \dontrun{
#' # Run bgmCompare on subset of the Boredom dataset
#' x = Boredom[Boredom$language == "fr", 2:6]
#' y = Boredom[Boredom$language != "fr", 2:6]
#'
#' fit = bgmCompare(x, y, chains = 2)
#'
#' # Posterior inclusion probabilities
#' summary(fit)$indicator
#'
#' # Bayesian model averaged main effects for the groups
#' coef(fit)$main_effects_groups
#'
#' # Bayesian model averaged pairwise effects for the groups
#' coef(fit)$pairwise_effects_groups
#' }
#'
#' @export
bgmCompare = function(
  x,
  y,
  group_indicator,
  difference_selection = TRUE,
  main_difference_selection = FALSE,
  variable_type = "ordinal",
  baseline_category,
  difference_scale = 1,
  difference_family = c("Normal", "Cauchy"),
  difference_prior = bernoulli_prior(0.5),
  difference_probability,
  interaction_prior = normal_prior(scale = 1),
  threshold_prior = beta_prime_prior(alpha = 0.5, beta = 0.5),
  iter = 2e3,
  warmup = 2e3,
  na_action = c("listwise", "impute"),
  update_method = c("nuts", "adaptive-metropolis"),
  target_accept,
  nuts_max_depth = 10,
  learn_mass_matrix = TRUE,
  chains = 4,
  cores = parallel::detectCores(),
  display_progress = c("per-chain", "total", "none"),
  seed = NULL,
  verbose = getOption("bgms.verbose", TRUE),
  progress_callback = NULL,
  # Deprecated prior arguments
  pairwise_scale,
  main_alpha,
  main_beta,
  beta_bernoulli_alpha,
  beta_bernoulli_beta,
  # Deprecated arguments (v0.1.6.0)
  main_difference_model,
  reference_category,
  main_difference_scale,
  pairwise_difference_scale,
  pairwise_difference_prior,
  main_difference_prior,
  pairwise_difference_probability,
  main_difference_probability,
  pairwise_beta_bernoulli_alpha,
  pairwise_beta_bernoulli_beta,
  main_beta_bernoulli_alpha,
  main_beta_bernoulli_beta,
  interaction_scale,
  threshold_alpha,
  threshold_beta,
  burnin,
  save,
  # Deprecated arguments (v0.2.0.0)
  standardize
) {
  # Set verbose option for internal functions, restore on exit
  old_verbose = getOption("bgms.verbose")
  options(bgms.verbose = verbose)
  on.exit(options(bgms.verbose = old_verbose), add = TRUE)

  # bgmCompare supports only NUTS and adaptive-Metropolis.
  update_method = match.arg(update_method)

  if(hasArg(main_difference_model)) {
    lifecycle::deprecate_warn("0.1.6.0", "bgmCompare(main_difference_model =)")
  }

  if(hasArg(reference_category)) {
    lifecycle::deprecate_warn("0.1.6.0", "bgmCompare(reference_category =)", "bgmCompare(baseline_category =)")
    if(!hasArg(baseline_category)) baseline_category = reference_category
  }

  if(hasArg(pairwise_difference_scale) || hasArg(main_difference_scale)) {
    lifecycle::deprecate_warn(
      "0.1.6.0", "bgmCompare(pairwise_difference_scale =, main_difference_scale =)",
      "bgmCompare(difference_scale =)"
    )
    if(!hasArg(difference_scale)) {
      difference_scale = if(!missing(pairwise_difference_scale)) pairwise_difference_scale else main_difference_scale
    }
  }

  if(hasArg(pairwise_difference_prior) || hasArg(main_difference_prior)) {
    lifecycle::deprecate_warn(
      "0.1.6.0",
      "bgmCompare(pairwise_difference_prior =, main_difference_prior =)",
      "bgmCompare(difference_prior =)"
    )
    if(!hasArg(difference_prior)) {
      difference_prior = if(!missing(pairwise_difference_prior)) pairwise_difference_prior else main_difference_prior
    }
  }

  if(hasArg(pairwise_difference_probability) || hasArg(main_difference_probability)) {
    lifecycle::deprecate_warn(
      "0.1.6.0",
      "bgmCompare(pairwise_difference_probability =, main_difference_probability =)",
      "bgmCompare(difference_probability =)"
    )
    if(!hasArg(difference_probability)) {
      difference_probability = if(!missing(pairwise_difference_probability)) pairwise_difference_probability else main_difference_probability
    }
  }

  if(hasArg(pairwise_beta_bernoulli_alpha) || hasArg(main_beta_bernoulli_alpha)) {
    lifecycle::deprecate_warn(
      "0.1.6.0",
      "bgmCompare(pairwise_beta_bernoulli_alpha =, main_beta_bernoulli_alpha =)",
      "bgmCompare(beta_bernoulli_alpha =)"
    )
    if(!hasArg(beta_bernoulli_alpha)) {
      beta_bernoulli_alpha = if(!missing(pairwise_beta_bernoulli_alpha)) pairwise_beta_bernoulli_alpha else main_beta_bernoulli_alpha
    }
  }

  if(hasArg(pairwise_beta_bernoulli_beta) || hasArg(main_beta_bernoulli_beta)) {
    lifecycle::deprecate_warn(
      "0.1.6.0",
      "bgmCompare(pairwise_beta_bernoulli_beta =, main_beta_bernoulli_beta =)",
      "bgmCompare(beta_bernoulli_beta =)"
    )
    if(!hasArg(beta_bernoulli_beta)) {
      beta_bernoulli_beta = if(!missing(pairwise_beta_bernoulli_beta)) pairwise_beta_bernoulli_beta else main_beta_bernoulli_beta
    }
  }

  if(hasArg(interaction_scale)) {
    lifecycle::deprecate_warn(
      "0.1.6.0", "bgmCompare(interaction_scale =)",
      "bgmCompare(interaction_prior =)"
    )
    if(!hasArg(pairwise_scale) &&
      identical(interaction_prior, normal_prior(scale = 1))) {
      interaction_prior = cauchy_prior(scale = interaction_scale)
    }
  }

  if(hasArg(pairwise_scale)) {
    lifecycle::deprecate_warn(
      "0.2.0", "bgmCompare(pairwise_scale =)",
      "bgmCompare(interaction_prior =)"
    )
    if(identical(interaction_prior, normal_prior(scale = 1))) {
      interaction_prior = cauchy_prior(scale = pairwise_scale)
    }
  }

  if(hasArg(threshold_alpha) || hasArg(threshold_beta)) {
    lifecycle::deprecate_warn(
      "0.1.6.0",
      "bgmCompare(threshold_alpha =, threshold_beta =)",
      "bgmCompare(threshold_prior =)"
    )
    if(identical(threshold_prior, beta_prime_prior(0.5, 0.5))) {
      ta = if(hasArg(threshold_alpha)) threshold_alpha else 0.5
      tb = if(hasArg(threshold_beta)) threshold_beta else 0.5
      threshold_prior = beta_prime_prior(alpha = ta, beta = tb)
    }
  }

  if(hasArg(main_alpha) || hasArg(main_beta)) {
    lifecycle::deprecate_warn(
      "0.2.0", "bgmCompare(main_alpha =)",
      "bgmCompare(threshold_prior =)"
    )
    if(identical(threshold_prior, beta_prime_prior(0.5, 0.5))) {
      ma = if(hasArg(main_alpha)) main_alpha else 0.5
      mb = if(hasArg(main_beta)) main_beta else 0.5
      threshold_prior = beta_prime_prior(alpha = ma, beta = mb)
    }
  }

  if(hasArg(burnin)) {
    lifecycle::deprecate_warn("0.1.6.0", "bgmCompare(burnin =)", "bgmCompare(warmup =)")
    if(!hasArg(warmup)) warmup = burnin
  }

  if(hasArg(save)) {
    lifecycle::deprecate_warn("0.1.6.0", "bgmCompare(save =)")
  }

  # --- Legacy deprecation: v0.2.0.0 removals ----------------------------------
  # standardize scaled each pair's baseline and difference prior by the product
  # of the two variables' maximum scores. FALSE, its old default, is what the
  # sampler does now, so it warns and proceeds; TRUE asks for an adjustment that
  # no longer exists, so it stops with the manual alternative.
  if(hasArg(standardize)) {
    if(isFALSE(standardize)) {
      lifecycle::deprecate_warn(
        "0.2.0", "bgmCompare(standardize =)",
        details = paste(
          "FALSE is what the sampler does, so the fit is unaffected; drop the",
          "argument."
        )
      )
    } else {
      lifecycle::deprecate_stop(
        "0.2.0", "bgmCompare(standardize = 'no longer supports TRUE')",
        details = paste(
          "Pairwise interactions are on the association scale and share one",
          "prior scale, so the per-pair adjustment by the maximum score",
          "product (scale * m_i * m_j) has been dropped. Set the scales",
          "directly, e.g. interaction_prior = normal_prior(scale =) for the",
          "baseline and difference_scale for the differences; there is no",
          "per-pair equivalent."
        )
      )
    }
  }

  # --- Handle difference_prior: accept both string (deprecated) and object ------
  if(is.character(difference_prior)) {
    lifecycle::deprecate_warn(
      "0.2.0", "bgmCompare(difference_prior = 'must be a prior object')",
      "bgmCompare(difference_prior = 'bernoulli_prior()')"
    )
    difference_prior_str = match.arg(difference_prior,
      choices = c("Bernoulli", "Beta-Bernoulli", "Stochastic-Block")
    )
    dp_prob = if(hasArg(difference_probability)) difference_probability else 0.5
    bba = if(hasArg(beta_bernoulli_alpha)) beta_bernoulli_alpha else 1
    bbb = if(hasArg(beta_bernoulli_beta)) beta_bernoulli_beta else 1

    difference_prior = switch(difference_prior_str,
      "Bernoulli" = bernoulli_prior(inclusion_probability = dp_prob),
      "Beta-Bernoulli" = beta_bernoulli_prior(alpha = bba, beta = bbb),
      "Stochastic-Block" = sbm_prior(alpha = bba, beta = bbb)
    )
  } else {
    if(hasArg(difference_probability)) {
      lifecycle::deprecate_warn(
        "0.2.0", "bgmCompare(difference_probability =)",
        "bgmCompare(difference_prior = 'bernoulli_prior()')"
      )
    }
  }

  # --- Unpack prior objects to flat parameters ---------------------------------
  ip = unpack_interaction_prior(interaction_prior)
  tp = unpack_threshold_prior(threshold_prior)
  difference_family = match.arg(difference_family)
  difference_prior_type = tolower(difference_family)

  # Unpack difference prior to flat params for bgm_spec. The prior needs the
  # variable count, so x is validated here rather than left to bgm_spec()
  # (data_check() is idempotent; the spec runs it again harmlessly).
  x = data_check(x, "x")
  num_variables = ncol(x)
  dp = unpack_indicator_prior(difference_prior, num_variables)

  # --- Build spec, sample, build output ----------------------------------------
  spec = bgm_spec(
    x = x,
    model_type = "compare",
    variable_type = variable_type,
    baseline_category = if(hasArg(baseline_category)) baseline_category else NULL,
    y = if(hasArg(y)) y else NULL,
    group_indicator = if(hasArg(group_indicator)) group_indicator else NULL,
    na_action = na_action,
    interaction_prior_type = ip$interaction_prior_type,
    pairwise_scale = ip$pairwise_scale,
    interaction_alpha = ip$interaction_alpha,
    interaction_beta = ip$interaction_beta,
    threshold_prior_type = tp$threshold_prior_type,
    main_alpha = tp$main_alpha,
    main_beta = tp$main_beta,
    threshold_scale = tp$threshold_scale,
    difference_selection = difference_selection,
    main_difference_selection = main_difference_selection,
    difference_prior = dp$edge_prior,
    difference_scale = difference_scale,
    difference_prior_type = difference_prior_type,
    difference_probability = dp$inclusion_probability,
    beta_bernoulli_alpha = dp$beta_bernoulli_alpha,
    beta_bernoulli_beta = dp$beta_bernoulli_beta,
    difference_beta_bernoulli_alpha_between = dp$beta_bernoulli_alpha_between,
    difference_beta_bernoulli_beta_between = dp$beta_bernoulli_beta_between,
    difference_dirichlet_alpha = dp$dirichlet_alpha,
    difference_lambda = dp$lambda,
    update_method = update_method,
    target_accept = if(hasArg(target_accept)) target_accept else NULL,
    iter = iter,
    warmup = warmup,
    nuts_max_depth = nuts_max_depth,
    learn_mass_matrix = learn_mass_matrix,
    chains = chains,
    cores = cores,
    seed = seed,
    display_progress = display_progress,
    verbose = verbose,
    progress_callback = progress_callback
  )

  raw = run_sampler(spec)
  output = build_output(spec, raw)

  return(output)
}

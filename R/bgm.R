#' @title Bayesian Estimation or Edge Selection for Markov Random Fields
#'
#' @description
#' The \code{bgm} function estimates the pseudoposterior distribution of the
#' parameters of a Markov Random Field (MRF) for binary, ordinal, continuous,
#' or mixed (discrete and continuous) variables. Optionally, it performs
#' Bayesian edge selection using discrete spike-and-slab priors to infer the
#' network structure.
#'
#' @details
#' Depending on the variable types, the model is an ordinal MRF, a Gaussian
#' graphical model (GGM), or a mixed MRF. Both regular ordinal variables and
#' Blume--Capel ordinal variables (with a baseline category) are supported.
#'
#' Edge selection uses spike-and-slab priors with Bernoulli, Beta-Bernoulli,
#' or Stochastic-Block priors on the edge inclusion indicators. Parameters are
#' sampled with NUTS (default) or adaptive Metropolis--Hastings, with a
#' multi-stage warmup schedule. Missing data can be handled via listwise
#' deletion or Gibbs imputation.
#'
#' For full details on model specification, prior choices, warmup, and
#' output interpretation, see the package website at
#' \url{https://bayesian-graphical-modelling-lab.github.io/bgms-docs/}.
#'
#' @seealso \code{vignette("intro", package = "bgms")} for a worked example.
#' @family model-fitting
#'
#' @param x A data frame or matrix with \code{n} rows and \code{p} columns.
#'   Columns may contain binary, ordinal, or continuous variables (see
#'   \code{variable_type}). Discrete variables are automatically recoded to
#'   non-negative integers (\code{0, 1, ..., m}); for regular ordinal
#'   variables, unobserved categories are collapsed, while Blume--Capel
#'   variables retain all categories. Continuous variables are column-centered
#'   internally so that the GGM likelihood is formulated with a zero-mean
#'   assumption.
#'
#' @param variable_type Character or character vector. Specifies the type of
#'   each variable in \code{x}. Allowed values: \code{"ordinal"},
#'   \code{"blume-capel"}, or \code{"continuous"}. A single string applies
#'   to all variables. A per-variable vector that mixes discrete
#'   (\code{"ordinal"} / \code{"blume-capel"}) and \code{"continuous"}
#'   types fits a mixed MRF. Binary variables are automatically treated as
#'   \code{"ordinal"}. Default: \code{"ordinal"}.
#'
#' @param baseline_category Integer or vector. Baseline category used in
#'   Blume--Capel variables. Can be a single integer (applied to all) or a
#'   vector of length \code{p}. Required if at least one variable is of type
#'   \code{"blume-capel"}.
#'
#' @param iter Integer. Number of post--burn-in iterations (per chain).
#'   Default: \code{2e3}.
#'
#' @param warmup Integer. Number of warmup iterations before collecting
#'   samples. Short warmups trigger progressive warnings (NUTS only); see
#'   \code{validate_sampler()} for the thresholds. Default: \code{2e3}.
#'
#' @param interaction_prior A prior specification object for pairwise
#'   interaction parameters, created by one of the prior constructor functions:
#'   \itemize{
#'     \item \code{\link{cauchy_prior}()}: Cauchy(0, scale) prior (default).
#'     \item \code{\link{normal_prior}()}: Normal(0, scale) prior.
#'     \item \code{\link{beta_prime_prior}()}: Beta-prime prior.
#'   }
#'   Default: \code{cauchy_prior(scale = 1)}.
#'
#' @param threshold_prior A prior specification object for threshold (main
#'   effect) parameters, created by one of the prior constructor functions:
#'   \itemize{
#'     \item \code{\link{beta_prime_prior}()}: Beta-prime prior (default).
#'     \item \code{\link{cauchy_prior}()}: Cauchy(0, scale) prior.
#'     \item \code{\link{normal_prior}()}: Normal(0, scale) prior.
#'   }
#'   Default: \code{beta_prime_prior(alpha = 0.5, beta = 0.5)}.
#'
#' @param means_prior A prior specification object for continuous variable
#'   means (mixed MRF models only), created by one of the prior constructor
#'   functions:
#'   \itemize{
#'     \item \code{\link{normal_prior}()}: Normal(0, scale) prior (default).
#'     \item \code{\link{cauchy_prior}()}: Cauchy(0, scale) prior.
#'     \item \code{\link{beta_prime_prior}()}: Beta-prime prior.
#'   }
#'   Only used when the model includes continuous variables. Ignored for
#'   pure ordinal or pure continuous (GGM) models.
#'   Default: \code{normal_prior(scale = 1)}.
#'
#' @param precision_scale_prior A prior specification object for the diagonal
#'   elements of the precision matrix, created by one of:
#'   \itemize{
#'     \item \code{\link{gamma_prior}()}: Gamma(shape, rate) prior (default).
#'     \item \code{\link{exponential_prior}()}: Exponential(rate) prior.
#'   }
#'   Only used for models with continuous variables (GGM and mixed MRF).
#'   Ignored for pure ordinal models.
#'   Default: \code{gamma_prior(shape = 1, rate = 1)}.
#'
#' @param delta Non-negative numeric, or \code{NULL} for the dimension-
#'   adaptive default. Determinant-tilt exponent on the continuous-block
#'   precision matrix \eqn{K} (GGM) or \eqn{K_{yy}} (mixed MRF):
#'   multiplies the prior by \eqn{|K|^{\delta}}, softly repelling the
#'   chain from the positive-definite cone boundary. \code{delta = NULL}
#'   (default) auto-resolves to \eqn{0.5 \log(p)} where \eqn{p} is the
#'   dimension of the continuous precision matrix (the number of
#'   variables for GGM, the number of continuous variables for mixed
#'   MRF). The rule is the simple form of the dimension-adaptive scaling
#'   \eqn{\delta(p) = c \log p} with \eqn{c \in (0.3, 0.6)} discussed in
#'   the companion paper on determinant-tilted spike-and-slab priors
#'   (Marsman et al., in preparation). Pass \code{delta = 0} for the
#'   untilted prior (the companion-paper baseline) or a non-negative
#'   numeric to override. Both NUTS and adaptive-Metropolis update paths
#'   apply the tilt. Not allowed for pure ordinal models (no precision
#'   matrix to tilt).
#'
#' @param pairwise_scale `r lifecycle::badge("deprecated")` Double.
#'   Scale of the Cauchy prior for pairwise
#'   interaction parameters. Use \code{interaction_prior} instead.
#'   Default: \code{1}.
#'
#' @param standardize Logical. If \code{TRUE}, the prior scale for each
#'   pairwise interaction is adjusted based on the range of response scores.
#'   Variables with more response categories have larger score products
#'   \eqn{x_i \cdot x_j}, which typically correspond to smaller interaction
#'   effects \eqn{\sigma_{ij}}. Without standardization, a fixed prior scale
#'   is relatively wide for these smaller effects, resulting in less shrinkage
#'   for high-category pairs and more shrinkage for low-category pairs.
#'   Standardization scales the prior proportionally to the maximum score
#'   product, ensuring equivalent relative shrinkage across all pairs.
#'   After internal recoding, regular ordinal variables have scores
#'   \eqn{0, 1, \ldots, m}. The adjusted scale for the interaction between
#'   variables \eqn{i} and \eqn{j} is \code{pairwise_scale * m_i * m_j},
#'   so that \code{pairwise_scale} itself applies to the unit interval case
#'   (binary variables where \eqn{m_i = m_j = 1}). For Blume-Capel variables
#'   with reference category \eqn{b}, scores are centered as
#'   \eqn{-b, \ldots, m-b}, and the adjustment uses the maximum absolute
#'   product of the score endpoints. For mixed pairs, ordinal variables use
#'   raw score endpoints \eqn{(0, m)} and Blume-Capel variables use centered
#'   score endpoints \eqn{(-b, m-b)}.
#'   Default: \code{FALSE}.
#'
#' @param main_alpha,main_beta `r lifecycle::badge("deprecated")` Double.
#'   Shape parameters of the beta-prime prior for threshold parameters.
#'   Use \code{threshold_prior} instead. Defaults: \code{main_alpha = 0.5}
#'   and \code{main_beta = 0.5}.
#'
#' @param edge_selection Logical. Whether to perform Bayesian edge selection.
#'   If \code{FALSE}, the model estimates all edges. Default: \code{TRUE}.
#'
#' @param edge_prior An edge prior specification object, or a character string
#'   (deprecated). Specifies the prior for edge inclusion.
#'   Preferred: pass an object from one of:
#'   \itemize{
#'     \item \code{\link{bernoulli_prior}()}: Fixed inclusion probability (default).
#'     \item \code{\link{beta_bernoulli_prior}()}: Beta-distributed inclusion.
#'     \item \code{\link{sbm_prior}()}: Stochastic Block Model.
#'   }
#'   Legacy character strings \code{"Bernoulli"}, \code{"Beta-Bernoulli"},
#'   \code{"Stochastic-Block"} are still accepted but deprecated.
#'   Default: \code{bernoulli_prior(0.5)}.
#'
#' @param inclusion_probability `r lifecycle::badge("deprecated")` Numeric
#'   scalar. Use \code{edge_prior = bernoulli_prior(inclusion_probability)}
#'   instead. Default: \code{0.5}.
#'
#' @param beta_bernoulli_alpha,beta_bernoulli_beta
#'   `r lifecycle::badge("deprecated")` Double. Use
#'   \code{edge_prior = beta_bernoulli_prior(alpha, beta)} instead.
#'   Defaults: \code{1}.
#'
#' @param beta_bernoulli_alpha_between,beta_bernoulli_beta_between
#'   `r lifecycle::badge("deprecated")` Double. Use
#'   \code{edge_prior = sbm_prior(alpha_between, beta_between)} instead.
#'   Defaults: \code{1}.
#'
#' @param dirichlet_alpha `r lifecycle::badge("deprecated")` Double. Use
#'   \code{edge_prior = sbm_prior(dirichlet_alpha = ...)} instead.
#'   Default: \code{1}.
#'
#' @param lambda `r lifecycle::badge("deprecated")` Double. Use
#'   \code{edge_prior = sbm_prior(lambda = ...)} instead.
#'   Default: \code{1}.
#'
#' @param na_action Character. Specifies missing data handling. Either
#'   \code{"listwise"} (drop rows with missing values) or \code{"impute"}
#'   (perform single imputation during sampling). Default: \code{"listwise"}.
#'
#' @param display_progress Character. Controls progress reporting during
#'   sampling. Options: \code{"per-chain"} (separate bar per chain),
#'   \code{"total"} (single combined bar), or \code{"none"} (no progress).
#'   Default: \code{"per-chain"}.
#'
#' @param progress_callback An optional R function with signature
#'   \code{function(completed, total)} that is called at regular intervals
#'   during sampling, where \code{completed} is the number of iterations
#'   completed across all chains and \code{total} is the total number of
#'   iterations. Useful for external front-ends (e.g., JASP) that supply
#'   their own progress reporting.
#'   When \code{NULL} (the default), no callback is invoked.
#'
#' @param verbose Logical. If \code{TRUE}, prints informational messages
#'   during data processing (e.g., missing data handling, variable recoding).
#'   Defaults to \code{getOption("bgms.verbose", TRUE)}. Set
#'   \code{options(bgms.verbose = FALSE)} to suppress messages globally.
#'
#' @param update_method Character. Specifies how the MCMC sampler updates
#'   the model parameters:
#'   \describe{
#'     \item{"adaptive-metropolis"}{Componentwise adaptive Metropolis--Hastings
#'       with Robbins--Monro proposal adaptation.}
#'     \item{"nuts"}{The No-U-Turn Sampler, a gradient-based sampler available
#'       for all variable types, including under edge selection. The Gaussian
#'       graphical model uses a free-element Cholesky parameterization that keeps
#'       the precision matrix positive-definite; the mixed model uses RATTLE
#'       constrained integration when excluded edges impose constraints.}
#'     \item{"gibbs"}{An exact Gibbs sampler for the Gaussian graphical model.
#'       Available only for all-continuous data on a fixed graph
#'       (\code{edge_selection = FALSE}) with a Normal (slab) interaction prior,
#'       a Gamma(shape = 1) scale prior on the precision diagonal, and
#'       \code{delta = 0}.}
#'   }
#'   Default: \code{"nuts"}.
#'
#' @param target_accept Numeric between 0 and 1. Target acceptance rate for
#'   the sampler. Defaults are set automatically if not supplied:
#'   \code{0.44} for adaptive Metropolis and \code{0.80} for NUTS.
#'
#' @param nuts_max_depth Integer. Maximum tree depth in NUTS. Must be positive.
#'   Default: \code{10}.
#'
#' @param learn_mass_matrix Logical. If \code{TRUE}, adapt a diagonal mass
#'   matrix during warmup (NUTS only). If \code{FALSE}, use the identity
#'   matrix. Default: \code{TRUE}.
#'
#' @param chains Integer. Number of parallel chains to run. Default: \code{4}.
#'
#' @param cores Integer. Number of CPU cores for parallel execution.
#'   Default: \code{parallel::detectCores()}.
#'
#' @param seed Optional integer. Random seed for reproducibility. Must be a
#'   single non-negative integer.
#'
#' @param interaction_scale,burnin,save,threshold_alpha,threshold_beta
#'   `r lifecycle::badge("deprecated")`
#'   Deprecated arguments as of \strong{bgms 0.1.6.0}.
#'   Use `interaction_prior`, `warmup`, and `threshold_prior` instead.
#'
#' @return
#' An S7 object of class \code{bgms} with posterior summaries, posterior mean
#' matrices, and access to raw MCMC draws. Its fields are accessible with
#' \code{$} and \code{[[} for backward compatibility, and it can be passed to
#' \code{print()}, \code{summary()}, and \code{coef()}.
#'
#' Main components include:
#' \itemize{
#'   \item \code{posterior_summary_main}: Data frame with posterior summaries
#'     (mean, sd, MCSE, ESS, Rhat) for main-effect parameters.
#'     For OMRF models these are category thresholds;
#'     for mixed MRF models these are discrete thresholds and
#'     continuous means. \code{NULL} for GGM models (no main effects).
#'   \item \code{posterior_summary_quadratic}: Data frame with posterior
#'     summaries for the residual variance parameters (GGM and mixed MRF).
#'     \code{NULL} for OMRF models.
#'   \item \code{posterior_summary_pairwise}: Data frame with posterior
#'     summaries for partial association parameters.
#'   \item \code{posterior_summary_indicator}: Data frame with posterior
#'     summaries for edge inclusion indicators (if \code{edge_selection = TRUE}).
#'
#'   \item \code{posterior_mean_main}: Posterior mean of main-effect
#'     parameters. \code{NULL} for GGM models. For OMRF: a matrix
#'     (p x max_categories) of category thresholds. For mixed MRF: a list
#'     with \code{$discrete} (threshold matrix) and \code{$continuous}
#'     (q x 1 matrix of means).
#'   \item \code{posterior_mean_pairwise}: Symmetric matrix of posterior
#'     mean partial associations (zero diagonal). For continuous variables
#'     these are unstandardized partial correlations; for discrete variables
#'     these are half the log adjacent-category odds ratio. Use
#'     [extract_precision()], [extract_partial_correlations()], or
#'     [extract_log_odds()] to convert to interpretable scales.
#'   \item \code{posterior_mean_residual_variance}: Named numeric vector of
#'     posterior mean residual variances \eqn{1/\Theta_{ii}}{1/Theta_ii}.
#'     Present for GGM and mixed MRF models; \code{NULL} for OMRF.
#'   \item \code{posterior_mean_indicator}: Symmetric matrix of posterior mean
#'     inclusion probabilities (if edge selection was enabled).
#'
#'   \item  Additional summaries returned when
#'     \code{edge_prior = "Stochastic-Block"}. For more details about this prior
#'     see \insertCite{SekulovskiEtAl_2025;textual}{bgms}.
#'    \itemize{
#'       \item \code{posterior_summary_pairwise_allocations}: Data frame with
#'       posterior summaries (mean, sd, MCSE, ESS, Rhat) for the pairwise
#'       cluster co-occurrence of the nodes. This serves to indicate
#'       whether the estimated posterior allocations,co-clustering matrix
#'       and posterior cluster probabilities (see blow) have converged.
#'       \item \code{posterior_mean_coclustering_matrix}: a symmetric matrix of
#'       pairwise proportions of occurrence of every variable. This matrix
#'       can be plotted to visually inspect the estimated number of clusters
#'       and visually inspect nodes that tend to switch clusters.
#'       \item \code{posterior_mean_allocations}: A vector with the posterior mean
#'       of the cluster allocations of the nodes. This is calculated using the method
#'       proposed in \insertCite{Dahl2009;textual}{bgms}.
#'       \item \code{posterior_mode_allocations}: A vector with the posterior
#'        mode of the cluster allocations of the nodes.
#'       \item \code{posterior_num_blocks}: A data frame with the estimated
#'       posterior inclusion probabilities for all the possible number of clusters.
#'       }
#'   \item \code{raw_samples}: A list of raw MCMC draws per chain:
#'     \describe{
#'       \item{\code{main}}{List of main effect samples.}
#'       \item{\code{pairwise}}{List of pairwise effect samples.}
#'       \item{\code{indicator}}{List of indicator samples
#'         (if edge selection enabled).}
#'       \item{\code{allocations}}{List of cluster allocations
#'         (if SBM prior used).}
#'       \item{\code{nchains}}{Number of chains.}
#'       \item{\code{niter}}{Number of post--warmup iterations per chain.}
#'       \item{\code{parameter_names}}{Named lists of parameter labels.}
#'     }
#'
#'   \item \code{arguments}: A list of function call arguments and metadata
#'     (e.g., number of variables, warmup, sampler settings, package version).
#' }
#'
#' The \code{summary()} method prints formatted posterior summaries, and
#' \code{coef()} extracts posterior mean matrices.
#'
#' NUTS diagnostics (tree depth, divergences, energy, E-BFMI) are included
#' in \code{fit$nuts_diag} if \code{update_method = "nuts"}.
#'
#' @references
#'   \insertAllCited{}
#'
#' @examples
#' \donttest{
#' # Run bgm on subset of the Wenchuan dataset
#' fit = bgm(x = Wenchuan[, 1:5], chains = 2)
#'
#' # Posterior inclusion probabilities
#' summary(fit)$indicator
#'
#' # Posterior pairwise effects
#' summary(fit)$pairwise
#' }
#'
#' @export
bgm = function(
  x,
  variable_type = "ordinal",
  baseline_category,
  iter = 2e3,
  warmup = 2e3,
  interaction_prior = cauchy_prior(scale = 1),
  threshold_prior = beta_prime_prior(alpha = 0.5, beta = 0.5),
  means_prior = normal_prior(scale = 1),
  precision_scale_prior = gamma_prior(shape = 1, rate = 1),
  delta = NULL,
  edge_selection = TRUE,
  edge_prior = bernoulli_prior(0.5),
  na_action = c("listwise", "impute"),
  update_method = c("nuts", "adaptive-metropolis", "gibbs"),
  target_accept,
  nuts_max_depth = 10,
  learn_mass_matrix = TRUE,
  chains = 4,
  cores = parallel::detectCores(),
  display_progress = c("per-chain", "total", "none"),
  seed = NULL,
  standardize = FALSE,
  verbose = getOption("bgms.verbose", TRUE),
  progress_callback = NULL,
  # Deprecated prior arguments (v0.1.6.0 and earlier)
  pairwise_scale,
  main_alpha,
  main_beta,
  inclusion_probability,
  beta_bernoulli_alpha,
  beta_bernoulli_beta,
  beta_bernoulli_alpha_between,
  beta_bernoulli_beta_between,
  dirichlet_alpha,
  lambda,
  # Deprecated arguments (v0.1.6.0)
  interaction_scale,
  burnin,
  save,
  threshold_alpha,
  threshold_beta
) {
  # Set verbose option for internal functions, restore on exit

  old_verbose = getOption("bgms.verbose")
  options(bgms.verbose = verbose)
  on.exit(options(bgms.verbose = old_verbose), add = TRUE)

  # --- Legacy deprecation: v0.1.6.0 renames -----------------------------------
  if(hasArg(interaction_scale)) {
    lifecycle::deprecate_warn(
      "0.1.6.0", "bgm(interaction_scale =)",
      "bgm(interaction_prior =)"
    )
    if(!hasArg(pairwise_scale) &&
      identical(interaction_prior, cauchy_prior(scale = 1))) {
      interaction_prior = cauchy_prior(scale = interaction_scale)
    }
  }

  if(hasArg(burnin)) {
    lifecycle::deprecate_warn("0.1.6.0", "bgm(burnin =)", "bgm(warmup =)")
    if(!hasArg(warmup)) warmup = burnin
  }

  if(hasArg(save)) {
    lifecycle::deprecate_warn("0.1.6.0", "bgm(save =)")
  }

  if(hasArg(threshold_alpha) || hasArg(threshold_beta)) {
    lifecycle::deprecate_warn(
      "0.1.6.0",
      "bgm(threshold_alpha =, threshold_beta =)",
      "bgm(threshold_prior =)"
    )
    if(identical(threshold_prior, beta_prime_prior(0.5, 0.5))) {
      ta = if(hasArg(threshold_alpha)) threshold_alpha else 0.5
      tb = if(hasArg(threshold_beta)) threshold_beta else 0.5
      threshold_prior = beta_prime_prior(alpha = ta, beta = tb)
    }
  }

  # --- Legacy deprecation: scalar prior parameters ----------------------------
  if(hasArg(pairwise_scale)) {
    lifecycle::deprecate_warn(
      "0.2.0", "bgm(pairwise_scale =)",
      "bgm(interaction_prior =)"
    )
    if(identical(interaction_prior, cauchy_prior(scale = 1))) {
      interaction_prior = cauchy_prior(scale = pairwise_scale)
    }
  }

  if(hasArg(main_alpha) || hasArg(main_beta)) {
    lifecycle::deprecate_warn(
      "0.2.0",
      "bgm(main_alpha =)",
      "bgm(threshold_prior =)"
    )
    if(identical(threshold_prior, beta_prime_prior(0.5, 0.5))) {
      ma = if(hasArg(main_alpha)) main_alpha else 0.5
      mb = if(hasArg(main_beta)) main_beta else 0.5
      threshold_prior = beta_prime_prior(alpha = ma, beta = mb)
    }
  }

  # Handle edge_prior: accept both string (deprecated) and object (new)
  if(is.character(edge_prior)) {
    lifecycle::deprecate_warn(
      "0.2.0", "bgm(edge_prior = 'must be a prior object')",
      "bgm(edge_prior = 'bernoulli_prior()')"
    )
    edge_prior_str = match.arg(edge_prior,
      choices = c("Bernoulli", "Beta-Bernoulli", "Stochastic-Block")
    )

    ip = if(hasArg(inclusion_probability)) inclusion_probability else 0.5
    bba = if(hasArg(beta_bernoulli_alpha)) beta_bernoulli_alpha else 1
    bbb = if(hasArg(beta_bernoulli_beta)) beta_bernoulli_beta else 1
    bbab = if(hasArg(beta_bernoulli_alpha_between)) beta_bernoulli_alpha_between else 1
    bbbb = if(hasArg(beta_bernoulli_beta_between)) beta_bernoulli_beta_between else 1
    da = if(hasArg(dirichlet_alpha)) dirichlet_alpha else 1
    lam = if(hasArg(lambda)) lambda else 1

    edge_prior = switch(edge_prior_str,
      "Bernoulli" = bernoulli_prior(inclusion_probability = ip),
      "Beta-Bernoulli" = beta_bernoulli_prior(alpha = bba, beta = bbb),
      "Stochastic-Block" = sbm_prior(
        alpha = bba, beta = bbb,
        alpha_between = bbab, beta_between = bbbb,
        dirichlet_alpha = da, lambda = lam
      )
    )
  } else {
    # Warn if loose edge params are also supplied alongside an object
    if(hasArg(inclusion_probability)) {
      lifecycle::deprecate_warn(
        "0.2.0", "bgm(inclusion_probability =)",
        "bgm(edge_prior = 'bernoulli_prior()')"
      )
    }
    if(hasArg(beta_bernoulli_alpha)) {
      lifecycle::deprecate_warn(
        "0.2.0", "bgm(beta_bernoulli_alpha =)",
        "bgm(edge_prior = 'beta_bernoulli_prior()')"
      )
    }
  }

  # --- Unpack prior objects to flat parameters for bgm_spec -------------------
  ip = unpack_interaction_prior(interaction_prior)
  tp = unpack_threshold_prior(threshold_prior)
  mp = unpack_parameter_prior(means_prior)
  sp = unpack_scale_prior(precision_scale_prior)

  # --- Build spec, sample, build output ----------------------------------------
  spec = bgm_spec(
    x = x,
    model_type = "omrf",
    variable_type = variable_type,
    baseline_category = if(hasArg(baseline_category)) baseline_category else 0L,
    na_action = na_action,
    interaction_prior_type = ip$interaction_prior_type,
    pairwise_scale = ip$pairwise_scale,
    interaction_alpha = ip$interaction_alpha,
    interaction_beta = ip$interaction_beta,
    threshold_prior_type = tp$threshold_prior_type,
    main_alpha = tp$main_alpha,
    main_beta = tp$main_beta,
    threshold_scale = tp$threshold_scale,
    means_prior_type = mp$prior_type,
    means_scale = mp$scale,
    means_alpha = mp$alpha,
    means_beta = mp$beta,
    scale_prior_type = sp$scale_prior_type,
    scale_shape = sp$scale_shape,
    scale_rate = sp$scale_rate,
    delta = delta,
    standardize = standardize,
    edge_selection = edge_selection,
    edge_prior = edge_prior,
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

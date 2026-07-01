# ==============================================================================
# bgm_spec: validated intermediate specification object
# ==============================================================================
#
# Central construction point for all bgm/bgmCompare models. Three layers:
#   bgm_spec()          --- user-facing: validates inputs, assembles sub-lists
#   new_bgm_spec()      --- low-level: type/presence assertions per field
#   validate_bgm_spec() --- cross-field invariant checks
#
# The result is an S3 list of class "bgm_spec" consumed by run_sampler()
# and build_output().
#
# The per-family spec builders (build_spec_*) live in build_spec.R, and the
# spec -> $arguments converters (build_arguments*) in build_arguments.R (S4 split).
# ==============================================================================


# ==============================================================================
# new_bgm_spec()  --- low-level constructor
# ==============================================================================
#
# Asserts presence and type of every field. Does NOT validate values
# (that's done upstream by the individual validators) or cross-field
# invariants (that's validate_bgm_spec).
# ==============================================================================
new_bgm_spec = function(model_type, data, variables, missing, prior,
                        sampler, precomputed = list()) {
  # --- top-level structure ---
  stopifnot(
    is.character(model_type), length(model_type) == 1L,
    model_type %in% c("ggm", "omrf", "compare", "mixed_mrf")
  )

  # --- data sub-list ---
  stopifnot(is.list(data))
  if(model_type == "mixed_mrf") {
    stopifnot(is.matrix(data$x_discrete))
    stopifnot(is.matrix(data$x_continuous))
  } else {
    stopifnot(is.matrix(data$x))
  }
  stopifnot(is.character(data$data_columnnames))
  stopifnot(is.integer(data$num_variables), length(data$num_variables) == 1L)
  stopifnot(is.integer(data$num_cases), length(data$num_cases) == 1L)

  if(model_type == "omrf" || model_type == "compare") {
    stopifnot(
      is.integer(data$num_categories),
      length(data$num_categories) == data$num_variables
    )
  }
  if(model_type == "mixed_mrf") {
    stopifnot(
      is.integer(data$num_categories),
      length(data$num_categories) == data$num_discrete
    )
  }


  if(model_type == "compare") {
    stopifnot(is.integer(data$group), length(data$group) == data$num_cases)
    stopifnot(is.integer(data$num_groups), length(data$num_groups) == 1L)
    stopifnot(is.matrix(data$group_indices))
    stopifnot(is.matrix(data$projection))
  }

  # --- variables sub-list ---
  stopifnot(is.list(variables))
  stopifnot(is.character(variables$variable_type))
  stopifnot(is.logical(variables$is_ordinal))
  stopifnot(is.logical(variables$is_continuous), length(variables$is_continuous) == 1L)
  stopifnot(is.integer(variables$baseline_category))

  # --- missing sub-list ---
  stopifnot(is.list(missing))
  stopifnot(
    is.character(missing$na_action), length(missing$na_action) == 1L,
    missing$na_action %in% c("listwise", "impute")
  )
  stopifnot(is.logical(missing$na_impute), length(missing$na_impute) == 1L)
  # missing_index can be NULL (no missing) or a matrix
  if(!is.null(missing$missing_index)) {
    stopifnot(is.matrix(missing$missing_index))
  }
  # mixed MRF uses separate indices for discrete and continuous
  if(!is.null(missing$missing_index_discrete)) {
    stopifnot(is.matrix(missing$missing_index_discrete))
  }
  if(!is.null(missing$missing_index_continuous)) {
    stopifnot(is.matrix(missing$missing_index_continuous))
  }

  # --- prior sub-list ---
  stopifnot(is.list(prior))
  # All model types now carry interaction_prior_type
  stopifnot(
    is.character(prior$interaction_prior_type),
    length(prior$interaction_prior_type) == 1L
  )
  stopifnot(is.numeric(prior$pairwise_scale), length(prior$pairwise_scale) == 1L)
  if(model_type %in% c("omrf", "compare", "mixed_mrf")) {
    stopifnot(
      is.character(prior$threshold_prior_type),
      length(prior$threshold_prior_type) == 1L
    )
    stopifnot(is.logical(prior$standardize), length(prior$standardize) == 1L)
  }
  if(model_type %in% c("omrf", "compare")) {
    stopifnot(is.matrix(prior$pairwise_scaling_factors))
  }
  if(model_type == "mixed_mrf") {
    stopifnot(is.logical(prior$standardize), length(prior$standardize) == 1L)
  }
  if(model_type %in% c("ggm", "omrf", "mixed_mrf")) {
    stopifnot(is.logical(prior$edge_selection), length(prior$edge_selection) == 1L)
    stopifnot(is.character(prior$edge_prior), length(prior$edge_prior) == 1L)
    stopifnot(is.matrix(prior$inclusion_probability))
  }
  if(model_type == "compare") {
    stopifnot(
      is.logical(prior$difference_selection),
      length(prior$difference_selection) == 1L
    )
    stopifnot(
      is.logical(prior$main_difference_selection),
      length(prior$main_difference_selection) == 1L
    )
    stopifnot(
      is.character(prior$difference_prior),
      length(prior$difference_prior) == 1L
    )
    stopifnot(
      is.numeric(prior$difference_scale),
      length(prior$difference_scale) == 1L
    )
    stopifnot(is.matrix(prior$inclusion_probability_difference))
  }

  # --- sampler sub-list ---
  stopifnot(is.list(sampler))
  stopifnot(is.character(sampler$update_method), length(sampler$update_method) == 1L)
  stopifnot(is.numeric(sampler$target_accept), length(sampler$target_accept) == 1L)
  stopifnot(is.integer(sampler$iter), length(sampler$iter) == 1L)
  stopifnot(is.integer(sampler$warmup), length(sampler$warmup) == 1L)
  stopifnot(is.integer(sampler$chains), length(sampler$chains) == 1L)
  stopifnot(is.integer(sampler$cores), length(sampler$cores) == 1L)
  stopifnot(is.integer(sampler$nuts_max_depth), length(sampler$nuts_max_depth) == 1L)
  stopifnot(is.logical(sampler$learn_mass_matrix), length(sampler$learn_mass_matrix) == 1L)
  stopifnot(is.integer(sampler$seed), length(sampler$seed) == 1L)
  stopifnot(is.integer(sampler$progress_type), length(sampler$progress_type) == 1L)

  # --- precomputed sub-list ---
  stopifnot(is.list(precomputed))

  structure(
    list(
      model_type  = model_type,
      data        = data,
      variables   = variables,
      missing     = missing,
      prior       = prior,
      sampler     = sampler,
      precomputed = precomputed
    ),
    class = "bgm_spec"
  )
}


# ==============================================================================
# validate_bgm_spec()  --- cross-field invariant checks
# ==============================================================================
validate_bgm_spec = function(spec) {
  mt = spec$model_type

  # GGM invariants
  if(mt == "ggm") {
    if(!isTRUE(spec$variables$is_continuous)) {
      stop("bgm_spec: model_type = 'ggm' requires is_continuous = TRUE.")
    }
  }

  # Compare invariants
  if(mt == "compare") {
    if(is.null(spec$data$group)) {
      stop("bgm_spec: model_type = 'compare' requires data$group.")
    }
    if(spec$data$num_groups < 2L) {
      stop("bgm_spec: model_type = 'compare' requires num_groups >= 2.")
    }
  }

  # Edge selection consistency
  if(mt %in% c("ggm", "omrf", "mixed_mrf")) {
    if(spec$prior$edge_selection && spec$prior$edge_prior == "Not Applicable") {
      stop("bgm_spec: edge_selection = TRUE but edge_prior = 'Not Applicable'.")
    }
  }

  # Scaling factors dimensions
  if(mt %in% c("omrf", "compare")) {
    nv = spec$data$num_variables
    sf = spec$prior$pairwise_scaling_factors
    if(nrow(sf) != nv || ncol(sf) != nv) {
      stop(
        "bgm_spec: pairwise_scaling_factors dimensions (",
        nrow(sf), "x", ncol(sf), ") don't match num_variables (", nv, ")."
      )
    }
  }

  # num_categories length (OMRF / compare)
  if(mt == "omrf" || mt == "compare") {
    if(length(spec$data$num_categories) != spec$data$num_variables) {
      stop("bgm_spec: num_categories length doesn't match num_variables.")
    }
  }
  if(mt == "mixed_mrf") {
    if(length(spec$data$num_categories) != spec$data$num_discrete) {
      stop("bgm_spec: num_categories length doesn't match num_discrete.")
    }
    allowed = c("adaptive-metropolis", "nuts")
    if(!(spec$sampler$update_method %in% allowed)) {
      stop(
        "bgm_spec: model_type = 'mixed_mrf' requires update_method in ",
        paste(sQuote(allowed), collapse = " or "), ". Got '",
        spec$sampler$update_method, "'."
      )
    }
  }

  invisible(spec)
}


# ==============================================================================
# bgm_spec()  --- user-facing constructor
# ==============================================================================
#
# Validates all user inputs via dedicated validators, assembles sub-lists,
# and passes through new_bgm_spec() and validate_bgm_spec().
#
# Parameters mirror the union of bgm() and bgmCompare() arguments.
# ==============================================================================
bgm_spec = function(x,
                    model_type = c("omrf", "ggm", "compare", "mixed_mrf"),
                    # Variable specification
                    variable_type = "ordinal",
                    baseline_category = 0L,
                    # Data (compare-specific)
                    y = NULL,
                    group_indicator = NULL,
                    # Missing data
                    na_action = c("listwise", "impute"),
                    # Priors (new: prior objects unpacked by bgm())
                    interaction_prior_type = "cauchy",
                    pairwise_scale = 1,
                    interaction_alpha = NA_real_,
                    interaction_beta = NA_real_,
                    threshold_prior_type = "beta-prime",
                    main_alpha = 0.5,
                    main_beta = 0.5,
                    threshold_scale = NA_real_,
                    means_prior_type = "normal",
                    means_scale = 1,
                    means_alpha = NA_real_,
                    means_beta = NA_real_,
                    scale_prior_type = "gamma",
                    scale_shape = 1,
                    scale_rate = 1,
                    delta = NULL,
                    standardize = FALSE,
                    edge_selection = TRUE,
                    edge_prior = bernoulli_prior(0.5),
                    # Legacy edge prior params (accepted for backward compat)
                    inclusion_probability = 0.5,
                    beta_bernoulli_alpha_between = 1,
                    beta_bernoulli_beta_between = 1,
                    dirichlet_alpha = 1,
                    lambda = 1,
                    # Priors (compare-specific)
                    difference_selection = TRUE,
                    main_difference_selection = FALSE,
                    difference_prior = c(
                      "Bernoulli", "Beta-Bernoulli", "Stochastic-Block"
                    ),
                    difference_scale = 1,
                    difference_probability = 0.5,
                    # Compare difference prior hyperparameters
                    beta_bernoulli_alpha = 1,
                    beta_bernoulli_beta = 1,
                    difference_beta_bernoulli_alpha_between = 1,
                    difference_beta_bernoulli_beta_between = 1,
                    difference_dirichlet_alpha = 1,
                    difference_lambda = 1,
                    # Sampler
                    update_method = c(
                      "nuts",
                      "adaptive-metropolis"
                    ),
                    target_accept = NULL,
                    iter = 10000L,
                    warmup = 1000L,
                    nuts_max_depth = 10L,
                    learn_mass_matrix = TRUE,
                    chains = 4L,
                    cores = parallel::detectCores(),
                    seed = NULL,
                    display_progress = c("per-chain", "total", "none"),
                    verbose = TRUE,
                    progress_callback = NULL) {
  model_type = match.arg(model_type)
  na_action = tryCatch(match.arg(na_action), error = function(e) {
    stop(paste0(
      "The na_action argument should be one of \"listwise\" or \"impute\", not \"",
      na_action, "\"."
    ), call. = FALSE)
  })

  # --- Data validation --------------------------------------------------------
  x = data_check(x, "x")
  data_columnnames = if(is.null(colnames(x))) {
    paste0("Variable ", seq_len(ncol(x)))
  } else {
    colnames(x)
  }
  num_variables = ncol(x)

  # --- Variable types ---------------------------------------------------------
  allow_continuous = (model_type != "compare")
  vt = validate_variable_types(
    variable_type    = variable_type,
    num_variables    = num_variables,
    allow_continuous = allow_continuous,
    allow_mixed      = (model_type != "compare"),
    caller           = if(model_type == "compare") "bgmCompare" else "bgm"
  )
  variable_type = vt$variable_type
  is_ordinal = vt$variable_bool
  is_continuous = vt$is_continuous
  is_mixed = vt$is_mixed

  # Resolve model_type if "omrf" default was kept but data is continuous
  if(model_type == "omrf" && is_continuous) {
    model_type = "ggm"
  }
  if(model_type == "omrf" && is_mixed) {
    model_type = "mixed_mrf"
  }

  # Auto-resolve delta = NULL to the dimension-adaptive default 0.5 * log(p),
  # where p is the dimension of the continuous precision matrix. For models
  # without a continuous block (omrf, compare) the tilt has no target, so
  # NULL resolves to 0.
  if(is.null(delta)) {
    delta = if(model_type == "ggm") {
      0.5 * log(max(num_variables, 1))
    } else if(model_type == "mixed_mrf") {
      0.5 * log(max(sum(!is_ordinal), 1))
    } else {
      0
    }
  }

  # Validate determinant-tilt exponent and reject for pure-ordinal models
  if(!is.numeric(delta) || length(delta) != 1L || is.na(delta) ||
    !is.finite(delta) || delta < 0) {
    stop("'delta' must be a single finite non-negative numeric, or NULL.")
  }
  if(delta > 0 && model_type %in% c("omrf", "compare")) {
    stop(
      "'delta' (determinant tilt) requires continuous variables; the ",
      "current model_type is '", model_type, "', which has no precision ",
      "matrix to tilt. Pass delta = 0 or use continuous data."
    )
  }

  # --- Sampler (needs is_continuous and edge_selection early) ------------------
  sampler = validate_sampler(
    update_method = update_method,
    target_accept = target_accept,
    iter = iter,
    warmup = warmup,
    nuts_max_depth = nuts_max_depth,
    learn_mass_matrix = learn_mass_matrix,
    chains = chains,
    cores = cores,
    seed = seed,
    display_progress = display_progress,
    is_continuous = is_continuous,
    edge_selection = if(model_type == "compare") FALSE else edge_selection,
    verbose = verbose,
    progress_callback = progress_callback
  )

  # --- Resolve edge prior object -----------------------------------------------
  if(inherits(edge_prior, "bgms_indicator_prior")) {
    ep_flat = unpack_indicator_prior(edge_prior, num_variables)
  } else if(is.character(edge_prior)) {
    # Legacy string path (tests and bgmCompare may call bgm_spec directly)
    edge_prior_str = match.arg(edge_prior,
      choices = c("Bernoulli", "Beta-Bernoulli", "Stochastic-Block")
    )
    ep_flat = validate_edge_prior(
      edge_selection = edge_selection, edge_prior = edge_prior_str,
      inclusion_probability = inclusion_probability,
      num_variables = num_variables,
      beta_bernoulli_alpha = beta_bernoulli_alpha,
      beta_bernoulli_beta = beta_bernoulli_beta,
      beta_bernoulli_alpha_between = beta_bernoulli_alpha_between,
      beta_bernoulli_beta_between = beta_bernoulli_beta_between,
      dirichlet_alpha = dirichlet_alpha, lambda = lambda
    )
    ep_flat$beta_bernoulli_alpha = beta_bernoulli_alpha
    ep_flat$beta_bernoulli_beta = beta_bernoulli_beta
    ep_flat$beta_bernoulli_alpha_between = beta_bernoulli_alpha_between
    ep_flat$beta_bernoulli_beta_between = beta_bernoulli_beta_between
    ep_flat$dirichlet_alpha = dirichlet_alpha
    ep_flat$lambda = lambda
  } else {
    stop(
      "'edge_prior' must be a bgms_indicator_prior object.",
      " Use bernoulli_prior(), beta_bernoulli_prior(), or sbm_prior()."
    )
  }
  # Override edge_selection if explicitly FALSE
  if(!edge_selection) {
    ep_flat$edge_selection = FALSE
    ep_flat$edge_prior = "Not Applicable"
    ep_flat$inclusion_probability = matrix(0.5, nrow = 1, ncol = 1)
  }

  # --- Build by model type ----------------------------------------------------
  if(model_type == "ggm") {
    spec = build_spec_ggm(
      x = x, data_columnnames = data_columnnames,
      num_variables = num_variables,
      variable_type = variable_type, is_ordinal = is_ordinal,
      is_continuous = is_continuous,
      baseline_category = as.integer(rep(0L, num_variables)),
      na_action = na_action, sampler = sampler,
      interaction_prior_type = interaction_prior_type,
      pairwise_scale = pairwise_scale,
      interaction_alpha = interaction_alpha,
      interaction_beta = interaction_beta,
      scale_prior_type = scale_prior_type,
      scale_shape = scale_shape,
      scale_rate = scale_rate,
      delta = delta,
      edge_prior_flat = ep_flat
    )
  } else if(model_type == "mixed_mrf") {
    spec = build_spec_mixed_mrf(
      x = x, data_columnnames = data_columnnames,
      num_variables = num_variables,
      variable_type = variable_type, is_ordinal = is_ordinal,
      baseline_category = baseline_category,
      na_action = na_action, sampler = sampler,
      interaction_prior_type = interaction_prior_type,
      pairwise_scale = pairwise_scale,
      interaction_alpha = interaction_alpha,
      interaction_beta = interaction_beta,
      threshold_prior_type = threshold_prior_type,
      main_alpha = main_alpha, main_beta = main_beta,
      threshold_scale = threshold_scale,
      means_prior_type = means_prior_type,
      means_scale = means_scale,
      means_alpha = means_alpha,
      means_beta = means_beta,
      scale_prior_type = scale_prior_type,
      scale_shape = scale_shape,
      scale_rate = scale_rate,
      delta = delta,
      standardize = standardize,
      edge_prior_flat = ep_flat
    )
  } else if(model_type == "omrf") {
    spec = build_spec_omrf(
      x = x, data_columnnames = data_columnnames,
      num_variables = num_variables,
      variable_type = variable_type, is_ordinal = is_ordinal,
      is_continuous = is_continuous,
      baseline_category = baseline_category,
      na_action = na_action, sampler = sampler,
      interaction_prior_type = interaction_prior_type,
      pairwise_scale = pairwise_scale,
      interaction_alpha = interaction_alpha,
      interaction_beta = interaction_beta,
      threshold_prior_type = threshold_prior_type,
      main_alpha = main_alpha, main_beta = main_beta,
      threshold_scale = threshold_scale,
      standardize = standardize,
      edge_prior_flat = ep_flat
    )
  } else {
    spec = build_spec_compare(
      x = x, y = y, group_indicator = group_indicator,
      data_columnnames = data_columnnames,
      num_variables = num_variables,
      variable_type = variable_type, is_ordinal = is_ordinal,
      is_continuous = is_continuous,
      baseline_category = baseline_category,
      na_action = na_action, sampler = sampler,
      interaction_prior_type = interaction_prior_type,
      pairwise_scale = pairwise_scale,
      interaction_alpha = interaction_alpha,
      interaction_beta = interaction_beta,
      threshold_prior_type = threshold_prior_type,
      main_alpha = main_alpha, main_beta = main_beta,
      threshold_scale = threshold_scale,
      standardize = standardize,
      difference_selection = difference_selection,
      main_difference_selection = main_difference_selection,
      difference_prior = difference_prior,
      difference_scale = difference_scale,
      difference_probability = difference_probability,
      beta_bernoulli_alpha = beta_bernoulli_alpha,
      beta_bernoulli_beta = beta_bernoulli_beta,
      beta_bernoulli_alpha_between = difference_beta_bernoulli_alpha_between,
      beta_bernoulli_beta_between = difference_beta_bernoulli_beta_between,
      dirichlet_alpha = difference_dirichlet_alpha,
      lambda = difference_lambda
    )
  }

  validate_bgm_spec(spec)
}


# ==============================================================================
# print.bgm_spec()  --- debugging summary
# ==============================================================================
#' @export
print.bgm_spec = function(x, ...) {
  s = x
  cat("bgm_spec object\n")
  cat("  model_type:", s$model_type, "\n")
  cat("  variables: ", s$data$num_variables, " (", s$data$num_cases, " cases)\n",
    sep = ""
  )
  cat(
    "  variable_type:",
    if(s$variables$is_continuous) {
      "continuous"
    } else {
      paste0(
        sum(s$variables$is_ordinal), " ordinal, ",
        sum(!s$variables$is_ordinal), " blume-capel"
      )
    },
    "\n"
  )
  cat("  sampler:", s$sampler$update_method,
    "(iter=", s$sampler$iter, ", warmup=", s$sampler$warmup,
    ", chains=", s$sampler$chains, ")\n",
    sep = ""
  )
  if(s$model_type %in% c("ggm", "omrf")) {
    cat(
      "  edge_selection:", s$prior$edge_selection,
      if(s$prior$edge_selection) paste0(" (", s$prior$edge_prior, ")"),
      "\n"
    )
  }
  if(s$model_type == "compare") {
    cat("  groups:", s$data$num_groups, "\n")
    cat(
      "  difference_selection:", s$prior$difference_selection,
      if(s$prior$difference_selection) paste0(" (", s$prior$difference_prior, ")"),
      "\n"
    )
  }
  cat(
    "  na_action:", s$missing$na_action,
    if(s$missing$na_impute) "(imputing)" else "(complete cases)", "\n"
  )
  invisible(s)
}

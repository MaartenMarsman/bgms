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
  # All model types carry interaction_prior_type
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
    stopifnot(
      is.character(prior$difference_prior_type),
      length(prior$difference_prior_type) == 1L,
      prior$difference_prior_type %in% c("cauchy", "normal")
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

  # Continuous-block scale prior must carry a resolved raw rate
  if(mt %in% c("ggm", "mixed_mrf")) {
    if(!is.null(spec$prior$scale_rate) && !is.finite(spec$prior$scale_rate)) {
      stop(
        "bgm_spec: prior$scale_rate is not finite; a standardized-frame ",
        "scale prior (eta) was not resolved to a raw rate."
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
# zratio_joint_realized_prior_notice()
# ==============================================================================
#
# Advisory notice for the joint precision-graph specification. Under that
# specification the graph marginal is pi(Gamma) * Z(Gamma), the edge prior
# reweighted by the per-graph normalizer of the tilted precision prior, so the
# realized edge-inclusion prior is not the nominal one whatever the edge prior
# is (a uniform Beta-Bernoulli at three variables realizes ~0.37; a fixed
# bernoulli_prior(0.5) at delta = 0 realizes ~0.27; the magnitudes are in
# ?bgm). Fires whenever a continuous precision block is under edge selection,
# with two wordings: a learned inclusion probability (Beta-Bernoulli, SBM) gets
# the corrected hyperparameter update, which is coherent with the joint model
# but does not restore the nominal prior; a fixed one has nothing to absorb the
# tilt. Advisory, not a warning: the joint specification is a modelling choice.
#
# @param precision_graph_prior  "joint" or "hierarchical".
# @param model_type             Model family; only ggm/mixed_mrf carry a tilt.
# @param edge_selection         Logical.
# @param edge_prior             Resolved edge-prior name (ep_flat$edge_prior).
# @param num_continuous         Number of continuous variables.
#
# Returns: invisible TRUE when the notice fired, FALSE otherwise.
# ==============================================================================
zratio_joint_realized_prior_notice = function(precision_graph_prior, model_type,
                                              edge_selection, edge_prior,
                                              num_continuous) {
  fires = identical(precision_graph_prior, "joint") &&
    model_type %in% c("ggm", "mixed_mrf") &&
    isTRUE(edge_selection) && num_continuous >= 2
  if(!fires || !isTRUE(getOption("bgms.verbose", TRUE))) {
    return(invisible(FALSE))
  }
  learned = edge_prior %in% c("Beta-Bernoulli", "Stochastic-Block")
  message(
    "Joint precision-graph specification: the realized edge-inclusion prior ",
    "is the edge prior reweighted by the per-graph normalizer, not the ",
    "nominal edge prior",
    if(learned) {
      paste0(
        " (the hyperparameter update is corrected, so it stays coherent with ",
        "the joint model, but the realized prior still differs). "
      )
    } else {
      paste0(
        " (the inclusion probability is fixed, so nothing absorbs the tilt ",
        "and no correction applies). "
      )
    },
    "Use extract_prior_inclusion_probabilities() to read the realized prior, ",
    "or precision_graph_prior = \"hierarchical\" to target the nominal one."
  )
  invisible(TRUE)
}


# ==============================================================================
# zratio_vacuous_spec_notice()
# ==============================================================================
#
# Advisory notice for a hierarchical specification with no continuous precision
# block: an ordinal model, or mixed data with fewer than two continuous
# variables. The two specifications differ only in how p(K | Gamma) is
# normalized across graphs, so with no K there is nothing for the argument to
# refer to and the fit is the same under either value.
#
# The fixed-graph case is deliberately silent. There the argument does refer to
# something -- the specifications coincide exactly, because with no between-
# model move there is no normalizer to compare across graphs -- and a message
# would report a difference that does not exist.
#
# The default case is silent for the same reason. Since F-010 the specification
# defaults to "hierarchical", so an ordinal fit carries the value without ever
# having asked for it; the notice reports on a request, and telling a user that
# an argument they never named has no effect is noise, not information.
#
# @param has_precision_block  Whether the model carries a continuous precision
#   block of at least two variables.
# @param explicit             Whether the caller named precision_graph_prior,
#   as opposed to inheriting the default.
#
# Returns: invisible TRUE when the notice fired, FALSE otherwise.
# ==============================================================================
zratio_vacuous_spec_notice = function(has_precision_block, explicit = TRUE) {
  if(has_precision_block || !isTRUE(explicit) ||
    !isTRUE(getOption("bgms.verbose", TRUE))) {
    return(invisible(FALSE))
  }
  message(
    "precision_graph_prior has no effect for this model: it normalizes the ",
    "continuous precision prior across graphs, and this model has no ",
    "continuous precision block. The fit is the same under either value."
  )
  invisible(TRUE)
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
                    baseline_category = NULL,
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
                    scale_prior_type = "exponential",
                    scale_shape = 1,
                    scale_rate = 1,
                    scale_eta = NA_real_,
                    delta = NULL,
                    edge_selection = TRUE,
                    edge_prior = bernoulli_prior(0.5),
                    # bgm() defaults to "hierarchical" since F-010 and resolves
                    # that default itself, passing a scalar. The default here
                    # serves the only other caller, bgmCompare(), which has no
                    # continuous precision block for the argument to refer to.
                    precision_graph_prior = c("joint", "hierarchical"),
                    # FALSE when the value was inherited from a default rather
                    # than named by the user; gates the vacuous-spec advisory.
                    precision_graph_prior_explicit = TRUE,
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
                    difference_prior_type = "normal",
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
                      "adaptive-metropolis",
                      "gibbs"
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
      0.5 * log(max(sum(variable_type == "continuous"), 1))
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

  # --- Hierarchical graph-prior spec eligibility --------------------------------
  # The two specifications differ only in how p(K | Gamma) is normalized across
  # graphs, so they differ only where a between-model move exists. Where none
  # does -- a fixed graph, or no continuous precision block -- the argument is
  # vacuous rather than wrong and is accepted; the fit follows the joint path,
  # which applies no correction on exactly those configurations either
  # (ggm_edge_prior_correction() returns NULL for both), so the two coincide.
  #
  # zratio_active is resolved here and nowhere else: run_sampler_*() and
  # build_output_*() read this flag rather than re-deriving eligibility, so the
  # engine, the surface build and the trust gauge cannot disagree about whether
  # the hierarchical machinery is in force.
  #
  # Vacuity is settled before the slab, and the order is load-bearing. The
  # Z-ratio constants are derived for a Normal or Cauchy slab, so a beta-prime
  # slab is rejected -- but only where the argument refers to something. With
  # no precision block there is no slab of the precision prior for that error
  # to be about, and reporting one would name the wrong cause.
  precision_graph_prior = match.arg(precision_graph_prior)
  num_continuous = sum(variable_type == "continuous")
  has_precision_block = model_type %in% c("ggm", "mixed_mrf") &&
    num_continuous >= 2
  zratio_active = precision_graph_prior == "hierarchical" &&
    has_precision_block && isTRUE(edge_selection)
  if(precision_graph_prior == "hierarchical") {
    # Rejected exactly when the request is meaningful and unsupported; zratio_
    # active is what "meaningful" means, so a vacuous cell never reaches this.
    if(zratio_active && !interaction_prior_type %in% c("normal", "cauchy")) {
      stop(sprintf(
        paste0(
          "precision_graph_prior = \"hierarchical\" supports a normal or Cauchy ",
          "interaction (slab) prior. Got %s_prior(). Use interaction_prior = ",
          "normal_prior() or cauchy_prior(), or keep precision_graph_prior = ",
          "\"joint\"."
        ),
        interaction_prior_type
      ))
    }
    zratio_vacuous_spec_notice(has_precision_block, precision_graph_prior_explicit)
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

  zratio_joint_realized_prior_notice(
    precision_graph_prior = precision_graph_prior,
    model_type = model_type,
    edge_selection = edge_selection,
    edge_prior = ep_flat$edge_prior,
    num_continuous = sum(variable_type == "continuous")
  )

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
      scale_eta = scale_eta,
      delta = delta,
      precision_graph_prior = precision_graph_prior,
      zratio_active = zratio_active,
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
      scale_eta = scale_eta,
      delta = delta,
      precision_graph_prior = precision_graph_prior,
      zratio_active = zratio_active,
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
      difference_selection = difference_selection,
      main_difference_selection = main_difference_selection,
      difference_prior = difference_prior,
      difference_scale = difference_scale,
      difference_prior_type = difference_prior_type,
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

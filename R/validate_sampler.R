# ==============================================================================
# Sampler validation
# ==============================================================================
#
# Pure validation functions for sampler-related arguments shared by
# bgm() and bgmCompare(). Each function takes input and returns
# validated output (or errors).
# ==============================================================================


# ------------------------------------------------------------------------------
# Generic input checkers
# ------------------------------------------------------------------------------
#
# Small reusable helpers used by validate_sampler() and other validators.
# ------------------------------------------------------------------------------

check_positive_integer = function(value, name) {
  if(!is.numeric(value) || abs(value - round(value)) > .Machine$double.eps || value <= 0) {
    stop(sprintf("Parameter `%s` must be a positive integer. Got: %s", name, value))
  }
}

check_non_negative_integer = function(value, name) {
  if(!is.numeric(value) || abs(value - round(value)) > .Machine$double.eps || value < 0) {
    stop(sprintf("Parameter `%s` must be a non-negative integer. Got: %s", name, value))
  }
}

check_logical = function(value, name) {
  value = as.logical(value)
  if(is.na(value)) {
    stop(sprintf("Parameter `%s` must be TRUE or FALSE. Got: %s", name, value))
  }
  return(value)
}

check_seed = function(seed) {
  if(is.null(seed)) {
    return(sample.int(.Machine$integer.max, 1L))
  }
  if(!is.numeric(seed) || length(seed) != 1 || is.na(seed) || seed < 0) {
    stop("Argument 'seed' must be a single non-negative integer.")
  }
  as.integer(seed)
}

progress_type_from_display_progress = function(display_progress = c("per-chain", "total", "none")) {
  if(is.logical(display_progress) && length(display_progress) == 1) {
    if(is.na(display_progress)) {
      stop("The display_progress argument must be a single logical value, but not NA.")
    }
    display_progress = if(display_progress) "per-chain" else "none"
  } else {
    display_progress = match.arg(display_progress)
  }
  return(if(display_progress == "per-chain") 2L else if(display_progress == "total") 1L else 0L)
}


# ------------------------------------------------------------------------------
# validate_sampler
# ------------------------------------------------------------------------------
#
# Validates and resolves all sampler-related arguments shared by bgm()
# and bgmCompare().
#
# @param update_method  Character vector: user-supplied value (full default
#   triple means "not explicitly chosen").
# @param target_accept  Numeric or NULL. NULL = user didn't provide it;
#   will be set to a method-specific default.
# @param iter  Integer: post-warmup iterations.
# @param warmup  Integer: warmup iterations.
# @param nuts_max_depth  Integer: max tree depth for NUTS.
# @param learn_mass_matrix  Logical: adapt diagonal mass matrix during warmup.
# @param chains  Integer: number of parallel chains.
# @param cores  Integer: number of CPU cores.
# @param seed  Integer or NULL.
# @param display_progress  Character or logical: progress display mode.
# @param is_continuous  Logical: TRUE for GGM model.
# @param edge_selection  Logical: affects warmup warning tiers.
# @param verbose  Logical: whether to emit warmup warnings.
# @param progress_callback  A function or NULL: an optional progress-reporting
#   callback forwarded unchanged into the returned list.
#
# Returns:
#   list(update_method, target_accept, iter, warmup,
#        nuts_max_depth, learn_mass_matrix, chains, cores, seed, progress_type,
#        progress_callback)
# ------------------------------------------------------------------------------
validate_sampler = function(update_method,
                            target_accept = NULL,
                            iter,
                            warmup,
                            nuts_max_depth = 10,
                            learn_mass_matrix = TRUE,
                            chains = 4,
                            cores = parallel::detectCores(),
                            seed = NULL,
                            display_progress = c("per-chain", "total", "none"),
                            is_continuous = FALSE,
                            edge_selection = FALSE,
                            verbose = TRUE,
                            progress_callback = NULL) {
  # --- update_method ----------------------------------------------------------
  update_method = match.arg(
    update_method,
    choices = c("nuts", "adaptive-metropolis", "gibbs")
  )

  # "gibbs" is the conjugate row-block Gibbs sampler for the Gaussian graphical
  # model only. With edge selection it adds or removes edges with a
  # full-conditional birth/death step, for both the Normal and Cauchy slabs.
  # The slab type is not gated here; the C++ side checks that the priors are
  # supported (a Normal or Cauchy slab with a Gamma diagonal).
  if(update_method == "gibbs" && !is_continuous) {
    stop(
      "update_method = \"gibbs\" is available only for the Gaussian ",
      "graphical model (all-continuous data)."
    )
  }

  # --- target_accept ----------------------------------------------------------
  if(!is.null(target_accept)) {
    target_accept = min(target_accept, 1 - sqrt(.Machine$double.eps))
    target_accept = max(target_accept, 0 + sqrt(.Machine$double.eps))
  } else {
    target_accept = switch(update_method,
      "adaptive-metropolis" = 0.44,
      "nuts"                = 0.80,
      # Exact draw: no acceptance target. Kept numeric (0.44, inert) so the
      # downstream numeric contract holds; unused by the Gibbs path.
      "gibbs"               = 0.44
    )
  }

  # --- iter / warmup ----------------------------------------------------------
  check_positive_integer(iter, "iter")
  check_non_negative_integer(warmup, "warmup")

  # --- warmup warnings --------------------------------------------------------
  if(verbose && update_method == "nuts") {
    if(edge_selection) {
      if(warmup < 50) {
        warning(
          "warmup = ", warmup,
          " is very short for edge selection. Consider >= 300."
        )
      } else if(warmup < 200) {
        warning(
          "warmup = ", warmup,
          ": proposal SD tuning skipped (needs >= 200). Consider >= 300."
        )
      } else if(warmup < 300) {
        warning(
          "warmup = ", warmup,
          ": limited proposal SD tuning. Consider >= 300."
        )
      }
    } else {
      if(warmup < 20) {
        warning(
          "warmup = ", warmup,
          ": no mass matrix estimation (needs >= 20)."
        )
      } else if(warmup < 150) {
        warning(
          "warmup = ", warmup,
          ": using proportional allocation (needs >= 150 for fixed buffers)."
        )
      }
    }
  }

  # --- nuts_max_depth ---------------------------------------------------------
  check_positive_integer(nuts_max_depth, "nuts_max_depth")
  nuts_max_depth = max(nuts_max_depth, 1L)

  # --- learn_mass_matrix ------------------------------------------------------
  learn_mass_matrix = check_logical(learn_mass_matrix, "learn_mass_matrix")

  # --- chains / cores ---------------------------------------------------------
  check_positive_integer(chains, "chains")
  check_positive_integer(cores, "cores")

  # --- seed -------------------------------------------------------------------
  seed = check_seed(seed)

  # --- display_progress -------------------------------------------------------
  progress_type = progress_type_from_display_progress(display_progress)

  list(
    update_method = update_method,
    target_accept = target_accept,
    iter = iter,
    warmup = warmup,
    nuts_max_depth = nuts_max_depth,
    learn_mass_matrix = learn_mass_matrix,
    chains = chains,
    cores = cores,
    seed = seed,
    progress_type = progress_type,
    progress_callback = progress_callback
  )
}

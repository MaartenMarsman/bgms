# ==============================================================================
# Internal Accessor Helpers for bgms / bgmCompare Fit Objects
# ==============================================================================
#
# These functions provide a single abstraction point for reading fields
# from bgms and bgmCompare fit objects. They handle both S7 (current) and
# legacy S3 list-based objects transparently.
#
# NOT for use in the `$.bgms`/`[[.bgms` compatibility methods (which use
# S7::prop() directly).
# ==============================================================================


# ------------------------------------------------------------------
# get_fit_cache
# ------------------------------------------------------------------
# Extracts the cache environment from a fit object.
#
# @param fit  A bgms or bgmCompare object (S7 or legacy S3).
#
# Returns: The cache environment, or NULL if absent.
# ------------------------------------------------------------------
get_fit_cache = function(fit) {
  if(inherits(fit, "S7_object")) {
    fit@cache
  } else {
    .subset2(fit, "cache")
  }
}


# ------------------------------------------------------------------
# get_fit_spec
# ------------------------------------------------------------------
# Extracts the embedded bgm_spec from a fit object.
#
# @param fit  A bgms or bgmCompare object (S7 or legacy S3).
#
# Returns: The bgm_spec, or NULL if absent.
# ------------------------------------------------------------------
get_fit_spec = function(fit) {
  if(inherits(fit, "S7_object")) {
    S7::prop(fit, ".bgm_spec")
  } else {
    .subset2(fit, ".bgm_spec")
  }
}


# ------------------------------------------------------------------
# get_raw_samples
# ------------------------------------------------------------------
# Extracts the raw_samples list from a fit object.
#
# @param fit  A bgms or bgmCompare object (S7 or legacy S3).
#
# Returns: A list with components main, pairwise, indicator (if present),
#   nchains, niter, parameter_names.
# ------------------------------------------------------------------
get_raw_samples = function(fit) {
  if(inherits(fit, "S7_object")) {
    fit@raw_samples
  } else {
    .subset2(fit, "raw_samples")
  }
}


# ------------------------------------------------------------------
# get_posterior_mean
# ------------------------------------------------------------------
# Extracts a named posterior_mean_* field from a fit object.
#
# @param fit    A bgms or bgmCompare object (S7 or legacy S3).
# @param field  The suffix after "posterior_mean_", e.g. "pairwise",
#   "main", "indicator", "residual_variance", "allocations",
#   "pairwise_baseline", "pairwise_differences".
#
# Returns: The posterior mean value (matrix, vector, or list), or NULL
#   if the field does not exist.
# ------------------------------------------------------------------
get_posterior_mean = function(fit, field) {
  if(inherits(fit, "S7_object")) {
    S7::prop(fit, paste0("posterior_mean_", field))
  } else {
    .subset2(fit, paste0("posterior_mean_", field))
  }
}

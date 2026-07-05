# ==============================================================================
# S7 Class Definitions for bgms and bgmCompare Fit Objects
# ==============================================================================
#
# These S7 classes define the structure of fit objects returned by bgm()
# and bgmCompare(). They replace the previous S3 list-based representation
# while preserving the same user-facing API via $ and [[ compatibility
# methods (see methods_bgms.R and methods_bgmcompare.R).
#
# Lazy summary computation:
# The posterior_summary_* properties use S7 custom getters that trigger
# ensure_summaries() on first access. Computed values are stored in
# the cache environment (reference semantics) and returned on subsequent
# access without recomputation.
#
# names(fit) contract:
# The .field_names property stores the set of names that names(fit)
# should return, matching the previous S3 list-based behavior where
# conditional fields only appear when present. This is computed during
# construction.
# ==============================================================================


# ------------------------------------------------------------------
# bgms S7 class
# ------------------------------------------------------------------

#' @importFrom S7 new_class new_property class_any class_character prop
bgms_class = new_class("bgms",
  package = NULL,
  properties = list(
    # --- Core (always present) ---
    arguments = new_property(class_any),
    raw_samples = new_property(class_any),
    cache = new_property(class_any),

    # --- Posterior means (set during construction, immutable) ---
    posterior_mean_main = new_property(class_any, default = NULL),
    posterior_mean_pairwise = new_property(class_any, default = NULL),
    posterior_mean_residual_variance = new_property(class_any, default = NULL),
    posterior_mean_indicator = new_property(class_any, default = NULL),
    posterior_mean_coclustering_matrix = new_property(class_any, default = NULL),
    posterior_mean_allocations = new_property(class_any, default = NULL),
    posterior_mode_allocations = new_property(class_any, default = NULL),
    posterior_num_blocks = new_property(class_any, default = NULL),

    # --- Pre-computed summaries (not lazy) ---
    posterior_summary_pairwise_allocations = new_property(class_any, default = NULL),

    # --- Sampled edge-prior inclusion parameter (per-chain list) ---
    inclusion_parameter_samples = new_property(class_any, default = NULL),

    # --- Lazy MCMC diagnostics (computed on first access via getter) ---
    posterior_summary_main = new_property(
      class = class_any,
      getter = function(self) {
        ensure_summaries(self)
        self@cache[["posterior_summary_main"]]
      }
    ),
    posterior_summary_pairwise = new_property(
      class = class_any,
      getter = function(self) {
        ensure_summaries(self)
        self@cache[["posterior_summary_pairwise"]]
      }
    ),
    posterior_summary_indicator = new_property(
      class = class_any,
      getter = function(self) {
        ensure_summaries(self)
        self@cache[["posterior_summary_indicator"]]
      }
    ),
    posterior_summary_quadratic = new_property(
      class = class_any,
      getter = function(self) {
        ensure_summaries(self)
        self@cache[["posterior_summary_quadratic"]]
      }
    ),

    # --- Optional ---
    nuts_diag = new_property(class_any, default = NULL),
    am_diag = new_property(class_any, default = NULL),
    zratio_diag = new_property(class_any, default = NULL),

    # --- easybgm compatibility (deprecated) ---
    indicator = new_property(class_any, default = NULL),
    interactions = new_property(class_any, default = NULL),
    thresholds = new_property(class_any, default = NULL),

    # --- Internal ---
    .bgm_spec = new_property(class_any, default = NULL),
    .field_names = new_property(class_character)
  )
)


# ------------------------------------------------------------------
# s3_list_to_bgms
# ------------------------------------------------------------------
# Converts an S3 list-based bgms fit (built incrementally in
# build_output_bgm / build_output_mixed_mrf) to the S7 bgms_class.
#
# @param results  A named list with class "bgms".
#
# Returns: A bgms_class S7 object.
# ------------------------------------------------------------------
s3_list_to_bgms = function(results) {
  bgms_class(
    arguments = .subset2(results, "arguments"),
    raw_samples = .subset2(results, "raw_samples"),
    cache = .subset2(results, "cache"),
    posterior_mean_main = .subset2(results, "posterior_mean_main"),
    posterior_mean_pairwise = .subset2(results, "posterior_mean_pairwise"),
    posterior_mean_residual_variance = .subset2(results, "posterior_mean_residual_variance"),
    posterior_mean_indicator = .subset2(results, "posterior_mean_indicator"),
    posterior_mean_coclustering_matrix = .subset2(results, "posterior_mean_coclustering_matrix"),
    posterior_mean_allocations = .subset2(results, "posterior_mean_allocations"),
    posterior_mode_allocations = .subset2(results, "posterior_mode_allocations"),
    posterior_num_blocks = .subset2(results, "posterior_num_blocks"),
    posterior_summary_pairwise_allocations = .subset2(results, "posterior_summary_pairwise_allocations"),
    inclusion_parameter_samples = .subset2(results, "inclusion_parameter_samples"),
    nuts_diag = .subset2(results, "nuts_diag"),
    am_diag = .subset2(results, "am_diag"),
    zratio_diag = .subset2(results, "zratio_diag"),
    indicator = .subset2(results, "indicator"),
    interactions = .subset2(results, "interactions"),
    thresholds = .subset2(results, "thresholds"),
    .bgm_spec = .subset2(results, ".bgm_spec"),
    .field_names = names(results)
  )
}


# ------------------------------------------------------------------
# bgmCompare S7 class
# ------------------------------------------------------------------

bgmCompare_class = new_class("bgmCompare",
  package = NULL,
  properties = list(
    # --- Core (always present) ---
    arguments = new_property(class_any),
    raw_samples = new_property(class_any),
    cache = new_property(class_any),

    # --- Posterior means (set during construction, immutable) ---
    posterior_mean_main_baseline = new_property(class_any, default = NULL),
    posterior_mean_pairwise_baseline = new_property(class_any, default = NULL),
    posterior_mean_main_differences = new_property(class_any, default = NULL),
    posterior_mean_pairwise_differences = new_property(class_any, default = NULL),

    # --- SBM difference-prior outputs (Stochastic-Block only) ---
    posterior_mean_coclustering_matrix = new_property(class_any, default = NULL),
    posterior_mean_allocations = new_property(class_any, default = NULL),
    posterior_mode_allocations = new_property(class_any, default = NULL),
    posterior_num_blocks = new_property(class_any, default = NULL),
    posterior_summary_pairwise_allocations = new_property(class_any, default = NULL),

    # --- Lazy MCMC diagnostics (computed on first access via getter) ---
    posterior_summary_main_baseline = new_property(
      class = class_any,
      getter = function(self) {
        ensure_summaries(self)
        self@cache[["posterior_summary_main_baseline"]]
      }
    ),
    posterior_summary_pairwise_baseline = new_property(
      class = class_any,
      getter = function(self) {
        ensure_summaries(self)
        self@cache[["posterior_summary_pairwise_baseline"]]
      }
    ),
    posterior_summary_main_differences = new_property(
      class = class_any,
      getter = function(self) {
        ensure_summaries(self)
        self@cache[["posterior_summary_main_differences"]]
      }
    ),
    posterior_summary_pairwise_differences = new_property(
      class = class_any,
      getter = function(self) {
        ensure_summaries(self)
        self@cache[["posterior_summary_pairwise_differences"]]
      }
    ),
    posterior_summary_indicator = new_property(
      class = class_any,
      getter = function(self) {
        ensure_summaries(self)
        self@cache[["posterior_summary_indicator"]]
      }
    ),

    # --- Optional ---
    nuts_diag = new_property(class_any, default = NULL),
    am_diag = new_property(class_any, default = NULL),

    # --- Internal ---
    .bgm_spec = new_property(class_any, default = NULL),
    .field_names = new_property(class_character)
  )
)


# ------------------------------------------------------------------
# s3_list_to_bgmCompare
# ------------------------------------------------------------------
# Converts an S3 list-based bgmCompare fit (built incrementally in
# build_output_compare) to the S7 bgmCompare_class.
#
# @param results  A named list with class "bgmCompare".
#
# Returns: A bgmCompare_class S7 object.
# ------------------------------------------------------------------
s3_list_to_bgmCompare = function(results) {
  bgmCompare_class(
    arguments = .subset2(results, "arguments"),
    raw_samples = .subset2(results, "raw_samples"),
    cache = .subset2(results, "cache"),
    posterior_mean_main_baseline = .subset2(results, "posterior_mean_main_baseline"),
    posterior_mean_pairwise_baseline = .subset2(results, "posterior_mean_pairwise_baseline"),
    posterior_mean_main_differences = .subset2(results, "posterior_mean_main_differences"),
    posterior_mean_pairwise_differences = .subset2(results, "posterior_mean_pairwise_differences"),
    posterior_mean_coclustering_matrix = .subset2(results, "posterior_mean_coclustering_matrix"),
    posterior_mean_allocations = .subset2(results, "posterior_mean_allocations"),
    posterior_mode_allocations = .subset2(results, "posterior_mode_allocations"),
    posterior_num_blocks = .subset2(results, "posterior_num_blocks"),
    posterior_summary_pairwise_allocations = .subset2(results, "posterior_summary_pairwise_allocations"),
    nuts_diag = .subset2(results, "nuts_diag"),
    am_diag = .subset2(results, "am_diag"),
    .bgm_spec = .subset2(results, ".bgm_spec"),
    .field_names = names(results)
  )
}

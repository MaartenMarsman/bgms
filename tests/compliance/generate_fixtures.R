# ==============================================================================
# Generate Bitwise Compliance Fixtures (baseline: CRAN bgms 0.1.6.3)
# ==============================================================================
#
# Installs bgms 0.1.6.3 from CRAN into an isolated library and runs all
# configurations. The resulting fixtures are the ground-truth reference for
# OMRF/Blume-Capel models (bgm + bgmCompare).
#
# These fixtures capture the LAST CRAN release before continuous variables
# were added. All future development must preserve bitwise-identical output
# for discrete (ordinal / binary / Blume-Capel) models.
#
# Usage:
#   Rscript tests/compliance/generate_fixtures.R
#
# Output:
#   tests/compliance/fixtures/  <U+2014> one .rds per configuration + manifest.rds
#
# ==============================================================================

library(callr)

fixture_dir = file.path("tests", "compliance", "fixtures")
dir.create(fixture_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Install CRAN 0.1.6.3 into a temporary library
# ==============================================================================

cran_lib = tempfile("bgms_cran_")
dir.create(cran_lib)

baseline_version = "0.1.6.3"

# The baseline must be 0.1.6.3 exactly. Take it from the current CRAN release
# while that still is 0.1.6.3; once a newer bgms is on CRAN, fall back to the
# pinned tarball in the CRAN Archive. Silently baselining against a newer
# release would make every comparison below vacuous.
cran_version = tryCatch(
  as.character(available.packages(repos = "https://cloud.r-project.org")["bgms", "Version"]),
  error = function(e) NA_character_
)

if(identical(cran_version, baseline_version)) {
  cat("Installing bgms", baseline_version, "from CRAN...\n")
  install.packages(
    "bgms",
    repos = "https://cloud.r-project.org",
    lib = cran_lib,
    quiet = TRUE
  )
} else {
  archive_url = sprintf(
    "https://cran.r-project.org/src/contrib/Archive/bgms/bgms_%s.tar.gz",
    baseline_version
  )
  cat(
    "CRAN currently serves", cran_version, "- installing pinned",
    baseline_version, "from the CRAN Archive...\n"
  )
  install.packages(archive_url, repos = NULL, type = "source", lib = cran_lib, quiet = TRUE)
}

installed_version = callr::r(
  function(lib_path) {
    as.character(packageVersion("bgms", lib.loc = lib_path))
  },
  args = list(lib_path = cran_lib)
)
cat("Installed version:", installed_version, "\n")

if(installed_version != baseline_version) {
  stop(sprintf(
    "Compliance baseline must be bgms %s, but %s was installed. Fixtures not generated.",
    baseline_version, installed_version
  ))
}

# ==============================================================================
# Configuration matrix <U+2014> OMRF/Blume-Capel bgm() configurations
# ==============================================================================
#
# Covers:
#   Datasets:       Wenchuan (ordinal 5-pt), ADHD (binary), Boredom (ordinal 7-pt)
#   Variable types: ordinal (default), blume-capel
#   Samplers:       NUTS, adaptive-metropolis
#   Edge priors:    Bernoulli, Beta-Bernoulli, Stochastic-Block
#   Edge selection: TRUE, FALSE
#   Missing data:   listwise, impute
#
# Short iterations (200 iter / 200 warmup / 2 chains) <U+2014> bitwise identity does
# not require convergence.
# ==============================================================================

bgm_configs = list(
  # --- Wenchuan ordinal ---
  list(
    id = "bgm_wenchuan_nuts_bernoulli",
    desc = "bgm: Wenchuan 6v, NUTS, Bernoulli, edge_sel",
    args = list(
      x = "wenchuan_small", iter = 200, warmup = 200, chains = 2,
      edge_selection = TRUE, edge_prior = "Bernoulli",
      update_method = "nuts", seed = 101, display_progress = "none"
    )
  ),
  list(
    id = "bgm_wenchuan_nuts_no_edgesel",
    desc = "bgm: Wenchuan 6v, NUTS, no edge selection",
    args = list(
      x = "wenchuan_small", iter = 200, warmup = 200, chains = 2,
      edge_selection = FALSE,
      update_method = "nuts", seed = 102, display_progress = "none"
    )
  ),
  list(
    id = "bgm_wenchuan_mh_bernoulli",
    desc = "bgm: Wenchuan 6v, MH, Bernoulli, edge_sel",
    args = list(
      x = "wenchuan_small", iter = 200, warmup = 200, chains = 2,
      edge_selection = TRUE, edge_prior = "Bernoulli",
      update_method = "adaptive-metropolis", seed = 103, display_progress = "none"
    )
  ),
  list(
    id = "bgm_wenchuan_nuts_betabern",
    desc = "bgm: Wenchuan 6v, NUTS, Beta-Bernoulli, edge_sel",
    args = list(
      x = "wenchuan_small", iter = 200, warmup = 200, chains = 2,
      edge_selection = TRUE, edge_prior = "Beta-Bernoulli",
      beta_bernoulli_alpha = 1, beta_bernoulli_beta = 1,
      update_method = "nuts", seed = 105, display_progress = "none"
    )
  ),
  list(
    id = "bgm_wenchuan_nuts_sbm",
    desc = "bgm: Wenchuan 6v, NUTS, Stochastic-Block, edge_sel",
    args = list(
      x = "wenchuan_small", iter = 200, warmup = 200, chains = 2,
      edge_selection = TRUE, edge_prior = "Stochastic-Block",
      update_method = "nuts", seed = 106, display_progress = "none"
    )
  ),
  list(
    id = "bgm_wenchuan_nuts_scaled_prior",
    desc = "bgm: Wenchuan 6v, NUTS, pairwise_scale=1.0",
    args = list(
      x = "wenchuan_small", iter = 200, warmup = 200, chains = 2,
      edge_selection = TRUE, edge_prior = "Bernoulli",
      pairwise_scale = 1.0,
      update_method = "nuts", seed = 107, display_progress = "none"
    )
  ),

  # --- Wenchuan Blume-Capel ---

  list(
    id = "bgm_wenchuan_nuts_blumecapel",
    desc = "bgm: Wenchuan 6v, NUTS, Blume-Capel, baseline=2",
    args = list(
      x = "wenchuan_small", iter = 200, warmup = 200, chains = 2,
      variable_type = "blume-capel", baseline_category = 2,
      edge_selection = TRUE, edge_prior = "Bernoulli",
      update_method = "nuts", seed = 401, display_progress = "none"
    )
  ),
  list(
    id = "bgm_wenchuan_mh_blumecapel",
    desc = "bgm: Wenchuan 6v, MH, Blume-Capel, baseline=2",
    args = list(
      x = "wenchuan_small", iter = 200, warmup = 200, chains = 2,
      variable_type = "blume-capel", baseline_category = 2,
      edge_selection = TRUE, edge_prior = "Bernoulli",
      update_method = "adaptive-metropolis", seed = 402, display_progress = "none"
    )
  ),
  list(
    id = "bgm_wenchuan_nuts_blumecapel_no_edgesel",
    desc = "bgm: Wenchuan 6v, NUTS, Blume-Capel, no edge_sel",
    args = list(
      x = "wenchuan_small", iter = 200, warmup = 200, chains = 2,
      variable_type = "blume-capel", baseline_category = 2,
      edge_selection = FALSE,
      update_method = "nuts", seed = 404, display_progress = "none"
    )
  ),
  list(
    id = "bgm_wenchuan_nuts_blumecapel_baseline1",
    desc = "bgm: Wenchuan 6v, NUTS, Blume-Capel, baseline=1",
    args = list(
      x = "wenchuan_small", iter = 200, warmup = 200, chains = 2,
      variable_type = "blume-capel", baseline_category = 1,
      edge_selection = TRUE, edge_prior = "Bernoulli",
      update_method = "nuts", seed = 405, display_progress = "none"
    )
  ),

  # --- ADHD binary ---

  list(
    id = "bgm_adhd_nuts_bernoulli",
    desc = "bgm: ADHD 6v, NUTS, Bernoulli, edge_sel",
    args = list(
      x = "adhd_small", iter = 200, warmup = 200, chains = 2,
      edge_selection = TRUE, edge_prior = "Bernoulli",
      update_method = "nuts", seed = 201, display_progress = "none"
    )
  ),
  list(
    id = "bgm_adhd_mh_bernoulli",
    desc = "bgm: ADHD 6v, MH, Bernoulli, edge_sel",
    args = list(
      x = "adhd_small", iter = 200, warmup = 200, chains = 2,
      edge_selection = TRUE, edge_prior = "Bernoulli",
      update_method = "adaptive-metropolis", seed = 202, display_progress = "none"
    )
  ),
  list(
    id = "bgm_adhd_nuts_no_edgesel",
    desc = "bgm: ADHD 6v, NUTS, no edge selection",
    args = list(
      x = "adhd_small", iter = 200, warmup = 200, chains = 2,
      edge_selection = FALSE,
      update_method = "nuts", seed = 203, display_progress = "none"
    )
  ),
  list(
    id = "bgm_adhd_nuts_sbm",
    desc = "bgm: ADHD 6v, NUTS, Stochastic-Block, edge_sel",
    args = list(
      x = "adhd_small", iter = 200, warmup = 200, chains = 2,
      edge_selection = TRUE, edge_prior = "Stochastic-Block",
      update_method = "nuts", seed = 204, display_progress = "none"
    )
  ),

  # --- Boredom ordinal 7-point ---

  list(
    id = "bgm_boredom_nuts_bernoulli",
    desc = "bgm: Boredom 6v, NUTS, Bernoulli, edge_sel",
    args = list(
      x = "boredom_small", iter = 200, warmup = 200, chains = 2,
      edge_selection = TRUE, edge_prior = "Bernoulli",
      update_method = "nuts", seed = 301, display_progress = "none"
    )
  ),
  list(
    id = "bgm_boredom_mh_betabern",
    desc = "bgm: Boredom 6v, MH, Beta-Bernoulli, edge_sel",
    args = list(
      x = "boredom_small", iter = 200, warmup = 200, chains = 2,
      edge_selection = TRUE, edge_prior = "Beta-Bernoulli",
      update_method = "adaptive-metropolis", seed = 302, display_progress = "none"
    )
  ),
  list(
    id = "bgm_boredom_nuts_no_edgesel",
    desc = "bgm: Boredom 6v, NUTS, no edge selection",
    args = list(
      x = "boredom_small", iter = 200, warmup = 200, chains = 2,
      edge_selection = FALSE,
      update_method = "nuts", seed = 304, display_progress = "none"
    )
  ),

  # --- Missing data imputation ---

  list(
    id = "bgm_wenchuan_nuts_impute",
    desc = "bgm: Wenchuan 6v +NAs, NUTS, impute",
    args = list(
      x = "wenchuan_na", iter = 200, warmup = 200, chains = 2,
      na_action = "impute",
      edge_selection = TRUE, edge_prior = "Bernoulli",
      update_method = "nuts", seed = 501, display_progress = "none"
    )
  ),
  list(
    id = "bgm_wenchuan_nuts_blumecapel_impute",
    desc = "bgm: Wenchuan 6v +NAs, NUTS, Blume-Capel, impute",
    args = list(
      x = "wenchuan_na", iter = 200, warmup = 200, chains = 2,
      variable_type = "blume-capel", baseline_category = 2,
      na_action = "impute",
      edge_selection = TRUE, edge_prior = "Bernoulli",
      update_method = "nuts", seed = 502, display_progress = "none"
    )
  )
)

# ==============================================================================
# bgmCompare() configurations
# ==============================================================================
#
# Uses x/y API (two separate data frames) since that is how 0.1.6.3 was
# primarily tested. Covers ordinal, binary, and Blume-Capel with
# difference_selection and main_difference_selection.
# ==============================================================================

compare_configs = list(
  # --- Wenchuan ordinal (split into two groups) ---
  list(
    id = "cmp_wenchuan_nuts_bernoulli",
    desc = "bgmCompare: Wenchuan, NUTS, Bernoulli, diff_sel",
    args = list(
      x = "wenchuan_g1", y = "wenchuan_g2",
      iter = 200, warmup = 200, chains = 2,
      difference_selection = TRUE,
      update_method = "nuts", seed = 601, display_progress = "none"
    )
  ),
  list(
    id = "cmp_wenchuan_mh_bernoulli",
    desc = "bgmCompare: Wenchuan, MH, Bernoulli, diff_sel",
    args = list(
      x = "wenchuan_g1", y = "wenchuan_g2",
      iter = 200, warmup = 200, chains = 2,
      difference_selection = TRUE,
      update_method = "adaptive-metropolis", seed = 602, display_progress = "none"
    )
  ),
  list(
    id = "cmp_wenchuan_nuts_no_diffsel",
    desc = "bgmCompare: Wenchuan, NUTS, no diff selection",
    args = list(
      x = "wenchuan_g1", y = "wenchuan_g2",
      iter = 200, warmup = 200, chains = 2,
      difference_selection = FALSE,
      update_method = "nuts", seed = 604, display_progress = "none"
    )
  ),
  list(
    id = "cmp_wenchuan_nuts_main_diffsel",
    desc = "bgmCompare: Wenchuan, NUTS, main+pairwise diff_sel",
    args = list(
      x = "wenchuan_g1", y = "wenchuan_g2",
      iter = 200, warmup = 200, chains = 2,
      difference_selection = TRUE,
      main_difference_selection = TRUE,
      update_method = "nuts", seed = 605, display_progress = "none"
    )
  ),
  list(
    id = "cmp_wenchuan_nuts_betabern",
    desc = "bgmCompare: Wenchuan, NUTS, Beta-Bernoulli, diff_sel",
    args = list(
      x = "wenchuan_g1", y = "wenchuan_g2",
      iter = 200, warmup = 200, chains = 2,
      difference_selection = TRUE,
      difference_prior = "Beta-Bernoulli",
      update_method = "nuts", seed = 606, display_progress = "none"
    )
  ),

  # --- Wenchuan Blume-Capel ---

  list(
    id = "cmp_wenchuan_nuts_blumecapel",
    desc = "bgmCompare: Wenchuan, NUTS, Blume-Capel, baseline=2, diff_sel",
    args = list(
      x = "wenchuan_g1", y = "wenchuan_g2",
      iter = 200, warmup = 200, chains = 2,
      variable_type = "blume-capel", baseline_category = 2,
      difference_selection = TRUE,
      update_method = "nuts", seed = 701, display_progress = "none"
    )
  ),
  list(
    id = "cmp_wenchuan_mh_blumecapel",
    desc = "bgmCompare: Wenchuan, MH, Blume-Capel, baseline=2, diff_sel",
    args = list(
      x = "wenchuan_g1", y = "wenchuan_g2",
      iter = 200, warmup = 200, chains = 2,
      variable_type = "blume-capel", baseline_category = 2,
      difference_selection = TRUE,
      update_method = "adaptive-metropolis", seed = 702, display_progress = "none"
    )
  ),
  list(
    id = "cmp_wenchuan_nuts_blumecapel_no_diffsel",
    desc = "bgmCompare: Wenchuan, NUTS, Blume-Capel, no diff_sel",
    args = list(
      x = "wenchuan_g1", y = "wenchuan_g2",
      iter = 200, warmup = 200, chains = 2,
      variable_type = "blume-capel", baseline_category = 2,
      difference_selection = FALSE,
      update_method = "nuts", seed = 703, display_progress = "none"
    )
  ),

  # --- ADHD binary ---

  list(
    id = "cmp_adhd_nuts_bernoulli",
    desc = "bgmCompare: ADHD, NUTS, Bernoulli, diff_sel",
    args = list(
      x = "adhd_g1", y = "adhd_g2",
      iter = 200, warmup = 200, chains = 2,
      difference_selection = TRUE,
      update_method = "nuts", seed = 801, display_progress = "none"
    )
  ),
  list(
    id = "cmp_adhd_mh_bernoulli",
    desc = "bgmCompare: ADHD, MH, Bernoulli, diff_sel",
    args = list(
      x = "adhd_g1", y = "adhd_g2",
      iter = 200, warmup = 200, chains = 2,
      difference_selection = TRUE,
      update_method = "adaptive-metropolis", seed = 802, display_progress = "none"
    )
  ),

  # --- Boredom ordinal 7-point ---

  list(
    id = "cmp_boredom_nuts_bernoulli",
    desc = "bgmCompare: Boredom, NUTS, Bernoulli, diff_sel",
    args = list(
      x = "boredom_g1", y = "boredom_g2",
      iter = 200, warmup = 200, chains = 2,
      difference_selection = TRUE,
      update_method = "nuts", seed = 901, display_progress = "none"
    )
  ),
  list(
    id = "cmp_boredom_mh_bernoulli",
    desc = "bgmCompare: Boredom, MH, Bernoulli, diff_sel",
    args = list(
      x = "boredom_g1", y = "boredom_g2",
      iter = 200, warmup = 200, chains = 2,
      difference_selection = TRUE,
      update_method = "adaptive-metropolis", seed = 902, display_progress = "none"
    )
  ),

  # --- Missing data imputation in bgmCompare ---

  list(
    id = "cmp_wenchuan_nuts_impute",
    desc = "bgmCompare: Wenchuan +NAs, NUTS, impute, diff_sel",
    args = list(
      x = "wenchuan_g1_na", y = "wenchuan_g2_na",
      iter = 200, warmup = 200, chains = 2,
      na_action = "impute",
      difference_selection = TRUE,
      update_method = "nuts", seed = 1001, display_progress = "none"
    )
  )
)

# ==============================================================================
# Dataset preparation
# ==============================================================================

prepare_datasets = function() {
  data(Wenchuan, package = "bgms", envir = environment())
  data(ADHD, package = "bgms", envir = environment())
  data(Boredom, package = "bgms", envir = environment())

  wenchuan_small = Wenchuan[, 1:6]
  adhd_small = ADHD[, 2:7]
  boredom_small = Boredom[, 2:7]

  # Wenchuan with NAs (deterministic injection)
  wenchuan_na = wenchuan_small
  set.seed(999)
  na_idx = sample(length(wenchuan_na), size = 20)
  wenchuan_na[na_idx] = NA

  # Split datasets for bgmCompare (deterministic halves)
  n_w = nrow(wenchuan_small)
  wenchuan_g1 = wenchuan_small[1:floor(n_w / 2), ]
  wenchuan_g2 = wenchuan_small[(floor(n_w / 2) + 1):n_w, ]

  n_a = nrow(adhd_small)
  adhd_g1 = adhd_small[1:floor(n_a / 2), ]
  adhd_g2 = adhd_small[(floor(n_a / 2) + 1):n_a, ]

  n_b = nrow(boredom_small)
  boredom_g1 = boredom_small[1:floor(n_b / 2), ]
  boredom_g2 = boredom_small[(floor(n_b / 2) + 1):n_b, ]

  # Wenchuan groups with NAs
  wenchuan_g1_na = wenchuan_g1
  wenchuan_g2_na = wenchuan_g2
  set.seed(998)
  na_idx1 = sample(length(wenchuan_g1_na), size = 10)
  wenchuan_g1_na[na_idx1] = NA
  set.seed(997)
  na_idx2 = sample(length(wenchuan_g2_na), size = 10)
  wenchuan_g2_na[na_idx2] = NA

  list(
    wenchuan_small = wenchuan_small,
    adhd_small     = adhd_small,
    boredom_small  = boredom_small,
    wenchuan_na    = wenchuan_na,
    wenchuan_g1    = wenchuan_g1,
    wenchuan_g2    = wenchuan_g2,
    adhd_g1        = adhd_g1,
    adhd_g2        = adhd_g2,
    boredom_g1     = boredom_g1,
    boredom_g2     = boredom_g2,
    wenchuan_g1_na = wenchuan_g1_na,
    wenchuan_g2_na = wenchuan_g2_na
  )
}

# ==============================================================================
# Fixture extraction helpers
# ==============================================================================

extract_bgm_fixture = function(fit, config) {
  list(
    id = config$id,
    desc = config$desc,
    posterior_summary_main = fit$posterior_summary_main,
    posterior_summary_pairwise = fit$posterior_summary_pairwise,
    posterior_summary_indicator = fit$posterior_summary_indicator,
    posterior_mean_main = fit$posterior_mean_main,
    posterior_mean_pairwise = fit$posterior_mean_pairwise,
    posterior_mean_indicator = fit$posterior_mean_indicator,
    raw_main_chain1 = fit$raw_samples$main[[1]],
    raw_pairwise_chain1 = fit$raw_samples$pairwise[[1]],
    raw_indicator_chain1 = if(!is.null(fit$raw_samples$indicator)) {
      fit$raw_samples$indicator[[1]]
    } else {
      NULL
    },
    nuts_diag = fit$nuts_diag,
    posterior_mean_coclustering_matrix = fit$posterior_mean_coclustering_matrix,
    posterior_mean_allocations = fit$posterior_mean_allocations,
    bgms_version = as.character(packageVersion("bgms"))
  )
}

extract_compare_fixture = function(fit, config) {
  list(
    id = config$id,
    desc = config$desc,
    posterior_summary_main_baseline = fit$posterior_summary_main_baseline,
    posterior_summary_pairwise_baseline = fit$posterior_summary_pairwise_baseline,
    posterior_summary_main_differences = fit$posterior_summary_main_differences,
    posterior_summary_pairwise_differences = fit$posterior_summary_pairwise_differences,
    posterior_summary_indicator = fit$posterior_summary_indicator,
    posterior_mean_main_baseline = fit$posterior_mean_main_baseline,
    posterior_mean_pairwise_baseline = fit$posterior_mean_pairwise_baseline,
    posterior_mean_main_differences = fit$posterior_mean_main_differences,
    posterior_mean_pairwise_differences = fit$posterior_mean_pairwise_differences,
    # bgmCompare never produced a posterior_mean_indicator field (0.1.6.3 sets
    # it only on the bgm path); difference inclusion probabilities are in
    # posterior_summary_indicator.
    raw_samples = fit$raw_samples,
    nuts_diag = fit$nuts_diag,
    bgms_version = as.character(packageVersion("bgms"))
  )
}

# ==============================================================================
# Generate fixtures via callr (isolated CRAN library)
# ==============================================================================

resolve_args = function(args, datasets) {
  resolved = args
  for(nm in names(resolved)) {
    if(is.character(resolved[[nm]]) && resolved[[nm]] %in% names(datasets)) {
      resolved[[nm]] = datasets[[resolved[[nm]]]]
    }
  }
  resolved
}

datasets = prepare_datasets()

# --- bgm fixtures ---

cat(sprintf("\nGenerating %d bgm fixtures...\n", length(bgm_configs)))
bgm_manifest = list()

for(config in bgm_configs) {
  cat(sprintf("  [%s] %s ... ", config$id, config$desc))
  resolved = resolve_args(config$args, datasets)

  result = tryCatch(
    {
      callr::r(
        function(args, lib_path, extract_fn) {
          .libPaths(c(lib_path, .libPaths()))
          library(bgms, lib.loc = lib_path)
          set.seed(args$seed)
          fit = do.call(bgm, args)
          extract_fn(fit, list(id = "tmp", desc = "tmp"))
        },
        args = list(args = resolved, lib_path = cran_lib, extract_fn = extract_bgm_fixture),
        show = FALSE
      )
    },
    error = function(e) {
      cat(sprintf("ERROR: %s\n", conditionMessage(e)))
      NULL
    }
  )

  if(!is.null(result)) {
    result$id = config$id
    result$desc = config$desc
    result$bgms_version = installed_version
    path = file.path(fixture_dir, paste0(config$id, ".rds"))
    saveRDS(result, path)
    bgm_manifest[[config$id]] = list(
      id = config$id, desc = config$desc, type = "bgm",
      file = paste0(config$id, ".rds")
    )
    cat("OK\n")
  }
}

# --- bgmCompare fixtures ---

cat(sprintf("\nGenerating %d bgmCompare fixtures...\n", length(compare_configs)))
compare_manifest = list()

for(config in compare_configs) {
  cat(sprintf("  [%s] %s ... ", config$id, config$desc))
  resolved = resolve_args(config$args, datasets)

  result = tryCatch(
    {
      callr::r(
        function(args, lib_path, extract_fn) {
          .libPaths(c(lib_path, .libPaths()))
          library(bgms, lib.loc = lib_path)
          set.seed(args$seed)
          fit = do.call(bgmCompare, args)
          extract_fn(fit, list(id = "tmp", desc = "tmp"))
        },
        args = list(args = resolved, lib_path = cran_lib, extract_fn = extract_compare_fixture),
        show = FALSE
      )
    },
    error = function(e) {
      cat(sprintf("ERROR: %s\n", conditionMessage(e)))
      NULL
    }
  )

  if(!is.null(result)) {
    result$id = config$id
    result$desc = config$desc
    result$bgms_version = installed_version
    path = file.path(fixture_dir, paste0(config$id, ".rds"))
    saveRDS(result, path)
    compare_manifest[[config$id]] = list(
      id = config$id, desc = config$desc, type = "bgmCompare",
      file = paste0(config$id, ".rds")
    )
    cat("OK\n")
  }
}

# ==============================================================================
# Save manifest
# ==============================================================================

manifest = c(bgm_manifest, compare_manifest)
saveRDS(manifest, file.path(fixture_dir, "manifest.rds"))

cat(sprintf(
  "\nDone. %d bgm + %d bgmCompare = %d total fixtures in %s\n",
  length(bgm_manifest), length(compare_manifest), length(manifest), fixture_dir
))
cat("Baseline version:", installed_version, "\n")
cat("R version:", R.version.string, "\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")

# Cleanup
unlink(cran_lib, recursive = TRUE)

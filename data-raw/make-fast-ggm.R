# ==============================================================================
# Generator for vignettes/fast-ggm.rds
# ==============================================================================
# The fast-GGM vignette reports wall-clock timings. Timings are machine facts,
# and a fit long enough to time is far too long to run inside R CMD check, so
# the runs happen here, once, and the vignette ships and prints the numbers they
# produced. This follows the pattern of data-raw/make-prior-sensitivity.R.
#
# Run from the package root:
#
#     Rscript data-raw/make-fast-ggm.R
#
# Every call below MUST stay identical to the call the vignette displays: that
# displayed call is the claim that these numbers came from it. The seeds are
# fixed, and each fit runs on its own, so the elapsed times are not contaminated
# by a neighbouring job.
#
# Machine budget: the runs use 4 chains on 4 threads (chains = 4, cores = 4),
# which is what the vignette states. Nothing here should be run in parallel with
# anything else.
#
# Total runtime is dominated by the 200-variable fit; the whole script is of the
# order of ten minutes on the machine recorded in the payload.
# ==============================================================================

library(bgms)

CHAINS = 4L
CORES = 4L

# ------------------------------------------------------------------------------
# A sparse Gaussian graphical model. The graph is Erdos-Renyi with an average
# degree of three; the precision matrix is made diagonally dominant so it is
# positive definite by construction, at any number of variables.
# ------------------------------------------------------------------------------
simulate_sparse_ggm = function(p, avg_degree = 3, weight = 0.25, seed) {
  set.seed(seed)
  K = matrix(0, p, p)
  pairs = which(upper.tri(K), arr.ind = TRUE)
  present = stats::runif(nrow(pairs)) < avg_degree / (p - 1)
  K[pairs[present, , drop = FALSE]] =
    sample(c(-1, 1), sum(present), replace = TRUE) * weight
  K = K + t(K)
  diag(K) = rowSums(abs(K)) + 1
  K
}

simulate_case = function(p, n, seed) {
  K = simulate_sparse_ggm(p, seed = seed)
  list(
    K = K,
    y = simulate_mrf(
      num_states = n, num_variables = p, pairwise = K,
      variable_type = "continuous", seed = seed
    ),
    present = K[upper.tri(K)] != 0
  )
}

edge_probabilities = function(fit) {
  pip = extract_posterior_inclusion_probabilities(fit)
  pip[upper.tri(pip)]
}

# Recovery against the graph the data came from, at the conventional cut. The
# payload carries these counts and never the inclusion-probability vector
# itself: at 200 variables that vector alone is 20,000 numbers, and the
# vignette quotes summaries of it rather than plotting it.
recovery = function(fit, present) {
  pip = edge_probabilities(fit)
  list(
    selected = sum(pip > 0.5),
    detected = sum(pip[present] > 0.5),
    false_positive = sum(pip[!present] > 0.5),
    max_rhat = max(summary(fit)$indicator[, "Rhat"], na.rm = TRUE)
  )
}

payload = list()

# ------------------------------------------------------------------------------
# Machine and software record. Every timing in the vignette is a fact about this
# machine and this build, and is worthless without it.
# ------------------------------------------------------------------------------
payload$machine = list(
  cpu = tryCatch(
    system("sysctl -n machdep.cpu.brand_string", intern = TRUE),
    error = function(e) NA_character_
  ),
  cores_total = parallel::detectCores(),
  os = paste(Sys.info()[["sysname"]], Sys.info()[["release"]]),
  r_version = R.version.string,
  platform = R.version$platform,
  bgms_version = as.character(utils::packageVersion("bgms")),
  chains = CHAINS,
  threads = CORES,
  date = as.character(Sys.Date())
)

# ==============================================================================
# 1. The mid-sized model, run both ways
# ==============================================================================
mid = simulate_case(p = 50, n = 1000, seed = 2026)

payload$mid = list(
  p = 50, n = 1000, edges = sum(mid$present),
  pairs = length(mid$present)
)

message("mid-sized (p = 50): default route")
t_default = system.time(
  fit_default <- bgm(mid$y,
    variable_type = "continuous",
    chains = CHAINS, cores = CORES, seed = 2026,
    display_progress = "none", verbose = FALSE
  )
)
payload$mid$default = c(
  list(elapsed = unname(t_default[["elapsed"]])),
  recovery(fit_default, mid$present)
)

message("mid-sized (p = 50): fast route")
t_fast = system.time(
  fit_fast <- bgm(mid$y,
    variable_type = "continuous",
    precision_graph_prior = "joint", update_method = "gibbs",
    chains = CHAINS, cores = CORES, seed = 2026,
    display_progress = "none", verbose = FALSE
  )
)
payload$mid$fast = c(
  list(elapsed = unname(t_fast[["elapsed"]])),
  recovery(fit_fast, mid$present)
)

payload$mid$speedup = payload$mid$default$elapsed / payload$mid$fast$elapsed
payload$mid$agreement = sum(
  (edge_probabilities(fit_default) > 0.5) ==
    (edge_probabilities(fit_fast) > 0.5)
)

# The two routes target different models, so the vignette shows what the choice
# does to the prior the analysis is actually run under. Six variables keeps the
# demonstration in seconds: reading the realized prior of a joint fit needs the
# correction table for that model size, and building one is a job of minutes.
message("realized prior inclusion probability (p = 6)")
small = simulate_case(p = 6, n = 300, seed = 1)
fit_small_joint = bgm(small$y,
  variable_type = "continuous",
  precision_graph_prior = "joint", update_method = "gibbs",
  edge_prior = bernoulli_prior(0.5),
  chains = 2, cores = 2, iter = 500, warmup = 500,
  seed = 1, display_progress = "none", verbose = FALSE
)
fit_small_hier = bgm(small$y,
  variable_type = "continuous",
  chains = 2, cores = 2, iter = 500, warmup = 500,
  seed = 1, display_progress = "none", verbose = FALSE
)
payload$prior = list(
  p = 6,
  nominal = 0.5,
  joint = extract_prior_inclusion_probabilities(fit_small_joint)[1, 2],
  hierarchical = extract_prior_inclusion_probabilities(fit_small_hier)[1, 2]
)

rm(fit_default, fit_fast, fit_small_joint, fit_small_hier)
gc()

# ==============================================================================
# 2. How the fast route scales, and the high-dimensional fit
# ==============================================================================
# The 200-variable fit is the vignette's high-dimensional example and also the
# last row of the scaling table, so it is run once and reported twice.
# ==============================================================================
grid = c(25, 50, 100, 150, 200)
scaling = data.frame(
  p = grid, pairs = grid * (grid - 1) / 2,
  edges = NA_integer_, elapsed = NA_real_
)

for(i in seq_along(grid)) {
  p = grid[i]
  message("scaling: p = ", p)
  case = simulate_case(p = p, n = 1000, seed = 2026)
  scaling$edges[i] = sum(case$present)
  tt = system.time(
    fit <- bgm(case$y,
      variable_type = "continuous",
      precision_graph_prior = "joint", update_method = "gibbs",
      chains = CHAINS, cores = CORES, seed = 2026,
      display_progress = "none", verbose = FALSE
    )
  )
  scaling$elapsed[i] = unname(tt[["elapsed"]])

  if(p == 200) {
    payload$high = c(
      list(
        p = p, n = 1000, pairs = scaling$pairs[i],
        edges = scaling$edges[i], elapsed = scaling$elapsed[i]
      ),
      recovery(fit, case$present)
    )
    # How decided the posterior is across all 19,900 pairs. The strongest
    # edges are uninformative to tabulate at this size -- there are hundreds
    # of them and they all sit at 1.000 -- so what the vignette reports is
    # the shape of the whole distribution.
    pip = edge_probabilities(fit)
    cuts = c(0, 0.05, 0.25, 0.75, 0.95, 1)
    payload$high$profile = data.frame(
      inclusion_probability = c(
        "below 0.05", "0.05 to 0.25", "0.25 to 0.75",
        "0.75 to 0.95", "above 0.95"
      ),
      pairs = as.integer(table(cut(pip, cuts, include.lowest = TRUE)))
    )
    # The first six pairs in the fit's own order, as they print.
    payload$high$head = utils::head(
      summary(fit)$indicator[, c("mean", "sd", "Rhat"), drop = FALSE], 6
    )
  }
  rm(fit, case)
  gc()
}
payload$scaling = scaling

# ------------------------------------------------------------------------------
# The payload ships with the package, so it holds the numbers the vignette
# quotes and nothing else: no chains, no draws, no per-edge vectors.
# ------------------------------------------------------------------------------
out = file.path("vignettes", "fast-ggm.rds")
saveRDS(payload, out, version = 2)

size_kb = file.size(out) / 1024
cat(sprintf("wrote %s (%.1f KB)\n", out, size_kb))
if(size_kb >= 100) {
  stop(
    "The vignette payload must stay in double-digit KB; got ",
    round(size_kb, 1), " KB. Thin it further before shipping."
  )
}

str(payload, max.level = 2)

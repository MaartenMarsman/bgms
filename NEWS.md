# bgms 0.2.0.0

## Breaking changes

* `update_method = "hamiltonian-mc"` has been removed. Use `update_method = "nuts"` instead. NUTS dynamically adapts trajectory length and is more reliable, especially with edge selection on GGM models.
* The `hmc_num_leapfrogs` argument has been removed along with pure HMC.

* Pairwise interaction parameters for ordinal MRFs are now stored on association scale (half the sigma scale used in 0.1.6.3). Code that interprets raw pairwise posterior samples or sets `pairwise_scale` explicitly will need adjustment.
* Default `pairwise_scale` changed from 2.5 to 1 to match the association-scale reparameterization.
* `extract_category_thresholds()` is deprecated in favor of `extract_main_effects()`, which covers category thresholds, continuous means, and precision diagonal entries.

## New features

* Gaussian graphical models (GGM): `bgm(x, variable_type = "continuous")` fits a GGM with Bayesian edge selection. Sampling uses NUTS on a free-element Cholesky (theta-space) parameterization of the precision matrix, which keeps the precision matrix positive-definite by construction; adaptive-metropolis is also available.
* GGM Gibbs sampler: `update_method = "gibbs"` fits a GGM with a conjugate row-by-row update of the precision matrix, with or without edge selection. It needs no step-size or proposal tuning and supports a Normal or Cauchy prior on the edges. Continuous data only.
* Gibbs warmup staging: with `update_method = "gibbs"` and edge selection, the first 15% of the warmup runs the full model so the precision matrix settles, and edge selection is active for the remaining 85%. Previously edge selection only started at the first retained iteration, so the graph's equilibration happened inside the retained samples. Both windows scale with the warmup budget; the warmup default is unchanged.
* Standardized-frame diagonal prior: `gamma_prior()` and `exponential_prior()` accept the rate as `eta`, the rate on the precision diagonal in the frame where the pairwise (slab) prior has unit scale. The raw rate is derived at fit time as `eta / s` for interaction-prior scale `s`, so the prior geometry is held fixed when the slab scale changes; at fixed `eta`, graph and partial-correlation inference is invariant to the slab scale. The default `precision_scale_prior` is now `gamma_prior(shape = 1, eta = 1)`, which is identical to the previous `gamma_prior(shape = 1, rate = 1)` at the default unit slab scale; with a non-unit slab scale the diagonal rate now co-scales. Applies to GGM and mixed MRF models. Two consequences of the new default: `sample_ggm_prior()`'s default diagonal rate resolves to `0.4` instead of `1` (its default interaction prior is `cauchy_prior(scale = 2.5)`), and a scale-free interaction prior (`beta_prime_prior()`) on a continuous model now requires specifying the diagonal prior with `rate` instead of `eta`.
* Mixed MRF models: `bgm()` accepts a per-variable `variable_type` vector that mixes `"ordinal"`, `"blume-capel"`, and `"continuous"` types to estimate networks with both discrete and continuous variables. `simulate.bgms()` and `predict.bgms()` also support mixed models.
* Missing data imputation: `na_action = "impute"` integrates over missing values during MCMC sampling for ordinal, continuous, and mixed models.
* `extract_precision()`: extract posterior precision matrix samples from GGM and mixed models.
* `extract_partial_correlations()`: extract posterior partial correlation samples from GGM and mixed models.
* `extract_log_odds()`: extract log-odds for discrete pairwise interactions.
* `extract_main_effects()`: extract main effect samples (category thresholds, continuous means, and precision diagonal).
* NUTS diagnostics now include the per-iteration mean Metropolis acceptance probability (`fit$nuts_diag$accept_prob`, paralleling Stan's `accept_stat__`) and a per-chain `mean_accept_prob` summary.
* `sample_ggm_prior()` accepts `update_method = "gibbs"` for the `spec = "joint"` prior chain, using the conjugate row and edge updates instead of adaptive Metropolis, and `edge_prior = "beta-bernoulli"` for sampling the inclusion probability under its Beta hyperprior.
* Corrected inclusion-probability updates for `beta_bernoulli_prior()` on continuous (GGM) models: under the determinant-tilted precision prior, the conjugate Beta update omits a normalizing-constant factor and biases the sampled inclusion probability toward sparsity (at 5 variables with a uniform hyperprior its prior mean lands near 0.37 instead of 0.5). The update now draws from the corrected conditional using a table built from the prior distribution at the first fit of a model configuration and cached on disk (`tools::R_user_dir("bgms", "cache")`; a one-time cost of the order of minutes, announced when `verbose = TRUE`). The sampled inclusion probability is returned per chain in `fit$inclusion_parameter_samples`.
* Corrected stochastic block model updates for `sbm_prior()` on continuous (GGM) models, using the same cached table: the block-probability draws, the block-allocation weights, and the new-cluster weight all carry the normalizing-constant correction. Without it the sampled partition collapses toward one block (at 20 variables the prior mean number of blocks lands near 1.0 instead of 1.87); with it a prior-only chain reproduces the model's partition prior (total variation 0.01-0.02 at 5-20 variables). `sample_ggm_prior()` accepts `edge_prior` objects (`bernoulli_prior()`, `beta_bernoulli_prior()`, `sbm_prior()`) and returns sampled allocations under the block model.
* The normalizing-constant correction extends to mixed models with `beta_bernoulli_prior()` or `sbm_prior()`: the determinant tilt acts on the continuous precision block, so the correction table is built for the continuous variables and the block-structure corrections read continuous-continuous edges only. With fewer than two continuous variables no edge is tilted and the plain conjugate updates apply unchanged; with exactly two, the single tilted pair supports the inclusion-probability correction but not the block-model slope curve, so `sbm_prior()` warns and keeps the plain conjugate updates there. Prior-only mixed chains reproduce the Beta hyperprior on the inclusion probability and the partition prior on the number of blocks; without the correction the inclusion probability biases toward sparsity and the partition toward fewer blocks.
* `extract_prior_inclusion_probabilities()`: prior edge-inclusion probabilities in the same matrix layout as `extract_posterior_inclusion_probabilities()`, for prior/posterior inclusion-odds computations. Under the joint spike-and-slab prior on a continuous block the graph marginal is reweighted by the per-graph normalizer (positive-definite-cone mass shaped by the determinant tilt), so continuous-continuous edges do not keep the edge-prior marginal at any `delta` — with a uniform Beta-Bernoulli hyperprior at 3 variables the prior edge probability is about 0.37 rather than 0.5, and a fixed `bernoulli_prior(0.5)` at `delta = 0` lands near 0.27. The values are read from the cached correction table (`bernoulli_prior()`, `beta_bernoulli_prior()`) or estimated by a prior-only chain with the fit's own correction settings (`sbm_prior()`, cached on the fit). Mixed models report per-edge-class values (discrete-discrete, continuous-continuous, cross); ordinal models use the analytic edge-prior marginals, including the exchangeable partition mixture for `sbm_prior()`.

## Other changes

* Fitted objects from `bgm()` and `bgmCompare()` are now S7 class objects (new dependency: `S7`). All existing `$`, `[[`, and `names()` access patterns continue to work. When an incompatible `easybgm` version is loaded, bgms returns plain S3 lists for backwards compatibility; this shim will be removed in a future release.
* Refactored the C++ backend: unified model hierarchy (`BaseModel` → `GGMModel` / `OMRFModel` / `MixedMRFModel`), shared NUTS/HMC infrastructure, and fused log-posterior and gradient computation.
* NUTS now uses Stan's multinomial candidate weighting (log-sum-exp of `H0 - h` per leaf, biased progressive sampling at the top level) in place of the Hoffman-Gelman slice variable. The two schemes target the same posterior; the multinomial variant produces lower-variance candidate selection and has been Stan's default since 2017. User-facing output is unchanged apart from the new `accept_prob` diagnostic.
* NUTS Stage-2 warmup windowing now matches Stan's `windowed_adaptation::compute_next_window`: when the window after the next would overshoot the Stage-3a boundary, the current next window is stretched to absorb the remaining Stage-2 budget instead of emitting a small trailing window. This eliminates a disruptive mass-matrix update + step-size reinit at the end of warmup and improves dual-averaging convergence.
* Dropped `coda` from Imports; ESS and R-hat are now computed in C++ with on-demand (lazy) evaluation, replacing the eager R-based computation from 0.1.6.3.
* `$` and `[[` accessors on fitted objects trigger lazy computation of MCMC diagnostics on first access.

## Bug fixes

* Fixed `extract_posterior_inclusion_probabilities()` on mixed MRF fits: the indicator samples are stored in block order (discrete-discrete, continuous-continuous, cross) over internally reordered variables, but the extractor filled the matrix in global pair order, misplacing the entries. It now maps them through the same block filler as `fit$posterior_mean_indicator`.
* Fixed an indicator bookkeeping asymmetry in mixed MRF models: cross (discrete-continuous) edge moves updated only the upper triangle of the internal edge-indicator matrix, while the stochastic block edge prior reads full columns. Under `edge_prior = sbm_prior()`, the block-allocation update for a discrete variable therefore saw its cross edges as always included. Discrete-discrete and continuous-continuous moves, other edge priors, and the reported indicator samples were not affected.
* Fixed the `delta = NULL` default for mixed MRF models: the determinant-tilt exponent applies to the continuous precision block, but the default counted blume-capel variables in its dimension. A mixed model with blume-capel variables now gets `0.5 * log(#continuous)` instead of `0.5 * log(#continuous + #blume-capel)`; models whose discrete variables are all ordinal are unchanged.
* A rank-1 Cholesky downdate that would leave the precision factor non-positive-definite is now detected and triggers a rebuild of the factor from the precision matrix. Previously the failure signal was never checked and a partially updated factor could silently corrupt the carried covariance in long GGM or mixed MRF edge-selection runs.
* Fixed compilation failure on Alpine/musl: `mrf_simulation.cpp` relied on a transitive include for `<tbb/global_control.h>` that is not available on all platforms.
* Fixed stale gradient cache after missing data imputation caused NUTS to use outdated cached values for leapfrog integration.
* Fixed stale observation transpose after missing data imputation caused the pairwise gradient to use stale data.
* Fixed NUTS acceptance probability: target_accept now correctly passed to lower-level NUTS functions.
* Fixed NUTS acceptance probability accumulation: the top-level trajectory loop overwrote the Metropolis contribution with the last subtree's value instead of summing across the full trajectory, biasing the signal used by dual-averaging step-size adaptation.

# bgms 0.1.6.3

## New features

* `extract_rhat()`: extract R-hat convergence diagnostics from fitted objects
* `extract_ess()`: extract effective sample size estimates from fitted objects
* `verbose` argument: control informational messages; set `options(bgms.verbose = FALSE)` to suppress globally
* `simulate_mrf()`: standalone MRF simulation with user-specified parameters
* `simulate.bgms()`: generate observations from fitted models (supports parallel processing)
* `predict.bgms()`: compute conditional probabilities P(X_j | X_{-j})
* `main_difference_selection` argument in `bgmCompare()`: control threshold difference selection
* `standardize` argument: scale Cauchy prior by response score range
* `baseline_category` now stored in fitted object for Blume-Capel simulation/prediction

## Bug fixes

* fixed matrix indexing for `posterior_mean_indicator`: now correctly maps C++ row-major order to R matrices (#77)
* fixed mass matrix adaptation: now correctly uses variance instead of precision
* fixed step size heuristic: re-runs after mass matrix updates, resamples momentum each iteration
* fixed E-BFMI diagnostic: now uses actual accepted trajectory momentum
* fixed Blume-Capel interaction: uses centered scores `(c - ref)` in pseudolikelihood denominator

## Other changes

* NUTS: implemented generalized U-turn criterion following Betancourt (2017) and STAN
* NUTS: fused log-posterior and gradient computation eliminates redundant probability evaluations
* bgmCompare: BLAS-vectorized gradient computation for improved performance
* expanded test suite: input validation, extractor functions, S3 methods, simulation, and numerical sanity tests
* improved warmup schedule: fixed buffers (75/25/50) with proportional fallback for short warmup
* edge selection warmup now within user budget: 85%/10%/5% split for stages 1-3a/3b/3c
* streamlined user messages: concise warnings, consolidated NUTS diagnostics
* E-BFMI threshold adjusted to 0.2 (standard)

## Deprecated

* `mrfSampler()` → use `simulate_mrf()`

# bgms 0.1.6.2

## New features

* added option to separately specify beta priors for the within- and between-cluster probabilities for the SBM prior.

## Other changes

* reparameterized the Blume-capel model to use (score-baseline) instead of score.
* implemented a new way to compute the denominators and probabilities. This made their computation both faster and more stable.
* refactored c++ code for better maintainability.
* removed the prepared_data field from bgm objects.

## Bug fixes

* fixed numerical problems with Blume-Capel variables using HMC and NUTS.
* fixed a reporting bug where category thresholds for ordinal variables with a single category were incorrectly expanded to two parameters, resulting in spurious NA values.

# bgms 0.1.6.1

## Other changes

* added extractor function for joint SBM output
* cleaned up documentation, and c++ files
* changed length of warmup phase I in warmup scheduler HMC / NUTS (15% → 7.5%)

## Bug fixes

* fixed a problem with warmup scheduling for adaptive-metropolis in bgmCompare()
* fixed stability problems with parallel sampling for bgm()
* fixed spurious output errors printing to console after user interrupt. 

# bgms 0.1.6.0

## New features

* added NUTS and HMC options for sampling `bgm()` and `bgmCompare()` models
* added support for running multiple chains in parallel
* added user interrupt handling for parallel sampling
* added Markov chain diagnostics (effective sample size and R-hat) for sampled parameters
* added `summary()`, `print()`, and `coef()` methods for fitted objects
* MCMC sampling in `bgm()` and `bgmCompare()` is now reproducible when a `seed` argument is specified

## Other changes

* improved progress bar for parallel sampling
* `summary()` now integrates the functionality of the old `summary_SBM()`
* removed options for modeling main differences; main differences are now always estimated or selected, equivalent to the previous `main_difference_model = "collapse"` setting

## Bug fixes

* fixed an out-of-bounds error in `bgmCompare()` when handling missing data
* fixed a bug in the SBM prior computation

## Deprecated

- In `bgm()`, the following arguments are deprecated:
  - `interaction_scale` → use `pairwise_scale`
  - `burnin` → use `warmup`
  - `save` → no longer needed (all outputs are returned by default)
  - `threshold_alpha`, `threshold_beta` → use `main_alpha`, `main_beta`

- In `bgmCompare()`, arguments related to difference models are deprecated:
  - `main_difference_model` (removed without replacement)
  - `reference_category` → use `baseline_category`
  - `pairwise_difference_*`, `main_difference_*` → use unified `difference_*` arguments
  - `pairwise_beta_bernoulli_*`, `main_beta_bernoulli_*` → use unified `beta_bernoulli_*` arguments
  - `interaction_scale` → use `pairwise_scale`
  - `threshold_alpha`, `threshold_beta` → use `main_alpha`, `main_beta`
  - `burnin` → use `warmup`
  - `save` → no longer needed

- Deprecated extractor functions:
  - `extract_edge_indicators()` → use `extract_indicators()`
  - `extract_pairwise_thresholds()` → use `extract_category_thresholds()`

- Deprecated object fields:
  - `$gamma` (pre-0.1.4) and `$indicator` (0.1.4–0.1.5) → replaced by `$raw_samples$indicator`
  - `$main_effects` (pre-0.1.4) and `$posterior_mean_main` (0.1.4–0.1.5) → replaced by `$raw_samples$main` (raw samples) and `$posterior_summary_main` (summaries)


# bgms 0.1.5.0 (GitHub only)

## New features

* The bgmCompare function now allows for network comparison for two or more groups.
* The new summary_sbm function can be used to summarize the output from the bgm function with the "Stochastic-Block" prior. 
* Two new data sets are included in the package: ADHD and Boredom.

## Other changes

* The bgm function with the "Stochastic-Block" prior can now also return the sampled allocations and block probabilities, and sample and return the number of blocks.
* The underlying R and c++ functions received a massive update to improve their efficiency and maintainance.
* Repository moved to the Bayesian Graphical Modelling Lab organization.
* Included custom c++ implementations for exp and log on Windows. 

## Bug fixes

* Fixed a bug in the bgmCompare function with selecting group differences of blume-capel parameters. Parameter differences that were not selected and should be fixed to zero were still updated.
* Fixed a bug in the bgmCompare function with handling the samples of blume-capel parameters. Output was not properly stored.
* Fixed a bug in the bgmCompare function with handling threshold estimation when missing categories and main_model = "Free". The sufficient statistics and number of categories were not computed correctly.
* Partially fixed a bug in which the bgms package is slower on Windows than on Linux or MacOS. This is because the computation of exp and log using the gcc compiler for Windows is really slow. With a custom c++ implementation, the speed is now closer to the speed achieved on Linux and MacOS.


# bgms 0.1.4.2

## Bug fixes
* fixed a bug with adjusting the variance of the proposal distributions
* fixed a bug with recoding data under the "collapse" condition

## Other changes
* when `selection = TRUE`, the burnin phase now runs `2 * burnin` iterations instead of `1 * burnin`. This ensures the chain starts with well-calibrated parameter values
* changed the maximum standard deviation of the adaptive proposal from 20 back to 2

# bgms 0.1.4.1

This is a minor release that adds some documentation and output bug fixes.

# bgms 0.1.4

## New features
* Comparing the category threshold and pairwise interaction parameters in two independent samples with bgmCompare().
* The Stochastic Block model is a new prior option for the network structure in bgm().

## Other changes
* Exported extractor functions to extract results from bgm objects in a safe way.
* Changed the maximum standard deviation of the adaptive proposal from 2 to 20.
* Some small bug fixes.

# bgms 0.1.3

## New features
* Added support for Bayesian estimation without edge selection to bgm().
* Added support for simulating data from a (mixed) binary, ordinal, and Blume-Capel MRF to mrfSampler()
* Added support for analyzing (mixed) binary, ordinal, and Blume-Capel variables to bgm()

## User level changes
* Removed support of optimization based functions, mple(), mppe(), and bgm.em()
* Removed support for the Unit-Information prior from bgm()
* Removed support to do non-adaptive Metropolis from bgm()
* Reduced file size when saving raw MCMC samples

# bgms 0.1.2

This is a minor release that adds some bug fixes.

# bgms 0.1.1

This is a minor release adding some new features and fixing some minor bugs.

## New features

* Missing data imputation for the bgm function. See the `na.action` option.
* Prior distributions for the network structure in the bgm function. See the `edge_prior` option.
* Adaptive Metropolis as an alternative to the current random walk Metropolis algorithm in the bgm function. See the `adaptive` option.

## User level changes

* Changed the default specification of the interaction prior from UnitInfo to Cauchy. See the `interaction_prior` option.
* Changed the default threshold hyperparameter specification from 1.0 to 0.5. See the `threshold_alpha` and `threshold_beta` options.
* Analysis output now uses the column names of the data.
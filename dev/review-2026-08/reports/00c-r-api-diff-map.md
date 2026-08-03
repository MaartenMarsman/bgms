# Report 00c — R layer & user-facing API diff map (internal sweep)

Produced 2026-08-01 by the review lead's mapping subagent (read-only sweep, no
brief). Basis: `cran-0.1.6.3` → `v0.2.0.0-rc1`. Findings folded into
FINDINGS.md; architecture sections feed MAINTAINERS.md.

---

R layer: 46 files / 21,988 LOC (was 15 files); `git diff --stat` over `R/` =
51 files changed, +19,539 / −4,585.

## 1. NAMESPACE delta

**23 new `export()`, 0 dropped.** Confirmed by set-difference of `export(`
lines. Every symbol exported at `cran-0.1.6.3` is still exported.

| Feature area | New exports |
|---|---|
| Prior constructors (`R/priors.R`) | `cauchy_prior`, `normal_prior`, `beta_prime_prior`, `gamma_prior`, `exponential_prior`, `bernoulli_prior`, `beta_bernoulli_prior`, `sbm_prior` |
| Prior-sensitivity / verdicts / calibration / plots / centrality | `prior_sensitivity_check`, `verdicts`, `calibration_check`, `plot_edge_posterior`, `extract_centrality` |
| Graph-prior sampling | `sample_graph_prior`, `sample_ggm_prior`, `sample_sbm_prior` |
| zratio / hierarchical-prior gauge | `summarize_zratio_gauge` |
| Extractors | `extract_inclusion_bf`, `extract_main_effects`, `extract_prior_inclusion_probabilities`, `extract_precision`, `extract_partial_correlations`, `extract_log_odds` |

**S3 methods added (35).** By class:

- `bgms` / `bgmCompare` (both): `extract_centrality`, `extract_inclusion_bf`, `extract_main_effects`, `prior_sensitivity_check`, `verdicts`, `$`, `[[`, `names`, `plot`
- `bgms` only: `calibration_check`, `extract_log_odds`, `extract_partial_correlations`, `extract_precision`, `extract_prior_inclusion_probabilities`
- `bgmCompare` only: `extract_sbm.bgmCompare` (new; `extract_sbm.bgms` pre-existed)
- New result classes: `print`/`plot` for `bgms_calibration`; `print`/`plot` for `bgms_prior_sensitivity`; `print.bgms_verdicts`; `summary`/`plot` for `bgms_centrality`
- New spec/prior print methods: `print.bgm_spec`, `print.bgms_parameter_prior`, `print.bgms_scale_prior`, `print.bgms_indicator_prior`

**S3 methods dropped (2):** `extract_category_thresholds.bgms` / `.bgmCompare` —
intentional; `extract_category_thresholds()` is now a deprecation shim
(`R/extractor_functions.R:1037`) forwarding to `extract_main_effects()`.

**Import changes:** `importFrom(coda, ...)` (4 entries) dropped; `importFrom(S7, ...)` (5),
11 more `stats` entries, `importFrom(utils, combn/head/packageVersion)` added.

**Asymmetry:** `calibration_check`, `extract_precision`,
`extract_partial_correlations`, `extract_log_odds`,
`extract_prior_inclusion_probabilities` have no `bgmCompare` method — calling
them on a `bgmCompare` fit gives a bare "no applicable method" error.
`extract_group_params` remains `bgmCompare`-only (pre-existing).

## 2. New user-facing functions

Legend: Rd = has `man/*.Rd`; Test = named in `tests/testthat/*.R`; Vig = referenced in `vignettes/`.

| Function (path:line) | Signature (defaults) | Purpose | Rd | Test | Vig |
|---|---|---|---|---|---|
| `cauchy_prior` `R/priors.R:49` | `(scale = 1)` | Cauchy slab | Y | 11 files | **no** |
| `normal_prior` `R/priors.R:93` | `(scale = 1)` | Normal slab | Y | 12 | prior-sensitivity |
| `beta_prime_prior` `R/priors.R:138` | `(alpha = 0.5, beta = 0.5)` | Beta-prime threshold prior | Y | 4 | **no** |
| `gamma_prior` `R/priors.R:208` | `(shape = 1, rate = NULL, eta = NULL)` | Gamma precision-diagonal prior | Y | 14 | diagnostics |
| `exponential_prior` `R/priors.R:297` | `(rate = NULL, eta = NULL)` | Exponential scale prior | Y | 2 | **no** |
| `bernoulli_prior` `R/priors.R:338` | `(inclusion_probability = 0.5)` | Fixed-p edge prior | Y | 9 | diagnostics |
| `beta_bernoulli_prior` `R/priors.R:375` | `(alpha = 1, beta = 1)` | Learned-p edge prior | Y | 13 | diagnostics |
| `sbm_prior` `R/priors.R:436` | `(alpha=1, beta=1, alpha_between=1, beta_between=1, dirichlet_alpha=1, lambda=1)` | Stochastic-block edge prior | Y | 11 | **no** |
| `prior_sensitivity_check` `R/prior_sensitivity.R:203` | `(bgms_object, anchors = c(0.4,0.63,1,1.6,2.5), evidence_threshold = 10, vary = c("auto","slab","slab-and-diagonal"), refit_sampler = "same-as-fit", iter = NULL, warmup = NULL, tolerance = 0.5, ess_floor = 400, include_preferred_scale = FALSE, cores = NULL, seed = 1L, keep_fits = FALSE, verbose = FALSE)` | Anchored sensitivity curve | Y | 453-line file | 2 |
| `verdicts` `R/verdicts.R:263` | `(bgms_object, evidence_threshold = 10, ...)` | Edge verdicts + MC fragility | Y | Y | 3 |
| `calibration_check` `R/calibration_check.R:277` | `(bgms_object, newdata = NULL, nrep = 200, probs = c(0.025,0.975), grid_size = 101, seed = NULL, ndraws = 500, ...)` | Posterior-predictive calibration | Y | Y | checking-your-model |
| `extract_centrality` `R/centrality.R:68` | `(bgms_object, measure = "strength", group = 1, ...)` | Draws × nodes centrality | Y | Y | checking-your-model |
| `plot_edge_posterior` `R/plot_bgms.R:534` | `(bgms_object, variable1, variable2, evidence_threshold = 10, binwidth = 0.01, ...)` | Single-edge posterior | Y | Y (qgraph-gated) | checking-your-model |
| `sample_graph_prior` `R/sample_graph_prior.R:230` | `(p, n_samples, edge_prior = bernoulli_prior(0.5), spec = c("hierarchical","joint"), interaction_prior = normal_prior(1), precision_scale_prior = exponential_prior(eta=1), delta = NULL, theta = NULL, allocations = NULL, block_probs = NULL, n_warmup = 2e3, seed = 1L, verbose = TRUE)` | Prior-only graph draws | Y | 1 file | **no** |
| `sample_ggm_prior` `R/sample_ggm_prior.R:190` | `(p, n_samples, n_warmup = 2e3, interaction_prior = normal_prior(1), precision_scale_prior = exponential_prior(eta=1), step_size = 0.1, max_depth = 10L, seed = 1L, verbose = TRUE, edge_indicators = NULL, delta = NULL, spec = c("conditional","joint","hierarchical"), edge_inclusion_prob = 0.5, update_method = c("adaptive-metropolis","gibbs"), edge_prior = NULL, apply_correction = TRUE, zratio_diagnostics = TRUE)` | Prior-only GGM chain | Y | 7 files | **no** |
| `sample_sbm_prior` `R/sample_graph_prior.R:379` | `(p, n_samples, edge_prior = sbm_prior(), seed = 1L)` | SBM hyperprior draws | Y | 59-line file | **no** |
| `summarize_zratio_gauge` `R/zratio_gauge.R:154` | `(chains, threshold = 0.01, verbose = TRUE, harm_inputs = NULL, harm_threshold = 0.01)` | Hierarchical trust gauge | Y | Y | mechanism only |
| `extract_inclusion_bf` `R/extractor_functions.R:419` | `(bgms_object, log = FALSE)` | RB inclusion BFs | Y | 2 | **no** |
| `extract_main_effects` `R/extractor_functions.R:932` | `(bgms_object)` | Replaces `extract_category_thresholds` | Y | 2 | **no** |
| `extract_precision` `R/extractor_functions.R:1682` | `(bgms_object)` | Posterior mean K | Y | 4 | intro |
| `extract_partial_correlations` `R/extractor_functions.R:1753` | `(bgms_object)` | Posterior mean partial r | Y | 1 | intro |
| `extract_log_odds` `R/extractor_functions.R:1806` | `(bgms_object)` | Pairwise log-odds | Y | 1 | **no** |
| `extract_prior_inclusion_probabilities` `R/extract_prior_inclusion_probabilities.R:443` | `(bgms_object, iter = 4000L, warmup = 1000L, recompute = FALSE)` | Prior PIPs (may run prior chain) | Y | 1 | **no** |

**Gaps.** All 23 have an Rd. 9/23 have no vignette mention (incl.
`extract_main_effects`, the replacement for a deprecated function). 25 Rd files
have no `\examples{}`, incl. 12 exported extractors. 3 exported Rd lack
`\value{}`: `extract_edge_indicators.Rd`, `extract_pairwise_thresholds.Rd`,
`mrfSampler.Rd` (all `\keyword{internal}`, but `mrfSampler` is a live export).
Thinnest tests: `sample_sbm_prior` (59 L), `extract_log_odds`,
`extract_partial_correlations`, `sample_graph_prior`, `summarize_zratio_gauge`,
`extract_prior_inclusion_probabilities` (1 file each).

## 3. Changes to pre-existing surface

### `bgm()` — `R/bgm.R:475` vs `cran-0.1.6.3:R/bgm.R:383`

Expected five — **all confirmed**: (a) pairwise prior Cauchy→`normal_prior(1)`
(`R/bgm.R:481`); (b) `precision_scale_prior = exponential_prior(eta = 1)` with
rate-XOR-eta validation (`R/priors.R:240`); (c) new `precision_graph_prior`,
default `"joint"` (`R/bgm.R:489`); (d) `update_method` gains `"gibbs"`,
GGM-only (`R/validate_sampler.R:116`); (e) `standardize` reinstated as bare
lifecycle-deprecated formal — FALSE warns, TRUE errors (`R/bgm.R:518/590/600`;
same in `R/bgmCompare.R:237/370/381`).

**FLAGGED — NOT on the expected list:**

1. **`iter`/`warmup` defaults 1e3 → 2e3** (`R/bgm.R:479-480`); same in `bgmCompare()`.
2. **NUTS `target_accept` default 0.60 → 0.80** (`R/validate_sampler.R:130`); bgmCompare 0.65 → 0.80.
3. **`hmc_num_leapfrogs` removed; `"hamiltonian-mc"` dropped** from `update_method` — silent argument removal.
4. **`means_prior = normal_prior(scale = 1)`** — new formal (`R/bgm.R:483`).
5. **`threshold_prior = beta_prime_prior(0.5, 0.5)`** — new formal, supersedes `main_alpha`/`main_beta`.
6. **`delta = NULL`** — new formal, determinant-tilt exponent, auto `0.5*log(p)`.
7. **`edge_prior` type change** character → prior object (character accepted via deprecation branch, `R/bgm.R:614`).
8. **`progress_callback = NULL`** — new formal.
9. Six ex-defaulted arguments now bare deprecated formals (`pairwise_scale`, `main_alpha`, `main_beta`, `inclusion_probability`, `beta_bernoulli_*`, `dirichlet_alpha`, `lambda`).

### `bgmCompare()` — `R/bgmCompare.R:185`

Same removals and 2e3 defaults. Additional: **`difference_family = c("Cauchy","Normal")`** new;
**`interaction_prior = cauchy_prior(1)` — differs from `bgm()`'s `normal_prior(1)`**
(deliberate per NEWS, but identically-named arguments now have different
defaults across the two entry points); `difference_prior` character→object;
`difference_probability` loses its default (lifecycle). `bgmCompare()` has
**no** `precision_graph_prior`, `means_prior`, `precision_scale_prior`,
`delta`, or `"gibbs"` — growing capability gap vs `bgm()`.

### Extractors

- **`extract_ess`** (`R/extractor_functions.R:1539`): `estimator = c("rb","mixt")`;
  default now RB. A 0.1.6.3-era call gets different numbers with no warning;
  pre-0.2.0.0 fits silently fall back to transition ESS.
- **`extract_posterior_inclusion_probabilities`**: `estimator = c("rb","raw")`,
  default `"rb"`; **bgmCompare + rb returns `NA` for unselected indicators** —
  new NA pattern in a previously all-numeric matrix (downstream risk: easybgm).
- **`summary.bgms`** return adds `quadratic` element + `main_label` attribute.

## 4. Priors and spec layer (architecture seed)

- **`R/priors.R` (739 L)** — 8 constructors returning tagged lists with classes
  `bgms_parameter_prior` / `bgms_scale_prior` / `bgms_indicator_prior`. Scale
  priors take `rate` XOR `eta` (`validate_rate_eta()` `:240`);
  `resolve_scale_rate()` (`:618`) maps standardized-frame `eta` → raw rate.
  Five `unpack_*()` helpers (`:557-682`) flatten prior objects into the scalar
  fields `bgm_spec()` expects — the single seam between object API and flat spec.
- **`R/bgm_spec.R` (719 L)** — `bgm_spec()` (`:330`): ~60 flat args in;
  validated list out (`$data $variables $missing $prior $sampler $model_type`).
  Owns hierarchical/joint eligibility + the two user notices
  (`zratio_joint_realized_prior_notice()` `:256`, `zratio_vacuous_spec_notice()` `:308`).
- **`R/build_spec.R` (700 L)** — per-family builders: `build_spec_ggm()` `:53`,
  `build_spec_omrf()` `:124`, `build_spec_mixed_mrf()` `:216`, `build_spec_compare()` `:389`.

**Flow of a `bgm()` call:** deprecation handling + prior unpacking
(`R/bgm.R:655-658`) → `bgm_spec()` (`:662`) → family `build_spec_*` (+
`validate_data/model/sampler.R`) → `run_sampler(spec)` (`R/run_sampler.R:27`)
switching on model type; GGM path resolves joint vs hierarchical and builds
zratio constants + Option-B surface (`R/run_sampler.R:80-98`) or the joint
correction table → `.Call` wrappers in `R/RcppExports.R` (`sample_ggm` :124,
`sample_omrf` :132, `sample_mixed_mrf` :128, `run_bgmCompare_parallel` :4) →
`build_output(spec, raw)` → `build_output_*` → `s3_list_to_bgms()` /
`s3_list_to_bgmCompare()` (`R/class_s7.R`) → S7 fit object.

## 5. Tests

- **77 `test-*.R` files, 24,803 LOC, 1,039 `test_that()` blocks.** Helpers:
  `helper-fixtures.R` 1,396 L (session-scoped `.test_cache`; **fits live models
  at first use** rather than loading RDS), `helper-validation.R` 433 L (mixed
  recovery harness), `helper-internals.R` 40 L (hoists 28 internals outside
  R CMD check), `setup.R` (verbose off; correction-table cache → tempdir).
- Fixtures: `tests/testthat/fixtures/` (340 K, zratio references + legacy/),
  `.Rbuildignore`d; gated by `has_golden_fixtures()`/existence checks.
  **`tests/compliance/` (2.9 MB, 33 RDS; weekly bitwise-vs-CRAN harness) is NOT
  `.Rbuildignore`d — ships in the tarball.**
- Skips: 155 `skip_on_cran`, 43 `skip_if`, 22 `skip_if_not`, 16
  `skip_if_not_installed`, 10 bare `skip()`. ~85% of blocks still run on CRAN;
  28 files have no skips (validators, gradients, spec construction, contracts,
  three zratio files). `skip_on_cran` clusters on SBC, correction-table
  identity, surface builds, recovery, plotting.
  **Risk: `helper-fixtures.R` fits real `bgm()` models unconditionally at first
  access, so CRAN pays for MCMC even in nominally cheap files.**
- **Genuinely thin/absent dedicated tests:** `R/build_output*.R` (4 files),
  `R/build_spec.R`, `R/class_s7.R`, `R/run_sampler.R`, `R/refit_engine.R`,
  `R/mcmc_summary.R`, `R/mcmc_summary_sbm.R`, `R/fit_accessors.R`,
  `R/compute_utils.R`, `R/diagnostics_am.R`, `R/diagnostics_nuts.R`,
  `R/validate_data.R`, `R/validate_model.R`, `R/zratio_surfaces.R`,
  `R/zratio_tables.R`, `R/datasets.R`, `R/zzz.R`.
- `tests/testthat/_problems/` present in working tree (unexplained);
  `test-zratio-surface-build.R` uncommitted-modified vs the tag.

## 6. Vignettes

| File | One-liner |
|---|---|
| `intro.Rmd` | Getting started; **untouched since 2025-07-01**, predates the 0.2 API (no prior constructors, no verdicts, no `precision_graph_prior`; only vignette without bibliography header) |
| `checking-your-model.Rmd` | verdicts / centrality / edge posterior / calibration |
| `diagnostics.Rmd` | R-hat/ESS/MCSE, spike-and-slab summaries, NUTS diagnostics, trust gauge (largest) |
| `prior-sensitivity.Rmd` | sensitivity workflow + tolerance band |
| `comparison.Rmd` | bgmCompare walkthrough (thinnest) |

- **The suspected wrong mixture-ESS gloss is NOT in `vignettes/diagnostics.Rmd`**
  — rewritten by `b42df13f` (RB columns described correctly; `n_eff_mixt`
  gone). **Residual soft gloss at `vignettes/diagnostics.Rmd:88`**: flip counts
  framed as "how much the chain explored the two states" — the weak survivor of
  the rejected reading.
- **Stale text lives on in `doc/diagnostics.Rmd:68-72`** (pre-build artifact;
  `doc/` is `.Rbuildignore`d and holds only 3 of 5 vignettes — a stale local
  build; delete).
- No vignette references removed/deprecated args (`standardize`,
  `pairwise_scale`, `hamiltonian-mc`, `hmc_num_leapfrogs`, `main_alpha`,
  `lambda =`, `extract_category_thresholds`) or states a Cauchy default. Clean.

## 7. DESCRIPTION delta

- Version 0.1.6.3 → 0.2.0.0; Title now "Bayesian Analysis of Graphical
  Models"; Description rewritten. **`Date: 2026-03-26` stale (F-003).**
- Imports: +`S7`, +`graphics`, +`grDevices`, +`parallel`, +`stats`, +`utils`;
  −`coda`. All new Imports genuinely used. Suggests: +`coda`, +`MASS`,
  +`withr` (**verify withr is used or drop**); −`ggplot2`. LinkingTo unchanged.
- **`RoxygenNote: 7.3.3` removed**, replaced by `Config/roxygen2/version: 8.0.0`
  + `Roxygen: list(markdown = TRUE)` — unconventional; Rd generated by a
  roxygen2 ahead of CRAN's line.
- **`SystemRequirements:` absent** while `configure` + `src/Makevars.in` use
  `include sources.mk` (GNU-make construct) and link RcppParallel/TBB — CRAN
  expects `SystemRequirements: GNU make`.
- Root artifacts: `..Rcheck/` present in working tree; `tests/compliance/`
  (2.9 MB) not `.Rbuildignore`d.

## 8. Class system (`R/class_s7.R`, 245 L)

Two S7 classes (`bgms_class` :29, `bgmCompare_class` :149), `package = NULL`,
properties all `class_any` except `.field_names`. **Lazy summaries** via
property getters calling `ensure_summaries(self)` mutating `self@cache` (a
reference-semantics environment) — reading a property has side effects.
**Construction is S3-first**: builders produce a plain S3 list; `s3_list_to_bgms()`
(`:105`) / `s3_list_to_bgmCompare()` (`:213`) convert field-by-field — a
hand-maintained mapping with no test file. All user-facing dispatch is S3 on
the class attribute; `$`, `[[`, `names` are hand-written shims
(`R/methods_bgms.R:303/326/351`, `R/methods_bgmcompare.R:347/370/395`);
`names(fit)` reads the stored `.field_names` vector, not the property set.
Three deprecated easybgm-compat properties (`indicator`, `interactions`,
`thresholds`) alongside canonical ones. Prior/spec objects are S3 lists — the
package carries S3 lists, S3-classed lists, and S7 objects simultaneously.

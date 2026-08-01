# bgms maintainer's guide

Living document, seeded 2026-08-01 from the release-review mapping sweeps
(reports 00a–00c). Audience: a maintainer with statistical expertise and
limited time. Grows as review reports land; the backlog at the end is the
post-release work queue.

## 1. The package in one paragraph

bgms fits Markov random fields to ordinal (`omrf`), continuous (`ggm`), and
mixed data, plus a multi-group comparison model (`bgmCompare`), with Bayesian
edge selection (spike-and-slab indicator variables) and inclusion Bayes
factors as the headline inference. An R layer (~22k LOC, 46 files) validates,
specs, and post-processes; a C++ core (~29k LOC, 87 files, Rcpp/Armadillo/
RcppParallel) runs the samplers. Since 0.1.6.3 the C++ core was effectively
rebuilt: model classes behind an abstract `BaseModel`, a shared chain runner,
multinomial NUTS, and a novel "zratio" subsystem approximating
normalizing-constant ratios for the hierarchical graph prior.

## 2. Architecture map

### Layers and the life of a `bgm()` call

```
user call: bgm() / bgmCompare()                        R/bgm.R, R/bgmCompare.R
  ├─ deprecation shims + prior-object unpacking        R/priors.R (unpack_*)
  ├─ bgm_spec(): flat args -> validated spec list      R/bgm_spec.R
  │    └─ per-family builders + validators             R/build_spec.R, R/validate_*.R
  ├─ run_sampler(spec): family dispatch                R/run_sampler.R
  │    ├─ hierarchical path: zratio constants +        R/zratio_tables.R,
  │    │  Option-B surface build (minutes-order!)      R/zratio_surfaces.R
  │    ├─ joint path: edge-prior correction table      R/correction_tables.R
  │    └─ .Call -> sample_ggm / sample_omrf /          src/sample_*.cpp,
  │       sample_mixed_mrf / run_bgmCompare_parallel   src/bgmCompare_interface.cpp
  │         └─ chain_runner: TBB fan-out over chain    src/mcmc/execution/
  │            clones; warmup schedule; samplers       src/mcmc/samplers|algorithms/
  │              └─ model subclass of BaseModel        src/models/{omrf,ggm,mixed}/
  │                 (bgmCompare: own worker, NOT on    src/models/bgmCompare/
  │                  chain_runner — the one exception)
  ├─ build_output() -> S3 list -> S7 fit object        R/build_output*.R, R/class_s7.R
  └─ post-processing on the fit: extractors, summary,  R/extractor_functions.R,
     verdicts/centrality/calibration/plots,            R/verdicts.R etc.
     prior_sensitivity_check (refits via bgm()!)       R/refit_engine.R, R/anchor_curve.R
```

### C++ subsystem inventory (status vs CRAN 0.1.6.3, risk)

| Subsystem | Status | Risk | One-liner |
|---|---|---|---|
| `src/mcmc/algorithms/` | rewritten | high | Multinomial NUTS (Stan-style, replaced slice NUTS), leapfrog+memoizer, RWM. No isolated test entry. |
| `src/mcmc/execution/` | new | high | chain_runner (TBB chain fan-out, warm starts, post-sampling gauge sweeps), warmup schedule (staged 85/10/5). |
| `src/mcmc/samplers/` | new | high | Sampler classes + dual-averaging and windowed diag-mass adaptation; controllers hold reference members. |
| `src/models/base_model.h` | new | med | Abstract model API; THREE parameter vectorizations (parameter/full/storage) — override consistently. |
| `src/models/omrf/` | rewritten | high | Ordinal/Blume-Capel pseudolikelihood; 3 hand-maintained caches (residual matrix, log-denominators, gradient). |
| `src/models/ggm/` (model+gradient) | new | high | Precision-matrix GGM; incremental Cholesky with drift detector; row-block Gibbs; free-element Cholesky NUTS param. |
| `src/models/ggm/zratio_*` | new | highest | Normalizing-constant-ratio engine: 3 routes (exact isolated, Option-B surfaces, saddle), block-Gibbs oracle, per-chain caches. |
| `src/models/mixed/` | new | high | Discrete+continuous MRF; densest caching (~12 caches with accept/reject adoption); 419-line gradient function. |
| `src/models/bgmCompare/` | heavily modified | high | Own TBB worker (45 args), new CompareSweepState cache + validity mask, thread_local scratch. NO test interface. |
| `src/priors/` | new+modified | high | Polymorphic priors; Beta-Bernoulli/SBM normalizing-constant corrections (grid-interpolated log C, mini-TI for SBM). |
| `src/math/` | new+renamed | med | chol rank-1 update/downdate (PD-loss flag), rank-2 det-lemma, log-sum-exp; custom exp/log unchanged (renamed). |
| `src/utils/` | modified | med | variable_helpers FAST/SAFE exp paths (EXP_BOUND=709); progress manager (atomics, main-thread R callbacks). |
| top-level `src/*.cpp` | new/modified | med-high | Thin entries per family; mrf_simulation (3 TBB workers), mcmc_diagnostics (AR-spectral ESS, split-Rhat, 2-state ESS). |

Unchanged since CRAN (deprioritize in any audit): `custom_explog.cpp` (byte-
identical rename), uncorrected MFM-SBM path (`sbm_edge_prior.cpp:17-407`),
rng engine, `sbm_edge_prior_interface.*`.

### R layer inventory

Fit construction: `bgm.R`, `bgmCompare.R`, `bgm_spec.R`, `build_spec.R`,
`build_arguments.R`, `validate_*.R`, `run_sampler.R`, `build_output*.R`,
`class_s7.R`. Post-processing: `extractor_functions.R` (all `extract_*`),
`mcmc_summary*.R`, `methods_*.R`. Checking layer (new 0.2.0.0): `verdicts.R`,
`centrality.R`, `calibration_check.R`, `plot_bgms.R`, `prior_sensitivity.R` +
`refit_engine.R` + `anchor_curve.R`. Hierarchical support: `zratio_tables.R`,
`zratio_surfaces.R`, `zratio_gauge.R`, `sample_graph_prior.R`,
`extract_prior_inclusion_probabilities.R`, `correction_tables.R`.

## 3. Conventions and invariants (the "laws")

1. **Association scale everywhere.** The linear predictor uses `2·omega·x`;
   pairwise parameters are stored on the association scale (half the 0.1.6.3
   sigma scale). Every path — bgm sampler, bgmCompare sampler, simulate,
   predict, mixed cross terms — must agree. A ~7s every-run cross-
   implementation guard defends this; it exists because bgmCompare silently
   diverged once (the 0.2.0.0 breaking fix). If you touch any likelihood or
   linear predictor: check the factor of 2 first.
2. **Priors travel as objects, run as scalars.** Users pass prior constructor
   objects (`normal_prior(1)`, ...); `unpack_*()` in `R/priors.R` flattens them
   into the spec; C++ rebuilds polymorphic priors from scalars. The seam is
   `unpack_*` — extend all three sides together. Scale priors use the
   standardized `eta` frame with `rate` XOR `eta` (`resolve_scale_rate()`).
3. **Spec is the single source of truth.** Everything the sampler needs is in
   the `bgm_spec` list; C++ interface files only translate. Policy decisions
   (MH targets, zratio deployment routes, gauge sweeps) are resolved in R at
   spec build; `zratio_engine_from_spec` is the one C++ reader, shared by ggm,
   mixed, and the test interface so routes cannot diverge.
4. **Worker threads never touch R.** Inside TBB workers: `std::runtime_error`,
   never `Rcpp::stop`; progress via atomic counters, R callbacks only from the
   main thread. RNG: per-chain xoshiro256++ streams; `metropolis_step_cached`
   exists specifically to preserve the RNG stream — reproducibility depends on
   draw-order discipline.
5. **Caches are adopted, not recomputed.** The models maintain incremental
   state (omrf residual matrix `= 2·X·pairwise`; ggm Cholesky + Sigma rows;
   mixed's ~12 caches; bgmCompare log-normalizer mask). Every accept/reject
   path must adopt or invalidate. Test interfaces pin the key invariants
   (`test_omrf_residual_invariant`, `sigma_rows_consistent_` drift detector,
   `mixed_indicator` both-triangle sync). If you add a move type, wire its
   cache repair AND extend the corresponding invariant test.
6. **Three parameter vectorizations** in `BaseModel` (parameter / full /
   storage dimension). Override all or none in a subclass.
7. **RB estimators are the reported numbers** since #182: PIPs/BFs come from
   Rao-Blackwellized alpha accumulators (log scale for saturation);
   `extract_ess` default is RB n_eff. The raw indicator draws still exist
   (`estimator="raw"/"mixt"`) — keep print/man wording aligned with which
   estimator a number comes from.
8. **Hierarchical path guard rails**: surfaces are validated on an anchor hull
   (cap 80, boundary-slope extrapolation beyond, counters `n_extrap`/
   `n_slope_floor`/`n_collapsed` reported per phase); the trust gauge is
   opt-in and post-hoc; sub-shape-0.5 additive zero-collapse is a documented
   limitation (C2 fix designed, not authorized). Don't relitigate decisions
   recorded in `dev/plans/active/2026-08-01_hier-followup_NOTE.md`.
9. **Never build in the Dropbox tree.** Stale `.o` + generated `Makevars`
   produce unloadable builds; `git archive` to a clean dir (F-009).
10. **Tests fit live models via `helper-fixtures.R`** (session-scoped cache,
    no stale RDS). Costly on CRAN — see F-030 before adding fixture users.
11. **SBC certifies the GGM path only.** The omrf, mixed, and bgmCompare
    paths target pseudolikelihood approximations, so simulation-based
    calibration against simulated data is not a valid check there — their
    correctness gates are the recovery and cross-validation suites. Never
    add an SBC-style test to a pseudolikelihood path, and never read one as
    evidence about it (MM, 2026-08-01).
12. **Anchor budgets: the end-to-end certificates are the authority, not the
    per-cell budget re-derivation (F-052 decision record, MM 2026-08-01).**
    The certification harness re-derives what anchor-budget multiplier each
    cell would need to reach reference precision and disagrees with the
    shipped 2× (`zratio_anchor_shape_multiplier()`) in 5 of 20 cells — it
    wants 4× — and in one cell (eta 2, common-neighbour, shape 0.5) parity
    is unmet even at 4×. ACCEPTED as shipped, because: the disagreement is
    banked (the gold bank was built at the same 2×; predates rc1), all five
    end-to-end route certificates PASS, and the trust gauge polices realized
    deployment harm at fit time — the operative guard in the stubborn cell.
    Do not raise the budgets casually: 4× means rebuilding every surface,
    rebuilding the gold bank, and recertifying, and still does not close the
    one cell. Revisit triggers: a route certificate fails, or the GGM
    paper's post-release calibration work reopens the derivation.

## 4. If you touch X, also check Y

| You touch | Also check |
|---|---|
| Any linear predictor / likelihood | factor-2 convention guard; simulate + predict + mixed cross terms; `test-tolerance`/regression suites |
| `unpack_*` or a prior constructor | `bgm_spec()` validation, C++ `create_parameter_prior`/`create_scale_prior`, `print.bgms_*_prior`, defaults table (downstream memo F-002) |
| `build_output_*` field set | `s3_list_to_bgms()`/`s3_list_to_bgmCompare()` field mapping AND `.field_names` (else the field is invisible to `names()`/`$`) |
| Edge-indicator moves (any model) | RB accumulators (#182), indicator-matrix symmetry tests, `extract_posterior_inclusion_probabilities` |
| zratio engine / surfaces | gold bank (`dev/validation/`), route certificates, extrapolation counters, `zratio_engine_from_spec` (shared ggm/mixed/test), disk-cache version key |
| Warmup schedule / adaptation | warm-start plumbing in `sample_omrf.cpp` (bypasses heuristics), `warmup_schedule_test_interface`, stage-boundary integer arithmetic |
| `extract_*` return shapes | easybgm + JASP (they consume these; see F-002 memo), `$`/`[[` shims, summary methods |
| NUTS internals | `accept_stat__`, dual-averaging targets (0.80 default; 0.44 forced for ggm-with-gibbs policy in `sample_ggm.cpp:85-89`), energy diagnostics |
| Vignettes / man pages | NEWS wording drift (the F-001 class of bug: prose vs computed quantity) |

## 5. Dragons (highest-risk corners)

1. **`zratio_engine.cpp`** — surface evaluation + extrapolation + hand-packed
   cache keys + raw RNG pointer + mutable counters. History of late-found
   defects (deploy gate inert, anchor alpha omission, stale pivot). Any change
   here needs the gold bank rerun + route certificates.
2. **`bgmCompare`** — own parallel worker (not chain_runner), thread_local
   scratch, hand-maintained cache validity mask, largest functions in the
   package, and NO test interface. Validate end-to-end (cross-path consistency
   vs `bgm()` on identical single-group data).
3. **Mixed-model caches** — a dozen incremental caches with adoption-on-accept;
   the densest state machine in the package.
4. **Adaptation controllers** — reference members (`&proposal_sd`, `&schedule`);
   lifetime must outlive the sampler; silent-bias territory.
5. **Edge-prior corrections** — grid-interpolated `log C(θ)` with
   state-dependent inverse-CDF window (a Markov draw, not iid — documented);
   SBM mini-thermodynamic-integration label moves.
6. **`helper-fixtures.R`** — fits real models at first access; innocuous test
   additions can silently add minutes to CRAN check time.
7. **S7/S3 hybrid** — dispatch is S3; S7 supplies storage + lazy getters with
   side effects (`ensure_summaries` mutates `@cache` env); `names()` is a
   stored vector, not the property set.

## 6. Test infrastructure and CI

- 77 test files, ~24.8k LOC, 1,039 `test_that` blocks. Tiers: every-run
  (~85% of blocks, incl. the 7s convention guard — deliberate, see scale-fix
  decisions); `skip_on_cran` for SBC/correction-identity/surface-build/
  recovery/plotting suites; golden fixtures gated by existence
  (`tests/testthat/fixtures/`, `.Rbuildignore`d).
- `tests/compliance/` = weekly bitwise-vs-CRAN-0.1.6.3 harness (2.9 MB RDS),
  run by `.github/workflows/weekly-compliance.yaml`. Must be `.Rbuildignore`d
  (F-017).
- CI: R-CMD-check + lint + test-coverage on push/PR (main, develop);
  nightly-validation Mon+Thu 03:00 UTC (`devtools::test`, slow tests);
  weekly-compliance Sun 05:00 UTC. pkgdown workflow retired (docs live in the
  bgms-docs repo).
- Gold references: `dev/validation/` (tracked in git deliberately — hours of
  compute) — zratio gold bank + route certificates + WP artifacts.
- Test-interface coverage gaps (00b §15): NUTS/adaptation/chain_runner/
  bgmCompare/mixed-kernels have no isolated entries — end-to-end only.

## 7. Build system

`configure` / `configure.win` generate `src/Makevars` (from `Makevars.in`) and
`src/sources.mk` (GNU-make `include` → needs `SystemRequirements: GNU make`,
F-027). Custom exp/log via `explog_macros.h` (`CUSTOM_EXP_LOG` switch; OpenLibM
on Windows). LinkingTo: Rcpp, RcppArmadillo, RcppParallel, dqrng, BH.
Generated `Makevars` in a used working tree goes stale — clean-export builds
only.

## 8. Maintenance backlog (post-release; ranked by risk-reduction per hour)

Reconciled with the July 2026 audit (00a §3): Phases 1–2 of that audit are
done and verified in git; Phase 3 partial; Phase 4 items are in FINDINGS as
release blockers (F-003, F-017, F-018). Remaining + new, ranked:

1. **Institutional memory** (F-016, decided: private): snapshot exists at
   `../bgms-audit-archive/dev-audit-plans-2026-08-01.tar.gz`; re-snapshot as
   the record grows. Optional hardening: push `dev/audit` + `dev/plans` to a
   private remote so a copy lives outside Dropbox.
2. **Test interfaces for the blind spots** (00b §15): bgmCompare first (it has
   none), then NUTS step / adaptation controllers. Each unlocks unit-level
   regression tests where today only end-to-end runs exist.
3. **Port bgmCompare onto BaseModel/chain_runner** (AUD-D2) — removes the
   duplicated parallel infrastructure and the thread_local hazard class.
   Large; schedule as its own project.
4. **CRAN-tier the test suite** (F-030 follow-up): explicit cheap-tier fixture
   (no live MCMC on CRAN), keep statistical suites nightly.
5. **Converter/dispatch hardening**: tests for `s3_list_to_*` + `.field_names`
   (today: hand-maintained, untested — 00c §8); thin/no-test R files
   (`build_output*`, `run_sampler`, `refit_engine`, `mcmc_summary*`).
6. **Docs debt**: intro.Rmd rewrite (F-026); Rd examples for the 12 extractor
   pages; vignette for `sample_*_prior` + `extract_prior_inclusion_
   probabilities` (9 exports have no vignette mention).
7. **Mixed gauge harm channel** (F-022) if not done for release.
8. **Debris sweep** (F-032): dead export attribute, commented-out harness,
   `doc/` stale build, `_problems/`, five absorbed branches (F-014), Readme.Rmd
   decision (F-007).
9. **Deferred perf items** (00a PERF-open): compare-imputation deltas, shared
   clones, O(q³) extraction floor + Phase-4 cross-sweep cache (SPH-2), opt-in
   tiers #23–27. Post-release; measured gains already banked.
10. **SLAB-2 conditional-means diagnostic** (F-020) if deferred at release —
    the known blind spot with a designed detector.

## 9. Reference points

- Frozen review target: tag `v0.2.0.0-rc1`; CRAN anchor: tag `cran-0.1.6.3`.
- Review record: `dev/review-2026-08/` (FINDINGS.md = master list).
- Decision records: `dev/audit/*.md`, `dev/plans/active/*_NOTE.md`,
  `plans/flagged-issues.md` (cross-repo tracker, outside this repo).
- Downstream consumers to notify on surface changes: easybgm, JASP, the
  bgms-docs site, the tutorial repo (defaults memo = F-002).

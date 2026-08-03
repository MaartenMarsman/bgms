# Report 00b — C++ core diff map (internal sweep)

Produced 2026-08-01 by the review lead's mapping subagent (read-only sweep).
Basis: `cran-0.1.6.3` → `v0.2.0.0-rc1`. Findings folded into FINDINGS.md;
subsystem sections feed MAINTAINERS.md.

**Scale:** src/ goes from ~11.4k tracked LOC (47 files) to ~28.9k LOC (87
files). Deleted wholesale: `src/bgm/` (2752 LOC), `src/bgm_interface.cpp`
(391), the flat `src/mcmc/*.{cpp,h}` layout, `src/utils/print_mutex.h`.
`src/RcppExports.cpp` 211 → 1029 LOC; 49 registered exports.

**Cross-cutting semantic change:** the factor-of-2 rest-score convention
(association scale) appears in both `src/mrf_prediction.cpp` (~:110/:140) and
`src/models/bgmCompare/bgmCompare_sampler.cpp:~100`. Same-signature R-visible
functions return different numbers than at CRAN (documented in NEWS Breaking
changes).

---

## 1. Top-level interface files

- **`sample_ggm.cpp` — 172 LOC — NEW.** Thin R entry: builds priors, constructs
  `GGMModel`, wires determinant tilt `delta`, the Z-ratio engine from spec,
  missing-data indices, edge prior + correction, calls `run_mcmc_sampler`.
  Risk **medium**: the single place where sampler-type→MH-target policy
  (`:85-89`, gibbs excluded, nuts forced to 0.44) and the zratio spec resolve.
  Review `:85-123`.
- **`sample_mixed.cpp` — 211 LOC — NEW.** Same shape for `MixedMRFModel`; four
  independently-defaulted prior families from an untyped R list; MH-target
  hardcode at `:132`. Risk medium.
- **`sample_omrf.cpp` — 172 LOC — NEW** (replaces `bgm_interface.cpp`).
  Adds warm start: per-chain `initial_parameters` / `initial_step_sizes` /
  `initial_inv_mass` (`:150-159`) bypassing step-size heuristic and windowed
  mass adaptation. Risk medium.
- **`bgmCompare_interface.cpp` — 487 LOC — substantially modified.**
  RcppParallel `GibbsCompareChainRunner` worker (`:94`, 189 lines) fanning
  chains over TBB; per-chain deep copies; `tbb::global_control` `:448`.
  **Risk high** — the only sampler path NOT migrated to
  `mcmc/execution/chain_runner`; 45 R arguments. Review `:94-283`, `:440-455`.
- **`zratio_interface.cpp` — 62 LOC — NEW.** `zratio_block_oracle_moments`:
  block-Gibbs oracle → weighted moments (S1, S2), called at fit setup from
  `R/zratio_surfaces.R` to generate Option-B surface anchors. **Risk high** —
  its output is the anchor data the deployed surfaces are fit to; an error
  propagates silently into every hierarchical fit.
- **`mrf_prediction.cpp` — 292 LOC (+200/−2).** Conditional predictive
  distributions (GGM `:24`, ordinal/BC `:72`, mixed `:181`). Risk medium (two
  functions new; third changed by ×2 convention).
- **`mrf_simulation.cpp` — 1147 LOC (+743/−66).** Gibbs simulators + three
  RcppParallel workers over posterior draws (`:324`, `:538`, `:915`), each with
  `tbb::global_control`. **Risk high.** Review `:699-846` `simulate_mixed_mrf`.
- **`mcmc_diagnostics.cpp` — 454 LOC — NEW.** ESS via AR spectral density
  (Levinson-Durbin + AIC, matching `coda::spectrum0.ar`), split-Rhat, two-state
  Markov ESS for indicators. Risk medium. Review `:29-242`, `:313-410`
  (`IndicatorESSWorker`, τ_int documented `:296-306`).
- `explog_interface.cpp` (28, moved), `sbm_edge_prior_interface.cpp/.h`
  (33/24, header rename only) — low risk.

## 2. `src/math` — 907 LOC (was 508)

- `custom_explog.cpp` (386) — **100% identical rename** of `e_exp.cpp`.
  Deprioritize.
- `explog_macros.h` (64, renamed from `explog_switch.h`) — **switch predicate
  changed meaning** (`:28-29`): explicit `-DCUSTOM_EXP_LOG=1` now wins on every
  platform (at baseline that combination disabled OpenLibM on Windows). Adds
  `MY_LOG1P`. 6 functional lines — worth a 5-minute read.
- `cholupdate.cpp`/`.h` (136/49) — **NEW**: Givens rank-1 update + hyperbolic
  downdate (ported from mgcv `mat.c`). Downdate can lose positive-definiteness
  (success flag exists for this). Review `chol_up` (~90 lines from `:30`).
- `cholesky_helpers.h` (140) — **NEW**: `get_log_det`,
  `compute_inv_submatrix_i`, rank-2 matrix-determinant-lemma log-ratio
  (`:48-95`) used in every GGM/mixed MH ratio.
- `log_sum_exp.h` (21) — NEW.

Risk medium overall. **Debris:** `cholupdate.cpp:123` carries a live
`// [[Rcpp::export]]` on `chol_update_arma` never scanned by
compileAttributes (subdirectory) — dead attribute, dead function weight.

## 3. `src/mcmc/algorithms` — 937 LOC — rewritten

- `nuts.cpp` (399) / `nuts.h` (61) — **the NUTS algorithm changed families**:
  slice-sampling NUTS (Hoffman–Gelman) → **multinomial NUTS** with log-sum-exp
  progressive/biased sampling and the Betancourt generalized U-turn criterion,
  following Stan's `base_nuts.hpp` (documented `nuts.cpp:15-26`).
- `leapfrog.cpp/.h` (81/145) — leapfrog + single-entry `Memoizer`.
- `hamiltonian_utils.cpp/.h` (82/75) — kinetic energy, heuristic initial step.
- `metropolis.cpp/.h` (45/49) — scalar RWM + `metropolis_step_cached`
  (RNG-stream-preserving).
- Deleted: `mcmc_hmc.{cpp,h}` — HMC update method gone.

**Risk high** — whole-algorithm replacement producing every posterior draw.
Review: `nuts.cpp:55-256` `build_tree` (201 lines), `:266-399` `nuts_step`,
`hamiltonian_utils.cpp` `heuristic_initial_step_size`. No isolated test entry
(end-to-end only).

## 4. `src/mcmc/execution` — 1227 LOC — NEW/rewritten

- `chain_runner.cpp` (439) / `.h` (147) — `run_mcmc_chain` (147 lines, `:68`),
  `MCMCChainRunner::operator()` (`:217`), `run_mcmc_sampler` (`:243`),
  `convert_results_to_list` (`:341`); TBB fan-out `:323-324`; post-sampling
  gauge-sweep block `:182-205`; `std::runtime_error`-not-`Rcpp::stop`
  discipline on worker threads (`:40-41`).
- `chain_result.h` (308) — sample storage + diagnostics incl. 15 zratio/gauge
  fields (`:99-126`).
- `warmup_schedule.h` (195) — stage 1/2/3a/3b/3c schedule; 85/10/5 split under
  selection + `learn_sd`; stage 3b skipped under 20 iterations; Gibbs settle
  fraction 0.15. Integer stage-boundary arithmetic; review `:34-195`.
- `sampler_config.h` (52), `step_result.h` (86).

**Risk high** (parallel chain clones, warm-start injection, post-sampling
gauge sweeps). No test interface for chain_runner/chain_result/conversion.

## 5. `src/mcmc/samplers` — 707 LOC — NEW/rewritten

- `sampler_base.h` (81) — abstract sampler (incl. warm step-size/inv-mass hooks).
- `nuts_sampler.h` (194), `metropolis_sampler.h` (66), `gibbs_sampler.h` (45).
- `nuts_adaptation.h` (257) — `DualAveraging` (`:17`),
  `DiagMassMatrixAccumulator` (`:87`; Welford + weak prior `prior_weight=5.0`,
  `prior_variance=1e-3` at `:113-114`), `NUTSAdaptationController` (`:139`).
- `metropolis_adaptation.h` (64) — Robbins-Monro proposal-SD controller;
  **holds `arma::mat& proposal_sd` as a reference member** (`:19`).

**Risk high** — adaptation is the classic source of silent bias; both
controllers store references whose lifetimes must outlive the sampler. Only
`WarmupSchedule` has a test interface; DualAveraging/mass-matrix/controllers
have none. Review `nuts_adaptation.h:139-257` (window boundaries, restarts).

## 6. `src/models/base_model.h` — 414 LOC — NEW

Pure-virtual model interface consumed by the whole mcmc/ stack; zratio
trust-gauge hooks (`:186-209`); `enum class ZRatioPhase` (`:17`). Risk medium:
**three distinct parameter vectorizations** (`parameter_dimension` /
`full_parameter_dimension` / `storage_dimension`) with defaulted fallbacks — a
subclass overriding one but not another silently mis-stores.

## 7. `src/models/ggm` — model + gradient — 3219 LOC — NEW

`ggm_model.h` (1000) / `ggm_model.cpp` (1323); `ggm_gradient.h/.cpp` (225/539);
`graph_constraint_structure.h` (132). GGM on the precision matrix; element-wise
adaptive Metropolis with incremental Cholesky update/downdate; exact row-block
Gibbs (Normal/Cauchy slab + Gamma diagonal only, gated `ggm_model.cpp:462`);
joint indicator+parameter edge MH with optional conjugate variant; NUTS via a
free-element Cholesky parameterization with per-column null-space constraints
and reverse-mode adjoint through stored Givens rotations.

**Risk high**: incremental Cholesky/covariance caching across accept/reject
(`sigma_rows_consistent_` drift detector `:1158`), rank-2 SMW update (`:372`),
mutable theta cache, four MH ratios that must each pick up the `delta·log|K|`
tilt. Review: `ggm_model.cpp:750-902` `update_edge_indicator_parameter_pair`;
`:510-653` `update_row_block_gibbs`; `ggm_gradient.cpp:163-538` (forward map,
logp/grad, adjoint).

## 8. `src/models/ggm` — zratio engine/law/gauge — 2444 LOC — NEW

`zratio_engine.cpp` (1190) / `.h` (593); `zratio_law.h` (501, **dormant**);
`zratio_gauge.h` (160). Computes per-edge normalizing-constant ratio
`J = Z(Γ−)/Z(Γ+)` for the hierarchical prior `p(K|Γ) = ρ_Γ(K)/Z(Γ)`. Three
routes by spec: (a) exact isolated-edge `log ψ0`; (b) Option-B absolute-moment
surfaces (bivariate quadratic in (log size, density), 9 monomials, fit offline
to block-Gibbs anchors; hull clamping + boundary-slope extrapolation); (c)
additive-counts two-moment saddle over cosine-transform tables. Block-Gibbs
oracle (Schur + rank-2 SMW + slice sampling on the pivot) supplies gold
references. Persistent caches: `cache_` (count key), `surf_cache_` (component
multiset), `comp_cache_`.

**Risk high — highest in the package.** Novel numerics; stateful per-chain
caches deep-copied on clone; raw `SafeRNG*` member (`zratio_engine.h:515`);
8 `mutable` counters mutated from `const` methods; hand-packed 64-bit cache
keys (`pack_count_key` `:454`, shifts 42/21); hand-rolled union-find
(`:533-544`). Review: `zratio_engine.cpp:408-542` `surface_logr_`;
`:193-238` `surface_eval_` + `:57-66` `saddle_ratio`; `:780-931` `gibbs_sweep_`
+ `:748-778` `smw_rank2_col_update_` + `:30-55` `slice_pivot_` (caps
`kSliceOutCap=20`, `kSliceShrinkCap=50`, `kPivotFloor=1e-12`).

`zratio_gauge.h`: in-chain trust gauge (D, per-pair SEs, MCSE, noise floor);
reached only via `chain_runner.cpp:182-205`; diagnostic-only (medium).

## 9. `src/models/mixed` — 3789 LOC — NEW

`mixed_mrf_model.cpp/.h` (1261/797), `mixed_mrf_metropolis.cpp` (978),
`mixed_mrf_gradient.cpp` (608), `mixed_mrf_likelihoods.cpp` (145). Joint MRF
over discrete + continuous: dd block, cross block, continuous precision Kyy
with tilt, means. **Densest caching layer in the package** (~a dozen caches
with accept/reject adoption via `adopt_kyy_proposal_caches` /
`adopt_cross_proposal_caches`; mutable scratch `mixed_mrf_model.h:478-482`,
`:510-528`); three edge-indicator move types keeping both triangles in sync.

**Risk high.** Review: `mixed_mrf_gradient.cpp:188-608` `logp_and_gradient`
(419 lines — largest function in src/); `mixed_mrf_metropolis.cpp:238-386`
(`log_ggm_ratio_edge`/`_diag`); `:698-977` (three `update_edge_indicator_*`).

## 10. `src/models/omrf` — 1746 LOC — rewrite of deleted `src/bgm/`

`omrf_model.cpp/.h` (1237/509). Git detects **no** rename similarity with the
old `bgm_sampler.cpp` — a rewrite into the BaseModel shape, not a port.
Ordinal/Blume-Capel pseudolikelihood; incrementally maintained
`residual_matrix_ = 2·X·pairwise` (`:302`); per-variable log-denominator cache;
gradient cache; AM moves; missing-data imputation.

**Risk high**: three interacting caches maintained by hand across every move;
`impute_missing` (`:1093`, 108 lines) mutates observations and must repair all
three. Review: `omrf_model.cpp:697-841` `logp_and_gradient`; `:925-985`
`update_edge_indicator`; `:597-644` log-denominator path.

## 11. `src/models/bgmCompare` — 4662 LOC (was 4287)

Moved from `src/bgmCompare/`. Non-comment changed lines:
`bgmCompare_sampler.cpp` 912, `bgmCompare_logp_and_grad.cpp` 433, helper 55.
New `bgmCompare_state.h` (85): `CompareSweepState` — per-group obs/pairwise/
residual + **per-variable × per-group log-normalizer cache with a validity
mask**. Cached MH steps (`:210`, `:288`); incremental rank-1 `pairwise_stats`
update during imputation; the ×2 rest-score change; `thread_local`
`LogZAndProbs`/`LogZScratch` at `bgmCompare_logp_and_grad.cpp:604-605`.

**Risk high**: hand-maintained validity mask; `thread_local` state inside a
TBB worker; largest functions in the package; **no test interface at all**.
Review: `bgmCompare_sampler.cpp:1612-1950` `run_gibbs_sampler_bgmCompare`;
`:1124-1353` `update_indicator_differences_metropolis_bgmcompare`;
`bgmCompare_logp_and_grad.cpp:275-763`.

## 12. `src/priors` — 1830 LOC (was 430)

- `parameter_prior.h` (203) — NEW polymorphic Cauchy/Normal/BetaPrime/GammaScale.
- `edge_prior.h` (361) — NEW Bernoulli / BetaBernoulli (correction hook `:102`)
  / StochasticBlock (`:158`); `attach_edge_prior_correction` (`:304`).
- `edge_prior_correction.h` (190) — NEW. Piecewise-linear interpolation of
  `log C(θ)` with linear extension (`:42-60`); **inverse-CDF draw on a
  401-point grid** centered at current θ, half-width 10 posterior SD
  (`:70-103`) — window depends on current state, so the draw is a Markov
  chain, not i.i.d. (documented at `edge_prior_correction_test_interface.cpp:38-40`).
- `sbm_edge_prior.cpp` (889, +575) / `.h` (187) — baseline MFM-SBM collapsed
  Gibbs retained (lines 17–407 ≈ CRAN); **new corrected branch** from `:408`:
  `compute_ce_sbm`, `degrees_ld_sbm`, mini thermodynamic integration
  (`miniti_node_sbm` `:491`, `miniti_removal_sbm` `:552`),
  `corrected_log_marginal_mfm_sbm` (`:613`),
  `block_allocations_mfm_sbm_corrected` (`:690`; **grows/shrinks `block_probs`
  in place**, `.h:79`), `draw_theta_local_density` (`:787`).

**Risk high.** Review: `edge_prior_correction.h:42-103`;
`sbm_edge_prior.cpp:491-611`; `:690-824`.

## 13. `src/rng` — `rng_utils.h` 218 LOC — lightly touched

Only functional addition: `rgamma` (boost gamma, rate→scale inversion), used by
row-block Gibbs and SBM corrected draws. Engine (xoshiro256++) unchanged. Low.

## 14. `src/utils` — 1471 LOC

- `variable_helpers.h/.cpp` (197/568) — split from a 598-line header. Ordinal +
  BC denominators/probabilities with FAST/SAFE branching at `EXP_BOUND = 709`;
  `*_into` scratch-buffer variants (`:339`, `:449`) commented as bit-identical
  to allocating versions. Risk medium (FAST power-chain path `:7-163` is the
  delicate part; raw `memptr()` scanning). No direct test entry.
- `progress_manager.cpp/.h` (450/177) — R callback support, `std::atomic`
  per-chain counters with explicit memory orders, main-thread-id gating,
  `poll()`. Risk medium. Debris: 45-line commented-out example harness
  (`:406-450`); the sole TODO in src/ (`:40`, cosmetic).
- `common_helpers.h` (79) — `hamiltonian_mc` enum removal only. Low.

## 15. Test interfaces (all 11 NEW)

| File | Exposes |
|---|---|
| `cholupdate_test_interface.cpp` | rank-1 downdate + PD flag |
| `edge_prior_correction_test_interface.cpp` | `log_C` interpolator; corrected θ draw as n-step chain |
| `ggm_gibbs_test_interface.cpp` | row-block Gibbs sweep in isolation |
| `ggm_gradient_interface.cpp` | GGM logp/gradient, forward map; `sample_ggm_prior_cpp` |
| `mixed_gradient_interface.cpp` | mixed logp_and_gradient in theta space |
| `mixed_indicator_test_interface.cpp` | indicator matrix both-triangle sync after n sweeps |
| `omrf_residual_test_interface.cpp` | `residual_matrix_ == 2·X·pairwise` invariant across warmup |
| `prior_test_interface.cpp` | parameter/scale prior logp/grad |
| `sbm_correction_test_interface.cpp` | `compute_ce_sbm`, mini-TI, corrected log-marginal |
| `warmup_schedule_test_interface.cpp` | stage boundaries + `selection_enabled` |
| `zratio_test_interface.cpp` | 9 entries: eval/reference/saddle/surface/spec/batch/gold/precompute/law |

**No test interface:** `mcmc/algorithms` (NUTS/leapfrog/metropolis),
`mcmc/samplers` (DualAveraging, mass matrix, both adaptation controllers),
`mcmc/execution` (chain_runner/chain_result/conversion), **`models/bgmCompare`
(none at all)**, `models/omrf` beyond the residual invariant, `models/mixed`
Metropolis kernels + cache adoption, `mrf_simulation` / `mrf_prediction` /
`mcmc_diagnostics` internals, `utils/variable_helpers` FAST/SAFE branches.

## 16. TODO/FIXME/HACK

Exactly **one** marker in all of src/: `progress_manager.cpp:40` (cosmetic).

## 17. Dead / transitional debris

1. `zratio_law.h` (501 LOC) — dormant by design ("DORMANT large-q insurance —
   NOT wired into the default surface build", `:12-21`); included only by the
   test interface; compiled into the CRAN binary.
2. `math/cholupdate.cpp:123` dead `// [[Rcpp::export]]` (`chol_update_arma`).
3. `progress_manager.cpp:406-450` commented-out RcppParallel example harness.
4. Working tree: 44 gitignored `.o` files; generated `src/Makevars` with a
   hard-coded absolute RcppParallel lib path (template `Makevars.in` is
   correct); verify tarball excludes.
5. Baseline debris removed: `src/RcppExports-9e54f986.o.tmp` (tracked at CRAN)
   is gone. Good.
6. No orphan headers.

## 18. Unchanged since CRAN (deprioritize)

`math/custom_explog.cpp` (identical rename), `custom_explog.h` /
`custom_arma_explog.h` (rename + doxygen), `sbm_edge_prior_interface.*`
(one include line), `sbm_edge_prior.cpp:17-407` (uncorrected MFM-SBM path),
`rng_utils.h` (except `rgamma`), `common_helpers.h` (except enum removal),
`explog_macros.h` (except 6 lines). Everything else is new, rewritten, or
substantially modified.

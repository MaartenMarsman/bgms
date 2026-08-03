# Brief 24 — hierarchical default, gauge harm wiring, the sweeps program, and the mixed degenerate guard (Opus agent)

You are the execution agent in the pre-CRAN review of the R package **bgms**.
You share no context with the review lead; everything you need is here. CODE
CHANGES authorized on a fix branch. Four maintainer decisions land together:
F-010 (GGM `precision_graph_prior` defaults to `"hierarchical"`), F-022
(mixed fits must report harm, not NA), F-103 (sweeps-first program, approved
verbatim "OK"), F-123 (degenerate-block guard, option (b)).

## Machine budget (standing rule — TWO other agents may be running)

~4 hardware threads; ONE fit at a time; heavy work last (task 3's
measurement, ~48 sequential fixture fits, is the heavy tail). Runtime is not
a grading criterion.

## Setup

- Repo (Dropbox; do NOT switch its branch or build in it):
  `/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms`
- Base on **`origin/develop`** (NOT the local ref) AT OR AFTER `1aefd2ff`:
  ```sh
  git -C /Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms \
      worktree add ~/bgms-review/wt-fix15 -b fix/hierarchical-gauge origin/develop
  cd ~/bgms-review/wt-fix15
  ```
- One commit per task; `fix(ggm):`/`fix(mixed):`/`test:`/`docs:` prefixes
  with the finding ID. No pushing; no attribution trailers.
- **HARD off-limits (two unmerged batches own them):**
  `R/methods_bgmcompare.R`, `tests/testthat/test-summary-marking.R`,
  `man/summary.bgmCompare.Rd`, `R/validate_data.R`, `R/bgmCompare.R`,
  `man/bgmCompare.Rd`, and every `tests/testthat/test-*ompare*.R` file.
  Nothing in your four tasks needs any of them. Also as always: `NEWS.md`
  (verbatim proposal in the report; the lead lands it), plot files,
  `vignettes/` (report-only list, task 1), `dev/review-2026-08/` except your
  report + assets.

## Tasks — 1, 2, 4 are small and sequential; 3 is the program

### 1. The default flip (F-010, maintainer-decided)

`bgm()` declares `precision_graph_prior = c("joint", "hierarchical")`
([R/bgm.R:488](R/bgm.R#L488)); `match.arg` makes `"joint"` the default.
Reorder to make **`"hierarchical"` the default**; `"joint"` stays fully
available. Then:
- Sweep roxygen and Rd for every place the default is named or implied
  (the `@param` at [R/bgm.R:170](R/bgm.R#L170) and anywhere else); update.
- `vignettes/` are off-limits: LIST in the report every vignette/pkgdown
  spot that names or assumes the joint default (the Phase-3 docs re-bake
  needs that list); do not edit them.
- Tests: expectations pinned to joint-default behaviour either (a) pass
  `precision_graph_prior = "joint"` explicitly if the test's INTENT is the
  joint path, or (b) re-derive under hierarchical if the intent is default
  behaviour — say which, per touched expectation. Snapshot re-records only
  with a stated reason (default-fit prints may gain gauge lines).
- Confirm `extract_arguments()` on a default GGM fit reports hierarchical.
- NEWS (report-only, VERBATIM proposal): the GGM never shipped, so this is
  FEATURE wording for the 0.1.6.3 reader — if the existing GGM feature
  bullet names the default or the prior composition, propose the adjusted
  bullet; never "changed from joint".

### 2. F-022 — mixed hierarchical fits must report harm, not NA

`R/build_output_mixed_mrf.R:305` calls
`summarize_zratio_gauge(zratio_chains, verbose = TRUE)` WITHOUT
`harm_inputs`, so `harm_pred` is permanently NA on mixed fits — and after
task 1 the gauge is part of the default mixed experience. The GGM builder
does it right: mirror [R/build_output_bgm.R:355](R/build_output_bgm.R#L355)
(the assembler for `harm_inputs` is documented at
[R/zratio_gauge.R:572](R/zratio_gauge.R#L572)). If the mixed path genuinely
cannot supply some ingredient, STOP on this task and report the limitation
plus a proposed Rd sentence instead of forcing it. Test: a seconds-scale
mixed hierarchical fit yields finite, non-NA `harm_pred` in its gauge
summary.

### 3. F-103 — the sweeps program (approved: sweeps-first)

The shipped harm alarm false-fires on ~25% of healthy negative-control fits
(report 13 §9b): `harm_pred` at the shipped configuration spans
0.0004–0.0125 across 12 seeds while `harm_threshold = 0.01` sits inside the
spread. `n_samples` does not shrink it (measured). The lever is the gauge's
own audit precision: `sample_ggm_prior()` hardwires
`gauge_sweeps = if(isTRUE(zratio_diagnostics)) 2L else 0L`
([R/sample_ggm_prior.R:330](R/sample_ggm_prior.R#L330)) while the deployed
path resolves `options(bgms.zratio_gauge_sweeps)` via
`zratio_gauge_sweeps()` ([R/run_sampler.R:90](R/run_sampler.R#L90), :242).
**`harm_threshold = 0.01` is a product constant — do NOT touch it.**

a. **Un-hardwire.** `sample_ggm_prior()`'s `2L` → the same
   `zratio_gauge_sweeps()` resolution the deployed path uses (preserving
   the `zratio_diagnostics`-off ⇒ 0 behaviour).
b. **Measure.** `harm_pred` for the negative control
   `biased_evidence_free_fit(2)`
   ([tests/testthat/test-zratio-gauge.R:204](tests/testthat/test-zratio-gauge.R#L204),
   the F-103 block at :243) across the SAME 12 seeds report 13 used, at
   `gauge_sweeps` ∈ {2, 4, 8, 16}, shipped `n_samples` (1200). One fit at a
   time. The sweeps-2 row must reproduce report 13's spread (sanity
   anchor). Report per sweep count: mean / sd / max / max÷threshold /
   number exceeding threshold / per-fit wall time.
c. **Choose the default.** The smallest sweep count whose 12-seed max sits
   below the threshold with REAL margin — target max ≤ 0.005 (half the
   threshold); state the achieved margin. **If even 16 sweeps leaves the
   max straddling the threshold, STOP: report the table, change no
   default — the threshold/margin question returns to the maintainer.**
d. **Set it.** The chosen count becomes `zratio_gauge_sweeps()`'s option
   default, with the measured per-fit runtime cost vs 2 sweeps stated in
   the roxygen that documents the option and in the report.
e. **Restore F-103's block to T1.** Remove the `skip_unless_certification`
   parking from the :243 block, re-tune its expectation with derivation,
   and PROVE the restored assertion passes 12/12 seeds at the new default
   (run exactly that loop; it is the point of the program).

### 4. F-123 (b) — mixed degenerate-block guard + loop-bound hygiene

Maintainer-decided (option (a) was rejected — do NOT make degenerate blocks
work). Two small changes:
- **Guard:** at the top of `build_spec_mixed_mrf()`
  ([R/build_spec.R:236](R/build_spec.R#L236) region), `stop()` when the
  discrete or continuous block is empty, message stating the mixed model
  requires at least one discrete and one continuous variable (pure-type
  data routes to the OMRF/GGM — this guard is defensive; no exported path
  reaches it). Unit test on the internal (both degenerate directions).
- **Hygiene:** totalize the seven unguarded upper-triangle bounds —
  `mixed_mrf_model.cpp` :851 (`p_`), :858 (`q_`) and
  `mixed_mrf_gradient.cpp` :31/:87/:146/:498 (`p_`), :570 (`q_`) — to the
  `i + 1 < n` form the same file already uses in the RB mirrors (:884,
  :898). Identical behaviour for all p, q ≥ 1: PROVE with one same-seed
  seconds-scale mixed fit before/after the change (identical draws).
- No NEWS (the mixed model never shipped).

## Verification gate

1. Both suite tiers: 0 failures / 0 warnings; every re-tuned expectation
   listed with its derivation and its (a)/(b) classification from task 1.
2. `R CMD check --as-cran` on a `git archive` tarball: 2 baseline NOTEs.
3. `devtools::document()` clean; NAMESPACE unchanged.
4. Task-3 table complete with runtimes; F-103 block restored to T1 and
   proven 12/12 at the chosen default (or the STOP branch taken and
   documented).
5. Task-4 same-seed identity check passes.
6. `extract_arguments()` default-fit check from task 1 passes.

## Deliverable

`/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms/dev/review-2026-08/reports/24-hierarchical-gauge.md`:
What was done (per task, commits) / task-3 sweeps table + chosen default +
margin + runtime cost / task-2 before/after gauge summary / task-1 touched
expectations with derivations + the vignette/docs-site re-bake list /
proposed NEWS wording VERBATIM / findings (severity-tagged) / open
questions. Copy to the Dropbox path; commit the report on the branch.

# Brief 20 — standalone batch: compliance harness, sparse-coding proof, extractor gap (Opus agent)

You are the execution agent in the pre-CRAN review of the R package **bgms**.
You share no context with the review lead; everything you need is here. CODE
CHANGES authorized on a fix branch. Three independent items, chosen because
they touch NOTHING that the three in-flight batches own.

## Machine budget (standing rule — THREE other agents may be running)

Cap your footprint at ~4 hardware threads; everything sequential; the only
non-trivial compute is one compliance-suite run and a handful of
seconds-scale fits. Runtime is not a grading criterion.

## Setup

- Repo (Dropbox; do NOT switch its branch or build in it):
  `/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms`
- Base on **`origin/develop`** (NOT the local ref) AT OR AFTER `8cab3c9b`:
  ```sh
  git -C /Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms \
      worktree add ~/bgms-review/wt-fix10 -b fix/standalone-batch origin/develop
  cd ~/bgms-review/wt-fix10
  ```
- One commit per item; `fix:`/`test:`/`docs:` prefixes with F-numbers. No
  pushing; no attribution trailers. Build/test in the worktree or exports.
- **HARD off-limits (other batches own them):** `NEWS.md`, `vignettes/`,
  `R/validate_data.R`, anything under `src/models/bgmCompare/`,
  `man/bgmCompare.Rd` and its roxygen source, every plot file
  (`R/plot_*.R`, `R/calibration_check.R`, `R/centrality.R`,
  `R/prior_sensitivity.R`), `R/verdicts.R`. Where an item below earns a
  NEWS clause, write the clause VERBATIM in your report — the lead lands
  NEWS at integration.

## Tasks

### 1. F-100 — the weekly-compliance harness is stale (major)

The scheduled compliance run on `main` fails with
`Error: Can't find property <bgms>@posterior_coclustering_matrix` — the
harness was never updated after the property rename; the current S7 class
carries `posterior_mean_coclustering_matrix`
(`tests/compliance/test_compliance.R:428,495`,
`tests/compliance/generate_fixtures.R:507`; the passing name is used in
`tests/testthat/test-bgmCompare.R:302`). S7's `@`/`$` errors on unknown
properties, so one stale name kills the run.

a. Sweep EVERY property name the compliance harness touches (both files,
   and any helper they source) against the current S7 class definitions —
   not just the one known offender. List every skew you find.
b. Fix the names; regenerate whatever fixtures the harness needs
   (`generate_fixtures.R` — seeded, sequential); run the compliance suite
   once against your tree and report its verdict line.
c. FORWARD-COMPATIBILITY FLAG, do not fix: another in-flight batch will
   add stored fields to compare fit objects (a category recode map and a
   per-group support table in the arguments). If the harness asserts an
   EXACT property/field set anywhere (would break on additions), name the
   assertion and line in your report so the lead can coordinate at that
   batch's merge — do not pre-change it.

### 2. F-109 — prove the sparse-coding fix, both ends (minor, evidence for a NEWS claim)

Context: at `cran-0.1.6.3`, `recode_data_for_prediction()`
(`R/simulate_predict.R:1238-1252` at the tag) recodes ordinal `newdata` by
MIN-SHIFT (`x - min(x)` when `min > 0`), so a SPARSE original coding — gaps,
e.g. values 1,2,4,5 — miscodes silently: the fit's own recode collapses to
0..3 while min-shift maps to 0,1,3,4. Current develop is believed to use a
stored-map recode. Your job is to turn "believed" into evidence:

a. CURRENT side: on your tree, fit a small ordinal model (seconds-scale:
   few variables, moderate n, seeded) on sparse-coded data (e.g. categories
   1,2,4,5), then (i) `simulate()` → `predict()` round trip and
   (ii) `predict()` on original-scale sparse `newdata` — verify both land
   on the correct internal categories (state HOW current code recodes:
   file:lines of the stored-map mechanism).
b. TAG side: export and install `cran-0.1.6.3` into a temp library
   (`git archive cran-0.1.6.3 | ...`; R CMD INSTALL to a scratch lib —
   minutes, sequential), run the SAME sparse-data calls, and show the
   miscode actually happening (the wrong category assignment, printed).
c. Add a regression test on the current behaviour (sparse coding round
   trip + original-scale predict), in the test file that owns
   simulate/predict behaviour — NOT in any compare test file.
d. In your report: the verbatim proposed NEWS Bug-fixes entry (0.1.6.3
   reader; the tag-side demo is its evidence). Do NOT edit NEWS.md.

### 3. F-114 — `extract_arguments()` has no `main_effect_indices` for bgmCompare fits (minor)

Surfaced when a diagnostic script errored on it (report 16 §7.iv). Read how
the field is populated for `bgm()` fits and what the compare fit stores
instead; then EITHER provide the field for compare fits (if the information
exists in the object — mirror the bgm shape) OR, if providing it is not a
contained change, document its absence in the extractor's Rd and make the
accessor fail with a helpful message instead of a bare NULL/error. State
which branch you took and why in the report; add the matching test.

## Verification gate

1. Full LOCAL default-tier suite: 0 failures / 0 warnings. (No src/
   changes expected; if item 3 forces one, run the slow tier too.)
2. The compliance suite runs GREEN against your tree (item 1b verdict
   line quoted).
3. The tag-side miscode demo output quoted in the report (item 2b).
4. `devtools::document()` clean; no NAMESPACE drift.

## Deliverable

`/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms/dev/review-2026-08/reports/20-standalone-batch.md`:
What was done (per item, commits) / Findings (severity-tagged — especially
any further property skews and the forward-compat assertions) / Evidence
(compliance verdict line, both sparse-coding runs, before/after for
item 3) / Proposed NEWS clause VERBATIM / Open questions.

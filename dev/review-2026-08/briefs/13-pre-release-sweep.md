# Brief 13 — pre-release sweep: plot follow-ups, guard closures, CI reds, vignette bake (Opus agent)

You are the execution agent in the pre-CRAN review of the R package **bgms**.
You share no context with the review lead; everything you need is here. CODE
CHANGES authorized on a fix branch. This is the last planned fix batch before
release mechanics: many small items, three open-ended diagnoses. Work the
tasks IN ORDER — cheap and closed first, the CI diagnoses last.

## Machine budget (standing rule — the release-gate batch is ALSO running)

Cap your footprint at ~4 hardware threads. ONE model fit at a time,
everything sequential. The 2-core CI reproductions (task 9) use `cores = 2`
by their nature; the tag-vs-rc1 fits (task 8) run one at a time. Runtime is
not a grading criterion — never parallelize to save wall-clock.

## Setup

- Repo (Dropbox; do NOT switch its branch or build in it):
  `/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms`
- Base on **`origin/develop`** (NOT the local ref) AT OR AFTER `1327c497`:
  ```sh
  git -C /Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms \
      worktree add ~/bgms-review/wt-fix11 -b fix/pre-release-sweep origin/develop
  cd ~/bgms-review/wt-fix11
  ```
- One commit per task; `fix:`/`test:`/`docs:`/`ci:` prefixes with F-numbers.
  No pushing; no attribution trailers. Build/test in the worktree or exports.
- **HARD off-limits (the release-gate batch owns them):** `R/validate_data.R`,
  anything under `src/models/bgmCompare/`, `man/bgmCompare.Rd` and its roxygen
  source, `vignettes/comparison.Rmd`. Also off-limits as always: `NEWS.md`
  (propose clauses VERBATIM in your report; the lead lands NEWS), and
  `dev/review-2026-08/` except your own report. `sv/*` and the `bgms-docs`
  repo are READ-ONLY.
- The other batch may merge into develop while you work. Before writing your
  report: `git fetch`, merge `origin/develop` into your branch, resolve, and
  re-run every test file you touched. Say in the report whether this happened.

## Tasks

### 1. F-105 — remove the `display_log_bf()` workaround (minor)

`format_log_bf()` (`R/verdicts.R`) now carries the −0.0 guard (merged, brief
17). The plot-side stopgap is therefore redundant: delete `display_log_bf()`
(`R/plot_style.R:419-430` and its header comment) and route its two call
sites (`R/plot_bgms.R:1288`, `:1352`) to `format_log_bf()`. The rendered
figures must be PIXEL-IDENTICAL — verify via the plot snapshots (no
re-records expected) and say so. If any pixel changes, STOP on this task and
report why instead of re-recording.

### 2. F-107 — remove `plot_edge_posterior(evidence_threshold =)` (minor)

The argument is inert (it only drove the deleted verdict word; the panel's
numbers don't depend on it) and `plot_edge_posterior()` never shipped —
0.1.6.3 exports no plot functions — so it is removed OUTRIGHT, no
deprecation cycle (maintainer-ratified never-shipped rule). Remove it from
the signature, the internal `verdicts()` row lookup (call at the default
threshold — the row's contents are threshold-free), the roxygen, and any
tests/examples passing it. No NEWS entry (never shipped). Check first that
nothing else began using the argument since the plot merge; if something
did, report it rather than forcing the removal.

### 3. F-115 — the default-device squeeze: documentation remedy only (minor)

The three-panel evidence display crowds at R's default 7×7 device
(pre-existing qgraph geometry; the review's renders use 13.5×5.2). Decided
remedies — DOCUMENTATION ONLY, no geometry changes (a coordinate-box fix was
already tried and declined; do not guess at qgraph internals):
a. A recommended-device-size sentence in `man/plot.bgms.Rd` and
   `man/plot.bgmCompare.Rd` (roxygen source): the three-panel display wants
   a wide device, with a concrete suggestion (e.g. `width = 13, height = 5`).
b. A one-line `bgms.verbose`-gated `message()` when the three-panel path
   draws on a device narrower than ~10 inches, naming the Rd advice. Follow
   the existing `bgms.verbose` message pattern; no message on the
   single-panel or groups displays.

### 4. F-066 residue — extend the slab-frame identity pin to GGM (test)

Report 10 proved the prior overlay needs no change of variable (slab frame =
extractor frame on every family) but the regression fixture pins this
ordinal-only. Extend it: a seconds-scale GGM fit asserting the extractor
frame equals the slab frame (the −0.5 double-application trap the review
caught before it shipped). Add mixed too if it is equally cheap; say which
you added. Place it with the existing frame-identity test.

### 5. F-019 — the GGM `n_ == 0` positive-definiteness question (major)

The July audit's last unresolved correctness item: `ggm_edge_move` and
`ggm_diag_move` guard their conjugate updates (`src/models/ggm/ggm_model.cpp:288`,
`:682`), but `update_edge_indicator_conjugate` (`:995`) has no visible
`n_ == 0` positive-definiteness guard. Settle it:
a. Targeted read: determine whether the unguarded path is reachable with
   `n_ == 0` (prior-only sampling, empty-data corners) and whether it can
   propose a non-PD state. Write the argument down either way.
b. Empirical: a long `sample_ggm_prior(update_method = "gibbs")` run
   (bounded — think minutes, not hours; state iterations and seed) watching
   for PD failures/NaNs.
c. If a guard is genuinely missing, add it MIRRORING the existing two (same
   shape, same tolerance), with a regression test. `src/models/ggm/` is in
   scope for this batch. If the path is unreachable or safe, close the
   finding with the written argument instead — do not add a guard for show.
d. If you fixed something real: does the defect exist at `cran-0.1.6.3`?
   Check, and if so propose the NEWS Bug-fixes clause (0.1.6.3 reader) in
   your report.

### 6. F-042 — delete the five dead golden-fixture tests (minor)

`test-simulate-predict-regression.R:225-332` (five blocks) resolve fixtures
under `dev/fixtures/scaffolding/`, which no longer exists and has no
generator in the repo — they can never run anywhere and report as skips
("golden fixtures not found"), which reads as deferred coverage when it is
dead code. Maintainer-ratified: DELETE the five blocks entirely (no stubs,
no skip placeholders). Sections 7 and 8 of that file (original-scale and
sparse-coding tests) are live — do not touch them. No NEWS (never
functional).

### 7. F-073 — close the four compare convention-guard holes (major, test coverage)

The every-run compare guard (`tests/testthat/test-bgmCompare.R:337-368`) has
lead-verified holes. Close all four:
a. A planted-difference scale pin: two groups with a KNOWN nonzero δ on
   named pairs (`difference_selection = FALSE`; reuse the guard's own
   construction; ~5 s budget), asserting the recovered difference matches
   the planted magnitude within a stated tolerance — this is the assertion
   a ×2 error in the difference parameterisation alone would fail, and
   today NOTHING in the every-run tier would catch it.
b. Make the existing guard two-sided: the current bounds pass a ×½
   divergence (`0.5·rmse(2·target)` only tests the doubling direction).
   Bound the error from both sides.
c. Assert group 2, not just `[, 1]`.
d. A numeric-convention test for `simulate.bgmCompare()` (its only
   group-difference test self-describes as soft): seeded simulated margins
   against `predict()` expectations at modest `nsim`, tolerance derived
   from the binomial MC error, in the file that owns simulate/predict
   behaviour.
State every tolerance's derivation in comments. These tests must FAIL on a
deliberately broken build — spot-prove (a) by injecting a ×2 locally,
running the test, reverting; quote the failure line in the report.

### 8. Tag-vs-rc1 ordinal estimate agreement (report-only validation)

The compliance harness compares structure only (fixtures are
machine-specific — accepted, F-116). What no artifact yet shows is ESTIMATE
agreement between the release candidate and CRAN 0.1.6.3 on the plain
ordinal model. Produce it:
a. Install `cran-0.1.6.3` into a scratch lib (`git archive cran-0.1.6.3 |
   tar -x` into a temp dir, `R CMD INSTALL` — the recipe report 20 §3.3
   used; reuse `~/bgms-review/lib10-tag` if it is still there).
b. Five seeded synthetic ordinal datasets (modest: e.g. p = 6, n = 500,
   4 categories, full support — category collapse must NOT be in play).
   Fit both versions on identical data with matched samplers/iterations,
   ONE fit at a time.
c. Compare posterior mean main effects and pairwise associations and the
   BF-10 presence verdicts. Tolerance: refit one seed twice on ONE version
   and use that run-to-run spread as the yardstick; agreement must sit
   within it. Table in the report: per seed, max |Δ| by block, verdict
   agreement count, and the yardstick.
No test, no shipped artifact — this is a one-off validation for the record.

### 9. The CI reds — diagnose IN ORDER, one at a time (major)

Shared context: all three surfaced when brief 12's tiering made blocks
actually run on CI. Local runs on this 15-core machine are green; the CI
runner has 2 cores and Linux BLAS. Documented mechanism to keep in mind:
oneTBB ≥ 2022 makes chains beyond the first follow different trajectories
across core counts (statistically equivalent, not bitwise).

a. **F-080 + F-081 together** (T2 blocks: `test-mixed-nuts.R` M.2F,
   `test-sbc-ggm.R` diagonal ranks). Reproduce locally at `cores = 2` with
   the T2 env set. If they fail at 2 cores: decide whether the assertion
   tests the property or the trajectory. A tolerance re-founding must be
   STATISTICALLY argued (e.g. from the assertion's own sampling
   distribution across seeds/core-counts), not widened until green — show
   the derivation. If they pass at 2 cores locally, say so and dig one
   level (BLAS? runner image?) before concluding; do not close as
   "unreproducible" without stating what WAS ruled out.
b. **F-103** (`test-zratio-gauge.R:243`, `harm_pred` 0.01071 ≥ 0.01000 on
   Linux only). Diagnose whether `harm_pred` drifts by platform (RNG
   path/BLAS in the assessment sweeps) or the fixture/threshold is
   mis-founded. Fix so the block can RETURN TO T1 — it is a named T1
   concern parked in T2 against the tier contract; restoring it is part of
   the fix. If the diagnosis turns into a harm-threshold POLICY question (a
   product constant, not a test constant), STOP that sub-task and flag for
   the maintainer.
c. **F-104** (`test-zratio-surface-build.R` socket block, T0, one observed
   error in four runs). First-class clue: the erroring run reported **4
   PSOCK worker nodes for a `cores = 2L` call** — `cores` did not reach the
   cluster constructor that time. Read the construction path end to end and
   find every way the worker count can diverge from the argument (option
   defaults, env vars, a racing default to `detectCores()`), then harden:
   plumb `cores` explicitly and assert `length(cl)` equals it inside the
   builder. Keep brief 12's visible-skip probe. If the original error
   cannot be reproduced, harden anyway and say plainly that the hardening
   is preventive, not proven causal.

### 10. Bake the prior-sensitivity vignette (build cost)

Precedent (maintainer): earlier vignettes ran analyses locally and shipped
the material as `.rds`. Apply it to `vignettes/prior-sensitivity.Rmd`:
a. First determine how the `ps` object's output currently gets into the
   built vignette (the F-086 fix made the shown report REAL — find whether
   it is a live chunk or pasted output).
b. Convert to the bake pattern: a TRACKED generator script produces the
   `ps` object once, locally, seeded; the vignette chunk shows the call
   un-evaluated and loads the `.rds`, printing the REAL object (F-086's
   honesty must survive — no pasted text).
c. Rails (maintainer's, hard): the `.rds` holds the ps object ONLY;
   double-digit KB at most (strip/thin if needed and say how); the
   generator script is tracked in the repo (and `.Rbuildignore`d if it
   would otherwise ship); `R CMD check --as-cran` stays at the 2 baseline
   NOTEs; vignette build time drops and you report before/after seconds.

### 11. OPTIONAL, only if it stays out of off-limits files — F-077

`arguments$baseline_category` stores raw (unshifted) values for ordinal
entries on mixed compare fits — inert today (`src/mrf_prediction.cpp:99-105`
reads it only in the Blume–Capel branch) but a trap for reference
implementations. Candidate: store `NA` for non-BC entries. Do this ONLY if
the write site is outside `R/validate_data.R` and `src/models/bgmCompare/`
(both off-limits — the release-gate batch owns them). Check the C++ read
sites tolerate `NA` first. If the write site is off-limits, write the
recommendation in the report and stop.

## Verification gate

1. Full LOCAL suite, default tier AND `BGMS_RUN_SLOW_TESTS=true` (src/
   changes expected from tasks 5/9): 0 failures / 0 warnings each. List
   every expectation you changed with its reason.
2. `R CMD check --as-cran` on a `git archive` tarball: the 2 baseline NOTEs
   only (the vignette task makes this non-optional).
3. Plot snapshots: NO re-records (task 1 is pixel-identical; tasks 2/3
   touch no drawn pixel). Any needed re-record = STOP and explain.
4. `devtools::document()` clean; NAMESPACE unchanged.
5. The task-7 planted-δ test proven to bite (the injected-×2 failure line
   quoted).
6. Post-fetch merge freshness confirmed (see Setup).

## Deliverable

`/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms/dev/review-2026-08/reports/13-pre-release-sweep.md`:
What was done (per task, commits) / Findings (severity-tagged — especially
the three CI diagnoses: mechanism, evidence, what was ruled out) / Evidence
(gate outputs; the tag-vs-rc1 table; before/after vignette build seconds;
the planted-δ failure line) / Proposed NEWS clauses VERBATIM (task 5d if
any) / Open questions.

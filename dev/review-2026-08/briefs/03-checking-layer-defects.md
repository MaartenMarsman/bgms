# Brief 03 — Checking-layer defect batch: reproduce, root-cause, fix (Opus agent)

You are the execution agent in the pre-CRAN review of the R package **bgms**.
You share no context with the review lead; everything you need is here. This
brief authorizes CODE CHANGES on a fix branch — the first brief that does.

## Context

MM's hands-on evaluation (report `dev/review-2026-08/reports/02-user-facing-checks.md`,
readable on branch `develop`) surfaced defects in the new user-facing checking
layer. The review lead cross-checked them against the code; root-cause
hypotheses below. Your job: reproduce, confirm or correct each hypothesis, fix,
add regression tests, and report.

- Repo (Dropbox, currently checked out on `main` — do NOT switch its branch):
  `/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms`
- Work in a git worktree OUTSIDE Dropbox so builds are clean and the user's
  checkout is untouched:
  ```sh
  git -C /Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms \
      worktree add ~/bgms-review/wt-fix -b fix/checking-layer-batch develop
  cd ~/bgms-review/wt-fix   # build/install/test from here only
  ```
- One branch, one commit per finding (`fix(verdicts): ...` style, reference the
  F-number). Do NOT push to GitHub; leave the branch local — the review lead
  merges after verifying your report. Do not add any attribution trailers to
  commits.
- House rules: no unrelated refactors; match surrounding style; every fix gets
  a regression test that fails before and passes after; preserve the two layer
  conventions (association scale ×2; RB estimators are the reported numbers).

## The findings (fix in this order)

### F-035 (major) — `calibration_check()` crashes on Blume-Capel and mixed fits

Repro:
```r
fit = bgm(Wenchuan, variable_type = "blume-capel", baseline_category = 1,
          chains = 2, seed = 1)
calibration_check(fit)
# Error in x[, j] <- levels_list[[j]][x[, j] + 1L] :
#   number of items to replace is not a multiple of replacement length
```
Also crashes on a mixed fit (6 continuous / 6 ordinal / 5 blume-capel columns).

Hypothesis (lead's code read): `fitted_observed_data()`'s `decode`
(`R/calibration_check.R:48-49`) indexes `arguments$category_levels` per column;
for blume-capel columns that entry is likely NULL or mis-keyed
(`R/calibration_check.R:363-364` builds `levels_list` and its names differently
for mixed) — `NULL[idx]` has length 0, which produces exactly this error.
Confirm, fix the level bookkeeping for BC (and mixed name alignment), and add
regression tests: calibration_check on (a) pure BC fit, (b) mixed fit with BC
columns, (c) ordinal fit with non-contiguous category scores if supported.

### F-036 (major) — `verdicts()` silently runs a 5000-iteration prior chain on hierarchical fits; interrupt still "returns results"

Repro: `bgm(Wenchuan, variable_type = "continuous",
precision_graph_prior = "hierarchical")`, then `verdicts(fit)` → minutes of
silence. MM interrupted at ~2 min and verdicts STILL printed a table.

Mechanism (verified by the lead): `extract_inclusion_bf()`
(`R/extractor_functions.R:466`) calls `extract_prior_inclusion_probabilities()`,
which for hierarchical specs runs `prior_only_chain_pips(spec, iter = 4000L,
warmup = 1000L)` (`R/extract_prior_inclusion_probabilities.R:134,299`).

Your tasks, in order:
1. **Correctness first**: determine what MM's interrupted call actually
   returned. Find the interrupt/tryCatch path — can an aborted prior chain
   yield partial prior PIPs (and therefore silently wrong Bayes factors and
   verdicts)? Test empirically: interrupt-simulating wrapper or reduced-iter
   comparison. If partial results are possible, make interruption ERROR
   cleanly instead of degrading.
2. Check the caching story: `recompute = FALSE` suggests the prior PIPs are
   cached on first computation — confirm where (fit object? env?), whether
   MM's second call was served from cache (explaining "then it does show
   results!"), and whether an interrupt can poison the cache.
3. UX: emit an upfront message ("computing prior inclusion probabilities for
   the hierarchical prior via a prior-only chain (~N iterations); one-time,
   cached") + progress; document the `iter`/`warmup` control in
   `?extract_prior_inclusion_probabilities` and reference it from `?verdicts`
   and `?extract_inclusion_bf`.
4. Regression tests: hierarchical fixture fit → first verdicts() call messages
   and caches; cached second call is instant; interrupt behavior per your
   fix.

### F-037 (major) — mixed-fit inclusion BFs anomalously small: DIAGNOSE, do not fix without sign-off

On the mixed fit above, MM found the edge 1–2 BF "very small" while the same
edge is decisively large on all-ordinal, all-BC, and all-continuous fits of the
same data. Protocol:
1. Controlled comparison: same columns, (a) `variable_type = "continuous"`
   (GGM) vs (b) a mixed fit with as many columns as possible declared
   continuous (if the mixed path requires a discrete column, use one). The
   shared continuous-block edges must have closely agreeing PIPs/BFs. Seeds
   fixed, chains = 4, defaults otherwise.
2. If they disagree: bisect the BF pipeline — raw indicator means vs RB PIPs
   vs prior-odds term (`extract_prior_inclusion_probabilities` vs
   `difference_prior_inclusion` paths), and the ×2 association-scale
   convention in the mixed cross terms (`src/models/mixed/`). Identify the
   first divergent quantity.
3. Report the numbers. If the root cause is an unambiguous mechanical bug
   (wrong index, wrong prior-odds source), fix it with a test. If it is
   anything with statistical judgment in it, STOP and write it up for MM.

### F-038 (minor) — `plot_edge_posterior()` title prints a 131-digit number

`R/plot_bgms.R` (title block in `plot_edge_posterior`): `bayes_factor = row$bf`
then `sprintf("%.1f", bayes_factor)` — a decisive edge (log10 BF ≈ 130) prints
the full non-scientific integer. Fix: display log10 BF (consistent with
`verdicts()` output), or switch to scientific notation past |log10 BF| ≈ 4.
Match whatever `print.bgms_verdicts` does after F-039. Test: title string for
a saturated-edge fixture stays short and parseable.

### F-039 (minor) — verdicts print must state the log10 base and boundaries

Classification is CORRECT (`build_verdicts`: `lthr = log10(evidence_threshold)`,
boundaries ±1 for threshold 10) but the printed header names only "10 / 0.1"
while the column is `log10_bf` — MM himself misread −1.474 as a natural log and
filed it as a classification bug. Amend the header to state the boundaries in
the displayed unit, e.g. "presence: log10 BF > 1; absence: log10 BF < −1
(threshold 10)". Snapshot-test the print.

### F-040 (minor) — calibration output UX

(a) Label `share_outside_band` units (proportion vs percent) in print and Rd.
(b) `plot.bgms_calibration`: add a `variables =` (and/or paging) argument so a
17-variable fit need not draw 17 panels at once. Tests for both.

### F-041 (note) — guard empty reductions in the sensitivity warm-start checks

`prior_sensitivity_check()` on a degenerate source fit leaks
`min()/max(): no non-missing arguments` warnings from reductions over
`wc$ebfmi_*` / `wc$var_ratio` (`R/prior_sensitivity.R`, warm-start convergence
checks). Guard the all-NA/empty case explicitly. While there: confirm the
refit-iteration floor (MM's iter=10 fit produced 1500-iteration refits) is
intended and documented; report, don't change.

### F-047 quick win — friendlier same-variable error

`plot_edge_posterior(fit, 1, 1)` errors correctly but tersely; make the message
name the resolved variable and suggest the fix. One-liner.

## Verification gate (all must pass before you write the report)

From the worktree, with a clean install of the branch:
1. The new/changed test files pass.
2. The full CRAN-mode suite passes (it is only ~70 s: run `testthat::test_dir`
   without `NOT_CRAN`).
3. `NOT_CRAN=true` tier for the touched areas passes (`test-verdicts.R`,
   `test-calibration-check.R`, `test-plot-methods.R`, `test-prior-sensitivity.R`,
   `test-prior-inclusion-probabilities.R`, `test-extractor-functions.R`).
4. `R CMD check --as-cran` on the branch tarball reports no new NOTES/WARNINGS
   vs report 01's baseline (2 NOTEs: stale Date, HTML tidy).

## Deliverable

Write `/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms/dev/review-2026-08/reports/03-checking-layer-defects.md`
(the path exists on develop; if the Dropbox checkout is on main, write the file
anyway — it is untracked there and the lead collects it):

1. **What was done** — per finding: reproduced? hypothesis confirmed/corrected;
   the fix (files, approach); commit SHA on `fix/checking-layer-batch`.
2. **Findings** — anything NEW you discovered, severity-tagged; explicitly
   answer the F-036 interrupt-correctness question and the F-037 verdict
   (bug vs legitimate, with numbers).
3. **Evidence** — before/after output for each repro, verification-gate
   results verbatim.
4. **Open questions** — anything needing MM (especially F-037 if statistical).

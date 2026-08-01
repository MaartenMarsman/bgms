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

**MM's directive (2026-08-01): the 2-minute chain is unexpected AND UNNEEDED.**
So the primary fix is not messaging — it is eliminating the chain where a
cheaper correct route exists. Your tasks, in order:
1. **Correctness first**: determine what MM's interrupted call actually
   returned. Find the interrupt/tryCatch path — can an aborted prior chain
   yield partial prior PIPs (and therefore silently wrong Bayes factors and
   verdicts)? Test empirically: interrupt-simulating wrapper or reduced-iter
   comparison. If partial results are possible, make interruption ERROR
   cleanly instead of degrading.
2. **Find the cheaper route.** PR #193 built machinery that computes/reports
   the REALIZED edge prior under the joint specification (see
   `tests/testthat/test-joint-realized-prior-notice.R` and the related
   `R/zratio_*`/spec code, plus `dev/plans/active/2026-07-31_hier-gating_NOTE.md`
   WP "realized-prior notice"). Establish whether the hierarchical spec's
   prior inclusion probabilities are obtainable from that machinery (or
   another closed/cheap form) instead of a 5000-iteration prior-only chain.
   Write up the mathematical claim explicitly (what quantity the realized-
   prior machinery yields, and why it equals — or does not equal — the prior
   PIP the BF needs); MM confirms validity before you switch the default
   route. If no cheap route survives scrutiny, fall back to: upfront message
   + progress + documented `iter`/`warmup` control.
3. Check the caching story either way: `recompute = FALSE` suggests the prior
   PIPs are cached on first computation — confirm where (fit object? env?),
   whether MM's second call was served from cache (explaining "then it does
   show results!"), and whether an interrupt can poison the cache. Does the
   cache survive `saveRDS`/`readRDS`?
4. Regression tests matching whichever route ships: correctness of the prior
   PIPs (against a long reference chain), caching, interrupt behavior.

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

### F-038 + F-039 (minor, DECIDED by MM) — display Bayes-factor evidence in NATURAL log, capped

MM's decision (2026-08-01): user-facing displays print the **natural** log
Bayes factor, labeled "log BF" — not log10 ("nobody does that"), not raw BF —
and large magnitudes are CAPPED to avoid fake accuracy, e.g. `log BF > 10,000`.

Background: classification is CORRECT today (`build_verdicts`:
`lthr = log10(evidence_threshold)`, boundaries ±1 for threshold 10) but the
`log10_bf` column invited a base misread — MM read −1.474 as natural log and
filed it as a classification bug. The unit itself is now decided: natural log.

Implement consistently:
1. `verdicts()`: the displayed evidence column becomes natural-log BF
   (`log_bf`), classification boundaries at `±log(evidence_threshold)`
   (≈ ±2.303 for 10); the printed header states them in the displayed unit
   ("presence: log BF > 2.30; absence: log BF < −2.30 (threshold 10)").
   The verdicts object is new in 0.2.0.0 — rename the stored `log10_bf`
   field to `log_bf` (natural) rather than carrying both; update Rd and any
   internal consumers (`R/plot_bgms.R` verdict tables, sensitivity internals
   that read verdict frames). Grep for `log10_bf` package-wide.
2. `plot_edge_posterior()` title: print `log BF = <x.x>` in natural log, and
   past a magnitude cap print an inequality (`log BF > 10,000`). No raw BFs,
   no 131-digit strings. Snapshot-test a saturated-edge title.
3. `extract_inclusion_bf(log = TRUE)` already returns natural log — say so
   in `?verdicts` so the printed unit and the extractor agree.
4. CONSISTENCY FLAG, do not change: `prior_sensitivity_check()` output and
   plot currently speak log10 ("0.41 log10 BF" noise line, curve scale). Its
   figures ship in the tutorial manuscript, so converting is NOT yours to
   decide — record in your report where its log10 usages live (file:line)
   and leave them; MM settles the package-wide unit with the guidelines
   terminology sync.

### F-018 (BLOCKER, fix authorized) — sampler path must honour CRAN's 2-core limit

Confirmed by MM ("CRAN rejects if you use more than 2 cores"). Report 01 f4
measured: 76 uncapped example fit calls sustain ~3.9 cores; the env var
`_R_CHECK_LIMIT_CORES_` is ignored on the sampler path
(`R/run_sampler.R:365` passes `s$cores` straight through), while
`R/correction_tables.R:151-160` (`normalize_builder_cores`) already implements
the correct guard.

Fix: apply the same guard on the sampler path — resolve the effective core
count where `run_sampler()` reads `s$cores` (cap to 2 when
`_R_CHECK_LIMIT_CORES_` is set and truthy, per `normalize_builder_cores`'s
exact semantics, Windows nuance included; keep `detectCores()` as the
interactive default). Also audit `simulate_predict.R:127` and any other
`cores =` consumer reaches the same guard. Verify empirically with brief 01's
probe (`~/bgms-review/probe-cores.R`): user/elapsed ratio for a default
`bgm()` call must drop to ≤ 2 under `_R_CHECK_LIMIT_CORES_=TRUE`, unchanged
without it. Add a regression test (env-var set → resolved cores ≤ 2).

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

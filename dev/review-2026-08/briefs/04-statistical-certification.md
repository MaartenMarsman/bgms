# Brief 04 — Statistical certification of rc1: slow-tier suites + zratio route certificates (Opus agent)

You are the execution agent in the pre-CRAN review of the R package **bgms**.
You share no context with the review lead. REPORT-ONLY brief: run, measure,
compare — change nothing.

## Context

Report 01 established that `R CMD check` and the standard test tiers are clean —
but ALL SBC certification, parameter recovery, and NUTS-vs-MH cross-validation
(8 files, 96 tests) sit behind the env var `BGMS_RUN_SLOW_TESTS` and only run in
a twice-weekly scheduled CI job (finding F-044). The hierarchical/zratio path
additionally has a bank of gold-reference validation scripts in
`dev/validation/` that were last run during PRs #193/#194. Neither has been
executed against the frozen release candidate as a unit. That is this brief.

- Frozen target: tag `v0.2.0.0-rc1`.
- Reuse brief 01's artifacts in `~/bgms-review/`: the clean source export
  `bgms-rc1/` and the installed library `lib/` (verify
  `packageVersion("bgms")` is 0.2.0.0 from `~/bgms-review/lib`; if the lib is
  missing, rebuild per brief 01 §Tasks 1+4).
- The Dropbox repo checkout (`.../Projecten/R/bgms`, currently on `main`) may
  be read for scripts but NEVER built in or written to. Copy anything you need
  to run into `~/bgms-review/certification/`.

## Task A — the dormant statistical suites, against rc1

From `~/bgms-review/bgms-rc1/`, with the rc1 library first on `.libPaths()`:

```sh
NOT_CRAN=true BGMS_RUN_SLOW_TESTS=true Rscript run-slow-suite.R
```

(where the runner mirrors brief 01's `run-suite.R`). This unlocks, at minimum:
`test-sbc-ggm.R`, `test-sbc-correction.R`, `test-mixed-nuts.R`,
`test-parameter-recovery-ggm.R`, `test-scaling-diagnostics.R`,
`test-validation-slow.R`, `test-zratio-cauchy.R`, `test-zratio-law.R`.

Report: per-file pass/fail/skip + runtime; EVERY failure verbatim; total wall
time. If a failure looks stochastic, rerun that file once with the same seed
settings and say whether it reproduces. These are the package's core
statistical-correctness gates — treat any failure as at least `major`.

**Scope caveat from MM (the maintainer): SBC is only a valid certification for
the GGM path.** The ordinal (omrf), mixed, and bgmCompare paths target
pseudolikelihood approximations, so simulation-based calibration against
simulated data does not certify them the same way — their gates are the
recovery/cross-validation suites instead. Interpret accordingly: do not read
any pseudolikelihood-path behavior as an "SBC failure", and if any test appears
to run SBC-style checks against a non-GGM path, FLAG it as a finding (the test
may be miscalibrated by construction) rather than interpreting its outcome.

## Task B — zratio route certificates and gold-bank spot checks

`dev/validation/` (present in the rc1 export) contains runnable certification
scripts with banked `.rds` references. Run at least:

- `zratio_deployed_route_cert.R` (deployed-route identity: `log_zratio == logR`)
- `zratio_isolated_route_cert.R`
- `zratio_anchor_gate.R`
- `zratio_shape_verdict.R`
- `zratio_collapse_reachability.R`

Method: copy each script plus whatever it reads into
`~/bgms-review/certification/`, read its header to get the intended invocation,
run against the rc1 install, and compare outputs to the banked `.rds` in
`dev/validation/` (read the banked files; NEVER overwrite them). Some scripts
were written PR-side and may need a path variable adjusted — record every
adaptation you make. If a script is not runnable as banked, say so rather than
improvising a replacement.

Report per script: what it certifies (one line), pass/deviation, and for
deviations the numbers side by side. The banked gold bank
(`dev/validation/zratio_gold_bank.md`) documents tolerances — judge against
those, not zero.

## Task C — MM's sensitivity drifter, reproduced properly

MM's quick pass (report 02) saw one large drifter (edge `intrusion-anger`) in
`plot(prior_sensitivity_check(fit))` on a Blume-Capel fit. Reproduce with a
publication-grade fit:

```r
fit = bgm(Wenchuan, variable_type = "blume-capel", baseline_category = 1,
          seed = 1)                       # defaults: 4 chains, iter 2000
sens = prior_sensitivity_check(fit, seed = 1)
```

Report: total runtime; the printed report verbatim; for `intrusion-anger`
specifically — its verdict at each anchor, whether the check classifies it
robust / within-noise / beyond-noise / not-certifiable, and whether the curve
shape matches MM's saved plot (`dev/review-2026-08/reports/assets/`, readable
on the Dropbox checkout). The question to answer: is that drift a stable
property of the data-prior combination (fine — the check is doing its job) or
an artifact of the quick fit (also fine) — or something unstable that
reproduces differently run to run (a finding).

## Task D — quantify the hierarchical prior-chain cost (feeds F-036)

On a hierarchical fit (`bgm(Wenchuan, variable_type = "continuous",
precision_graph_prior = "hierarchical", seed = 1)`):

- Time the FIRST `verdicts(fit)` call (wall clock) — this triggers the silent
  5000-iteration prior-only chain (F-036).
- Time the SECOND call — confirm whether caching makes it instant, and where
  the cache lives (does it survive `saveRDS`/`readRDS` of the fit?).

Report both timings and the cache observation. Do not fix anything — brief 03
owns the fix; your numbers calibrate its messaging.

## Deliverable

Write `/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms/dev/review-2026-08/reports/04-statistical-certification.md`
(write it even though the path is untracked on the `main` checkout — the lead
collects it):

1. **What was done** — environment, library provenance, exact commands.
2. **Findings** — numbered, severity-tagged (blocker / major / minor / note);
   a suite failure or certificate deviation is at least major.
3. **Evidence** — verbatim failures/deviations, per-file and per-script tables,
   timings; artifact paths under `~/bgms-review/certification/`.
4. **Open questions** — anything unrunnable, ambiguous tolerances, anything
   needing MM's statistical judgment.

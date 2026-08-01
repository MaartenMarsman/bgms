# Brief 17 — docs + small-defect fix batch (Opus agent)

You are the execution agent in the pre-CRAN review of the R package **bgms**.
You share no context with the review lead; everything you need is here. CODE
CHANGES authorized on a fix branch. This batch lands every fix from the
report-15 docs sweep plus four small code defects found alongside them.
Report 15 (`dev/review-2026-08/reports/15-docs-accuracy-sweep.md`) is
REQUIRED READING — it carries the reproduction for every item below.

## Machine budget (standing rule)

This machine has 15 cores and is SHARED (another agent may be running a long
fit batch). Cap your total footprint at ~6 hardware threads; run the test
suite and vignette builds SEQUENTIALLY, one process at a time. Runtime is
not a grading criterion.

## Setup

- Repo (Dropbox; do NOT switch its checked-out branch or build in it):
  `/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms`
- **Base on `develop` AT OR AFTER the brief-12 merge** (the maintainer
  launches you only once it has landed; verify with
  `git log --oneline -5 develop` that a nightly-respec merge is present —
  if it is not, STOP and say so in your report rather than proceeding):
  ```sh
  git -C /Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms \
      worktree add ~/bgms-review/wt-fix7 -b fix/docs-batch develop
  cd ~/bgms-review/wt-fix7
  ```
- One commit per item or small coherent group; `docs:`/`fix:` prefixes with
  F-numbers. No pushing; no attribution trailers. Build/test only in
  exports or the worktree, never the Dropbox tree.
- `R CMD build` needs pandoc on PATH (F-045): use RStudio's bundled pandoc
  as report 15 did.

## Tasks

### Vignette corrections

1. **F-087 (major) — the inverted Bayes-factor reading**,
   `vignettes/diagnostics.Rmd:146,152-153`. The chunk prints
   `BF_10 = 24.59`; the prose reads it as "little evidence for inclusion"
   and 1/BF = 0.041 as "strong evidence for the absence". Rewrite the prose
   to read the actual numbers correctly (this edge is `presence`, log BF
   3.20 — agree with `verdicts()`), and keep the 1/BF sentence only as the
   correct illustration of reading evidence for absence (i.e., describe
   what a small BF_10 *would* mean, or pick an edge that actually shows
   absence — your choice, but the worked numbers and the words must agree).
2. **F-086 (major) — the mock sensitivity report**,
   `vignettes/prior-sensitivity.Rmd:62-103`. Lead decision: make the chunk
   LIVE (report 15 measured ~33 s for `bgm(Wenchuan[, 1:6], chains = 2)` +
   `prior_sensitivity_check(fit)`) so the shown output IS the shipped print
   method's output and can never drift again. Seed it. This also resolves
   **F-090** (the prose naming a nonexistent "Details line" and
   "not certifiable" wording — rewrite those sentences against the live
   output: the footer is `Method:`/`Refits:`/`Noise:`, the wording is
   "too noisy to assess").
3. **F-088** — `prior-sensitivity.Rmd:143`: the lingering
   `\log_{10}` becomes natural log, matching its neighbours.
4. **F-089** — `diagnostics.Rmd:58`: the summary carries `n0->1`/`n1->0`
   only; stop naming `n0->0`/`n1->1` as reported columns.
5. **F-092** — `diagnostics.Rmd:84-86`: "NA only for … constant to double
   precision" is too narrow; add the second route (variance below the
   autocovariance kernel's numerical floor — `R/mcmc_summary.R:263-268`),
   ideally using report 15's observed example (sd 2.2e-11, NA mcse/n_eff,
   finite Rhat).
6. **F-094** — `diagnostics.Rmd:64-66`: reword the full-chain-ESS sentence;
   `mean`/`sd`/`Rhat` use the full chain but the reported `n_eff` is the
   composite RB estimator, and `R/mcmc_summary.R:334-336` deliberately
   rejects the raw-chain route. The sentence must not justify something the
   package chooses not to do.
7. **F-002 (one sentence)** — `vignettes/comparison.Rmd`: state the
   asymmetry: `bgm()` defaults to `interaction_prior = normal_prior(1)`
   while `bgmCompare()` defaults to `cauchy_prior(1)`, so comparing a
   compare fit against separate `bgm()` fits at stated defaults compares
   different models; matching priors explicitly is the remedy.
8. **F-026 (minimum only)** — `vignettes/intro.Rmd`: (a) the variable-type
   sentence must include the mixed discrete+continuous case; (b) make the
   GGM demo chunk runnable as written (define its data or use a shipped
   dataset; it may stay `eval = FALSE`, but a reader pasting it must not
   hit an undefined object); (c) add a short "what's new in 0.2.0" pointer
   paragraph naming verdicts(), the prior constructors, GGM/mixed models,
   and the checking vignettes. Do NOT rewrite the vignette beyond this —
   the full rewrite is the maintainer's own pass.

### Code fixes

9. **F-091 — `print.bgms_verdicts` subset error**, `R/verdicts.R:451-470`.
   `subset(v, fragile)` errors because `[.data.frame` drops the
   `evidence_threshold` attribute while all five display columns survive;
   the guard tests columns only. Lead decision: extend the guard to ALSO
   test the attribute (`is.null(attr(x, "evidence_threshold"))` → plain
   data-frame print fallback). Do NOT add a `[.bgms_verdicts` method (that
   is a recorded backlog item). Regression test: `print(subset(v, fragile))`
   succeeds. Ride-alongs in the same print: "(1 indicators)" pluralization
   and the `%g` Bayes-factor formatting ("0.0333333") — make both read
   cleanly.
10. **F-093 — wrong order-of-magnitude in runtime messages**,
    `R/zratio_surfaces.R:545-546` and `R/zratio_gauge.R:392-395`: 0.003 /
    0.00028 = 10.7 — say "roughly a tenth of" (matching the vignette), not
    "two orders below".
11. **F-095 — chain-boundary transition fabrication**,
    `src/mcmc_diagnostics.cpp:352-371`: the pooled scan carries `prev`
    across chain boundaries, counting `nchains − 1` spurious transitions
    per parameter. Reset per chain (`int start = 1;` with `prev` taken from
    each chain's own first draw). Add a unit test pinning: a 2-chain
    indicator that is all-zeros in chain 1 and all-ones in chain 2 reports
    `n0->1 = 0`. Re-record any snapshots whose transition counts shift, and
    LIST every re-recorded snapshot in your report with its before/after.
12. **F-031 (two examples)** — add `\examples{}` to `extract_inclusion_bf`
    and `extract_prior_inclusion_probabilities`, following the sibling
    extractors' pattern exactly (same fixture scale, same
    `\donttest` discipline — examples run on CRAN, mind F-018: no
    uncapped cores).

### Harmonization check (from the report-14 integration note)

13. **F-049 Rd sentence** — brief 12 landed a sentence in
    `?prior_sensitivity_check` Details about extrapolation error and anchor
    exactness. Verify it states exactness the way the code does: exactness
    lives on the anchor fits' own RB statistics (per-anchor verdict
    columns, chosen-scale quantities), NOT on the curve row at the anchor's
    scale — the pooled curve is pooled everywhere
    (`R/anchor_curve.R:207-211`). If the sentence says "the anchors are
    exact" without that distinction, sharpen it to match the NEWS wording
    brief 14 merged. One sentence; cite both files in your report.

### Explicitly OUT of scope

- The "NUTS issues: Warmup may be incomplete" notices in three flagship
  examples: deliberate honesty on short demo fits — leave them.
- The docs-site / tutorial copies of the F-087 paragraph: out-of-repo,
  tracked in the maintainer's plans folder.
- The full intro.Rmd rewrite (F-026), the wider 25-file examples sweep
  (F-031), a `[.bgms_verdicts` method (F-091 backlog).

## Verification gate

1. Full LOCAL suite (default tier, no slow env vars): 0 failures,
   0 warnings; plus the slow tier IF `src/` changed (it does — F-095):
   run with `BGMS_RUN_SLOW_TESTS=true` too and report both.
2. All five vignettes render clean (0 errors / 0 warnings / 0 messages),
   sequentially; report the new prior-sensitivity render time (the live
   chunk must keep the total vignette budget sane — state it).
3. `R CMD check --as-cran` on a `git archive` tarball: the 2 baseline
   NOTEs only.
4. The F-087 chunk's printed numbers and its prose agree — quote both in
   the report. The F-086 chunk's output is live — confirm no pasted block
   remains.

## Deliverable

Write
`/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms/dev/review-2026-08/reports/17-docs-fix-batch.md`:
What was done (per item, with commits) / Findings (anything new the fixes
surfaced) / Evidence (gate outputs, before/after for F-087 and F-086,
snapshot re-record list) / Open questions.

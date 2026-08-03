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
- **Base on `develop` AT OR AFTER the brief-12 merge `f24ad3c8`** (verify
  with `git log --oneline -5 develop` that it is present — if it is not,
  STOP and say so in your report rather than proceeding):
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

13. **F-049 Rd sentence — verify only.** The lead checked the merged
    sentence at 12's integration: it correctly places exactness on the
    anchor fits' own RB statistics with the between-anchor curve
    importance-reweighted (`R/prior_sensitivity.R:71-78` /
    `man/prior_sensitivity_check.Rd:144-151`), matching
    `R/anchor_curve.R:207-211`. Your job: cross-check that wording against
    the NEWS sentence brief 14 merged and edit ONLY if the two disagree on
    where exactness lives. One line in your report either way.

### NEWS amendments from the maintainer's ratification read (F-096)

14. Lead-verified facts at `cran-0.1.6.3`, which you build on rather than
    re-derive: (i) BOTH released samplers use the old scale with NO factor 2 —
    compare: `src/bgmCompare/bgmCompare_sampler.cpp:118-139`
    (`category * rest_score` off plain `observations · effects`); bgm:
    `src/bgm/bgm_logp_and_grad.cpp:178`; simulator: `src/mrf_simulation.cpp:98`
    — i.e. the released pair was MUTUALLY CONSISTENT, and the
    "compare = 2× bgm" mismatch existed only mid-development; (ii) the tag's
    NAMESPACE exports `predict`/`simulate` methods for both classes, plus
    `simulate_mrf` and `mrfSampler`. Amendments, all in NEWS.md:
    a. REWRITE the bgmCompare association-scale entry for the 0.1.6.3
       reader: both `bgm()` and `bgmCompare()` moved to the association
       scale together in 0.2.0.0; compare pairwise values are ~half their
       0.1.6.3 values (mirroring the bgm() entry). KEEP refit-not-rescale,
       the prior-tightening consequence, and the scale-contingency caveat.
       DROP the omega-vs-2omega internal-mismatch narrative.
    b. DROP the "predict()/simulate() on a bgmCompare fit were wrong"
       entry outright (dev-only defect; at 0.1.6.3 sampler and predictor
       agreed). The Wenchuan numbers go with it.
    c. The `simulate_mrf()` bullet STAYS as merged (it WAS exported at
       0.1.6.3 — verified).
    d. The hamiltonian-mc removal bullet gains one clause cross-referencing
       the GGM `update_method = "gibbs"` addition.
    e. The `standardize` bullet gains the maintainer's rationale: a
       forthcoming g-prior formulation of the interaction prior addresses
       standardization in the model itself; the per-pair max-score
       adjustment is retired in favour of that direction.
    f. The compare RB-NA bullet gains the maintainer-requested
       clarification: an always-included edge indicator in `bgm()` is still
       PROPOSED every iteration, so RB evaluates its conditional odds and
       reports a finite Bayes factor; an unselected main-difference
       indicator is never proposed at all, so no RB quantity exists — NA is
       the honest value (`R/extractor_functions.R:342`: NA ⇔ zero visits).
       Mirror the clause in the Rd where the NA is documented if a natural
       slot exists; say so in the report if not.
    g. Apply the never-shipped rule UNIFORMLY to the Bug-fixes section: for
       EACH entry, establish at `cran-0.1.6.3` whether the defect could
       reach a user there. Mixed-model and GGM-only fixes (PIP block order,
       cross-indicator asymmetry, `delta = NULL` mixed default, Cholesky
       downdate) — confirm their features are absent at the tag and DROP.
       Verify the ambiguous ones AT THE TAG before deciding: the
       category-scale recode fix, the Alpine/musl include, the two
       imputation-cache fixes, and the two NUTS fixes (`target_accept`
       pass-through; acceptance accumulation — NUTS EXISTED at 0.1.6.3, so
       these may be genuine user-facing fixes: keep with 0.1.6.3 framing if
       the tag shows the defect). Your report carries a per-entry verdict
       table with tag citations. Re-run the markdown-NEWS parse gate from
       report 14 after all edits.

### Archive-and-remove the dormant analytic law (F-099) + self-contained comments (F-097)

15. The maintainer has decided (2026-08-02): the DORMANT analytic law
    leaves the package and lives on its own documented branch. Verified
    seams: the only include in the tree is `src/zratio_test_interface.cpp:14`;
    nothing deployed touches it. Execute in this order:
    a. **Archive first.** From your fix branch BEFORE any removal, create
       `archive/zratio-analytic-law`. On it, one commit: rewrite the three
       companion-referencing comments in `src/models/ggm/zratio_law.h`
       self-contained (`:24` — describe what the engine IS: a
       deterministic, tableless CN anchor engine, one self-consistent
       spectral solve per (size, density) cell yielding S1, S2 — without
       the port framing; `:47` — "composite trapezoid weights on a
       non-uniform grid"; `:61` — state the notation map as a fact: "eta
       here is the tilt rate often written beta; t2 = 2*eta*sigma^2"), and
       add a short archival note at the top of the header: why archived
       (maintainer decision, 2026-08-02, review F-099), the re-wiring
       condition (the large-q crossover the banner already describes), and
       that regenerating `fixtures/zratio_law_reference.rds` requires the
       companion R implementation (private). Keep the DORMANT banner. Do
       NOT push the branch — the lead pushes it at integration.
    b. **Remove from the package** (on the fix branch): delete
       `src/models/ggm/zratio_law.h`, the `zratio_law_moments` block in
       `src/zratio_test_interface.cpp` (include line + comment + function),
       `tests/testthat/test-zratio-law.R`,
       `tests/testthat/fixtures/zratio_law_reference.rds`, and
       `tests/testthat/fixtures/make_zratio_law_reference.R`. Regenerate
       RcppExports (Rcpp::compileAttributes). Verify with a tree-wide grep
       that no reference to `zratio_law` survives in `src/`, `R/`, or
       `tests/` (dev/ record files are fine).
    c. **Fix the fourth F-097 site in-package**: `src/models/ggm/
       zratio_engine.h:229` — state what the reference IS (the
       per-component Monte-Carlo oracle evaluation of the same
       decomposition) with no companion referent.
    d. **One MAINTAINERS line** (architecture section): the analytic law
       is archived on `archive/zratio-analytic-law`, with the re-wiring
       condition.
    e. **No NEWS entry** — never deployed, never user-visible (the
       never-shipped rule applies to code too).
    f. In your report: the tarball-size and test-time deltas, and note
       that `test-zratio-law.R` vanishing from the tier classification is
       AUTHORIZED by F-099 (not a silent drop).

### CI hygiene (F-102)

16. Add `paths-ignore: ['dev/**']` to the `on: push` triggers of FOUR
    workflows: the three r-lib ones (`.github/workflows/` lint /
    R-CMD-check / test-coverage — the exact filenames as found) AND the
    new `fast-checks.yaml` (its `push:` trigger only — leave its
    `pull_request:` trigger exactly as is). Review-record commits must
    stop launching package CI: dev/ files cannot change any test outcome.
    Do not touch `nightly-validation.yaml` or `weekly-certification.yaml`
    (schedule-triggered; no push trigger to filter).

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

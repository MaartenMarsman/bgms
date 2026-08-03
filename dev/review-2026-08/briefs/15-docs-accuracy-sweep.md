# Brief 15 — vignette + man-page accuracy sweep (Opus agent)

You are the execution agent in the pre-CRAN review of the R package **bgms**.
You share no context with the review lead; everything you need is here.
**ANALYSIS ONLY — change nothing.** Discrepancies you find become findings in
your report; the fixes land in a later authorized batch (another agent is
editing Rd/roxygen in parallel — write access here would collide).

## Machine budget (standing rule)

This machine has 15 cores and is SHARED. Cap your total footprint at ~6
hardware threads; SEQUENCE the vignette renders (one at a time); verification
snippets stay tiny (p ≤ 6, seconds — the vignettes' own code is the upper
bound of what you run). Runtime is not a grading criterion.

## Setup

- Repo (Dropbox; do NOT switch its checked-out branch or build in it):
  `/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms`
- Work from a clean export of `develop` AT OR AFTER `5b850410`; build and
  install ONCE into a private library (this is your oracle):
  ```sh
  mkdir -p ~/bgms-review/val15 && cd ~/bgms-review/val15
  git -C /Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms \
      archive develop --prefix=bgms-val/ | tar -x
  R CMD build bgms-val && R CMD INSTALL --library=~/bgms-review/lib-val15 bgms_0.2.0.0.tar.gz
  ```
- REQUIRED READING first: `dev/review-2026-08/FINDINGS.md` — the docs rows
  **F-005, F-006, F-026, F-031** are KNOWN; do not re-discover them (but DO
  verify F-006's substance yourself, below). Rows F-049/F-057 have Rd lines
  landing in a parallel brief; skip those two surfaces entirely.

## Tasks

1. **Vignette render.** All vignettes in `vignettes/` (`intro`, `comparison`,
   `diagnostics`, `checking-your-model`, `prior-sensitivity`), rendered
   sequentially against the installed build. Any error/warning during render
   is a finding.
2. **Claim-by-claim accuracy.** For each vignette, extract every checkable
   factual claim — default values, argument names, output semantics, printed
   conventions, estimator descriptions — and verify it against the installed
   build (run it) or the source (cite file:line). Produce one table per
   vignette: claim / how checked / verdict. Specifics that MUST be in the
   sweep:
   - the natural-log convention: every Bayes-factor mention should read as
     ln-based ("log BF"); any lingering log10 phrasing or stale threshold
     numbers (e.g. old `1.15/6.91`-style caps mismatching printed output) is
     a finding;
   - the F-006 surface: the diagnostics vignette's mixture-ESS description,
     verified against `src/mcmc_diagnostics.cpp:375-387` — report what the
     code computes vs what the vignette says (the lead holds a corrected
     reading; your independent read cross-checks it);
   - stated defaults vs actual formals (`iter`, `warmup`, `target_accept`,
     `interaction_prior` — these changed in 0.2.0; also
     `bgm()` = `normal_prior(1)` vs `bgmCompare()` = `cauchy_prior(1)`, an
     asymmetry users must not learn the hard way);
   - every code chunk's shown output vs what the installed build actually
     prints today.
3. **Man-page spot-checks.** The 23 new exports (list in
   `dev/review-2026-08/reports/00c-r-api-diff-map.md`): run each export's
   examples against the installed build; check the Description/Details
   claims match observed behavior; check documented defaults equal actual
   formals. Table: export / example runs? / claims verdict.
4. **Cross-consistency.** Where a vignette, a man page, and printed runtime
   output describe the same thing (verdict categories, evidence thresholds,
   fragility caveats, group numbering), flag any disagreement among the
   three — cite all sides.

## Verification gate (for the report)

Every claim table states how the check was run (snippet or file:line). Every
finding carries a minimal reproduction. Coverage statement at the end: which
man pages were NOT checked (anything outside the 23 new exports you didn't
reach), listed explicitly — silent partial coverage is a FAIL.

## Deliverable

Write
`/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms/dev/review-2026-08/reports/15-docs-accuracy-sweep.md`:
What was done / Findings (severity-tagged; new discrepancies only, known
rows cross-referenced) / Evidence (the tables) / Open questions.

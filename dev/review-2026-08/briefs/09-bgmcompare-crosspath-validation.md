# Brief 09 — bgmCompare cross-path consistency validation (Opus agent)

You are the execution agent in the pre-CRAN review of the R package **bgms**.
You share no context with the review lead; everything you need is here.
**ANALYSIS ONLY — no package code changes.** Defects you find become findings
in your report; fixes come later. This brief is the numbers-level validation
of the bgmCompare path — the review's strategy hunts cross-path convention
divergence (the 0.2.0.0 breaking fix exists because bgmCompare once silently
used `omega·x` where every other path used `2·omega·x`), and bgmCompare is
the one sampler with no C++ test interface, so end-to-end statistical
validation is its primary defense.

## Setup

- Repo (Dropbox; do NOT switch its checked-out branch or build in it):
  `/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms`
- Work from a clean export of `develop` AT OR AFTER merge `4969f843` (all
  compare fixes, the calibration method, and the nats conversion are in —
  validating older code wastes the run):
  ```sh
  mkdir -p ~/bgms-review/val09 && cd ~/bgms-review/val09
  git -C /Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms \
      archive develop --prefix=bgms-val/ | tar -x
  R CMD build bgms-val && R CMD INSTALL --library=~/bgms-review/lib-val09 bgms_0.2.0.0.tar.gz
  ```
- Context on the worktrees/exports: reports 03–06 under `dev/review-2026-08/
  reports/` (readable via `git show develop:<path>`). Report 06 fixed the
  compare-path defects this brief validates on top of.
- Seeds fixed and reported for every run. Runtimes noted. Where a check is
  Monte-Carlo, derive its tolerance from measured run-to-run noise (fit the
  same cell twice with different seeds) and state the derivation — no bare
  magic thresholds.

## Items

### 1. Split-halves identity (no true differences ⇒ no difference evidence)

Take `Wenchuan` complete cases, randomly split rows into two "groups"
(seeded), fit `bgmCompare(x, group_indicator = split, seed = ...)` at
defaults with enough iterations for stable RB estimates. There are no true
group differences, so:

- Difference verdicts: no `presence` verdicts beyond a false-positive budget
  you derive from the prior (state it); the bulk should be absence/undecided.
- Report the full difference-PIP distribution (max, count > 0.5) and the
  difference posterior means (should straddle 0).
- Repeat over ≥ 3 split seeds to show it is not one lucky split.

### 2. Planted-difference recovery

Simulate two groups from known MRFs that differ in a controlled way: take a
fitted Wenchuan-like parameter set, plant differences on ~4 pairwise
associations (magnitudes spanning small→large ON THE ASSOCIATION SCALE — mind
the ×2 convention; a planted difference of δ in the stored parameter means
2·δ·x in the linear predictor), simulate n per group at two sizes (e.g. 400
and 2000) with `simulate_mrf()`, fit `bgmCompare`, and report:

- Difference-PIP separation: planted vs unplanted edges (an ROC-style summary
  or the two PIP distributions).
- Posterior difference means vs planted truth (bias, and whether truth sits
  inside 95% intervals ~95% of the time across edges).
- How recovery scales with n and with magnitude. This is the operating-
  characteristics picture the review cites at sign-off.

### 3. Per-group estimates vs separate single-group fits (THE convention check)

On the split-halves data (item 1's fixed split), also fit each half separately
with `bgm()` at matched priors. From the compare fit, reconstruct each
group's pairwise associations (baseline + group offset, per the fit's own
extractors). Compare, per group:

- Correlation and max |Δ| between compare-reconstructed and separate-fit
  posterior means of pairwise associations, against a tolerance derived from
  seed-to-seed refit noise of `bgm()` itself.
- CRITICALLY: test for a SYSTEMATIC ×2 (or ×½) factor — regress compare
  estimates on separate estimates and report the slope with its interval.
  A slope near 2 or 0.5 is the resident defect class; near 1 is health.
- Note where disagreement is EXPECTED (difference-selection shrinkage when
  `difference_selection = TRUE` pools groups): run this item with
  `difference_selection = FALSE` for the clean comparison, and once with
  defaults to show/quantify the pooling effect.

### 4. Convention-guard audit (read, don't plant)

Read `tests/testthat/test-bgmCompare.R` around line 337 (the ~7 s every-run
cross-implementation guard, review law #1). Report: what exactly it pins
(which quantities, which paths, what tolerance), whether it would catch a
reintroduced ×2 divergence on (a) the pairwise likelihood, (b) simulate/
predict, (c) the mixed cross terms, and any coverage hole (e.g. a path the
guard never exercises). Reading and reasoning only.

### 5. Blume-Capel end-to-end on the compare path

One compare fit containing a Blume-Capel variable (Wenchuan column with a
sensible baseline, or simulated): fit, `verdicts()`, `predict()` (this is the
F-068 regression surface — report 06 fixed baseline mis-centering up to 0.23
probability), `calibration_check()` per group, and `simulate()` round-trip
sanity (simulated margins vs model expectations). Report anything that looks
off; the F-068 manual-reference test pins the convention, so disagreement
here would be new.

## Verification gate (for the report, since no code changes)

Every item reports: seeds, runtimes, derived tolerances with their
derivations, and a pass/fail verdict per check. Any FAIL gets a minimal
reproducible snippet and a severity-tagged finding.

## Deliverable

Write
`/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms/dev/review-2026-08/reports/09-crosspath-validation.md`:
What was done / Findings (severity-tagged; item-3 slope FIRST) / Evidence
(tables + any figures to `reports/assets/`) / Open questions.

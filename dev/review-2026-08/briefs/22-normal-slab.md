# Brief 22 — compare slab default Cauchy → Normal, the ridge picture, and the zero-support test evaluation (Opus agent)

You are the execution agent in the pre-CRAN review of the R package **bgms**.
You share no context with the review lead; everything you need is here. CODE
CHANGES authorized on a fix branch. Maintainer decision F-119: bgmCompare's
slab default moves from Cauchy to Normal, mirroring `bgm()`. Two
maintainer-requested evaluation artifacts ride along.

## Machine budget (standing rule — TWO other agents may be running)

~4 hardware threads; ONE model fit at a time; heavy work last (task 5 is the
heavy tail). Runtime is not a grading criterion.

## Setup

- Repo (Dropbox; do NOT switch its branch or build in it):
  `/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms`
- Base on **`origin/develop`** (NOT the local ref) AT OR AFTER `c2fb48e1`:
  ```sh
  git -C /Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms \
      worktree add ~/bgms-review/wt-fix14 -b fix/normal-slab origin/develop
  cd ~/bgms-review/wt-fix14
  ```
- One commit per task; `fix(compare):`/`test:`/`docs:` prefixes with F-119.
  No pushing; no attribution trailers.
- **HARD off-limits (an unmerged batch owns them):** `R/methods_bgmcompare.R`,
  `tests/testthat/test-summary-marking.R`, `man/summary.bgmCompare.Rd`,
  `R/validate_data.R`. Also as always: `NEWS.md` (verbatim proposal in the
  report; the lead lands it), plot files, `dev/review-2026-08/` except your
  report + assets. `~/bgms-review/val16/` is READ-ONLY input.

## Tasks

### 1. The switch (F-119)

`bgmCompare()` currently defaults `interaction_prior = cauchy_prior(scale = 1)`
(`R/bgmCompare.R:242`); `bgm()` defaults `normal_prior(scale = 1)`. Before
changing anything, MAP which parameter families draw from this slab —
baseline pairwise, pairwise differences, main-effect THRESHOLD differences —
with file:line citations; if main differences are governed by a different
prior object, the maintainer's intent covers the DIFFERENCE priors, so find
and flip that too, and state exactly what governs what. Then: default becomes
the Normal mirroring `bgm()`; `cauchy_prior()` stays fully available; roxygen
and Rd sweep for every place the default is named.

### 2. Tests under the new default

Full suite, both tiers. Every expectation that legitimately shifts is
re-tuned WITH its derivation stated (the compare convention guard's planted-δ
tolerances were derived, not guessed — re-derive, don't inflate). List every
touched expectation in the report. No snapshot re-record without a stated
reason.

### 3. Operating-characteristics transfer (3 seeds)

Reuse `~/bgms-review/val16/out/` truths and the saved `bgm()` baselines
(READ-ONLY; do not refit bgm). Refit compare on seeds 01–03 under the new
default, one at a time; compare group slopes, noise sd, δ = 0.40 recovery,
detection, false presences against report 19's fixed-build table. Escalate to
the full ten seeds ONLY if any number moves beyond the seed spread report 16
established; otherwise three suffice and say so.

### 4. The ridge visualization (maintainer request — report asset)

The identification geometry of a structural-zero cell, drawn EXACTLY, not
schematically. Construct the adverse fit (four-level ordinal, supports
{0,1,2} vs {1,2,3}, seeded, seconds-scale). For the zero-support cell, grid-
evaluate the TRUE log-posterior over (overall threshold, group difference)
via the `bgmCompare_test_logp_and_gradient()` hook, all other parameters held
at their posterior means from a reference fit. Three same-axes contour
panels:
a. the likelihood surface alone (compute the prior term analytically in R
   and subtract it from the hook's log-posterior) — the one-sided ridge;
b. the posterior under `cauchy_prior(scale = 1)`;
c. the posterior under `normal_prior(scale = 1)` — the regularization made
   visible.
PNG(s) + the generating script under
`dev/review-2026-08/reports/assets/` (`f119-ridge-*`). Paper-grade axes and
labels; the maintainer will judge the figure and may reuse it for teaching.

### 5. Difference-test evaluation on zero-support cells (maintainer request)

How do the main-difference inclusion verdicts behave when a cell is empty?
Three sub-studies, sequential, thin grids — report every dropped cell:
a. TRUE shift: the non-observing group's size n₁ ∈ {100, 400, 1600} × the
   category's rate in the observing group {common ~0.25, rare ~0.05}; one
   fit per cell; report the main-difference inclusion BF's direction and
   magnitude per cell, under the NEW default.
b. SAMPLING zero (the over-call risk): both groups share one distribution
   with a rare category (~0.02), n per group ∈ {50, 100}, ~20 seeded
   replicates each; report how often difference selection calls the empty
   cell's difference SUPPORTED at BF 10 — that rate is the finding.
c. Slab sensitivity: Cauchy vs Normal vs scale on the (a) fits — use
   `prior_sensitivity_check()` reweighting IF it covers main-difference
   verdicts (verify and say); refit only where it does not.

### 6. NEWS (report-only)

0.1.6.3 shipped bgmCompare WITH Cauchy priors, so this default change is
visible to the tag reader — it earns a real entry, maintainer-directed. Find
how `bgm()`'s own Cauchy→Normal default switch is worded in NEWS.md and
propose (VERBATIM, in the report) a matching one-line entry for bgmCompare,
including that results under defaults change and `cauchy_prior()` remains
available. Do not edit NEWS.md.

## Verification gate

1. Both suite tiers: 0 failures / 0 warnings; every re-tuned expectation
   listed with derivation.
2. `R CMD check --as-cran` on a `git archive` tarball: 2 baseline NOTEs.
3. `devtools::document()` clean; NAMESPACE unchanged.
4. Task-3 table complete (3 or 10 seeds, stated which and why).
5. Ridge assets render and are committed; scripts alongside.

## Deliverable

`/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms/dev/review-2026-08/reports/22-normal-slab.md`:
What was done (per task, commits) / The slab-coverage map (task 1) / OC
transfer table vs report 19 / The two evaluation artifacts with their
readings / Proposed NEWS entry VERBATIM / Findings (severity-tagged) / Open
questions. Copy to the Dropbox path; commit the report on the branch.

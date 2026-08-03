# Brief 23 — mixed relational certification: mixed vs GGM / OMRF (Opus agent)

You are the execution agent in the pre-CRAN review of the R package **bgms**.
You share no context with the review lead; everything you need is here.
**REPORT-ONLY batch: no source, test, or doc changes.** The mixed MRF is well
certified internally (likelihood/gradient hook, prior-only correction
identities, deliberate SBC scoping) but has never been checked AGAINST its
sibling models. That relation is this batch. If you find a divergence, you
characterize it — you do not fix it.

## Machine budget (standing rule — TWO other agents are running)

Cap your footprint at ~4 hardware threads; ONE model fit at a time; every
heavy step strictly sequential. The recovery study (task 4) is the heavy
tail — run it last, and size it so the whole batch's fitting stays around an
hour. Runtime is not a grading criterion.

## Setup

- Repo (Dropbox; do NOT switch its branch or build in it):
  `/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms`
- Base on **`origin/develop`** (NOT the local ref) AT OR AFTER `1393da20`:
  ```sh
  git -C /Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms \
      worktree add ~/bgms-review/wt-val23 -b review/mixed-relational-cert origin/develop
  cd ~/bgms-review/wt-val23
  ```
- Scratch outputs under `~/bgms-review/val23/` (fits, tables, seeds — keep
  them; the lead may re-derive). The ONLY tracked changes on your branch are
  the report and its assets. No pushing; no attribution trailers.
- Everything else is off-limits: `R/`, `src/`, `tests/`, `man/`,
  `vignettes/`, `NEWS.md`, workflows, and `dev/review-2026-08/` beyond your
  own report.

## Tasks — in this order

### 1. The routing map (decides what kind of claim the rest makes)

Read how `bgm()` dispatches to the OMRF / GGM / mixed machinery by
`variable_type`, with file:line citations. Answer precisely: can a USER
cause the mixed path to run on pure-type data (all-continuous or
all-ordinal), or does pure-type data always route to the pure sampler? If
the mixed path is user-reachable on pure data, tasks 2–3 are user-facing
consistency checks; if not, they are internal-machinery certification —
say which, in one sentence, at the top of the report.
To drive the mixed machinery directly where the API will not route to it,
mirror how the suite does it — `test-mixed-correction.R` builds
`sample_mixed_mrf` inputs by hand (`run_mixed_prior_chain()` shows the
shape).

### 2. All-ordinal reduction: mixed vs OMRF (the sharp check)

Same model, same inference type → the two paths must agree to Monte Carlo
error. Construct one all-ordinal dataset (modest: e.g. p = 6, 4 categories,
full support, n = 1000, seeded). Fit it through the OMRF path and through
the mixed machinery with IDENTICAL priors (document the mapping; the
default normal slab on both). Two comparisons:
a. `edge_selection = FALSE` (or its internal equivalent): posterior mean
   thresholds and pairwise associations, compared THROUGH THE PUBLIC
   EXTRACTORS on both sides so parameter frames align by construction.
b. Selection on (matched indicator priors): posterior inclusion
   probabilities.
Yardstick: refit the SAME data on ONE path with a different seed; that
run-to-run spread is the agreement tolerance. Report max |Δ| per block
against it. **If the paths diverge systematically — a stable slope away
from 1, a constant offset — STOP this task and characterize it fully**
(which parameter family, which direction, does it look like a factor):
that is a first-class finding of exactly the cross-path class this review
hunts.

### 3. All-continuous, calibrated: mixed vs GGM

These are NOT expected to agree exactly — the GGM path fits the exact
Gaussian likelihood, the mixed sampler a pseudolikelihood — so the check is
calibrated, not sharp:
a. One continuous dataset at n = 500 and one at n = 5000 (same truth,
   seeded). Fit both paths at both n.
b. Through the extractors, regress the mixed estimates on the GGM estimates
   across edges: report slope and intercept at each n. A pseudo-vs-exact
   finite-sample discrepancy SHRINKS as n grows; a convention bug (the
   omega-vs-2·omega class) is a slope away from 1 that is STABLE across n.
   Say which pattern you see.
c. Report max |Δ| per block at both n, against the same refit yardstick.

### 4. Cross-block recovery on genuinely mixed data (the no-sibling surface)

Cross-block edges (continuous–ordinal) have no sibling model to reduce to;
only planted truth certifies them. Design: 3 continuous + 3 ordinal
variables, a known truth with nonzero edges in ALL THREE blocks
(continuous–continuous, ordinal–ordinal, cross) and known zeros in each,
moderate n, 5 seeds, fits strictly one at a time. Report, per block family:
recovery slope (estimate vs truth), noise sd, detection of the planted
nonzero edges, false presences on the planted zeros. The cross block is the
result; the pure blocks anchor it against tasks 2–3.

### 5. OPTIONAL, only if the estimand map is clean — mgm overlap

An external cross-check against `mgm` was planned early and never run. Its
nodewise parameterization limits comparable estimands; if edge-detection
agreement on one mixed dataset is cleanly comparable, run it and report;
otherwise SKIP and say so explicitly (no silent caps).

## Verification gate

No suite runs (no code changed). The gate is internal to the report:
1. Every comparison carries its yardstick (refit spread), stated next to
   the number it judges.
2. Seeds, dimensions, and prior mappings stated for every fit.
3. Task 1's routing answer present with file:line citations.
4. Anything skipped or down-sized is named with its reason.

## Deliverable

`/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms/dev/review-2026-08/reports/23-mixed-relational-cert.md`:
Routing map / task 2 tables + verdict (agree to MC error: yes/no) / task 3
slopes at both n + the shrinking-vs-stable read / task 4 recovery table per
block / findings (severity-tagged; a systematic divergence is first-class)
/ Open questions. Copy to the Dropbox path as usual; branch committed, not
pushed.

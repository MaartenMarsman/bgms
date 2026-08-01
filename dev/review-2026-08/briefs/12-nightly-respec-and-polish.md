# Brief 12 — nightly respecification + release-polish batch (Opus agent)

You are the execution agent in the pre-CRAN review of the R package **bgms**.
You share no context with the review lead; everything you need is here. CODE
CHANGES authorized on a fix branch. Part A restructures the CI test tiers to a
maintainer-approved specification; Part B lands four small decided fixes.

## Setup

- Repo (Dropbox; do NOT switch its checked-out branch or build in it):
  `/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms`
- Base on `develop` AT OR AFTER `4969f843`:
  ```sh
  git -C /Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms \
      worktree add ~/bgms-review/wt-fix5 -b fix/nightly-respec develop
  cd ~/bgms-review/wt-fix5
  ```
  If brief 10's branch (`feat/edge-posterior-panel`) has merged by the time you
  finish, merge develop into your branch before the gates.
- One commit per item, `test:`/`ci:`/`fix:`/`docs:` prefixes, F-numbers
  referenced. No pushing; no attribution trailers.
- REQUIRED READING first:
  `dev/plans/backlog/2026-07-29_test-suite-tiering-audit_AUDIT.md` (untracked;
  read from the Dropbox path above). Its classifying principle — the CLASS OF
  BUG a test catches, not its runtime — governs every decision in Part A. Its
  T0 (~90 s every-run local split) already shipped as #185 (`c1e52c8e`); your
  job is the tier ABOVE it.

## Part A — the nightly respecification (F-071, closes F-044's branch gap)

Background: the current "nightly" tier (`BGMS_RUN_SLOW_TESTS=true`, Mon+Thu)
has outgrown CI — the 2026-08-01 dispatched run on develop was CANCELLED at
`timeout-minutes: 120` with no verdict (run 30709427193; locally the same
suite takes ~12 min on 15 cores). The maintainer's direction: "we specify and
condense the nightly tests" — a real respecification, not a timeout patch.

### The approved tier contract (implement exactly; flag disagreements)

| Tier | Cadence | Content (by bug class) | Budget | Enforcement |
|---|---|---|---|---|
| T0 every-run | local + (D1) PR | contracts, validation, unit guards, product surface | ~90 s local | shipped (#185) — do not touch |
| T1 nightly heartbeat | daily 03:00 UTC, **develop** | curated CALIBRATION subset: graph-law identities, prior-chain identities, gauge detector, surface-vs-gold single cells, RB saturation, concordance smokes — "does the settled math still hold tonight" | **≤ 60 min wall on the 2-core runner**; `timeout-minutes: 90` | `BGMS_RUN_SLOW_TESTS=true` |
| T2 weekly certification | Sunday 03:00 UTC, **develop** | the heavy Monte-Carlo machinery: SBC (1000-fit suites), full parameter-recovery sweeps, full mixed-nuts condition grids, n=2e6 MC channels, refit cross-validations (incl. the F-049 gate) | `timeout-minutes: 360` | `BGMS_RUN_CERTIFICATION=true` (new; T2 workflow sets BOTH vars) |

Scheduled runs move to **develop** (the branch where change happens); `main`
is verified manually at the release re-merge. The existing weekly-compliance
bitwise workflow is UNTOUCHED.

MECHANICS you must design around, not "fix": GitHub fires `schedule:` crons
from the DEFAULT branch's copy of the workflow file. So (a) the way a
scheduled run tests develop is an explicit `ref: develop` in the checkout
step, NOT the file's location; (b) the new schedules only go LIVE when this
batch reaches `main` at the release re-merge — until then, main's old
Mon+Thu nightly keeps firing and cancelling (known, pre-announced noise);
(c) `weekly-certification.yaml` is not even dispatchable until it exists on
main, so its budget is proven at its first real Sunday run post-re-merge —
state this in the report rather than working around it.

### Tasks

1. **Measure with real 2-core numbers.** Harvest per-file wall times from the
   cancelled run's log (`gh run view 30709427193 --log` — it executed ~2 h of
   files before the kill, which is most of the suite). For files it never
   reached, measure locally and scale by the observed local→CI ratio of files
   present in both (state the ratio and its spread).
2. **Classify every current slow-tier block** into T1 or T2 by the audit's
   bug-class principle FIRST, budget second. Produce the full classification
   table (file / block / class / CI-seconds / tier) in your report. Where
   principle and budget conflict, or a block is genuinely equivocal, put it in
   T2 and flag it for the maintainer. Mind the audit's shared-fixture warning
   (the cached Cauchy cell: co-located consumers must move together or the
   cost just relocates).
3. **Implement the gates**: `BGMS_RUN_CERTIFICATION` skips for T2 blocks
   (helper in the test setup mirroring the existing slow-tier gate), keeping
   `BGMS_RUN_SLOW_TESTS` as the T1 gate. A T2 block skipped under
   nightly-only env must say so in its skip message.
4. **Rewrite the workflows**: `nightly-validation.yaml` → explicit
   `ref: develop` checkout, daily cron, `timeout-minutes: 90`, T1 env only,
   informative job name; new `weekly-certification.yaml` → Sunday cron,
   `ref: develop`, `timeout-minutes: 360`, both env vars. Document the tier
   contract in a comment header in each file AND as a new subsection of
   `dev/review-2026-08/MAINTAINERS.md` §6.
5. **Prove T1's budget**: `workflow_dispatch` the rewritten nightly on your
   branch ref and report the wall time and totals. SCOPED EXCEPTION to the
   no-pushing rule, for this task only: you MAY push the fix branch itself
   (`git push origin fix/nightly-respec` — NEVER develop, never main),
   because dispatch needs the ref to exist on the remote; this works for
   `nightly-validation.yaml` since that file already exists on the default
   branch. Expected: completes well under 60 min, ZERO failures (Part B
   item 7 re-founds the one red test and moves it to T2). If it exceeds
   budget, re-cut and re-run — the budget is the spec.
6. **D1 — APPROVED (maintainer, 2026-08-01): the fast tier gates every
   push/PR.** Add `fast-checks.yaml`: on push to develop and on pull_request,
   run the ~90 s T0 (plain `devtools::test()`, no slow env vars) plus
   `R CMD check --no-manual --no-vignettes` if it fits the runner budget
   (drop the check half if it pushes the job past ~15 min — the test tier is
   the point). Push/PR triggers use the workflow file on the pushed branch,
   so this activates the moment it lands on develop — no re-merge wait.

## Part B — decided release polish

7. **F-049 gate re-founding** (maintainer-approved shape). Keep the test's
   construction; fix the fragile denominator: the noise estimate pools over
   FOUR independent refit pairs (deterministic seeds), gate stays
   `gap < 4 × pooled_noise`. Verify: the shipped construction passes on ≥ 3
   seeds including the previously-failing one, and the test carries a comment
   citing the 20-seed study (report 06) and the known ~0.01-pip bound.
   Classify this test T2. Plus the documentation half: one sentence in
   `?prior_sensitivity_check` Details — between-anchor curve values are
   importance-reweighted and can deviate up to ~0.01 in inclusion probability
   at the extrapolation ends; the anchors themselves (including the 1× anchor,
   the user's own fit) are exact.
8. **F-072 group labels.** Store `group_labels` (character) in the compare
   fit's arguments at build — the F-068 `baseline_category` precedent: one
   line in `build_arguments_compare()` + the field whitelist; labels =
   `as.character(unique_g)` (first-appearance order, matching the group
   numbering) on the indicator path, `c("x", "y")` on the x/y path. Then
   surface them on HUMAN displays only — extractor colnames (`group1`,
   `group2`) are consumed by easybgm/JASP and MUST NOT change:
   - `print`/`summary` compare headers gain one mapping line:
     "groups: 1 = en (n = 496), 2 = fr (n = 490)"
     (`R/methods_bgmcompare.R:185,212` vicinity);
   - plot group-panel titles "group 1 (en)" (`R/plot_bgms.R:509`);
   - calibration panel titles likewise (`R/calibration_check.R:828`);
   - centrality labels likewise (`R/centrality.R:173,178`).
   Old fits without the stored field degrade gracefully to bare numbers.
   Tests: label propagation on a Boredom fit (fr = 1 — the data is fr-first),
   x/y path, and a legacy-arguments fit.
9. **F-057 residue**: one Rd paragraph in `?bgmCompare` documenting how groups
   get their numbers (first appearance in `group_indicator`) and that
   original labels are carried into output displays (item 8).
10. **F-004**: the README development-install line gains the `@develop` ref
    (`install_github("Bayesian-Graphical-Modelling-Lab/bgms@develop")` form).

## Verification gate

1. Full LOCAL suite with BOTH env vars set (= everything): ZERO failures —
   the F-049 re-founding removes the last red. Zero warnings.
2. The dispatched T1 CI run from task 5: within budget, totals reported.
3. CRAN-mode suite: 0 failures. `R CMD check --as-cran` on a `git archive`
   tarball: the 2 baseline NOTEs only.
4. The classification table covers 100% of previously-slow blocks (none
   silently dropped — a block that vanishes from both tiers is a FAIL).

## Deliverable

Write
`/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms/dev/review-2026-08/reports/12-nightly-respec.md`:
What was done / The classification table + tier budgets measured / Findings
(severity-tagged) / Evidence (CI run links, gate outputs) / Open questions
(equivocal classifications flagged for the maintainer).

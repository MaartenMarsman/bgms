# Brief 05 — Test-repair and release-hygiene batch (Opus agent)

You are the execution agent in the pre-CRAN review of the R package **bgms**.
You share no context with the review lead; everything you need is here. CODE
CHANGES authorized on a fix branch. Deadline pressure: the twice-weekly nightly
(Mon 2026-08-03 03:00 UTC) will run the slow test tier against current code and
go red on the stale assertions below — this batch should land first.

## Setup

- Repo (Dropbox; do NOT switch its checked-out branch or build in it):
  `/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms`
- Base your branch on `develop` AT OR AFTER commit `adf87013` (the merged
  checking-layer fix batch — your work depends on it):
  ```sh
  git -C /Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms \
      worktree add ~/bgms-review/wt-fix2 -b fix/test-repair-batch develop
  cd ~/bgms-review/wt-fix2
  ```
- One commit per item, `test:`/`fix:`/`build:` prefixes, reference the
  F-number. No pushing; no attribution trailers. Regression discipline: every
  behavioral change gets a test that fails before, passes after.
- Context reports (readable via `git show develop:<path>` or on this worktree):
  `dev/review-2026-08/reports/03-checking-layer-defects.md` and
  `04-statistical-certification.md`.

## Items

### 1. F-048 — repair the three stale slow-tier assertions (test-only)

Report 04 §04-1 has the full diagnosis. The shipped surface-deployment fence is
`[0.5, 10]` (`R/zratio_surfaces.R:225-226`, set by PR #193), and three
assertions still encode the pre-#193 fence:

- `tests/testthat/test-zratio-surface-build.R:149-157` expects
  `zratio_build_surfaces()` to return NULL at shapes 2.5 and 5 — both INSIDE
  the shipped fence. Rewrite to assert the shipped contract: shapes inside
  [0.5, 10] (test at least 2.5 and 5) build non-NULL surfaces; shapes outside
  (test at least 12, and one below 0.5 if the code path allows) return NULL.
- `tests/testthat/test-zratio-gauge.R:213` ("a known-biased evidence-free fit
  fires the harm flag") — its fixture premise says shape 2 falls back to the
  additive kernel; under the shipped fence shape 2 deploys the surface, the
  projected distortion drops 5-fold (report 04 measured amplification 17.19,
  harm_pred 0.00794 < threshold 0.01), and the flag correctly stays down. Two
  acceptable repairs, your judgment which: (a) re-point the test at what it
  now demonstrates (the deployed kernel keeps a previously-harmful fixture
  under threshold — assert harm_pred < threshold and flag FALSE, with a
  comment explaining the kernel change), AND/OR (b) construct a fixture that
  is genuinely harmful under the DEPLOYED kernel so the flag's firing path
  stays tested. (b) is preferred if a harmful cell is reachable (consider a
  fence-outside shape where the additive kernel still serves); document
  whichever you do.
- DO NOT touch `tests/testthat/test-prior-sensitivity.R:432` (the 4×-noise
  gate failure) — that is finding F-049, a pending statistical decision for
  MM. Leave it failing in the slow tier and say so in your report.

### 2. F-050 — migrate the slow tier off deprecated argument names (test-only)

1152 deprecation warnings, all from four files (report 04 §04-4):
`test-sbc-ggm.R` (1000× `pairwise_scale`), `test-parameter-recovery-ggm.R`
(100×), `test-mixed-nuts.R` (19× `pairwise_scale` + 19× `main_alpha`),
`test-scaling-diagnostics.R` (9× + 5×). Replace with the 0.2.0 forms
(`interaction_prior = cauchy_prior(scale = ...)` / `normal_prior(...)`,
`threshold_prior = beta_prime_prior(...)`) preserving the numeric settings
exactly — these are certification suites; the tested cells must not move.
Acceptance: the slow tier emits ZERO deprecation warnings.

### 3. F-054 — version the correction-table disk cache key

`ggm_ctable_v1_...rds` omits the package version; the zratio surface cache
(`zratio_surf_v2_0.2.0.0_...`) includes it, and report 03 §4.2 observed real
cross-version contamination under the shared key (max |edens| diff 0.042
between rc1-built and develop-built tables). Add the package version to the
ctable key following the surface cache's convention (`R/correction_tables.R`,
key construction near the cache read/write). Old-key files are simply never
read again — do not migrate or delete them. Test: the key embeds
`packageVersion("bgms")`; a file under a wrong-version key is ignored.

### 4. F-055 — guard the q ≤ 3 hierarchical crash

`bgm(variable_type = "continuous", precision_graph_prior = "hierarchical")`
errors at q = 2 or 3: the bipartite anchor grid's smallest size is 4, so
`zratio_anchor_grids(cap ≤ 3)` filters to 0 rows and
`zratio_build_surfaces()`'s `bip_jobs$fam = "bip"` (`R/zratio_surfaces.R:411`)
fails on the 0-row frame (report 03 §4.1). Fix shape, per the review lead's
recommendation (flag in your report if anything in the code argues otherwise):
a 0-row bipartite job table yields a NULL bip family — bipartite bridges need
2+2 nodes, so the family is genuinely empty at q ≤ 3 and the engine's existing
NULL-family handling serves those fits through its other routes. Smoke tests:
q = 2, 3, 4, 5 hierarchical fits run to completion (small iter counts), and
q = 2, 3 produce no bip surface while q = 4, 5 do.

### 5. F-017 + hardening — `.Rbuildignore` lines

Add `^tests/compliance$` and `^tests/fixtures$` (keep the existing
`^tests/testthat/fixtures$` line), plus `^\.git$` (hardens tarball builds from
git worktrees, where `.git` is a file — report 03 §6). Acceptance: build the
tarball and list it — `tests/compliance/` (2.85 MB) and `tests/fixtures/` are
gone; nothing else changed vs the baseline listing.

### 6. F-051 — one docs line on the sensitivity plot's cap

`?plot.bgms_prior_sensitivity` documents the pip clamp; add the reading
consequence (report 04 §04-5): curve values are capped at log10 BF = 6, so
flat segments at 6.0 are display clamps, not evidence plateaus, and the
uncapped values live in `$edges$chosen_scale_log10_bf`; note that
`$edges$saturated` means zero-flip, not clamped. Keep it to a few sentences in
the existing Details; do NOT change any computation or unit (the log10 → nats
re-expression is a separate pending decision).

## Verification gate

From a clean install of the branch:
1. Changed test files pass.
2. **The full slow tier is green except the one documented F-049 failure**:
   `NOT_CRAN=true BGMS_RUN_SLOW_TESTS=true` full suite — report the totals;
   expected: 0 failures other than `test-prior-sensitivity.R:432`, 0
   deprecation warnings.
3. Full CRAN-mode suite: 0 failures.
4. `R CMD check --as-cran` on a `git archive` tarball: same 2 NOTEs as report
   01's baseline, and confirm the tarball listing for item 5.

## Deliverable

Write `/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms/dev/review-2026-08/reports/05-test-repair-batch.md`
(write it even if the checkout is on main — the lead collects it): What was
done per item (with commit SHAs) / Findings (anything new, severity-tagged) /
Evidence (gate outputs verbatim, tarball listing delta) / Open questions.

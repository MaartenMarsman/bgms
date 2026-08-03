# Brief 19 — the F-075 fix: category collapse in bgmCompare (Opus agent)

You are the execution agent in the pre-CRAN review of the R package **bgms**.
You share no context with the review lead; everything you need is here. CODE
CHANGES authorized on a fix branch. This batch fixes the review's RELEASE
GATE. Required reading, in this order:
`dev/review-2026-08/reports/16-f075-resolution.md` (the diagnosis — its §2
and §5 are the spec's foundation) and the maintainer's decisions below.

## The defect, in two sentences

`collapse_categories_across_groups()` (`R/validate_data.R:384-425`) keeps
only ordinal categories observed in EVERY group and silently folds the rest
downward, so groups that differ in location — the ordinary group-comparison
case — lose well-observed categories (a five-level variable became binary
on the F-075 data; its pairwise parameters then ran away). The defect
SHIPPED: `cran-0.1.6.3`'s `R/data_utils.R:388-394` has the identical loop,
so this is a true user-facing fix, not dev history.

## The maintainer's decisions (final; do not relitigate)

1. **Union semantics.** Ordinal variables retain the UNION of categories
   observed across groups. Values never observed in ANY group (true gaps)
   are still collapsed — that is the "unused categories" behaviour the Rd
   always claimed. A category unobserved in one group is a structural zero
   for that group: its group-specific threshold difference is data-free and
   sits at the prior. That is accepted and must be made VISIBLE (task 2).
2. **Blume–Capel stays exempt, deliberately** — BC scores are meaningful;
   collapsing/renumbering would destroy them. Document the rationale where
   the exemption is described.
3. **The compare log-posterior test hook lands in this batch** (task 6).

## Machine budget (standing rule — a second agent may be running)

~6 hardware threads total; ONE fit at a time (`cores = 4`); all heavy
steps sequential. The re-measure (task 7) is the only heavy part: 10
compare fits, run strictly one after another. Runtime is not graded.

## Setup

- Repo (Dropbox; do NOT switch its branch or build in it):
  `/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms`
- Base on **`origin/develop`** (NOT the local ref) AT OR AFTER `d7c2fe39`;
  verify report 16 is present in the tree, STOP if not:
  ```sh
  git -C /Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms \
      worktree add ~/bgms-review/wt-fix9 -b fix/category-collapse origin/develop
  cd ~/bgms-review/wt-fix9
  ```
- One commit per task or coherent group; `fix(compare):`/`test(compare):`/
  `docs(compare):` prefixes with F-numbers. No pushing; no attribution
  trailers. Build/test in the worktree or exports only.
- Do NOT touch plot files, NEWS entries other than your own additions, or
  anything another batch owns (`R/plot_*.R`, `R/calibration_check.R`,
  `R/centrality.R`, `R/prior_sensitivity.R`).

## Tasks

1. **The semantics fix (F-075).** Rewrite the ordinal branch of
   `collapse_categories_across_groups()`: recode onto the contiguous union
   of observed-in-any-group values; collapse only true gaps; BC branch
   untouched. Multi-group: union over all K groups.
   **Before flipping the switch, audit the downstream consumers.** The
   compare sampler was built against intersection-collapsed data; map every
   consumer of `num_categories` / the recoded matrix in
   `src/models/bgmCompare/` and the R-side prep, and check each survives a
   retained category with ZERO observations in some group: sufficient
   statistics / per-group count tables (zero cells), threshold start
   values (no `log(0)`/−Inf), adaptive proposal scales, and any indexing
   that assumed every category is observed everywhere. If something
   genuinely cannot handle a structural zero, STOP and report the options
   rather than forcing it — that finding would itself be first-class.
2. **Always-warn + record (F-110).** Two conditions, two volumes:
   `message()` when true gaps are renumbered (benign, informative);
   `warning()` when any RETAINED category has zero observations in some
   group — naming the variable, the category, the group, and the
   consequence in plain language (that group's threshold for the category
   is informed by the prior alone). Store the full recode map and the
   per-group support table in the fitted object's arguments so
   `summary()`/`print()` can surface them later. No jargon in either text.
3. **Documentation (F-111).** `man/bgmCompare.Rd` (the roxygen source):
   replace the false "unused categories are collapsed" sentence with the
   union semantics, the structural-zero consequence, and the BC-exemption
   rationale (decision 2). Matching note in `vignettes/comparison.Rmd`
   (one short paragraph; the vignette otherwise belongs to another batch's
   merged work — touch only your addition).
4. **NEWS.** One Bug-fixes entry, written for the 0.1.6.3 reader (the
   defect shipped there — cite nothing internal): what bgmCompare did
   (merged categories not observed in every group, silently), what it does
   now (retains all observed categories; warns when a group lacks support),
   and the consequence: any prior group comparison where groups differed in
   observed support should be refit. Re-run the markdown-NEWS parse gate
   after editing.
5. **Regression tests.** (a) The report-16 reproducer as a unit test of the
   NEW semantics: four-category variable, 40 obs/category, supports
   {0,1,2} vs {1,2,3} — union keeps four levels, nothing merged, the
   structural-zero warning fires; plus a true-gap case (a value observed in
   NO group collapses, message fires); plus a BC case (untouched).
   (b) Estimate survival: a SMALL group-differing-support simulation
   (choose dimensions so the fit runs in seconds at T0 scale, or gate it
   T1 if it genuinely needs more) asserting the recovered pairwise
   parameters match two matched `bgm()` fits within a derived tolerance —
   the assertion `test-collapse-categories.R` never had. Update that
   file's nine mechanics cases to the new semantics.
6. **The test hook (F-112).** Add `bgmCompare_test_logp_and_gradient()`
   as an `[[Rcpp::export]]` mirroring the GGM/mixed hooks
   (`ggm_test_logp_and_gradient` / `mixed_test_logp_and_gradient` — copy
   their interface shape). One unit test: compare's log-posterior and
   gradient against an independent R-side reference (finite-difference
   gradient check) on a tiny dataset, INCLUDING a group-differing-support
   case under the new semantics.
7. **Post-fix re-measure (the gate's second half).** Re-run report 16's
   operating-characteristics table on the SAME ten seeds. Reuse everything
   reusable from `~/bgms-review/val16/out/` (READ-ONLY): the truth objects
   and data seeds, and the saved `bgm()` baselines in `A_seed01..10.rds`
   (`bgm()` is untouched by this fix — do NOT refit it). Ten compare fits
   on the fixed build, one at a time. Report the table side by side with
   report 16's §5 defective-behaviour table: group slopes (expect ≈
   bgm()'s), noise sd, max null-difference error, δ = 0.40 recovery,
   detection, false presences — and the HEALTHY-EDGE slope explicitly: if
   the ~4–5% healthy-edge offset vanishes on fixed code, say so in its own
   sentence (it settles F-074's simulation side as collapse contamination).
8. **Real-data support scan (F-113).** Zero fits: for each reachable
   multi-group dataset — start with the ADHD comparison the docs site
   bakes (`~/Dropbox/Projecten/R/bgms-docs/fixtures/bake.R` shows where its
   data comes from; site repo READ-ONLY) and the shipped Boredom (known
   negative, include as the control row) — tabulate per-variable observed
   support by group and whether the new warning would fire. One table in
   the report: dataset, variables with differing support, would-warn.

## Verification gate

1. Full LOCAL suite, default tier AND `BGMS_RUN_SLOW_TESTS=true` (src/
   changed): 0 failures / 0 warnings each; list any test whose expectation
   you changed, with the semantics reason.
2. `R CMD check --as-cran` on a `git archive` tarball: the 2 baseline
   NOTEs only.
3. The task-5 tests pass on the fixed build AND the reproducer test FAILS
   on the pre-fix build (run it once against an unfixed export to prove
   the test bites — state the failure line).
4. The re-measure table complete, ten of ten seeds.

## Deliverable

`/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms/dev/review-2026-08/reports/19-category-collapse-fix.md`:
What was done (per task, commits) / The downstream-consumer audit (task 1
— what you checked and what each does with a zero-support cell) / The
re-measure table vs report 16's / Findings (severity-tagged; anything the
union semantics exposed) / Evidence (gate outputs, warning texts verbatim,
the pre-fix reproducer failure) / Open questions.

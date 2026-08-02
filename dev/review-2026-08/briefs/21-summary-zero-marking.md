# Brief 21 — mark prior-only rows in bgmCompare summaries (Opus agent, micro-batch)

You are the execution agent in the pre-CRAN review of the R package **bgms**.
You share no context with the review lead; everything you need is here. CODE
CHANGES authorized on a fix branch. One small, closed task: maintainer
decision F-117.

## Context, in three sentences

bgmCompare now retains every ordinal category any group observes (union
semantics, just merged). A retained category with ZERO observations in some
group makes that group's category-threshold difference prior-driven — a
classed warning fires at fit time and the per-group counts are stored on the
fit (`extract_arguments(fit)$category_support`) — but `summary()`/`print()`
still show those rows unmarked, and the reader of a saved fit never saw the
warning. The maintainer has decided: MARK the rows in the printed output.

## Machine budget

Trivial batch — another agent is running. ~4 threads, one seconds-scale fit
total, everything sequential.

## Setup

- Repo (Dropbox; do NOT switch its branch or build in it):
  `/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms`
- Base on **`origin/develop`** (NOT the local ref) AT OR AFTER `acc19428`:
  ```sh
  git -C /Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms \
      worktree add ~/bgms-review/wt-fix12 -b fix/summary-zero-marking origin/develop
  cd ~/bgms-review/wt-fix12
  ```
- One or two commits; `fix(compare):` prefix with F-117. No pushing; no
  attribution trailers.
- **HARD off-limits (another batch owns them):** every plot file and plot
  Rd/test, `tests/testthat/test-bgmCompare.R`,
  `tests/testthat/test-simulate-predict-regression.R`,
  `tests/testthat/test-methods.R`, anything under `src/`,
  `.github/workflows/`, `vignettes/`, `R/validate_data.R`, `NEWS.md`
  (report-only), `dev/review-2026-08/` except your report. Your surface is
  `R/methods_bgmcompare.R`, the Rd that documents compare summary/print
  output, and a NEW test file.

## The task

1. **Mark the rows.** Wherever compare main-effect difference rows are
   rendered (`R/methods_bgmcompare.R` — `summary()` builds `main_diff` at
   ~:144, `print()` shows it at ~:234), mark every row whose
   group-by-category cell in `category_support` has ZERO observations:
   an asterisk on the row (label or a trailing marker — pick what the table
   structure makes cleanest and say why in the report), plus ONE footnote
   line under the block:
   `* no observations in this group for this category; the estimate reflects the prior, not the data`
   The footnote prints whenever any rendered row carries the mark. Map rows
   to support cells structurally (indices used to build the rows), not by
   parsing row-name strings.
2. **Old fits degrade gracefully.** A fit whose arguments carry no
   `category_support` prints exactly as today — no mark, no footnote, no
   error.
3. **Rd.** One short paragraph in the Rd that documents the summary/print
   output, naming the mark and pointing to `category_support`. A cross-ref
   sentence in `man/bgmCompare.Rd`'s roxygen is allowed if natural; nothing
   more there.
4. **Tests, in a NEW file** `tests/testthat/test-summary-marking.R` (a new
   file keeps you clear of files the other batch edits): (a) snapshot of the
   marked summary on a seconds-scale adverse fit — four-level variable,
   supports {0,1,2} vs {1,2,3}, seeded, small iter; expect the classed
   `bgms_group_support_warning` explicitly rather than blanket-suppressing;
   (b) a shared-support fit prints unmarked and footnote-free; (c) the
   graceful old-fit case (strip `category_support` from the arguments).

## Verification gate

1. Full LOCAL default-tier suite: 0 failures / 0 warnings. No existing
   snapshot should change (existing snapshot fits share support); if one
   does, STOP and explain rather than re-record.
2. `devtools::document()` clean; NAMESPACE unchanged.

## Deliverable

`/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms/dev/review-2026-08/reports/21-summary-zero-marking.md`:
What was done (commits) / **BEFORE and AFTER `summary()` console output,
verbatim** — the maintainer judges the printed format from this report /
Evidence (gate outputs) / Open questions.

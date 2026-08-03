# Brief 14 — NEWS.md verification and reconstruction for 0.2.0.0 (Opus agent)

You are the execution agent in the pre-CRAN review of the R package **bgms**.
You share no context with the review lead; everything you need is here. CODE
CHANGES authorized on a fix branch touching **NEWS.md only** (plus your
report). This clears release blocker F-001 and the verification debt of
F-023/F-029.

## Machine budget (standing rule)

This machine has 15 cores and is SHARED (the maintainer's interactive work
plus other agents). Cap your total footprint at ~6 hardware threads; one
process at a time; sequence anything heavy. Runtime is not a grading
criterion. This brief needs essentially no compute — do not run the test
suite or R CMD check (you change no code); verification fits, if any, stay
tiny (p ≤ 5, seconds).

## Setup

- Repo (Dropbox; do NOT switch its checked-out branch or build in it):
  `/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms`
- Base on `develop` AT OR AFTER `5b850410`:
  ```sh
  git -C /Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms \
      worktree add ~/bgms-review/wt-fix6 -b fix/news-0.2.0.0 develop
  cd ~/bgms-review/wt-fix6
  ```
- One commit (or few, `docs:` prefix). No pushing; no attribution trailers.
- REQUIRED READING first: `dev/review-2026-08/FINDINGS.md` rows **F-001,
  F-023, F-029** (your mandate), and skim `dev/review-2026-08/reports/00c-*`
  (the R/API diff map — it lists the 23 new exports and the default changes).

## Context

CRAN has **0.1.6.3**; this release is **0.2.0.0** (233 commits, tag
`cran-0.1.6.3` marks the anchor). NEWS.md already carries substantial
accumulated text, but its top section is framed as "Changes since the
0.2.0.0 development build" — written for people tracking develop. The
released file must speak to a **0.1.6.3 user**: one coherent
`# bgms 0.2.0.0` section describing the change from the CRAN version.
Accuracy law: a NEWS claim that misdescribes the code is a release defect
(that is exactly what F-001 is — the prior-sensitivity entry misdescribes
the mechanism).

## Tasks

1. **Inventory.** `git log --oneline cran-0.1.6.3..HEAD`, the merged-PR
   titles it contains, the 23 new exports (report 00c), and the F-029 list
   of silent behavior changes. This is your completeness checklist.
2. **Verify every existing NEWS claim** against the code at HEAD. Produce a
   claims table in your report: entry / verified-against (file:line or a
   snippet you ran) / verdict (accurate, corrected, removed). For F-001
   specifically: reread `man/prior_sensitivity_check.Rd` and
   `R/anchor_curve.R` (`assemble_curve()`), then fix the sensitivity-check
   entry so it describes the actual mechanism (anchored refits + importance
   reweighting between anchors — the anchors exact, the between-anchor curve
   reweighted).
3. **Restructure for the 0.1.6.3 reader.** One `# bgms 0.2.0.0` section,
   ordered by user impact:
   - **Breaking changes** first. At minimum: the bgmCompare pairwise-scale
     semantic fix (old compare fits used `omega·x` where every other path
     used `2·omega·x`; stored fits should be refit — state the user action);
     `extract_ess()` indicator semantics now Rao-Blackwellized;
     `simulate_mrf()`'s factor-2 input-scale change (the F-023 sentence: a
     user-supplied `pairwise` matrix means something different now);
     defaults: `iter`/`warmup` 1e3→2e3, NUTS `target_accept` → 0.80,
     `interaction_prior` → `normal_prior(1)` in `bgm()`/`sample_ggm_prior()`
     (bgmCompare keeps Cauchy — say so); `hamiltonian-mc` and
     `hmc_num_leapfrogs` removed; `standardize` deprecated.
   - **New features**: GGM and mixed MRF model classes, the bgmCompare
     rewrite, RB inclusion machinery, `prior_sensitivity_check()`,
     `calibration_check()` (incl. the bgmCompare method), the hierarchical
     precision-graph prior + trust gauge, the plot redesign, new extractors.
   - **Fixes** and **Deprecations** last.
   Every F-029 item gets verified line-by-line coverage — a silent behavior
   change with no NEWS sentence is a FAIL of this brief.
4. **Every entry answers three questions** for a user: what changed, who is
   affected, what (if anything) they should do. Natural-log convention
   throughout ("log BF" means ln). No marketing tone.
5. **Nothing invented.** Every entry maps to a commit, PR, or FINDINGS row;
   the mapping table goes in your report. Where you cannot establish an
   entry's truth from code + record, put it in Open questions rather than
   guessing — the maintainer authored the intent and reviews your draft.

## Verification gate

1. The claims table covers 100% of pre-existing NEWS entries.
2. The F-029 checklist: every listed change has an entry, cited.
3. `tools:::news2Rd` (or `utils::news()` on the installed structure) parses
   the file without error — run against the file alone, no package build.

## Deliverable

Write
`/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms/dev/review-2026-08/reports/14-news-verification.md`:
What was done / Findings (severity-tagged — any pre-existing NEWS claim
found FALSE is a finding) / Evidence (the claims + mapping tables) / Open
questions (entries needing the maintainer's intent).

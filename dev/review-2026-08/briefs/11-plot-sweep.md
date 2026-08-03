# Brief 11 — package-wide plot restyle + compare evidence parity (Opus agent)

You are the execution agent in the pre-CRAN review of the R package **bgms**.
You share no context with the review lead; everything you need is here. CODE
CHANGES authorized on a fix branch. This brief finishes the plot program a
previous brief started: it restyled the single-edge posterior panel onto a
shared style module (`R/plot_style.R`); you now carry that language across
every remaining plot and close the maintainer's ratified change requests.

The maintainer will judge BEFORE/AFTER renders of EVERY figure before any
of this merges — that is a standing requirement at his request. Your job is
to make those renders and the code behind them; the merge decision is his.

## Machine budget (standing rule — you are the THIRD concurrent agent)

This machine has 15 cores, SHARED: one agent is running a long fit batch,
another a docs+check batch. Cap your total footprint at ~4 hardware
threads. Everything sequential: one fit at a time, one suite run at a time.
Reuse the smallest viable fits (the test fixtures in
`tests/testthat/helper-fixtures.R`, or shipped-data fits at low `iter`) and
seed every fit so renders are reproducible. Runtime is not a grading
criterion.

## Setup

- Repo (Dropbox; do NOT switch its checked-out branch or build in it):
  `/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms`
- **Base on `develop` AT OR AFTER `f24ad3c8`** (the brief-12 merge — it
  rewrote the compare-panel titles you will restyle; verify with
  `git log --oneline -5 develop`, STOP if absent):
  ```sh
  git -C /Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms \
      worktree add ~/bgms-review/wt-fix8 -b fix/plot-sweep develop
  cd ~/bgms-review/wt-fix8
  ```
- One commit per coherent item, `feat(plot):`/`fix(plot):`/`style(plot):`
  prefixes with F-numbers. No pushing; no attribution trailers.

## HARD parallel-safety constraints (another agent's batch is in flight)

DO NOT TOUCH: `NEWS.md`, anything under `vignettes/`, `R/verdicts.R`,
anything under `src/`, `.github/workflows/`, `test-mcmc-diagnostics.R`.
Where your changes deserve a NEWS clause (they do — at least the compare
evidence display), write the proposed clause VERBATIM in your report; the
lead lands NEWS at integration. In `R/prior_sensitivity.R` you own ONLY
`plot.bgms_prior_sensitivity` (`:1245` on) — do not edit the roxygen
Details block near the top of the file (another batch owns its wording).

## The visual language (ratified; do not relitigate)

JASP/compendium feel, distilled not ported (base graphics here; JASP's
`jaspGraphs` holds their working constants; the R Graph Compendium is the
broader reference). Concretely, via `R/plot_style.R` helpers everywhere:
typography and muted ink from the module; offset axes; no plot box; NO big
main titles; wheels/rings carry probabilities (`probability_wheel()`,
`plot_style.R:214`); Bayes-factor text is natural log ("log BF"), never
log10; muted single-line captions where a plot needs context. The
maintainer owns this direction explicitly ("I am happy with the drift. I
own it") — deviations from older docs' presentation rules are DECIDED.

## Tasks

1. **(F-067) Restyle every remaining plot onto the style module.** The five
   entry points: `plot.bgms` (`R/plot_bgms.R:128`), `plot.bgmCompare`
   (`R/plot_bgms.R:328`), `plot.bgms_calibration`
   (`R/calibration_check.R:764`), `plot.bgms_centrality`
   (`R/centrality.R:307`), `plot.bgms_prior_sensitivity`
   (`R/prior_sensitivity.R:1245`). The single-edge panel
   (`plot_edge_posterior`) is already on the module — it is your reference,
   not your target. Remove headline `main` titles (the sensitivity plot's
   big title is the explicitly named offender); carry any needed context in
   the module's caption/subtitle idiom.
2. **(MM ratified) The edge panel drops the verdict-word annotation in
   every case** — no "evidence of absence", no "undecided" wording on the
   panel; the numbers and the picture carry it. (The network legends at
   `plot_bgms.R:187` and `:445` also print "undecided" as a line-type key:
   do NOT change those unilaterally, but flag the consistency question in
   your report and render both variants of ONE network figure so the
   maintainer can choose with his eyes.)
3. **(MM ratified) "PIP" leaves the rendered text.** Every printed "PIP"
   (e.g. `plot_bgms.R:879`, composed per `plot_style.R:335`) becomes a
   self-explanatory inclusion label — lead proposal `"P(included)"`; keep
   it one short token, same composition idiom with the relation operator.
4. **(MM ratified) Three decimals via the style constant**:
   `format_probability()` (`plot_style.R:342-353`) moves from two decimals
   to three, no-leading-zero style preserved, ALL call sites move with the
   constant. Leave `format_log_bf()`'s precision as is unless the renders
   argue otherwise — if you change it, show before/after and say why.
5. **(F-079, MM requirement) `plot.bgmCompare` shows difference Bayes
   factors with the SAME evidence conventions `plot.bgms` uses.** Read
   `plot.bgms`'s evidence display first (all its `type`s), then give the
   compare method the mirror-image display for difference evidence: same
   thresholds (natural-log, ln 10 presence), same formatting helpers, same
   visual grammar. This is the network/overview level; a single-DIFFERENCE-
   edge panel is explicitly OUT (recorded post-release item) — do not build
   it. In your report, state in one paragraph exactly what "parity" ended
   up meaning in code.
6. **(report-10 deviation 8) Compare node rings: qgraph pie →
   `probability_wheel()`** — the compare network's per-node qgraph `pie`
   ring is replaced by the module's wheel so probability marks look the
   same on every figure in the package.
7. **(F-070) Sensitivity-plot right-margin edge labels clip at the device
   edge at default width** (pre-existing, visible in old and new renders).
   Fix properly: measure the widest label and reserve the margin, or place
   labels inside — no hard-coded magic width. Show the fix at default
   device size.

## Renders for the maintainer (the deliverable that gates the merge)

For EVERY figure the package can draw — the five methods (each meaningful
`type`/variant) plus the single-edge panel in its presence / absence /
undecided cases — write a BEFORE (develop) and AFTER (your branch) PNG at
identical size and data to
`dev/review-2026-08/renders/11/<figure>-{before,after}.png`, seeded so they
re-render byte-stable. Add `renders/11/INDEX.md`: one row per figure —
what changed and which task did it. Small, honest fits are fine; the
figures must show real structure (some presence, some absence, some
undecided; for compare, a real difference and a null).

## Verification gate

1. Full LOCAL default-tier suite (no slow env vars), sequential: 0
   failures, 0 warnings. List every plot-related snapshot you re-record
   with a one-line reason; touch no other snapshots.
2. `devtools::document()` clean; no unrelated NAMESPACE/Rd drift.
3. Every render pair regenerates from a fresh session via a single script
   you commit at `dev/review-2026-08/renders/11/make_renders.R`.
4. Roxygen for the five methods matches the new behaviour (no stale
   "PIP"/verdict-word/title wording left in their Rd text).

## Deliverable

Write
`/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms/dev/review-2026-08/reports/11-plot-sweep.md`:
What was done (per task, with commits) / The parity paragraph (task 5) /
Proposed NEWS clauses VERBATIM / Findings (anything the restyle surfaced
in plot logic — severity-tagged) / Evidence (gate outputs, snapshot
re-record list) / Open questions (incl. the legend-word consistency
question from task 2).

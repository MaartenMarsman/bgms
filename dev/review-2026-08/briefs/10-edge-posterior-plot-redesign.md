# Brief 10 — edge-posterior panel redesign, JASP standard (Opus agent)

You are the execution agent in the pre-CRAN review of the R package **bgms**.
You share no context with the review lead; everything you need is here. CODE
CHANGES authorized on a fix branch. This is a maintainer-requested redesign of
`plot_edge_posterior()` to publication quality — his reference is the JASP
prior-and-posterior figure in the R Graph Compendium:
https://www.shinyapps.org/apps/RGraphCompendium/index.php#prior-and-posterior
(base-R code is on the page). His verdict on the shipped panel: "Our plot
really is not **that** good." The target: a figure a reader could lift into a
research paper.

## Setup

- Repo (Dropbox; do NOT switch its checked-out branch or build in it):
  `/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms`
- Base on `develop` AT OR AFTER the brief-06 merge (this brief deliberately
  waits for it — both touch `R/plot_bgms.R`):
  ```sh
  git -C /Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms \
      worktree add ~/bgms-review/wt-fix4 -b feat/edge-posterior-panel develop
  cd ~/bgms-review/wt-fix4
  ```
- Conventional commits, F-066 referenced. No pushing; no attribution trailers.
- Current implementation: `R/plot_bgms.R`, `plot_edge_posterior()` (~line 535)
  — a stem at zero (posterior spike probability) plus a slab density scaled to
  probability-per-bin. Honest, but no prior, no interval, unnumbered y axis.

## The design

One panel, base graphics only (no new dependencies), compendium conventions
throughout: large fonts, offset axes, no box, needless ink omitted.

1. **Densities, conditional framing.** Solid accent curve with a light fill =
   the posterior density of the edge weight CONDITIONAL on inclusion (a real
   density integrating to 1 — the y axis gets its numbers back, labeled
   "Density"). Dashed grey curve = the slab prior density, computed
   ANALYTICALLY from the fit's own interaction prior (read the spec:
   `cauchy_prior`/`normal_prior` and its scale from `extract_arguments()`; do
   not hardcode the default, which changed in 0.2.0). Legend "Posterior /
   Prior" top-left, compendium style.
2. **The inclusion mass moves to a probability wheel.** Top corner inset: a
   small wheel (pie) filled by the POSTERIOR inclusion probability — accent
   share = included, pale share = excluded — with "PIP = .87"-style text
   beside it. Draw it with `cos()`/`sin()` + `polygon()` as the compendium
   does. This deliberately matches the package convention forming in the
   compare plots (F-062d: node rings filled by pip), one visual language:
   wheels/rings carry probabilities, text carries log BF.
3. **Evidence text.** Next to the wheel: "log BF = 4.31" via the existing
   `format_log_bf()` (natural log, capped — the panel title convention from
   `edge_panel_title()`; keep the title as is or simplify now that the wheel
   carries part of the story, your judgment).
4. **Estimate text.** Top right, JASP style: "median = 0.31" and
   "95% CI [0.12, 0.48]" of the CONDITIONAL posterior (quantiles of the
   nonzero draws).
5. **One honest caption line** (small, muted, replacing the current mtext):
   the slab is shown conditional on inclusion; absence probability is the
   wheel's pale share.

**The dots rule (maintainer-set): the panel shows whichever estimator the
model licenses.** The JASP figure reads the BF off the prior/posterior
ordinates at zero (Savage-Dickey). Whether that transfers depends on the fit:

- `edge_selection = TRUE` (default): bgms BFs are Rao-Blackwellized INDICATOR
  BFs from the spike-and-slab, not density ratios — Savage-Dickey dots would
  visually assert an estimator the package does not use. NO dots; the PIP
  wheel + printed RB log BF carry the evidence.
- `edge_selection = FALSE`: no indicator exists, the posterior is continuous,
  and Savage-Dickey IS the licensed estimator. Draw the JASP figure exactly:
  both grey ordinate dots at zero; BF = prior ordinate (analytic, from the
  fit's interaction prior) / posterior ordinate (from the draws — use
  `stats::density(bw = "SJ")` interpolated at 0; JASP's own implementations
  use logspline, but do NOT add a dependency — document the estimator choice
  in the Rd); print it as natural-log "log BF" via `format_log_bf()`; wheel
  filled by BF/(1+BF) with JASP's data|H1 / data|H0 labeling (that fill is a
  posterior probability at equal prior odds, so the wheels-carry-
  probabilities convention holds); median/CI as usual.

## Edge cases (all tested)

- PIP ≈ 1 (saturated): wheel effectively full, log BF prints as the capped
  form ("log BF > ..."); CI/median unaffected.
- PIP ≈ 0 / fewer than 2 nonzero draws: today this `stop()`s. Instead: draw
  the PRIOR density, the wheel (nearly all pale), and the log BF, with the
  caption noting no included draws — a decisive-absence edge deserves a
  figure too, not an error.
- Fits without edge selection: the Savage-Dickey case in the dots rule above
  — full posterior density, prior, both ordinate dots, SD log BF, BF-filled
  wheel, median/CI. Snapshot this case specifically, including one where the
  posterior ordinate at 0 is near-zero (decisive) and one where it exceeds
  the prior ordinate (evidence for absence).
- Blume-Capel / mixed fits: weights are continuous; nothing special, but
  include one in the snapshots.

## API

- `binwidth` loses its meaning under the conditional framing: lifecycle-
  deprecate it (warn-and-ignore, following the package's deprecation
  pattern), do not silently drop it.
- No other signature changes; `evidence_threshold` keeps feeding the verdict
  in the title.

## The style becomes a module (this panel is the reference implementation)

The maintainer has decided the JASP/compendium feel goes PACKAGE-WIDE (brief
11 will restyle every other plot). So do not inline the styling: extract it as
roxygen-documented internal helpers — a new `R/plot_style.R` — written as if
they are the package's plotting law:

- `bgms_panel_par()` — typography, margins, offset-axis (eps) conventions,
  no box, large fonts, ink/muted/accent colors from `mover_palette()`.
- `probability_wheel(x, y, prob, radius, labels)` — the `cos()`/`sin()` +
  `polygon()` wheel, usable by this panel, the compare node rings (F-062d),
  and anything else that shows a probability.
- `annotation_block()` — the median/CI/log-BF text placement.
- A NO-BIG-TITLES convention: JASP plots carry annotations, not headline
  titles. State it in the helper docs; this panel may keep only a compact
  title (or none — your judgment against the reference figures).

Consult the compendium beyond the one figure, and the JASP source
(github.com/jasp-stats — `jaspGraphs` encodes their working constants: font
sizes, axis break logic, expansion factors) for the conventions; distill,
don't port (they are ggplot2, we are base graphics).

## Verification gate

1. Snapshot tests in the existing `_snaps/plot-methods` pattern: a decisive
   edge, an undecided edge, a saturated edge, a decisive-absence edge, a
   no-selection fit, one non-ordinal fit.
2. Visual proof for the maintainer: render BEFORE/AFTER PNGs of the same
   three Wenchuan edges (decisive / undecided / absent) to
   `dev/review-2026-08/reports/assets/edge-panel-{before,after}-*.png`.
3. Full CRAN-mode suite: 0 failures; `R CMD check --as-cran` on a
   `git archive` tarball: the 2 baseline NOTEs only.
4. The Rd example still runs; Rd regenerated.

## Deliverable

Write
`/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms/dev/review-2026-08/reports/10-edge-posterior-plot.md`:
What was done / Findings (severity-tagged) / Evidence (the before/after PNGs,
gate outputs) / Open questions. The maintainer reviews the AFTER figures
before this merges — flag anywhere you deviated from the design and why.

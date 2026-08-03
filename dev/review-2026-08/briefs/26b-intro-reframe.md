# Brief 26b — reframe the vignette set under the maintainer's terminology rulings (Opus agent)

You are the execution agent in the pre-CRAN review of the R package **bgms**.
You share no context with the review lead; everything you need is here. This
is the REMEDIATION round of brief 26: the first intro rewrite was executed
cleanly (gates green) but the maintainer REJECTED the rendered result on its
frame and terminology. His verdict, verbatim (2026-08-03):

> "The vignettes are misaligned with the documentation, tutorial, and website,
> in that they only consider the discrete variable model a Markov random
> field. ALL MODELS IN BGMS ARE MARKOV RANDOM FIELDS. Also we DO NOT CALL THEM
> network models BUT graphical models. Continuing reading the material, I am
> unhappy. Why discuss problems with regularization. This package is not about
> them. The idea for spike and slab priors will never land. Really, this
> material is terrible."

He then chose FIX over revert. Your job is the fix.

## The four rulings — HARD CONSTRAINTS on every sentence you write

1. **Every model in bgms is a Markov random field.** The Gaussian graphical
   model is the Gaussian MEMBER of that one family. Never present "MRFs for
   discrete data" beside "GGMs for continuous data" as different kinds of
   thing.
2. **The models are GRAPHICAL models. Never "network models."**
3. **No regularization-contrast framing.** Do not motivate the package by the
   problems of penalized/regularized estimation. The package stands on what it
   does.
4. **Do not lead with spike-and-slab prior machinery.** The front door leads
   with what the user gets — evidence per edge, the three verdicts. The
   spike-and-slab prior is named once, in the Priors section, where priors are
   the topic.

The canonical frame already exists in the package and tutorial — inherit it:

- `DESCRIPTION`: "Bayesian Analysis of Graphical Models"; "the variable types
  determine the model: an ordinal Markov random field for discrete data, a
  Gaussian graphical model for continuous data, or a mixed Markov random field
  combining both."
- Tutorial manuscript (read-only,
  `~/Library/CloudStorage/Dropbox/projecten/SV/BayesianNetworkTutorial/manuscript.tex:222`):
  "They are called Markov random fields, and their central idea is simple: Two
  variables are connected by an edge when they remain associated after we
  account for every other variable in the network."

**The maintainer-approved opening.** The vignette MUST open with this text
(verbatim up to typographic details; if you believe a word must change, keep
the change minimal and flag it in the report):

> **bgms** is for Bayesian analysis of graphical models: Markov random
> fields, in which two variables are connected by an edge when they remain
> associated after accounting for every other variable in the model. Every
> model in the package is a Markov random field — `variable_type` selects the
> member: the ordinal Markov random field (with binary variables as a special
> case), the Blume–Capel model, the Gaussian graphical model for continuous
> variables, and the mixed Markov random field that joins discrete and
> continuous variables in one model.
>
> The graph is treated as unknown, not fixed: the analysis returns a
> posterior distribution over graphs, and for every pair of variables an
> inclusion Bayes factor — the factor by which the data shift the odds that
> the edge is present. Because a Bayes factor can point both ways, every pair
> ends in one of three states — evidence of presence, evidence of absence, or
> undecided — and `verdicts()` reports which.

On "network" as a plain noun: the rulings govern what the MODELS are called.
The plot API's own `type = "network"` argument and established phrases for
the drawn object are not banned — but in prose PREFER "graph" or "model",
and your report must list EVERY retained "network" instance across all five
vignettes with its sentence, marked keep/changed, for the maintainer's eye.

## Machine budget (standing rule)

~4 hardware threads; ONE fit/knit at a time; check `ps aux | grep "[e]xec/R"`
first and wait out any long job. Another review agent may be active — do not
touch `dev/review-2026-08/reports/25*` in any tree.

## Setup

- Repo (Dropbox; do NOT switch its branch or build in it):
  `/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms`
- Base on **`origin/develop`** (NOT the local ref) AT OR AFTER `d9b5a4dd`:
  ```sh
  git -C /Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms \
      worktree add ~/bgms-review/wt-fix17 -b fix/intro-vignette-2 origin/develop
  cd ~/bgms-review/wt-fix17
  ```
- Starting material for intro: the PARKED round-1 branch (do not modify it):
  `git show fix/intro-vignette:vignettes/intro.Rmd > vignettes/intro.Rmd`
  — its machinery (nine-item worked example, defaults table, verdicts/plot
  flow, GGM/mixed and bgmCompare sections, bibliography header) SURVIVES;
  its frame and terminology do not.
- Your surface: `vignettes/intro.Rmd`, `vignettes/checking-your-model.Rmd`,
  `vignettes/comparison.Rmd`, plus terminology-only edits in
  `vignettes/diagnostics.Rmd` / `vignettes/prior-sensitivity.Rmd` if your
  sweep (task 4) finds violations there. EVERYTHING else is off-limits:
  `R/`, `src/`, `tests/`, `man/`, `NEWS.md` (verbatim proposals in the
  report; the lead lands it), `vignettes/refs.bib`, workflows,
  `dev/review-2026-08/` except your report+assets.
- Small commits by concern; `docs(vignette):` prefix with the finding id
  (F-026 / F-127 / F-128). No pushing; no attribution trailers.

## Tasks

### 1. Reframe `intro.Rmd` (F-026)

The approved opening replaces the old "What the package answers" +
regularization paragraph. Then bring every remaining sentence under the four
rulings: no "network model(s)"; no discrete-MRF-vs-GGM split (the models
section presents ONE family and its members); the regularization contrast
deleted; spike-and-slab named exactly once, in Priors. Keep the round-1
machinery listed above unless a ruling forces a change. Keep the show-uncapped
/ run-capped chunk pair, `fig.width = 12` for the plot, the bibliography
header, and the website/easybgm/NEWS pointers.

### 2. `checking-your-model.Rmd` (F-127 — release-visible)

(a) The plot chunk at line ~80 (`fig.height = 5.5, fig.width = 6`) →
`fig.width = 12, fig.height = 4.5, out.width = "100%"`. The shipped HTML
currently prints `plot.bgms`'s narrow-device advisory INTO the vignette and
the three panel titles collide; your rebuilt HTML must show that chunk's
output clean (no advisory) — quote it in the report.
(b) The paragraph beneath (line ~85) still describes the pre-0.2.0.0
single-network display ("dotted grey edges are undecided; edges with evidence
of absence are not drawn at all"). Rewrite it to describe the three-panel
display, in ruled terminology.
(c) Line ~205: "a network model" → ruled terminology. Then read the WHOLE
file for the same two diseases (old-display descriptions, model naming) and
fix what you find, listing each change.

### 3. `comparison.Rmd` (F-128)

Replace the hand-rolled `qgraph` block built from
`coef(fit)$pairwise_effects_groups` with `plot(fit, type = "groups")` (the
supported display: each group's own network, shared layout; see
`?plot.bgmCompare`). Add the bibliography header the other vignettes carry
(`bibliography: refs.bib`, `csl: apa.csl`, `link-citations: TRUE`) — cite
only keys already in `refs.bib`; the .bib file itself is off-limits. Then a
full terminology pass over the file.

### 4. Terminology sweep, all five vignettes

`grep -in "network"` across `vignettes/*.Rmd`: fix every model-naming
violation wherever it sits; produce the keep/changed table for every
remaining instance. Gate: `grep -ci "network model" vignettes/*.Rmd` = 0
everywhere.

### 5. NEWS proposals (report-only)

Rewrite round 1's proposed intro line under the rulings (it belongs under
`## Other changes` — the vignette existed at 0.1.6.3, tag-visible). State
whether the `checking-your-model.Rmd` and `comparison.Rmd` changes need NEWS
touches, given each vignette's tag status (check whether each existed at
`cran-0.1.6.3`), and propose exact wording for any that do. The lead lands
all NEWS.

## Verification gate

1. `R CMD check --as-cran` on a `git archive` tarball: 2 baseline NOTEs only
   (stale Date + local HTML Tidy; a third examples-timing NOTE can appear
   under `--run-donttest` — known, not yours).
2. Full default-tier suite: 0 failures / 0 warnings (nothing should move —
   say so).
3. Terminology gates: "network model" count 0; the retained-"network" table
   complete; the approved opening present verbatim.
4. All touched vignettes' rendered HTML read end to end — every chunk runs,
   nothing leaks; the F-127 advisory is GONE; intro knit time measured
   (budget ~30 s) and the full-set vignette rebuild time from the check
   reported.

## Deliverable

`/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms/dev/review-2026-08/reports/26b-intro-reframe.md`:
per-file change map / the retained-"network" table / verbatim clean plot
output for F-127 / knit + rebuild times / NEWS proposals VERBATIM / findings /
open questions. **Copy the rendered HTML of every touched vignette from the
built tarball's `inst/doc/` to
`dev/review-2026-08/reports/assets/f026b-<vignette>.html`** — the maintainer
judges the renders before merge, same as round 1. Copy the report to the
Dropbox path; commit report + assets on the branch.

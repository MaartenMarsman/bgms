# Brief 26 — rewrite the front-door vignette against the 0.2.0.0 package (Opus agent)

You are the execution agent in the pre-CRAN review of the R package **bgms**.
You share no context with the review lead; everything you need is here. CODE
CHANGES authorized on a fix branch, scoped to ONE file plus its build assets.
Maintainer decision (F-026, 2026-08-03, verbatim): "Do it now. The package
has changed dramatically, it needs to align with the package, website, and
tutorial."

## Machine budget (standing rule — ONE other agent may be running)

~4 hardware threads; ONE fit at a time; check `ps aux | grep "[e]xec/R"`
before anything that fits or knits and wait if a long job is mid-flight.
Your own fits are seconds-scale; the vignette build must stay light.

## Setup

- Repo (Dropbox; do NOT switch its branch or build in it):
  `/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms`
- Base on **`origin/develop`** (NOT the local ref) AT OR AFTER `f0b07fd5`:
  ```sh
  git -C /Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms \
      worktree add ~/bgms-review/wt-fix16 -b fix/intro-vignette origin/develop
  cd ~/bgms-review/wt-fix16
  ```
- Your surface: `vignettes/intro.Rmd`, optionally one small baked `.rds`
  under `vignettes/` with its generator under `dev/vignette-data/` (the
  prior-sensitivity vignette set that precedent). EVERYTHING else is
  off-limits: `R/`, `src/`, `tests/`, `man/`, other vignettes, `NEWS.md`
  (verbatim proposal in the report; the lead lands it), workflows,
  `dev/review-2026-08/` except your report.
- One or two commits; `docs(vignette):` prefix with F-026. No pushing; no
  attribution trailers.

## Alignment sources — read before writing

1. **The package is authoritative** — current roxygen/Rd and NEWS.md are the
   record of what 0.2.0.0 is. Where the old vignette contradicts them, the
   package wins. The maintainer is the design authority over all of it.
2. The website strategy (READ-ONLY, outside the repo):
   `/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/plans/tutorial-website/bgms-website-strategy.md`
   — the vignette is the front door; the site carries the deep teaching. Do
   not duplicate the site's job; point to it where it exists.
3. The tutorial repo (READ-ONLY, tone/terminology reference):
   `~/Library/CloudStorage/Dropbox/projecten/SV/BayesianNetworkTutorial`
   — match its terminology for shared concepts (inclusion Bayes factors,
   evidence categories); do not import its length.

## What the rewrite must cover (the old file predates ALL of this)

- The prior constructor API (`normal_prior()`, `cauchy_prior()`,
  `beta_prime_prior()`, edge priors) and the CURRENT defaults — Normal slab
  in both `bgm()` and `bgmCompare()`, `precision_graph_prior =
  "hierarchical"` for continuous data (one honest sentence on the trust
  gauge; details belong to the docs).
- What the package fits now: ordinal MRF, GGM (continuous), mixed MRF,
  Blume-Capel — one paragraph each at most, with `variable_type` doing the
  routing.
- Inference the package's way: `verdicts()`, inclusion Bayes factors on the
  natural-log display convention, the evidence/absence/undecided
  trichotomy, `summary()`.
- The three-panel edge evidence plot as the flagship display (and the
  qgraph-in-Suggests note); `bgmCompare()` gets a short section with a
  pointer to its own vignette; the easybgm package named for extended
  plotting/summary workflows.
- A worked example, seconds-scale, seeded: small ordinal dataset (a shipped
  dataset if suitable), fit → summary → verdicts → plot. If any chunk needs
  more than a few seconds, bake it (`.rds` + generator script), matching the
  prior-sensitivity vignette's pattern.
- Bibliography header parity with the other vignettes (it is the only one
  without; use `inst/REFERENCES.bib` entries where claims need them).

## Build discipline (F-045 — measured, stated in the report)

- `R CMD build` knits all vignettes live: keep intro's knit under ~30 s.
  State the measured knit time in the report.
- Preserve the show-uncapped/run-capped chunk pattern used by the other
  vignettes for any fit the reader sees.

## Verification gate

1. `R CMD check --as-cran` on a `git archive` tarball: baseline NOTEs only
   (2 plain; a third examples-timing NOTE appears under `--run-donttest` —
   known, not yours).
2. Full default-tier suite: 0 failures / 0 warnings (nothing should move —
   say so).
3. The rendered HTML reads clean end to end; every code chunk runs; every
   claim about the package checked against the current Rd/NEWS (list any
   place the old vignette said something now false — those are the
   findings).
4. Knit time measured and under budget.

## Deliverable

`/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms/dev/review-2026-08/reports/26-intro-vignette.md`:
What was done (commits) / the old-vs-new content map (what was dropped, what
was added, and why) / knit time / proposed NEWS line VERBATIM (the vignette
existed at 0.1.6.3, so its rewrite is tag-visible) / findings / open
questions. **Copy the rendered `intro.html` to
`dev/review-2026-08/reports/assets/f026-intro.html`** — the maintainer
judges the rendered vignette before merge. Copy the report to the Dropbox
path; commit report + asset on the branch.

# Brief 02 — Hands-on evaluation of the new user-facing checks (MM)

Time box: ~2–3 hours at the console. You are evaluating as a *user* — a
methodologically sophisticated one — not reading code. Record findings as you
go, each with a severity (blocker / major / minor / note).

## Orientation — where this layer sits

0.2.0.0 adds a user-facing "checking and reporting" layer that post-processes
fit objects; none of it runs during sampling. The components and their call
paths:

- **`verdicts(fit)`** (`R/verdicts.R`) — three-category edge classification
  (presence / absence / undecided) from inclusion BFs; consumes the extractors
  (`extract_inclusion_bf()`, RB-based posterior inclusion probabilities from
  PR #182). Print method `print.bgms_verdicts`.
- **`extract_centrality(fit)`** (`R/centrality.R`) — posterior centrality
  measures with `summary()` and `plot()` methods; consumes posterior draws of
  the pairwise parameters.
- **`calibration_check(fit)`** (`R/calibration_check.R`) — model-calibration
  check with print/plot methods; consumes the simulate/predict machinery. Its
  scope was deliberately narrowed for continuous data — decision record in
  `dev/audit/2026-07-30-ppc-continuous-scope-decision.md`.
- **Network plots** (`R/plot_bgms.R`) — `plot.bgms`, `plot.bgmCompare`,
  `plot_edge_posterior`; the evidence-split display convention comes from the
  docs-site A6 work.
- **`prior_sensitivity_check(fit)`** (`R/prior_sensitivity.R`, PR #183) — the
  anchored sensitivity curve: `refit_engine.R` refits the model at scaled slab
  values (the 1× anchor is your own fit), `anchor_curve.R` assembles the
  continuous inclusion-BF curve by importance reweighting between anchors and
  inverse-variance pooling of anchors above an ESS floor. Print/plot methods.
  Design record: `dev/audit/2026-07-28-anchored-curve-spec.md`.

API design rationale for the whole layer:
`dev/audit/2026-07-30-user-facing-checks-api-proposal.md`.

## Setup

Use a **clean install of the frozen target** — not the Dropbox working tree
(stale `.o` files make it produce unloadable builds). Easiest: after the Opus
agent finishes brief 01, `.libPaths(c("~/bgms-review/lib", .libPaths()))` gives
you its rc1 install. Otherwise replicate its recipe (git archive
`v0.2.0.0-rc1` to a clean dir, `R CMD build` + `R CMD INSTALL`).

Data: `Wenchuan` (shipped with the package), your usual fit settings. Fit once
at defaults, keep the fit object for all components. Note fit runtime.

## Questions to answer (per component, in order)

1. **Defaults**: run each function with no arguments beyond the fit. Is the
   default output something you would put in a paper or show a PhD student
   without caveats? If not, what exactly is off?
2. **Wording vs computation**: does every printed sentence match what is
   actually computed? Check each print method against its man page
   (`?verdicts`, `?calibration_check`, `?prior_sensitivity_check`, ...). We
   already know NEWS.md misdescribes the sensitivity pooling (finding F-001) —
   look for the same class of drift in the print methods themselves.
3. **Terminology**: do category names, evidence thresholds, and prior names
   match the tutorial/guidelines vocabulary (presence / absence / undecided;
   the 3/10/30 discussion)? List every mismatch — this is the
   decide-once-with-NS terminology set.
4. **Sensitivity check specifically**: at defaults on the Wenchuan fit —
   runtime; does the stability headline match what the curve shows; does the
   1×-anchor verdict set equal your reported analysis exactly (the design
   property that makes the check defensible); is the plot the one you would
   ship in the tutorial?
5. **Centrality**: which measures, on what scale, are the uncertainty
   intervals interpretable? Anything statistically misleading?
6. **Calibration**: what does it actually check, is that what the name
   promises, and does the output tell a user what to *do*?
7. **Plots**: publication-ready at defaults (labels, legends, sizing)? Is the
   evidence-split network consistent with the tutorial's three-panel figures?
8. **Failure behavior**: run `verdicts()` and `prior_sensitivity_check()` on a
   deliberately bad fit (e.g. 100 iterations, no warmup) — are the
   errors/warnings helpful or misleading?
9. **The embarrassment test**: anything here that a sharp reviewer of the
   tutorial paper could use against the package?

## Deliverable

Write your report to `dev/review-2026-08/reports/02-user-facing-checks.md`:

1. **What was done** — fit settings, versions, which functions exercised.
2. **Findings** — numbered, severity-tagged (blocker / major / minor / note).
3. **Evidence** — console excerpts / saved plots (drop images in
   `dev/review-2026-08/reports/assets/` if useful).
4. **Open questions** — anything needing a code-level answer; these become
   your next code-read briefs.

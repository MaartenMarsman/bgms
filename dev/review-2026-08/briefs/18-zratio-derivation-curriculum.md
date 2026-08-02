# Brief 18 — assemble the zratio derivation curriculum (Opus agent; F-098)

You are the execution agent in the pre-CRAN review of the R package **bgms**.
You share no context with the review lead; everything you need is here.
**READ AND WRITE ONLY** — you produce ONE new document plus your report; you
change nothing in the package, and you treat every source directory below as
READ-ONLY (several are live manuscript projects).

## Why this brief exists

The maintainer must be able to verify the hierarchical-prior z-ratio
machinery's mathematics before release. He can follow derivations when they
are in front of him — but the derivations are scattered across two manuscript
projects, several side projects, and a stack of dev notes, and were never
assembled into one readable chain. Your deliverable is that chain: a
self-contained tutorial-grade document of the DEPLOYED z-ratio path, every
step derived, sourced, or explicitly flagged as ungrounded.

## Machine budget

Effectively none: reading and writing. Tiny numeric illustrations in R are
fine (seconds, single process). No fits, no builds.

## Sources (read-only)

**The mathematics (the ground truth):**
- `~/Dropbox/Projecten/sv/ggm_paper/` — the GGM paper; its
  `analytic-correction-companion.tex/.pdf` IS "the companion" that
  `src/models/ggm/zratio_law.h` and `zratio_engine.h` comments reference;
  its `R/` holds the companion implementation the C++ was ported from.
- `~/Dropbox/Projecten/sv/Z/` — the normalizing-constant project
  (`manuscript.tex`, `estimate_ggm_w_tilt.pdf`, `notes/`).
- Related SV projects as needed: `Z_follow_up/`, `z_graph_prior/`,
  `spikeslab/` (`spike_slab_w_tilt.pdf` — the tilted spike-and-slab prior),
  `zising/`, `degord/`.

**The engineering record (bgms repo, Dropbox
`~/Dropbox/Projecten/R/bgms/`):**
- `dev/plans/backlog/2026-07-02_hier-zratio-approximation_PLAN.md` (the
  original approximation plan), `dev/plans/active/
  2026-07-13_zratio-calibration-redesign_DESIGN.md`,
  `2026-07-22_zratio-surface-b-migration_PLAN.md` (Option B),
  `2026-07-14_zratio-calibration-and-trust-gauge_NOTE.md`,
  `2026-07-06_zratio-standardized-cell-cauchy_ANALYSIS.md`;
  `dev/audit/2026-07-23-pr172-zratio-surface-review.md`;
  `dev/validation/zratio_gold_bank.md`.
- The code: `src/models/ggm/zratio_engine.{h,cpp}`, `zratio_gauge.h`,
  `src/zratio_interface.cpp` (the block-Gibbs oracle
  `zratio_block_oracle_moments` — the DEPLOYED anchor source),
  `R/zratio_surfaces.R`, `R/zratio_tables.R`, `R/zratio_gauge.R`; tests
  `tests/testthat/test-zratio-*.R`.

## Scope: the DEPLOYED chain, in this order

1. **The model and the object.** The hierarchical specification
   p(Γ) p(K | Γ) with p(K | Γ) = ρ_Γ(K) / Z(Γ); why every between-graph
   move needs J = Z(Γ⁻)/Z(Γ⁺); the standardized cell and its invariance
   (the engine header's conventions block, `zratio_engine.h:128-131`).
2. **Locality.** Why J depends only on the toggled edge's mediating block
   (common neighbours + 2-hop bridges between exclusive neighbour sets) —
   the theorem and its proof or proof sketch, from the Z/ggm_paper sources.
3. **The two-moment saddle.** How J is evaluated from per-channel moment
   constants (S1, S2) of the tilted prior's pair integrals via the saddle;
   the additive composition over components; where the constants come from
   (`R/zratio_tables.R`, build at fit time) and the certified shape band
   [0.5, 20].
4. **The Option-B absolute-moment surface.** The alpha = 1 correction:
   log S as the 9-monomial bivariate quadratic in (log size, density), fit
   to block-Gibbs oracle anchors; hull clamps, the +/- 0.1 log-moment
   clamp, boundary-slope extension past size_hi, the additive fallback
   below size_min (exact there); why the MC oracle beat the analytic law
   as anchor source (cost + parity on gold — the `zratio_law.h` header
   states this).
5. **Routing.** Which requests are served by surface / additive saddle /
   isolated-edge route (shape > 10); the exact routing predicate and how it
   compares to the certification predicate.
6. **The trust gauge.** What the two assessment sweeps replay, what
   flip_rate and harm_pred measure, what the audit can and cannot resolve
   (coherent error vs rare edge-specific failures), the fixed per-chain cap
   as a sampling design.
7. **Appendix, one page max:** the dormant analytic law (`zratio_law.h`) —
   what it is, that it is validated (test-zratio-law.R) but NOT deployed
   and not in the paper, and when it would become load-bearing (the
   large-q crossover its header describes).

## Rules — these are the brief

- Every mathematical step is (a) DERIVED in place at working-statistician
  level (algebra shown), (b) QUOTED/adapted from a named source with file +
  section/equation citation, or (c) **FLAGGED as ungrounded** in a gap list
  that is a first-class deliverable. DO NOT invent derivations — a
  plausible-but-wrong step is this brief's failure mode and poisons the
  maintainer's read that builds on it. When a source derives something only
  for a special case, SAY SO.
- Every mathematical object gets a correspondence row: math symbol ↔ source
  notation ↔ code name ↔ file:line. (Watch the renamings the code comments
  hint at: the engine's `eta` is the companion's `beta`; `t2 = 2·beta·sigma²`.)
- Every approximation carries its error story AND where the review's
  existing evidence checks it (report 04's five certificates, the gold bank,
  the gauge's runtime audit) — the reader should always know what is proved,
  what is measured, and what is merely plausible.
- Audience: the maintainer — a Bayesian statistician expert in MRFs and
  graphical models who does not read C++ fluently. Algebra in full;
  numerics summarized; no code walkthroughs beyond the correspondence
  tables.
- Length: whatever the chain needs; expected order 15-30 pages of markdown
  with LaTeX math.

## Deliverable

1. `dev/review-2026-08/orientation/zratio-derivation-curriculum.md` — the
   document.
2. `dev/review-2026-08/reports/18-zratio-curriculum.md` — What was
   assembled from where (per section: derived / sourced-from / flagged) /
   **The gap list** (every ungrounded step, prominently) / Findings
   (anything you found where source, code, and dev-notes disagree — cite
   all sides) / Open questions.

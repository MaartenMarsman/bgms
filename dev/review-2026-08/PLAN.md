# Phased review plan — 0.2.0.0 (working document)

Review lead maintains this; updated as reports land. Brief ledger at the end.

## Target and scale

`cran-0.1.6.3` → `v0.2.0.0-rc1`: 233 commits, 429 files, ~+80k/−20k lines.
Package code: R +19.5k/−4.6k (51 files, 14 new); C++ reorganized (old
`src/bgm/`, `src/bgmCompare/` deleted; new `src/models/` tree ~16k lines,
`src/mcmc/` +2.8k, `src/priors/` +1.5k); tests +25.8k lines (123 files);
23 new exports, 0 removed; 6 vignettes (2 new).

## Strategy

1. **Verify known findings first.** `dev/audit/` (July 2026 audit + PR reviews
   + specs) and `plans/flagged-issues.md` are the starting inventory (report
   00a); rediscovery is waste.
2. **Risk-based depth.** Depth per subsystem ∝ novelty × statistical criticality
   × (defects already found there). The hierarchical/zratio path had late
   significant defects → deepest treatment.
3. **Hunt the resident defect class.** The bgmCompare scale bug (`omega*x` vs
   `2*omega*x` — one path diverging from the package-wide convention) is the
   signature failure mode: cross-path convention divergence (scales,
   parameterizations, prior framings). Every validation brief includes
   cross-path consistency checks, not just single-path correctness.
4. **Three lanes in parallel.** Opus agent executes (checks, suites,
   simulations); MM evaluates user-facing behavior and reads math-critical
   code (curriculum-sequenced); review lead triages, cross-checks, synthesizes.

## Phase 0 — Freeze and baseline (2026-08-01 → ~08-03)

- [x] Merge develop→main, tag `v0.2.0.0-rc1`; tag CRAN anchor `cran-0.1.6.3`.
- [x] Verify all six unmerged branches are absorbed or deliberately
      post-release (F-014); freeze integrity confirmed.
- [x] 00a prior-audit inventory; 00b C++ diff map; 00c R/API diff map
      (internal sweeps, done 2026-08-01 — see reports/, findings folded into
      FINDINGS.md F-017..F-034; MAINTAINERS.md seeded).
- [ ] **Brief 01 (Opus)**: baseline — clean-export build, CRAN-tarball anchor
      verification (F-012), `R CMD check --as-cran`, full suite with and
      without CRAN skips.
- [ ] **Brief 02 (MM)**: hands-on user evaluation of the new checking layer
      (verdicts, centrality, calibration, plots, sensitivity check).
- Exit: check/test status known; anchor verified; 00a findings folded into
  FINDINGS.md; Phase 1 briefs finalized against the risk map.

## Phase 1 — Correctness of the new statistical machinery

Provisional subsystem ranking (refined when 00a–00c land):

| # | Subsystem | Why this depth | Planned work |
|---|---|---|---|
| 1 | Hierarchical precision-graph prior / zratio path (`src/models/ggm/zratio_*`, `R/zratio_*.R`, `R/sample_graph_prior.R`) | Newest, most intricate (surface approximation of normalizing-constant ratios, anchor hulls, trust gauge, isolated-edge routing); late defects found here; ~2.5k C++ + ~2.3k R lines | Opus: rerun `dev/validation/` gold bank + route certificates against rc1; targeted SBC. MM: math read of `zratio_law.h` + gauge. Lead: cross-check vs PR #172/#193/#194 review docs |
| 2 | bgmCompare path (`src/models/bgmCompare/`, sampler rewritten, +749 lines) | Site of the scale bug; breaking semantic change ships this release | Opus: bgm-vs-bgmCompare cross-path consistency on identical single-group data; recovery study; verify the ~7s convention guard tier |
| 3 | GGM + mixed samplers (`src/models/ggm`, `src/models/mixed`, `src/mcmc/`) | New model classes (headline feature); NUTS/Gibbs restructuring; HMC removed, hamiltonian_utils new | Opus: SBC + parameter recovery; comparison vs BGGM (GGM) and mgm (mixed) where estimands overlap |
| 4 | RB inclusion machinery (PR #182; extractors, `src/mcmc_diagnostics.cpp`) | Changes every reported PIP/BF number; semantic break in `extract_ess()` | Opus: RB vs raw-indicator agreement on long reference runs. MM: estimator math read |
| 5 | Sensitivity check internals (`R/refit_engine.R`, `R/anchor_curve.R`) | Flagship adoption asset; pooling/reweighting math; NEWS already misdescribes it (F-001) | MM: pooling math read (after 02). Opus: small numerical check — curve vs brute-force refits at off-anchor scales |
| 6 | Ordinal MRF (`src/models/omrf/`, ~435 lines changed) + priors (`src/priors/`, SBM +552, edge-prior correction new) | Existing CRAN functionality touched; regression risk for current users | Opus: 0.1.6.3-vs-rc1 posterior agreement on fixed seeds/data (the true regression test); SBM correction tests review |
| 7 | Diagnostics (rhat df adjustment, ESS variants) | Fix docs exist (07-26/07-27); verify implemented-as-documented | Lead + targeted read |

Each row produces: an Opus validation brief, findings in FINDINGS.md, an
architecture-map section in MAINTAINERS.md, and (rows 1, 4, 5 minimum) an MM
code-read brief.

## Phase 2 — Surface, docs, and adequacy

- API coherence review of all 23 new exports (lead, from 00c).
- NEWS.md reconstruction 0.1.6.3→0.2.0.0 driven from the commit/PR log (MM
  brief; folds in F-001).
- Vignette accuracy sweep (F-006 mixture-ESS gloss + full pass), man-page
  spot-checks against behavior.
- Defaults-freeze memo to easybgm/JASP/docs/tutorial (F-002).
- Test-adequacy verdict from brief 01's skip analysis: what does CRAN actually
  exercise; nightly-validation workflow coverage.

## Phase 3 — CRAN mechanics and ship

- Resolve all blockers; land fixes on develop; re-merge develop→main;
  re-verify (full `--as-cran` + targeted re-tests of every fixed finding).
- Submission-batch commit: Date bump (F-003), README install line (F-004),
  sensitivity-check roxygen framing sentence (F-005).
- win-builder + mac-builder; reverse-dependency check (easybgm at minimum);
  `cran-comments.md`; bump-and-tag; submit.
- Post-acceptance: the docs-site launch chain (bgms-docs steps in
  `plans/flagged-issues.md`) — outside this repo, tracked there.

## MM curriculum (code reads, sequenced)

Goal: by the end, MM has read every high-risk statistical component with
orientation, cumulatively covering the architecture.

1. **02** New checking layer, hands-on (no code) — the API surface. *(issued)*
2. zratio law + trust gauge: the surface approximation's math and its guard
   rails (after 00b lands, with the Opus gold-bank rerun as companion).
3. RB inclusion estimator + mixture-ESS (with `src/mcmc_diagnostics.cpp`).
4. Sensitivity pooling: `assemble_curve()` + `refit_engine` importance
   reweighting.
5. Spec→sampler flow: one `bgm()` call traced end-to-end (R spec build →
   `.Call` → chain runner → model), closing the architecture loop.

## Brief ledger

| Brief | Assignee | Status | Report |
|---|---|---|---|
| 01 baseline check | Opus | **issued 2026-08-01** | pending |
| 02 user-facing checks | MM | **issued 2026-08-01** | pending |
| 03 hierarchical/zratio validation | Opus | draft after 00a/00b | — |
| 04 bgmCompare cross-path consistency | Opus | draft after 00b | — |
| 05+ | per Phase 1 table | — | — |

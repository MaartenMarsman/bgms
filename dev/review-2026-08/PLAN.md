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
- [x] **Brief 01 (Opus)**: DONE 2026-08-01 (report 01). Check clean (2 known
      NOTEs), suite 0 failures in every tier, determinism verified, anchor
      faithful (F-012 closed), CRAN test time a non-issue (F-030 closed).
      Corrections: examples NOT core-capped → F-018 upgraded to blocker;
      F-017 regraded major.
- [x] **Brief 02 (MM)**: DONE 2026-08-01 (report 02). New defects F-035..F-041,
      F-046/F-047 (calibration crash on BC/mixed; silent prior chain in
      verdicts on hierarchical fits; mixed-BF anomaly; title/print UX).
      Lead cross-checked: the reported "classification boundary bug" is a
      log10-vs-ln misread — classification verified correct (F-039).
- Exit: **met 2026-08-01.** Phase 1 opens with briefs 03 (defect batch) and
  04 (statistical certification).

## Phase 1 — Correctness of the new statistical machinery

Provisional subsystem ranking (refined when 00a–00c land):

| # | Subsystem | Why this depth | Planned work |
|---|---|---|---|
| 1 | Hierarchical precision-graph prior / zratio path (`src/models/ggm/zratio_*`, `R/zratio_*.R`, `R/sample_graph_prior.R`) | Newest, most intricate (surface approximation of normalizing-constant ratios, anchor hulls, trust gauge, isolated-edge routing); late defects found here; ~2.5k C++ + ~2.3k R lines | Opus: rerun `dev/validation/` gold bank + route certificates against rc1; targeted SBC. MM: math read of `zratio_law.h` + gauge. Lead: cross-check vs PR #172/#193/#194 review docs |
| 2 | bgmCompare path (`src/models/bgmCompare/`, sampler rewritten, +749 lines) | Site of the scale bug; breaking semantic change ships this release | Opus: bgm-vs-bgmCompare cross-path consistency on identical single-group data; recovery study; verify the ~7s convention guard tier. **Report 09 (2026-08-01): convention question CLOSED — no ×2/×½ anywhere (contrast slope 1.014, predict to 3e-16 with real power). Open residue: F-074 level offset, F-075 magnitude anomaly (control running), F-073 guard holes → brief 13** |
| 3 | GGM + mixed samplers (`src/models/ggm`, `src/models/mixed`, `src/mcmc/`) | New model classes (headline feature); NUTS/Gibbs restructuring; HMC removed, hamiltonian_utils new | Opus: SBC + parameter recovery; comparison vs BGGM (GGM) and mgm (mixed) where estimands overlap |
| 4 | RB inclusion machinery (PR #182; extractors, `src/mcmc_diagnostics.cpp`) | Changes every reported PIP/BF number; semantic break in `extract_ess()` | Opus: RB vs raw-indicator agreement on long reference runs. MM: estimator math read |
| 5 | Sensitivity check internals (`R/refit_engine.R`, `R/anchor_curve.R`) | Flagship adoption asset; pooling/reweighting math; NEWS already misdescribes it (F-001) | MM: pooling math read (after 02). Opus: small numerical check — curve vs brute-force refits at off-anchor scales |
| 6 | Ordinal MRF (`src/models/omrf/`, ~435 lines changed) + priors (`src/priors/`, SBM +552, edge-prior correction new) | Existing CRAN functionality touched; regression risk for current users | Opus: 0.1.6.3-vs-rc1 posterior agreement on fixed seeds/data (the true regression test); SBM correction tests review |
| 7 | Diagnostics (rhat df adjustment, ESS variants) | Fix docs exist (07-26/07-27); verify implemented-as-documented | DONE 2026-08-01, lead read: classic split-Rhat with the coda df-adjustment intentionally omitted (`src/mcmc_diagnostics.cpp:152-161` states the rationale); matches the fix doc; backlog item 39's branch SHA is post-squash but the content is in develop |

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
- The re-merge ACTIVATES the new CI schedules (crons fire from main's files):
  the first real Sunday run is the T2 budget proof (~70 min measured blocks +
  the 11-fit F-049 block vs `timeout-minutes: 360` — ample on paper, unproven;
  watch it) and un-reds the weekly-compliance harness once F-100's fix lands.
- Submission-batch commit: Date bump (F-003), sensitivity-check roxygen
  framing sentence (F-005). (F-004's README line landed in brief 12.)
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

## Brief conventions (standing; every future brief carries these)

- **Machine budget (MM directive, 2026-08-01).** The execution machine has 15
  cores and is SHARED — MM's interactive work plus up to three agents at once
  have produced 20 concurrent processes and full contention. Every brief's
  setup section must state: cap your total footprint at ~6 hardware threads;
  one fit at a time (`cores = 4` max); no concurrent sweeps, background R
  sessions, or parallel builds beyond that; SEQUENCE heavy steps. Runtime is
  not a grading criterion — contended runs waste more wall-clock than
  parallelism saves. Compute beyond what a check's power requires is waste:
  prefer design over brute force, and flag (don't run) any step estimated
  over ~30 min of busy cores.
- Reports to `dev/review-2026-08/reports/NN-<topic>.md` (What was done /
  Findings severity-tagged / Evidence / Open questions); seeds, runtimes, and
  derived tolerances stated. No pushing (unless the brief grants a scoped
  exception), no attribution trailers, never build in the Dropbox tree.
- **Steering as relay blocks (MM directive, 2026-08-02).** MM relays
  verbatim: ALL agent-directed content — briefs, mid-flight steering,
  follow-up questions, scope changes — is delivered to MM as a
  self-contained fenced block addressed to the agent, ready to paste.
  Never prose addressed to MM that he must excerpt.
- **Base worktrees on `origin/develop`, named explicitly.** The Dropbox
  checkout's LOCAL `develop` ref lags (two agents hit it, 23 commits stale
  at brief 11); briefs must say `origin/develop` in the worktree command
  and any base-verification check.
- **Ledger rows stay ONE line** (MM trimmed the overgrown rows,
  2026-08-02). Running state, verdicts, and evidence live in FINDINGS.md
  rows and reports; the ledger answers only who/what/status/report.

## Brief ledger

| Brief | Assignee | Status | Report |
|---|---|---|---|
| 01 baseline check | Opus | done | `reports/01-baseline-check.md` |
| 02 user-facing checks | MM | done | `reports/02-user-facing-checks.md` |
| 09–17 (validation, plots, tiering, NEWS, docs sweeps, F-075 program) | Opus/MM | done and merged — full state in FINDINGS rows + reports | `reports/09..17-*.md` |
| 11 plot restyle + evidence displays | Opus | round 3 in flight — state in F-067/F-079/F-115; merge gated on MM's delta look | branch `fix/plot-sweep` |
| 13 post-11 batch (CI reds, regression check, plot follow-ups, bake) | Opus | to draft after 11 merges — scope carried by FINDINGS "brief 13" routing tags | — |
| 03 checking-layer defect batch | Opus | done — **merged to develop `adf87013`** (MM confirmed the F-036 math, 2026-08-01) | `reports/03-checking-layer-defects.md` |
| 04 statistical certification | Opus | done — all 5 zratio certificates PASS vs the gold bank; 4 deterministic slow-tier failures, all test-side (F-048/F-049); drifter resolved (quick-fit artifact) | `reports/04-statistical-certification.md` |
| 05 test-repair + release-hygiene batch (F-048 stale fences, F-050 deprecation sweep, F-054 cache-key version, F-055 q≤3 guard, F-017 .Rbuildignore, F-051 docs line, `^\.git$` hardening) | Opus | done — **merged to develop `b04dbd06`** (all gates green; slow tier 4→1 failures, 1152→0 warnings; the remaining red is F-049 by design). New: F-058 (.Rbuildignore live regex), F-059 (additive band = uncertified band). Correction from the report: the scheduled nightly runs MAIN, so Monday is red regardless — see FINDINGS details | `reports/05-test-repair-batch.md` |
| 06 bgmCompare defect + units batch (F-061 sensitivity misalignment FIRST, F-056/F-057 input handling, F-060 verdicts print, F-021 method (adapt the parked patch), F-063 centrality removal, F-062a–d plots, nats-everywhere conversion, F-065 startup silence) | Opus | done — **merged to develop `4969f843`** (11 commits; full calibration method shipped, not the stub; nats invariance exact; F-049 layout hypothesis refuted + 20-seed diagnostic delivered — disposition with MM; new F-068 predict-BC fix in-batch, F-069/F-070 opened) | `reports/06-bgmcompare-defect-batch.md` |
| 07 MM: bgmCompare user pass (difference verdicts, group/difference plots, difference-scale sensitivity trace; ~1 h) | MM | done — report in 2026-08-01. Yields F-060..F-065, confirms F-057; decisions: natural log EVERYWHERE (closes the F-039 residue), defer F-021 with a 0.2.0 stub, drop difference centrality (F-063) | `reports/07-bgmcompare-user-pass.md` |
| 08 MM curriculum step 2: zratio machinery math read (carries F-059) | MM | **REISSUED as v2, 2026-08-02** — v1 stalled and the stall was diagnostic: v1's reading order led with the DORMANT `zratio_law.h` (lead error; its own header says not-deployed, not-in-paper) and the derivations were never assembled anywhere readable (F-098; companion comment leakage = F-097). v2: engine header `:111-131` readable TODAY; the full read runs on brief 18's curriculum after lead verification, INTERACTIVELY (MM pastes anything that loses him; confusion = documentation findings). F-059 rides unchanged | — |
| 18 zratio derivation curriculum assembly (F-098): ONE self-contained tutorial of the DEPLOYED chain (model → locality → two-moment saddle → Option-B surface + block-Gibbs oracle → routing → gauge; dormant law = 1-page appendix) from the REAL sources — `sv/ggm_paper` (the analytic-correction companion), `sv/Z`, related SV projects, dev/plans + dev/audit notes; derive-source-or-FLAG rules (invented math = failure mode); math↔source↔code correspondence tables; error stories tied to report-04 certificates | Opus | **issued 2026-08-02** — `briefs/18-zratio-derivation-curriculum.md`; read/write only, zero machine, launchable NOW in parallel; lead verifies the document before MM reads it. TIMING = MM-OWNED (2026-08-02): he launches when he has time; blocks only his own brief-08 read, nothing else (his ruling). If 08 lands post-release, F-059 resolves by its no-cost default: accept the additive-band coverage asymmetry as a documented property. Lead stops nudging. | — |
| 19 **F-075 fix batch (release gate)**: `collapse_categories_across_groups()` recode semantics per MM's decision (union vs interim adjacent-merge; BC-exemption question); F-110 always-warn + store the map in the fit; F-111 Rd + comparison-vignette correction; NEWS behaviour-change entry with refit guidance; regression tests — the 8-line reproducer as a unit test + an estimate-survival assertion vs matched bgm() fits; POST-FIX re-measure of report 16's operating-characteristics table (same 10 seeds) incl. the healthy-edge slope (resolves the 09-3 side of F-074); F-113 real-data support-difference scan (ADHD fixture first); F-112 compare logp test hook IF MM pulls it forward | Opus | **issued 2026-08-02** — `briefs/19-category-collapse-fix.md`. MM decisions IN: (1) UNION semantics; (2) BC exemption DELIBERATE (scores meaningful — Rd documents rationale); (3) F-112 hook pulled into the batch. Lead tag-verified the defect SHIPPED (0.1.6.3 `data_utils.R:388-394`, identical loop) ⇒ true Bug-fixes NEWS entry with refit guidance. Key risk called out: downstream-consumer audit before flipping semantics (zero-support cells in sufficient stats / starts / proposals). Re-measure reuses 16's truth + bgm baselines; 10 sequential compare fits. | handed to MM 2026-08-02 |
| 20 standalone zero-contention batch: F-100 compliance property sweep + fixture regen + green run (forward-compat assertions flagged for 19's merge); F-109 sparse-coding proof BOTH ends (current stored-map verified + tag miscode demonstrated at runtime) + regression test + NEWS clause proposed in report; F-114 extract_arguments main_effect_indices for compare fits (provide or document+helpful-error). NEWS/vignettes/compare-model/plot files OFF-LIMITS. | Opus | **issued 2026-08-02** — `briefs/20-standalone-batch.md`; parallel-safe with 11/18/19 by construction; ~4 threads | handed to MM 2026-08-02 |

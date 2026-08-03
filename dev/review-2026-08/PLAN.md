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
- Defaults-freeze memo to easybgm/JASP/docs/tutorial (F-002) — draft ready
  (report 25 §4.3 + change list §4.1); MM sends AT submission time (bgms
  first, notify at submission, not after — report 25 §4.4). Notification
  only, ZERO WAIT: bgms submits on its own schedule (MM 2026-08-03); the
  easybgm adaptation is Nikola's, in-house, with report §4.1 as his work
  order; easybgm 0.5.0 must FOLLOW bgms 0.2.0.0 to CRAN, never precede it.
  EASYBGM LIMB DISCHARGED 2026-08-03: standalone handoff
  `reports/25b-easybgm-handoff.md` written and RELAYED TO NIKOLA by MM.
  Residual audiences: JASP, docs site, tutorial (+ MM's call on a courtesy
  note to the easybgm CRAN maintainer of record at submission).
- Test-adequacy verdict from brief 01's skip analysis: what does CRAN actually
  exercise; nightly-validation workflow coverage.

## Phase 3 — CRAN mechanics and ship

- Unrouted-tail disposition sweep (lead): **DONE 2026-08-03** — every
  remaining open row dispositioned in FINDINGS (routed to the submission
  batch / Phase-3 / wrap-up, backlogged post-release, or accepted; F-006
  closed; F-080/081 open-watch for the nightly proof). MM ruled 2026-08-03 on
  three of the four: F-007 dropped (done), F-026 rewrite NOW (brief 26),
  F-046 status quo + graceful message + easybgm pointer (done). F-020
  RESOLVED ship-later (MM 2026-08-03, "Ship later indeed") — caveat landed
  `fbb0adaa` (`extract_rhat` @details), implementation = flagged-issues
  item 16. All four rulings executed; sweep CLOSED.

- Resolve all blockers; land fixes on develop; re-merge develop→main;
  re-verify (full `--as-cran` + targeted re-tests of every fixed finding).
- The re-merge ACTIVATES the new CI schedules (crons fire from main's files):
  the first real Sunday run is the T2 budget proof (~70 min measured blocks +
  the 11-fit F-049 block vs `timeout-minutes: 360` — ample on paper, unproven;
  watch it) and un-reds the weekly-compliance harness once F-100's fix lands.
  It is ALSO the live proof of the F-125 joint-pins (eight T2 blocks pinned at
  22's integration after the hierarchical flip re-pointed them) and of the two
  restored-at-0.02 gauge harm blocks (F-103).
- Check-gate baseline for future briefs/agents: 2 NOTEs plain `--as-cran`;
  **3 NOTEs under `--run-donttest`** (the third is the F-018 examples-timing
  NOTE, verified base-vs-branch identical in report 22 §7).
- Phase-3 entry gate (post-sign-off): run 1 came back 2-baseline + a
  SELF-INFLICTED third NOTE (the lead's wrapper wrote `build.log` into the
  build dir and `R CMD build` swept it into the tarball — harness lesson:
  build from a pristine archive dir, logs outside); rerun killed for machine
  contention (MM's own suite run measured 286 s against the agents' idle
  221.7 s — same job, contention only). Post-sign-off hygiene from MM's
  console observations landed first: F-129 (verbose leak) + F-130 (RATTLE
  retirement); final gate run includes both. **FINAL GATE GREEN 2026-08-03
  (tree `f5b823d1`): Status: 2 NOTEs, both baseline (stale Date / local HTML
  Tidy); testthat `[108s/86s]` — the idle-machine number, confirming the
  earlier 118 s was contention; donttest examples `[485s/248s]` OK, NO
  timing NOTE; vignettes `[93s/50s]` OK. Phase-3 entry condition met; next
  step = develop→main re-merge (MM's word).**
- Submission-batch commit: Date bump (F-003), sensitivity-check roxygen
  framing sentence (F-005). (F-004's README line landed in brief 12.)
- win-builder + mac-builder; reverse-dependency check (easybgm at minimum —
  dress rehearsal DONE, report 25 §3.4: `--as-cran` status IDENTICAL on both
  builds (1W/2N, none bgms-caused), runtime the only delta; re-run against
  the final tarball at rc); `cran-comments.md`; bump-and-tag; submit.
- F-121 residue (MM lifted 2026-08-03): one containerized ASan run (instrumented R,
  e.g. `wch1/r-debug` clang-ASAN) — the only defect class the diagnosis leaves
  unexamined; seed-docs portability sentence (both fits) lands at 22's
  integration (lead); afterwards clear `~/bgms-review/f121/` incl.
  `obj-fix12/` (only surviving copy of the anomalous objects).
- Compliance baseline (F-116): fixtures are machine-specific; structure-only
  stays green everywhere. Optional pre-release: regenerate the canonical set on
  `ubuntu-latest` (needs a regenerate job on `weekly-compliance.yaml`); required
  before any bitwise restoration, which also meets the dead `identical()` at
  `tests/compliance/test_compliance.R:606`.
- Wrap-up housekeeping (MM, 2026-08-03): clear `~/bgms-review/` — merged
  worktrees and local `fix/*` branches, exports, temp libs, `val16/` outputs —
  only after every batch that reads them has merged.
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
| 11 plot restyle + evidence displays | Opus | DONE — merged `f3f18941`, NEWS `66f7d17f`; record in F-067 | branch `fix/plot-sweep` |
| 13 pre-release sweep | Opus | DONE — merged `8bfa6148` (11 task commits, both tiers 0/0/0, --as-cran 2 NOTEs, vignette baked at 9.4 KB); F-103 escalated to MM, F-080/081 await CI proof, F-122 opened; records in FINDINGS | reports/13-pre-release-sweep.md |
| 03 checking-layer defect batch | Opus | done — **merged to develop `adf87013`** (MM confirmed the F-036 math, 2026-08-01) | `reports/03-checking-layer-defects.md` |
| 04 statistical certification | Opus | done — all 5 zratio certificates PASS vs the gold bank; 4 deterministic slow-tier failures, all test-side (F-048/F-049); drifter resolved (quick-fit artifact) | `reports/04-statistical-certification.md` |
| 05 test-repair + release-hygiene batch (F-048 stale fences, F-050 deprecation sweep, F-054 cache-key version, F-055 q≤3 guard, F-017 .Rbuildignore, F-051 docs line, `^\.git$` hardening) | Opus | done — **merged to develop `b04dbd06`** (all gates green; slow tier 4→1 failures, 1152→0 warnings; the remaining red is F-049 by design). New: F-058 (.Rbuildignore live regex), F-059 (additive band = uncertified band). Correction from the report: the scheduled nightly runs MAIN, so Monday is red regardless — see FINDINGS details | `reports/05-test-repair-batch.md` |
| 06 bgmCompare defect + units batch (F-061 sensitivity misalignment FIRST, F-056/F-057 input handling, F-060 verdicts print, F-021 method (adapt the parked patch), F-063 centrality removal, F-062a–d plots, nats-everywhere conversion, F-065 startup silence) | Opus | done — **merged to develop `4969f843`** (11 commits; full calibration method shipped, not the stub; nats invariance exact; F-049 layout hypothesis refuted + 20-seed diagnostic delivered — disposition with MM; new F-068 predict-BC fix in-batch, F-069/F-070 opened) | `reports/06-bgmcompare-defect-batch.md` |
| 07 MM: bgmCompare user pass (difference verdicts, group/difference plots, difference-scale sensitivity trace; ~1 h) | MM | done — report in 2026-08-01. Yields F-060..F-065, confirms F-057; decisions: natural log EVERYWHERE (closes the F-039 residue), defer F-021 with a 0.2.0 stub, drop difference centrality (F-063) | `reports/07-bgmcompare-user-pass.md` |
| 08 MM curriculum step 2: zratio machinery math read (carries F-059) | MM | **NON-BLOCKING for release (MM 2026-08-03: "i will do it, but lets not make it a blocking change")**; REISSUED as v2, 2026-08-02 — v1 stalled and the stall was diagnostic: v1's reading order led with the DORMANT `zratio_law.h` (lead error; its own header says not-deployed, not-in-paper) and the derivations were never assembled anywhere readable (F-098; companion comment leakage = F-097). v2: engine header `:111-131` readable TODAY; the full read runs on brief 18's curriculum after lead verification, INTERACTIVELY (MM pastes anything that loses him; confusion = documentation findings). F-059 rides unchanged | — |
| 18 zratio derivation curriculum assembly (F-098) — **NON-BLOCKING for release (MM 2026-08-03)**: ONE self-contained tutorial of the DEPLOYED chain (model → locality → two-moment saddle → Option-B surface + block-Gibbs oracle → routing → gauge; dormant law = 1-page appendix) from the REAL sources — `sv/ggm_paper` (the analytic-correction companion), `sv/Z`, related SV projects, dev/plans + dev/audit notes; derive-source-or-FLAG rules (invented math = failure mode); math↔source↔code correspondence tables; error stories tied to report-04 certificates | Opus | **issued 2026-08-02** — `briefs/18-zratio-derivation-curriculum.md`; read/write only, zero machine, launchable NOW in parallel; lead verifies the document before MM reads it. TIMING = MM-OWNED (2026-08-02): he launches when he has time; blocks only his own brief-08 read, nothing else (his ruling). If 08 lands post-release, F-059 resolves by its no-cost default: accept the additive-band coverage asymmetry as a documented property. Lead stops nudging. | — |
| 19 F-075 fix batch (release gate: F-075/110/111/112/113) | Opus | DONE — merged `02c609c3`; GATE CLEARED, OCs at bgm parity 10/10 seeds; F-074 settled; F-117/F-118 opened; record in F-075 | reports/19-category-collapse-fix.md |
| 20 standalone batch (F-100/F-109/F-114) | Opus | DONE — merged `a1c0e768`, NEWS `121890e5`; F-116 opened (fixture provenance); records in FINDINGS | reports/20-standalone-batch.md |
| 21 mark prior-only compare summary rows (F-117) + F-121 diagnosis | Opus | **DONE 2026-08-03** — merged `20e0e1d7` (widened mark, 906-char warning, BC pin, build-invariant snapshot; lead-verified green on a second build); report + addendum landed at integration; F-121 mechanism: legal FP compilation variance (-O0 stowaway build), NOT UB — lift pending MM | report 21 + addendum |
| 26 intro-vignette rewrite against 0.2.0.0 (F-026, MM-ordered) | Opus | **executed 2026-08-03, RENDER REJECTED — NOT MERGED** — branch `fix/intro-vignette` (`36559f84`+`e3aec725`) parked; gates green (2 baseline NOTEs, 8518/0/0/0, knit 16.6 s) but MM rejected the frame: all models are MRFs, GRAPHICAL not network models, no regularization contrast, no spike-and-slab-led pedagogy (F-026 row has the verbatim + the four rulings). MM ruled FIX 2026-08-03 ("Ok fix") — brief 26b executed and MERGED `7bc67695`; MM rewrote the intro opening himself and SIGNED OFF ("Done. And signed off.", landed `cd2013a3`); NEWS + report-26b micros landed `a75052aa`. **CLOSED — Phase 2 COMPLETE 2026-08-03** | report 26; report 26b; MM sign-off 2026-08-03 |
| 25b easybgm handoff for Nikola (MM-ordered follow-on to 25) | Opus | **DONE 2026-08-03, RELAYED TO NIKOLA by MM** — `reports/25b-easybgm-handoff.md` (799 lines, self-contained: S7-native change, MCSE_BF fix, moved numbers, runtime one-argument fix, new model classes, two pre-existing bugs, release plan with the `bgms (>= 0.2.0.0)` pin, reproduction); same round strengthened report 25 in place (N4 build-verified, E10/E12 archaeology, F-025-7 cross-backend table) — landed at lead integration, spot-checks green | landed post-relay |
| 25 easybgm compatibility report (shim exercise vs CRAN easybgm, delta enumeration, revdep dress rehearsal, submission-ordering answer) | Opus | **DONE 2026-08-03** — merged `23bf0234` + lead fixes `fbb0adaa` (B1/B2 + F-020 caveat); 0 BREAKS: easybgm 0.4.0 suite 0/0 on both builds (lead re-ran, 185/0/0 reproduced), 14 workflows 0 errors on new, `--as-cran` status IDENTICAL; ordering = bgms FIRST with notification; only revdep-visible change is runtime (check 83→196 s; two heavy examples ~6.4× elapsed → flagged item 17); memo draft §4.3 ready for MM at submission; findings registered as F-126 | report 25 |
| 22 compare slab Cauchy → Normal (F-119) + ridge visualization + zero-support difference-test evaluation | Opus | **DONE 2026-08-03** — merged `436a7a89`; both defaults flipped (baseline + difference families), shim guards re-pointed (F-119-a), tautological default test fixed (F-119-e); ridge 7.5 → 2.9 (MM judges figures); 0/16 sampling-zero over-calls; escalation declined on the paired contrast (lead accepted; full-ten completion running, F-119-b); NEWS three edits + seed sentence landed; T2 joint-pins (F-125) ride the same integration | report 22 |
| 23 mixed relational certification (vs OMRF sharp, vs GGM calibrated, cross-block planted recovery) | Opus | **DONE 2026-08-03** — merged `0e94e62e`; mixed NOT user-reachable on pure data (internal-machinery cert); sharp reduction blocked by F-123 (degenerate-block spin/error, first-class); substitutes certify mixed≡OMRF (slope 1.0036, ΔPIP ≤ 0.0096, residual = companion by measurement) and mixed-vs-GGM shrinking 1.00239→1.00037 (no convention bug); cross-block 15/15, 0/30, slope 1.039; mgm 9/9 cross; F-124 opened (ord-ord slope 1.164, Q1) | report 23 |
| 24 precision_graph_prior default → "hierarchical" (F-010) + F-022 harm wiring + F-103 sweeps program + F-123 (b) guard/hygiene | Opus | **DONE 2026-08-03** — merged `d3d8dbcf`; flip shipped with request-aware advisory; mixed harm channel live (pool-aware gain, GGM byte-identical); F-103 STOP branch: lever spent, dispersion is the fit's not the audit's, threshold question back to MM with decision inputs; F-123 guard + 7 bounds (8 residual bounds = lead follow-up); NEWS N1–N4 landed (parse 33); gate T0/T1 0F/0W + as-cran 2 baseline NOTEs | report 24 |

# Master findings list — 0.2.0.0 pre-CRAN review

Maintained by the review lead. One entry per deduplicated finding.
Severity: blocker / major / minor / note. Status: open / fixed-claimed /
fixed-verified / decided / wontfix. Every fixed finding needs a verification
(commit, test, or report reference) before it becomes fixed-verified.

Seeded 2026-08-01 from the cross-repo tracker (`plans/flagged-issues.md`),
MM's A2 review log (`plans/a2-review-findings.md`), and repo state checks.
Findings from the `dev/audit/` inventory and the mapping sweeps are folded in
as reports 00a–00c and later land.

| ID | Severity | Status | Finding | Source |
|---|---|---|---|---|
| F-001 | blocker | open | NEWS.md misdescribes the sensitivity check: says each curve point is "read from the anchor with the highest importance ESS"; `assemble_curve()` and the man page do inverse-variance pooling of every anchor clearing the ESS floor. NEWS ships with the release. | flagged-issues #1 |
| F-002 | major | open | Defaults-freeze memo to downstream (easybgm, JASP, docs site, tutorial) not yet sent: bgm() interaction prior Cauchy→Normal(1); precision_scale_prior raw-rate→exponential-eta; new precision_graph_prior="joint"; update_method gains "gibbs"; standardize= lifecycle-deprecated; extract_ess() indicators now RB n_eff (transition ESS via estimator="mixt"); extract_inclusion_bf() gained log=. First two change numbers downstream callers quote. | flagged-issues #2 |
| F-003 | blocker | open | DESCRIPTION `Date: 2026-03-26` is stale (guaranteed CRAN NOTE). Bump in the submission commit, not before. | flagged-issues #5 |
| F-004 | minor | open | README development-install line lacks `@develop` ref — plain `install_github` call installs main. Fix in the submission batch. | flagged-issues #5 (A2 pass) |
| F-005 | minor | open | `prior_sensitivity_check()` roxygen Details needs one framing sentence before the bare `\insertRef{bartos2026}`: Bartos et al.'s reweighting identity is used locally around refit anchors; the anchored construction is the package's own. MM approves wording; regenerate Rd. | flagged-issues #5 |
| F-006 | minor | open (downgraded) | The wrong mixture-ESS gloss is ALREADY FIXED in `vignettes/diagnostics.Rmd` source (rewritten by `b42df13f`, #190; verified in report 00c). Residue: (a) soft gloss at `vignettes/diagnostics.Rmd:88` — flip counts framed as "how much the chain explored the two states", the weak survivor of the rejected reading; (b) stale pre-build artifact `doc/` dir still carries the old text (`doc/diagnostics.Rmd:68-72`; `.Rbuildignore`d, but delete to prevent confusion). | a2 #1; report 00c §6 |
| F-007 | minor | open | `Readme.Rmd` stale vs hand-edited `README.md`; knitting would regress the README. Decide: drop the Rmd (recommended) or re-sync. | flagged-issues #12 |
| F-008 | note | open | `.github/copilot-instructions.md` gap after pkgdown removal: "register a new exported function" procedure should point at the bgms-docs reference workflow. | flagged-issues #13 |
| F-009 | note | decided | Build hygiene: never build/check in the Dropbox tree (stale `src/*.o` against old RcppParallel + generated `src/Makevars` → unloadable `bgms.so`; Dropbox sync resurrects files mid-render). Standing practice: `git archive` to a clean dir outside Dropbox for every build. | flagged-issues #16, #14; A5 §2 |
| F-010 | note | decided | `precision_graph_prior` default STAYS "joint" for 0.2.0.0 (hierarchical path: trust gauge off by default + post-hoc only, silent clamping beyond size-44/22 anchor hulls, minutes-order surface-build pause, continuous/mixed-only). Revisit post-release with the GGM paper; flip would be a breaking change. Decision record: flagged-issues #2b, PR #172 review. | flagged-issues #2b |
| F-011 | note | fixed-verified | pkgdown decommissioned before any release tag: `.github/workflows/pkgdown.yaml` and `_pkgdown.yml` absent from rc1 (verified in merge diff, 2026-08-01). Required because the workflow also fired on published releases. | flagged-issues #4/2; verified in rc1 |
| F-012 | major | open | In-repo CRAN anchor (`cran-0.1.6.3` = `18e660a2`) not yet verified against the actual CRAN 0.1.6.3 tarball. All review diffs assume they match. Verification is task 2 of brief 01. | review lead, 2026-08-01 |
| F-013 | note | open | `develop` is pushed and in sync (after the 2026-08-01 history repair), but `main` (= rc1 merge) and the tags `v0.2.0.0-rc1` / `cran-0.1.6.3` are still local-only. Push main + both tags when MM decides to publish the rc. | review lead, 2026-08-01 |
| F-014 | note | open | Branch hygiene: six unmerged local branches. Verified 2026-08-01 that all release-relevant content is in rc1: `fix/bgmcompare-association-scale` and `feature/indicator-diagnostics-cleanup` are patch-equivalent in develop (`git cherry` = 0); `feature/hier-followup`, `feature/hier-gating`, `feature/user-facing-checks` residual diffs are pre-squash/pre-refactor states (develop strictly newer — checked hunks). `feat/graphical-g-prior` (May 2026, +14k/−37k vs develop) is the deliberate post-release GGM line — keep. Delete the five absorbed branches after MM confirms; retire the stale `wt-ufc` worktree on user-facing-checks. | review lead, 2026-08-01 |
| F-015 | note | fixed-verified | The 2026-07-31 scale-fix landing decisions (`dev/audit/2026-07-31-scale-fix-landing-decisions.md`) were executed: (1) fix landed as its own commit, patch-equivalent in develop; (3) the scale-contingency sentence is in `vignettes/comparison.Rmd:70-75`; NEWS.md "Breaking changes" documents the omega*x vs 2*omega*x divergence and the ~½ rescaling of reported pairwise effects. Still to verify: (2) the ~7s cross-implementation guard runs in the every-run test tier (check in report 01's test listing). | review lead, 2026-08-01 |
| F-017 | blocker | open | `.Rbuildignore` gaps: `tests/compliance/` (~2.9 MB weekly bitwise-vs-CRAN harness), `tests/fixtures/`, and stale `tests/testthat/_problems/` all ship in the tarball. Add `^tests/compliance$`, `^tests/fixtures$`, `^tests/testthat/_problems$`; delete `_problems/`. Confirm empirically from brief 01's built tarball listing. | 00a AUD-R1; 00c §5,7 |
| F-018 | major | open | `cores = parallel::detectCores()` remains the signature default in `R/bgm.R:495`, `R/bgmCompare.R:207`, `R/simulate_predict.R:127` — a known CRAN-reviewer trigger and a surprise on user machines. Examples now capped at 2, but verify tests/vignettes stay ≤2 workers under CRAN settings and decide: keep with policy comment, or default to a capped value. | 00a AUD-R2 |
| F-019 | major | open | AUD-H6 (July audit) is the ONLY Phase-1 correctness item with no visible landing: no n_==0 positive-definiteness guard inside `update_edge_indicator_conjugate` (`src/models/ggm/ggm_model.cpp:995`), though `ggm_edge_move`/`ggm_diag_move` have guards (`:288`, `:682`). Settle by targeted read + a long `sample_ggm_prior(update_method="gibbs")` run. | 00a Top-10 #4 |
| F-020 | major | open | SLAB-2: no cross-chain conditional-means check on slab weights in `R/mcmc_summary.R` — the one detector with measured full power for weight multimodality where indicator diagnostics are structurally blind (slab-diagnostics brief). Decide: implement for 0.2.0.0, or ship-later with a NEWS/docs caveat. MM decision. | 00a SLAB-2 |
| F-021 | minor | open | `calibration_check.bgmCompare` exists as a ready 19 KB patch (`dev/audit/2026-07-30-b4-calibration-bgmcompare.patch`) but was never applied; NAMESPACE has only the `bgms` method. Decide apply-or-drop for 0.2.0.0 (also closes part of the extractor asymmetry in 00c §1). | 00a PAR-B4 |
| F-022 | minor | open | Mixed hierarchical fits: `harm_pred` is permanently NA — `R/build_output_mixed_mrf.R:305` calls `summarize_zratio_gauge()` without `harm_inputs`. Wire the inputs or state the per-model limitation in the Rd/vignette. | 00a PR172-open |
| F-023 | minor | open | NEWS lacks a `simulate_mrf()`-specific sentence for the factor-2 input-scale change (a user-supplied `pairwise` matrix means something different now). One sentence under Breaking changes; fold into the NEWS reconstruction. | 00a AUD-R3 |
| F-024 | note | fixed-verified | The uncommitted 10-line edit to `tests/testthat/test-zratio-surface-build.R` was committed on develop as "test(zratio): stop paying for constants and anchors the surface-build claims do not use" (2026-08-01, now `dd43696a` after the history repair). Note: this test change postdates rc1 — the rc1 tag intentionally does not contain it. | 00a; resolved 2026-08-01 |
| F-025 | minor | open | Verify shipped wording on the two decided-open hierarchical limitations: (a) sub-shape-0.5 additive zero-collapse — NEWS/`?bgm` text vs the measured 12–32-variable boundary (`dev/plans/active/2026-08-01_hier-followup_NOTE.md` Disposition); (b) the eta<2 gate deviation record lives only in ggm_paper's `deployment-approach.md` — confirm it was updated or copy the decision in-repo. | 00a HF-1/PR172-B2 |
| F-026 | minor | open | `vignettes/intro.Rmd` untouched since 2025-07-01: predates the 0.2 API (no prior constructors, no verdicts, no `precision_graph_prior`), only vignette without a bibliography header. The front-door vignette misrepresents the package. | 00c §6 |
| F-027 | major | open | DESCRIPTION has no `SystemRequirements:` while the package ships `configure` + `src/Makevars.in` using `include sources.mk` (GNU-make construct) and links RcppParallel/TBB. Add `SystemRequirements: GNU make` (confirm against brief 01 check output). | 00c §7 |
| F-028 | minor | open | `RoxygenNote:` removed from DESCRIPTION, replaced by `Config/roxygen2/version: 8.0.0` (ahead of CRAN's roxygen2 line). Unconventional; confirm the Rd regenerate cleanly with a released roxygen2 and restore the conventional field. | 00c §7 |
| F-029 | major | open | Silent behavior changes for existing users whose NEWS coverage must be verified line-by-line during NEWS reconstruction: `iter`/`warmup` defaults 1e3→2e3; NUTS `target_accept` 0.60/0.65→0.80; `hamiltonian-mc` + `hmc_num_leapfrogs` removed; new formals `means_prior`, `threshold_prior`, `delta`, `progress_callback`, `difference_family`; `extract_ess` default now RB with no warning on new fits; prediction/simulation on the ×2 association-scale convention; bgmCompare RB PIPs return NA for unselected indicators (breaks all-numeric assumption downstream — add to the F-002 memo). | 00c §3; 00b |
| F-030 | major | open | CRAN check-time risk: ~85% of the 1,039 test blocks run under CRAN settings, and `tests/testthat/helper-fixtures.R` fits real `bgm()` models at first access even in nominally cheap files. Quantify total CRAN-mode test time in brief 01; if over ~10 min the suite needs a CRAN-tier gate. | 00c §5 |
| F-031 | minor | open | Rd polish before submission: `mrfSampler` is a live export with `\keyword{internal}` and no `\value` (CRAN checks exported Rd for \value); 25 Rd files lack `\examples{}` incl. 12 exported extractors; `withr` sits in Suggests — confirm used or drop. | 00c §2,7 |
| F-032 | note | open | C++ debris (backlog, not release): dead `// [[Rcpp::export]]` at `src/math/cholupdate.cpp:123` (never scanned; function dead weight); 45-line commented-out harness at `src/utils/progress_manager.cpp:406-450`; `src/models/ggm/zratio_law.h` (501 LOC) dormant-by-design but compiled into the shipped binary (documented insurance — keep, but note). | 00b §17 |
| F-033 | note | open | Engineering-risk flags for Phase-1 expert reads (no confirmed bug): adaptation controllers hold reference members (`src/mcmc/samplers/metropolis_adaptation.h:19` `arma::mat&`; `NUTSAdaptationController` holds `&schedule`) — lifetime discipline to verify; `base_model.h` has three parameter vectorizations with defaulted fallbacks (silent mis-storage risk); `thread_local` scratch inside the bgmCompare TBB worker (`bgmCompare_logp_and_grad.cpp:604-605`). | 00b §5,6,11 |
| F-034 | note | open | bgmCompare is the only sampler path NOT on `mcmc/execution/chain_runner` and the only one with NO test interface — its new `CompareSweepState` caching is reachable only end-to-end. Elevates Phase-1 row 2 (cross-path consistency validation is the primary defense); the BaseModel port is backlog AUD-D2. | 00b §1,15; 00a AUD-D2 |

## Details and verification notes

**F-001** — one-line fix in NEWS.md; verify by rereading `man/prior_sensitivity_check.Rd`
and `R/anchor_curve.R::assemble_curve()` after the edit. Fold into the NEWS
reconstruction task (PI) rather than fixing twice.

**F-002** — not a code change. Producing the memo is a review deliverable; the
five bgm() deltas were code-confirmed by BGMS-1 (see its defaults table) and are
re-verified in report 00c. Send after the defaults table is final.

**F-006** — the same wrong gloss existed in the docs site (routed to brief A6
there); the vignette is the in-package copy. Corrected reading is drafted in
`plans/a2-review-findings.md` and grounded in `src/mcmc_diagnostics.cpp:375-387`.

**F-012** — if the tarball diff shows real differences, every "unchanged since
CRAN" claim in the diff map must be rechecked against the true baseline.

**F-016 (major, open, maintainability+process)** — `dev/audit/` (the July 2026
audit, all PR reviews, specs, decision records — the package's institutional
memory) is **not tracked by git** (`.gitignore:25 dev/*`; only `dev/validation/`
and now `dev/review-2026-08/` are allowlisted). It exists only in the Dropbox
working tree: one bad sync from loss, invisible to fresh clones, unversioned.
MM to decide: allowlist `!dev/audit/` (contents become public on push) or move
to a private archive. Until decided, treat the Dropbox copy as the single
fragile source. Source: review lead, 2026-08-01.

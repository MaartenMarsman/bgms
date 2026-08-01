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
| F-006 | major | open | `vignettes/diagnostics.Rmd` carries the wrong mixture-ESS gloss ("frequent switching = precise / well explored"). Grounded refutation and corrected wording in `plans/a2-review-findings.md` finding 1 (statistic normalizes away PIP-driven flip frequency; penalizes stickiness relative to the chain's own p). Vignette ships with the release. Grep all vignettes for the same gloss. | a2 #1; flagged-issues #5 |
| F-007 | minor | open | `Readme.Rmd` stale vs hand-edited `README.md`; knitting would regress the README. Decide: drop the Rmd (recommended) or re-sync. | flagged-issues #12 |
| F-008 | note | open | `.github/copilot-instructions.md` gap after pkgdown removal: "register a new exported function" procedure should point at the bgms-docs reference workflow. | flagged-issues #13 |
| F-009 | note | decided | Build hygiene: never build/check in the Dropbox tree (stale `src/*.o` against old RcppParallel + generated `src/Makevars` → unloadable `bgms.so`; Dropbox sync resurrects files mid-render). Standing practice: `git archive` to a clean dir outside Dropbox for every build. | flagged-issues #16, #14; A5 §2 |
| F-010 | note | decided | `precision_graph_prior` default STAYS "joint" for 0.2.0.0 (hierarchical path: trust gauge off by default + post-hoc only, silent clamping beyond size-44/22 anchor hulls, minutes-order surface-build pause, continuous/mixed-only). Revisit post-release with the GGM paper; flip would be a breaking change. Decision record: flagged-issues #2b, PR #172 review. | flagged-issues #2b |
| F-011 | note | fixed-verified | pkgdown decommissioned before any release tag: `.github/workflows/pkgdown.yaml` and `_pkgdown.yml` absent from rc1 (verified in merge diff, 2026-08-01). Required because the workflow also fired on published releases. | flagged-issues #4/2; verified in rc1 |
| F-012 | major | open | In-repo CRAN anchor (`cran-0.1.6.3` = `18e660a2`) not yet verified against the actual CRAN 0.1.6.3 tarball. All review diffs assume they match. Verification is task 2 of brief 01. | review lead, 2026-08-01 |
| F-013 | note | open | Local `develop` is 1 commit ahead of `origin/develop` (`50617ae2`); rc1 contains it. Push develop + main + tags together when MM decides to publish the rc; until then all work is local. | review lead, 2026-08-01 |
| F-014 | note | open | Branch hygiene: six unmerged local branches. Verified 2026-08-01 that all release-relevant content is in rc1: `fix/bgmcompare-association-scale` and `feature/indicator-diagnostics-cleanup` are patch-equivalent in develop (`git cherry` = 0); `feature/hier-followup`, `feature/hier-gating`, `feature/user-facing-checks` residual diffs are pre-squash/pre-refactor states (develop strictly newer — checked hunks). `feat/graphical-g-prior` (May 2026, +14k/−37k vs develop) is the deliberate post-release GGM line — keep. Delete the five absorbed branches after MM confirms; retire the stale `wt-ufc` worktree on user-facing-checks. | review lead, 2026-08-01 |
| F-015 | note | fixed-verified | The 2026-07-31 scale-fix landing decisions (`dev/audit/2026-07-31-scale-fix-landing-decisions.md`) were executed: (1) fix landed as its own commit, patch-equivalent in develop; (3) the scale-contingency sentence is in `vignettes/comparison.Rmd:70-75`; NEWS.md "Breaking changes" documents the omega*x vs 2*omega*x divergence and the ~½ rescaling of reported pairwise effects. Still to verify: (2) the ~7s cross-implementation guard runs in the every-run test tier (check in report 01's test listing). | review lead, 2026-08-01 |

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

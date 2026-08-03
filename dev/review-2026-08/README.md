# August 2026 pre-CRAN release review (0.2.0.0)

Review of the merged `develop`→`main` state before submitting 0.2.0.0 to CRAN.
CRAN currently has 0.1.6.3.

## Fixed reference points

| Ref | Meaning |
|---|---|
| `cran-0.1.6.3` (tag, `18e660a2`) | Last CRAN release state (in-repo anchor; verified against the CRAN tarball in report 01) |
| `v0.2.0.0-rc1` (tag, `4f5bddab`) | **Frozen evaluation target.** The develop→main merge of 2026-08-01. All checks run against this tag. |

Rules of engagement:

- The evaluation target stays frozen. Fixes land on `develop` via short-lived
  branches; `main` does not move until the evaluation is complete.
- When the evaluation is done: merge `develop`→`main` again, run the lighter
  re-verification pass (full `R CMD check --as-cran`, targeted re-tests of every
  fixed finding), then submit.
- Never build/check inside the Dropbox working tree (stale `src/*.o`, generated
  `Makevars`, Dropbox sync interference — see FINDINGS F-009). Always
  `git archive <tag>` to a clean directory outside Dropbox.

## Team

- **Review lead** (planning session): triage, briefs, synthesis,
  FINDINGS.md, MAINTAINERS.md, final report.
- **Opus agent** (execution session, no shared context): everything that runs
  code — checks, test suite, simulations, recovery studies, benchmarks.
- **MM (PI)**: user-perspective evaluation of new user-facing functionality,
  NEWS.md reconstruction, targeted code reads where statistical judgment is
  decisive. Code-read briefs are sequenced as a curriculum (see PLAN.md).

## Protocol

- Briefs live in `briefs/NN-<topic>.md`, numbered in order of issue. Every brief
  is self-contained (goal, repo/branch/commit, exact commands, deliverable).
- Every executed brief produces `reports/NN-<topic>.md` with the fixed
  structure: *What was done / Findings (each with severity: blocker, major,
  minor, note) / Evidence (commands, output, paths) / Open questions.*
- Reports `00a`–`00c` are the review lead's internal mapping sweeps (no brief).
- `FINDINGS.md` is the deduplicated master list. Only the review lead edits it.
  Severity: **blocker** (must fix before CRAN) / **major** / **minor** / **note**.
  Status: **open** / **fixed-claimed** / **fixed-verified** / **decided** / **wontfix**.
- Maintainability findings are NOT release blockers unless they are also
  correctness or CRAN-policy issues; they go to the backlog in MAINTAINERS.md.

## Deliverables (definition of done)

1. FINDINGS.md with all blockers resolved and verified.
2. Final review report (what changed since 0.1.6.3, what was checked and how,
   residual risks).
3. Clean NEWS.md covering 0.1.6.3 → release, PI-approved.
4. Executed CRAN checklist: `R CMD check --as-cran` clean, win-builder /
   mac-builder, reverse dependencies (easybgm!), Date bump, cran-comments.md.
5. MAINTAINERS.md: architecture map, conventions, dragons, prioritized
   maintenance backlog reconciled with the July 2026 audit (`dev/audit/`).

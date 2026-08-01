# Brief 01 — Baseline verification of the frozen release candidate (Opus agent)

You are the execution agent in a pre-CRAN release review of the R package
**bgms** (Bayesian analysis of graphical models: Rcpp/RcppParallel C++ core, S7/S3
R layer, testthat suite). You share no context with the review lead; everything
you need is in this brief.

## Goal

Establish the executable baseline for the frozen evaluation target: does the
release candidate build, check, and test cleanly — and is our in-repo notion of
"what is on CRAN" actually what is on CRAN?

## Repo and reference points

- Repo (local clone, Dropbox): `/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms`
- Evaluation target: tag `v0.2.0.0-rc1` (commit `4f5bddab`, on `main`)
- CRAN anchor: tag `cran-0.1.6.3` (commit `18e660a2`)
- **Never build or check inside the Dropbox tree** (stale `src/*.o` compiled
  against an older RcppParallel, stale generated `src/Makevars`, Dropbox sync
  interference). Work in `~/bgms-review/` (create it; outside Dropbox).

## Tasks

### 1. Clean exports

```sh
mkdir -p ~/bgms-review
cd ~/bgms-review
git -C /Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms archive v0.2.0.0-rc1 --prefix=bgms-rc1/ | tar -x
git -C /Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms archive cran-0.1.6.3 --prefix=bgms-cran-anchor/ | tar -x
```

### 2. CRAN-anchor honesty check

Download the published CRAN source of 0.1.6.3
(`https://cran.r-project.org/src/contrib/bgms_0.1.6.3.tar.gz`; if 404, check
`https://cran.r-project.org/src/contrib/Archive/bgms/`). Unpack and
`diff -r` its `bgms/` against `~/bgms-review/bgms-cran-anchor/`, ignoring files
generated at build/submission time (`MD5`, `Meta/`, `build/`, `inst/doc/`,
anything `.Rbuildignore`d such as `dev/`, `.github/`). Report **every remaining
difference verbatim**. Purpose: all review diffs use the in-repo anchor; if it
diverges from the real CRAN tarball, the whole diff map has a false baseline
(finding F-012).

### 3. Build and `R CMD check --as-cran` the rc

```sh
cd ~/bgms-review
R CMD build bgms-rc1                      # vignettes on; note build time
R CMD check --as-cran bgms_0.2.0.0.tar.gz
```

Report: R version, platform, compiler (`R CMD config CXX` + version); every
ERROR / WARNING / NOTE **verbatim**; check duration; installed package size;
any compiler warnings during compilation of `src/` (capture the install log —
`00install.out` in the check dir); vignette build times and whether any vignette
runs long MCMC at build time.

### 4. Full test suite, CRAN skips disabled

Install the built tarball into a fresh library and run the full suite with
`NOT_CRAN=true` so `skip_on_cran()` tests run:

```sh
mkdir -p ~/bgms-review/lib
R CMD INSTALL --library=~/bgms-review/lib bgms_0.2.0.0.tar.gz
cd ~/bgms-review/bgms-rc1
NOT_CRAN=true Rscript -e '.libPaths(c("~/bgms-review/lib", .libPaths()));
  library(testthat); library(bgms);
  res <- test_dir("tests/testthat", package = "bgms", reporter = "summary");
  saveRDS(as.data.frame(res), "~/bgms-review/test-results-rc1.rds")'
```

(Adapt mechanics if needed — e.g. `testthat::test_local()` — but the
requirements stand.) Report: totals (pass / fail / warning / skip); **every
failure and every warning verbatim** with test file and name; every skip with
its stated reason, grouped (which tests never run anywhere?); total runtime and
the 10 slowest test files; whether any test is stochastic-flaky (rerun failures
once to distinguish deterministic from flaky, and say which).

### 5. What does CRAN actually exercise?

Rerun the suite WITHOUT `NOT_CRAN` (i.e. as CRAN would). Report the same totals
and the delta: which test files effectively vanish under CRAN settings. This
feeds the review's test-adequacy assessment.

## Deliverable

Write your report to
`/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms/dev/review-2026-08/reports/01-baseline-check.md`
with exactly this structure:

1. **What was done** — commands run, environment (R/compiler/platform versions).
2. **Findings** — numbered, each with severity: `blocker` (must fix before
   CRAN), `major`, `minor`, or `note`.
3. **Evidence** — verbatim output for every finding (trim only obvious
   repetition), paths to logs/artifacts under `~/bgms-review/`.
4. **Open questions** — anything you could not settle, with what would settle it.

Leave all artifacts (tarball, check directory, test logs) in `~/bgms-review/`
and list their paths at the end of the report. Do not modify anything in the
repo except writing this one report file. Do not push, commit, or fix anything —
report only.

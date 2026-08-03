# Report 01 — Baseline verification of the frozen release candidate

Agent: brief 01 (Opus). Date: 2026-08-01.
Target: tag `v0.2.0.0-rc1` = commit `4f5bddab` (verified).
CRAN anchor: tag `cran-0.1.6.3` = commit `18e660a2` (verified).

**Headline.** The release candidate builds, checks and tests cleanly. `R CMD
check --as-cran` returns **Status: 2 NOTEs**, no ERROR and no WARNING; both
NOTEs are known and neither is a code defect. The full test suite passes
**0 failures / 0 errors / 0 warnings** in every configuration run (in-check,
`NOT_CRAN=true`, and CRAN settings). The in-repo CRAN anchor is **faithful to
the published CRAN tarball** — F-012 can be closed. The substantive findings
are about what the release *ships* and what it *does not exercise*, not about
whether it works.

---

## 1. What was done

### 1.1 Environment

| | |
|---|---|
| R | 4.6.0 (2026-04-24), platform `aarch64-apple-darwin23` |
| OS | macOS Tahoe 26.4.1 (Darwin 25.4.0), arm64 |
| Cores | `parallel::detectCores()` = **15** |
| C++ compiler used for the package | **Apple clang 21.0.0 (clang-2100.1.1.101)**, SDK `MacOSX26.5.sdk` |
| R itself was compiled by | Apple clang 17.0.0 (clang-1700.3.19.1), GNU Fortran (GCC) 14.2.0 |
| `R CMD config CXX` | `clang++ -arch arm64 -std=gnu++20` |
| `R CMD config CXXFLAGS` | `-falign-functions=64 -Wall -g -O2` |
| C++ standard requested by the package | `CXX_STD = CXX20` (`src/Makevars.in`) — **unchanged from CRAN 0.1.6.3** |
| pandoc | none on `PATH`; used the RStudio bundle via `RSTUDIO_PANDOC` (see 1.3) |

All work was done in `~/bgms-review/` (outside the Dropbox tree), from
`git archive` exports. Nothing was built or checked inside the working copy.

### 1.2 Commands

```sh
mkdir -p ~/bgms-review && cd ~/bgms-review
git -C <repo> archive v0.2.0.0-rc1  --prefix=bgms-rc1/         | tar -x
git -C <repo> archive cran-0.1.6.3  --prefix=bgms-cran-anchor/ | tar -x

# task 2 — CRAN anchor honesty
curl -sO https://cran.r-project.org/src/contrib/bgms_0.1.6.3.tar.gz
tar -xzf bgms_0.1.6.3.tar.gz -C cran-tarball
diff -r cran-tarball/bgms bgms-cran-anchor            > anchor-diff-full.txt
Rscript dcf-compare.R                                 > dcf-compare.txt   # normalized DCF field compare

# task 3 — build + check
export RSTUDIO_PANDOC=/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64
env PATH="$RSTUDIO_PANDOC:$PATH" R CMD build bgms-rc1                  > build-rc1.log
env PATH="$RSTUDIO_PANDOC:$PATH" _R_CHECK_CRAN_INCOMING_REMOTE_=false \
    R CMD check --as-cran bgms_0.2.0.0.tar.gz                          > check-rc1.log

# task 4/5 — test suite, both tiers, into a fresh library
R CMD INSTALL --library=~/bgms-review/lib bgms_0.2.0.0.tar.gz
NOT_CRAN=true Rscript run-suite.R notcran     > suite-notcran.log
env -u NOT_CRAN Rscript run-suite.R cran      > suite-cran.log
NOT_CRAN=true Rscript run-suite.R notcran2    > suite-notcran2.log   # determinism re-run
Rscript delta.R                               > delta.txt
```

`run-suite.R` loads the **installed** package from `~/bgms-review/lib` and runs
`testthat::test_dir("~/bgms-review/bgms-rc1/tests/testthat", package = "bgms")`,
so both tiers see the same files and the only variable is `NOT_CRAN`.

### 1.3 One deviation worth recording

The first `R CMD build` **failed outright** — no system `pandoc`, so all five
vignettes failed to knit and the build aborted (`Error: Vignette re-building
failed.`). It succeeded only after putting the RStudio-bundled pandoc on
`PATH`. This is an environment property of this machine, not a package defect,
but it means *a plain `R CMD build` in this checkout does not work* without
`RSTUDIO_PANDOC` set. See finding 11.

---

## 2. Findings

| # | Severity | Summary | Maps to |
|---|---|---|---|
| 1 | note | In-repo CRAN anchor is faithful to the published tarball — **F-012 resolved** | F-012 |
| 2 | blocker | `Date` field over a month old → guaranteed CRAN NOTE (confirmed empirically) | F-003 |
| 3 | major | `tests/compliance/` (2.85 MB, 66% of the tarball) ships and **never executes** under `R CMD check` | F-017 |
| 4 | major | Examples are **not** capped at 2 cores: 76 uncapped fit calls across 26 Rd files run 4 chains; `_R_CHECK_LIMIT_CORES_` is ignored on the sampler path | F-018 |
| 5 | minor | 5 tests can never run anywhere — their fixture directory is gitignored and absent | new |
| 6 | minor | 3 more tests can never run *from the shipped tarball* — fixtures are `.Rbuildignore`d | new |
| 7 | note | 8 test files are 100% dormant without `BGMS_RUN_SLOW_TESTS`, including **all SBC certification** | new |
| 8 | note | CRAN-mode test time is 67 s — comfortably inside budget; **F-030's concern does not materialise** on wall clock | F-030 |
| 9 | note | Zero compiler warnings across 43 translation units at `-Wall` | new |
| 10 | note | `R CMD check` HTML-Tidy NOTE is an environment artifact, not a package issue | new |
| 11 | note | `R CMD build` fails without pandoc on `PATH`; all 5 vignettes run live MCMC at build time (no `eval=FALSE`, no cache) | new |
| 12 | note | F-015 item (2) verified: the association-scale guard runs in the every-run tier | F-015 |

Detailed evidence for each is in §3.

---

## 3. Evidence

### Finding 1 — CRAN anchor is faithful (note; **closes F-012**)

`bgms_0.1.6.3.tar.gz` is live on the CRAN main contrib directory (HTTP 200,
741,191 bytes, `last-modified: Sat, 14 Feb 2026 16:50:13 GMT`), not in Archive.

`diff -r` between the published tarball's `bgms/` and `~/bgms-review/bgms-cran-anchor/`
produced 33 differing files. **Every one of them is a build-time artifact.**

*Structural differences — all expected:*

```
Only in bgms-cran-anchor: .Rbuildignore      # dev-only
Only in bgms-cran-anchor: .github            # .Rbuildignore'd
Only in bgms-cran-anchor: .gitignore         # dev-only
Only in bgms-cran-anchor: Meta               # local build artifact, .Rbuildignore'd
Only in bgms-cran-anchor: Readme.Rmd         # .Rbuildignore'd
Only in bgms-cran-anchor: _pkgdown.yml       # .Rbuildignore'd
Only in bgms-cran-anchor: bgms.Rproj         # .Rbuildignore'd
Only in bgms-cran-anchor: dev                # .Rbuildignore'd
Only in bgms-cran-anchor/tests/testthat: fixtures   # .Rbuildignore'd
Only in bgms-cran-anchor/vignettes: .gitignore
Only in cran-tarball/bgms: MD5               # generated at submission
Only in cran-tarball/bgms: build             # generated at build
Only in cran-tarball/bgms/inst: doc          # generated at build
```

*Content differences — 32 of 33 are a missing trailing newline only.* I
re-hashed every differing file with trailing newlines stripped; 32 of 33 hash
identically. `R CMD build` appends the final newline that the repo files lack:

```
diff -r cran-tarball/bgms/src/bgm/bgm_helper.cpp bgms-cran-anchor/src/bgm/bgm_helper.cpp
251c251
< }
---
> }
\ No newline at end of file
```

The affected files (all `src/`, all newline-only):
`bgm/bgm_helper.cpp`, `bgm/bgm_logp_and_grad.{cpp,h}`, `bgm/bgm_sampler.{cpp,h}`,
`bgmCompare/bgmCompare_helper.{cpp,h}`, `bgmCompare/bgmCompare_logp_and_grad.{cpp,h}`,
`bgmCompare/bgmCompare_sampler.{cpp,h}`, `bgmCompare_interface.cpp`,
`bgm_interface.cpp`, `math/custom_exp.cpp`, `math/explog_switch.h`,
`mcmc/mcmc_adaptation.h`, `mcmc/mcmc_hmc.h`, `mcmc/mcmc_leapfrog.{cpp,h}`,
`mcmc/mcmc_memoization.h`, `mcmc/mcmc_nuts.{cpp,h}`, `mcmc/mcmc_rwm.h`,
`mcmc/mcmc_utils.h`, `priors/sbm_edge_prior.{cpp,h}`,
`sbm_edge_prior_interface.{cpp,h}`, `utils/common_helpers.h`,
`utils/print_mutex.h`, `utils/progress_manager.h`, `utils/variable_helpers.h`.

*The 33rd is `DESCRIPTION`, and it too is build-time only.* A normalized DCF
field-by-field comparison (collapsing whitespace, since `R CMD build` re-wraps
continuation lines) shows **no field whose value differs**:

```
Fields only in CRAN tarball: NeedsCompilation, Packaged, Author, Repository, Date/Publication
Fields only in repo anchor :
--- comparison complete ---
```

All five extra fields are added by `R CMD build`/CRAN at submission.

**Conclusion: `R/`, `man/`, `NAMESPACE`, `data/`, `src/`, `tests/` and
`vignettes/` are byte-identical modulo trailing newlines. The in-repo anchor is
a sound baseline; every review diff built on it stands.**

Artifacts: `anchor-diff-full.txt`, `anchor-diff-classified.txt`, `dcf-compare.txt`.

---

### Finding 2 — stale `Date` (blocker; confirms F-003)

`R CMD check --as-cran`, verbatim:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Maarten Marsman <m.marsman@uva.nl>’

The Date field is over a month old.
```

`DESCRIPTION` has `Date: 2026-03-26`. This is one of the two NOTEs and the only
one attributable to the package. One-line fix in the submission commit, as
F-003 already prescribes. Confirmed empirically here.

---

### Finding 3 — `tests/compliance/` ships 2.85 MB of never-executed payload (major; confirms F-017, with corrections)

Size breakdown of `tests/` inside the built tarball (uncompressed member bytes):

```
     2.85 MB    37 files  tests/compliance
     0.85 MB    82 files  tests/testthat
     0.01 MB     2 files  tests/fixtures
     0.00 MB     1 files  tests/testthat.R
--- tarball is 4.31 MB compressed ---
```

`tests/compliance/` is **66% of the shipped `tests/` payload**, dominated by 35
`.rds` bitwise fixtures (largest: `cmp_boredom_nuts_bernoulli.rds` 277 KB,
`cmp_wenchuan_nuts_no_diffsel.rds` 247 KB).

**It never runs under `R CMD check`.** `R CMD check` executes only `*.R` files
at the top level of `tests/`; it does not descend into subdirectories. The check
log confirms only one test file ran:

```
* checking tests ...
  Running ‘testthat.R’ [81s/65s]
 [81s/65s] OK
```

`tests/compliance/test_compliance.R` and `generate_fixtures.R` are both copied
into the check directory and both ignored. The harness is driven separately by
`.github/workflows/weekly-compliance.yaml` (`Rscript tests/compliance/test_compliance.R`),
which runs from a git checkout and does not need the tarball copy.

`tests/fixtures/` also ships, but it is only
`tests/fixtures/generate_legacy_fixtures.R` — an orphan generator (12 KB in the
tree, 2 entries in the tarball).

**Two corrections to F-017 as written:**

1. `tests/testthat/_problems/` **does not exist in rc1** — not in the tarball
   and not in the `v0.2.0.0-rc1` tree. That part of F-017 is already resolved;
   nothing to delete.
2. The bloat does **not** trip any CRAN check on this platform. Both size checks
   passed:
   ```
   * checking installed package size ... OK
   ```
   Installed size is 4.6 MB (`libs/bgms.so` 2.1 MB, `doc/` 880 KB, `R/` 700 KB),
   under the 5 MB NOTE threshold; the tarball is 4.31 MB. So this is a hygiene
   and honesty issue — shipping a retired bitwise harness that no user can run —
   rather than a submission blocker. The lead may want to re-grade F-017 from
   `blocker` to `major` on that basis; it remains worth fixing, and
   `^tests/compliance$` + `^tests/fixtures$` in `.Rbuildignore` would recover
   ~2.86 MB.

---

### Finding 4 — examples are not capped at 2 cores (major; **corrects F-018**)

F-018 states "Examples now capped at 2". **They are not.** At tag
`v0.2.0.0-rc1`:

- `bgm()` signature: `cores = parallel::detectCores()`, `chains = 4`
  (`R/bgm.R:495`, `R/bgm.R:494`). Same for `bgmCompare()` (`R/bgmCompare.R:207`)
  and `simulate_predict.R:127`.
- **76 fit calls across 26 Rd files** call `bgm()`/`bgmCompare()` with no
  `chains` argument, so they run 4 chains. Only `man/simulate.bgms.Rd:96`
  passes `cores = 2`.

Per-file counts of uncapped example fit calls:

```
bgm.Rd 1, bgmCompare.Rd 1, calibration_check.Rd 1, coef.bgmCompare.Rd 4,
coef.bgms.Rd 6, extract_centrality.Rd 1, extract_log_odds.Rd 2,
extract_main_effects.Rd 2, extract_partial_correlations.Rd 2,
extract_precision.Rd 2, plot.bgmCompare.Rd 4, plot.bgms.Rd 5,
plot.bgms_calibration.Rd 1, plot.bgms_centrality.Rd 1, plot_edge_posterior.Rd 5,
predict.bgmCompare.Rd 1, predict.bgms.Rd 2, print.bgmCompare.Rd 4,
print.bgms.Rd 6, simulate.bgmCompare.Rd 1, simulate.bgms.Rd 2,
simulate_mrf.Rd 2, summary.bgmCompare.Rd 4, summary.bgms.Rd 6,
summary.bgms_centrality.Rd 1, verdicts.Rd 5
```

**Measured, not inferred.** A default `bgm()` call — exactly the form in
`?calibration_check` — sustains ~4 concurrent workers, and setting the env var
`R CMD check` uses to signal its core limit changes nothing:

```
########## DEFAULT (no check limit) ##########
detectCores()           : 15
_R_CHECK_LIMIT_CORES_   : ‘’
formals(bgm)$cores      : parallel::detectCores()
formals(bgm)$chains     : 4
   user  system elapsed
 28.944   0.175   7.395
user/elapsed ratio      : 3.91

########## WITH _R_CHECK_LIMIT_CORES_=TRUE ##########
detectCores()           : 15
_R_CHECK_LIMIT_CORES_   : ‘TRUE’
   user  system elapsed
 29.176   0.173   7.447
user/elapsed ratio      : 3.92
```

The check's own example pass shows the same at aggregate:

```
* checking examples with --run-donttest ... [500s/173s] OK
```

500 s CPU against 173 s elapsed is a **2.9× mean** parallel ratio over the whole
examples pass, peaking at 4.

The package *does* have the right guard, but only on one path.
`R/correction_tables.R:151-160` honours the limit:

```r
normalize_builder_cores = function(cores) {
  cores = max(1L, as.integer(cores))
  if(identical(.Platform$OS.type, "windows")) return(1L)
  check_limit = Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if(nzchar(check_limit) && !identical(tolower(check_limit), "false")) {
    cores = min(cores, 2L)
  }
  min(cores, parallel::detectCores())
}
```

The sampler path does not: `R/run_sampler.R:365` passes `nThreads = s$cores`
straight through from the signature default.

**Scope — the other two contexts are already clean, so this is examples-only:**

- *Vignettes* use a show-uncapped / run-capped pattern and are compliant. Each
  displays an unadorned call under `eval=FALSE` and runs a capped one under
  `include=FALSE`, e.g. `intro.Rmd`:
  ```
  ```{r, eval=FALSE}
  fit = bgm(data, seed = 1234)
  ```
  ```{r, include=FALSE}
  fit = bgm(data, seed = 1234, chains = 2, display_progress = "none", verbose = FALSE)
  ```
  ```
  Same in `checking-your-model.Rmd`, `comparison.Rmd`, `diagnostics.Rmd`.
  Measured: `* checking re-building of vignette outputs ... [78s/34s]` = 2.3×.
- *Tests* do not sustain more than ~2 workers: 188.9 s CPU / 149.7 s elapsed
  (`NOT_CRAN=true`) and 81.8 s / 67.4 s (CRAN settings) — ratios of 1.26 and
  1.21.

CRAN's Repository Policy is explicit that a package must never use more than two
cores simultaneously in its checks. This is a known reviewer trigger; the
examples are where it bites.

---

### Finding 5 — 5 tests can never run anywhere (minor; new)

`test-simulate-predict-regression.R:212-222` resolves its fixtures to
`<pkg>/dev/fixtures/scaffolding/`:

```r
golden_fixture_path = function(id) {
  fixture_dir = file.path(
    testthat::test_path(), "..", "..", "dev", "fixtures", "scaffolding"
  )
  file.path(fixture_dir, paste0(id, ".rds"))
}
has_golden_fixtures = function() file.exists(golden_fixture_path("manifest"))
```

That directory does not exist in the `v0.2.0.0-rc1` tree, does not exist in the
working copy, and is not in git — `dev/fixtures` was removed by commit
`d55f3916` ("chore: untrack dev/ directory from git"), and `.gitignore:25`
(`dev/*`, with only `!dev/validation/` re-included) keeps it out. There is no
generator script for it in the repo.

So these five skip in *every* configuration — locally, in PR CI, in nightly CI
(where `BGMS_RUN_SLOW_TESTS=true` is set), and on CRAN:

```
• golden fixtures not found (5): 'test-simulate-predict-regression.R:225:3',
  'test-simulate-predict-regression.R:252:3',
  'test-simulate-predict-regression.R:279:3',
  'test-simulate-predict-regression.R:305:3',
  'test-simulate-predict-regression.R:332:3'
```

They report as skips, which reads as "deliberately deferred" rather than "dead".
Either restore the fixtures with a tracked generator, or delete the block.

---

### Finding 6 — 3 more tests cannot run from the shipped tarball (minor; new)

`tests/testthat/fixtures` is `.Rbuildignore`d (`^tests/testthat/fixtures`), so
340 KB of fixtures present in git are stripped from the tarball. Consequence
under `R CMD check`, from the in-check run:

```
• zratio_reference.rds not generated (3): 'test-zratio-engine.R:8:3',
  'test-zratio-engine.R:44:3', 'test-zratio-engine.R:97:3'
```

Those same tests **run fine** from the source tree — my standalone runs executed
64 assertions in `test-zratio-engine.R` because
`tests/testthat/fixtures/zratio_reference.rds` (62 KB) is present there. The
same mechanism silently disables the legacy-fit fixture tests in
`test-extractor-functions.R:829,834` (`tests/testthat/fixtures/legacy/`, 16
frozen fits back to v0.1.3).

This is a defensible trade (they are regeneration-backed developer fixtures,
and `dev/validation/regenerate_zratio_reference.R` exists), but it means the
shipped tarball's self-test is weaker than the repo's, and the skip message
("not generated") misdescribes the cause on CRAN — the file was generated, then
excluded.

---

### Finding 7 — 8 test files are entirely dormant by default (note; new)

Independently of `NOT_CRAN`, a second gate `BGMS_RUN_SLOW_TESTS=true` governs 96
tests. It is defined ad hoc in 13 test files (e.g. `test-sbc-ggm.R:26-31`) rather
than in a helper. Files where **no test runs at all** even with `NOT_CRAN=true`:

```
             test-mixed-nuts.R      11 skipped
 test-parameter-recovery-ggm.R       2
         test-sbc-correction.R       6
                test-sbc-ggm.R       6
    test-scaling-diagnostics.R       9
        test-validation-slow.R       3
          test-zratio-cauchy.R       5
             test-zratio-law.R       4
```

This covers **all SBC certification** (`test-sbc-ggm.R`, `test-sbc-correction.R`),
GGM/mixed NUTS-vs-MH correctness, and parameter recovery — i.e. the package's
core statistical-correctness gates. They are not dead: `.github/workflows/nightly-validation.yaml:14`
sets `BGMS_RUN_SLOW_TESTS: true` and runs on a `cron: '0 3 * * 1,4'` schedule
(Mondays and Thursdays). But they run **twice a week on a schedule**, never on a
PR, never in `R CMD check`, and never for a user who installs from source. Worth
stating explicitly in the review's test-adequacy assessment; it is a design
choice, not a defect.

---

### Finding 8 — what CRAN actually exercises (note; answers F-030)

Same harness, same directory, only `NOT_CRAN` varies:

```
                         NOT_CRAN       CRAN      delta
passed                       7954       7294       -660
failed                          0          0          0
error                           0          0          0
warning                         0          0          0
skipped                       103        221       +118
elapsed (s)                 149.7       67.4      -82.2
```

CRAN settings drop **660 assertions (8.3%)** and **118 additional tests** across
34 files. The `skip_on_cran()` skips number 162 in the in-check run.

**Files that vanish entirely under CRAN settings** (ran with `NOT_CRAN=true`,
nothing runs under CRAN):

```
                        file passed.notcran skipped.cran
           test-centrality.R             37            6
     test-mixed-correction.R             23            9
 test-hier-zratio-identity.R              3            9
```

**Files losing the most coverage:**

```
                                 file  notcran   cran  lost
           test-extractor-functions.R      649    419   230
             test-calibration-check.R       66      6    60
                      test-verdicts.R       58     13    45
                           test-bgm.R      161    124    37
                    test-centrality.R       37      0    37
                  test-plot-methods.R       55     27    28
             test-prior-sensitivity.R       58     30    28
 test-prior-inclusion-probabilities.R       30      4    26
              test-mixed-correction.R       23      0    23
                 test-regressions-2.R       31      9    22
                 test-bgm-hier-spec.R       21      1    20
                     test-ggm-gibbs.R       20      4    16
```

**F-030's concern does not materialise on wall clock.** F-030 asked whether
CRAN-mode test time exceeds ~10 minutes and needs a CRAN-tier gate. Measured:
**67.4 s** standalone, and **81 s CPU / 65 s elapsed** inside `R CMD check`. That
is an order of magnitude inside budget. Caveat: CRAN runs on ~2 cores, but the
CRAN-mode CPU/elapsed ratio is only 1.21, so a 2-core farm machine would land
around 80–90 s, still comfortable. **No CRAN-tier gate is needed for time.**

10 slowest test files, CRAN settings (seconds):

```
   test-rb-inclusion-probabilities.R 12.14
 test-zratio-isolated-edge-routing.R  8.96
                test-zratio-engine.R  6.17
                   test-bgmCompare.R  5.44
         test-zratio-surface-build.R  4.42
              test-prior-interface.R  3.78
                          test-bgm.R  2.69
     test-zratio-surface-extension.R  2.18
            test-prior-sensitivity.R  1.92
                      test-methods.R  1.84
```

10 slowest with `NOT_CRAN=true`:

```
         test-zratio-surface-build.R 19.80
   test-rb-inclusion-probabilities.R 12.21
            test-prior-sensitivity.R 10.27
                   test-bgmCompare.R  9.22
 test-zratio-isolated-edge-routing.R  9.11
                          test-bgm.R  8.89
                test-bgm-hier-spec.R  6.81
                test-zratio-engine.R  6.30
                   test-centrality.R  4.63
                test-bb-correction.R  4.58
```

**Failures, warnings, flakiness.** There were **no failures, errors or warnings
to report** in any of the four suite executions (in-check; `NOT_CRAN=true`;
CRAN settings; `NOT_CRAN=true` repeat). Since nothing failed, there was nothing
to re-run to separate deterministic from flaky. To probe determinism anyway I
ran the `NOT_CRAN=true` tier a second time:

```
                 run 1    run 2
passed            7954     7954
failed               0        0
warning              0        0
skipped            103      103
elapsed (s)      149.7    152.5
```

Not just equal in total — **all 77 files matched pass-for-pass, skip-for-skip**
(`files differing between the two runs: 0`). No evidence of stochastic
flakiness anywhere in the suite. Every MCMC-touching test that runs is
seed-pinned. Elapsed differed by 1.9%, which is scheduler noise.

Caveat on scope: this establishes run-to-run determinism on one machine with one
build, not cross-platform or cross-core-count reproducibility — and NEWS already
records that bitwise serial/parallel identity does **not** hold on Windows under
RcppParallel >= 6.0.0.

---

### Finding 9 — clean compilation (note; new)

43 translation units compiled at `-Wall -g -O2` with Apple clang 21.0.0.
**Zero warnings, zero notes.** Filtering `00install.out` for
`warning:|Warning:|note:` returns nothing; the log contains only the 43 compile
commands plus:

```
* installing *source* package ‘bgms’ ...
** this is package ‘bgms’ version ‘0.2.0.0’
** using staged installation
** libs
specified C++20
using C++ compiler: ‘Apple clang version 21.0.0 (clang-2100.1.1.101)’
using C++20
using SDK: ‘MacOSX26.5.sdk’
...
** testing if installed package can be loaded from temporary location
** checking absolute paths in shared objects and dynamic libraries
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (bgms)
```

`R CMD check` also passed every compiled-code stage: `checking compilation flags
in Makevars ... OK`, `checking pragmas in C/C++ headers and code ... OK`,
`checking compiled code ... OK`, `checking for GNU extensions in Makefiles ... OK`.

Note on **F-027** (proposes adding `SystemRequirements: GNU make`): the check's
`* checking for GNU extensions in Makefiles ... OK` did **not** flag
`src/Makevars.in`'s `include sources.mk`. That is not a refutation — `include`
is portable across make implementations, whereas GNU-specific *functions* are
what the check looks for — but the empirical trigger F-027 expected is absent
here. Flagging for the lead rather than resolving it.

---

### Finding 10 — the second NOTE is environmental (note; new)

```
* checking HTML version of manual ... NOTE
Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy.
Please obtain a recent version of HTML Tidy by downloading a binary
release or compiling the source code from <https://www.html-tidy.org/>.
```

macOS ships an ancient `tidy`. This NOTE says the validation was *skipped*, not
that it failed; it will not appear on CRAN's machines. No action.

For completeness, the full check verdict and cost:

```
Status: 2 NOTEs
real 361.64   user 729.62   sys 15.63
```

6.0 minutes wall clock on 15 cores, 730 s CPU. Stage costs:

```
* checking whether package ‘bgms’ can be installed ... [67s/68s] OK
* checking examples with --run-donttest ...           [500s/173s] OK
  Running ‘testthat.R’                                 [81s/65s] OK
* checking re-building of vignette outputs ...         [78s/34s] OK
```

Examples dominate at 500 s CPU. On a 2-core CRAN machine that pass alone would
run ~4 minutes; capping example chains (finding 4) would cut it and fix the
policy issue at the same time.

Slowest examples (elapsed, s):

```
prior_sensitivity_check  21.99   (user 43.27)
predict.bgmCompare       20.60   (user 40.79)
plot.bgmCompare           9.70   (user 35.45)
plot.bgms_calibration     8.04   (user 30.76)
plot.bgms                 7.89   (user 30.55)
extract_centrality        7.86   (user 30.31)
calibration_check         7.82   (user 30.09)
plot.bgms_centrality      7.52   (user 29.42)
plot_edge_posterior       7.43   (user 29.06)
bgm                       6.85   (user 13.59)
predict.bgms              6.67   (user 13.28)
```

31 example topics total: 367.9 s user, 123.1 s elapsed (first pass, before
`--run-donttest`).

---

### Finding 11 — build requires pandoc; vignettes run live MCMC (note; new)

`R CMD build` **fails** without pandoc:

```
--- re-building ‘comparison.Rmd’ using rmarkdown
Error: processing vignette 'comparison.Rmd' failed with diagnostics:
Pandoc is required to build R Markdown vignettes but not available. Please make sure it is installed.
--- failed re-building ‘comparison.Rmd’
...
SUMMARY: processing the following files failed:
  ‘checking-your-model.Rmd’ ‘comparison.Rmd’ ‘diagnostics.Rmd’
  ‘intro.Rmd’ ‘prior-sensitivity.Rmd’
Error: Vignette re-building failed.
Execution halted
```

With `RSTUDIO_PANDOC` on `PATH` the build succeeded in **1341.8 s wall clock
(22.4 min)**. Most of that is the package compile that `R CMD build` performs to
process help pages, not the vignettes: `R CMD check` re-built all five vignettes
in 78 s CPU / 34 s elapsed.

All five vignettes **do** run live MCMC at build time. None sets `eval = FALSE`
globally or uses a cache; every `opts_chunk$set()` call is limited to `collapse`,
`comment`, `fig.width`, `fig.height`. They are cheap because each capped fit uses
`chains = 2` and `Wenchuan[, 1:5]`-scale data, which is the right pattern —
worth keeping in mind before anyone enlarges a vignette example.

---

### Finding 12 — F-015 item (2) verified (note)

F-015 left open: "verify the ~7 s cross-implementation guard runs in the
every-run test tier". It does. `test-bgmCompare.R:337` ("bgmCompare pairwise
effects are on the association scale") carries **no** `skip_on_cran()` and ran in
both tiers:

```
under NOT_CRAN=true:   test-bgmCompare.R  ... passed 2  skipped FALSE  real 2.207
under CRAN settings:   test-bgmCompare.R  ... passed 2  skipped FALSE  real 2.194
```

Measured cost is 2.2 s, not ~7 s. Its companion
`test-mixed-mrf-simulate-predict.R` ("sample_mixed_mrf_gibbs: continuous marginal
SD matches association scale") also runs in both tiers (0.03 s).

Two further association-scale assertions **do** skip on CRAN, which the lead may
want to note: `test-regressions-2.R` ("GGM selection summary is a mixture summary
on the association scale") and `test-regressions.R` ("GGM summary() pairwise
means match coef() (association scale)").

---

## 4. Open questions

1. **F-017 severity.** The 2.85 MB `tests/compliance/` payload trips no CRAN
   check on this platform (installed size 4.6 MB, tarball 4.31 MB — both under
   threshold, `checking installed package size ... OK`). Should it stay a
   `blocker`? *Settled by:* a decision from the lead; the measurement is done.
   Note `tests/testthat/_problems/` no longer exists, so that clause is already
   satisfied.

2. **Finding 4 — fix shape.** Capping `chains` in the 76 example calls, capping
   `cores` there, or teaching the sampler path to honour `_R_CHECK_LIMIT_CORES_`
   the way `normalize_builder_cores()` already does, are three different fixes
   with different user-visible consequences. *Settled by:* MM's decision; F-018
   already frames it as "keep with policy comment, or default to a capped value".
   I did not change anything.

3. **Finding 5 — golden fixtures.** Were `dev/fixtures/scaffolding/*.rds`
   deliberately retired along with the rest of `dev/`, or lost by accident when
   `d55f3916` untracked `dev/`? *Settled by:* whether anyone still wants that
   cross-check against the pre-`bgm_spec()` pipeline. If yes it needs a tracked
   generator; if no, delete the five test blocks.

4. **F-027 (`SystemRequirements: GNU make`).** `checking for GNU extensions in
   Makefiles` passed, so I could not confirm F-027's expected trigger. *Settled
   by:* a check on a platform with non-GNU make, or an r-hub/win-builder run —
   I only had macOS/arm64 here.

5. **Single-platform result.** Everything above is macOS 26.4 / arm64 /
   Apple clang 21 / R 4.6.0. The RcppParallel 6.0.0 Windows reproducibility
   caveat in NEWS is untested here by construction. *Settled by:* win-builder
   and r-hub runs before submission.

---

## 5. Artifacts

All under `~/bgms-review/` (i.e. `/Users/maartenmarsman/bgms-review/`), left in
place. Nothing in the repo was modified except this report file.

| Path | What |
|---|---|
| `bgms-rc1/` | clean `git archive` export of `v0.2.0.0-rc1` |
| `bgms-cran-anchor/` | clean export of `cran-0.1.6.3` |
| `cran-tarball/bgms/` | unpacked published CRAN `bgms_0.1.6.3.tar.gz` |
| `bgms_0.1.6.3.tar.gz` | as downloaded from CRAN (741,191 bytes) |
| `bgms_0.2.0.0.tar.gz` | the built rc tarball (4,524,133 bytes) |
| `bgms.Rcheck/` | full check directory (94 MB) |
| `bgms.Rcheck/00check.log` | check log |
| `bgms.Rcheck/00install.out` | compile log — the source for finding 9 |
| `bgms.Rcheck/bgms-Ex.timings` | per-topic example timings |
| `bgms.Rcheck/tests/testthat.Rout` | in-check test output, full 224-skip listing |
| `lib/` | fresh library holding the installed rc (4.5 MB) |
| `env-info.txt` | R / platform / compiler versions |
| `build-rc1.log` | `R CMD build` log (successful run) |
| `check-rc1.log` | `R CMD check --as-cran` log + timing |
| `install-lib.log` | `R CMD INSTALL` into the fresh library |
| `anchor-diff-full.txt` | raw `diff -r` CRAN tarball vs in-repo anchor |
| `anchor-diff-classified.txt` | per-file newline-only vs genuine classification |
| `dcf-compare.R`, `dcf-compare.txt` | normalized DESCRIPTION field comparison |
| `run-suite.R` | suite runner (installed pkg, source-tree tests) |
| `analyse-results.R` | failure/warning/skip/timing extractor |
| `delta.R`, `delta.txt`, `delta-byfile.csv` | NOT_CRAN vs CRAN comparison |
| `probe-cores.R` | the finding-4 measurement |
| `suite-notcran.log`, `test-results-notcran.rds`, `analysis-notcran.txt` | `NOT_CRAN=true` run |
| `suite-cran.log`, `test-results-cran.rds`, `analysis-cran.txt` | CRAN-settings run |
| `suite-notcran2.log`, `test-results-notcran2.rds` | determinism re-run |
| `skips-*.csv`, `byfile-*.csv` | per-test skip reasons and per-file totals |

**Repo state note.** At the time of writing, the working copy is checked out on
`main`, not `develop`. `dev/review-2026-08/` is tracked on `develop` only, so I
created `dev/review-2026-08/reports/` and wrote this file there as an untracked
path (`dev/*` is gitignored on both branches, and `develop` does not track
`01-baseline-check.md`, so switching back to `develop` will restore the tracked
review files alongside this one without conflict).

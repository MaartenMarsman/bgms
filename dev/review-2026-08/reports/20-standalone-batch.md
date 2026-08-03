# Report 20 — standalone batch: compliance harness, sparse-coding proof, extractor gap

Brief 20. Code changes authorized. Branch `fix/standalone-batch`, based on
`origin/develop` = `3317f54a` (contains `8cab3c9b`), 4 commits — **not pushed**.
Worktree `~/bgms-review/wt-fix10`; the Dropbox tree's branch was never switched
and nothing was built in it. No off-limits file was touched: the branch changes
`tests/compliance/{test_compliance.R,generate_fixtures.R}` + its 33 fixture
files, `tests/testthat/test-{simulate-predict-regression,extractor-functions}.R`,
`R/extractor_functions.R`, `man/extract_arguments.Rd` and this report.

`origin/develop` advanced to `53b5a73c` while this batch ran (the lead's own
`FINDINGS.md` / `PLAN.md` commits). None of those files are mine; the branch
still rebases cleanly onto them as far as file overlap goes.

**Headline.** All three items landed. F-100 was worse than briefed: the known
stale name was one of **three** skews, the other two masked behind it because
the first one aborts the script before any compare fixture is reached. The
harness is now green at 32/32. F-109 is proven at both ends — the tag-side
miscode is reproduced with the wrong category assignment printed. F-114 took the
"provide the field" branch; the brief's premise that `bgm()` fits populate
`main_effect_indices` is incorrect (they do not, and do not need to), which is
stated in §2.

Two findings beyond the brief are worth the lead's eye: the compliance fixtures
are **not reproducible across machines** (§2, finding 5), and the harness's
`migrate_fixture()` rename table mapped every name to itself (§2, finding 4).

---

## 1. What was done

| # | Item | Commit | Type |
|---|---|---|---|
| 1 | F-100 compliance harness | `4386f545` | `fix:` |
| 2 | F-109 sparse-coding proof + regression test | `aa9cb9cb` | `test:` |
| 3 | F-114 `main_effect_indices` for compare fits | `83782eaa` | `fix:` |

### Item 1 — F-100, `4386f545`

**(a) The sweep.** Every property name the two harness files read off a fit
object, checked against `S7::props()` of the live `bgms_class` and
`bgmCompare_class`. Neither file sources a helper; `library(bgms)` and
`readRDS()` are their only external reads, so the sweep is complete over the
harness.

| Name read | On | Status |
|---|---|---|
| `posterior_summary_{main,pairwise,indicator}` | bgm | ok |
| `posterior_mean_{main,pairwise,indicator}` | bgm | ok |
| `raw_samples$`{`main`,`pairwise`,`indicator`} | bgm | ok |
| `nuts_diag` + `$treedepth,$divergent,$energy,$ebfmi` | both | ok |
| `posterior_mean_allocations` | bgm | ok |
| **`posterior_coclustering_matrix`** | bgm | **SKEW 1** → `posterior_mean_coclustering_matrix` |
| `posterior_summary_{main,pairwise}_{baseline,differences}` | compare | ok |
| `posterior_summary_indicator` | compare | ok |
| `posterior_mean_{main,pairwise}_baseline` | compare | ok |
| `posterior_mean_{main,pairwise}_differences` | compare | present now, **absent in the 0.1.6.3 baseline — SKEW 3** |
| **`posterior_mean_indicator`** | compare | **SKEW 2** — exists on neither side |
| `raw_samples` (whole list) | compare | ok |

**(b) The fixes, and the run.** Names corrected in both files; fixtures
regenerated against CRAN 0.1.6.3 (`Rscript tests/compliance/generate_fixtures.R`,
32/32 OK, sequential); suite run once against this tree. Verdict line:

```
=== Results: 32 PASS, 0 FAIL, 0 ERROR, 0 SKIP (of 32) ===
All fixtures match (structure-only pending association-scale fixture regeneration).
```

(exit status 0; log `scratchpad/compliance2.log`.)

Also in this commit, all consequences of the sweep rather than separate ideas:
`check_structure()` now checks the coclustering matrix and the allocations
vector so the renamed property is actually exercised (it was extracted but never
compared, since every config is `structure_only`); `migrate_fixture()` removed
(finding 4); `generate_fixtures.R` pins the 0.1.6.3 baseline instead of
installing whatever CRAN currently serves (finding 6); header notes 11 and 12
added.

**(c) Forward-compatibility flag — nothing to coordinate.** I looked for an
assertion over an exact property or field set. **There is none, and in
particular the harness never reads `fit$arguments` at all** — so the other
batch's category recode map and per-group support table can be added to the
compare arguments without touching this harness. Concretely:

- `compare_fields()` and `check_structure()` iterate explicit allow-lists
  (`test_compliance.R:498-513`, `637-652`); unlisted fields are ignored.
- `strip_diag_cols()` uses `intersect()` on column names
  (`test_compliance.R:494`), so extra summary columns are tolerated.
- Nothing calls `names(fit)`, `setdiff()` over a field set, or compares
  `.field_names`.

The one strict whole-object comparison is `identical(exp_val, act_val)` at
**`tests/compliance/test_compliance.R:606`**, which the compare path reaches for
the `raw_samples` list. It would break on any addition *to `raw_samples`* — not
to `arguments` — and it is currently unreachable, because
`structure_only_ids = names(all_configs)` routes every config to
`check_structure()` instead. Flagging it only so that whoever eventually
re-enables bitwise comparison (header note 8) knows it is there. **Not
pre-changed.**

### Item 2 — F-109, `aa9cb9cb`

Evidence in §3. No source change was needed: current develop already recodes
through the stored map. The regression test is section 8 of
`tests/testthat/test-simulate-predict-regression.R` — the file that owns
end-to-end simulate/predict behaviour, alongside the existing section 7
original-scale tests. No compare test file touched.

### Item 3 — F-114, `83782eaa`

**Branch taken: provide the field.** The information is in the object:
`R/build_spec.R:580-592` builds the per-variable main-effect block layout and
`R/build_output_compare.R:46` files it under `cache$main_effect_indices`. Only
`$arguments` never carried it, so `extract_arguments()` returned `NULL`.

**Correcting the brief's premise.** `bgm()` fits do **not** populate
`main_effect_indices` either — the field is built only by the bgmCompare spec
path, and there is no bgm shape to mirror. Nor should there be:
`extract_main_effects.bgms()` already returns one row per variable, whereas
`extract_main_effects.bgmCompare()` returns a flat baseline block whose
per-variable widths vary (one column per category for ordinal, two for
Blume-Capel) and from which categories unobserved in a group are dropped — which
is exactly why the index map is needed there and nowhere else. This is the same
parameter-dropping that broke the first log-PL run in report 16 §3.6.

I surfaced it in `extract_arguments.bgmCompare()` from the cache rather than
storing it into `$arguments` at build time. Reasons: the value is derived, not
an argument the user passed; `$arguments` is where callers look; and it keeps
this change clear of the compare argument builder that another batch is
extending. Fits whose cache predates the field fall through unchanged (no error,
field simply absent). Documented in `man/extract_arguments.Rd`.

---

## 2. Findings

**1. `posterior_coclustering_matrix` → `posterior_mean_coclustering_matrix`
(major, fixed).** `test_compliance.R:428`, `generate_fixtures.R:507` at base.
S7 errors on unknown properties, so this aborted the run at the *first* bgm
fixture — the extraction is outside the `tryCatch` that wraps model fitting, so
it is a hard stop, not a per-config error.

The name was **never** correct. CRAN 0.1.6.3 also stored
`posterior_mean_coclustering_matrix` (`output_utils.R:142` at the tag); only its
roxygen block said `posterior_coclustering_matrix` (`R/bgm.R:332` at the tag),
and the current `R/bgm.R:424` has since been corrected. The harness author
copied the name from the then-wrong documentation. Consequence: on the fixture
side `$posterior_coclustering_matrix` silently returned `NULL` from the S3 list,
so **every committed SBM fixture had a NULL coclustering matrix** — verified on
the pre-change `bgm_wenchuan_nuts_sbm.rds`. The regenerated fixture now carries
a 6×6 matrix.

**2. `posterior_mean_indicator` on the compare path (major, fixed).**
`test_compliance.R:444` and `generate_fixtures.R:526` at base. `bgmCompare` has
no such property in the current S7 class, and 0.1.6.3 set the field only on the
bgm path (`output_utils.R:130`, inside the bgm branch). So it errored against
S7 and was silently `NULL` against every committed fixture. Removed from both
extractors and both field lists; difference inclusion probabilities were already
compared through `posterior_summary_indicator`, so no coverage is lost.

**3. `posterior_mean_{main,pairwise}_differences` have no 0.1.6.3 counterpart
(major, handled).** Only visible once findings 1 and 2 were fixed and the compare
configs became reachable for the first time: **all 13 failed**, each with two
`one is NULL, the other is not` mismatches. CRAN 0.1.6.3's compare output
builder stores posterior means for the **baselines only**
(`output_utils.R:344-386` at the tag); the two `_differences` means are new.
Handled with an explicit `new_since_baseline` skip plus header note 11, rather
than by deleting the fields — they become live comparisons automatically once
fixtures are regenerated against a baseline that has them.

**4. `migrate_fixture()` was a no-op that documented a rename that does not
exist (low, fixed).** Its table mapped `posterior_mean_pairwise` →
`posterior_mean_pairwise` and likewise for the two `_baseline`/`_differences`
variants, and its comment read "PR #84 renamed posterior_mean_pairwise ->
posterior_mean_pairwise". A global search-and-replace of the old name had
clobbered the left-hand side of its own migration table. Confirmed dead:
0.1.6.3 already used `posterior_mean_pairwise` (`output_utils.R:122` at the
tag), so nothing needs migrating. Function and call site removed.

**5. The compliance fixtures are not reproducible across machines (major,
measured, NOT fixed — lead's call).** Regenerating with the same seeds and the
same 0.1.6.3 source does not reproduce the committed fixtures: the chains differ
from iteration 1 (`max |diff|` on `posterior_mean_main` = 0.66; first three
draws `0.4259 -2.0982 -4.7967` committed vs `0.3724 -2.0456 -5.3576`
regenerated). Diagnosis:

- **Not `cores`.** `bgm()` at 0.1.6.3 defaults `cores = parallel::detectCores()`
  and passes it to the parallel sampler as `nThreads`, which looked like the
  obvious culprit. Measured: `cores = 1`, `2` and `detectCores()` give
  bit-identical draws, and `cores = 1` twice is identical. The RNG stream does
  not depend on thread count. **I therefore did not pin `cores`.**
- **Not the data.** `Wenchuan`, `ADHD` and `Boredom` are `identical()` between
  0.1.6.3 and 0.2.0.0.
- **Toolchain/floating point.** A source build of the `cran-0.1.6.3` tag and the
  CRAN binary agree bit-for-bit *on this machine*, and both differ from the
  committed fixtures. Consistent with header note 7: an FP perturbation in the
  gradient cascades through the fixed-step leapfrog integrator.

Consequence: "bitwise compliance" only holds within one machine. This is
harmless today — structure-only comparison is machine-independent, which is why
the suite is green — but it means the committed fixtures now encode *my*
machine's 0.1.6.3 behaviour rather than the previous generator's, and it is a
hard blocker on re-enabling bitwise comparison (header note 8) as long as
fixtures are generated on a developer machine and checked on `ubuntu-latest`.
**Recommendation for the lead:** regenerate the canonical fixtures on the runner
that executes the check (`weekly-compliance.yaml` already has
`workflow_dispatch`), and commit those. Recorded as header note 12 so it is not
rediscovered.

**6. `generate_fixtures.R` did not pin its baseline (medium, fixed).** It ran
`install.packages("bgms", repos = cloud)`, i.e. whatever CRAN currently serves,
and only `cat()`-ed a warning on a version mismatch, tagging the fixtures with
the actual version and continuing. CRAN is still on 0.1.6.3 today, so the run
above was correct — but the moment 0.2.0.0 reaches CRAN, a regeneration would
silently baseline the package against itself and make every comparison vacuous.
Now: resolve the CRAN version first, fall back to the pinned CRAN Archive
tarball if it is not 0.1.6.3, and `stop()` rather than warn if the installed
version is wrong.

**7. Deprecation warning inside the harness (informational, not fixed).** The
bgm configs pass `edge_prior = "Bernoulli"` as a string, which now emits *"The
`edge_prior` argument of `bgm()` must be a prior object as of bgms 0.2.0"*. The
harness wraps fits in `suppressWarnings()`, so it is invisible and does not
affect the verdict; and the string form is deliberate here, since the fixtures
must be generated through 0.1.6.3's API. Noted only so it is not mistaken for a
harness defect later.

**8. `extract_arguments()` returned no `main_effect_indices` for compare fits
(minor, fixed).** See §1 item 3. Severity is minor for users but it is a real
dead end for anyone reconstructing per-variable main effects from a compare fit,
which is precisely what report 16 §7.iv hit.

---

## 3. Evidence

### 3.1 Compliance verdict line (item 1b)

```
=== Results: 32 PASS, 0 FAIL, 0 ERROR, 0 SKIP (of 32) ===
All fixtures match (structure-only pending association-scale fixture regeneration).
```

Before the fix, on the same tree:

```
Error: Can't find property <bgms>@posterior_coclustering_matrix
```

Intermediate state, after fixing skews 1 and 2 only (this is what finding 3 looks
like when it first becomes reachable):

```
=== Results: 19 PASS, 13 FAIL, 0 ERROR, 0 SKIP (of 32) ===
  cmp_wenchuan_nuts_bernoulli:
      posterior_mean_main_differences: one is NULL, the other is not
      posterior_mean_pairwise_differences: one is NULL, the other is not
  ... (all 13 compare configs, identically)
```

### 3.2 F-109 current side (item 2a)

Data: n = 400, 4 ordinal variables, categories coded **1, 2, 4, 5** (gap at 3,
minimum > 0), seeded; `bgm(iter = 300, warmup = 300, chains = 2, seed = 42)`.

The mechanism: `recode_data_for_prediction()`
(**`R/simulate_predict.R:1046-1092`**) maps each value through the stored map
`arguments$category_levels` — for OMRF an unnamed sorted vector of the original
values, recoded category = position − 1 (line 1070) — instead of subtracting the
column minimum, which survives only as the no-map fallback (lines 1081-1088).
The inverse for `simulate()` is `recode_simulated_to_original()`
(**`R/simulate_predict.R:1110-1134`**, OMRF branch at line 1130). Call sites:
`R/simulate_predict.R:673` (bgm) and `:949` (bgmCompare).

```
--- stored recode map (arguments$category_levels) ---
[[1]]
[1] 1 2 4 5
...
num_categories: 3 3 3 3

=== (ii) predict() on ORIGINAL-SCALE sparse newdata ===
newdata (original scale): 1 2 4 5
recoded to internal categories: 0 1 2 3
EXPECTED (fit's own recode)   : 0 1 2 3
min-shift would have given    : 0 1 3 4   <- 0.1.6.3 behaviour
predict() ok; prob columns per variable: 4 4 4 4
all row-probabilities sum to 1: TRUE

=== (i) simulate() -> predict() round trip ===
simulate() returns values:
[1] 1 2 4 5
on the ORIGINAL sparse scale (1,2,4,5): TRUE
round-tripped internal codes: 0 1 2 3  (expected 0 1 2 3)
predict() on simulate() output: no warning, ncol per var = 4 4 4 4
```

### 3.3 F-109 tag side (item 2b) — the miscode, printed

`git archive cran-0.1.6.3 | tar -x` → `R CMD INSTALL` into
`~/bgms-review/lib10-tag`; same data, same seed. The tag's recode is
`x - min(x)` with no map (**`R/simulate_predict.R:1238-1252`** at the tag).
`newdata` spans the full observed range in every column, so the per-column
minimum is 1, exactly as in the training data.

```
bgms version: 0.1.6.3

training data values      : 1 2 4 5
fit num_categories        : 3 3 3 3   -> internal categories 0..3
fit stores a recode map?  : no (category_levels absent)

=== THE MISCODE ===
newdata (original scale), one row per category:
     V1 V2 V3 V4
[1,]  1  1  1  1
[2,]  2  2  2  2
[3,]  4  4  4  4
[4,]  5  5  5  5

0.1.6.3 min-shift recode -> internal categories:
     V1 V2 V3 V4
[1,]  0  0  0  0
[2,]  1  1  1  1
[3,]  3  3  3  3
[4,]  4  4  4  4

CORRECT recode (the fit collapsed 1,2,4,5 to 0,1,2,3):
     V1 V2 V3 V4
[1,]  0  0  0  0
[2,]  1  1  1  1
[3,]  2  2  2  2
[4,]  3  3  3  3

per-value mapping  original -> 0.1.6.3  vs  correct
   1 -> 0   (correct: 0)
   2 -> 1   (correct: 1)
   4 -> 3   (correct: 2)   *** MISCODED ***
   5 -> 4   (correct: 3)   *** MISCODED ***

category 4 does not exist in a 4-category model (valid: 0..3)

=== consequence for predict() ===
predict() returned without error or warning.
control newdata min-shifts to the intended codes 0,1,2,3: TRUE

P(V1 = cat | V2,V3,V4) for each of the 4 rows:

-- MISCODED (0.1.6.3 on original-scale sparse newdata) --
       cat_0   cat_1   cat_2   cat_3
[1,] 0.17955 0.25389 0.26944 0.29712
[2,] 0.21754 0.27312 0.25735 0.25198
[3,] 0.30371 0.30061 0.22331 0.17237
[4,] 0.35017 0.30774 0.20298 0.13911

-- CORRECT (same conditioning categories) --
       cat_0   cat_1   cat_2   cat_3
[1,] 0.17955 0.25389 0.26944 0.29712
[2,] 0.21754 0.27312 0.25735 0.25198
[3,] 0.25919 0.28893 0.24173 0.21015
[4,] 0.30371 0.30061 0.22331 0.17237

max |difference| : 0.04646
rows that differ  : 3 4

=== 0.1.6.3 simulate() scale ===
simulate() returns values: 0 1 2 3
  training data was on scale: 1 2 4 5
  -> simulate() output is NOT on the original scale
```

Read the two tables together: the MISCODED row 3 (original value **4**, coded 3)
returns **exactly** the CORRECT row 4 distribution — the off-by-one across the
gap, made visible. MISCODED row 4 (original value **5**, coded 4) conditions on
a category that does not exist in a four-category model and still returns a
number, silently. The last block is the second half of the same bug: 0.1.6.3's
`simulate()` emits internal codes, so its own round trip does not return to the
scale the user supplied.

### 3.4 F-114 before / after

On the cached `bgmCompare` test fit (4 binary variables, 2 groups):

```
BEFORE (stored arguments, unchanged by this batch): NULL
AFTER  (extract_arguments):
     [,1] [,2]
[1,]    0    0
[2,]    1    1
[3,]    2    2
[4,]    3    3
num_categories: 1 1 1 1
baseline block width: 4  total main cols: 8
```

Four binary variables → one column each, so `mei[nrow, 2] + 1 = 4`, which equals
the baseline block width of `extract_main_effects()`; the remaining 4 of the 8
`raw_samples$main` columns are the differences. `$arguments` itself is
deliberately unchanged, so `fit$arguments$main_effect_indices` is still `NULL`
(§6 question 3). `bgm()` fits are untouched: they never had the field and the
accessor's bgm method is unmodified.

### 3.5 Gates

| Gate | Result |
|---|---|
| Local default-tier suite (`devtools::test()`) | **0 failures, 0 warnings, 0 errors**; 101 skips, exit 0 |
| Compliance suite against this tree | `32 PASS, 0 FAIL, 0 ERROR, 0 SKIP`, exit 0 |
| Tag-side miscode demo | quoted in §3.3 |
| `devtools::document()` | clean, writes only `extract_arguments.Rd` |
| NAMESPACE drift | none (`git diff NAMESPACE` empty) |



---

## 4. Proposed NEWS clause (VERBATIM — `NEWS.md` not edited)

Under **Bug fixes**:

```markdown
* `predict()` now recodes `newdata` through the category map stored with the
  fit rather than by subtracting each column's minimum. Under the old rule a
  sparse category coding -- one with gaps, such as values 1, 2, 4, 5 -- was
  miscoded silently: the fit collapses those values to categories 0, 1, 2, 3,
  but the minimum shift mapped them to 0, 1, 3, 4, so every value above a gap
  was attributed to the wrong category and the largest value fell outside the
  fitted category range entirely, with no error and no warning. `simulate()`
  correspondingly returns ordinal data on the original category scale, so the
  `simulate()` -> `predict()` round trip stays on the scale the model was
  fitted on.
```

Optional companion entry, if the lead wants F-114 in NEWS:

```markdown
* `extract_arguments()` now reports `main_effect_indices` for `bgmCompare()`
  fits: the per-variable column blocks of the baseline main-effect parameters.
  The layout was computed during fitting but kept internal, so there was no
  supported way to map main-effect parameters back to variables.
```

F-100 is test infrastructure and needs no NEWS entry.

---

## 5. Gate result

`devtools::test()` on this tree, default tier (neither `BGMS_RUN_SLOW_TESTS` nor
`BGMS_RUN_CERTIFICATION` set), exit code 0:

```
══ DONE ════════════════════════════════════════════════════════════════════════
```

**0 failures, 0 warnings, 0 errors.** 101 skips, every one of them a declared
tier gate or a missing optional fixture: T1 nightly blocks ("Set
BGMS_RUN_SLOW_TESTS=true …"), T2 certification blocks ("Weekly certification
tier (T2) …"), and the five `golden fixtures not found` skips in
`test-simulate-predict-regression.R` that are equally skipped on `origin/develop`.
No skip is attributable to this batch. Both new tests run and pass.

One thing to be aware of when re-running pieces of this by hand: running
`test-extractor-functions.R` through `testthat::test_file()` against an
*installed* bgms reports two errors (`indicator_ess_column`,
`compute_Vn_mfm_sbm` "could not find function"). Those are pre-existing and
environmental — the file reaches for package internals, which are visible under
`devtools::test()`'s `load_all()` and not under a plain `library(bgms)`. They do
not appear in the gate run above.

No `src/` change was made by any of the three items, so the slow tier was not
required by the brief's gate.

---

## 6. Open questions for the lead

1. **Fixture provenance (finding 5).** The committed fixtures are now mine,
   generated on this machine; the previous ones were someone else's. Under
   structure-only comparison this makes no difference and CI stays green. Do you
   want the canonical set regenerated on `ubuntu-latest` via
   `weekly-compliance.yaml`'s `workflow_dispatch` before the release, so the
   baseline lives on the machine that checks it? That is a prerequisite for ever
   restoring bitwise comparison (header note 8).

2. **`raw_samples` strict comparison (item 1c).** `test_compliance.R:606` compares
   the whole compare-side `raw_samples` list with `identical()`. It is dead code
   today because every config is `structure_only`. If a batch adds an element to
   `raw_samples` (not to `arguments` — that is safe), this is where it would
   surface when bitwise comparison is switched back on. Left untouched per the
   brief.

3. **F-114 placement.** I attached `main_effect_indices` in
   `extract_arguments.bgmCompare()` rather than storing it in `$arguments` at
   build time, specifically to stay out of the compare argument builder another
   batch is extending. If you would rather have `fit$arguments$main_effect_indices`
   and `extract_arguments()` agree exactly, moving it is a two-line change in
   `build_output_compare.R` — better done at that batch's merge than now.

4. **Compare-side coverage gap, noted not acted on.** With
   `posterior_mean_indicator` gone, the compare configs compare four posterior-mean
   fields, five summaries and `raw_samples`. Two of those means
   (`*_differences`) are skipped against the 0.1.6.3 baseline. If you want the
   difference means genuinely covered before release, that needs a second,
   current-version fixture baseline — a different job from this harness, which is
   specifically a 0.1.6.3 compliance check.

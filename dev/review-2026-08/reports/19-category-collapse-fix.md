# Report 19 — the F-075 fix: category collapse in bgmCompare (Opus agent)

Branch `fix/category-collapse`, worktree `~/bgms-review/wt-fix9`, based on
`origin/develop` @ `e76a1e82` (report 16 present; `d7c2fe39` an ancestor).
Seven commits, no pushes, no attribution trailers. Everything measured on
exports/installs under `~/bgms-review/val19/`.

**Verdict in one line.** The release gate is closed: `bgmCompare()` now retains
the union of the categories observed across groups, the group-1 slope drops
from 1.222 to ≈ `bgm()`'s, and the ~4–5% healthy-edge offset that F-074 was
chasing **vanishes** — it was collapse contamination, exactly as report 16's §4
hypothesised.

---

## 1. What was done

| task | commit | what |
|---|---|---|
| 1 semantics + audit | `c3b315a6` | ordinal branch of `collapse_categories_across_groups()` rewritten to the union; downstream-consumer audit in §2 |
| 2 warn + record | `c3b315a6`, `8b6b8ee2` | `warning()` per empty group-by-category cell; `message()` for dropped values; `category_support` in the fit's arguments |
| 3 documentation | `4fa472bd` | new *Categories across groups* Rd section; `vignettes/comparison.Rmd` paragraph |
| 4 NEWS | `4fa472bd` | Bug-fixes entry written for the 0.1.6.3 reader, with refit guidance |
| 5 tests | `adc7c1fd`, `f706848c`, `bc64f7f4` | nine mechanics cases rewritten; report-16 reproducer, true-gap, BC, and `bgmCompare()`-surface cases; T1 estimate-survival case; warning classed and pinned at incidental sites |
| 6 test hook | `adc7c1fd` | `bgmCompare_test_logp_and_gradient()` + `test-bgmCompare-gradient.R` |
| 7 re-measure | — | §4 |
| 8 support scan | — | §5 |

Diffstat against the merge base: 21 files, +1298 / −119 (13 package files
plus eight test files touched only to pin the new warning, §6). Nothing under
`R/plot_*.R`, `R/calibration_check.R`, `R/centrality.R`,
`R/prior_sensitivity.R`, and no NEWS entry but my own was touched.
`origin/develop` advanced by two commits during the batch
(`02be39b4`, `8cab3c9b`); both are review bookkeeping in
`dev/review-2026-08/`, no `R/` or `src/` change, so the base stands.

### The semantics, precisely

`reformat_ordinal_data()` already maps the **pooled** sorted unique values onto
contiguous 0-based codes, so a true gap is closed before the cross-group pass
ever runs. What remained for `collapse_categories_across_groups()` was the
intersection rule, and that is what came out. The ordinal branch now:

* recodes onto the contiguous union of the values observed in **any** group —
  on input from `reformat_ordinal_data()` this is the identity map, and it is
  written to be correct standalone, so a direct call with non-contiguous codes
  still closes gaps and only gaps;
* sets `num_categories[node] = (number of observed values) − 1`;
* builds a `(num_categories + 1) × K` per-group count table and returns it as
  `category_support`;
* raises the two conditions of task 2.

Blume–Capel variables remain exempt, and the rationale is now written down in
the function header and in the Rd: their two parameters are functions of the
numeric category *score*, so renumbering would change the model rather than
relabel it, and an unobserved score is still a meaningful point on the scale.

---

## 2. The downstream-consumer audit (task 1)

The compare sampler was built against intersection-collapsed data, so before
flipping the switch I mapped every consumer of `num_categories` and of the
recoded matrix and asked what each does with a retained category that has zero
observations in some group. Nothing had to be forced; the summary is that a
structural zero is a *statistical* exposure, not a numerical one.

| consumer | zero-support behaviour | verdict |
|---|---|---|
| `compute_counts_per_category()` (`R/compute_utils.R:24`) | `sum(x == category)` returns `0L`; the cell is a legitimate zero | safe |
| `compute_blume_capel_stats()`, `compute_pairwise_stats()` | never indexed by category | not applicable |
| `main_effect_indices` / `pairwise_effect_indices` (`R/build_spec.R:584`) | sized from `num_categories`, which grows; layout is dense over categories, not over *observed* categories | safe, and now larger |
| `log_pseudoposterior()` data term (`bgmCompare_logp_and_grad.cpp:575`) | `counts(c,v) * mu` contributes exactly 0 | safe |
| the normalizer `compute_logZ_and_probs_ordinal_into()` | runs over all `K+1` categories regardless of counts; no `log(0)` anywhere in `src/models/bgmCompare/` | safe |
| gradient, MAIN block (`:626`) | `grad -= sum_col_s`, a strictly negative pull on an empty cell — see below | **identified exposure** |
| `log_pseudolikelihood_ratio_main()` (`:1157`) | multiplies the count, so an empty cell contributes 0 to the ratio | safe |
| threshold start values (`bgmCompare_sampler.cpp:1657`) | `arma::mat main_effects(num_main, num_groups, fill::zeros)` — no empirical-log start anywhere in the compare path | safe |
| adaptive proposal SDs (`:999`) | Robbins–Monro on the realised acceptance rate; adapts to a wide parameter rather than assuming a narrow one | safe |
| NUTS step size / mass matrix | generic over the flat vector | safe |
| `build_output_compare.R`, `generate_param_names_bgmCompare()` | loop `seq_len(num_categories[v])`; a category name is emitted whether or not it was observed | safe (and correct — the parameter exists) |
| `predict.bgmCompare()`, `category_levels` recode map | built from `unique(cbind(x_original, x_recoded))`, so it now records a one-to-one map instead of a many-to-one one | safe, and more informative |

**The one exposure, stated plainly.** For a category `c` that group `g` never
observes, the pseudolikelihood gradient with respect to that group's threshold
is `0 − Σ_i P_i(c)`, which is negative everywhere: the likelihood pushes the
threshold toward −∞ and only the prior stops it. The Cauchy difference prior
does stop it — its gradient falls off as `2/|x|` while the likelihood term
falls off exponentially — so the posterior is proper and the sampler is well
behaved, but the parameter settles far out with a wide posterior. Measured on
a deliberately adverse fit (`out/04_audit_surface.R`, three variables, one
structural zero per group): the two affected threshold differences have
posterior means 5.89 (sd 2.16) and 14.27 (sd 4.67), against 0.31 / −0.50 /
0.62 / 0.12 for the unaffected ones. Nothing non-finite appeared anywhere.

That is why the warning is worded the way it is, and why it says *"expect a
large and very uncertain number there, and do not read it as evidence of a
group difference."* It is a statement about identification, not about a bug.

**End-to-end survival check** (`out/02_audit_zerocell.R`; p = 5, n = 400/group,
4 chains, seven structural-zero cells): all pairwise and main-effect output
finite; pairwise baseline max R̂ 1.001 / min ESS 3109, differences max R̂ 1.005
/ min ESS 639, main block max R̂ 1.007 / min ESS 1404. One NUTS advisory
("energy not stationary in chain 4"), consistent with the heavy tail on a
structural-zero threshold. Compare's group pairwise estimates matched a matched pair of
`bgm()` fits to within 0.14 (group 1) and 0.05 (group 2), on posterior SDs of
about the same size.

**The whole user-facing surface was exercised on a structural-zero fit**
(`out/04_audit_surface.R`): `print()`, `summary()`, `coef()`,
`extract_main_effects()`, `extract_group_params()`, `verdicts()`,
`extract_arguments()`, `predict()` — all clean, no non-finite values in the
`main`, `pairwise`, `main_diff` or `pairwise_diff` summary means.

Nothing in the audit required stopping.

---

## 3. Warning, message, and the recorded support (task 2)

Two conditions, two volumes:

* **`warning()`**, always (not gated on `bgms.verbose`), listing every
  `variable, category, group` triple with no observations — up to ten, then
  "… and N more" — followed by a plain-language statement of the consequence.
  No jargon: no "structural zero", no "identification", no "posterior".
* **`message()`**, gated on `bgms.verbose` like the rest of the data-cleaning
  output, when a category value that no group used is dropped.

The message needed a second commit. Written inside
`collapse_categories_across_groups()` it can never reach a `bgmCompare()`
user, because `reformat_ordinal_data()` runs first and has already renumbered
the values — verified directly (`bgm_spec()` on values {0,1,3} produces
`category_levels` `0→0, 1→1, 3→2` with the collapse function seeing contiguous
input). `8b6b8ee2` moves the user-facing message to `build_spec_compare()`,
which still holds the supplied values in the recode map, and distinguishes a
dropped value from a scale that merely starts above zero — so Boredom (1–7)
and Wenchuan (1–5) stay silent. The in-function message is kept as a guard for
direct calls.

Recorded on the fit: `extract_arguments(fit)$category_support` is a list with
one `(num_categories + 1) × K` integer matrix per ordinal variable
(`NULL` for Blume–Capel), dimnamed `category 0…` × `group 1…`. The recode map
itself was already stored as `category_levels` and is unchanged in shape.

---

## 4. The re-measure (task 7)

Report 16's Phase A re-run on the **same ten data seeds and the same fit
seeds**, on the fixed build (`lib-run19`, frozen at `8b6b8ee2`). `OM`, `MAIN`,
the four planted differences and the construction are byte-identical — reused
from `~/bgms-review/val09/out/item2.rds` as report 16 did. The `bgm()`
baselines are **not refitted**: `bgm()` is untouched by this fix, so `e1`/`e2`
are read out of report 16's saved `A_seed01..10.rds`. Ten compare fits, one
at a time, `chains = 4, cores = 4`. All ten seeds reported.

### 4.1 The operating-characteristics table, side by side

| quantity | report 16 §5 (defective) | **fixed** | two `bgm()` |
|---|---|---|---|
| group-1 pairwise slope vs truth | 1.222 (sd 0.109) | **0.998 (sd 0.029)** | 0.993 (sd 0.025) |
| group-1 slope range | [1.112, 1.418] | **[0.940, 1.028]** | [0.940, 1.020] |
| group-2 pairwise slope vs truth | 1.148 (sd 0.281) | **1.012 (sd 0.033)** | 1.003 (sd 0.064) |
| group-2 slope range | [0.953, **1.926**] | **[0.958, 1.049]** | [0.876, 1.101] |
| slope, collapsed variable's edges excluded | 1.043 / 1.027 | **0.994 / 0.996** | 0.986 / 0.991 |
| slope on the `avoidact` edges alone | 1.313 (level) | **1.002 / 1.027** | 1.001 / 1.015 |
| noise sd on true-zero differences | 0.048 (range 0.007–0.111) | **0.019 (range 0.004–0.044)** | 0.045 |
| noise ratio compare / `bgm()` | 1.11 (range 0.15–2.68) | **0.38 (range 0.10–0.76)** | — |
| max error on a true-zero difference | up to 0.68 | **0.159** | 0.173 |
| δ = 0.40 recovery | 1.59× | **1.10×** | 1.10× |
| δ = 0.40 detection | 10/10 seeds | **10/10 seeds** | — |
| false presences, null edges | 3.4% (14/410) | **2.2% (9/410)** | — |
| worse group is group 1 | 8 of 10 (p = 0.109) | **4 of 10 (p = 0.754)** | — |

Every ordinal variable now keeps all five categories on every seed
(`num_categories` = `4 4 4 4 4 4 4 4 4 4`, against report 16's
`4 3 3 3 3 3 1 4 3 4` where `avoidact` was binary). Each seed raises exactly
one group-support warning, covering 3–9 empty group-by-category cells.

### 4.2 Per seed, none omitted

| seed | cmp g1 (was) | cmp g2 (was) | bgm g1 | bgm g2 | noise cmp (was) | noise bgm | max null err |
|---|---|---|---|---|---|---|---|
| 1 | 0.964 (1.270) | 1.015 (**1.926**) | 0.963 | 1.034 | 0.015 (0.111) | 0.041 | 0.071 |
| 2 | 1.026 (1.112) | 1.024 (1.107) | 1.019 | 1.034 | 0.013 (0.016) | 0.052 | 0.067 |
| 3 | 0.940 (1.331) | 1.020 (1.102) | 0.940 | 1.012 | 0.044 (0.058) | 0.057 | 0.150 |
| 4 | 1.023 (1.121) | 1.044 (1.131) | 0.997 | 1.101 | 0.004 (0.007) | 0.045 | 0.021 |
| 5 | 0.989 (1.125) | 0.960 (0.953) | 0.993 | 0.876 | 0.025 (0.043) | 0.054 | 0.091 |
| 6 | 1.028 (1.273) | 1.033 (1.095) | 1.020 | 1.016 | 0.017 (0.037) | 0.039 | 0.071 |
| 7 | 1.005 (1.299) | 1.029 (1.071) | 1.012 | 1.009 | 0.008 (0.060) | 0.036 | 0.037 |
| 8 | 0.983 (**1.418**) | 0.986 (0.985) | 0.992 | 0.947 | 0.010 (0.096) | 0.042 | 0.041 |
| 9 | 1.007 (1.145) | 0.958 (0.984) | 1.000 | 0.947 | 0.028 (0.037) | 0.047 | 0.159 |
| 10 | 1.012 (1.130) | 1.049 (1.120) | 0.996 | 1.051 | 0.012 (0.016) | 0.035 | 0.072 |

Convergence: max R̂ on the difference block 1.020–1.113, min ESS 76–175; max
R̂ on the baseline block 1.010–1.083. Comparable to report 16's run
(which reached 1.119).

### 4.3 The healthy-edge slope — F-074's simulation side

**The ~4–5% healthy-edge offset is gone.** On the fixed build, with the nine
`avoidact` edges dropped, the compare level slope is **0.994** (group 1) and
**0.996** (group 2), against `bgm()`'s 0.986 and 0.991 — an offset of
**+0.8% and +0.5%**, against report 16's +4.3% / +1.043 on the same metric
and the same seeds. Report 16 §4 offered contamination by the collapsed
variable as the lead hypothesis for that offset, having withdrawn the
induced-prior-width account; the fix removes the contamination and the offset
goes with it. **F-074's simulation side is settled: it was collapse
contamination, not a prior-width or level effect.** The two cells report 16
cut to test it — the null split in the degenerate regime and the λ
dose–response — are no longer needed for this question.

The corroborating number is §4.1's row for the `avoidact` edges alone: those
edges carried a level slope of 1.313 pre-fix and now sit at 1.002 / 1.027,
against `bgm()`'s 1.001 / 1.015. The variable that was destroyed is now
estimated as well as `bgm()` estimates it.

### 4.4 Runtime

Nine of the ten fits ran on an awake machine and took **511.6–612.6 s, mean
557.7 s** (total 5019 s), against report 16's 471–621 s for the same cells.
The union semantics carry no material runtime cost, even though the parameter
vector grows (every variable now has four free thresholds rather than one to
four). Seed 1 logs 35,307 s: the laptop was closed during that fit. That is a
wall-clock artifact, confirmed by direct measurement — 9 h 49 m elapsed
against 40 min 05 s of CPU time consumed, with the process accruing 239 CPU-s
per 60 s of wall while awake (`out/TIMING_NOTE.txt`). Report 16 recorded the
same artifact for its own seed 3.

---

## 5. Real-data support scan (task 8)

Zero fits. For each reachable multi-group dataset, per-variable observed
support by group, after the listwise deletion and the pooled ordinal recode
that `bgmCompare()` would apply. `out/40_support_scan.R`, full output in
`out/support_scan.log`.

| dataset | variables | with differing support | levels the old rule would have destroyed | would warn |
|---|---|---|---|---|
| ADHD, clinical vs control (the comparison `bgms-docs` bakes) | 18 | 0 | 0 | no |
| Boredom, English vs French (shipped two-group; the control row) | 8 | 0 | 0 | no |
| Wenchuan, random null split | 17 | 0 | 0 | no |
| Boredom, random 3-way null split | 8 | 0 | 0 | no |
| Boredom, split at the median total score | 8 | **2** | **2** | **YES** |

Reading. The docs-site ADHD comparison could never have been hit: every ADHD
variable is binary, so union and intersection coincide. Boredom's real
en/fr split and both null splits are clean, which reproduces report 16's
Boredom negative and its Wenchuan observation, and extends it to three groups —
the case where an intersection rule compounds fastest.

The one positive row is the informative one. Splitting Boredom at its median
total score is an ordinary thing for a user to do, and it is exactly the
regime the defect punished: two groups that differ in level. Two of eight
7-point variables (`loose_ends`, `keep_interest`) then have a category one
group never uses, and the pre-fix code would have merged it away. So the
answer to report 16's open question 3 — "how often do real group comparisons
have a variable whose observed support differs by group" — is: **not on the
data bgms ships with its own natural grouping, but readily as soon as the
groups differ in location**, which is what a group comparison is usually for.

---

## 6. Verification gate

**1. Full local suite, both tiers — clean.**

| tier | files | tests | pass | fail | error | warning | skip | elapsed |
|---|---|---|---|---|---|---|---|---|
| default (`NOT_CRAN=true`) | 78 | 1203 | 8320 | 0 | 0 | 0 | 102 | 182.6 s |
| `BGMS_RUN_SLOW_TESTS=true` | 78 | 1203 | 8524 | 0 | 0 | 0 | 66 | 289.3 s |

Tests whose expectation I changed, with the reason:

* **`test-collapse-categories.R`, all nine mechanics cases.** Rewritten for
  the union semantics. Cases 2, 3, 4, 5 and 6 asserted merged categories and
  the "Only one value was observed" error on disjoint supports; under the
  union nothing merges and disjoint supports are legal. These changed because
  the semantics changed, not because they were wrong about the old code — the
  file says so at the top.
* **`test-build-arguments.R`, "Compare build_arguments: all expected field
  names present".** `category_support` added to the whitelist. It is a new
  argument on the fitted object (task 2).
* **14 tests across seven files gained `without_support_warning()`.** They are
  about reproducibility, priors, imputation-row mapping, group-indicator
  coercion, and object structure, and their fixtures are small enough
  (25 rows per group of a five-point scale) that a category is missing in one
  group by construction. The warning is correct there; it is simply not what
  those tests are about. The muffle is class-scoped
  (`bgms_group_support_warning`), so any *other* new warning still fails the
  run. No test that is about support behaviour was muffled — those assert on
  the warning directly.

**2. `R CMD check --as-cran` on a `git archive` tarball — the 2 baseline
NOTEs only.**

```
Status: 2 NOTEs
* checking CRAN incoming feasibility ... NOTE   (The Date field is over a month old)
* checking HTML version of manual ... NOTE      ('tidy' not recent enough)
```

Both are the known baseline (identical to `~/bgms-review/check-fix4.log`).
Examples, `--run-donttest`, `testthat.R`, and **re-building of vignette
outputs** all OK. Log: `gate/check19full.log`, tree `gate/checkdir2/`.

**3. The new tests bite on the pre-fix build.** The task-5 tests pass on the
fixed build and fail on an unfixed `origin/develop` export
(`lib-prefix`, built from `git archive origin/develop`):
**34 failures, 5 errors, 31 passes** (`out/prefix-bite.log`). The two that
matter:

* the report-16 reproducer, `test-collapse-categories.R:335` —
  `Expected result$num_categories[1] to equal 3. Differences: [1] 1 - 3 == -2`
  (the four-category variable comes back binary), and `:337` —
  `Expected result$x to be identical to x. Differences: 200/480 mismatches`;
* the estimate-survival case, `:584` —
  `Expected max(abs(z(cmp1, draws1))) < 3. Actual comparison: 5.0 >= 3.0`,
  with `:585` at 3.81 and `:588`/`:589` (mean |z|) at 2.0 and 1.76.

**4. The re-measure table is complete, ten of ten seeds** (§4.2).

---

## 7. Findings

### 19-1 · `resolved` · F-075 is fixed and the gate's second half is closed

§4. Both group slopes land on `bgm()`'s, the group-1 systematic inflation is
gone, the group asymmetry is gone (worse group 4/10 vs 8/10), and detection is
undamaged. Nothing in the fixed table is worse than report 16's.

### 19-2 · `resolved` · F-074's simulation side was collapse contamination

§4.3. The healthy-edge offset falls from ~+4.3% to +0.8% / +0.5% on the same
seeds and the same metric. This confirms report 16 §4's lead hypothesis and
retires the two cells it had cut to test it.

### 19-3 · `note` · a structural zero is pinned by the prior alone, and the number it produces is large

§2. This is the accepted cost of the maintainer's decision 1, not a defect,
and the whole point of the warning. Worth stating in one place for the record:
a group-by-category threshold difference with no supporting observations
settles where the Cauchy difference prior stops it. Measured posterior means of
5.9 (sd 2.2) and 14.3 (sd 4.7) on a deliberately adverse three-variable fit.
Anyone reading `summary()`'s `main_diff` block on a fit that warned will see
those numbers; the warning tells them not to interpret them, and
`category_support` tells them which rows are affected. If the maintainer wants
them suppressed rather than explained, that is a separate decision — the
honest options are to hard-zero those differences, to shrink them with a
category-specific prior, or to mark them in the summary output. I did not do
any of that, because none was authorised and all three change reported
quantities.

### 19-4 · `minor` · the true-gap message could not reach a `bgmCompare()` user

§3. Task 2's `message()`, written where the brief located it, was dead code
from the entry point: `reformat_ordinal_data()` closes gaps first. Fixed in
`8b6b8ee2` by reporting from `build_spec_compare()`. **The same blind spot
still exists for `bgm()`**: an ordinal variable with an unused interior
category is silently renumbered there too, and nothing says so. That is a
shared-path change touching `bgm()`'s output, outside this batch, and I left
it alone — flagging it as a candidate finding rather than acting on it.

### 19-5 · `note` · the intersection rule made structural zeros impossible

Worth recording because it explains why the compare sampler had never been
exercised in this regime: under the old rule every retained category was by
construction observed in every group, so a zero cell in
`counts_per_category` could not occur. The union rule makes them ordinary. The
audit (§2) and the new gradient hook (task 6, which includes a
group-differing-support case) are what now cover it.

### 19-6 · `note` · the group-support warning is sample-size sensitive

It fired in 14 existing tests with 20–25 rows per group. That is correct
behaviour — small groups genuinely lack support — but it means users fitting
small comparisons will see it routinely. The class
`bgms_group_support_warning` is now part of the API so it can be handled
programmatically; whether the Rd should say so explicitly is a maintainer
call.

---

## 8. Open questions

1. **19-3's presentation.** Should `summary()` mark, or suppress, the
   main-effect difference rows whose group-by-category cell is empty? The data
   to do it is now on the object (`category_support`). Maintainer's call; it
   changes printed output.
2. **F-074's real-data side.** §4.3 settles the simulation side. Report 16's
   09-2 null-split observation and anything F-074 claims about real data are
   untouched by this batch.
3. **`bgm()`'s silent gap renumbering** (19-4). Same class of defect as F-111,
   different entry point, shared function.
4. **Three or more groups.** Report 16's open question 4 asked whether the
   damage grows with K. It is now moot for the defect — the union does not
   compound — but the *structural-zero* exposure does grow with K, since each
   extra group is another chance for an empty cell. The three-group case is
   covered mechanically (`test-collapse-categories.R`) and in the gradient hook,
   but not in a fitted-estimate sense.
5. **`nuts_max_depth` under structural zeros.** The audit fit raised one
   "energy not stationary" advisory and the re-measure's min ESS dipped to 76
   on seed 8. Neither is a failure, but heavy-tailed threshold differences are
   the kind of geometry that eventually costs tree depth. Not measured
   systematically here.

---

## 9. Artifacts

`~/bgms-review/val19/`:

* `lib-fix19/` — rolling install of the fix branch; `lib-run19/` — the frozen
  copy the re-measure ran against (`RUN_COMMIT` = `8b6b8ee2`).
* `lib-prefix/`, `prefix/` — `git archive origin/develop` build used for the
  bite check.
* `out/lib19.R`, `out/20_phaseA_fixed.R`, `out/30_digest_fixed.R` — the
  re-measure, seeds and construction identical to report 16's `20_phaseA.R`.
* `out/F_seed01..10.rds`, `out/digestF.rds`, `out/phaseA_fixed.log`.
* `out/02_audit_zerocell.R`, `out/04_audit_surface.R` — the task-1 audit.
* `out/40_support_scan.R`, `out/support_scan.log` — the task-8 scan.
* `out/prefix-bite.log` — the new tests run against the unfixed export.

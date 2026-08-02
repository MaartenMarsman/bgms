# Report 19 — the F-075 fix: category collapse in bgmCompare (Opus agent)

Branch `fix/category-collapse`, worktree `~/bgms-review/wt-fix9`, based on
`origin/develop` @ `e76a1e82` (report 16 present; `d7c2fe39` an ancestor).
Five commits, no pushes, no attribution trailers. Everything measured on
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
| 5 tests | `adc7c1fd`, `f706848c` | nine mechanics cases rewritten; report-16 reproducer, true-gap, BC, and `bgmCompare()`-surface cases; T1 estimate-survival case |
| 6 test hook | `adc7c1fd` | `bgmCompare_test_logp_and_gradient()` + `test-bgmCompare-gradient.R` |
| 7 re-measure | — | §4 |
| 8 support scan | — | §5 |

Diffstat against the merge base: 13 files, +1200 / −98. Nothing under
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
structural-zero threshold. Compare's group pairwise estimates matched matched
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

*(filled in below)*

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

*(filled in below)*

---

## 7. Findings

*(filled in below)*

---

## 8. Open questions

*(filled in below)*

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

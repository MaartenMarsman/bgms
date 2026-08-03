# Report 22 — the Normal slab default (F-119), the ridge picture, and the zero-support test evaluation

Branch `fix/normal-slab`, based on `origin/develop` at `f60101b6` (at or after
`c2fb48e1`). Worktree `~/bgms-review/wt-fix14`. Nothing pushed. Off-limits files
(`R/methods_bgmcompare.R`, `tests/testthat/test-summary-marking.R`,
`man/summary.bgmCompare.Rd`, `R/validate_data.R`, `NEWS.md`, plot files,
`dev/review-2026-08/` beyond this report and its assets) were not touched;
`~/bgms-review/val16/` and `~/bgms-review/val19/` were read only.

---

## 0. What was done

| task | outcome | commit |
|---|---|---|
| 1 — the switch | `interaction_prior` → `normal_prior(scale = 1)` and `difference_family` → `"Normal"`, after mapping which parameter families each governs (§1). Two internal spec defaults flipped with them; one legacy fallback deliberately left. Roxygen, Rd and the comparison vignette swept. One latent defect caught (F-119-a). | `71c508fb` |
| 2 — tests | Both tiers clean. Two expectations touched, each with its derivation; one of them was a tautology that never tested what it claimed. The planted-δ guard re-measured and **not** widened. | `44dc9c2f` |
| 3 — OC transfer | Seeds 01–03 refitted at the new defaults against report 19's own fits on the same seeds. Paired differences ≤ 6% of the ten-seed spread on every slope and noise quantity. **Three seeds; no escalation**, with the reasoning in §3.3. | (analysis only, `~/bgms-review/val22/`) |
| 4 — ridge figure | Exact log-posterior over a structural-zero cell through `bgmCompare_test_logp_and_gradient()`, prior subtraction verified to 1.4e-12. 95% width in the unidentified direction: open → 7.5 (Cauchy) → 2.9 (Normal). | `1c837714` |
| 5 — zero-support evaluation | (a) 6/6 cells fitted, BF decisive and n₁-scaling on a true shift, one conservative cell; (b) **0 of 16 sampling zeros over-called** at BF 10; (c) `prior_sensitivity_check()` coverage verified, no verdict flips across family or scale. | (analysis only, `~/bgms-review/val22/`) |
| 6 — NEWS | Three edits proposed verbatim (§6): the new entry, plus two shipped statements this change makes false. **`NEWS.md` not edited.** | — |

Headline: **the switch is inert for the operating characteristics report 19
measured, and decisive for the one thing it was meant to fix** — the width of
the unidentified direction at a cell one group never observed.

---

## 1. The slab-coverage map (task 1, before any change)

`bgmCompare()` prices four families of real-valued parameters. Three prior
objects govern them, and `interaction_prior` — the argument the decision names —
governs only one.

| parameter family | governed by | where the prior is applied |
|---|---|---|
| baseline pairwise interactions | `interaction_prior` | `bgmCompare_logp_and_grad.cpp:746` (`interaction_prior.logp`), gradient at `:482`, `:749`; single-parameter Metropolis at `:1082` |
| **pairwise-interaction differences** | `difference_family` + `difference_scale` | `bgmCompare_logp_and_grad.cpp:754`, gradient `:489`, `:757`; Metropolis `:1084` |
| **main-effect (threshold) differences** | `difference_family` + `difference_scale` | `bgmCompare_logp_and_grad.cpp:714` (ordinal) and `:731` (Blume–Capel), gradient `:448`, `:466`, `:470`; per-parameter Metropolis `:917`, `:920` |
| baseline thresholds | `threshold_prior` | `bgmCompare_logp_and_grad.cpp:706`, `:723`, `:908`, `:911` |

Column 0 of `main_effects` / `pairwise_effects` is the overall (baseline) value;
columns 1..K−1 are the group-contrast differences. That is the whole of the
split: `interaction_prior` touches column 0 of the pairwise block and nothing
else.

**So main differences were indeed governed by a different object.** The
difference slab is not a prior *object* at all but a pair of loose arguments,
`difference_family = c("Cauchy", "Normal")` (`R/bgmCompare.R:239` before the
change) and `difference_scale = 1`, unpacked at `R/bgmCompare.R:469-470` into
`difference_prior_type = tolower(difference_family)` and handed to
`create_parameter_prior()` at `src/bgmCompare_interface.cpp:407`. One object
serves both difference families — pairwise and threshold — because
`bgmCompare()` gives them one scale and one family.

Under `difference_selection = TRUE` (the default) this difference prior **is the
slab of the spike-and-slab**. Reading F-119 as covering "the slab default", it
therefore covers `difference_family` as squarely as `interaction_prior`, and per
the brief both were flipped.

Two further sites carry the same default and were flipped with it, both internal
formal defaults that `bgmCompare()` always overrides but that would otherwise
hand a caller who omits the argument the old family:
`R/build_spec.R:401` and `R/bgm_spec.R:374` (`difference_prior_type = "cauchy"`
→ `"normal"`). They are inert for `bgm()`, whose spec path never reads the field.

One site was deliberately **not** flipped: `R/anchor_curve.R:103`,
`spec$prior$difference_prior_type %||% "cauchy"`. `build_spec_compare()` always
populates that field, so the fallback fires only for a spec that predates it —
i.e. a fit from the Cauchy era. Flipping it would misreport those fits.

### What the switch actually is

| | before | after |
|---|---|---|
| `bgmCompare(interaction_prior =)` | `cauchy_prior(scale = 1)` | `normal_prior(scale = 1)` |
| `bgmCompare(difference_family =)` | `c("Cauchy", "Normal")` | `c("Normal", "Cauchy")` |
| `bgm(interaction_prior =)` | `normal_prior(scale = 1)` | unchanged |

`cauchy_prior()` and `difference_family = "Cauchy"` remain fully available and
are unchanged in behaviour.

### A defect the switch would have introduced, caught and fixed

`bgmCompare()`'s deprecation shims for `interaction_scale=` and
`pairwise_scale=` guard on `identical(interaction_prior, cauchy_prior(scale = 1))`
(`R/bgmCompare.R:361`, `:371` before the change) — that is, "the user did not
override the default". Moving the default without moving the guard would have
left both guards permanently false, silently ignoring the deprecated arguments
instead of honouring them. Both now guard on `normal_prior(scale = 1)`, exactly
as `bgm()` does at `R/bgm.R:533`, `:566`. The legacy arguments still construct a
*Cauchy* at the requested scale, preserving 0.1.6.3 behaviour, which mirrors
`bgm()` and is what `NEWS.md:123` documents for that path.

### Documentation swept

* `R/bgmCompare.R` roxygen: `interaction_prior` (now names the Normal as the
  default and states it governs the baseline pairwise interactions **only**),
  `difference_family` (now names the Normal as the default and states it governs
  both difference families and is the spike-and-slab's slab), `pairwise_scale`
  (reworded so the deprecated argument's Cauchy is not read as the default),
  and the `standardize` deprecation message.
* `man/bgmCompare.Rd` regenerated by `devtools::document()`. NAMESPACE unchanged.
* `vignettes/comparison.Rmd`: the paragraph that told readers the two entry
  points ship different baseline priors was true and is now false. Rewritten to
  say they agree, and to say where the difference parameters get their prior
  instead.

Commit: `71c508fb` `fix(compare): default slab Cauchy -> Normal, mirroring bgm() (F-119)`.

---

## 2. Tests under the new default (task 2)

Commit: `44dc9c2f` `test: re-tune the two compare expectations the Normal slab default moves (F-119)`.

### Every touched expectation, with its derivation

**(1) `tests/testthat/test-prior-sensitivity.R:329`** — the suite's only failure
under the new default.

```
- expect_equal(d$family, "cauchy")
+ expect_equal(d$family, "normal")
```

Derivation: `anchor_draws()` on a compare fit returns
`tolower(spec$prior$difference_prior_type)` (`R/anchor_curve.R:103`), and that
field is `tolower(difference_family)` (`R/bgmCompare.R:470`). The fit in this
test is at defaults, so the value is now `"normal"`. Exact string equality — no
tolerance is involved and none was widened.

**(2) `tests/testthat/test-regressions-2.R`, "difference_family selects the
difference prior independently"** — did **not** fail, and that is the finding.

The block claimed to check the default:

```r
# The default is Cauchy.
default = run("Cauchy")
expect_equal(default$raw_samples$pairwise[[1]], cauchy$raw_samples$pairwise[[1]])
```

`run()` took the family as an argument, so `default` was an explicit `"Cauchy"`
run compared against another explicit `"Cauchy"` run. The assertion was a
tautology and would have passed whatever the default was. It has been rewritten
so `run()` forwards `...`, the default fit omits `difference_family`, and the
claim is two-sided:

```r
default = run()
expect_equal(default$raw_samples$pairwise[[1]], normal$raw_samples$pairwise[[1]])
expect_false(isTRUE(all.equal(default$raw_samples$pairwise[[1]],
                              cauchy$raw_samples$pairwise[[1]])))
```

Same seed, same sampler path, so the equality is exact and needs no tolerance.
The same two-sided idiom was added for `interaction_prior`, because
`extract_arguments()` does not surface the slab family for compare fits and the
draws are therefore the only public evidence. The file's header comment was
updated from "defaulting to Cauchy" to "defaulting to Normal (F-119)".

### Not re-tuned, and why — the planted-δ guard

`tests/testthat/test-bgmCompare.R:394` ("bgmCompare recovers a planted group
difference at its planted size") passed unchanged. Its bounds were derived, not
guessed — from a 12-seed spread that reached rmse 0.089, max absolute error
0.125 and slope in [0.813, 1.167] — so the question is whether the switch moves
the realised statistics enough to eat that margin. Re-measured at the test's own
data and fit seeds across all four slab combinations:

| difference slab | interaction slab | rmse | max abs err | slope |
|---|---|---|---|---|
| Cauchy | cauchy | 0.0274 | 0.0378 | 0.9764 |
| Cauchy | normal | 0.0266 | 0.0366 | 0.9762 |
| Normal | cauchy | 0.0300 | 0.0403 | 0.9822 |
| **Normal** | **normal** (new default) | **0.0290** | **0.0377** | **0.9827** |

The switch moves rmse by 0.0016 — about 2% of the spread the bounds were built
from — and moves the slope *towards* 1. The bounds (rmse < 0.15, max < 0.20,
slope in [0.7, 1.4]) keep their full separation from what a mis-scaled
parameterization would produce (slope 0.5 or 2, rmse 0.204 or 0.408). Widening
them would have been unjustified.

**No snapshots were re-recorded.**

### Suite results

Both tiers, run against a fresh install of this branch into a private library
(`~/bgms-review/wt-fix14-lib`), `test_dir()` over the worktree's own
`tests/testthat`:

| tier | tests | pass | fail | error | warning | skip | elapsed |
|---|---|---|---|---|---|---|---|
| `NOT_CRAN=true` (final, after the re-tunings) | 1214 | 8448 | **0** | **0** | **0** | 97 | 609 s |
| CRAN settings (`env -u NOT_CRAN`) | 1214 | 7611 | **0** | **0** | **0** | 253 | 202 s |

The first `NOT_CRAN=true` run, before the re-tunings, had exactly one failure —
`test-prior-sensitivity.R:329` — and no warnings. The pass count rises 8444 →
8448 across the two `NOT_CRAN` runs: +1 from that failure becoming a pass, +3
from the assertions added to `test-regressions-2.R`.

A third gate exists and was **not** run: `BGMS_RUN_CERTIFICATION=true` (the
weekly T2 tier), which holds the SBC certification and the Z-ratio Monte-Carlo
channels. It is outside the brief's "both tiers", but it is the tier most
exposed to a slab change and the lead may want it before merge.

---

## 3. Operating-characteristics transfer (task 3)

Report 16's Phase A construction, unchanged: the same truths from
`~/bgms-review/val09/out/item2.rds`, the same data seeds, the same fit seeds,
the same `iter = 2000, warmup = 2000, chains = 4, cores = 4`. The only edit to
report 19's `20_phaseA_fixed.R` is that `interaction_prior = cauchy_prior(1)` is
dropped, so the fit runs at the shipped defaults — which is the change under
test. `bgm()` is **not** refitted: `e1`/`e2` are read from report 16's saved
`A_seed01..03.rds`, exactly as report 19 did. Seeds 01–03, none omitted.

Scripts and outputs: `~/bgms-review/val22/out/` (`20_phaseA_normal.R`,
`30_digest_normal.R`, `31_paired.R`, `N_seed01..03.rds`, `digestN.rds`,
`paired.rds`, and the three logs). `val16`, `val19` and `val09` were read only.

### 3.1 Per seed, paired against report 19's own fits on the same seeds

| seed | g1 slope N / C | g2 slope N / C | noise sd N / C | δ=0.40 N / C | false pres. N / C | max R̂ diff | min ESS | s |
|---|---|---|---|---|---|---|---|---|
| 1 | 0.964 / 0.964 | 1.012 / 1.015 | 0.0155 / 0.0146 | 0.552 / 0.568 | 1 / 1 | 1.035 | 170 | 2234 |
| 2 | 1.026 / 1.026 | 1.024 / 1.024 | 0.0150 / 0.0127 | 0.417 / 0.416 | 1 / 0 | 1.031 | 151 | 2235 |
| 3 | 0.939 / 0.940 | 1.019 / 1.020 | 0.0449 / 0.0438 | 0.431 / 0.437 | 4 / 4 | 1.035 | 216 | 2489 |

N = the Normal slab (this run); C = the Cauchy slab (report 19's `F_seed*.rds`,
same seeds). Runtimes are wall clock on a contended machine and are not
comparable to report 19's.

### 3.2 The report-19 table rows

| quantity | **Normal (new default)** | Cauchy, same 3 seeds | report 19, all 10 |
|---|---|---|---|
| group-1 pairwise slope | **0.977 (sd 0.045)** | 0.977 (sd 0.045) | 0.998 (sd 0.029) |
| group-2 pairwise slope | **1.018 (sd 0.006)** | 1.020 (sd 0.004) | 1.012 (sd 0.033) |
| slope, `avoidact` edges excluded | **0.991 / 0.983** | 0.992 / 0.983 | 0.994 / 0.996 |
| slope, `avoidact` edges alone | **0.960 / 1.054** | 0.959 / 1.056 | 1.002 / 1.027 |
| noise sd on true-zero differences | **0.025 (0.015–0.045)** | 0.024 (0.013–0.044) | 0.019 (0.004–0.044) |
| noise ratio compare / `bgm()` | **0.48** | 0.45 | 0.38 (0.10–0.76) |
| max error on a true-zero difference | **0.150** | 0.150 | 0.159 |
| δ = 0.40 recovery | **1.17×** | 1.18× | 1.10× |
| δ = 0.40 detection | **3/3 seeds** | 3/3 seeds | 10/10 seeds |
| false presences, null edges | **6/123 (4.9%)** | 5/123 (4.1%) | 9/410 (2.2%) |

The Normal column and the Cauchy column are two fits of the *same three data
sets*, so their difference is the slab and nothing else. Where the Normal column
differs from report 19's ten-seed column, the Cauchy column differs the same
way — that gap is the three-seed subset, not the switch. Seed 3 is the worst of
report 19's ten on several rows and is one of the three.

### 3.3 Escalation: three seeds suffice

The pre-registered mechanical rule — every three-seed value inside report 19's
ten-seed range — **fired**, on `g1`, `g1_kk`, `noise` and `d010`. It fires on
excursions of 0.0011 or less:

| flagged | Normal, seeds 1–3 | Cauchy, same seeds | report 19 band (10) | excursion N / C | slab moves it |
|---|---|---|---|---|---|
| `g1` | 0.9645 0.9645¹ 0.9393 | 0.9639 1.0264 0.9396 | [0.9396, 1.0282] | 0.0003 / 0.0000 | 0.0005 |
| `g1_kk` | 0.9481 1.0401 0.8907 | 0.9469 1.0404 0.8910 | [0.8910, 1.0627] | 0.0002 / 0.0000 | 0.0012 |
| `noise` | 0.0155 0.0150 0.0449 | 0.0146 0.0127 0.0438 | [0.0043, 0.0438] | 0.0011 / 0.0000 | 0.0023 |
| `d010` | 0.1297 0.0709 0.0011 | 0.0998 0.0637 0.0011 | [0.0011, 0.1780] | 0.0000 / 0.0000 | 0.0299 |

¹ seed 2's `g1` is 1.0262; the table lists seeds 1, 2, 3 in order.

**The rule is the wrong test here and I am not escalating on it.** Seeds 1–3 are
three of the same ten that *define* the band, and seed 3 sets the band edge on
`g1`, `g1_kk` and `noise`. A containment test can then fail by an arbitrarily
small amount for a reason that has nothing to do with the slab. The paired
contrast is the test the question actually asks, and it is decisive: the largest
move the switch produces in any continuous quantity, expressed as a fraction of
the between-seed width report 19 measured over ten seeds, is

| | max abs. paired diff | r19 ten-seed width | fraction |
|---|---|---|---|
| group-1 slope | 0.0005 | 0.0886 | 0.6% |
| group-2 slope | 0.0028 | 0.0911 | 3.1% |
| noise sd | 0.0023 | 0.0396 | 5.9% |
| δ = 0.40 recovery | 0.0155 | 0.2728 | 5.7% |
| δ = 0.10 recovery | 0.0299 | 0.1770 | **16.9%** |
| max PIP on a null edge | 0.0978 | 0.6954 | 14.1% |
| false presences (count) | 1 | 4 | 25.0% |

Against **report 16's** spread, which is what the brief names and which is much
wider (group-1 slope [1.112, 1.418], noise 0.007–0.111), nothing is close to
outside. Three seeds suffice, and the full ten were not run.

Two directional effects are small but real and worth recording rather than
rounding away:

* **δ = 0.10 recovery improves.** Mean 0.067 under the Normal against 0.055
  under the Cauchy, on the same three data sets — the Normal shrinks a small
  true difference *less*. δ = 0.05 moves the same way (0.006 vs 0.004), δ = 0.20
  and δ = 0.40 are flat.
* **One extra false presence**, on seed 2: 6/123 against 5/123. One event in a
  count whose per-seed range over report 19's ten seeds is 0–4. Not separable
  from noise at three seeds; flagged, not concluded (F-119-b below).

Convergence is comparable to report 19's (max R̂ on the difference block
1.031–1.035 against its 1.020–1.113; min ESS 151–216).

---

## 4. The ridge visualization (task 4)

Assets, committed at `1c837714`:
`dev/review-2026-08/reports/assets/f119-ridge-panels.png`,
`f119-ridge-profile.png`, `f119-ridge-values.rds`, and the generating script
`f119-ridge.R` (self-contained, `Rscript f119-ridge.R`, ~80 s).

### 4.1 The construction and why the surface is the true one

Two groups, four variables, 250 persons each, four-level ordinal. Both groups
are drawn from the same three-level MRF and group 2's codes are shifted up by
one, so group 1 lives on {0,1,2} and group 2 on {1,2,3} exactly. Under the union
semantics all four categories are retained (`num_categories` = 3 3 3 3, i.e.
three non-baseline categories each), and category 3 of V1 is observed 70 times
in group 2 and **0 times in group 1**. Coupling 0.10 and thresholds (−0.1, −1.2)
are chosen so no category is starved: group 1's V1 runs 73/106/71/0.

Group *g*'s threshold at category 3 is `mu_3 + proj[g] * delta_3` with
`proj = (−0.5, +0.5)`. Group 2's data pin `mu_3 + 0.5 * delta_3`. Group 1's pin
nothing: no observation ever lands in category 3, so the pseudolikelihood
improves monotonically as `mu_3 − 0.5 * delta_3` is driven down. That is the
ridge, and it runs to infinity.

The surface is evaluated on a 201 × 201 grid through the package's own
`bgmCompare_test_logp_and_gradient()`, every other parameter frozen at its
posterior mean from a defaults reference fit (`difference_selection = FALSE`, so
every difference is active and the flat vector is full length). The hook returns
log-likelihood plus log-prior; the log-prior is recomputed in R from the closed
forms in `src/priors/parameter_prior.h` and subtracted.

**That subtraction is verified, not assumed.** The recovered log-pseudolikelihood
cannot depend on which slab the hook was called with, so the hook is called at
the same parameter vector under three family combinations:

```
normal/normal  -2105.7337562930
cauchy/cauchy  -2105.7337562930   (delta  1.4e-12)
cauchy/normal  -2105.7337562930   (delta -4.5e-13)
```

Agreement to 1.4e-12 on a value of order 2×10³. The script `stopifnot()`s this.

### 4.2 What the figure shows

Three same-axes, same-fill-scale, same-contour-level panels over
(`mu_3`, `delta_3`): (a) the pseudolikelihood alone, no prior at all, hence
improper; (b) the posterior under a Cauchy(0, 1) slab on the difference; (c) the
posterior under a Normal(0, 1) slab. (b) and (c) both carry the shipped
beta-prime(0.5, 0.5) threshold prior on `mu_3`, which F-119 does not change, so
they differ from each other **only** in the slab.

The reading, as the width in `delta_3` of the 95% highest-density region
(`qchisq(.95, 2)/2` = 2.996 log units below the maximum):

| panel | mode (`mu_3`, `delta_3`) | 95% region in `delta_3` | width |
|---|---|---|---|
| (a) likelihood alone | −14.74, +23.92 (at the grid edge) | runs off the grid | ≥ 25.0 |
| (b) posterior, Cauchy(0, 1) | −3.82, +2.00 | [+0.40, +7.92] | **7.52** |
| (c) posterior, Normal(0, 1) | −3.61, +1.52 | [+0.24, +3.12] | **2.88** |

Walking the ridge axis itself (group 2's identified combination held at its
fitted value, so only the empty cell's own threshold moves) gives the same
statement in one dimension: 95% widths of 24.96 / 7.79 / 2.88.

The Cauchy's logarithmic tail bends the ridge back but leaves an arm still
inside the 95% cut at `delta_3` = 8 — a group difference of eight logits on a
category one group never observed. The Normal's quadratic tail closes it by 3.
**That factor of 2.6 in the width of the unidentified direction is what F-119
buys.** It is the clearest statement of the change I can make.

---

## 5. Difference-test evaluation on zero-support cells (task 5)

Shared construction (`~/bgms-review/val22/out/50_zero_common.R`): four
variables, four-level ordinal, coupling 0.10, thresholds (−0.1, −1.2) shared,
V2–V4 carrying an unremarkable top category (~23%) identically in both groups,
V1 carrying the cell under study. All fits at the **new default** with
`difference_selection = TRUE` and `main_difference_selection = TRUE` — the
latter is not the shipped default and is required, because without it main
differences are always in the model and there is no inclusion Bayes factor to
report at all. Every fit `iter = 1500, warmup = 1500, chains = 2, seed` fixed.

### 5.0 A structural fact that governs the whole of task 5

**The main-difference indicator is per variable, not per category.**
`compare_anchor_draws()` builds `main_owner = rep(rep(main_indicator, times =
block), times = num_contrasts)` with `block = num_categories`
(`R/anchor_curve.R`), so one indicator gates a variable's *entire* block of
threshold differences. Confirmed on a fitted object: 10 indicator columns for
p = 4 (6 pairwise + 4 main), against 18 difference parameters.

Consequences the maintainer should know:

* The empty cell's difference **cannot be selected on its own**. Every number
  below is a joint test of "does this variable have any threshold difference",
  of which the empty cell is one of three components.
* The test therefore pays a dimension penalty: two well-identified, genuinely
  null thresholds dilute the evidence from the one that shifted. This is the
  most likely explanation of the conservative `rare / n1 = 100` cell in 5a.

### 5.1 (a) TRUE shift — the non-observing group genuinely cannot produce the category

The non-observing group's generating third threshold is set to −50, so its
category-3 probability is zero by construction; the observing group (n = 1000)
carries it at a tuned rate. Realised rates 0.2532 (target 0.25) and 0.0512
(target 0.05). **Zero cells dropped.**

| category rate in the observing group | n₁ (non-observing) | count in observing group | PIP | log BF | BF | verdict |
|---|---|---|---|---|---|---|
| common, 0.253 | 100 | 242 | 1.0000 | **+27.07** | 5.7e+11 | presence |
| common, 0.253 | 400 | 230 | 1.0000 | **+118.29** | 2.4e+51 | presence |
| common, 0.253 | 1600 | 228 | 1.0000 | **+271.70** | 9.9e+117 | presence |
| rare, 0.051 | 100 | 59 | 0.2798 | **−0.95** | 0.39 | undecided |
| rare, 0.051 | 400 | 47 | 0.9996 | **+7.79** | 2.4e+03 | presence |
| rare, 0.051 | 1600 | 44 | 1.0000 | **+35.84** | 3.7e+15 | presence |

Direction and magnitude: the BF points at *presence* in five of six cells and
grows close to linearly in n₁ — roughly 0.17 log units per non-observing
observation when the category is common, 0.023 when it is rare. The control
variables V2–V4, which have no planted difference, point the other way in all
eighteen control readings — every one negative, the largest magnitude 6.13
(`common / n₁=1600`).

The one **undecided** cell is the informative one. With a rare category and a
small non-observing group, seeing 0 of 100 where the observing group's rate is
5.1% is a p ≈ 0.005 event, yet the block test returns BF 0.39 — mild evidence
*against* a difference. Two things drive that: the block dilution of §5.0, and
the fact that under a Normal(0,1) slab a shift large enough to zero out a
category is far into the slab's tail, so the marginal likelihood under the
alternative is penalised. **The test is conservative here, not liberal.**

### 5.2 (b) SAMPLING zero — the over-call risk

Both groups drawn from one distribution whose V1 top category is rare (tuned
rate 0.0203, target 0.02); n ∈ {50, 100} per group; 20 seeded replicates each.
There is no difference to find. Every replicate is classified and none hidden:

| n per group | sampling zero (one group empty) | category absent in both (**dropped**) | no zero |
|---|---|---|---|
| 50 | 8 | **2** | 10 |
| 100 | 8 | **0** | 12 |

The two dropped cells are `n=50 rep 11` and `n=50 rep 15`, where neither group
used the category, so it is removed from the model entirely and no zero-support
cell exists.

**The finding: 0 of 16 sampling zeros were called supported at BF 10.**

| n | sampling zeros | called supported at BF 10 | median log BF | max log BF | max BF |
|---|---|---|---|---|---|
| 50 | 8 | **0 (0%)** | −0.67 | +0.07 | 1.1 |
| 100 | 8 | **0 (0%)** | −0.92 | +1.17 | 3.2 |

All 16 land on *undecided*. The largest Bayes factor any sampling zero produced
was 3.2, against a threshold of 10. The 22 no-zero replicates behave the same
way (0 supported, median log BF −1.39 at n = 50 and −1.90 at n = 100), so the
empty cell is not pulling the verdict even slightly. **Under the new default,
this design shows no over-call risk from a sampling zero.**

The honest limit: 16 sampling zeros bound the over-call rate loosely. Zero out
of 16 is consistent with a true rate up to ~17% at 95% confidence. The result is
"no sign of over-calling", not "over-calling is ruled out".

### 5.3 (c) Slab sensitivity — family, and scale

**Does `prior_sensitivity_check()` cover main-difference verdicts? Yes,
verified.** Not taken from the documentation (which does claim it, at
`R/prior_sensitivity.R:283-287`) but checked on a fitted object:
`compare_anchor_draws()` returns 10 indicator columns of which 4 are
main-difference, and `anchor_reweight()` returns a reweighted inclusion
probability for all 10. So the **scale** axis needs no refit beyond the check's
own anchors.

**What the reweighting cannot do is change the family.** `anchor_log_weights()`
takes `draws$family` from the fit and only ever moves the scale `s` within it
(`R/anchor_curve.R`). The Cauchy-vs-Normal contrast was therefore done by
refitting, same data, same seed.

#### Family axis (refits)

| cell | PIP N | log BF N | PIP C | log BF C | N − C | verdict flips? |
|---|---|---|---|---|---|---|
| common / n₁=100 | 1.0000 | +27.07 | 1.0000 | +30.23 | −3.16 | no |
| common / n₁=400 | 1.0000 | +118.29 | 1.0000 | +176.08 | −57.79 | no |
| common / n₁=1600 | 1.0000 | +271.70 | 1.0000 | +407.61 | −135.91 | no |
| rare / n₁=100 | 0.2798 | −0.95 | 0.1867 | −1.47 | +0.53 | no |
| rare / n₁=400 | 0.9996 | +7.79 | 0.9999 | +11.79 | −4.01 | no |
| rare / n₁=1600 | 1.0000 | +35.84 | 1.0000 | +41.50 | −5.66 | no |

**No verdict flips in any cell.** But the direction is systematic: the Normal
returns *less* evidence than the Cauchy in five of six cells, and the gap widens
with the strength of the evidence (−3 at log BF 27, −136 at log BF 272). This is
the expected consequence of the geometry in §4 — a structural-zero shift wants
to be enormous, the Cauchy's heavy tail assigns it far more prior mass than the
Normal's, and the marginal likelihood under the alternative follows. **A
zero-support difference earns less evidence under the new default than it did
under the old one.** Here that costs nothing, because the surviving evidence is
still overwhelming and the one undecided cell stays undecided under both. It is
worth knowing that the effect exists and its sign.

#### Scale axis (the package's anchored curve, anchors 0.5×/1×/2×, thin per the brief)

The V1 main-difference verdict is **constant across the whole swept range in all
six cells** — `presence` in five, `undecided` in the rare/n₁=100 cell — with
wobble q95 between 0.20 and 0.78 log units.

In four of the six the curve sits flat at log BF **+13.82**. That is the
package's documented display cap, `ln(10^6)`, from clamping the pooled inclusion
probability at 1 − 1e-6 (`R/prior_sensitivity.R:609`, documented at `:1227`).
It means the verdict is saturated across the sweep, not that the evidence is
literally 13.82. The two unsaturated cells run [−1.75, −0.55] (rare/n₁=100,
undecided throughout) and [+5.00, +13.82] (rare/n₁=400, presence throughout).

Two diagnostics from the thin grid, reported rather than suppressed:

* `common / n₁=100`: the 0.5× refit failed the convergence gate (edge-inclusion
  R̂ above 1.01), so the curve does not cover the scales nearest 0.5×.
* `rare / n₁=1600`: the package warned that the 0.5× and 1× anchor radii do not
  overlap, leaving a gap in the curve, and asked for another anchor.

Both are consequences of the deliberately thin three-anchor grid, not of the
fits. Neither affects the verdicts, which are read from the exact per-anchor
statistics rather than from the interpolated curve.

---

## 6. Proposed NEWS entry (task 6, report-only — `NEWS.md` was not edited)

`bgm()`'s own Cauchy→Normal switch is worded in the **Changed defaults.**
section, `NEWS.md:57-67`. Matching it needs **three** edits, not one: the new
entry, plus two existing statements that this change makes false.

**(1) New bullet, to follow the `bgm()` interaction-prior bullet in
**Changed defaults.** (`NEWS.md`, after line 67):**

```
* `bgmCompare()` now defaults to `interaction_prior = normal_prior(scale = 1)`
  for the baseline pairwise interactions and to `difference_family = "Normal"`
  for the group differences, both of which were Cauchy in 0.1.6.3 and through
  0.2.0.0's development; the baseline default now matches `bgm()`, so the two
  entry points no longer ship different priors for an identically named
  argument. Results under the defaults change: a fit at the new defaults is a
  different model from a fit at the old ones, and difference verdicts on
  categories one group never observed move the most, because a Normal slab
  bounds an unidentified threshold difference far more tightly than a Cauchy
  one. `cauchy_prior(scale = 1)` and `difference_family = "Cauchy"` remain
  fully available and reproduce the previous behaviour.
```

**(2) Required correction to the existing `bgm()` bullet (`NEWS.md:59-62`),
which now states the opposite of what ships.** Replace

```
  the same default, so the GGM prior chain matches `bgm()`. `bgmCompare()`
  keeps a Cauchy default, `cauchy_prior(scale = 1)` — the two entry points
  therefore ship different defaults for an identically named argument, which
  matters if you compare a compare fit against separate `bgm()` fits at stated
  defaults. Pass `interaction_prior = cauchy_prior(scale = 2.5)` to `bgm()` for
```

with

```
  the same default, so the GGM prior chain matches `bgm()`. `bgmCompare()`
  takes the same default for its baseline pairwise interactions, so a compare
  fit and separate `bgm()` fits at stated defaults now price a baseline
  interaction alike. Pass `interaction_prior = cauchy_prior(scale = 2.5)` to `bgm()` for
```

**(3) Required correction to the `difference_family` bullet (`NEWS.md:247-249`),
whose parenthesis is now wrong.** Replace

```
* `bgmCompare(difference_family =)` chooses the family of the prior on the
  group differences, `"Cauchy"` (the default, and 0.1.6.3's fixed behaviour) or
  `"Normal"`.
```

with

```
* `bgmCompare(difference_family =)` chooses the family of the prior on the
  group differences, `"Normal"` (the default) or `"Cauchy"` (0.1.6.3's fixed
  behaviour). It governs the pairwise-interaction differences and the
  main-effect threshold differences alike, and under
  `difference_selection = TRUE` it is the slab of the spike-and-slab.
```

---

## 7. Verification gate

| # | gate | result |
|---|---|---|
| 1 | both suite tiers, 0 failures / 0 warnings | **pass** (§2) — `NOT_CRAN=true` 1214 tests, 8448 pass, 0/0/0, 97 skip; CRAN settings 1214 tests, 7611 pass, 0/0/0, 253 skip. Every re-tuned expectation listed with its derivation in §2. |
| 2 | `R CMD check --as-cran` on a `git archive` tarball, 2 baseline NOTEs | **pass, with the baseline restated as 3 on this machine** — see below |
| 3 | `devtools::document()` clean; NAMESPACE unchanged | **pass** — `roxygenise()` rewrote only `man/bgmCompare.Rd`; `git diff NAMESPACE` is empty |
| 4 | task-3 table complete, 3 or 10 seeds stated with reason | **pass** — 3 seeds, none omitted; escalation declined with the paired evidence in §3.3 |
| 5 | ridge assets render and are committed, scripts alongside | **pass** — `1c837714`, four files under `reports/assets/f119-ridge-*` |

### Gate 2 in full

Tarball from `git archive HEAD`, built **with** vignettes (pandoc 3.8.3 from the
RStudio bundle), then `R CMD check --as-cran --run-donttest`:

```
Status: 3 NOTEs
* checking CRAN incoming feasibility ... NOTE   (The Date field is over a month old)
* checking examples ................... NOTE   (examples with CPU time > 5s)
* checking HTML version of manual ..... NOTE   ('tidy' not recent enough)
```

The first and third are the two the brief names. **The third — the examples
timing NOTE — is not mine, and I verified that rather than asserting it.** The
identical build-and-check was run on a `git archive origin/develop` tarball:

```
base   (origin/develop):  Status: 3 NOTEs   examples [483s/247s]
branch (fix/normal-slab): Status: 3 NOTEs   examples [483s/246s]
```

`diff` of the two logs' check-step outcomes shows **no difference in any step's
verdict** — the only differences anywhere are install/example/vignette timings
of 1–3 seconds. The eight examples over the 5 s threshold are the same eight in
both, and half of them (`plot.bgms`, `extract_centrality`, `calibration_check`,
`plot.bgms_calibration`) are `bgm()` examples this change does not touch. The
NOTE is the known F-018 condition (examples run 4 chains uncapped) surfacing on
this machine under `--run-donttest`; report 19 did not pass that flag.

**The branch adds no NOTE, no WARNING and no ERROR relative to its base.**
Examples, `--run-donttest`, `testthat.R` and re-building of vignette outputs all
OK. Logs: `~/bgms-review/f119chk/check.log` (branch),
`~/bgms-review/f119base/check.log` (base).

---

## 8. Findings

| id | severity | finding |
|---|---|---|
| **F-119** | — | **Implemented.** `bgmCompare()` defaults to `normal_prior(scale = 1)` for the baseline pairwise interactions and `difference_family = "Normal"` for both difference families. `cauchy_prior()` and `"Cauchy"` remain available. Commit `71c508fb`. |
| **F-119-a** | **medium** | **A naive implementation of F-119 would have silently broken two deprecated arguments.** The `interaction_scale=` and `pairwise_scale=` shims guard on `identical(interaction_prior, cauchy_prior(scale = 1))` — "the user did not override the default". Moving the default without moving the guard leaves both permanently false and the deprecated arguments quietly ignored, with no warning and no error. Fixed here (`R/bgmCompare.R:361`, `:371`). **Generalisable: any future default flip must move its `identical()` guards with it.** No test covered this; the rewritten `test-regressions-2.R` block now would not catch it either, because it tests the default rather than the shim. Worth a dedicated test. |
| **F-119-b** | **low, watch** | Three-seed OC transfer shows one extra false presence under the Normal (6/123 vs 5/123 on the same data, seed 2), and δ = 0.10 recovery improves (0.067 vs 0.055). Both are inside the between-seed spread and neither is separable from noise at three seeds. If the lead wants either resolved, it needs the full ten. |
| **F-119-c** | **medium** | **The main-effect difference indicator is per variable, not per category** (`R/anchor_curve.R`, `main_owner = rep(rep(main_indicator, times = block), ...)`). One indicator gates a variable's entire block of threshold differences, so a structural-zero category's difference cannot be selected on its own and every main-difference Bayes factor is a block test. Nothing in `?bgmCompare` or the vignette says so, and `main_difference_selection`'s documentation reads as though selection were per difference. This is a documentation gap, not a bug, but it changes how every main-difference verdict should be read. |
| **F-119-d** | **low** | **`difference_prior_type` names two different priors in one file.** In `src/bgmCompare_interface.cpp`, the worker member `difference_prior_type` (line 103) carries the *indicator* prior family (`"Beta-Bernoulli"`, `"Stochastic-Block"`), while `difference_prior_type_str` (line 392) carries the *slab* family (`"cauchy"`, `"normal"`). They are unrelated priors with near-identical names, passed adjacently. A rename would cost nothing and remove a real trap for whoever next touches this path. |
| **F-119-e** | **low** | **Pre-existing test defect, fixed.** `test-regressions-2.R`'s "The default is Cauchy" assertion passed `difference_family = "Cauchy"` explicitly and compared it to another explicit `"Cauchy"` run, so it asserted nothing about the default and did not fail when the default moved. Rewritten two-sided (`44dc9c2f`). |
| **F-119-f** | **medium (docs)** | Two shipped statements were made false by this change and are corrected or proposed: `vignettes/comparison.Rmd`'s "the two entry points ship different defaults" paragraph (**fixed on this branch**), and `NEWS.md:59-62` plus `NEWS.md:247-249` (**proposed verbatim in §6; not edited, per the brief**). The lead must land the NEWS edits or the tag ships a NEWS that contradicts the code. |
| **F-119-g** | **info** | A zero-support difference earns systematically **less** evidence under the Normal slab than under the Cauchy — up to 136 log units in the strongest cell — with no verdict flips in six cells. Expected from the geometry (§4): the shift a structural zero wants is far into the slab's tail, and the Normal penalises it more. Recorded because it is a real behavioural change, and covered by the proposed NEWS wording. |
| **F-119-i** | **info** | The brief's "2 baseline NOTEs" is stale under `--run-donttest`: on this machine both `origin/develop` and this branch return **3**, the third being the examples-timing NOTE (the known F-018 condition — examples run 4 chains uncapped). Verified by running the identical build-and-check on a base tarball; see §7. Worth correcting in the standing gate text so the next agent does not spend a check cycle on it. |
| **F-119-h** | **info** | `R/anchor_curve.R:103`'s `%||% "cauchy"` fallback was deliberately **not** flipped. `build_spec_compare()` always populates the field, so the fallback fires only for specs predating it — Cauchy-era fits. Flipping it would misreport them. |

---

## 9. Open questions

1. **Should `main_difference_selection = TRUE` be reachable per category?** F-119-c
   makes every main-difference verdict a block test over a variable's whole
   threshold-difference block. For the zero-support case this is the difference
   between "does this variable differ" and "is this empty category's threshold
   different", and users comparing groups with unequal category support will want
   the second. Out of scope here; it is a modelling decision, not a defect.
2. **Is `difference_scale = 1` still right under the Normal?** `NEWS.md:40-42`
   already records that the `difference_scale` default is "still under study"
   under the association-scale parameterization. The switch changes what that
   scale means: the same number now buys a much tighter bound on an unidentified
   direction (§4, 95% width 7.5 → 2.9). The scale sweep in §5.3 found no verdict
   moving over 0.5×–2×, so nothing is urgent, but the calibration question is now
   coupled to the family.
3. **Is the 5b bound tight enough?** 0 of 16 sampling zeros over-called is
   consistent with a true over-call rate up to ~17% at 95% confidence. If the
   maintainer wants a usable upper bound rather than "no sign of it", this needs
   roughly 60–100 sampling-zero replicates, which at these fit sizes is a few
   hours.
4. **Does the switch interact with the normalizer correction table?**
   `NEWS.md:64-67` records that under the joint precision-graph specification the
   correction table is keyed on the interaction prior, so a new default cell
   builds and caches a fresh table once. `bgmCompare()` has no precision block, so
   this should not apply — but `R/correction_tables.R:481` gates on
   `interaction_prior_type %in% c("cauchy", "normal")` and I did not exercise the
   compare path through it. Flagging rather than asserting.
5. **Should the deprecated `pairwise_scale=` path still produce a Cauchy?** It
   does, mirroring `bgm()` and preserving 0.1.6.3 behaviour, which is defensible.
   But it now means `bgmCompare(pairwise_scale = 1)` and
   `bgmCompare(interaction_prior = normal_prior(1))` differ, where before they
   agreed. `NEWS.md:123` documents this for `bgm()`; the same sentence is not
   there for `bgmCompare()`.

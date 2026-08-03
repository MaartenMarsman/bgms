# Report 23 — mixed relational certification: mixed vs GGM / OMRF (Opus agent)

Report only; no package code, test, or documentation was changed.

**Headline (task 1 answer, one sentence):** a user cannot cause the mixed
path to run on pure-type data — `bgm()` resolves `model_type` from the
declared variable types and pure-type data always routes to the pure
sampler — so tasks 2 and 3 are **internal-machinery certification**, not
user-facing consistency checks.

**Second headline:** the literal tasks 2 and 3 could not be run. Driving the
mixed machinery on pure-type data does not produce a divergent number; it
produces an error (all-continuous, and all-ordinal without selection) or an
**unbounded spin that never returns** (all-ordinal with selection). That is
finding **F1**, characterized in full in §5. Tasks 2 and 3 were then run on
the nearest reachable configuration — the same pure block plus one
independent companion variable of the other type — and are labelled
throughout as substitutes.

---

## 0. Environment

```sh
git -C <dropbox>/bgms worktree add ~/bgms-review/wt-val23 \
    -b review/mixed-relational-cert origin/develop      # 07813a28
R CMD INSTALL --library=~/bgms-review/lib-val23 ~/bgms-review/wt-val23
```

`origin/develop` at **`07813a28`**, verified to contain `1393da20`
(`git merge-base --is-ancestor 1393da20 origin/develop`). R 4.6.0, bgms
0.2.0.0, darwin/aarch64 (Apple M5 Pro). Nothing was built in the Dropbox
tree and its branch was not switched.

Every fit uses `chains = 4, cores = 4`, `iter = 2000, warmup = 2000`,
`update_method = "nuts"`, and the `bgm()` default priors, one fit at a time
(standing budget: two other agents on the machine). Scripts, logs and RDS
outputs are in `~/bgms-review/val23/`: `harness.R`, `probe.R`, `trace.R`,
`task2.R`, `task3.R`, `task4.R`, `task5.R`, and the matching `.log`/`.rds`.

### Prior mapping (identical on every path, by construction)

The mixed sampler is not reachable through `bgm()` on pure-type data, so
`harness.R::fit_path()` reproduces `bgm()`'s body verbatim with one change:
`model_type` is passed through instead of being pinned to `"omrf"` at
[bgm.R:664](R/bgm.R#L664). Everything downstream — `bgm_spec()`,
`run_sampler()`, `build_output()`, the extractors — is the package's own
wiring. The priors are therefore not hand-matched across paths; they are the
same `bgm()` defaults resolved by the same `unpack_*_prior()` calls:

| prior | value | reaches |
|---|---|---|
| `interaction_prior` | `normal_prior(scale = 1)` | slab on all pairwise, all paths |
| `threshold_prior` | `beta_prime_prior(0.5, 0.5)` | OMRF + mixed discrete block |
| `means_prior` | `normal_prior(scale = 1)` | mixed continuous means |
| `precision_scale_prior` | `exponential_prior(eta = 1)` | GGM + mixed precision diagonal |
| `edge_prior` | `bernoulli_prior(0.5)` | selection runs only |
| `delta` | `NULL` → resolved (below) | GGM + mixed |
| `precision_graph_prior` | `"joint"` | GGM + mixed |

`threshold_prior_type` is read from the same `inputFromR` field with the same
default on both discrete paths
([sample_omrf.cpp:67-68](src/sample_omrf.cpp#L67-L68),
[sample_mixed.cpp:83-84](src/sample_mixed.cpp#L83-L84)) and handed to the
same `create_parameter_prior()`. `delta = NULL` resolves to
`0.5 * log(num_variables)` for `ggm` and `0.5 * log(num_continuous)` for
`mixed_mrf` ([bgm_spec.R:443-451](R/bgm_spec.R#L443-L451)); in the task 3
comparison both sides therefore get `0.5 * log(6)`, and in task 2 both get
`0` (`num_continuous = 1` → `log(1)`).

---

## 1. The routing map

`bgm()` never chooses a model family itself. It calls `bgm_spec()` with
`model_type = "omrf"` hardcoded ([bgm.R:662-664](R/bgm.R#L662-L664)) and
lets the spec builder re-resolve it:

1. [bgm_spec.R:419-429](R/bgm_spec.R#L419-L429) — `validate_variable_types()`
   is called with `allow_continuous = TRUE, allow_mixed = TRUE` (both are
   `model_type != "compare"`), returning `is_continuous` and `is_mixed`.
2. [bgm_spec.R:431-437](R/bgm_spec.R#L431-L437) — the only re-resolution:
   ```r
   if(model_type == "omrf" && is_continuous) model_type = "ggm"
   if(model_type == "omrf" && is_mixed)      model_type = "mixed_mrf"
   ```
3. [bgm_spec.R:571-619](R/bgm_spec.R#L571-L619) — dispatch to
   `build_spec_ggm()` / `build_spec_mixed_mrf()` / `build_spec_omrf()`.
4. [run_sampler.R:27-36](R/run_sampler.R#L27-L36) — `switch(spec$model_type,
   ggm = …, omrf = …, mixed_mrf = …)`, reaching `sample_mixed_mrf()` at
   [run_sampler.R:279](R/run_sampler.R#L279).

So everything turns on `is_mixed`, set in
[validate_model.R:74-76](R/validate_model.R#L74-L76):

```r
has_continuous = any(variable_type == "continuous")
has_discrete   = any(variable_type %in% c("ordinal", "blume-capel"))
is_mixed       = has_continuous && has_discrete
```

Three things close the question:

* **The single-string branch cannot produce `is_mixed`.** When
  `variable_type` has length 1 it is replicated to all variables
  ([validate_model.R:42-64](R/validate_model.R#L42-L64)) and `is_mixed` is
  never touched from its `FALSE` initialization
  ([validate_model.R:40](R/validate_model.R#L40)). `variable_type =
  "ordinal"` → OMRF; `"continuous"` → GGM.
* **The vector branch requires one of each.** `is_mixed` is the conjunction
  above, so an all-`"continuous"` vector gives `is_continuous = TRUE`
  ([validate_model.R:126-134](R/validate_model.R#L126-L134)) → GGM, and an
  all-discrete vector gives both flags `FALSE` → OMRF. A vector mixing
  `"continuous"` with anything not in
  `{ordinal, blume-capel, continuous}` is rejected
  ([validate_model.R:78-86](R/validate_model.R#L78-L86)).
* **There is no second door.** `bgm_spec` is not exported (it is absent from
  `NAMESPACE`, which exports `bgm`, `bgmCompare` and the extractors), and it
  has exactly two callers in the package: `bgm()` with `"omrf"`
  ([bgm.R:664](R/bgm.R#L664)) and `bgmCompare()` with `"compare"`
  ([bgmCompare.R:480](R/bgmCompare.R#L480)). The other route into
  `sample_mixed_mrf()`,
  [extract_prior_inclusion_probabilities.R:239](R/extract_prior_inclusion_probabilities.R#L239),
  takes `num_discrete` / `num_continuous` from an already-fitted mixed
  object, so it inherits that fit's block sizes and cannot manufacture a
  degenerate one.

**Answer.** Pure-type data always routes to the pure sampler; the mixed path
is reachable only with at least one discrete *and* at least one continuous
variable. Tasks 2 and 3 are therefore internal-machinery certification.
Corollary, used throughout §5: `p ≥ 1` and `q ≥ 1` hold for every mixed fit
a user can produce, so **F1 is not user-reachable**.

---

## 2. Task 2 — all-ordinal reduction: mixed vs OMRF

### 2.0 The sharp check could not be run

The brief's construction — one all-ordinal dataset put through the OMRF path
and through the mixed machinery — requires the mixed model at `q = 0`. It
does not run. With `edge_selection = FALSE` it raises
`Col::subvec(): indices out of bounds or incorrectly used` from inside
`sample_mixed_mrf()`; with `edge_selection = TRUE` it never returns. Both are
finding **F1** (§5.1), characterized there in full. Per the brief's
stop-and-characterize rule this is reported as the first-class result of
task 2; what follows is a labelled substitute, not the sharp check.

### 2.1 The substitute and what it can still catch

Nearest reachable configuration: the **same** six ordinal variables plus one
continuous companion that is independent noise, which is the minimal mixed
model the machinery accepts. The ordinal block of that fit is compared with
the OMRF fit on the same six columns.

The substitute is weaker than the sharp check in exactly one way: the mixed
fit carries six extra cross parameters that the OMRF fit does not have, and
their estimates are not exactly zero at finite `n`, so the ordinal block is
perturbed by an amount that has nothing to do with cross-path consistency.
It is not weaker in the way that matters for this review: a convention
divergence in the discrete block — a factor, a constant offset — would
survive the substitution untouched, because the companion perturbs the
ordinal estimates by ~10⁻² and a convention bug moves them by 50–100%.
§2.4 measures the confound rather than assuming it is small.

**Design.** `p = 6` ordinal, 4 categories, `n = 1000`, data seed
`20230801`, threshold ramp `-1.0·c` (with zero thresholds this coupling
strength saturates the top category and one variable loses support
entirely; the ramp restores full support, verified 4/4/4/4/4/4). Planted
graph `V1-V2 = 0.40, V2-V3 = -0.35, V3-V4 = 0.30, V4-V5 = -0.25,
V1-V6 = 0.30`, all other pairs zero. Companion `C1 ~ N(0,1)` from seed
`20230802`, `max |r(C1, V·)| = 0.0275`. Fit seeds: `4001` (both paths),
`4002` (yardstick refit, OMRF path, same data).

**Yardstick** = the OMRF refit spread at seed 4002 vs seed 4001, stated next
to every number it judges.

### 2.2 (a) `edge_selection = FALSE` — posterior means

| edge | OMRF (4001) | OMRF (4002) | mixed (4001) |
|---|---|---|---|
| V1-V2 |  0.39059 |  0.39071 |  0.39282 |
| V1-V3 | -0.02753 | -0.02775 | -0.03035 |
| V2-V3 | -0.32269 | -0.32264 | -0.32141 |
| V1-V4 |  0.02917 |  0.02905 |  0.03160 |
| V2-V4 |  0.01932 |  0.01948 |  0.01793 |
| V3-V4 |  0.35475 |  0.35499 |  0.35689 |
| V1-V5 | -0.00112 | -0.00129 |  0.00084 |
| V2-V5 | -0.01563 | -0.01539 | -0.01664 |
| V3-V5 | -0.08214 | -0.08173 | -0.08062 |
| V4-V5 | -0.28460 | -0.28474 | -0.28609 |
| V1-V6 |  0.22686 |  0.22693 |  0.22818 |
| V2-V6 |  0.02165 |  0.02169 |  0.02150 |
| V3-V6 | -0.01316 | -0.01315 | -0.01279 |
| V4-V6 | -0.01824 | -0.01824 | -0.01893 |
| V5-V6 | -0.01941 | -0.01942 | -0.02032 |

| block | max \|Δ\| mixed vs OMRF | yardstick (OMRF refit) | ratio | slope | intercept |
|---|---|---|---|---|---|
| pairwise (15) | **0.00282** | 0.00041 | 6.9× | 1.00362 | +0.00026 |
| thresholds (18) | **0.02721** | 0.00505 | 5.4× | 1.00811 | +0.00633 |

Thresholds, both paths (`extract_main_effects()`; rows are variables,
columns categories 1–3):

| var | OMRF cat1 / cat2 / cat3 | mixed cat1 / cat2 / cat3 |
|---|---|---|
| V1 | -0.8492 / -1.6537 / -2.2581 | -0.8477 / -1.6590 / -2.2835 |
| V2 | -1.2106 / -2.1858 / -3.3389 | -1.2142 / -2.1966 / -3.3593 |
| V3 | -0.9327 / -1.7425 / -2.7379 | -0.9308 / -1.7435 / -2.7495 |
| V4 | -0.9931 / -2.5381 / -3.7206 | -0.9971 / -2.5492 / -3.7478 |
| V5 | -0.8672 / -1.5924 / -2.8643 | -0.8703 / -1.5975 / -2.8768 |
| V6 | -0.6326 / -1.2870 / -1.8906 | -0.6335 / -1.2919 / -1.9030 |

The six companion cross estimates, planted at exactly zero, come back at
`-0.0375, 0.0241, -0.0258, 0.0276, 0.0205, 0.0103` — the confound, and it is
an order of magnitude larger than the ordinal-block gap it has to explain.

### 2.3 (b) `edge_selection = TRUE` — posterior inclusion probabilities

`bernoulli_prior(0.5)` on both paths (the `bgm()` default), RB estimator.

| edge | planted | OMRF (4001) | OMRF (4002) | mixed (4001) |
|---|---|---|---|---|
| V1-V2 | 0.40 | 1.00000 | 1.00000 | 1.00000 |
| V1-V3 | 0 | 0.04106 | 0.03858 | 0.03798 |
| V2-V3 | -0.35 | 1.00000 | 1.00000 | 1.00000 |
| V1-V4 | 0 | 0.03876 | 0.03800 | 0.04196 |
| V2-V4 | 0 | 0.04573 | 0.04642 | 0.05529 |
| V3-V4 | 0.30 | 1.00000 | 1.00000 | 1.00000 |
| V1-V5 | 0 | 0.02922 | 0.02877 | 0.03209 |
| V2-V5 | 0 | 0.02630 | 0.02884 | 0.02665 |
| V3-V5 | 0 | 0.15699 | 0.15726 | 0.15754 |
| V4-V5 | -0.25 | 1.00000 | 1.00000 | 1.00000 |
| V1-V6 | 0.30 | 1.00000 | 1.00000 | 1.00000 |
| V2-V6 | 0 | 0.05322 | 0.05611 | 0.05446 |
| V3-V6 | 0 | 0.06718 | 0.06565 | 0.06179 |
| V4-V6 | 0 | 0.02083 | 0.02095 | 0.02126 |
| V5-V6 | 0 | 0.02930 | 0.02855 | 0.02687 |

| quantity | max \|Δ\| mixed vs OMRF | yardstick | ratio |
|---|---|---|---|
| PIP (15 ordinal edges) | **0.00956** | 0.00289 | 3.3× |

All five planted edges are selected at PIP = 1.0000 on both paths; the ten
planted zeros agree to ≤0.0096 everywhere. The companion's own cross PIPs
max at 0.0305, i.e. the extra node is correctly read as disconnected.

### 2.4 Is the residual gap a path difference or the companion?

The gaps in §2.2 sit 3–7× above the pure Monte-Carlo yardstick, so
"agrees to MC error" cannot be asserted from those numbers alone. The
confound probe (`task2b_confound.R`) separates the two explanations: refit
the mixed model on the **same** ordinal data at the **same** fit seed, with
the companion column redrawn from a different noise seed. A cross-path
convention difference is a property of the paths and must reproduce; a
companion artefact must move.

| quantity | companion draw 1 | companion draw 2 |
|---|---|---|
| max \|mixed − OMRF\|, pairwise | 0.00282 | 0.00796 |
| max \|mixed − OMRF\|, thresholds | 0.02721 | 0.06885 |
| max \|companion cross estimate\| | 0.03749 | 0.07241 |

* per-edge gap correlation across the two draws: **r = 0.442** (pairwise),
  **r = 0.076** (thresholds);
* `max |gap₁ − gap₂| = 0.00514` on the pairwise block — a fixed convention
  offset would give ≈ 0;
* the gap magnitude tracks the companion's spurious cross magnitude (both
  roughly double from draw 1 to draw 2).

The residual is the companion, not the paths.

### 2.5 Verdict

**Sharp check (as briefed): not runnable — finding F1.**

**Substitute: no systematic divergence.** Slope 1.0036 on pairwise means,
1.0081 on thresholds, intercepts +0.0003 / +0.0063, PIP agreement ≤0.0096
with identical selection decisions on all 15 edges. Nothing resembling a
factor (0.5 or 2), and no constant offset — §2.4 shows what offset there is
does not reproduce across companion draws. The residual excess over the
pure-MC yardstick is attributed, by measurement, to the extra companion
node the substitute is forced to carry. The discrete block of the mixed
sampler and the OMRF sampler are the same estimator to the resolution this
design can reach (~3 × 10⁻³ on associations, ~3 × 10⁻² on thresholds).

---

## 3. Task 3 — all-continuous, calibrated: mixed vs GGM

### 3.0 Substitution, again

All-continuous data through the mixed machinery needs `p = 0`, which fails
in R before reaching the sampler (finding **F1**, §5.1). Same substitution
as task 2, mirrored: six continuous variables plus one **ordinal**
companion (3 categories, drawn independently), compared against the GGM fit
on the same six columns.

`edge_selection = FALSE` throughout. A spike-and-slab shrinks estimates
toward zero by an amount that depends on the model dimension, and the two
sides do not have the same dimension (the mixed side carries the companion
and its six cross edges); that shrinkage would land directly on the slope
this task exists to read. Turning selection off removes it. This is a
deliberate deviation from an unstated default, recorded here.

**Design.** `Q = 6` continuous; truth `K` with
`K12 = -0.35, K23 = 0.30, K34 = -0.25, K45 = 0.20, K16 = -0.30`, unit
diagonal, positive definite; association truth `-K/2`. Data seed
`20230802`, `n ∈ {500, 5000}` from the same truth. Fit seeds `5001`
(both paths), `5002` (GGM yardstick refit). `delta` resolves to
`0.5·log(6)` on both paths (§0).

### 3.1 Slope and intercept across edges

| n | slope (mixed ~ GGM) | intercept | max \|Δ\| assoc | yardstick (GGM refit) | ratio |
|---|---|---|---|---|---|
| 500  | **1.00239** | +1.72 × 10⁻⁴ | 0.00133 | 0.00023 | 5.8× |
| 5000 | **1.00037** | +2.6 × 10⁻⁵ | 0.00017 | 0.00008 | 2.1× |

Per-block max \|Δ\| through `extract_precision()`, against the same refit
yardstick:

| n | block | max \|Δ\| | yardstick | ratio |
|---|---|---|---|---|
| 500  | precision off-diagonal | 0.00267 | 0.00046 | 5.8× |
| 500  | precision diagonal     | 0.00520 | 0.00099 | 5.3× |
| 5000 | precision off-diagonal | 0.00035 | 0.00016 | 2.2× |
| 5000 | precision diagonal     | 0.00054 | 0.00025 | 2.2× |

Both paths at `n = 5000` (association scale, `extract_pairwise_interactions()`):

| edge | truth | GGM (5001) | GGM (5002) | mixed (5001) |
|---|---|---|---|---|
| Y1-Y2 |  0.175 |  0.186229 |  0.186288 |  0.186399 |
| Y1-Y3 |  0     |  0.002654 |  0.002678 |  0.002787 |
| Y2-Y3 | -0.150 | -0.150416 | -0.150442 | -0.150409 |
| Y1-Y4 |  0     |  0.003352 |  0.003273 |  0.003262 |
| Y2-Y4 |  0     |  0.003158 |  0.003189 |  0.003210 |
| Y3-Y4 |  0.125 |  0.122639 |  0.122593 |  0.122649 |
| Y1-Y5 |  0     |  0.009674 |  0.009618 |  0.009554 |
| Y2-Y5 |  0     | -0.008528 | -0.008481 | -0.008635 |
| Y3-Y5 |  0     | -0.003467 | -0.003515 | -0.003461 |
| Y4-Y5 | -0.100 | -0.105369 | -0.105303 | -0.105359 |
| Y1-Y6 |  0.150 |  0.168624 |  0.168655 |  0.168727 |
| Y2-Y6 |  0     |  0.000617 |  0.000547 |  0.000646 |
| Y3-Y6 |  0     |  0.002610 |  0.002612 |  0.002644 |
| Y4-Y6 |  0     | -0.012569 | -0.012503 | -0.012510 |
| Y5-Y6 |  0     | -0.000408 | -0.000368 | -0.000234 |

### 3.2 Shrinking or stable? — the read

**Shrinking, unambiguously.**

| quantity | n = 500 | n = 5000 | factor | implied rate |
|---|---|---|---|---|
| slope − 1 | 2.39 × 10⁻³ | 3.7 × 10⁻⁴ | 6.5× | ≈ n^-0.81 |
| intercept | 1.72 × 10⁻⁴ | 2.6 × 10⁻⁵ | 6.6× | ≈ n^-0.82 |
| max \|Δ\| assoc | 1.33 × 10⁻³ | 1.7 × 10⁻⁴ | 7.8× | ≈ n^-0.89 |
| max \|companion cross\| | 3.56 × 10⁻² | 1.18 × 10⁻² | 3.0× | ≈ n^-0.48 (√n) |

A convention bug of the `omega`-vs-`2·omega` class is a slope at 0.5 or 2.0
that does not move with `n`. The observed slope is 1.0024 at `n = 500` and
1.0004 at `n = 5000` — within 0.24% and 0.04% of 1, and closing at roughly
`1/n`. That is a finite-sample pseudolikelihood-vs-exact-likelihood
discrepancy, plus the residue of the companion confound (bottom row, which
closes at the slower √n rate expected of a spurious estimate).

The two paths also agree with each other far better than either agrees with
the truth, which is the right ordering: at `n = 500` the regression of
estimates on truth has slope 1.16338 (mixed) and 1.16009 (GGM), and at
`n = 5000` 1.04920 and 1.04881 — a shared bias from the shared prior and
determinant tilt, ~500× larger than the between-path difference and
identical on both sides to 3 decimal places. Whatever the two paths do to
the truth, they do together.

### 3.3 Verdict

Calibrated agreement holds; the discrepancy is finite-sample, not
structural. No convention bug in the continuous block. The
`extract_precision()` round-trip (GGM stores raw samples on the precision
scale and converts at
[extractor_functions.R:873-875](R/extractor_functions.R#L873-L875); the mixed model
stores associations directly and converts at
[extractor_functions.R:1748](R/extractor_functions.R#L1748)) lands both
paths on the same scale — which is the specific thing that would have gone
wrong, and did not.

## 4. Task 4 — cross-block recovery on genuinely mixed data

This is the only task that runs through the **public** route: the data are
genuinely mixed, so `bgm()` reaches the mixed sampler by itself and no
`fit_path()` forcing is involved.

**Design.** 3 ordinal (4 categories) + 3 continuous, `n = 1000`, 5 data
seeds `9101…9105`, fits strictly one at a time. Data generated by the
package's own mixed Gibbs sampler,
`bgms:::sample_mixed_mrf_gibbs(..., iter = 500)`, so the planted values are
in the sampler's own parameterization (association scale throughout;
`pairwise_cont = -K/2`).

| block | planted nonzero | planted zero |
|---|---|---|
| ord-ord (3) | `X1-X2 = 0.40`, `X2-X3 = -0.30` | `X1-X3` |
| cont-cont (3) | `Y1-Y2 = 0.175`, `Y2-Y3 = -0.150` (`K12 = -0.35`, `K23 = 0.30`) | `Y1-Y3` |
| cross (9) | `X1-Y1 = 0.30`, `X2-Y2 = -0.25`, `X3-Y3 = 0.20` | the other six |

Two passes per seed: `edge_selection = FALSE` for an unshrunk recovery
estimand, `edge_selection = TRUE` (`bernoulli_prior(0.5)`) for detection.
Ordinal support was 4/4/4 at seeds 9101–9103 and 3/4/4 at seeds 9104–9105
(one variable lost its top category); noted, not corrected — the fit
recodes to the observed levels on both passes.

### 4.1 Per-block summary (5 seeds pooled, 15 fits worth of edges per block family)

| block | recovery slope | intercept | noise sd | detected (PIP > 0.5) | false presence | max PIP on a planted zero |
|---|---|---|---|---|---|---|
| ord-ord   | 1.164 | +0.0011 | 0.0573 | **10/10** | **0/5** | 0.093 |
| cont-cont | 0.993 | −0.0010 | 0.0174 | **10/10** | **0/5** | 0.022 |
| **cross** | **1.039** | **−0.0075** | **0.0423** | **15/15** | **0/30** | **0.455** |

Slopes from the `edge_selection = FALSE` means; detection from the
`edge_selection = TRUE` pass. With selection on, the pooled slopes are
1.171 / 0.981 / 1.018 — i.e. the spike-and-slab costs essentially nothing
here, because every planted edge sits at PIP = 1.

### 4.2 Per-edge detail (mean over the 5 seeds)

| edge | block | truth | mean est (sel=F) | sd over seeds | PIP range |
|---|---|---|---|---|---|
| X1-X2 | ord-ord   |  0.400 |  0.4691 | 0.0703 | 1.000 – 1.000 |
| X1-X3 | ord-ord   |  0     | −0.0047 | 0.0730 | 0.051 – 0.093 |
| X2-X3 | ord-ord   | −0.300 | −0.3448 | 0.0339 | 1.000 – 1.000 |
| Y1-Y2 | cont-cont |  0.175 |  0.1714 | 0.0160 | 1.000 – 1.000 |
| Y1-Y3 | cont-cont |  0     |  0.0018 | 0.0136 | 0.015 – 0.022 |
| Y2-Y3 | cont-cont | −0.150 | −0.1515 | 0.0245 | 1.000 – 1.000 |
| X1-Y1 | cross     |  0.300 |  0.2960 | 0.0288 | 1.000 – 1.000 |
| X1-Y2 | cross     |  0     |  0.0009 | 0.0694 | 0.043 – 0.455 |
| X1-Y3 | cross     |  0     | −0.0200 | 0.0786 | 0.049 – 0.247 |
| X2-Y1 | cross     |  0     |  0.0084 | 0.0535 | 0.044 – 0.200 |
| X2-Y2 | cross     | −0.250 | −0.2761 | 0.0224 | 1.000 – 1.000 |
| X2-Y3 | cross     |  0     | −0.0129 | 0.0396 | 0.032 – 0.100 |
| X3-Y1 | cross     |  0     | −0.0058 | 0.0170 | 0.024 – 0.049 |
| X3-Y2 | cross     |  0     |  0.0000 | 0.0258 | 0.023 – 0.119 |
| X3-Y3 | cross     |  0.200 |  0.2020 | 0.0335 | 1.000 – 1.000 |

Per-seed cross-block slopes: 1.043 / 1.061 / 1.022 / 1.025 / 1.046 — stable
across seeds, residual sd 0.018 – 0.052.

### 4.3 Reading the cross block against the anchors

The cross block is the result; the pure blocks tell us how to read it.

* **Slope.** 1.039, sitting between the continuous anchor (0.993) and the
  ordinal anchor (1.164) — which is what a block spanning the two ought to
  do. Nothing like a factor.
* **Intercept.** −0.0075, small against a noise sd of 0.042.
* **Noise.** sd 0.042, again between the two anchors (0.017 / 0.057) and
  closer to the noisier ordinal side, as expected of a parameter with one
  discrete foot.
* **Detection.** All 15 planted cross edges (3 edges × 5 seeds) at PIP =
  1.000; none of the 30 planted cross zeros exceeds 0.5. The single closest
  call is `X1-Y2` at 0.455 on seed 9105 — no false presence, but the
  smallest margin anywhere in the study, and worth knowing that the cross
  block is where the margin is thinnest.
* **The ordinal anchor's slope is not a mixed-path artefact.** The ord-ord
  slope of 1.164 is nominally ~3 se above 1 (per-seed slopes 1.055 / 1.165 /
  1.145 / 1.103 / 1.351, sd 0.111). But task 2 pins the mixed discrete block
  to the OMRF discrete block at 0.4%, so this belongs to the ordinal MRF
  estimator, not to the mixed path. See open question Q1.

**Verdict.** The cross block recovers its planted truth on the same terms as
the two blocks that do have siblings: slope within 4% of 1, no offset,
perfect detection, no false presence. This is the certification the cross
block could not otherwise get.

---

## 4b. Task 5 (optional) — mgm overlap: RUN, detection only

The estimand map is **not** clean in magnitude, and this is stated rather
than worked around: `mgm` 1.2.15 is nodewise-regularized (glmnet, EBIC,
`lambdaGam = 0.25`, AND rule) and collapses the `(m−1)` parameters of a
categorical node into a single non-negative edge weight, whereas bgms fits a
joint ordinal MRF with one signed association per pair. Weights and
associations are therefore incomparable. The **edge set** is a shared
estimand, so the check was run as detection agreement and nothing else.

One dataset: task 4 seed 9101 (3 ordinal + 3 continuous, `n = 1000`),
ordinal columns shifted to 1..4 for `type = "c"`, continuous as `"g"`.
bgms edge = PIP > 0.5; mgm edge = nonzero `wadj`.

| edge | block | truth | bgms PIP | bgms | mgm w | mgm |
|---|---|---|---|---|---|---|
| X1-X2 | ord-ord | 0.400 | 1.000 | ✓ | 0.0000 | ✗ |
| X1-X3 | ord-ord | 0 | 0.085 | ✗ | 0.0000 | ✗ |
| X2-X3 | ord-ord | −0.300 | 1.000 | ✓ | 2.1579 | ✓ |
| Y1-Y2 | cont-cont | 0.175 | 1.000 | ✓ | 0.3090 | ✓ |
| Y1-Y3 | cont-cont | 0 | 0.022 | ✗ | 0.0000 | ✗ |
| Y2-Y3 | cont-cont | −0.150 | 1.000 | ✓ | 0.2750 | ✓ |
| X1-Y1 | cross | 0.300 | 1.000 | ✓ | 0.1030 | ✓ |
| X1-Y2 | cross | 0 | 0.062 | ✗ | 0.0000 | ✗ |
| X1-Y3 | cross | 0 | 0.247 | ✗ | 0.0000 | ✗ |
| X2-Y1 | cross | 0 | 0.060 | ✗ | 0.0000 | ✗ |
| X2-Y2 | cross | −0.250 | 1.000 | ✓ | 0.9075 | ✓ |
| X2-Y3 | cross | 0 | 0.045 | ✗ | 0.0000 | ✗ |
| X3-Y1 | cross | 0 | 0.027 | ✗ | 0.0000 | ✗ |
| X3-Y2 | cross | 0 | 0.023 | ✗ | 0.0000 | ✗ |
| X3-Y3 | cross | 0.200 | 1.000 | ✓ | 0.5927 | ✓ |

| block | concordance | bgms TP / FP | mgm TP / FP |
|---|---|---|---|
| ord-ord   | 2/3 | 2/2 · 0/1 | 1/2 · 0/1 |
| cont-cont | 3/3 | 2/2 · 0/1 | 2/2 · 0/1 |
| **cross** | **9/9** | **3/3 · 0/6** | **3/3 · 0/6** |
| total     | 14/15 | 7/7 · 0/8 | 6/7 · 0/8 |

**Read.** On the cross block — the one with no sibling and therefore the
point of this whole batch — the two packages agree on all nine edges and
both recover the planted structure exactly. The single disagreement is
`X1-X2`, a planted ordinal edge at 0.40 that bgms detects at PIP = 1.000 and
mgm's EBIC-tuned LASSO drops to exactly zero; with one dataset and a
regularization path involved this is a miss by mgm, not evidence about bgms,
and no attempt was made to tune `lambdaGam` until it agreed.

Nothing was down-sized or skipped in this task; it is one dataset by design,
because that is all a detection-only comparison supports.

## 5. Findings

### 5.1 F1 — the mixed model is undefined at a degenerate block: hangs on `q = 0` with selection, errors otherwise (first-class; not user-reachable)

**Severity: moderate.** Not reachable through any exported entry point
(§1), so no user can trigger it on 0.2.0.0 as shipped. First-class because
it is exactly the cross-path class this review hunts — the mixed machinery
does not reduce to its siblings, it refuses to run there — and because one
of the three failure modes is a **silent non-terminating spin**, the worst
failure shape a sampler can have. It is also what blocked tasks 2 and 3 as
briefed.

**Reproduction.** `probe.R`, one configuration per invocation, `p = 4`,
`n = 200`, `iter = warmup = 100`, `chains = 1`:

| configuration | `edge_selection` | result | wall |
|---|---|---|---|
| mixed, all-ordinal (`q = 0`) | FALSE | error `Col::subvec(): indices out of bounds or incorrectly used` | 0.0 s |
| mixed, all-ordinal (`q = 0`) | TRUE  | **no return** | killed at 150 s |
| mixed, all-continuous (`p = 0`) | FALSE | error `subscript out of bounds` | 0.0 s |
| mixed, all-continuous (`p = 0`) | TRUE  | error `subscript out of bounds` | 0.0 s |
| mixed, 3 ordinal + 1 continuous | TRUE | OK | 0.3 s |
| mixed, 1 ordinal + 3 continuous | TRUE | OK | 0.1 s |

The last two rows matter: every configuration a user can actually reach
(`p ≥ 1`, `q ≥ 1`) is fine.

**Mechanism 1 — the hang (`q = 0`, selection on).** A `sample(1)` of the
stopped process put 2387 of 2387 stack samples in
`MixedMRFModel::get_vectorized_indicator_parameters()` at
[mixed_mrf_model.cpp:858-859](src/models/mixed/mixed_mrf_model.cpp#L858-L859):

```cpp
// 2. Upper-triangle of Gyy
for(size_t i = 0; i < q_ - 1; ++i) {
    for(size_t j = i + 1; j < q_; ++j) {
        out(idx++) = gyy(i, j);
```

`q_` is `size_t` ([mixed_mrf_model.h:355](src/models/mixed/mixed_mrf_model.h#L355)),
so at `q_ = 0` the bound `q_ - 1` wraps to `SIZE_MAX` and the outer loop is
asked for ~1.8 × 10¹⁹ iterations. The inner loop `j < q_` is empty for every
`i`, so **nothing is written out of bounds** — it is a pure no-op spin: no
crash, no memory growth, no progress, 100% of one core indefinitely. At
~10⁹ iterations/s it would finish in roughly 600 years. The same wrap sits
on the `Gxx` loop one block above
([:851](src/models/mixed/mixed_mrf_model.cpp#L851)) for `p_ = 0`.

Reached from
[chain_runner.cpp:150](src/mcmc/execution/chain_runner.cpp#L150) and
[chain_runner.cpp:287](src/mcmc/execution/chain_runner.cpp#L287), both
guarded by `config.edge_selection` — which is why turning selection off
changes the failure mode instead of curing it.

**The fix pattern already exists in the same file**, which is what makes
this an oversight rather than a design decision:

* `get_vectorized_rb_alpha()`
  ([:875](src/models/mixed/mixed_mrf_model.cpp#L875), loops at
  [:882-884](src/models/mixed/mixed_mrf_model.cpp#L882-L884)) and
  `get_vectorized_rb_pregamma()`
  ([:892](src/models/mixed/mixed_mrf_model.cpp#L892), loops at
  [:896-899](src/models/mixed/mixed_mrf_model.cpp#L896-L899)) — whose own
  comment says they "mirror `get_vectorized_indicator_parameters()`
  exactly" — use the underflow-safe form
  `for(size_t i = 0; i + 1 < q_; ++i)`.
* The two adaptive-Metropolis sweeps at
  [:1129](src/models/mixed/mixed_mrf_model.cpp#L1129) and
  [:1201](src/models/mixed/mixed_mrf_model.cpp#L1201) wrap their `q_ - 1`
  loop in `if(q_ >= 2) { … }`.

So of the four places that traverse the `Gyy` upper triangle, three are
protected and one is not. `mixed_mrf_gradient.cpp` has the same unguarded
shape at [:570](src/models/mixed/mixed_mrf_gradient.cpp#L570) (`q_ - 1`) and
[:31, :87, :146, :498](src/models/mixed/mixed_mrf_gradient.cpp#L31)
(`p_ - 1`).

**Mechanism 2 — the error (`q = 0`, selection off).** With selection off the
indicator vectorizer is never called and the run reaches the parameter
vectorizer instead, where the zero-width Cholesky block gives
`subvec(idx, idx + 0 - 1)`:

* [mixed_mrf_model.cpp:630](src/models/mixed/mixed_mrf_model.cpp#L630) —
  `out.subvec(idx, idx + chol_constraint_structure_.active_dim - 1) = theta_yy_;`
* [mixed_mrf_model.cpp:775](src/models/mixed/mixed_mrf_model.cpp#L775) —
  `arma::vec theta_yy = params.subvec(idx, idx + chol_dim - 1);`

R traceback: `fit_path → run_sampler → run_sampler_mixed_mrf →
sample_mixed_mrf` (`trace.R q0`).

**Mechanism 3 — the error (`p = 0`).** This one never reaches C++. R
traceback (`trace.R p0`):

```
Error in x[, node] : subscript out of bounds
build_spec_mixed_mrf -> reformat_ordinal_data -> sort -> unique
```

`build_spec_mixed_mrf()` splits the data at
[build_spec.R:236-243](R/build_spec.R#L236-L243) and hands the zero-column
discrete part to `reformat_ordinal_data()` at
[build_spec.R:310-314](R/build_spec.R#L310-L314), which indexes a column
that is not there.

**What it costs this review.** The sharp all-ordinal reduction — the single
strongest cross-path test available for the mixed sampler, and the reason
the OMRF sibling exists as a reference — cannot be run at all without
touching source, which this batch may not do. Tasks 2 and 3 fall back to
substitutes (§2.1, §3.0) that carry a measurable confound.

**Suggested resolution** (not applied; report-only batch). Either make the
three sites total — `i + 1 < q_`, matching the sibling functions — so the
degenerate blocks work and the sharp reduction becomes testable, or reject
`p == 0 || q == 0` at the top of `build_spec_mixed_mrf()` with a clear
message so the machinery documents its own domain. The first is
substantially more valuable: it turns an untestable relation into a
testable one. Either way the loop bounds should be fixed, because a
non-terminating spin is not an acceptable failure mode even on an
unreachable path.

### 5.2 F2 — `bgm()`'s `cores` default oversubscribes; unrelated to this batch but seen throughout (informational)

**Severity: informational.** `bgm(cores = parallel::detectCores())`
([bgm.R:495](R/bgm.R#L495)) with `chains = 4` requests every hardware
thread on the machine for four chains. Every fit in this report passed
`cores = 4` explicitly to stay inside the batch's machine budget. Noted
only because it is visible from this seat; it is not a finding about the
mixed model and is presumably already known to the review.

### 5.3 Non-findings — things checked that came back clean

* **No `omega`-vs-`2·omega` convention divergence anywhere.** Task 3's
  slope of mixed on GGM is 1.0024 / 1.0004 and shrinks with `n`; task 2's
  slope of mixed on OMRF is 1.0036 with an intercept of +0.0003. The
  `extract_precision()` round-trip lands both continuous paths on the same
  scale (§3.3).
* **Prior plumbing is shared, not parallel.** Both discrete paths read
  `threshold_prior_type` from the same field with the same `"beta-prime"`
  default and hand it to the same `create_parameter_prior()`
  ([sample_omrf.cpp:67-68](src/sample_omrf.cpp#L67-L68),
  [sample_mixed.cpp:83-84](src/sample_mixed.cpp#L83-L84)).
* **`delta` resolves consistently.** `0.5·log(num_variables)` for GGM and
  `0.5·log(num_continuous)` for mixed
  ([bgm_spec.R:443-451](R/bgm_spec.R#L443-L451)) coincide whenever the two
  are being compared on the same continuous block.
* **Selection decisions match across paths.** All 15 task-2 edges get the
  same call from OMRF and from the mixed sampler, with `max |ΔPIP| =
  0.0096`.
* **No degenerate-block risk in reachable territory.** `p ≥ 1, q ≥ 1` is
  guaranteed by the routing (§1), and the two minimal mixed
  configurations (3+1 and 1+3) fit without incident.

---

## 6. Open questions

**Q1 — is the ordinal block's ~15% upward slope at `n = 1000, p = 3` real,
and if so what is it?** Task 4's ord-ord recovery slope is 1.164 (per-seed
1.055 / 1.165 / 1.145 / 1.103 / 1.351, sd 0.111, so ≈3 se above 1); both
planted ordinal edges come back inflated in the same direction
(`0.400 → 0.469`, `−0.300 → −0.345`), while the continuous block in the same
fits is unbiased (0.993). Three things are already known about it:

* it is **not** a mixed-path artefact — task 2 pins the mixed discrete block
  to the OMRF discrete block at 0.4%, so both paths would carry it equally;
* it is **not** visible on the task-2 ordinal design (`p = 6`, `n = 1000`,
  OMRF): recovery slope 0.9876, though with wide per-edge scatter (per-edge
  ratios 0.76 – 1.18, residual sd 0.037);
* it runs **against** the prior, which shrinks toward zero, so it is not
  slab shrinkage.

Candidates not separated here: finite-sample bias of the ordinal
pseudolikelihood at `p = 3`; an interaction with the cross-block coupling
(the ordinal variables in task 4 are each tied to a continuous one);
or simply the slope being carried by two distinct truth values. Resolving it
needs an ordinal-only sweep across `n` and `p`, which is outside this
batch's remit — it is a property of the OMRF estimator, not of the
mixed/OMRF relation this batch certifies.

**Q2 — should the degenerate block be supported or rejected?** F1 offers a
choice (§5.1). Supporting it has a concrete payoff beyond tidiness: it makes
the sharp all-ordinal reduction runnable, and that is the strongest
cross-path test the mixed sampler can be given. If it is instead rejected at
the R layer, the sharp reduction stays permanently untestable and this
batch's tasks 2–3 remain substitutes forever.

**Q3 — the thinnest margin is in the cross block.** `X1-Y2` reaches
PIP = 0.455 on a planted zero at seed 9105 (§4.2), against a maximum of
0.093 anywhere in the pure blocks. Five seeds is not enough to say whether
the cross block's null distribution has a heavier tail than the pure blocks'
or whether this is one draw. A larger null sweep on cross edges would settle
it.

---

## 7. Verification gate

The brief's gate is internal to this report (no suite runs — no code
changed).

1. **Every comparison carries its yardstick.** Task 2: §2.2 table (0.00041
   pairwise, 0.00505 thresholds), §2.3 (0.00289 PIP) — all from the same-path
   refit at seed 4002. Task 3: §3.1, both tables, yardstick column at each
   `n` from the GGM refit at seed 5002. Task 4 has no cross-path claim to
   judge, so its yardstick is the between-seed sd, reported per edge in
   §4.2 and per block in §4.1.
2. **Seeds, dimensions and prior mappings stated for every fit.** Priors:
   §0 table, one mapping for all paths. Dimensions and seeds: §2.1 (p = 6,
   4 categories, n = 1000, data 20230801, companion 20230802, fits
   4001/4002; confound draw 2 from data seed 20230803), §3.0 (Q = 6,
   n ∈ {500, 5000}, data 20230802, fits 5001/5002), §4 (3+3, n = 1000, seeds
   9101–9105), §4b (task-4 seed 9101, mgm seeded at 9101).
3. **Task 1's routing answer with file:line citations.** §1, ten citations,
   plus the `NAMESPACE` negative.
4. **Anything skipped or down-sized is named with its reason.**
   * Tasks 2 and 3 **as briefed** were not run — blocked by F1 (§2.0, §3.0).
     Both were replaced by named substitutes whose weakness is measured
     rather than assumed (§2.4).
   * Task 3 uses `edge_selection = FALSE` throughout, a deliberate deviation
     recorded with its reason in §3.0.
   * Task 5 was **run**, not skipped, but restricted to edge detection; the
     magnitude comparison is not available and §4b says why.
   * Task 4 seeds 9104 and 9105 have one ordinal variable at 3 of 4
     categories (§4); not corrected, not excluded.
   * Nothing else was reduced. All 5 task-4 seeds ran both passes; no
     top-N truncation or sampling anywhere.

**Deliverable verdicts, restated compactly:**

| task | verdict |
|---|---|
| 1 routing | mixed path **not** user-reachable on pure data → tasks 2–3 are internal-machinery certification |
| 2 mixed vs OMRF | sharp check **not runnable** (F1); substitute shows **no systematic divergence** — slope 1.0036, PIP agreement ≤0.0096, residual gap demonstrated to be the companion |
| 3 mixed vs GGM | slope 1.00239 at `n = 500`, 1.00037 at `n = 5000` — **shrinking**, so finite-sample pseudo-vs-exact, **not** a convention bug |
| 4 cross-block recovery | slope 1.039, noise sd 0.042, **15/15** detected, **0/30** false presences; bracketed by the pure-block anchors |
| 5 mgm overlap | run, detection only: **9/9** agreement on the cross block, 14/15 overall |

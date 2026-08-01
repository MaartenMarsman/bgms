# Report 04 — Statistical certification of rc1

Slow-tier suites, zratio route certificates, sensitivity reproduction, and the
hierarchical prior-chain cost. **Report only — nothing was changed.**

Agent: brief 04. Date: 2026-08-01.

---

## 1. What was done

### Environment

| | |
|---|---|
| Frozen target | tag `v0.2.0.0-rc1` (= `4f5bddab`, merged to `main` 2026-08-01 12:52) |
| Source | `~/bgms-review/bgms-rc1/` (clean export from brief 01) |
| Library | `~/bgms-review/lib/bgms`, `packageVersion("bgms")` = **0.2.0.0** — verified |
| R | 4.6.0 (2026-04-24), aarch64-apple-darwin23 |
| Machine | Darwin 25.4.0, macOS 26.4.1, 15 cores, Apple clang 21.0.0 |
| Working dir | `~/bgms-review/certification/` |

Every R invocation in this report used `R_LIBS=~/bgms-review/lib` or an explicit
`.libPaths()` prepend, so the rc1 install resolved first. Verified once
explicitly:

```
/Users/maartenmarsman/bgms-review/lib
/Library/Frameworks/R.framework/Versions/4.6/Resources/library
ver: 0.2.0.0
path: /Users/maartenmarsman/bgms-review/lib/bgms
```

### Commands

```sh
# Task A
cd ~/bgms-review
NOT_CRAN=true BGMS_RUN_SLOW_TESTS=true Rscript run-slow-suite.R slow
NOT_CRAN=true BGMS_RUN_SLOW_TESTS=true Rscript certification/rerun-failing.R

# Task B (from ~/bgms-review/certification/, which holds a copy of dev/validation/)
R_LIBS=~/bgms-review/lib Rscript dev/validation/zratio_deployed_route_cert.R --cores=6
R_LIBS=~/bgms-review/lib Rscript dev/validation/zratio_isolated_route_cert.R --cores=12
R_LIBS=~/bgms-review/lib Rscript dev/validation/zratio_anchor_gate.R --cores=12
R_LIBS=~/bgms-review/lib Rscript dev/validation/zratio_shape_verdict.R --shape=10 --cores=12
R_LIBS=~/bgms-review/lib Rscript dev/validation/zratio_collapse_reachability.R

# Tasks C, D
Rscript certification/task-c-sensitivity.R ; Rscript certification/task-c-analyse.R
Rscript certification/task-d-hier-cost.R   ; Rscript certification/task-d2-chain-trigger.R
```

`run-slow-suite.R` mirrors brief 01's `run-suite.R` exactly; only the env gate differs.

### Protection of the gold bank

The validation scripts write their results back into `dev/validation/`. To keep
the banked references pristine I copied `dev/validation/` into
`~/bgms-review/certification/dev/validation/` and ran everything with
`~/bgms-review/certification/` as the working directory, so all writes landed on
the copy. Pristine copies of all 13 banked `.rds` were kept in
`~/bgms-review/certification/banked/` for comparison.

Verified after the run:

- all 13 banked `.rds` in `~/bgms-review/bgms-rc1/dev/validation/` byte-identical
  to the pre-run backup;
- `git status --porcelain dev/validation/` on the Dropbox checkout: **empty**.

The Dropbox repo was read only — never built in, never written to.

### Adaptations made to the banked scripts (complete list)

| Script | Adaptation | Why |
|---|---|---|
| `zratio_shape_verdict.R` | `suppressMessages(devtools::load_all(".", quiet = TRUE))` → `suppressMessages(library(bgms))` | Script assumes a source tree; the brief requires running against the rc1 **install**. Same source, compiled at -O2. |
| `zratio_collapse_reachability.R` | same substitution | same |
| `zratio_deployed_route_cert.R` | ran with `--cores=6` instead of the documented `--cores=12` | Ran alongside another job. Affects surface-build scheduling only; see finding 04-7 for the measured impact (3.9e-8). |
| `zratio_anchor_gate.R` | ran with default `--cost_shapes=1,10`; the bank holds four shapes (1, 2, 5, 10) | The script's own default. Cost arm coverage is narrower than the bank — noted, not a deviation. |
| all | run from `~/bgms-review/certification/` rather than the repo root | protects the bank (above) |

No script was rewritten or replaced. All five requested scripts were runnable as
banked.

The gold bank itself (`zratio_gold_bank.rds`, 470 cells) was pure lookup: **0 new
cells computed, 0 banked values changed**. Every gold figure below is the banked
one, so all deviations reported are on the surface/deployed side.

---

## 2. Findings

Severity per the brief: a suite failure or certificate deviation is at least `major`.

---

### 04-1 · `major` · Three slow-tier assertions encode a fence that PR #193 replaced — and the tier that would have caught it does not run on that branch

Three of the four slow-tier failures share **one** root cause, which I traced and
then confirmed causally.

PR #193 (`cb8e81c1`, "de-gate, and deploy the surface across its validated shape
range") set the surface-deployment fence to:

```r
# R/zratio_surfaces.R:225-226
.zratio_surface_shape_lo = 0.5
.zratio_surface_shape_hi = 10
```

Measured against the rc1 install:

```
shape 0.5   build -> surface        shape 3     build -> surface
shape 1     build -> surface        shape 5     build -> surface
shape 2     build -> surface        shape 10    build -> surface
shape 2.5   build -> surface        shape 12    build -> NULL (fenced)
```

The same commit **wrote** this assertion, which asserts the opposite for two
shapes that are inside its own new fence:

```r
# tests/testthat/test-zratio-surface-build.R:149-157
for(shape in c(2.5, 5)) {
  zc = bgms:::zratio_constants(0.5 * log(12), 2, alpha = shape)
  expect_null(bgms:::zratio_build_surfaces(zc, max_size = 8L, cores = 1L))
}
```

It is contradictory as authored. It shipped because it sits behind
`BGMS_RUN_SLOW_TESTS` and never executed.

The third failure, `test-zratio-gauge.R:213`, states the same stale premise in
its own comment — *"the non-unit Gamma-diagonal fence (shape = 2), where the
surface is not deployed and the engine falls back to the additive-counts
saddle"* — which the shipped fence falsifies. I confirmed this is the cause
rather than inferring it, by re-running the test's exact fixture under both
fences:

```
shipped fence [0.5, 10]            amplification  17.19 | harm_pred 0.00794 | threshold 0.01 | harm_flag FALSE
fence narrowed to alpha == 1 only  amplification   5.20 | harm_pred 0.03937 | threshold 0.01 | harm_flag TRUE
```

So: shape 2 now deploys the surface, the more accurate kernel cuts the projected
distortion 5-fold, and the flag correctly stays down. **The harm channel is not
broken — the fixture is stale.** It still fires under the kernel it was written
for.

**Why this was not caught.** `nightly-validation.yaml` is the only job that sets
`BGMS_RUN_SLOW_TESTS`, and it is `schedule:` (Mon/Thu 03:00 UTC) plus
`workflow_dispatch` only — **no `push`, no `pull_request` trigger**. It therefore
runs on `main`. PRs #193 and #194 landed on `main` only today (2026-08-01 12:52,
in the release-evaluation merge); the last nightly ran 2026-07-30 and was green
because it never saw this code. The next scheduled run (Mon 2026-08-03) is the
first that would go red.

This is finding F-044 materialising concretely: a core statistical gate asserting
the opposite of shipped, intended behaviour, merged and released without ever
executing.

**Not a product defect.** In all three cases the shipped code behaves as designed
and documented; the assertions are stale. Fixing them is test work, not sampler
work. But rc1 currently ships with a red slow tier.

---

### 04-2 · `major` · `test-prior-sensitivity.R:432` fails on a self-calibrating tolerance derived from a 2-sample noise estimate

Separate root cause from 04-1.

```
── 1. Failure ('test-prior-sensitivity.R:432:3'): the difference-scale reweighting reproduces a refit at that scale
Expected `max(abs(rw$pip[2, pairwise] - pip_of(f2)[pairwise]))` < `4 * noise`.
Actual comparison: 0.0255 >= 0.0202
Difference: 0.0053 >= 0
```

Overshoot is 26% — the observed gap is 5.05× the noise estimate where the gate
allows 4×. Deterministic: reproduces to the same digits on rerun (fixed seeds
11/12/13).

The weak point is the threshold's construction. `noise` is
`max|pip_of(f2) - pip_of(f2b)|` — a **two-refit** estimate of run-to-run spread,
used as the scale for a max-over-edges statistic. A 2-sample spread estimate is
itself extremely noisy, so the gate's own threshold has large sampling variance;
`4 * noise` here evaluated to 0.0202, and a modestly luckier seed pair would have
passed. This is fragile by construction independent of whether the underlying
reweighting is sound.

This is a `bgmCompare` difference-selection test — a pseudolikelihood path. Under
MM's scope caveat that is exactly where the recovery/cross-validation suites are
the correct gate rather than SBC, and this test *is* such a cross-validation
(reweighting vs. an independent refit). So the failure is in scope and meaningful;
what I cannot say from one run is whether the reweighting is genuinely off by
~0.005 pip or whether the tolerance is simply too tight. See Open questions.

---

### 04-3 · `note` (positive) · No SBC test is miscalibrated against a pseudolikelihood path

The brief asked me to flag any test running SBC-style checks against a non-GGM
path. I checked all five files mentioning SBC. **Nothing to flag** — the suite
already encodes MM's caveat correctly:

- `test-sbc-ggm.R` — all six SBC tests use `variable_type = "continuous"`
  exclusively. Entirely within the valid-certification scope.
- `test-sbc-correction.R` — the only non-GGM cells are **prior-only** identities.
  The file header states the reasoning explicitly:

  > *"The with-data cells are restricted to the all-continuous models because
  > those fit the exact Gaussian likelihood; the mixed sampler fits a
  > pseudolikelihood, whose pseudo-posterior is not calibrated in the exact-SBC
  > sense (rank uniformity fails by construction, and a sandwich/Godambe
  > recalibration is out of scope for a package test). With no data the
  > pseudolikelihood factor drops out, so the prior-only chain is exact and gates
  > the mixed correction sharply."*

- `test-ggm-nuts.R`, `test-regressions.R`, `test-sample-ggm-prior.R` — mention SBC
  only in prose/TODO comments; no SBC assertions.

No pseudolikelihood-path behaviour was read as an SBC failure anywhere in this
report.

---

### 04-4 · `minor` · The entire slow tier is unmigrated to the 0.2.0 argument names

All **1152** warnings raised by the slow tier are deprecation notices — zero
substantive warnings.

| File | Warnings | Message |
|---|---:|---|
| `test-sbc-ggm.R` | 1000 | `pairwise_scale` deprecated → use `interaction_prior` |
| `test-parameter-recovery-ggm.R` | 100 | `pairwise_scale` deprecated |
| `test-mixed-nuts.R` | 38 | `pairwise_scale` (19) + `main_alpha` (19) |
| `test-scaling-diagnostics.R` | 14 | `pairwise_scale` (9) + `main_alpha` (5) |

The fast tier emits **0** warnings, so this is confined to the dormant tier —
consistent with 04-1's diagnosis that this code has not been exercised since the
0.2.0 deprecations landed. Cosmetic, but it means 1152 lines of noise will bury
any real warning the next time the nightly runs.

---

### 04-5 · `minor` · `plot()` and `print()`/`$edges` report different Bayes factors for the same scale point

Found while doing Task C. At the chosen scale (1x — the original fit, the same
quantity by both routes):

- `$log10_bf[chosen_index, ]` (what `plot()` draws) clamps pip at `1 - 1e-6`,
  capping log10 BF at **6.0**;
- `$edges$chosen_scale_log10_bf` (what `print()` and the table report) is the
  uncapped Rao-Blackwellised value.

Measured on the Task C fit (136 edges):

| | |
|---|---:|
| edges where the two differ by > 0.01 | **93** |
| edges where they differ by > 0.5 | **25** |
| max abs difference | **8.26 log10 units** |
| curve points pinned at the 6.0 cap at 1x | 20 |
| …of those, `$edges$saturated == FALSE` | 11 |

Examples:

```
            edge curve_1x table_1x      pip saturated
 intrusion-flash        6  9.27280 1.000000     FALSE
 intrusion-anger        6  4.60857 0.999975     FALSE
     sleep-anger        6 14.25562 1.000000     FALSE
   sleep-startle        6 13.53463 1.000000     FALSE
```

Two mitigating facts, both verified: **no verdict changes** (0/136
disagreements between the two routes, at any threshold), and the cap *is*
documented in `?plot.bgms_prior_sensitivity`. So this is presentation, not
inference.

Still worth a line in the docs, because it has a concrete reading consequence:
the flat-at-6.0 plateaus in MM's saved figure are clamp artifacts, not evidence
plateaus, and `$edges$saturated` does not mark them (`saturated` means
`zeroflip` — the indicator never flipped — which is a different concept and is
correctly `FALSE` for these edges).

---

### 04-6 · `minor` · F-036 is wrong as stated: the default hierarchical fit does not trigger the prior-only chain, and the real trigger costs ~1.4 s

Detail and numbers in §5 (Task D). Summary: the 5000-iteration chain is real
(`iter = 4000L` + `warmup = 1000L`) but is **not** reached by
`bgm(..., precision_graph_prior = "hierarchical")`. The trigger is the **edge
prior**, not the precision-graph prior. Brief 03 should not describe this as a
hierarchical-fit cost.

---

### 04-7 · `note` · The anchor gate reproduces exactly — including a standing disagreement between the gate and the shipped multiplier

The certificate **passes** (bit-identical reproduction of the bank, §4). But
reproducing it re-exposes a pre-existing disagreement that is in the bank too, so
it is not rc1 drift:

In **5 of 20** cells the gate's re-derived multiplier exceeds the shipped
`zratio_anchor_shape_multiplier()`, i.e. the shipped anchor budget is **lower**
than the gate says parity needs:

| eta | family | alpha | gate resolves | shipped |
|---:|---|---:|---:|---:|
| 1 | cn | 3 | 4× | 2× |
| 2 | cn | 3 | 4× | 2× |
| 2 | cn | 5 | 4× | 2× |
| 2 | cn | 10 | 4× | 2× |
| 2 | bip | 10 | 4× | 2× |

And in **1** cell (eta 2, cn, shape 0.5) parity is *not met even at the 4× hard
cap* (`resolved = NA`), where shipped is already 4×.

Both states are identical in the banked run, so this is a known standing
condition, not a regression. Flagging it because it is a live statistical
question the release arguably ought to answer. Needs MM — see Open questions.

---

## 3. Evidence — Task A: the dormant suites against rc1

### Totals

| | fast tier (brief 01) | **slow tier (this run)** |
|---|---:|---:|
| files | 77 | 77 |
| tests | 1129 | 1129 |
| **pass** | 7954 | **8650** |
| **fail** | 0 | **4** |
| error | 0 | **0** |
| warning | 0 | 1152 |
| skip | 103 | **7** |
| wall clock | 149.7 s | **710.4 s** (`real 711.48`) |

Unlocking the tier adds **696 assertions** and **4.7×** runtime. Skips drop
103 → 7; the 7 remaining are 2 documented `test-regressions.R` skips plus 5
others unrelated to the slow gate. **All 23 slow-gated files ran with 0 skips.**

### Per-file, the 23 slow-gated files

Sorted by runtime. `fast_skip` = how many of that file's tests were skipped in
brief 01's fast run, i.e. what the gate was hiding.

| file | pass | fail | err | warn | skip | secs | fast_skip |
|---|---:|---:|---:|---:|---:|---:|---:|
| `test-validation-slow.R` | 3 | 0 | 0 | 0 | 0 | 122.7 | 3 |
| `test-prior-sensitivity.R` | 83 | **1** | 0 | 0 | 0 | 118.9 | 5 |
| `test-mixed-nuts.R` | 264 | 0 | 0 | 38 | 0 | 50.0 | 11 |
| `test-hier-zratio-identity.R` | 39 | 0 | 0 | 0 | 0 | 44.9 | 8 |
| `test-sbc-correction.R` | 16 | 0 | 0 | 0 | 0 | 38.8 | 6 |
| `test-scaling-diagnostics.R` | 45 | 0 | 0 | 14 | 0 | 33.7 | 9 |
| `test-sbc-ggm.R` | 12 | 0 | 0 | 1000 | 0 | 32.3 | 6 |
| `test-zratio-extrapolation-notice.R` | 14 | 0 | 0 | 0 | 0 | 31.0 | 1 |
| `test-zratio-pivot-slice.R` | 19 | 0 | 0 | 0 | 0 | 27.3 | 4 |
| `test-rb-inclusion-probabilities.R` | 68 | 0 | 0 | 0 | 0 | 24.7 | 2 |
| `test-zratio-surface-build.R` | 25 | **2** | 0 | 0 | 0 | 24.2 | 3 |
| `test-bgm-hier-spec.R` | 47 | 0 | 0 | 0 | 0 | 14.8 | 4 |
| `test-mixed-correction.R` | 38 | 0 | 0 | 0 | 0 | 11.6 | 4 |
| `test-zratio-law.R` | 49 | 0 | 0 | 0 | 0 | 10.8 | 4 |
| `test-zratio-gauge.R` | 50 | **1** | 0 | 0 | 0 | 9.2 | 1 |
| `test-zratio-gamma-shape.R` | 40 | 0 | 0 | 0 | 0 | 8.3 | 5 |
| `test-sbm-correction.R` | 39 | 0 | 0 | 0 | 0 | 7.3 | 2 |
| `test-parameter-recovery-ggm.R` | 32 | 0 | 0 | 100 | 0 | 6.0 | 2 |
| `test-bb-correction.R` | 19 | 0 | 0 | 0 | 0 | 4.8 | 2 |
| `test-zratio-surface-extension.R` | 12 | 0 | 0 | 0 | 0 | 2.1 | 1 |
| `test-zratio-cauchy.R` | 16 | 0 | 0 | 0 | 0 | 1.3 | 5 |
| `test-ggm-nuts.R` | 86 | 0 | 0 | 0 | 0 | 0.5 | 6 |
| `test-rattle-edge-selection.R` | 13 | 0 | 0 | 0 | 0 | 0.2 | 2 |
| **totals** | **1029** | **4** | **0** | **1152** | **0** | **625.4** | **96** |

All 8 files the brief named ran, none skipped.

Two caveats on the runtime column, both checked:

1. Per-file seconds sum to 625.4 vs 710.4 s wall — the ~85 s difference is work
   outside `test_that()` blocks, which `testthat` does not attribute per test.
2. `test-ggm-nuts.R` at 0.5 s contradicts its own header ("take several
   minutes"), so I verified it rather than assuming. It genuinely does the work:
   86 real assertions, 0 skips, and a standalone timing of one of its fits — a
   2-chain × 5000-iteration NUTS fit at p=4, n=200 — is **0.04 s** on this
   machine. The header's cost estimate is simply stale for current hardware. The
   gate is real; no finding.

### The four failures, verbatim

```
══ Failed ══════════════════════════════════════════════════════════════════════
── 1. Failure ('test-prior-sensitivity.R:432:3'): the difference-scale reweighti
Expected `max(abs(rw$pip[2, pairwise] - pip_of(f2)[pairwise]))` < `4 * noise`.
Actual comparison: 0.0255 >= 0.0202
Difference: 0.0053 >= 0

── 2. Failure ('test-zratio-gauge.R:213:3'): a known-biased evidence-free fit fi
Expected `pc$harm_flag` to be TRUE.
Differences:
`actual`:   FALSE
`expected`: TRUE


── 3. Failure ('test-zratio-surface-build.R:156:5'): the build fences shapes out
Expected `bgms:::zratio_build_surfaces(zc, max_size = 8L, cores = 1L)` to be NULL.
Differences:
`actual` is a list
`expected` is NULL


── 4. Failure ('test-zratio-surface-build.R:156:5'): the build fences shapes out
Expected `bgms:::zratio_build_surfaces(zc, max_size = 8L, cores = 1L)` to be NULL.
Differences:
`actual` is a list
`expected` is NULL
```

Failures 3 and 4 are the same assertion at `shape = 2.5` and `shape = 5`.

### Reproducibility

All three files rerun individually, same env, same seeds
(`certification/rerun-failing.log`):

| file | pass | fail | reproduces? |
|---|---:|---:|---|
| `test-prior-sensitivity.R` | 83 | 1 | **yes — identical numbers** (`0.0255 >= 0.0202`) |
| `test-zratio-gauge.R` | 50 | 1 | **yes — identical** |
| `test-zratio-surface-build.R` | 25 | 2 | **yes — identical** |

**None of the four is stochastic.** All are deterministic properties of rc1.

---

## 4. Evidence — Task B: zratio route certificates

Tolerances are the ones the scripts and the bank encode, not zero. The gold bank
`sd` column bounds what any claim can resolve; each script carries its own
acceptance constant (`1e-12` identity, `2.84e-4` mediation bound, `0.003` nat
in-hull envelope, `60 s` build cost).

| # | script | certifies | verdict |
|---|---|---|---|
| 1 | `zratio_deployed_route_cert.R` | the hot path a fit takes (`log_zratio`) is the same function as the scored surface path (`logR`), and the composition hits gold | **PASS** |
| 2 | `zratio_isolated_route_cert.R` | above shape 10 the deployed value is exactly `log(psi0)`, and the mediation it drops is inside the recorded bound | **PASS** (bit-identical to bank) |
| 3 | `zratio_anchor_gate.R` | anchor Monte-Carlo error at non-unit shape matches the shape-1 reference; cap-80 build stays under 60 s | **PASS** (bit-identical to bank; see 04-7) |
| 4 | `zratio_shape_verdict.R` | the surface tracks block-Gibbs gold inside the 0.003-nat envelope at shape 10, with shape 2 as a banked control | **PASS** |
| 5 | `zratio_collapse_reachability.R` | whether the additive zero-collapse is reachable by a real fit at `delta = 0.5*log(p)` | **PASS** (bit-identical to bank) |

### 1. Deployed-route certificate — PASS

Arm 1 (identity, `log_zratio == logR`), 36 cells across shapes {0.5, 1, 2, 3, 5,
10} × {cn, bip} × sizes {10, 24, 42}:

```
identity worst 0 over 36 cells -> PASS       (limit 1e-12)
```

Exactly zero, not merely under tolerance. Matches the bank (also 0).

Arm 2 (end-to-end gold spot per shape, eta 2, cn, size 42):

| alpha | gold | gold_sd | deployed | err_deployed | err_additive |
|---:|---:|---:|---:|---:|---:|
| 0.5 | 0.3679 | 2.09e-04 | 0.3678 | 5.94e-05 | 0.36788 |
| 1.0 | 0.3431 | 4.44e-04 | 0.3426 | 5.31e-04 | 0.34310 |
| 2.0 | 0.2895 | 2.73e-04 | 0.2883 | 1.24e-03 | 0.17328 |
| 3.0 | 0.2353 | 1.81e-04 | 0.2350 | 3.28e-04 | 0.00558 |
| 5.0 | 0.1469 | 5.58e-05 | 0.1467 | 2.40e-04 | 0.02257 |
| 10.0 | 0.0477 | 9.39e-07 | 0.0477 | 4.77e-06 | 0.04773 |

Against the bank: gold identical (lookup, as designed);
**max |deployed_new − deployed_banked| = 3.9e-8**, four orders of magnitude below
the errors being reported. That residual is the `--cores=6` vs `--cores=12`
surface-build scheduling difference. No deviation.

Runtime 113 s.

### 2. Isolated-route certificate — PASS, bit-identical

```
identity worst 0 over 48 cells -> PASS (exact)
arm 3 route counters: every evaluation isolated, none predicted, none additive: PASS
worst deployed error 0.000218 against the recorded shape-10 bound 0.000284 -> inside
worst additive error  0.0335 (154x the deployed one)
gold sd worst 8.53e-07; deployed error / gold sd worst 4102.9
```

Comparison to bank: `identity` frame `identical()` TRUE; `spots` frame
`all.equal()` TRUE; max |deployed diff| = **0**; max |gold diff| = **0**.

Runtime 25 s.

### 3. Anchor gate — PASS, bit-identical (with the standing 04-7 disagreement)

All 20 gate cells resolve to exactly the banked multiplier, and the `ratio_1x`
values match to displayed precision — expected, since `anchor_spread()` uses
fixed seeds 5001–5008.

```
shipped constant unchanged between bank and now: TRUE
resolved identical bank vs new: TRUE
cells whose RESOLVED multiplier moved since the bank: 0
cells where parity was NOT met even at the 4x cap — banked: 1  this run: 1
```

Cap-80 build cost, serial, installed build:

| alpha | banked | this run | within 60 s |
|---:|---:|---:|---|
| 1 | 24.5 s | **26.0 s** | yes |
| 10 | 43.4 s | **47.5 s** | yes |

Both inside the standing ≤ 60 s acceptance; +6% and +9% vs the bank, consistent
with different hardware. Coverage note: the bank holds four cost shapes
(1, 2, 5, 10), I ran the script default (1, 10).

Runtime 124 s.

### 4. Shape verdict at shape 10 — PASS

All 8 in-hull cells inside the 0.003-nat envelope; the reference moves in every
cell (so every row discriminates — the script's standing power gate):

| alpha | eta | family | banked max | this run max | headroom vs 0.003 |
|---:|---:|---|---:|---:|---:|
| 2 | 1 | bip/cn | 5.68e-04 | 5.68e-04 | 5× |
| 2 | 2 | bip/cn | 1.24e-03 | 1.24e-03 | 2× |
| 10 | 1 | bip/cn | 6.59e-07 | 2.03e-06 | 1474× |
| 10 | 2 | bip/cn | 4.73e-06 | 4.77e-06 | 628× |

```
banked in-hull max surface_max : 0.001242389
new    in-hull max surface_max : 0.001242389
gold identical (bank lookup)   : TRUE
max |new - banked| surface_max over all 44 rows: 1.95e-06  (0.065% of the envelope)
rows where |new-banked| > 1e-5: 0
```

The banked control shape 2 reproduced to the digit. The alpha 10 / eta 1 rows
differ from the bank by ratios up to 109× — reported honestly, but the absolute
values are 1e-8 to 2e-6, i.e. ~1500× below the envelope, so this is Monte-Carlo
noise on an error that is effectively zero, not a deviation.

Runtime 62 s.

### 5. Collapse reachability — PASS, bit-identical

```
identical: TRUE      all.equal: TRUE
max abs diff k_star: 0
reachable agreement: TRUE
reachable at default delta in 17 of 40 (q, shape, eta) cells   (bank: 17 of 40)
```

The answer the script exists to give is unchanged: at `delta = 0.5*log(p)` the
additive zero-collapse **is** reachable for `q >= 20` in the shallow-shape cells
and from `q = 30` upward generally, and is unreachable for `q <= 16`.

Runtime 209 s.

---

## 5. Evidence — Task C: MM's sensitivity drifter, publication-grade

```r
fit  = bgm(Wenchuan, variable_type = "blume-capel", baseline_category = 1, seed = 1)
sens = prior_sensitivity_check(fit, seed = 1)
```

**Runtime: fit 37.9 s + check 84.6 s = 122.6 s total.** (Model resolved as
`omrf`; 5 NUTS refits warm-started from the fit.)

### Printed report, verbatim

```
Prior sensitivity check: are the edge verdicts robust to the slab scale?
Bayes-factor curve from 0.4x to 2.5x the chosen scale (anchors at 0.4x, 0.63x,
1x, 1.6x, 2.5x; the 1x anchor is the original fit); 136 edges.

92 of 136 verdicts hold across the whole 0.4x-2.5x range; the exceptions are
named below.

  robust (same verdict at every scale)     91
  changed, within run-to-run noise          1
  changed, beyond run-to-run noise         32
  not certifiable (too noisy to assess)    12

32 edges' verdicts genuinely depend on the scale:
  edge               0.4x       0.63x      1x         1.6x       2.5x
  intrusion-amnesia  undecided  absence    absence    absence    absence
  intrusion-lossint  undecided  absence    absence    absence    absence
  intrusion-numb     undecided  absence    absence    absence    absence
  intrusion-hyper    undecided  absence    absence    absence    absence
  dreams-avoidact    undecided  absence    absence    absence    absence
  flash-upset        undecided  absence    absence    absence    absence
  flash-avoidth      undecided  undecided  undecided  undecided  absence
  flash-anger        presence   presence   presence   presence   undecided
  flash-hyper        undecided  absence    absence    absence    absence
  upset-avoidact     undecided  undecided  undecided  absence    absence
  ...and 22 more; see $edges.

12 edges are too noisy to assess: intrusion-startle, dreams-upset, dreams-numb,
flash-sleep, upset-concen, physior-avoidact, avoidth-hyper, avoidact-lossint,
avoidact-startle, amnesia-lossint, and 2 more (see $edges).
Their Bayes factor sits within Monte Carlo error of an evidence threshold, or
their chains disagree on the verdict, at the chosen scale itself; a rerun with
a fresh seed could flip them without any prior change. Run more iterations to
settle these verdicts before reading their sensitivity.

Verdict counts by scale (at the anchors):
            0.4x 0.63x 1x 1.6x 2.5x
  presence    38    38 37   37   34
  undecided   67    47 38   30   29
  absence     31    51 61   69   73
More absence at wider scales is expected: a wider slab strengthens evidence
against borderline edges.

Note: the chosen scale (1) is much wider than the estimated interactions (about
0.195 [0.161, 0.236]); absence verdicts in particular depend on this choice.

Method:  43-point curve from 5 anchor fits (0.4x to 2.5x the chosen scale),
         joined by importance reweighting; the 1x anchor is the original fit.
         Points with reweighting effective sample size below 400 are not shown.
Refits:  5 nuts refits, warm-started from the original fit, 84 s total.
Noise:   two identical refits at 1.6x differed by up to 0.31 log10 BF across
         threshold-relevant edges; verdict moves smaller than that are reported
         as run-to-run noise, not prior sensitivity.
         2 edges saturated in one of the two and were left out of that spread.
See ?prior_sensitivity_check for the full construction.
```

### `intrusion-anger` specifically

**It does not appear in the mover list.** On the publication-grade fit it is
classified **robust**:

```
chosen-scale pip        : 0.9999754
chosen-scale log10 BF   : 4.608568
chosen-scale mcse       : 0.3255321
chosen-scale verdict    : presence
mover                   : stable
insufficient            : FALSE (noisy FALSE / disagree FALSE)
saturated               : FALSE
stability lower / upper : 0.4 / 2.5
CLASSIFICATION          : ROBUST
```

Per-anchor verdict — **presence at every anchor**:

| anchor | verdict | log10 BF | mcse |
|---:|---|---:|---:|
| 0.40 | presence | 6.000 | 0.0000 |
| 0.63 | presence | 5.627 | 0.2119 |
| 1.00 | presence | 6.000 | 0.0000 |
| 1.60 | presence | 5.601 | 0.2139 |
| 2.50 | presence | 2.136 | 0.9942 |

So: it never crosses a threshold, `stability_lower/upper` spans the full
0.4×–2.5× range, and it is neither within-noise, beyond-noise, nor
not-certifiable — it is **robust**.

### Curve shape vs MM's saved plot

MM's figure (`assets/sensitivity.pdf`, read on the Dropbox checkout — the blue
highlighted curve) runs roughly **6.0 → 6.0 → 4.7 → 1.5 → 1.0**, with the cliff
between 1× and 1.6×, ending in the undecided band.

The publication-grade fit runs **6.0 → 5.63 → 6.0 → 5.60 → 2.14**.

**The shape matches qualitatively — a high plateau, then one sharp cliff — but
the cliff has moved and shrunk.** From the full 43-point curve, the drop is
localised between multipliers 2.179 and 2.281:

```
 2.0814 4.961      1.8991 5.169
 2.1790 4.866   <- last point before the cliff
 2.2811 2.232   <- cliff
 2.3880 2.182
 2.5000 2.136
```

On the quick fit the same cliff sat between 1× and 1.6× and bottomed near
log10 BF 1.0 (undecided). On the good fit it sits at ~2.2× and bottoms at 2.14,
comfortably inside presence.

### Is it stable, or does it move run to run?

I ran a second `prior_sensitivity_check(fit, seed = 2)` on the **same** fit
(87.3 s) to separate the two.

| anchor | run 1 (seed 1) | run 2 (seed 2) | diff | verdict 1 | verdict 2 |
|---:|---:|---:|---:|---|---|
| 0.40 | 6.000 | 5.077 | −0.923 | presence | presence |
| 0.63 | 5.627 | 6.000 | +0.373 | presence | presence |
| 1.00 | 6.000 | 6.000 | 0.000 | presence | presence |
| 1.60 | 5.601 | 5.628 | +0.027 | presence | presence |
| 2.50 | 2.136 | 1.775 | −0.361 | presence | presence |

Classification is **ROBUST in both runs**. Headline counts barely move
(robust 91→90, within 1→1, beyond 32→33, not-certifiable 12→12). The cliff
reproduces at the same location (between 2.179× and 2.281×) with a similar floor
(2.14 vs 1.78).

**Answer to the brief's question: it is a stable property of the data–prior
combination, and the check is doing its job.** The 2.5× decline reproduces across
seeds in location and direction; only its depth wobbles by ~0.36 log10. What was
an artifact of the quick fit is the *severity* — the quick fit put the cliff at
1×–1.6× and pushed the edge into `undecided`, which is what made it look like a
large drifter. On a publication-grade fit the edge never leaves `presence`.
Nothing here is unstable run-to-run; **no finding on the drift itself.**

Two incidental notes, both verified: the low-scale end (0.4×–0.75×) is the noisy
region — run 1 sits at the cap while run 2 dips to 5.08 then recovers, and the
1× point is seed-independent by construction (it is the original fit, not a
refit), which is why it is identical across runs. And the flat-at-6.0 segments
are the clamp described in finding 04-5, not evidence plateaus.

---

## 6. Evidence — Task D: the hierarchical prior-chain cost (feeds F-036)

```r
fit = bgm(Wenchuan, variable_type = "continuous",
          precision_graph_prior = "hierarchical", seed = 1)
```

### The requested timings

| | |
|---|---:|
| fit | 6.16 s |
| **first `verdicts(fit)`** | **0.26 s** |
| **second `verdicts(fit)`** | **0.016 s** |
| third `verdicts(fit)` | 0.019 s |
| `verdicts()` after `saveRDS`/`readRDS` of the verdicted fit | 0.014 s |
| `verdicts()` on a reloaded fit that was **never** verdicted | 0.11 s |

### The cache

- It lives on **`fit$cache`, an R `environment`** (`bgms:::get_fit_cache()`;
  `is.environment(cache)` TRUE), keyed `prior_inclusion_class_values`.
- Before the first `verdicts()`: `edge_names, edge_selection, is_continuous,
  model_type, names_main, raw, summaries_computed`.
- After: the same **plus** `posterior_summary_indicator`,
  `posterior_summary_main`, `posterior_summary_pairwise`,
  `posterior_summary_quadratic`, `prior_inclusion_class_values`.
- **It does survive `saveRDS`/`readRDS`.** R serialises environments by value, so
  the reloaded fit carries the estimate: `prior_inclusion_class_values` present
  after `readRDS` = TRUE, values `all.equal()` to the original = TRUE, and
  `verdicts()` on the reloaded fit returns in 0.014 s with identical results.
- Caching gives a **~16×** speedup here (0.26 s → 0.016 s), and the saved fit is
  28.9 MB.

### The important correction for brief 03

**The 5000-iteration chain did not run in any of the above.** It is real —
`extract_prior_inclusion_probabilities()` defaults are `iter = 4000L` +
`warmup = 1000L` = 5000 — but the default hierarchical fit never reaches it.

Reading `prior_pip_class_values()` (`R/extract_prior_inclusion_probabilities.R`),
`prior_only_chain_pips()` is called only when either:

1. the edge prior is **Stochastic-Block** with a joint continuous block, or
2. `prior_pip_table()` returns `NULL` — i.e. the slab family has no tabulated
   correction.

For Bernoulli / Beta-Bernoulli with a **normal or cauchy** slab the code reads
the tabulated `edens` curve instead. I confirmed both slab families return a
table, so route 2 is unreachable through any documented `bgm()` slab argument.

Measured across four configurations:

| config | 1st `verdicts()` | 2nd | chain fires? |
|---|---:|---:|---|
| hierarchical / Bernoulli / normal slab *(the brief's config)* | **0.13 s** | 0.012 s | no — table lookup |
| hierarchical / Beta-Bernoulli | **0.11 s** | 0.012 s | no — table lookup |
| **hierarchical / Stochastic-Block** | **1.38 s** | 0.012 s | **yes** |
| **joint / Stochastic-Block** | **1.39 s** | 0.012 s | **yes** |

So, for brief 03's messaging:

- **The trigger is the edge prior, not the precision-graph prior.** Hierarchical
  and joint pay identically (1.38 vs 1.39 s); switching the precision-graph prior
  changes nothing. Describing this as a hierarchical-fit cost is wrong.
- **The measured cost is ~1.4 s at p = 26**, not a multi-second stall — and it is
  paid once, then cached, and the cache survives serialisation.
- It **is** silent (no message, no progress bar) and it **is** a 5000-iteration
  MCMC chain, so the transparency part of F-036 stands. The cost part should be
  restated.

(Cost will scale with `p`; Wenchuan at p = 26 is the only size measured.)

---

## 7. Open questions

1. **`test-prior-sensitivity.R:432` — tolerance or substance?** (finding 04-2)
   The gate is `4 × noise` where `noise` is a **two-refit** spread estimate; it
   failed at 5.05× noise. I cannot tell from one seed triple whether the
   difference-scale reweighting is genuinely biased by ~0.005 pip or the
   threshold is simply too tight to be stable. Needs MM: what is the intended
   resolution of this gate, and should `noise` come from more than two refits?
   A cheap diagnostic would be to run the triple over ~20 seed sets and look at
   the distribution of the ratio.

2. **The anchor-gate multiplier disagreement (finding 04-7) — is it accepted?**
   The gate re-derives 4× where `zratio_anchor_shape_multiplier()` ships 2× in
   five cells, and cannot reach parity at all in one (eta 2, cn, shape 0.5). This
   is identical in the bank, so it predates rc1 and someone presumably decided it
   was tolerable. But nothing in `dev/validation/` records that decision, and it
   is a live statistical question for a release: are the shipped anchor budgets
   knowingly below what the gate asks for? Purely MM's call.

3. **Ambiguous tolerance: the gold bank documents none directly.**
   `zratio_gold_bank.md` is the reference table; its only stated tolerance
   language is that the `sd` column *"bounds what any claim against these
   references can resolve."* The actual acceptance constants live in the scripts
   (`1e-12`, `2.84e-4`, `0.003`, `60 s`). I judged against those. If the lead
   expected a documented tolerance table in the `.md`, it is not there.

4. **Should the nightly gate PRs?** (finding 04-1) `nightly-validation.yaml` has
   no `push`/`pull_request` trigger, so branch work is never checked against the
   slow tier — which is exactly how three contradictory assertions merged. This
   is a process decision, not something I can settle. Related: the next scheduled
   run (Mon 2026-08-03 03:00 UTC) will be the first red one; if rc1 ships before
   then, the red will surface post-release.

5. **`test-ggm-nuts.R`'s "several minutes" header is stale** (0.6 s measured, work
   verified real). Not a finding, but if the tier is ever re-balanced by cost,
   that header will mislead.

---

## 8. Artifacts

All under `~/bgms-review/certification/`.

| path | contents |
|---|---|
| `suite-slow.log` | full slow-tier run, 274 KB, includes the verbatim failure and warning blocks |
| `test-results-slow.rds` | complete `testthat` result object for the slow run |
| `rerun-failing.log` / `.rds` | the three failing files rerun individually |
| `cert-deployed-route.log` | Task B certificate 1 |
| `cert-isolated-route.log` | Task B certificate 2 |
| `cert-anchor-gate.log` | Task B certificate 3 |
| `cert-shape-verdict-10.log` | Task B certificate 4 |
| `cert-collapse-reachability.log` | Task B certificate 5 |
| `banked/` | pristine copies of all 13 banked `.rds`, for comparison |
| `dev/validation/` | the working copy the scripts wrote to (originals untouched) |
| `cmp-verdict10.R`, `cmp-anchor-gate.R` | bank-vs-run comparison scripts |
| `task-c.log`, `task-c-analyse.log` | Task C, both seeds |
| `task-c-fit1.rds`, `task-c-run2.rds` | Task C fit and both sensitivity objects |
| `mm-sensitivity.pdf` / `.png` | MM's saved plot, copied and rasterised for comparison |
| `task-d.log`, `task-d2.log` | Task D timings |
| `task-d-timings.rds`, `task-d2-timings.rds` | Task D numbers |
| `slow-file-table.R` | per-file table generator |
| `~/bgms-review/run-slow-suite.R` | the runner (mirrors brief 01's `run-suite.R`) |

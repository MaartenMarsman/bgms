# Report 12 — nightly respecification + release-polish batch

Branch `fix/nightly-respec`, based on `origin/develop` @ `8dcdcf9d`, with
`origin/develop` @ `6ed36b44` (brief 10) merged in before the final gates.

> **Base note.** The brief says "develop AT OR AFTER `4969f843`". The Dropbox
> checkout's local `develop` was stale at `289d028d`, which does *not* contain
> `4969f843`. The worktree is based on `origin/develop` @ `8dcdcf9d`, which
> contains it plus the brief-12 amendment commits.

---

## The headline

**The old nightly never measured its own cost.** It died at
`timeout-minutes: 120` in the same place every run — part-way through the
`zratio-*` tail — so nine test files had never once been reached. Their runtime
was unknown, and so were their failures. Everything below follows from finally
getting a run to the end of the suite.

The respecified nightly (T1) now runs **green in 43.1 min** against its
60-minute budget, and the tier split is **38 T1 / 60 T2 = 98** — every
previously-slow block accounted for, none dropped.

---

## What was done

| # | Item | Commit |
|---|---|---|
| A3 | `BGMS_RUN_CERTIFICATION` gate (`tests/testthat/helper-tiers.R`) | `bde60827` |
| A2 | Re-tier the 98 slow-gated blocks | `5c2e8a95` |
| A4/A6 | `nightly-validation.yaml` rewrite; new `weekly-certification.yaml`; new `fast-checks.yaml` (D1) | `3f97da84` |
| A4 | Tier contract in `MAINTAINERS.md` §6 | `e63dec34` |
| B8 | F-072 group labels stored and surfaced on human displays | `f4c683a2` |
| B9 | F-057 residue: `?bgmCompare` group-numbering paragraph | `b6bb498c` |
| B10 | F-004: README development install gains `@develop` | `b6f237b5` |
| B7 | F-049 gate re-founding + `?prior_sensitivity_check` sentence | `8ef96121` |
| — | Merge brief 10 (`6ed36b44`) | `edfcfc26` |
| A5 | **T1 re-cut to budget** — 8 blocks T1→T2 | `cd457608` |
| B7 | F-049: pool the anchor side as well | `0f488b46` |
| — | **Two latent CI defects** the re-cut exposed | `2af3c68c` |

---

## Tier budgets, measured

All CI figures are from the 2-core `ubuntu-latest` runner.

| Tier | Cadence | Blocks | Measured | Budget |
|---|---|--:|---|---|
| T0 every run | push to develop + every PR | — | ~90 s of tests (`fast-checks.yaml`) | ~90 s |
| **T1 nightly** | daily 03:00 UTC, develop | 38 gated | **43.1 min wall, `success`** (run [30735988783](https://github.com/Bayesian-Graphical-Modelling-Lab/bgms/actions/runs/30735988783): setup 1.1 min + tests 41.7 min; `Duration: 2297.8 s`; `FAIL 0 \| WARN 0 \| SKIP 67 \| PASS 8460`) | ≤ 60 min, `timeout-minutes: 90` |
| **T2 weekly** | Sunday 03:00 UTC, develop | 60 gated | not yet run — see *Open questions* | `timeout-minutes: 360` |

Four T1 runs were needed:

| run | commit | wall | outcome |
|---|---|--:|---|
| [30715557912](https://github.com/Bayesian-Graphical-Modelling-Lab/bgms/actions/runs/30715557912) | `e63dec34` | 90.3 min | **cancelled at timeout** — over budget; first run ever to reach the `zratio` tail |
| [30719120127](https://github.com/Bayesian-Graphical-Modelling-Lab/bgms/actions/runs/30719120127) | `cd457608` | 37.6 min | in budget; 2 failures (finding 3) |
| [30720543250](https://github.com/Bayesian-Graphical-Modelling-Lab/bgms/actions/runs/30720543250) | `2af3c68c` | 59.8 min | green, but 12 s inside the budget — the package install cost 10.2 min of setup |
| [**30735988783**](https://github.com/Bayesian-Graphical-Modelling-Lab/bgms/actions/runs/30735988783) | `87b3eee1` | **43.1 min** | **green, in budget, ~17 min margin** |

**Runner variance is large and has to be budgeted for.** Runs 2 and 4 ran
identical test work and reported `Duration` 1993.5 s and 2297.8 s; run 3
reported 2716.9 s. That is a 36 % spread on the same suite, which is why 59.8
min was not accepted as "inside 60".

### How the first cut went wrong

The estimate from the cancelled develop run's per-file totals was ~45 min. It
came in at ~89 min. The whole discrepancy is the nine files that run had never
reached — chiefly `bgmCompare`'s full-defaults Boredom fit (828 s) and
`zratio-extrapolation-notice`'s size-cap fit (448 s). **A cancelled run gives
you a floor, not a cost.**

### The re-cut (8 blocks, T1 → T2)

| block | CI s | why |
|---|--:|---|
| `bgmCompare` — full-defaults Boredom fit | ~828 | **cost demotion.** A product-surface check by class, but a quarter of the budget by itself. The cheap end of the same surface stays T0 (the F-072 label test fits the same data at `iter = 50`). |
| `zratio-extrapolation-notice` — size-cap fit | ~448 | **cost demotion.** A heavy surface case, not a heartbeat. |
| `zratio-surface-build` — 3 build blocks | ~135 | contract-conformant: T1 keeps surface-vs-gold **single cells**; these are builds. |
| `zratio-gauge` — p=16 detector (2 blocks) | ~90 | **not cost — it fails on CI.** See finding 2. |
| `zratio-surface-extension` — block-Gibbs gold | ~31 | oracle machinery, not a heartbeat. |

Kept in T1 under the maintainer's ~2–3 min allowance: `zratio-law`'s fidelity
probe and two solver guards (~99 s), and `zratio-pivot-slice`'s bitwise pin and
two backstop guards (~86 s). The gauge-detector half of that allowance could
**not** be taken — it is red.

### Shared-fixture check

The audit's warning was respected: the marginal-Cauchy cell is still built
inside T1 by `hier-zratio-identity`'s two Cauchy identities and
`zratio-cauchy`'s cache-cell identity, so moving `zratio-surface-build`'s
Cauchy block (a *surface* cell, a different object) strands nothing and
relocates no cost.

---

## The classification table

Every one of the 98 previously-slow blocks, by the class of bug it catches.
`CI s` is the block's share of its file's CI time, apportioned by local block
share; `*` marks a block measured locally (serial, idle) and scaled by the
observed local→CI ratio (median 2.7, range 1.8–5.8).

#### T1 — 38 blocks (11.5 min of gated-block CI time)

| file | L | block | bug class | CI s |
|---|--:|---|---|--:|
| `hier-zratio-identity` | 156 | Cauchy graph law holds at dense high q (coupling regime) | graph-law identity | 132 |
| `mixed-correction` | 194 | corrected mixed prior chain returns the MFM partition prior | prior-chain identity | 85 |
| `bb-correction` | 85 | corrected prior-only chain returns the Beta hyperprior | prior-chain identity | 71 |
| `rb-inclusion-probabilities` | 127 | RB is interior for Monte-Carlo-saturated edges and never worse tha | RB saturation boundary | 58 |
| `hier-zratio-identity` | 24 | hierarchical graph marginal holds at gamma shapes | graph-law identity | 54 |
| `zratio-law` | 17 | law reproduces the companion eval_mu_law reference (fidelity) | graph-law identity | 50 |
| `sbm-correction` | 364 | corrected prior-only chain returns the MFM partition prior | prior-chain identity | 39 |
| `zratio-pivot-slice` | 99 | the pivot slice never exhausts its shrinkage budget | backstop guard | 37 |
| `hier-zratio-identity` | 135 | hierarchical graph marginal holds for the Cauchy slab | graph-law identity | 25 |
| `prior-sensitivity` | 385 | prior_sensitivity_check traces bgmCompare difference verdicts | compare surface trace | 24 |
| `bgm-hier-spec` | 248 | the hierarchical spec accepts a Cauchy slab on every update method | hierarchical engine breadth | 22 |
| `hier-zratio-identity` | 71 | hierarchical graph marginal matches Bernoulli(p); joint does not | graph-law identity | 13 |
| `zratio-pivot-slice` | 129 | the pivot slice handles a non-log-concave conditional | backstop guard | 10 |
| `bgm-hier-spec` | 273 | bgm fits the hierarchical spec and attaches the trust gauge | gauge attachment wiring | 10 |
| `hier-zratio-identity` | 190 | hierarchical BB identity: theta ~ Beta(a, b), PIP = a/(a+b) | graph-law identity | 10 |
| `mixed-correction` | 217 | uncorrected mixed prior chain under-segments the partition | prior-chain identity | 7 |
| `prior-sensitivity` | 195 | the chosen-scale verdict uses the per-edge prior odds, not 1/2 | fixed-bug regression guard | 7 |
| `zratio-law` | 80 | the gate rejects an unconverged (hard-corner) solve | solver gate guard | 6 |
| `zratio-pivot-slice` | 41 | the exponential shape is bit-identical to the pre-split kernel | bitwise kernel pin | 5 |
| `zratio-cauchy` | 96 | Cauchy constants live in their own cache cell | cache-cell identity | 4 |
| `ggm-nuts` | 198 | NUTS edge selection recovers true graph structure (p=6) | NUTS-vs-MH concordance | 3 |
| `ggm-nuts` | 164 | NUTS posterior means agree on p=6 tridiagonal | NUTS-vs-MH concordance | 2 |
| `mixed-correction` | 154 | corrected mixed prior chain returns the Beta hyperprior on theta | prior-chain identity | 2 |
| `rattle-edge-selection` | 85 | NUTS RATTLE interaction estimates agree with MH (p=5) | NUTS-vs-MH concordance | 2 |
| `sbm-correction` | 394 | uncorrected prior-only chain misses the partition prior | prior-chain identity | 2 |
| `ggm-nuts` | 51 | NUTS and MH posteriors agree on means and variances (p=4) | NUTS-vs-MH concordance | 2 |
| `zratio-law` | 95 | the law solve is deterministic | determinism guard | 2 |
| `bgm-hier-spec` | 337 | mixed data supports the hierarchical spec on the continuous block | hierarchical engine breadth | 2 |
| `mixed-correction` | 175 | uncorrected mixed prior chain biases theta toward sparsity | prior-chain identity | 1 |
| `ggm-nuts` | 109 | NUTS and MH 95% credible intervals overlap (p=4) | NUTS-vs-MH concordance | 1 |
| `ggm-nuts` | 384 | NUTS and MH posteriors agree under tilt (p=4, delta=1) | NUTS-vs-MH concordance | 1 |
| `hier-zratio-identity` | 115 | hierarchical graph marginal holds at a non-unit slab scale | graph-law identity | 1 |
| `bgm-hier-spec` | 152 | a vacuous hierarchical fit round-trips through a refit | refit round-trip wiring | 1 |
| `ggm-nuts` | 362 | NUTS diagnostics are clean for well-specified model | NUTS-vs-MH concordance | 0 |
| `rattle-edge-selection` | 141 | NUTS diagnostics are reasonable with edge selection (p=4) | NUTS-vs-MH concordance | 0 |
| `bb-correction` | 106 | uncorrected prior-only chain misses the hyperprior | prior-chain identity | 0 |
| `rb-inclusion-probabilities` | 206 | extract_inclusion_bf is finite for saturated edges and matches the | RB saturation boundary | 0 |
| `hier-zratio-identity` | 51 | gamma-shape constants build in the standardized cell | cache-cell identity | 0 |

#### T2 — 60 blocks (70.0 min of gated-block CI time)

| file | L | block | bug class | CI s |
|---|--:|---|---|--:|
| `zratio-extrapolation-notice` | 78 | the retained split reaches R from a fit past the size cap | surface size-cap notice | 448* |
| `validation-slow` | 43 | mixed MRF parameter recovery: cor > 0.8 (small network) | parameter-recovery / agreement sweep | 428 |
| `validation-slow` | 78 | MH vs NUTS posterior agreement: cor > 0.95 | parameter-recovery / agreement sweep | 305 |
| `sbc-correction` | 503 | SBC: corrected GGM SBM ranks are uniform | SBC replicate suite | 250 |
| `validation-slow` | 116 | estimate-simulate-re-estimate cycle: cor > 0.7 (mixed MRF) | parameter-recovery / agreement sweep | 195 |
| `scaling-diagnostics` | 281 | S.M3: Mixed NUTS healthy at p=7, q=5 with edge selection | NUTS health condition grid | 176 |
| `sbc-correction` | 249 | SBC: corrected GGM beta-bernoulli ranks are uniform (eta = 1) | SBC replicate suite | 149 |
| `sbc-correction` | 269 | SBC: corrected GGM beta-bernoulli ranks are uniform (eta = 3) | SBC replicate suite | 140 |
| `sbc-correction` | 259 | SBC: corrected GGM beta-bernoulli ranks are uniform (eta = 2) | SBC replicate suite | 129 |
| `prior-sensitivity` | 468 | the difference-scale reweighting reproduces a refit at that scale | refit cross-validation | 121 |
| `scaling-diagnostics` | 263 | S.M2: Mixed NUTS healthy at p=5, q=3 with edge selection | NUTS health condition grid | 121 |
| `scaling-diagnostics` | 304 | S.M4: Mixed NUTS healthy at p=5, q=3, marginal PL | NUTS health condition grid | 116 |
| `hier-zratio-identity` | 208 | hierarchical graph law across cells (slow battery) | graph-law condition grid | 113 |
| `bgmCompare` | 378 | the shipped data's own language column works as the group indicato | full-defaults product surface | 91 |
| `sbc-ggm` | 621 | SBC: GGM Gibbs produces uniform ranks at a gamma-shape diagonal | SBC replicate suite | 90 |
| `scaling-diagnostics` | 325 | S.M5: Mixed NUTS survives near-singular Kyy | NUTS health condition grid | 86 |
| `parameter-recovery-ggm` | 165 | PR.2: GGM parameter recovery, sparse graph with edge selection (p= | parameter-recovery sweep | 78 |
| `mixed-nuts` | 220 | M.2C: NUTS vs MH agree (conditional PL, ES, grouped) | NUTS-vs-MH condition grid | 73 |
| `mixed-nuts` | 265 | M.2D: NUTS vs MH agree (marginal PL, ES, grouped) | NUTS-vs-MH condition grid | 72 |
| `sbc-ggm` | 107 | SBC: GGM NUTS produces uniform ranks (p=3, no edge selection) | SBC replicate suite | 71 |
| `mixed-nuts` | 309 | M.2E: NUTS vs MH agree (conditional PL, ES, interleaved) | NUTS-vs-MH condition grid | 69 |
| `parameter-recovery-ggm` | 130 | PR.1: GGM parameter recovery, dense graph (p=5) | parameter-recovery sweep | 63 |
| `mixed-nuts` | 147 | M.2A: NUTS vs MH agree (conditional PL, no ES, grouped) | NUTS-vs-MH condition grid | 62 |
| `mixed-nuts` | 602 | M.2T: NUTS vs MH agree under tilt (conditional PL, no ES, delta=1) | NUTS-vs-MH condition grid | 62 |
| `sbc-ggm` | 433 | SBC: GGM NUTS produces uniform ranks under tilt (p=3, delta=1) | SBC replicate suite | 57 |
| `mixed-nuts` | 386 | M.2F: NUTS vs MH agree (conditional PL, no ES, interleaved) | NUTS-vs-MH condition grid | 55 |
| `sbc-ggm` | 521 | SBC: GGM joint-spec produces uniform ranks (p=5, edge selection) | SBC replicate suite | 52 |
| `sbc-correction` | 348 | SBC: corrected mixed prior chain returns Beta(1, 1) at the normal  | SBC replicate suite | 51 |
| `mixed-nuts` | 440 | M.2G: NUTS vs MH main effects agree (conditional PL, grouped) | NUTS-vs-MH condition grid | 50 |
| `mixed-nuts` | 185 | M.2B: NUTS vs MH agree (marginal PL, no ES, grouped) | NUTS-vs-MH condition grid | 50 |
| `scaling-diagnostics` | 245 | S.M1: Mixed NUTS healthy at p=3, q=2, no edge selection | NUTS health condition grid | 44 |
| `sbc-ggm` | 315 | SBC: GGM MH produces uniform diagonal ranks (p=3, edge selection) | SBC replicate suite | 42 |
| `prior-sensitivity` | 216 | warm-started short refits agree with cold full refits within wobbl | refit cross-validation | 39 |
| `sbc-ggm` | 182 | SBC: GGM MH produces uniform ranks (p=3, no edge selection) | SBC replicate suite | 37 |
| `zratio-surface-extension` | 77 | the extension beats the freeze against block-Gibbs gold | surface-vs-gold cell | 31* |
| `scaling-diagnostics` | 226 | S.G4: GGM NUTS healthy at p=15 with edge selection | NUTS health condition grid | 26 |
| `mixed-nuts` | 576 | M.2J: NUTS diagnostics are clean for mixed MRF | NUTS-vs-MH condition grid | 21 |
| `zratio-surface-build` | 149 | the build fences shapes outside the validated range | build fence guard | 16 |
| `mixed-nuts` | 503 | M.2H: coef/summary/simulate/predict work on mixed NUTS fit | NUTS-vs-MH condition grid | 16 |
| `mixed-nuts` | 547 | M.2I: simulate/predict preserve interleaved column order | NUTS-vs-MH condition grid | 15 |
| `zratio-law` | 49 | certified law cells match the all-MC oracle within noise | all-MC oracle cells | 14 |
| `zratio-gamma-shape` | 155 | clique-2 moments at alpha != 1 match an importance-sampled referen | MC channel reference | 12* |
| `zratio-pivot-slice` | 73 | the slice path reproduces the conjugate law it replaces | conjugate ground-truth oracle | 11 |
| `zratio-gamma-shape` | 81 | isolated-edge ratio psi0 at alpha != 1 matches direct Monte Carlo | MC channel reference | 11* |
| `scaling-diagnostics` | 210 | S.G3: GGM NUTS healthy at p=10 with edge selection | NUTS health condition grid | 8 |
| `scaling-diagnostics` | 194 | S.G2: GGM NUTS healthy at p=10, no edge selection | NUTS health condition grid | 6 |
| `sbc-correction` | 365 | SBC: uncorrected mixed prior chain misses the hyperprior at the no | SBC replicate suite | 6 |
| `zratio-gauge` | 235 | the deployed surface holds that fixture under the harm threshold | gauge detector | 5 |
| `zratio-surface-build` | 279 | the Cauchy slab builds and deploys its own surface cell | surface build/deploy wiring | 4 |
| `zratio-surface-build` | 37 | the built surface tracks the gold oracle far tighter than additive | surface-vs-gold cell | 3 |
| `zratio-gauge` | 215 | a known-biased evidence-free fit fires the harm channel | gauge detector | 2 |
| `zratio-gamma-shape` | 59 | spike ratio at alpha != 1 matches direct Monte Carlo | MC channel reference | 2* |
| `prior-sensitivity` | 255 | prior_sensitivity_check runs for GGM and mixed fits (cold refits) | refit cross-validation | 1 |
| `scaling-diagnostics` | 178 | S.G1: GGM NUTS healthy at p=5, no edge selection | NUTS health condition grid | 1 |
| `zratio-gamma-shape` | 128 | bridge channel at alpha != 1 matches direct Monte Carlo | MC channel reference | 1* |
| `zratio-gamma-shape` | 108 | node channel at alpha != 1 matches direct Monte Carlo | MC channel reference | 1* |
| `zratio-cauchy` | 47 | Cauchy bridge channel matches direct Monte Carlo | MC channel reference | 0 |
| `zratio-cauchy` | 77 | Cauchy isolated-edge ratio psi0 matches direct Monte Carlo | MC channel reference | 0 |
| `zratio-cauchy` | 25 | Cauchy node channel matches direct Monte Carlo | MC channel reference | 0 |
| `zratio-cauchy` | 111 | Cauchy exact Monte-Carlo evaluation tracks the additive prediction | MC channel reference | 0 |

---

## Findings

### 1 (high) — the cancelled run had three failures, not one; the tail hid more

The brief expected F-049 to be "the one red test". The 2026-08-01 develop run
carried three, and the first full run surfaced two more that no CI run had ever
executed:

| test | discovered | now |
|---|---|---|
| `mixed-nuts` M.2F | 2026-08-01 develop run | T2 — **F-080** |
| `sbc-ggm` diagonal ranks | 2026-08-01 develop run | T2 — **F-081** |
| F-049 reweighting gate | 2026-08-01 develop run | T2, re-founded here |
| `zratio-gauge` harm threshold | **run 30715557912** | T2 — **needs an owner** |
| `zratio-surface-build` PSOCK workers | **run 30719120127** | T0 — fixed here |
| `mcmc-diagnostics` RB-Rhat pin | **run 30719120127** | T0 — fixed here |

Per the lead's directive, F-080/F-081 were not chased.

### 2 (high) — `zratio-gauge` fails on Linux CI and passes locally

```
Failure ('test-zratio-gauge.R:243:3'): the deployed surface holds that fixture under the harm threshold
Expected `pc$harm_pred` < `f$zratio_diagnostics$harm_threshold`.
Actual comparison: 0.01071 >= 0.01000
```

A 7 % overshoot of the harm threshold; `harm_flag` flips to `TRUE`. It passes
on macOS. **This is not a cost demotion** — at ~1.5 min it fits the nightly
comfortably and the tier contract names the gauge detector a T1 concern. It is
in T2 only so the nightly can be green. It is the same species as F-080/F-081
and needs the same treatment: an owner and a decision on whether the threshold
or the fixture is wrong.

### 3 (high) — two more failures behind the old timeout; one fixed, one unexplained

Run 2 was the first in-budget run and it was red on two blocks neither new nor
previously reachable.

**`test-mcmc-diagnostics.R:434` — fixed.** The RB-Rhat identity was pinned at
`1e-8`. The two sides sum the same autocovariances in different orders, so the
last digits follow the machine's BLAS: the Linux runner lands at 3.1e-8 while
the pin holds on macOS. Widened to `1e-6` — still four orders tighter than the
df adjustment or wrong split the identity guards against, and the same
reasoning `test-zratio-pivot-slice.R` already records for its own pins.

**`test-zratio-surface-build.R:275` — intermittent, mechanism NOT established.**
It errored in run 2 with `4 nodes produced errors; first error: could not find
function "zratio_anchor_cn"`, then passed in runs 3 and 4. I first read this as
the documented cause — the block's PSOCK workers load the *installed* `bgms`
namespace (its own comment says so, and `zratio_build_surfaces()` ships a
closure calling that internal to `parLapplyLB`) — and added `local::.` to the
workflows. **That diagnosis does not hold up**: run 4 passed without
`local::.`, on a runner where the dependency step installs only `bgms-deps` and
no `bgms` or `easybgm` appears anywhere in the log. Two sub-hypotheses were
checked and are dead: a stale CRAN `bgms` pulled in transitively (nothing in
Suggests would do it, and no version string appears), and the block skipping
(it reports 3 skips in every run — the three T2 blocks — so it ran).

Also unexplained: run 2 reported **4** worker nodes, while the test asks for
`cores = 2L`.

What ships is the budget-preserving shape: the workflows install only
`devtools`, and the block asks a PSOCK worker whether `bgms` is loadable and
skips if not. That probe cannot turn a real failure into a pass — only an
error into a visible skip — so it is safe, but it is **not** a demonstrated
fix. Treat this as an intermittent needing an owner, alongside F-080/F-081.
Had `local::.` been kept, it would have cost ~10 min of every nightly.

### 4 (high) — the brief's prescribed F-049 fix does not work; pooling the anchor does

The brief specified pooling the *denominator* over four refit pairs. Measured
over 8 seed bases, that leaves the gate seed-fragile. The residual variance is
in the term neither the brief's fix nor my first extension pooled: **the
reweighting prediction still came from one anchor fit and carried its Monte
Carlo error whole.** Under the lead's authorization to pool the anchor side:

| construction | ratios over 8 bases | median | max | over ×4 |
|---|---|--:|--:|--:|
| 1 anchor, 1 refit (original) | 0.61 1.45 1.48 2.70 3.31 3.33 4.58 4.94 | 3.00 | 4.94 | **2/8** |
| 1 anchor, 8 refits | 0.78 1.00 1.53 2.01 2.89 3.52 3.87 5.44 | 2.45 | 5.44 | **1/8** |
| **3 anchors, 8 refits (shipped)** | 0.41 0.47 0.92 1.41 1.69 1.81 2.12 2.83 | **1.55** | **2.83** | **0/8** |

Target met: 8/8 green at the maintainer's ×4, with 29 % headroom on the worst
base. The pooled gap over those bases is 0.0018–0.0088 pip, inside the ~0.01
bound `?prior_sensitivity_check` now documents. Cost: 11 fits per run, in T2.

### 5 (medium) — `ggm-nuts`'s header was stale by two orders of magnitude

It claimed "Gated behind `BGMS_RUN_SLOW_TESTS` because they take several
minutes." Measured: **12.7 s** for the whole file. This decided its tier — by
bug class it is a NUTS-vs-MH agreement check like `mixed-nuts` (T2), but at 13 s
it is the concordance *smoke* the contract puts in T1. Header corrected.

### 6 (low) — the brief's illustrative mapping line has the groups swapped

The brief gives `"groups: 1 = en (n = 496), 2 = fr (n = 490)"` and, in the same
sentence, `(fr = 1 — the data is fr-first)`. `Boredom$language[1]` is `"fr"`, so
first-appearance numbering makes **fr group 1 (n = 490)**. Implemented as
first-appearance; the shipped line reads
`groups: 1 = fr (n = 490), 2 = en (n = 496)`.

### 7 (low, deviation) — D1's optional `R CMD check` half was skipped

`fast-checks.yaml` runs the T0 tier only. `R-CMD-check.yaml` already runs
`R CMD check` on the *same* push/PR triggers across five platforms; a second
copy would add cost and no signal. What was missing, and what this file adds,
is a fast pass/fail on the tests alone ahead of that matrix. Recorded in the
workflow header.

---

## Evidence

### Verification gate

| # | Gate | Result |
|---|---|---|
| 1 | Full **local** suite, `BGMS_RUN_SLOW_TESTS=true BGMS_RUN_CERTIFICATION=true` (= everything) | **`FAIL 0 \| WARN 0 \| SKIP 7 \| PASS 8942`** — zero failures, zero warnings |
| 2 | Dispatched T1 CI run, within budget, totals reported | **43.1 min wall, `success`** — run [30735988783](https://github.com/Bayesian-Graphical-Modelling-Lab/bgms/actions/runs/30735988783), `Duration: 2297.8 s`, `FAIL 0 \| WARN 0 \| SKIP 67 \| PASS 8460` |
| 3 | CRAN-mode suite; `R CMD check --as-cran` on a `git archive` tarball | **`FAIL 0 \| WARN 0 \| SKIP 254 \| PASS 7418`**; **`Status: 2 NOTEs`** — stale `Date` (F-003) and the environmental HTML-Tidy NOTE, i.e. the two baselines only |
| 4 | Classification covers 100 % of previously-slow blocks | **38 T1 + 60 T2 = 98.** Census script re-run against the final tree; no block appears in neither tier |

The 7 local skips are all pre-existing and unrelated to tiering: 5 × "golden
fixtures not found" (`tests/testthat/fixtures/`, `.Rbuildignore`d and absent
from a fresh worktree), 1 three-block-recovery case covered elsewhere, 1
warmup-efficiency case with no posterior claim to assert.

### F-049, all 8 seed bases

| base | identity gap | min ESS | pooled noise | gap | ratio |
|--:|--:|--:|--:|--:|--:|
| 11 | 0.0036 | 3231 | 0.00633 | 0.0030 | **0.47** |
| 31 | 0.0040 | 2291 | 0.00242 | 0.0051 | **2.12** |
| 51 | 0.0045 | 3106 | 0.00466 | 0.0079 | **1.69** |
| 71 | 0.0044 | 3388 | 0.00435 | 0.0061 | **1.41** |
| 91 | 0.0035 | 2825 | 0.00312 | 0.0088 | **2.83** |
| 111 | 0.0035 | 2885 | 0.00666 | 0.0061 | **0.92** |
| 131 | 0.0037 | 3491 | 0.00268 | 0.0049 | **1.81** |
| 151 | 0.0028 | 2615 | 0.00439 | 0.0018 | **0.41** |

Gate is `ratio < 4`: worst base 2.83, 29 % headroom. Identity gap ≤ 0.0045
against its 0.02 assertion; min ESS 2291 against its 400 assertion.

### Reproducing the measurements

- CI per-file times: `gh run view <id> --log`, per-file `[Ns]` from testthat's
  progress reporter (files run serially — no `Config/testthat/parallel`).
- Local per-block times: `testthat::test_file(reporter = ListReporter)` per
  file. **`NOT_CRAN=true` must be set** — without it `skip_on_cran()` silently
  skips most blocks and every timing is wrong. Local runs were contended (up to
  8 concurrent R processes on 15 cores), so they are used only as *within-file
  proportions* to apportion a file's CI total; the nine never-reached files
  were measured serially on an idle machine and scaled by the observed
  local→CI ratio (median 2.7, range 1.8–5.8).

---

## Open questions for the maintainer

1. **`zratio-gauge`'s harm-threshold failure** (finding 2) needs an owner. It is
   in T2 for greenness, not on merit; by the contract it belongs in T1.
1b. **`zratio-surface-build`'s PSOCK block** (finding 3) is an unexplained
   intermittent. It is in T0, so it gates every push and PR through
   `fast-checks.yaml`; if it recurs there it will be noisy.
2. **T2's budget is unproven.** `weekly-certification.yaml` is not dispatchable
   until it reaches `main`, because GitHub only exposes a workflow's triggers
   from the default branch. Its budget will be proven at its first real Sunday
   run. The gated T2 blocks alone total ~70 min of measured CI time, and the
   F-049 block grew from 3 fits to 11, so a 360-minute timeout has ample room —
   but that is an estimate, not a measurement.
3. **Equivocal classifications**, flagged per the brief:
   - `bgmCompare`'s full-defaults fit and `zratio-extrapolation-notice`'s
     size-cap fit are in T2 on **cost**, against their bug class. If the nightly
     budget is ever raised, these are the first two to bring back.
   - `ggm-nuts` (6 blocks) is the same bug class as `mixed-nuts` (T2) but 40×
     cheaper; it is in T1 as a smoke. If you would rather the tier be decided by
     class alone, it moves to T2 at a cost of 13 s.
   - `hier-zratio-identity`'s ten-cell battery went T2 as a condition grid while
     its single-cell siblings stayed T1. That is the contract's line, but the
     battery is the same identity.
4. **Schedule activation.** The new crons do not fire until this branch reaches
   `main` at the release re-merge; until then `main`'s Mon+Thu nightly keeps
   firing and cancelling. Known and expected.

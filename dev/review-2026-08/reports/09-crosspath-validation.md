# Report 09 — bgmCompare cross-path consistency validation (Opus agent)

Analysis only; no package code was changed.

## 1. What was done

### Environment

The local Dropbox `develop` did **not** contain merge `4969f843` (it was 29
commits behind), so the export was cut from `origin/develop` at
**`1bb1c7ca6ef80606d34a3eb51bff9182101e3ed9`** — which does contain
`4969f843` (verified with `git merge-base --is-ancestor`), plus the two
review commits that follow it. Nothing was built in the Dropbox tree.

```sh
mkdir -p ~/bgms-review/val09 && cd ~/bgms-review/val09
git -C <dropbox>/bgms archive origin/develop --prefix=bgms-val/ | tar -x
R CMD build bgms-val --no-build-vignettes --no-manual
R CMD INSTALL --library=~/bgms-review/lib-val09 bgms_0.2.0.0.tar.gz
```

R 4.6.0, bgms 0.2.0.0, darwin/arm64 (M5 Pro, 15 cores). Every fit below uses
`chains = 4, cores = 4` and, unless stated, `iter = 2000, warmup = 2000`
(package defaults). Scripts and raw logs live in `~/bgms-review/val09/out/`.

### Runtimes

| run | cell | wall |
|---|---|---|
| item 1 | 3 split-halves compare fits, p = 17, n = 344 | 67 / 69 / 79 s |
| item 1 | seed-repeat of split 101 | 80 s |
| item 2 | 4 planted cells (n = 400 sel T/F, n = 2000 sel T/F) | 184 / 160 / 575 / 1291 s |
| item 2 | `simulate_mrf()` mixing + round-trip checks | 130 / 190 s |
| item 2 | control: 2 × `bgm()` at n = 4000/group | 2367 s |
| item 2 | matched-n control: 2 × `bgm()` at n = 2000/group | 850 s |
| item 3 | 4 × `bgm()` half-fits + 1 full-data `bgm()` | 121 + 57 s |
| item 3 | 2 compare fits (sel on / off) | 78 / 79 s |
| item 4 | the every-run guard, reproduced verbatim | **2.3 s** |
| item 5 | compare Blume-Capel fit, p = 6 | 15 s |

### Derived tolerances

Every threshold below is measured, not assumed.

* **`bgm()` seed-to-seed refit noise** (item 3): three fits of the same half
  at seeds 11 / 13 / 14, all 136 pairwise means compared pairwise —
  `max|Δ| = 0.00211`, `sd = 0.00047`. Tolerance for compare-vs-separate
  agreement set at **3 × max = 0.0063**.
* **Compare sampler-seed noise on PIPs** (item 1): the same split refitted at
  seed 777 — `max|Δpip| = 0.0884`, `rms = 0.0169`; on difference posterior
  means `max|Δ| = 0.0181`.
* **Null noise floor on difference estimates** (item 2): the 41 unplanted
  edges, whose truth is exactly 0 — `rmse = 0.262` at n = 400/group,
  `0.111` at n = 2000/group. This is the yardstick the planted δ have to
  clear, and three of the four do not.
* **False-positive budget** (item 1): `difference_prior = bernoulli_prior(0.5)`
  gives prior odds 1, so `evidence_threshold = 10` calls presence at
  `PIP > 10/11 = 0.9091`. The Bayes-factor bound `P(BF > 10 | H0) ≤ 1/10`
  gives **≤ 13.6 of 136** pairwise differences.

---

## 2. Findings

### 09-1 · `note` (positive, and the headline) · No ×2 or ×½ convention divergence anywhere on the compare path

The resident defect class is absent. Regressing bgmCompare's reconstructed
per-group pairwise associations on matched separate `bgm()` fits of the same
rows (`difference_selection = FALSE`, priors matched, split seed 101):

| quantity | slope through origin | 95% CI | distance from 2 | distance from ½ |
|---|---|---|---|---|
| group 1 | **1.026** | [1.019, 1.034] | −0.974 | +0.526 |
| group 2 | **1.100** | [1.090, 1.110] | −0.900 | +0.600 |
| contrast (g1 − g2) | **1.014** | [0.994, 1.034] | −0.986 | +0.514 |
| level ((g1+g2)/2) | **1.064** | [1.060, 1.069] | −0.936 | +0.564 |

The contrast row is the one that matters most and is the cleanest: the
**difference parameter is on exactly the same scale as the separate-fit
contrast**, slope 1.014 with an interval that contains 1 and whose nearest
edge sits 25 half-widths from 0.5. Correlations are 0.9935–0.9997.
Independently, item 5 reproduces
`predict.bgmCompare()` from a hand-written linear predictor carrying
`2 · omega · x` to **3 × 10⁻¹⁶** on all six variables in both groups, mixed
Blume-Capel and ordinal; the same reference with the rest score multiplied by
2 or by ½ disagrees by 0.79 and 0.66 probability. Item 2's simulate/predict
round-trip and item 1's null all agree.

### 09-2 · `minor` · The compare fit's pairwise **level** sits 6–9 % above matched `bgm()` fits, and pooling does not explain it

Compare-vs-separate `max|Δ|` is 0.029 (group 1) and 0.082 (group 2), against
the derived 0.0063 tolerance — 4.6× and 13× the `bgm()` refit noise. So the
literal item-3 agreement check **FAILS**, and it fails in a specific,
narrow way: the *scale* is right (09-1) but the *magnitude* runs high.

The obvious explanation — the compare fit's shared baseline pools both halves
and is therefore less prior-shrunk than a half-sized `bgm()` — is **refuted by
the control**. Against one `bgm()` fit of all 344 rows:

| | slope through origin | max abs deviation |
|---|---|---|
| compare level ~ full-data `bgm()` | 1.088 [1.068, 1.109] | 0.0535 |
| mean of the two half `bgm()` fits ~ full-data `bgm()` | 1.022 [1.002, 1.042] | 0.0389 |

If pooling were the mechanism, the compare level would sit *closer* to the
full-data fit than the half-fits do. It sits further. The two half-fits are
themselves 2.2 % high against the full-data fit, so the direction is shared
and the compare path roughly quadruples it. I have not established the
mechanism and am not going to guess at one; the observation is that a
bgmCompare group estimate is a few percent larger in magnitude than the
`bgm()` estimate a user would get on the same rows, consistently in one
direction. At the magnitudes involved (0.08 at the worst of 136 edges on
n = 172 per group) this changes no verdict, but it is the kind of systematic
offset the review exists to surface. **Needs a maintainer call** on whether
it is expected behaviour of the shared-baseline parameterisation.

### 09-3 · `major` · On identical simulated data at identical n, bgmCompare's per-group pairwise estimates are ~4× noisier than two separate `bgm()` fits, and group 2's slope reads 1.93 where `bgm()` reads 1.03

This is item 2's operating-characteristics result and the one I would not
sign off without a maintainer reading. Truth is planted on the association
scale as `group2 − group1 = δ` (the fit's own projection is `[−0.5, +0.5]`,
so the stored difference parameter *is* the group contrast — verified to
machine precision against `extract_group_params()` in every cell).

At n = 2000 per group, `difference_selection = TRUE` (defaults):

| planted edge | δ | posterior mean | 95 % interval | covers δ | PIP | verdict |
|---|---|---|---|---|---|---|
| intrusion-upset | 0.05 | 0.000 | [0.000, 0.000] | no | 0.019 | absence |
| dreams-avoidact | 0.10 | 0.402 | [0.000, 0.897] | yes | 0.780 | undecided |
| upset-avoidact | 0.20 | 0.019 | [−0.267, 0.405] | yes | 0.208 | undecided |
| avoidact-lossint | **0.40** | **1.212** | **[0.697, 1.732]** | **no** | 1.000 | presence |

The δ = 0.40 row is the concerning one: the detection is right (PIP 1.000,
presence, the only presence among the four), but the estimate is **3.0× the
truth** and the interval misses it entirely, and it does *not* improve from
n = 400 to n = 2000 (0.855 → 1.212, and 1.276 with selection off at n = 400).
Ratios across the four planted edges are 0.0, 4.0, 0.09, 3.0 — **not a
constant factor**, so this is not the ×2 class of 09-1; it is a
magnitude/identifiability problem, not a units problem.

The mitigating context, which is substantial: the estimator's own noise floor
on a difference whose truth is exactly 0 is `rmse = 0.111` at n = 2000/group
and `0.262` at n = 400/group. Three of my four planted δ (0.05, 0.10, 0.20)
sit **at or below that floor**, so their non-recovery is a power statement
about my design, not about the code. Only δ = 0.40 clears it — and that one
is detected but over-estimated 3×. Aggregate coverage of the planted truth
across all 45 edges is 0.933 (42/45) at n = 2000 and 0.911 at n = 400 with
selection on, close to nominal, because 41 of the 45 truths are 0 and those
are covered 40/41.

**What I ruled out, and what is left.** I chased four explanations that
would have made this a design artefact rather than a package problem. Three
are dead:

1. **`simulate_mrf()` mixing failure** — refuted. Converged by `iter = 250`;
   two draws from the same omega differ only by sampling noise that halves
   as n grows (§3.2).
2. **A broken `bgm()` → `simulate_mrf()` round trip** (i.e. my planted truth
   not being the truth) — refuted. Simulated correlations reproduce real
   Wenchuan's to a slope of 0.974, mean abs error 0.014.
3. **An ill-conditioned near-critical regime** — refuted, and this was my own
   wrong turn: I read the distribution's sharp sensitivity to omega as
   ill-conditioning, when it means the opposite. `bgm()` recovers omega from
   data simulated at exactly this operating point with **slope 1.0007
   [0.975, 1.027], rmse 0.0123**.
4. **A fixed units defect** — refuted by 09-1's four independent checks, and
   independently by the fact that which group is inflated *flips* between the
   n = 400 and n = 2000 cells, which no fixed convention error can do.

**The bgmCompare-free control settles the design question: the planted truth
is fully recoverable.** Two independent `bgm()` fits of the same planted
`OM1` / `OM2`, n = 4000 per group (2367 s):

```
bgm() group 1 ~ planted OM1        slope 0.996 [0.964, 1.028]  cor 0.993  rmse 0.0152
bgm() group 2 ~ planted OM2        slope 0.968 [0.920, 1.016]  cor 0.982  rmse 0.0242
bgm() contrast ~ planted delta     slope 0.883 [0.774, 0.992]  cor 0.924  rmse 0.0259
```

Per planted edge, `bgm()` reads the contrast off almost exactly, and its
noise floor on the 41 true-zero contrasts is small:

| planted edge | truth δ | `bgm()` contrast, n = 4000/gp | bgmCompare difference, n = 2000/gp |
|---|---|---|---|
| intrusion-upset | 0.05 | **0.064** | 0.000 |
| dreams-avoidact | 0.10 | **0.094** | 0.402 |
| upset-avoidact | 0.20 | **0.196** | 0.019 |
| avoidact-lossint | 0.40 | **0.339** | 1.212 |
| unplanted (truth 0) | — | sd **0.0255**, max **0.0648** | rmse **0.111**, max **0.677** |

So every design-side explanation is exhausted: the simulation mixes, the
round trip is faithful, the regime is informative, and an independent
single-group fit recovers all four planted δ to within 0.06.

**The matched-n control removes the last confound and confirms the gap.**
`bgm()` fitted separately to the *byte-identical* n = 2000/group datasets
from item 2's own cell (same data seeds 2500 / 2501, 850 s):

| quantity, n = 2000/group | two independent `bgm()` fits | bgmCompare (sel = T) | bgmCompare (sel = F) |
|---|---|---|---|
| group 1 slope vs OM1 | **0.963** [0.921, 1.005] | 1.270 [1.159, 1.382] | 1.295 [1.170, 1.420] |
| group 1 rmse | **0.0207** | 0.0657 | 0.0729 |
| group 2 slope vs OM2 | **1.034** [0.931, 1.136] | **1.926** [1.599, 2.253] | **1.986** [1.636, 2.337] |
| group 2 rmse | **0.0511** | 0.2144 | 0.2291 |
| unplanted contrast noise (truth 0) | sd **0.0414**, max **0.1104** | — | sd 0.1559, max 0.8299 |

On identical rows, at identical n, the compare path is **3.8× noisier** on
true-zero differences (sd 0.156 vs 0.041, max 0.830 vs 0.110) and **4.2×
worse in rmse** on group 2, whose slope sits at 1.93 [1.60, 2.25] where two
separate fits give 1.03 [0.93, 1.14] — non-overlapping intervals. That is the
finding, and it is not attributable to sample size, to the simulation, to the
regime, or to a units convention.

**Two honest qualifications.** First, `bgm()` is not itself clean on the
*contrast* at this n: its contrast slope is 1.500 [1.314, 1.687] (per-edge
−0.001, 0.097, 0.288, 0.629 against truth 0.05, 0.10, 0.20, 0.40), against
0.883 at n = 4000. So some contrast inflation at n = 2000/group is a
small-sample property both paths share, and only the *group-level* slopes and
the noise inflation are compare-specific. Second, and more important, **this
rests on one simulated dataset per cell.** The group-1/group-2 asymmetry
flips between the n = 400 and n = 2000 cells, which says run-to-run
variability is large; a single seed is not a calibration study. Replication
across data seeds is the prerequisite before anyone scopes a fix.

Item 1 is genuine evidence on the other side: on *real* Wenchuan rows split
at random, the compare path produced difference means with sd 0.046–0.058
and near-nominal interval behaviour. Whatever §3.2 is showing does not
reproduce on real data at n = 172 per group.

Independently of which way that resolves, the *calibration* question stands:
whether a Cauchy(0, `difference_scale = 1`) slab on a parameter whose
realistic scale is ~0.1 is the right default is exactly the open question
`verdicts()` already prints a caveat about.

### 09-4 · `minor` · The every-run convention guard is blind to a ×½ divergence, and to every difference-parameter convention

`tests/testthat/test-bgmCompare.R:337`, run verbatim (2.3 s, seed 1234):
group 1 recovers `(0.626, −0.079, 0.586)` against a truth of
`(0.5, 0, 0.45)`. Applying the two assertions to hypothetically diverged
estimates:

| divergence applied to the estimate | `rmse < 0.2` | `rmse < 0.5·rmse(2·target)` | guard |
|---|---|---|---|
| as shipped (×1) | 0.116 PASS | 0.143 PASS | passes |
| **×2 (the historical defect)** | 0.609 FAIL | FAIL | **catches it** |
| **×½ (a double-applied correction)** | **0.143 PASS** | **0.265 PASS** | **misses it** |

Idealised (estimate exactly ½ the truth) the margin is just as comfortable:
rmse 0.194 against a 0.2 threshold, and 0.194 against a 0.291 relative bound.
Both assertions pass a ×½ divergence with room to spare. The guard is
one-sided by construction: the `0.5 * rmse(2 * target)` comparison only ever
tests the doubling direction.

Four further coverage holes, all read-only:

1. **Group 2 is never asserted.** The guard reads
   `pairwise_effects_groups[, 1]` only. In the observed run group 2 came out
   `(0.513, 0.125, 0.317)` — visibly worse than group 1 against the same
   truth, and entirely unchecked. Since both groups are drawn from the same
   omega, a sign error in the projection's second row would leave every
   assertion intact.
2. **No every-run test constrains the scale of a *nonzero* group
   difference.** The guard's two groups share one omega, so its difference
   parameters are noise around zero (observed: −0.113, 0.203, −0.269), and
   `difference_selection = FALSE` removes the indicators. A ×2 introduced in
   the difference parameterisation alone would pass the whole every-run tier.
   This is the largest hole, and it is precisely the parameterisation the
   0.2.0.0 breaking change touched.
3. **The anchor is relative, not absolute.** The guard pins bgmCompare
   against `simulate_mrf()`. A change flipping *both* conventions together
   passes. (`simulate_mrf()`'s own convention is pinned independently — hand
   computations in `test-mixed-mrf-likelihood.R` — so the composition is
   sound today; the point is that this test alone does not anchor it.)
4. **Main effects are not pinned at all** by the guard; only pairwise.

On the brief's three specific questions:

* **(a) pairwise likelihood** — caught in the ×2 direction, missed in ×½,
  and only for group 1's baseline.
* **(b) simulate / predict** — the guard does not touch either.
  `predict.bgmCompare()` is well covered by the F-068 test at line 387, which
  is a full manual reference carrying `2 * pw` and would catch a ×2 in either
  direction — but it is Blume-Capel-only and group-1-only, and it compares
  predict against `extract_group_params()`, so a divergence living upstream
  in the sampler moves both sides together and survives.
  **`simulate.bgmCompare()` has no numeric-convention test at all.** Its
  coverage is input validation (`test-input-validation.R:433+`) plus purely
  structural assertions in `test-methods.R:670-745` — dimensions, column
  names, integrality, in-range category values, and seed reproducibility.
  The one test that could have compared groups
  (`"simulate.bgmCompare produces different results for different groups"`,
  line 721) asserts only `is.matrix()` and equal dimensions; its own comment
  calls it "a soft test — we just verify they can be different". Nothing
  anywhere pins a simulated margin to a model expectation. Item 5 fills that
  gap empirically (max |simulated margin − predicted expectation| = 0.0083
  over 12 variable × group cells at nsim = 20000, against a binomial MC error
  of ~0.003); there is no test pinning it.
* **(c) mixed cross terms** — **not reachable on the compare path.**
  `R/bgm_spec.R:418` sets `allow_continuous = (model_type != "compare")`, so
  `bgmCompare()` rejects continuous variables outright. No compare test can
  cover mixed cross terms and none needs to; that convention is pinned in
  `test-mixed-mrf-likelihood.R` against hand computations, on the `bgm()`
  path only.

### 09-5 · `note` · The compare path passes the null cleanly, but its strongest null false positives are very strong

Item 1 is a pass: 4 / 3 / 2 presence verdicts out of 136 pairwise differences
across three split seeds (2.9 %, 2.2 %, 1.5 %), all inside the ≤ 13.6 budget;
posterior mean differences straddle zero (mean −0.001 to −0.002, sign split
69/76/69 of 136 against a binomial expectation of 68 ± 5.8); 3/2/2 of 136
95 % intervals exclude zero, below the nominal 5 %.

The decisive diagnostic is that these presences are **not structural**: no
edge reaches presence in more than one split, and cross-split PIP
correlations are 0.025, 0.037, 0.009 — indistinguishable from zero. The
machinery is not systematically manufacturing differences.

What is worth recording is the *strength* of the tail. On split 103,
`sleep-startle` reached PIP 0.9999998, **log BF = 15.41**, non-fragile, on
data with no true difference whatsoever; split 101's `upset-numb` reached
log BF 7.60. Nine of 408 null edge-fits crossed the presence threshold and
the largest crossed it by 13 nats. None was flagged fragile. This is an
operating-characteristic note tied directly to the `difference_scale`
calibration caveat `verdicts()` already prints, not a demonstrated defect.

### 09-6 · `note` · `arguments$baseline_category` stores a meaningless value for ordinal variables

On a mixed fit (`variable_type = c("blume-capel", "blume-capel", rep("ordinal", 4))`,
`baseline_category = 2`) the stored vector is `1 1 2 2 2 2`: the two
Blume-Capel entries are correctly shifted to the 0-based scale, and the four
ordinal entries carry the raw input value 2, which is not a baseline of
anything. It is inert — `src/mrf_prediction.cpp:99-105` reads
`baseline_category[vertex]` only inside the `blume-capel` branch, and item 5's
manual reference (which centres ordinal neighbours at 0) matches `predict()`
to 3 × 10⁻¹⁶, confirming the ordinal entries are never used. Recording it
because it is a live trap for anyone writing a reference implementation
against `extract_arguments()`.

### 09-7 · `note` (positive) · Blume-Capel end-to-end on the compare path is clean

`verdicts()`, `predict()`, `calibration_check()` per group, and `simulate()`
all behave. Details in §3.5; the F-068 surface shows no regression, and
per-group calibration is indistinguishable from `bgm()` fitted to the same
halves.

---

## 3. Evidence

### 3.1 Item 1 — split-halves identity

`Wenchuan` complete cases (n = 344, p = 17, all 5-category), split
50/50 by `set.seed(ss); rep(1:2, length.out = n)[sample(n)]` for
ss ∈ {101, 102, 103}; `bgmCompare(W, group_indicator = grp, seed = 2026)` at
defaults (`difference_selection = TRUE`, `main_difference_selection = FALSE`,
`difference_prior = bernoulli_prior(0.5)`, `difference_family = "Cauchy"`,
`difference_scale = 1`). Convergence: max indicator Rhat 1.028 / 1.032 /
1.031, min n_eff 310 / 309 / 326.

| | split 101 | split 102 | split 103 |
|---|---|---|---|
| presence | 4 | 3 | 2 |
| undecided | 54 | 49 | 62 |
| absence | 78 | 84 | 72 |
| max PIP | 0.9995 | 0.9942 | 1.0000 |
| # PIP > 0.5 | 16 | 11 | 11 |
| median PIP | 0.073 | 0.067 | 0.080 |
| difference means, range | [−0.332, 0.222] | [−0.233, 0.285] | [−0.264, 0.129] |
| difference means, mean | −0.00107 | −0.00046 | −0.00203 |
| # positive (of 136) | 69 | 76 | 69 |
| # 95 % CIs excluding 0 | 3 | 2 | 2 |

The nine edges that reached presence in any split, showing that none repeats:

```
                             split101 split102 split103
intrusion-physior              0.2974   0.0779   0.9880
intrusion-concen               0.9872   0.0951   0.0555
intrusion-startle              0.0685   0.9889   0.1759
dreams-lossint                 0.9869   0.0615   0.2359
dreams-numb                    0.6457   0.9942   0.2828
flash-numb                     0.0640   0.9653   0.7460
upset-numb                     0.9995   0.0510   0.6008
avoidth-avoidact               0.9553   0.0555   0.2067
sleep-startle                  0.1155   0.0376   1.0000

edges reaching presence in more than one split: 0
cross-split PIP correlation: 101-102 0.025, 101-103 0.037, 102-103 0.009
```

Presence verdicts with their evidence (all non-fragile):

```
split 101: intrusion-concen 4.34 | dreams-lossint 4.33 | upset-numb 7.60 | avoidth-avoidact 3.06
split 102: intrusion-startle 4.49 | dreams-numb 5.14 | flash-numb 3.33
split 103: intrusion-physior 4.41 | sleep-startle 15.41      (natural-log BF)
```

Verdict: **PASS** on the budget and on the cross-split independence
diagnostic; the log-BF-15 tail is 09-5.

### 3.2 Item 2 — planted-difference recovery

Truth: `bgm(Wenchuan[, 1:10], edge_selection = FALSE, interaction_prior = cauchy_prior(1), seed = 31)`
gives baseline omega (range [−0.062, 0.631], median 0.050) and a 10 × 4
threshold matrix. Differences planted on pairs 3, 14, 27, 41 as
`OM1 = OM − δ/2`, `OM2 = OM + δ/2` with δ = (0.05, 0.10, 0.20, 0.40).
Data from `simulate_mrf(n, 10, num_categories = 4, ..., iter = 1000)` — which
runs one independent Gibbs chain per row (`src/mrf_simulation.cpp:88-98`), so
rows are iid, not an autocorrelated trace.

**The design is sound.** Simulated margins reproduce the real Wenchuan
margins closely (e.g. `intrusion` 0.066/0.393/0.210/0.222/0.108 simulated
against 0.078/0.378/0.224/0.212/0.108 observed), so the truth is a realistic
network, not a degenerate one.

Posterior mean of the difference parameter, all four cells:

| planted | δ | n=400 sel=T | n=400 sel=F | n=2000 sel=T | n=2000 sel=F |
|---|---|---|---|---|---|
| intrusion-upset | 0.05 | −0.015 | −0.222 | 0.000 | 0.011 |
| dreams-avoidact | 0.10 | 0.642 | 0.777 | 0.402 | 0.537 |
| upset-avoidact | 0.20 | 0.017 | 0.160 | 0.019 | −0.040 |
| avoidact-lossint | 0.40 | 0.855 | 1.276 | **1.212** | **1.042** |
| unplanted rmse (truth 0) | — | 0.185 | 0.262 | 0.111 | 0.157 |
| unplanted max abs mean | — | 1.092 | 1.279 | 0.677 | 0.830 |
| coverage, all 45 edges | — | 0.911 | 0.800 | 0.933 | **0.689** |
| AUC planted vs unplanted | — | 0.768 | — | 0.720 | — |

Selection *off* is worse, not better, in every cell, so the inflation is not
a spike-and-slab selection artefact. The n = 2000 selection-off cell is the
worst: 95 % intervals cover the truth on only 31 of 45 edges (0.689), and
28 of the 41 edges whose truth is *exactly zero* — a prior centred at zero
should over-cover a zero truth, not under-cover it by 27 points.

Decomposing each cell into its two group-level estimates against the planted
`OM1` / `OM2` (all four cells confirm `dm == g2 − g1` to ≤ 9 × 10⁻¹⁵, so the
projection convention is exact):

| cell | group 1 slope vs OM1 | group 2 slope vs OM2 |
|---|---|---|
| n=400 sel=T | 1.802 [1.450, 2.155] | 0.977 [0.793, 1.161] |
| n=400 sel=F | 1.968 [1.473, 2.464] | 0.983 [0.763, 1.203] |
| n=2000 sel=T | 1.270 [1.159, 1.382] | 1.926 [1.599, 2.253] |
| n=2000 sel=F | 1.295 [1.170, 1.420] | 1.986 [1.636, 2.337] |

These are incoherent as a convention story: which group is inflated *flips*
between the n = 400 and n = 2000 cells, and the two groups of a single fit
differ by non-overlapping intervals. A units defect cannot do that.

**`simulate_mrf()` mixing is ruled out as the explanation.** Two draws from
the *same* omega differ only by sampling noise, and that noise halves as n
grows the way iid draws must:

```
A. two samples, same omega, different seeds (truth: identical)
   n = 2000 : max |cor(a) - cor(b)| = 0.0558  (mean 0.0145)
   n = 8000 : max |cor(a) - cor(b)| = 0.0237  (mean 0.0093)

B. Gibbs convergence, same seed, increasing iter (truth: stabilises)
   iter =   250 : mean cor 0.4923      iter =  4000 : mean cor 0.4878
   iter =  1000 : mean cor 0.4859      iter = 16000 : mean cor 0.4887
```

The chain is converged by `iter = 250`; the shipped `iter = 1000` is ample.
What part C does show is that the planted δ propagate: at n = 2000/group the
two datasets differ by max |Δcor| = 0.261 against a same-omega reference of
0.046, and the largest marginal shift sits on an *unplanted* pair (planted
edges shift by 0.138, 0.021, 0.070, 0.233). That is expected — the model is
parameterised conditionally, so a conditional change on four edges moves
marginal correlations everywhere — and it does not invalidate the planted
truth in model space.

**The `bgm()` → `simulate_mrf()` round trip is clean on the pairwise term**,
which both validates the planted truth and adds a fourth independent
confirmation of the association-scale convention on a path the brief did not
ask about:

```
real Wenchuan pairwise correlations : mean 0.5042  range [0.366, 0.802]
simulated (n = 20000) from the fit  : mean 0.4901  range [0.350, 0.792]
mean |simulated - real| = 0.0144 ; max = 0.0380 ; slope sim ~ real = 0.9735
```

The operating point is a sensitive one — rescaling omega moves the
distribution sharply:

```
  omega x 0.5 : mean cor 0.0313   (real 0.5042)   mean |diff| 0.4729
  omega x 1.0 : mean cor 0.4890   (real 0.5042)   mean |diff| 0.0153
  omega x 2.0 : mean cor 0.0007   (real 0.5042)   mean |diff| 0.5036
```

I first read that sensitivity as ill-conditioning that would invalidate the
design. **That reading is wrong, and the next check refutes it.** Sharp
dependence of the distribution on the parameters means the data are highly
*informative* about them — large Fisher information, so estimation is
easier, not harder. And indeed `bgm()` recovers omega from data simulated at
exactly this operating point essentially perfectly:

```
bgm() refit on its own simulated data, n = 4000:
  slope = 1.0007 [0.9748, 1.0267]   cor 0.9950   rmse 0.0123
```

So the regime is fine, the truth is valid, and a single-group fit of the same
network at comparable n lands on slope 1.000 with rmse 0.012 — against
bgmCompare's group slopes of 1.270 and 1.926 and rmse 0.066 and 0.214 in the
n = 2000 cell. That gap now points at the compare path rather than at the
design, which is the opposite of my earlier reading in this same section.

The last confound — the self-recovery used n = 4000 in one homogeneous
sample, while the compare cell used n = 2000 per group — is removed by the
matched-n control tabulated in 09-3: on the byte-identical n = 2000/group
datasets, two separate `bgm()` fits reach group slopes 0.963 and 1.034 with
rmse 0.021 and 0.051 and an unplanted-contrast noise floor of sd 0.041,
against bgmCompare's 1.270 / 1.926, rmse 0.066 / 0.214, and sd 0.156.

Verdict: **detection PASS at δ = 0.40** (PIP 1.000, the only presence among
the four planted edges at n = 2000); **magnitude FAIL, and compare-specific**
— see 09-3. AUC of 0.72–0.77 is weak, but three of four planted δ sit under
the estimator's own noise floor, so the AUC understates the achievable
separation. All magnitude conclusions rest on one dataset per cell and need
seed replication before they are acted on.

### 3.3 Item 3 — per-group estimates vs separate single-group fits

Fixed split (seed 101, n = 172 per half). Priors matched explicitly, because
**the two entry points ship different pairwise defaults**: `bgm()` defaults to
`normal_prior(scale = 1)` (`R/bgm.R:481`) and `bgmCompare()` to
`cauchy_prior(scale = 1)`. Both sides were pinned to
`interaction_prior = cauchy_prior(1), threshold_prior = beta_prime_prior(0.5, 0.5)`
and `bgm(edge_selection = FALSE)` so no spike-and-slab acts on the pairwise
level on either side.

```
bgm() seed-to-seed noise, half 1, seeds 11/13/14 (3 pairs, 136 edges each):
  max |delta| = 0.00211   sd = 0.00047   rms = 0.00047
  => tolerance 3 x max = 0.0063

compare difference_selection = FALSE vs separate bgm():
  group 1: cor 0.9991  max|delta| 0.0289  mean delta +0.00290
           slope 1.0264 [1.0192, 1.0336]  (with intercept 1.0222 [1.0149, 1.0295])
  group 2: cor 0.9984  max|delta| 0.0818  mean delta +0.00580
           slope 1.1003 [1.0904, 1.1101]  (with intercept 1.0949 [1.0845, 1.1053])

decomposition:
  level    ((g1+g2)/2)  slope 1.0643 [1.0597, 1.0689]  cor 0.9997
  contrast ( g1-g2   )  slope 1.0138 [0.9935, 1.0341]  cor 0.9935
```

Verdict: **slope PASS** (09-1) — 1.014 on the contrast, 1.03/1.10 on the
groups, no interval anywhere near 2 or 0.5. **max|Δ| FAIL** against the
0.0063 tolerance (09-2), located in the level, not the scale.

**Pooling effect of `difference_selection = TRUE`**, as the brief asked:

```
contrast slope vs separate fits:  sel=FALSE 1.0138   sel=TRUE 0.4025 [0.3506, 0.4543]
sd of the between-group contrast: sel=FALSE 0.11669  sel=TRUE 0.05769  (factor 0.494)
max |g1 - g2|:                    sel=FALSE 0.4231   sel=TRUE 0.3327
separate-fit contrast (bgm h1 - bgm h2): sd 0.11423, max 0.3879
```

With selection off the compare contrast reproduces the independent-fit
contrast almost exactly (slope 1.014, sd 0.1167 vs 0.1142). Turning selection
on halves it — on null data, where halving is the correct behaviour. Note
also that the *level* comes back to 1.00 under selection
(slope 0.9964 [0.9822, 1.0105]), i.e. the 09-2 level offset is specific to
the unselected fit.

### 3.4 Item 4 — convention-guard audit

Covered in full in 09-4. What the guard pins, restated compactly:

* **Quantity**: `extract_group_params(fit)$pairwise_effects_groups[, 1]` —
  group 1's three reconstructed pairwise associations, posterior means only.
* **Path**: `bgmCompare()` sampler → `.compute_group_param_matrices()` →
  `extract_group_params()`, against `simulate_mrf()`'s ordinal Gibbs kernel
  as the data-generating reference.
* **Setup**: p = 3, 3-category, n = 400 per group, both groups from the same
  `omega = (0.5, 0.0, 0.45)`, `difference_selection = FALSE`, 1 chain,
  `iter = 600 / warmup = 300`, seed 1234.
* **Tolerance**: `rmse(estimate, omega) < 0.2` and
  `rmse(estimate, omega) < 0.5 * rmse(estimate, 2 * omega)`.
* **Tier**: no `skip_on_cran()`, no env gate — it runs on every test
  invocation including CRAN. Measured 2.3 s here (the brief's ~7 s is the
  right order).

### 3.5 Item 5 — Blume-Capel end-to-end

`Wenchuan[, 1:6]` complete cases, `variable_type = c("blume-capel",
"blume-capel", "ordinal" × 4)`, `baseline_category = 2`, random split
(seed 555), `seed = 606`. Stored `baseline_category` `1 1 2 2 2 2`,
`blume_capel_shift` `1 1 NA NA NA NA` (09-6).

**`predict()` against a hand-written reference** — the sampler's own
convention, `exp(main_lin·(c−r) + main_quad·(c−r)² + (c−r)·rest)` for
Blume-Capel and `exp(threshold_c + c·rest)` for ordinal, with
`rest = Σ 2·pairwise·(x_v − centre_v)`, centre = the baseline for Blume-Capel
neighbours and 0 for ordinal ones (`src/mrf_prediction.cpp:99-105`),
reconstructed per group from `extract_group_params()`:

```
  group 1 intrusion  (BC): 2.776e-16      group 2 intrusion  (BC): 3.331e-16
  group 1 dreams     (BC): 3.331e-16      group 2 dreams     (BC): 2.220e-16
  group 1 flash      (ord): 3.331e-16     group 2 flash      (ord): 4.441e-16
  group 1 upset      (ord): 5.551e-16     group 2 upset      (ord): 4.441e-16
  group 1 physior    (ord): 3.331e-16     group 2 physior    (ord): 3.331e-16
  group 1 avoidth    (ord): 2.776e-16     group 2 avoidth    (ord): 2.220e-16

sensitivity of that reference to a scale error on the rest score:
  multiplier 1.0 -> 0.0000    2.0 -> 0.7949    0.5 -> 0.6575
```

Machine precision on both groups and both variable kinds, and the check has
enormous power against the defect class it exists to catch. This extends the
F-068 test (group 1, all-Blume-Capel) to group 2 and to mixed types; no
regression.

**`calibration_check()` per group** runs clean and produces one curve per
variable per group (12 panels, 1212 grid rows). Against a control —
`bgm()` fitted separately to each half with the same variable types — the
compare fit is indistinguishable:

```
                      mean_dev range     share_outside_band range
bgmCompare, per group  [0.0389, 0.0691]   [0.000, 0.168]
bgm() half 1           [0.0334, 0.0652]   [0.000, 0.178]
bgm() half 2           [0.0398, 0.0594]   [0.000, 0.208]
```

**`simulate()` round-trip**, nsim = 20000 per group, comparing the simulated
category margins against the mean predicted category probabilities from
`predict()` on that same simulated data (two different C++ kernels —
`sample_bcomrf_gibbs` and `compute_conditional_probs` — reading the same
reconstructed parameters):

```
  g1: intrusion 0.0034  dreams 0.0083  flash 0.0047  upset 0.0058  physior 0.0033  avoidth 0.0049
  g2: intrusion 0.0039  dreams 0.0054  flash 0.0026  upset 0.0059  physior 0.0019  avoidth 0.0049
```

Max 0.0083 against a binomial Monte-Carlo error of ~0.003 at nsim = 20000 —
consistent within 3 MC errors on every cell. `verdicts()` prints correctly,
excludes the unselected main rows from its counts (the F-060 fix), and
carries both the fragility and scale-contingency caveats.

Verdict: **PASS** on all four surfaces.

### 3.6 Figure

`assets/crosspath_09.pdf` — four panels: item 1's null PIP distributions
across the three splits with the presence threshold marked; the split-101 vs
split-102 PIP scatter showing the absence of cross-split structure; item 3's
compare-vs-separate scatter with the y = x, y = 2x and y = x/2 references;
and the four slope estimates with their intervals against those references.

---

## 4. Open questions

1. **09-3 — now measured, and it needs a maintainer decision.** On
   byte-identical data at identical n, the compare path is ~3.8× noisier on
   true-zero differences and ~4.2× worse in rmse on group 2 than two separate
   `bgm()` fits, with group 2's slope at 1.93 [1.60, 2.25] against `bgm()`'s
   1.03 [0.93, 1.14]. Five alternative explanations were tested and refuted
   (simulation mixing, round-trip fidelity, regime conditioning, truth
   recoverability, units), including one wrong turn of my own that this
   report records rather than hides.

   **What I did not do, and would do next, in order:**
   * **Replicate across data seeds.** Everything above is one simulated
     dataset per cell, and the group-1/group-2 asymmetry flips between the
     n = 400 and n = 2000 cells, so run-to-run variability is clearly large.
     Ten data seeds at n = 2000/group would turn this from an observation
     into an operating characteristic. This is the gate before any fix is
     scoped — I would not act on the single-seed number.
   * **Sweep `difference_scale`.** If the inflation shrinks as the scale
     narrows from 1 toward the realistic ~0.1, the default is the lever; if
     not, the sampler is.
   * **Check whether it reaches real data.** Item 1 says it does not at
     n = 172/group on Wenchuan splits, where difference means had sd
     0.046–0.058. Whether there is an n or an effect-size regime where it
     does is the question that decides release impact.
2. **09-2.** Is the compare path's 6–9 % magnitude excess on the pairwise
   level over matched `bgm()` fits expected from the shared-baseline
   parameterisation? Pooling is refuted as the mechanism; I have no
   replacement.
3. **09-4.** Which of the four coverage holes are worth closing before CRAN?
   My reading of the risk ordering: (i) no every-run test pins a *nonzero*
   difference's scale — this is the parameterisation the breaking change
   touched and it is unguarded; (ii) `simulate.bgmCompare()` has no numeric
   test; (iii) the ×½ blind spot, fixable by tightening `rmse < 0.2` to
   `< 0.15` or adding the symmetric `rmse(target) < 0.5 * rmse(0.5 * target)`
   assertion; (iv) group 2 never asserted, a one-line addition. All four are
   cheap; (i) is the only one that needs a design decision (what nonzero δ,
   at what n, inside a 7-second budget).
4. **09-5.** Should the difference-verdict false-positive tail (log BF 15.4
   on a null split, non-fragile) be characterised properly before release?
   Three splits is not a calibration study. It bears on brief 07's question
   about validating the fragility flag for difference indicators, which the
   same runs would answer.
5. **Documented default divergence** (not a finding, but it bit this brief):
   `bgm()` ships `interaction_prior = normal_prior(1)`, `bgmCompare()` ships
   `cauchy_prior(1)`. Any user comparing a compare fit against separate
   `bgm()` fits at *stated* defaults is comparing different models. Worth one
   sentence in the comparison vignette; it is already on F-002's
   defaults-freeze memo list for `bgm()` but the compare-side asymmetry is
   not called out.

## 5. Artifacts

* Scripts and logs: `~/bgms-review/val09/out/` (`02_item1.R` … `16_apples.R`,
  one `.log` per run, `item1.rds` / `item2.rds` / `item3.rds` /
  `control.rds`). The item-2 diagnostic chain, in the order it was run and
  should be re-read: `04_item2.R` (the four cells) → `11_item2b.R`
  (group-level decomposition) → `14_simcheck.R` (mixing) → `15_roundtrip.R`
  (round-trip fidelity and `bgm()` self-recovery) → `12_control.R`
  (n = 4000/group) → `16_apples.R` (matched n = 2000/group, decisive).
* Export SHA: `1bb1c7ca6ef80606d34a3eb51bff9182101e3ed9`
  (`~/bgms-review/val09/EXPORT_SHA.txt`).
* Figure: `dev/review-2026-08/reports/assets/crosspath_09.pdf`.

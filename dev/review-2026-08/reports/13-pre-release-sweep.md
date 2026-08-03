# Report 13 — pre-release sweep: plot follow-ups, guard closures, CI reds, vignette bake

Branch `fix/pre-release-sweep`, based on `origin/develop` at `3a30c0ea`, merged
forward to `44413b56` and then to `d4f76a48` before this report (see **Merge
freshness**). Eleven commits, all gate items green.

Two headline items for the lead:

* **F-019 is a real defect, not a theoretical one.** The conjugate GGM edge move
  proposes and accepts positive-definiteness violations under prior-only
  sampling; the guard is in with a regression test that bites. It does **not**
  exist at `cran-0.1.6.3` — that version has no GGM model at all — so there is
  no 0.1.6.3-facing NEWS clause. Details in §5.
* **F-103 is stopped and flagged, not fixed.** The diagnosis turned into the
  harm threshold's false-alarm rate, which is a product constant. The block is
  therefore **not** restored to T1. Evidence and the decision needed in §9b.

One more item wants a decision: task 8's tag-vs-rc comparison found that
**0.1.6.3 fails its own convergence diagnostics on three of five synthetic
ordinal datasets** where the release candidate is clean on all five (§8).

---

## What was done

| # | Finding | Outcome | Commit |
|---|---------|---------|--------|
| 1 | F-105 | `display_log_bf()` deleted, both call sites route to `format_log_bf()` | `b0c4b431` |
| 2 | F-107 | `plot_edge_posterior(evidence_threshold =)` removed outright | `11eb07e1` |
| 3 | F-115 | Device-size Rd sections + gated narrow-device advisory | `82a20f56` |
| 4 | F-066 | Slab-frame identity pinned on GGM **and** mixed | `46fffbf8` |
| 5 | F-019 | **Defect confirmed**; guard added, regression test added | `c6a9c40c` |
| 6 | F-042 | Five dead golden-fixture blocks deleted | `f32bdb30` |
| 7 | F-073 | All four convention-guard holes closed | `b5c6d85e`, `+ per-group support fix` |
| 8 | — | Tag-vs-rc estimate agreement (report-only) | *(no commit)* |
| 9a | F-080, F-081 | Not reproducible at 2 cores; one platform dependence removed | `9f3101bb` |
| 9b | F-103 | **Stopped — maintainer decision needed**; block left in T2 | *(no commit)* |
| 9c | F-104 | Worker count settled once + cluster read-back | `030d830a` |
| 10 | — | Prior-sensitivity vignette baked: 42.3 s → 0.1 s | `2e1d1838` |
| 11 | F-077 | **Declined** — provenance runs through an off-limits file | *(no commit)* |

---

### 1. F-105 — `display_log_bf()` removed (minor)

`format_log_bf()` already collapses a Bayes factor that rounds to nothing
(`shown = round(log_bf, 1); if (shown == 0) shown = 0`), so the plot-side
wrapper that pre-zeroed `|log BF| < 0.05` was mapping a value the formatter
already collapsed. Every input reaches the same string: finite values below the
cap round identically, `NA` and `±Inf` bypass the wrapper's `is.finite()` test
and hit the same branches.

**Pixel identity confirmed.** The plot snapshots are the panel-description
snapshots (`tests/testthat/_snaps/plot-methods.md`), which include the rendered
`log BF = …` strings. `test-plot-methods.R` runs 0 failures with **no
snapshot re-record** — `tests/testthat/_snaps/` is byte-unchanged on the whole
branch. Nothing to stop on.

### 2. F-107 — `plot_edge_posterior(evidence_threshold =)` removed (minor)

Removed from the signature, the roxygen, the `check_evidence_threshold()` call,
and `edge_selection_evidence()`, which now calls `verdicts()` at its default and
reads a row whose contents are threshold-free.

**Checked first, as instructed:** nothing else in the tree had begun using the
argument. The only two references were the argument-validation test
(`test-plot-methods.R:46`) and one direct internal call
(`test-plot-methods.R:477`), both updated. `vignettes/checking-your-model.Rmd`
and `dev/review-2026-08/renders/11/make_renders.R` call the function without it.
No NEWS entry — never shipped.

### 3. F-115 — the default-device squeeze, documentation remedy only (minor)

No geometry touched.

(a) A `\strong{Device size.}` section in `?plot.bgms` and `?plot.bgmCompare`
(roxygen source in `R/plot_bgms.R`), naming the concrete suggestion
`width = 13, height = 5` with a copyable `dev.new()` / `pdf()` line, and saying
which displays *don't* need it.

(b) `warn_narrow_device()`, called only from `draw_evidence_panels()` — the
three-panel path — when `par("din")[1] < 10`. It follows the `bgms.verbose`
pattern and names the right topic per fit type, derived from `unit$kind` so no
signature changed. Three tests cover it: the message fires on a 7-inch device
and names `?plot.bgms`; a 13-inch device is silent; `bgms.verbose = FALSE` is
silent; and the single-panel and `type = "groups"` displays are silent.

### 4. F-066 residue — slab-frame identity extended (test)

**Both GGM and mixed were added** — both are cached seconds-scale fixtures, so
mixed was free. Placed with the existing ordinal pin in `test-plot-methods.R`.

The GGM block is the one with content: it asserts `extract_pairwise_interactions()
== anchor_draws()$theta` *and* `== -0.5 * raw$pairwise`, plus
`expect_false(isTRUE(all.equal(extracted, raw)))` so the equality is a claim
rather than a tautology. That is the −0.5 double-application trap: applying the
factor in the extractor and again on the way to the prior breaks the first
equality. Mixed is added because it lays its pairwise draws out by block, so the
agreement has to hold column for column in that layout.

### 5. F-019 — the GGM `n_ == 0` positive-definiteness question (**major — defect confirmed**)

#### (a) Targeted read

**Reachability: yes.** `sample_ggm_prior(spec = "joint"/"hierarchical")` builds
`inputFromR` with `n = 0L, suf_stat = matrix(0, p, p)`
(`R/sample_ggm_prior.R:290-292`), and `GibbsSampler::initialize()`
(`src/mcmc/samplers/gibbs_sampler.h:34`) calls
`model.set_conjugate_edge_proposal(true)`. So `update_method = "gibbs"` reaches
`update_edge_indicator_conjugate` with no data.

**In exact arithmetic the move is PD-safe.** The cofactor constants give
`K_ij = c₂ + c₃·φ` and `K_jj = c₄ + φ²` with `c₄ = K_jj − φ_old²`; that is the
Roverato parameterization, in which only one entry of the Cholesky factor moves
and every diagonal entry of the factor is untouched, so any real `φ` keeps `K`
positive definite. Both branches lie on that curve — the spike the delete
proposes is `constrained_diagonal(0)`.

**But the constants are not exact.** They are read from `covariance_matrix_` and
`log_det_precision_`, which are maintained by SMW rank-1/rank-2 updates within a
sweep and drift. And unlike `ggm_edge_move`/`ggm_diag_move`, this ratio carries
**no determinant term at all** — the tilt and the likelihood determinant cancel
by design (`ggm_model.cpp:996-1001`) — so there is nothing in it to veto a
drifted proposal. With data, `n_ · logdet` → −∞ as the proposal approaches the
cone boundary and the likelihood rejects it; with `n_ == 0` nothing does. This is
the mechanism the header comment on `proposal_is_positive_definite_()`
(`ggm_model.h:720-728`) already documents — it had simply never been applied to
the conjugate move.

#### (b) Empirical

`sample_ggm_prior(spec = "joint", update_method = "gibbs",
edge_inclusion_prob = 0.5)`, 4000 samples after 1000 warmup, seeds 1–60 per
cell. `delta = 0` is the adversarial setting: the determinant tilt is exactly the
term that repels the chain from the cone boundary.

| p | delta | `gibbs` (conjugate, unguarded) | `adaptive-metropolis` (guarded) |
|---|-------|-------------------------------|---------------------------------|
| 8 | 0 | 0/60 | 0/60 |
| 8 | auto | 0/60 | 0/60 |
| 15 | 0 | **2/60** (seeds 31, 35) | 0/60 |
| 15 | auto | 0/60 | 0/60 |
| 25 | 0 | **5/60** (seeds 2, 9, 10, 33, 36) | 0/60 |
| 25 | auto | 0/60 | 0/60 |

Failure mode: `chol(): decomposition failed`, thrown by `refresh_cholesky()`.
A longer probe (20 000 draws, seeds 101–105 / 201–203) reproduced the same throw
at p = 15 and otherwise found no non-PD accepted draw and no non-finite value;
records are in `~/bgms-review/f019/`.

**With data: 0 failures in 300 runs** — `bgm(variable_type = "continuous",
update_method = "gibbs")` at p ∈ {15, 25} × n ∈ {20, 50, 200} × delta ∈ {0, auto}
× 25 seeds. That is why the guard is scoped to `n_ == 0`, mirroring the existing
two exactly rather than guarding unconditionally.

**Attribution.** A temporary instrumented build marked each move before it ran
and evaluated the *proposed* matrix directly. Six of the seven throws are inside
`update_edge_indicator_conjugate`, and in every one the proposal itself was
already outside the cone and was accepted:

```
F019 THROW in=EDGE_CONJ i=1  j=2   pre_min_eig= 5.545026e-06  proposal_min_eig=-1.427075e-04  cond=1.196e+05
F019 THROW in=EDGE_CONJ i=11 j=12  pre_min_eig= 1.193669e-09  proposal_min_eig=-9.784721e-09  cond=1.592e+09
F019 THROW in=EDGE_CONJ i=20 j=22  pre_min_eig= 8.693696e-06  proposal_min_eig=-1.652327e-06  cond=1.029e+07
F019 THROW in=EDGE_CONJ i=1  j=24  pre_min_eig= 2.387359e-06  proposal_min_eig=-2.499444e-06  cond=6.639e+06
F019 THROW in=EDGE_CONJ i=6  j=21  pre_min_eig= 4.071810e-06  proposal_min_eig=-2.378084e-06  cond=8.274e+06
F019 THROW in=EDGE_CONJ i=9  j=23  pre_min_eig= 2.097743e-06  proposal_min_eig=-3.058359e-06  cond=5.554e+06
F019 THROW in=ROW_GIBBS  i=20 j=-1 pre_min_eig= 1.403617e-02  proposal_min_eig=nan  now_min_eig=-5.901838e-10  cond=3.590e+10
```

In each conjugate case it was the chain's *first* non-PD proposal
(`nonpd_proposals_so_far=1`), so the guard costs nothing in acceptance terms —
it refuses exactly the move that breaks the run.

#### (c) The guard

`edge_proposal_is_positive_definite_(i, j, kij, kjj)` fills
`precision_proposal_` from the three changed entries and defers to the existing
`proposal_is_positive_definite_()` — same `arma::chol` test, same tolerance, same
shape. Guarded on both branches under `n_ == 0`. A refused move returns
`accept_prob = 0`, which the Rao-Blackwellized estimator reads as "would not have
added" for a birth and "stays included" for a death — correct, since a non-PD
state carries no prior mass.

**After the guard**, the same 60-seed sweep: 0/60 everywhere except p = 25,
delta = 0, seed 33 — the one `ROW_GIBBS` case, exactly as the attribution
predicted.

**Regression test** (`test-sample-ggm-prior.R`, ~2.4 s, `skip_on_cran()`): the
four reproducing (p, seed) pairs must complete and every retained draw must be
positive definite. Proven to bite against the unguarded build:

```
── 1. Failure ('test-sample-ggm-prior.R:362:5'): the conjugate edge move keeps prior-only chains inside the PD cone
p = 15, seed = 31: chol(): decomposition failed
```

**Cost.** An O(p²) copy and an O(p³) factorization per prior-only edge move —
the same price the existing two guards already pay, and off the data path
entirely. Prior-only Gibbs sampling at `delta = 0`, matched configurations:
p = 8 / 20 000 draws 0.30 s → 0.44 s; p = 15 / 20 000 draws 1.5 s → 2.67 s;
p = 25 / 10 000 draws 3.5 s → 8.83 s. Roughly 1.5× at p = 8 rising to 2.5× at
p = 25. I mirrored the existing guards as instructed rather than substituting a
cheaper rank-2 determinant test; if that cost matters on the SBC reference path,
say so and it can be revisited.

#### (d) Does it exist at `cran-0.1.6.3`?

**No.** `git ls-tree cran-0.1.6.3 -- src/` has no `models/ggm/` at all — the tag
ships only `bgm/`, `bgmCompare/`, `math/` and `mcmc/`. `update_edge_indicator_conjugate`,
`sample_ggm_prior()` and the Gibbs update method are all new in 0.2.0.0.

**No NEWS clause is proposed** (see *Proposed NEWS clauses* below).

#### New finding — the row-block Gibbs draw can lose definiteness too

The seventh throw is a different path: `update_row_block_gibbs()`, an *exact*
conjugate draw with no proposal to reject, taking `min eig` from +1.40e-02 to
−5.90e-10 on a matrix of condition number 3.6e10. That is a rounding-level loss,
not an acceptance error, so a PD guard is the wrong instrument — rejecting an
exact Gibbs draw would break the invariant distribution. The right remedy is
probably for `refresh_cholesky()` to fall back rather than throw when `K` is
numerically singular. **Out of scope for this brief; opening as a new finding.**
Reproducer: `sample_ggm_prior(p = 25, n_samples = 4000, n_warmup = 1000,
seed = 33, spec = "joint", update_method = "gibbs", delta = 0)`.

### 6. F-042 — the five dead golden-fixture blocks deleted (minor)

Deleted outright per the maintainer ruling: no stubs, no skip placeholders. The
two resolver helpers (`golden_fixture_path()`, `has_golden_fixtures()`) had no
other caller and went with them, as did the file header's entry for the group.
`dev/fixtures/` does not exist in the repo, confirmed.

Sections 7 and 8 untouched. Surviving section numbers keep their identities
rather than being renumbered, so the brief's references to "sections 7 and 8"
still resolve. The file now runs **0 failures, 0 skips** (previously five skips
reading "golden fixtures not found"). No NEWS — never functional.

### 7. F-073 — the four compare convention-guard holes (major, test coverage)

**(a) Planted-difference scale pin.** New block in `test-bgmCompare.R`, ~4 s. The
two group matrices *swap* their two nonzero pairs — `(0.60, 0, 0.10)` against
`(0.10, 0, 0.60)` — so the planted difference is large (±0.5) while both groups
stay in a well-identified coupling range. Planting a large difference by
inflating one group instead was tried and rejected: at ω ≈ 1.2 on binary
variables the posterior widens faster than the difference grows, and the
read-back gets noisier than the thing being measured (worst-case rmse 0.58
against a planted 0.7 — measured, not assumed).

Assertions: rmse of the recovered difference < 0.15, worst per-pair error < 0.20,
regression slope of recovered on planted in [0.7, 1.4], and the qualitative
signs. **Tolerance derivation** (stated in the test's comments): over 12
alternative data/fit seeds the rmse reached 0.089, no pair was off by more than
0.125, and the slope stayed within [0.813, 1.167].

**(b) Two-sided.** The old `0.5 * rmse(2 * target)` tests only the doubling
direction, and at that data size the halving direction is *not resolvable by rmse
at all*: `target` and `0.5 * target` are 0.194 apart while the run's own error
reaches 0.19. The scale is now read as the slope of the estimate on the target,
which is bounded on both sides by construction — band [0.6, 1.6], from a measured
slope range of [0.851, 1.369] over 12 seeds, excluding both 0.5 and 2.

**(c) Group 2** is asserted: the guard now loops over both groups, which are
drawn from the same ω and so must both recover it.

**(d) `simulate.bgmCompare()` numeric convention**, in
`test-simulate-predict-regression.R` (the file that owns simulate/predict). Data
simulated from a group, scored back through `predict()` with the **same** group,
must reproduce its own category margins — `E[1{X_v = c}] = E[P(X_v = c | X_-v)]`
holds for any fit that puts the same parameters into both paths. Band: the
simulated margin's binomial Monte-Carlo error `sqrt(π(1−π)/nsim)` at four
standard errors, floored at 0.01 so a near-certain category does not get a
vanishing band from its own vanishing variance; the predicted margin averages
conditional probabilities over the same rows and is far tighter, so it
contributes little. Observed same-group residuals 0.0006–0.0114 against MC
standard errors 0.0022–0.0112. Cross-group residuals are 0.03–0.69, and the test
requires at least a tenfold gap, so the identity cannot be satisfied by a fit
that ignores `group`. Predicted columns are mapped to original values through the
fit's own `category_levels`, so collapse cannot turn this into a test of the
coding.

**Spot-proof.** Injecting a factor of 2 into the projection expansion in
`.compute_group_param_matrices()` (`R/extractor_functions.R`), rebuilding and
running `test-bgmCompare.R`:

```
── 4. Failure ('test-bgmCompare.R:442:3'): bgmCompare recovers a planted group difference at its planted size
Expected `sqrt(mean((recovered - planted)^2))` < 0.15.
Actual comparison: 0.39 >= 0.15
── 6. Failure ('test-bgmCompare.R:447:3'): bgmCompare recovers a planted group difference at its planted size
Expected `slope` < 1.4.
Actual comparison: 1.95 >= 1.40
```

The injection was reverted and the build restored before anything else ran.

**Incidental finding — the old guard was seed-fragile.** Its `0.5` factor is met
by its own seed but by only 8 of 12 alternatives (worst observed ratio 0.69) and
its absolute bound of 0.2 was met at 0.190 in the worst case. The retained
doubling bound is therefore re-founded at 0.8 from the same measured spread. The
guard was passing because it was pinned to one seed, not because it had margin.

### 8. Tag-vs-rc1 ordinal estimate agreement (report-only)

Five seeded synthetic ordinal datasets, p = 6, n = 500, 4 categories, **full
support enforced** (generator seeds that produced any variable with a missing
level were rejected and redrawn — 11 of 16 were rejected, so category collapse is
provably not in play). `cran-0.1.6.3` from `~/bgms-review/lib10-tag`. Matched
specification: the tag's own defaults (Cauchy(2.5) slab, beta-prime(0.5, 0.5)
thresholds, Bernoulli(0.5) edges), asked of the rc explicitly since the rc's
defaults changed; `update_method = "nuts"`, `iter = 2000`, `warmup = 2000`,
`chains = 4`, `cores = 2`, one fit at a time.

**Yardstick** — the rc, dataset 1, three sampler seeds (777 / 881 / 882):

| pair | max \|Δ\| main | max \|Δ\| pairwise | max \|Δ\| PIP |
|---|---|---|---|
| 881 vs 777 | 0.0297 | 0.0015 | 0.0072 |
| 882 vs 777 | 0.0345 | 0.0012 | 0.0058 |
| 882 vs 881 | 0.0200 | 0.0009 | 0.0101 |
| **yardstick (max)** | **0.0345** | **0.0015** | **0.0101** |

Verdict agreement across the three reruns: 15/15 every pair.

**The comparison:**

| dataset | max \|Δ\| main | max \|Δ\| pairwise | max \|Δ\| PIP | verdicts agreeing | tag s | rc s |
|---|---|---|---|---|---|---|
| 1 | 0.6899 | 0.2042 | 0.5067 | 15/15 | 11.3 | 31.6 |
| 2 | 0.0378 | 0.3672 | 0.0331 | 15/15 | 15.7 | 24.4 |
| 3 | 0.7704 | 0.1348 | 0.7816 | 14/15 | 2.6 | 17.4 |
| 4 | 1.2119 | 0.6709 | 0.9425 | 14/15 | 2.6 | 33.0 |
| 5 | 0.1088 | 0.4704 | 0.0848 | 15/15 | 16.6 | 29.7 |

**Agreement does not sit within the yardstick** — the ratios run 1.1× to 443×.
That result cannot be read as an rc regression, because on three of the five
datasets the tag's chains are not a posterior sample at all:

| dataset | 0.1.6.3 max R̂ | 0.1.6.3 min ESS | 0.1.6.3 rmse vs truth | 0.2.0.0 max R̂ | 0.2.0.0 min ESS | 0.2.0.0 rmse vs truth |
|---|---|---|---|---|---|---|
| 1 | **12.259** | **1.7** | 0.0794 | 1.024 | 70.1 | **0.0235** |
| 2 | 1.024 | 194.7 | 0.1589 | 1.004 | 581.8 | **0.0107** |
| 3 | **Inf** | **0.0** | 0.0774 | 1.034 | 117.3 | **0.0288** |
| 4 | **Inf** | **0.0** | 0.3119 | 1.037 | 63.5 | **0.1307** |
| 5 | 1.011 | 381.7 | 0.2432 | 1.007 | 295.6 | **0.0391** |

The rc converges on all five and is closer to the generating truth on all five,
by 1.9× to 14.9×. The tag's fast fits (2.6 s on datasets 3 and 4) are the
degenerate ones. On dataset 4 the tag reports every category threshold of
variable 1 as exactly 0.000 and PIPs spread over 0.377–1.000 with no near-zero,
against the rc's 0.005–1.000.

**Reading.** The release candidate's ordinal estimates are not merely compatible
with 0.1.6.3 — where they differ, the evidence points at 0.1.6.3. Two caveats,
stated rather than buried: (i) five synthetic datasets at one size and one
sampler configuration is a spot check, not a certification; (ii) the PIP
comparison may be confounded by an estimator change (the rc's inclusion
probability may be Rao-Blackwellized where the tag's is a sample proportion),
which the run-to-run yardstick cannot see. Neither caveat touches the R̂/ESS
result, which is about the tag alone.

Artifacts: `~/bgms-review/tagvsrc/` (data, per-fit summaries, scripts).

### 9a. F-080 + F-081 — the 2-core CI reds (major)

**Neither reproduces at `cores = 2` on this machine, and the documented oneTBB
mechanism is ruled out for both rather than merely untested.**

**F-080 (`test-mixed-nuts.R` M.2F).** Run at `cores = 2` and `cores = 15`, the
block's numbers are **bit-identical** — all ten pairwise z values, all ten SD
ratios, all ten CI overlaps agree to the printed precision. So the core count
changes nothing here at all. And every assertion clears by more than an order of
magnitude:

| assertion | limit | observed |
|---|---|---|
| pairwise \|z\| | < 3 | max **0.143** |
| SD ratio | (0.6, 1.7) | **0.955 – 1.042** |
| CI overlap | > 0.7 | min **0.966** |
| main \|z\| | < 3 | max **0.156** |
| `total_divergences` | == 0 | **0** in 18 runs (3 data seeds × 6 sampler seeds) |

Statistically equivalent trajectories cannot move `z` from 0.14 to 3, so **this
is not a tolerance question and no tolerance was touched.** The one assertion
with no margin by construction — the exact-zero divergence count — was zero every
time it was sampled.

**Ruled out:** the core count (bit-identical), tolerance margin (20×), and the
divergence assertion (0/18). **Not ruled out, and untestable from here:** the
runner's BLAS and image. That is where the remaining explanation must live.

**F-081 (`test-sbc-ggm.R` diagonal ranks).** The block fits `chains = 1`, so the
oneTBB multi-chain mechanism **cannot apply to it by construction**; `cores = 2`
and `cores = 15` again agree bit for bit. Across eight independent realizations
of the block's own null (seed offsets 0–7 at `cores = 2`) it passed every time:
no KS failure against a limit of one, and a global chi-squared p between **0.078
and 0.723** against a limit of **0.001**. Moving that statistic to failure means
a 78-fold move, not drift.

**One real platform dependence was found and removed rather than tolerated.** The
block drew its data with `MASS::mvrnorm()`, which decomposes Σ with `eigen()`.
An eigendecomposition is not unique — eigenvector signs, and the ordering of
near-equal eigenvalues, are whatever the platform's LAPACK returns — so *the same
seed produced different data on a different LAPACK*, and every "fixed" SBC
realization in that file was in fact re-randomized per platform. The Cholesky
factor of a positive-definite matrix is unique, so the draw now goes through
`chol()` (all six sites in the file). Re-verified across the same eight offsets:
still no KS failure, chi-squared p between 0.030 and 0.886. The whole file at T2
with `cores = 2`: **0 failures, 0 errors, 0 warnings, 0 skips**.

**Honest labelling:** without the CI logs it cannot be shown that the LAPACK
route is what bit. This is preventive. What it does buy is that a future red in
that file is reproducible from its seed, which the old draw could not promise.

### 9b. F-103 — STOPPED, maintainer decision needed

**This became a harm-threshold policy question, so I stopped, as instructed.**

The block asserts `harm_pred < harm_threshold` for
`biased_evidence_free_fit(2)` — the *negative control* on a cell inside the
surface's validated deployment range. `harm_threshold` is 0.01, a product
constant.

Measuring the fixture's own dispersion (only the seed changes; everything else is
the shipped configuration):

| `n_samples` | mean `harm_pred` | sd | max | max / threshold | exceeds threshold |
|---|---|---|---|---|---|
| 1 200 (shipped) | 0.00483 | 0.00407 | 0.01253 | 1.25 | **1 / 12** |
| 4 800 | 0.00779 | 0.00606 | 0.01449 | 1.45 | **6 / 12** |
| 12 000 | 0.00318 | 0.00427 | 0.01298 | 1.30 | **2 / 12** |

At the shipped size the twelve seeds give 0.00036, 0.00076, 0.00095, 0.00111,
0.00203, 0.00336, 0.00435, 0.00662, **0.00794 (the fixture's seed 7)**, 0.00876,
0.00921, **0.01253**. The Linux value of 0.01071 is an unremarkable draw from
that distribution — indistinguishable from what seeds 2, 8 and 12 produce here.

So **this is not platform drift.** The fixture was chosen on one realization that
happened to land at 79 % of the threshold with no margin, and its statistic's own
spread straddles the threshold. Raising the sample size does **not** shrink the
spread, so it is not simply Monte-Carlo-limited in `n_samples`.

**Why that is a product question, not a test question.** Across 36 runs
(12 seeds × 3 sample sizes) the deployed configuration fires the harm alarm on
**9 of 36** healthy fits at a cell the surface is supposed to cover — a ~25 %
false-alarm rate on the negative control. Either the threshold is set below the
statistic's own dispersion, or the gauge's estimator is noisier than the
threshold assumes. Both are decisions about the shipped diagnostic's behaviour,
not about a test constant.

**Consequences for the tier contract, stated plainly:** I have **not** restored
the block to T1. Restoring it would put a roughly one-in-four flaky assertion on
every push, which is worse than the contract violation it would fix. The block
stays in T2 with `skip_unless_certification()` until the threshold question is
settled.

**A lever for whoever settles it.** `harm_pred = |mean(m_e s_e)| · A` is computed
from the gauge's audit stream, and `sample_ggm_prior()` hardwires
`gauge_sweeps = 2L` (`R/sample_ggm_prior.R:330`) while the deployed path reads
`options(bgms.zratio_gauge_sweeps)` via `zratio_gauge_sweeps()`. Raising the
audit sweeps is the parameter that shrinks the *gauge's* own error, as distinct
from the chain's — which is what `n_samples` moves and what the table above shows
does not help. Making the prior sampler honour the same option would let this be
measured without another code change. I did not make that change, because it
touches the shipped diagnostic and the answer belongs to the maintainer.

### 9c. F-104 — the PSOCK worker count (major)

**The original error did not reproduce.** The hardening is therefore preventive,
not proven causal — stated plainly, as the brief allows.

Reading the construction path end to end, three places let the number that
decides differ from the number that was passed:

1. **`zratio_surface_build_cores()` defaulted its argument to `1L`.** A caller
   that said nothing got `getOption("bgms.zratio_surface_cores", 1L)` — the
   *option*, not the call, decided the width. `sample_ggm_prior.R:337` was
   exactly such a caller. This is the clearest "`cores` did not reach the
   constructor" route: nothing the caller intended was ever in play. The default
   is removed; all three call sites now state what they mean.
2. **`normalize_parallel_cores()` ended in `min(cores, detectCores())`**, and
   `detectCores()` is documented to return `NA` when it cannot tell. `min(k, NA)`
   is `NA`, which then travels to a cluster constructor as neither a width nor an
   error. `NA`, empty, non-numeric and non-positive inputs now all resolve to 1.
3. **`zratio_build_surfaces()` read `cores` raw**, so the branch choice, the
   verbose announcement and the constructor each re-derived it. It now normalizes
   once at entry and everything below reads that one value.

And the read-back the brief asked for: after `makePSOCKcluster()` the builder
compares `length(cl)` with the width it asked for and stops with both numbers
named. `makePSOCKcluster`'s own default is `getOption("mc.cores", 2L)`, so a
width that fails to arrive is otherwise entirely silent — the build simply runs
at a width nobody chose, and under `R CMD check` possibly wider than the two
workers CRAN policy allows.

Tests added: the resolver's degenerate inputs including a mocked `NA` from
`detectCores()`; and the read-back firing on a stand-in cluster of the wrong
length (`"asked for 2 PSOCK worker(s) and got 1"`). **Brief 12's visible-skip
probe is untouched.**

### 10. Prior-sensitivity vignette baked

**(a) How the output got in:** a **live chunk** (`vignettes/prior-sensitivity.Rmd`
lines 62–68) that ran `bgm()` and `prior_sensitivity_check()` on every build —
not pasted output. F-086's fix made the report real by making it live.

**(b) The bake.** `dev/vignette-data/make-prior-sensitivity.R` runs the analysis
once, locally, with the seeds and arguments the vignette displays. The vignette
now shows that call `eval = FALSE` and then, in an `echo = FALSE` chunk, reads
the `.rds` and **prints the object**. Verified in the knitted output: the shown
call is un-evaluated, and the report beneath it begins
`#> Prior sensitivity check: are the edge verdicts robust to the slab scale?` —
`print()` on a genuine `bgms_prior_sensitivity` object. F-086's honesty survives:
nothing is transcribed.

**(c) Rails.**

* **Payload:** `vignettes/prior-sensitivity-ps.rds` is **9.2 KB** and holds the
  `ps` object and nothing else. **No thinning was needed** —
  `prior_sensitivity_check()` keeps refits only under `keep_fits`, so no chains
  were ever attached. The generator sets `ps$fits = NULL` as belt-and-braces and
  **refuses to write a file that reaches three digits of KB**, so the rail is
  enforced rather than merely met once.
* **Generator tracked, not shipped:** `.gitignore` gains
  `!dev/vignette-data/` (with a comment saying why, following the existing
  `!dev/validation/` pattern); `^dev$` in `.Rbuildignore` already keeps it out of
  the tarball.
* **Build time: 42.3 s → 0.1 s** (`knitr::knit()` on the vignette, which isolates
  the chunk-evaluation cost that changed; pandoc is not on this machine's
  non-interactive PATH, so the full `render()` was not the instrument).
* **`R CMD check --as-cran`:** see *Verification gate* below.

### 11. F-077 — DECLINED, recommendation instead

**Not done, and the brief's own condition is why.** The write of
`arguments$baseline_category` is `R/build_arguments.R:141`, which is in scope —
but the value it writes comes from `spec$variables$baseline_category`, built in
`R/build_spec.R` from output produced in **`R/validate_data.R`** (lines 358, 428),
which the release-gate batch owns. Changing what is stored without touching that
provenance would mean two disagreeing copies of the same field.

The C++ read sites make it worse rather than better. The finding names
`src/mrf_prediction.cpp:99-105` as the only reader, and that one is indeed
BC-only — but `baseline_category` is also converted to `arma::ivec` wholesale in
`src/sample_mixed.cpp:70`, `src/mixed_gradient_interface.cpp:16`,
`src/mrf_simulation.cpp:50` and `src/bgmCompare_interface.cpp:109/149/373`. An
`NA_integer_` becomes `INT_MIN` in an `arma::ivec`, and none of those sites has
any NA handling: a path that touched it outside a Blume-Capel branch would
produce silent garbage rather than an error. Two of those files are under
`src/models/bgmCompare/`'s owner.

**Recommendation for a later batch:** store `NA` for non-Blume-Capel entries at
the point the field is *built* (`R/validate_data.R` / `R/build_spec.R`, one
place), not at the `arguments` layer; and before doing so, add an explicit
`is_finite` / BC-branch assertion at each `arma::ivec` read site so a future
unconditional read fails loudly instead of computing on `INT_MIN`. Doing this the
week of release mechanics, across two batches' files, is not worth the trap it
would close — the field is inert today.

---

## Verification gate

All six items pass. Every number below is from the merged branch.

### 1. Full local suite, both tiers — PASS

| tier | files | tests | pass | **fail** | **error** | **warning** | skip | elapsed |
|---|---|---|---|---|---|---|---|---|
| default (`NOT_CRAN=true`) | 78 | 1214 | 8445 | **0** | **0** | **0** | 97 | 204.6 s |
| `BGMS_RUN_SLOW_TESTS=true` | 78 | 1214 | 8649 | **0** | **0** | **0** | 61 | 322.3 s |

A first default-tier run showed **1 warning**, and it was fixed rather than
tolerated: the new `simulate.bgmCompare()` block's data had a category observed
in one group but not the other, so the union-category machinery merged from
develop warned about a prior-set threshold — true, and nothing to do with what
the block tests. New data seeds give every group every category, and the block
now **asserts** that before fitting rather than suppressing the warning
(`test-simulate-predict-regression.R`). Both tiers rerun clean afterwards.

**Every expectation I changed, and why:**

| where | change | reason |
|---|---|---|
| `test-plot-methods.R` | dropped `expect_error(plot_edge_posterior(…, evidence_threshold = 0.5), "greater than 1")` | the argument no longer exists (F-107, task 2) |
| `test-plot-methods.R` | `edge_selection_evidence(fit, label, 10)` → `(fit, label)` | signature change (F-107) |
| `test-simulate-predict-regression.R` | five golden-fixture blocks and their two helpers deleted | maintainer-ratified (F-042, task 6) |
| `test-bgmCompare.R` | `expect_lt(rmse(target), 0.5 * rmse(2 * target))` → `0.8 * …`, **on both groups** | **this is the one loosening.** The 0.5 factor was met by the block's own seed but by only 8 of 12 alternatives (worst ratio 0.69), so it was passing by seed rather than by margin. It is re-founded at 0.8 from the measured spread, and more than offset by the two-sided slope band and the new planted-difference block, both of which the old guard had nothing corresponding to. |
| `test-sbc-ggm.R` | `MASS::mvrnorm()` → `rmvnorm_chol()` at all six sites | data draw only — **no assertion, tolerance or limit changed** (F-081, task 9a) |

No other expectation was weakened, re-founded, or removed.

### 2. `R CMD check --as-cran` on a `git archive` tarball — PASS

`git archive HEAD` → `R CMD build` (vignettes built) → `R CMD check --as-cran`
under `_R_CHECK_LIMIT_CORES_=TRUE`:

```
Status: 2 NOTEs
```

Both are the baseline pair, unchanged in text from `check-rc1.log`:

* `checking CRAN incoming feasibility ... NOTE` — "The Date field is over a
  month old."
* `checking HTML version of manual ... NOTE` — "'tidy' doesn't look like recent
  enough HTML Tidy."

The vignette task made this non-optional, and it is clean:
`checking re-building of vignette outputs ... [105s/57s] OK` (that 105 s is all
five vignettes; prior-sensitivity alone was 42.3 s of it before the bake),
`checking package vignettes ... OK`, `checking files in 'vignettes' ... OK`.
Tarball contents confirm the rails: `bgms/vignettes/prior-sensitivity-ps.rds`
ships, and `dev/vignette-data/` does not appear at all.

### 3. Plot snapshots — PASS, no re-records

`git diff 3a30c0ea..HEAD -- tests/testthat/_snaps/` is **empty** over the whole
branch. Task 1 is pixel-identical as required, and tasks 2 and 3 touch no drawn
pixel. Nothing to stop on.

### 4. `devtools::document()` — PASS

Run on the merged tree; `git status` is clean afterwards, so the committed `man/`
is current. `git diff 3a30c0ea..HEAD -- NAMESPACE` is **empty** — NAMESPACE
unchanged.

### 5. The task-7 planted-δ test proven to bite — PASS

Quoted in §7 above:
`Expected sqrt(mean((recovered - planted)^2)) < 0.15. Actual comparison: 0.39 >= 0.15`,
with the slope at 1.95 against its 1.4 bound. Injection reverted and the build
restored before anything else ran.

### 6. Merge freshness — PASS, and it happened twice

`origin/develop` moved **twice** while this batch was in flight.

* First merge: `3a30c0ea` → **`44413b56`** (18 commits), which brought the
  category-collapse and compare-gradient work — including `test-bgmCompare.R`
  and `test-simulate-predict-regression.R`, both of which this batch also edits.
  **No conflicts.** The package was rebuilt and re-documented after it, and every
  test file this batch touched was re-run: `test-plot-methods.R`,
  `test-bgmCompare.R`, `test-simulate-predict-regression.R`,
  `test-sample-ggm-prior.R`, `test-sbc-ggm.R` (full file at T2, `cores = 2`),
  `test-zratio-surface-build.R` — all via the two full-suite runs above, which
  postdate the merge.
* Second merge: `44413b56` → **`d4f76a48`**, a two-line docs-only advance to
  `FINDINGS.md` and `PLAN.md`. No code, so no re-run was warranted.

The first merge is also what surfaced the suite warning in item 1 — the
union-category semantics did not exist when the new block was written.

---

## Proposed NEWS clauses

**None.**

Task 5(d) asked for a Bug-fixes clause aimed at the 0.1.6.3 reader if the F-019
defect exists at the tag. It does not: `cran-0.1.6.3` has no
`src/models/ggm/` and no `sample_ggm_prior()`, so the affected code path is new
in 0.2.0.0 and has never shipped. A 0.1.6.3 reader has nothing to be warned
about.

If the lead nonetheless wants the guard recorded for a 0.2.0.0 reader — it is a
correctness fix in an unreleased path, so this is a judgement call, not a
requirement — the clause would be, verbatim:

> * The conjugate edge update used by `update_method = "gibbs"` could accept a
>   proposal outside the positive-definite cone when sampling with no data
>   (`sample_ggm_prior()`), which aborted the chain with
>   `"chol(): decomposition failed"`. Such proposals are now rejected, as they
>   already were on the Metropolis update paths.

No NEWS for F-042 (never functional) or F-107 (never shipped), per the brief.

---

## Open questions

1. **F-103 (blocking a tier-contract item).** Is the harm threshold of 0.01
   right, given that the negative control's own statistic exceeds it on 9 of 36
   healthy runs? Until that is answered the block cannot go back to T1. The
   audit-sweep lever in §9b is the cheapest way to find out whether the
   estimator or the constant is the problem.
2. **The row-block Gibbs definiteness loss (new finding, §5).** Should
   `refresh_cholesky()` fall back instead of throwing when `K` is numerically
   singular? Reproducer given.
3. **F-019's guard cost.** 2.5× on prior-only Gibbs at p = 25. I mirrored the
   existing guards as instructed; a rank-2 determinant test would be far cheaper
   and equally sound. Worth it, or not?
4. **Task 8's convergence result.** 0.1.6.3 fails R̂/ESS on three of five modest
   synthetic ordinal datasets. Is that already known? If not, it is a stronger
   argument for the release than any agreement table would have been, and it may
   deserve a line in the release notes.
5. **F-077.** Confirm the recommendation in §11 is the shape a later batch should
   take, or overrule it.
6. **F-080 / F-081.** Both are closed here as "not reproducible on macOS at
   2 cores, with the core-count mechanism ruled out". If the CI logs still exist,
   the specific failing assertion for M.2F would settle whether the divergence
   count or an agreement bound was the one that fired — I could not tell from
   here, and the margins say it should have been neither.

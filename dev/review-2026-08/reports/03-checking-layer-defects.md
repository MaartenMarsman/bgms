# Report 03 — Checking-layer defect batch: reproduce, root-cause, fix

Agent: brief 03 (Opus). Date: 2026-08-01.
Branch: `fix/checking-layer-batch`, 8 commits off `develop` (`abe0312c`), local only.
Worktree: `~/bgms-review/wt-fix` (outside Dropbox). The Dropbox checkout was
not touched and stayed on `main`.

**Headline.** All nine findings reproduced and fixed. Two of them turned out to
be worse than the brief assumed, and both are **correctness** defects rather
than the UX defects they were filed as:

- **F-037 is a blocker, not a diagnosis task.** On a mixed fit, `verdicts()`
  labelled its rows from one permutation of the edges and filled its numbers
  from another. On MM's 17-variable Wenchuan fit **all 136 rows were
  mismatched**, and the edge MM flagged — intrusion–dreams, true log BF
  **+122.5** — was reported as **evidence of absence at log BF −1.33**. The
  same misalignment put the network plot's edge weights on the wrong pairs.
  Unambiguous mechanical bug; fixed with a test.
- **F-036 is a wrong number, not just a slow one.** For
  `precision_graph_prior = "hierarchical"`, `extract_prior_inclusion_probabilities()`
  applied the *joint* specification's determinant tilt. On MM's fit it returned
  a prior inclusion probability of **0.4219 where the specification says 0.5**,
  shifting **every** inclusion Bayes factor by **+0.315 nats**. The 2-minute
  wait was the package building a correction table it had no use for.

The brief's mechanism hypothesis for F-036 was **incorrect** in a way that
matters: no prior-only chain runs. See §2.2.

The interrupt question has a clean answer: an interrupt during that route is
**swallowed, not partially applied**. The call runs to completion and returns
the correct table. There is no path to partial prior PIPs. See §2.3.

Three findings are new; one of them is a **crash on a documented
configuration** (§4.1).

---

## 1. What was done, per finding

| Finding | Commit | Reproduced | Hypothesis |
|---|---|---|---|
| F-035 | `0131bcdc` | yes | **confirmed** |
| F-038 + F-039 | `20b56996` | yes | n/a (decided) |
| F-037 | `3b42a67c` | yes | **corrected** — mechanical index bug, blocker |
| F-018 | `a6c0725d` | yes | **confirmed** |
| F-036 | `ff29d177` | yes | **corrected** — no chain involved; wrong quantity |
| F-041 | `0972109a` | yes | **confirmed** |
| F-040 | `aff555f7` | yes | partly pre-existing (`variables =` already shipped) |
| F-047 | `1aa527b4` | yes | n/a |

21 files, +788/−100. Every fix has a regression test that fails on `develop`
and passes on the branch.

---

## 2. The three substantive findings

### 2.1 F-035 — Blume-Capel category bookkeeping (major, confirmed)

**Reproduced** exactly as reported, on a pure Blume-Capel fit and on a mixed
fit with Blume-Capel columns:

```
Error in x[, j] <- levels_list[[j]][x[, j] + 1L] :
  number of items to replace is not a multiple of replacement length
```

**Hypothesis confirmed, with one correction.** `reformat_ordinal_data()`
(`R/validate_data.R:246-322`) records a recode map in `category_levels` only
for *regular ordinal* columns; a Blume-Capel column records an additive shift
in `blume_capel_shift` instead and leaves `category_levels[[j]]` **NULL**.
`NULL[idx]` has length 0, which is the error. Measured on the repro:

```
pure BC:   category_levels = list(NULL, NULL, NULL, NULL, NULL, NULL)
           blume_capel_shift = 1 1 1 1 1 1
mixed:     category_levels = list(1:5, 1:5, NULL, NULL)
           blume_capel_shift = NA NA 1 1
```

The brief also suspected mixed **name** mis-alignment. That part is **not** a
defect: on the mixed path both vectors are indexed over the discrete columns
alone, and `data_columnnames_discrete` has the matching length. The names line
up; only the NULL entries break.

There was a **second, silent** occurrence of the same mistake that the crash
masked. `calibration_check.bgms()` looked the observed category up with
`match(newdata[, v], levels_list[[v]])`; for a Blume-Capel column that is
`match(x, NULL)` → all `NA`, i.e. a silently wrong calibration curve rather
than an error. Both sites are fixed.

**Fix** (`R/calibration_check.R`). `fitted_observed_data()` now decodes through
`recode_simulated_to_original()` — the inverse `simulate()` already uses, which
handles both bookkeeping forms — and a new internal `discrete_category_index()`
reads the shift where there is no recode map, erroring explicitly if a fit
carries neither.

**Tests** (`test-calibration-check.R`): pure BC fit, mixed fit with BC columns,
and an ordinal fit with non-contiguous scores (1, 3, 5, 7, 9 …) — the third
case the brief asked for; it is supported and now covered. Each asserts that
`fitted_observed_data()` round-trips to the input data.

Before/after: §3.1.

### 2.2 F-036 — hierarchical fits were served the wrong prior (major; mechanism corrected)

**The brief's mechanism is wrong.** `prior_only_chain_pips()` is never reached
on this fit. With the package defaults (`edge_prior = bernoulli_prior(0.5)`,
`interaction_prior = normal_prior(scale = 1)`) the route is:

```
verdicts() → extract_inclusion_bf() → extract_prior_inclusion_probabilities()
           → prior_pip_class_values()  [has_joint_block = TRUE]
           → prior_pip_table()         → ggm_correction_table()
           → build_ggm_correction_table()   ← the minutes
```

That is the **tilted prior sweep** for the joint specification, not a
prior-only chain. Measured cold, at q = 17: **68.1 s and 77.9 s** on two clean
runs (`~/bgms-review/repro/r06-clean-build.R`). MM's own cache file
`ggm_ctable_v1_q17_delta1.4166067_...rds` is dated **Aug 1 14:32**, which is
that build.

**The wait was not the real defect. The number was.** The `bgm()` documentation
already states the mathematics (`R/bgm.R:190`):

> `"hierarchical"` — the hierarchical specification \(p(\Gamma)\,p(K\mid\Gamma)\)
> with \(p(K\mid\Gamma)\) normalized per graph, **so the graph marginal is
> exactly the edge prior** \(\pi(\Gamma)\).

and PR #193's own test asserts it (`test-joint-realized-prior-notice.R`):

```r
# Hierarchical targets the nominal edge prior, so there is nothing to report.
expect_no_message(notice(spec = "hierarchical"))
```

**The mathematical claim, explicitly.** The prior inclusion Bayes factor needs
\(\Pr(\gamma_{ij}=1)\) under the fit's own prior, i.e. the graph marginal of
\(p(K,\Gamma)\) after integrating \(K\) out.

- *Joint* (`precision_graph_prior = "joint"`): \(p(K,\Gamma) \propto
  \rho_\Gamma(K)\,\pi(\Gamma)\) is un-normalized, so
  \(p(\Gamma) \propto \pi(\Gamma)\,Z(\Gamma)\) with
  \(Z(\Gamma)=\int\rho_\Gamma(K)\,dK\). The tilt does not cancel and the
  realized edge density must be measured — this is what the correction table's
  `edens` curve is, and it is correct there.
- *Hierarchical*: \(p(K\mid\Gamma) = \rho_\Gamma(K)/Z(\Gamma)\) integrates to 1
  **for every** \(\Gamma\), so
  \(p(\Gamma) = \pi(\Gamma)\int p(K\mid\Gamma)\,dK = \pi(\Gamma)\) exactly.
  The prior inclusion probability is the edge prior's own marginal in **every**
  edge class, continuous–continuous included: \(\theta_0\) for
  `bernoulli_prior()`, \(\alpha/(\alpha+\beta)\) for `beta_bernoulli_prior()`,
  and the exchangeable partition mixture for `sbm_prior()`. This is exact and
  closed-form; there is nothing to tabulate and no chain to run.

  *Caveat, stated for the record:* the deployed hierarchical sampler evaluates
  \(Z(\Gamma')/Z(\Gamma)\) with the Option-B surface, so the **realized**
  chain sits within the surface's approximation error of that target
  (~0.003 nats in range; audited by `summarize_zratio_gauge()`). The closed
  form is the prior of the *specified* model — the same standard the joint
  route's table meets, which is itself a Monte-Carlo estimate of its own tilt.

**Measured impact of the old behaviour** on MM's fit (17 continuous variables,
Bernoulli(0.5), δ = 1.4167):

| | prior PIP | prior log-odds |
|---|---|---|
| specification (hierarchical) | **0.5** | 0.000 |
| what 0.2.0.0 returned | 0.4219 | −0.3148 |

Every inclusion Bayes factor on a hierarchical fit was inflated by
**exp(0.3148) = 1.37×**, i.e. **+0.315 nats** (+0.137 on the log10 scale) —
a uniform shift of every edge toward "presence", biggest where it matters
most, at the boundary.

**Fix** (`R/extract_prior_inclusion_probabilities.R`). New
`prior_graph_is_tilted()`; `has_joint_block` is now `num_cont >= 2 &&
prior_graph_is_tilted(spec)`, and the heterogeneous-Bernoulli pass-through
follows the same predicate. A hierarchical fit returns the closed form
instantly, builds no table and runs no chain. Where a table *is* genuinely
needed (joint spec), `prior_pip_table()` now passes `verbose` and
`show_progress = TRUE`, so the build announces itself and draws the sampler's
own progress bar rather than pausing in silence — the fallback the brief asked
for, kept for the case that still needs it.

**Timing, after:** 0.001 s (was 68 s cold). `elapsed < 5 s` is asserted in the
test, run against an empty cache directory.

> **Gate.** The brief asked for MM's confirmation of the mathematical claim
> before switching the default route. The switch is implemented as its own
> commit, `ff29d177`, and the branch is not pushed. If MM rejects the claim,
> dropping that one commit reverts the route and leaves the other seven fixes
> intact. My reading is that the claim is not really new: it is the
> specification `bgm()` documents and PR #193 already tests for.

**Caching (brief item 3).** Answered empirically
(`~/bgms-review/repro/r03-f036.R`):

- Cached on the fit's diagnostics environment as
  `cache$prior_inclusion_class_values` (`R/extract_prior_inclusion_probabilities.R:499-509`).
  Confirmed absent before the first call, present after.
- Second call is served from it (0.000 s). This is what MM's "then it does
  show results!" is, once the disk table exists.
- It **survives `saveRDS`/`readRDS`** — verified, and now asserted in a test.
  (Environments serialize.)
- **An interrupt cannot poison it.** Nothing is written until the value is
  computed; the assignment is the last statement.
- The correction table has a *second*, disk-level cache at
  `tools::R_user_dir("bgms", "cache")`, which raises a separate concern —
  see §4.2.

### 2.3 The interrupt-correctness question (F-036 item 1) — answered

**No, an interrupted call cannot yield partial prior PIPs.** The interrupt is
ignored outright.

Experiment (`~/bgms-review/repro/r04-interrupt.R`): a q = 17 hierarchical fit,
an emptied cache directory, `verdicts()` wrapped in
`tryCatch(interrupt = ...)`, and `SIGINT` sent from outside 20 s into a build
that takes 69 s.

```
t0= 2026-08-01 15:52:06
READY                       <- SIGINT sent at ~15:52:26
t1= 2026-08-01 15:53:15     <- 69 s later
=== RESULT CLASS: bgms_verdicts/data.frame ===
!!! verdicts() RETURNED A TABLE AFTER THE INTERRUPT !!!
cache files left behind: 2
```

No `interrupt` condition was raised, no error, no truncation. And the table it
produced is **bit-for-bit identical** to a clean build:

```
clean A vs clean B    max|raw edens diff|: 0
clean A vs interrupt  max|raw edens diff|: 0
```

(Two clean runs of the identical script also agree exactly, so the sweep is
deterministic and the comparison is meaningful.)

**Mechanism.** `Rcpp::checkUserInterrupt()` appears only in
`src/mrf_simulation.cpp`. No sampler translation unit calls it, so a pending
`SIGINT` is not acted on until the C++ chain returns — and by then the sweep
cell has completed and the loop moves on. So: **no correctness risk, but
Ctrl+C genuinely does not stop a bgms computation.** Logged as a new finding
(§4.3) rather than fixed here — adding interrupt checks to the sampler loop is
a C++ change well outside this brief.

The upshot for the release: no silently wrong Bayes factors came out of MM's
interrupted call. What MM saw was the correct table, arriving late.

### 2.4 F-037 — mixed-fit verdicts were attributed to the wrong edges (BLOCKER)

**Verdict: unambiguous mechanical bug.** No statistical judgment involved, so
it is fixed here per the brief's protocol.

**Step 1 — controlled comparison.** Same 17 Wenchuan columns, chains = 4,
seed = 1, defaults otherwise. Edge 1–2 (intrusion–dreams), from the
**matrices** (which were always right):

| fit | log BF | posterior PIP | prior PIP |
|---|---|---|---|
| `variable_type = "continuous"` (GGM, joint) | 65.45 | 1.000 | 0.4219 |
| 16 continuous + 1 ordinal (mixed) | 87.79 | 1.000 | 0.4200 |
| MM's mixed (6 cont / 6 ord / 5 BC) | 122.47 | 1.000 | 0.3952 |
| all ordinal | 297.03 | 1.000 | 0.5000 |
| all Blume-Capel | 304.31 | 1.000 | 0.5000 |

Shared continuous block, GGM vs the 16-continuous mixed fit, all 120
continuous–continuous pairs: **median |Δ log BF| = 0.103**, max |Δ PIP| = 0.141,
cor(PIP) = 0.997. The blocks agree; the modelling is fine.

**Step 2 — where it diverges.** `verdicts.bgms()` built its index as the
row-major upper triangle of the *variable* order:

```r
idx = which(upper.tri(matrix(0, num_variables, num_variables)), arr.ind = TRUE)
idx = idx[order(idx[, "row"], idx[, "col"]), , drop = FALSE]
log_bf = extract_inclusion_bf(bgms_object, log = TRUE)[idx]
pip    = bgms_object$posterior_mean_indicator[idx]
build_verdicts(parameter = raw$parameter_names$indicator, ...)
```

but a **mixed fit lays its indicators out by block** — discrete–discrete, then
continuous–continuous, then cross (`fill_mixed_symmetric()`,
`build_output_mixed_mrf.R:120-152`). Measured on MM's fit:

```
row-major upper-tri labels[1:3]:  intrusion-dreams  intrusion-flash  intrusion-upset
raw parameter_names$indicator:    avoidact-amnesia  avoidact-lossint avoidact-distant
IDENTICAL ORDER? FALSE
n mismatched positions: 136 of 136
```

The `parameter` column came from one permutation, every numeric column from
the other. Consequence:

```
verdicts(mm_mixed)["intrusion-dreams", ]
  pip 0.0443   bf 0.0463   log10_bf -1.334   verdict ABSENCE
true value from the matrix: log BF +122.47, RB pip 1.000
```

That is exactly MM's "very small … must be a bug". It is.

Everything else raw-side is already in the block layout — the pairwise draws,
`posterior_summary_indicator` rows, and `parameter_names$indicator` all agree
(verified). `verdicts()` was the single point of failure. A GGM or ordinal fit
is unaffected (the two orders coincide; verified).

**Second casualty.** `plot.bgms(type = "network")` builds `pairs` as the
row-major upper triangle but takes `weight = colMeans(extract_pairwise_interactions(x))`,
which is in the **block** layout — so on a mixed fit every drawn edge carried
another edge's weight. The fix repairs both.

Also mis-attributed, and now repaired as a side effect: `mcse` and `draws`
enter `build_verdicts()` in block order, so the **fragility flag** was computed
from another edge's standard errors on mixed fits.

**Fix** (`R/verdicts.R`, `R/plot_bgms.R`). New internal
`indicator_pair_index(fit, num_variables)` reads the positions off the fit's
own layout (by filling `fill_mixed_symmetric()` with draw positions and reading
them back). `verdicts.bgms()` and `plot.bgms()` both use it. The documented
contract — "one row per indicator, in the order of the fit's raw indicator
draws" — is preserved; it is now actually honoured.

**Tests**: an interleaved-type mixed fit where the two layouts genuinely differ,
asserting every reported `log_bf`/`pip` equals the matrix cell its own row
names; plus a single-type fit asserting the row-major layout is unchanged; plus
a mixed-network test asserting the weights line up with the pairs they are
drawn at.

Before/after: §3.3.

---

## 3. Evidence — before and after

### 3.1 F-035

Before (develop, `~/bgms-review/repro/r01-f035.R`):

```
=== (a) pure Blume-Capel fit ===
category_levels: [[1]] NULL ... [[6]] NULL
blume_capel_shift: 1 1 1 1 1 1
[1] "number of items to replace is not a multiple of replacement length"

=== (b) mixed fit with BC columns ===
category_levels: [[1]] 1 2 3 4 5  [[2]] 1 2 3 4 5  [[3]] NULL  [[4]] NULL
blume_capel_shift: NA NA 1 1
[1] "number of items to replace is not a multiple of replacement length"
```

After (branch):

```
--- bc6 ---
Calibration of the conditional predictions, 95% consistency band from 20 resamples:
  variable kind mean_dev max_dev share_outside_band
 intrusion  pav    0.055   0.160              0.297
    dreams  pav    0.050   0.151              0.188
 ...
--- mixed6 ---
  variable kind mean_dev max_dev share_outside_band
     flash  pav    0.038   0.140              0.089
   physior  pav    0.030   0.119              0.050
     upset  pav    0.037   0.112              0.010
   avoidth  pav    0.037   0.108              0.158
    dreams  pit    0.035   0.103              0.337
 intrusion  pit    0.025   0.067              0.297
--- ord6_noncontig ---   (scores 1,3,5,7,9)
  variable kind mean_dev max_dev share_outside_band
    dreams  pav    0.046   0.166              0.208
 ...
```

### 3.2 F-036

```
BEFORE   hier prior pip: 0.4219387   (68 s cold, silent)
         joint prior pip: 0.4219387

AFTER    hier prior pip: 0.5         (0.001 s, no table built)
         joint prior pip: 0.4219387  (unchanged, announced when built)
```

### 3.3 F-037 / F-038 / F-039

Before, MM's mixed fit:

```
        parameter   pip   log10_bf   verdict
 intrusion-dreams 0.044     -1.334   absence      <- true log BF +122.5
```

After:

```
Edge verdicts at an inclusion Bayes factor of 10 (and 0.1 for absence):
presence: log BF > 2.30; absence: log BF < -2.30

  presence 21 | undecided 56 | absence 59   (136 indicators)

        parameter   pip log_bf   verdict fragile
 avoidact-amnesia 0.653  0.634 undecided   FALSE
 avoidact-lossint 0.829  1.578 undecided   FALSE
 avoidact-distant 0.774  1.234 undecided   FALSE
    avoidact-numb 0.058 -2.794   absence   FALSE
 ...

7 verdicts are Monte-Carlo fragile: ...
```

and `intrusion-dreams` now reads `pip 1.000, log_bf 122.47, presence`.

Panel titles (F-039), the string MM saw replaced:

```
mm_mixed   intrusion-dreams | presence, log BF = 122.5
ord_all    intrusion-dreams | presence, log BF = 297.0
ggm_joint  intrusion-dreams | presence, log BF = 65.5
```

and past the cap: `intrusion-dreams | presence, log BF > 10,000`
(snapshot-tested, `tests/testthat/_snaps/plot-methods.md`).

MM's original complaint case is worth spelling out: the BC fit's
`intrusion-dreams` had `log10_bf = 130.9`, so the old title printed
`sprintf("%.1f", 10^130.9)` — a **131-digit** number.

### 3.4 F-038 classification — what actually changed

The classification was, as the brief says, already correct. MM read −1.474
against the wrong boundary because the column was log10 while everything
around it was natural log. Now:

- stored column: `log_bf`, natural log (`log10_bf` is gone package-wide —
  grepped; the only remaining `log10_bf` symbols are
  `prior_sensitivity_check()`'s, deliberately untouched, §5)
- boundaries: ±`log(evidence_threshold)`, printed in the header
- `boundary_distance()` no longer divides the standard errors by `log(10)` —
  they arrive on the logit (natural log-odds) scale and are now directly
  commensurable. Distances, and therefore the fragility flag, are unchanged in
  value; the conversion cancelled on both sides.
- `bf = exp(log_bf)` retained
- the calibration-study figure in `?verdicts` and the vignette is restated in
  the new unit: 0.25 log10 → **0.58** natural log
- `?verdicts` now states that `extract_inclusion_bf(log = TRUE)` returns the
  same natural-log scale

MM's example row re-reads correctly: log BF −1.474 is **undecided** at
threshold 10 (boundary −2.30), which is what MM expected all along.

### 3.5 F-018 — measured

Brief 01's probe, against the branch build, unchanged script:

```
=== WITHOUT _R_CHECK_LIMIT_CORES_ ===
detectCores(): 15
user 143.067  elapsed 36.406     user/elapsed ratio: 3.93

=== WITH _R_CHECK_LIMIT_CORES_=TRUE ===
detectCores(): 15
user 135.092  elapsed 68.094     user/elapsed ratio: 1.98
```

3.93 → **1.98** under the check, and 3.93 unchanged without it — matching
report 01's uncapped baseline exactly.

Applied at every site that resolves a worker count: `validate_sampler()` (so
`spec$sampler$cores` is honest and every downstream reader inherits it),
`simulate()`/`predict()` (`simulate_predict.R:135`), `refit_engine.R:110`
(which bypasses `validate_sampler()`), and `zratio_surface_build_cores()`.
`normalize_builder_cores()` keeps its extra Windows no-forking rule layered on
the shared predicate.

**One deliberate deviation from the brief**, flagged for MM: the brief asked
for `normalize_builder_cores()`'s semantics "Windows nuance included". Its
Windows rule is `return(1L)` — correct for `mclapply`, which cannot fork there.
The sampler's workers are **threads**, not forked sessions, so applying that
rule would serialize all Windows sampling for no reason. `normalize_parallel_cores()`
therefore carries the env-var cap and the `detectCores()` clamp but not the
fork rule; `normalize_builder_cores()` is now that function plus the fork rule,
so the env-var semantics have one definition.

### 3.6 F-040, F-041, F-047

F-040(a): print method and `@return` now say "the proportion of the curve that
does not, on a 0 to 1 scale (0.3 is 30% of the curve, not 0.3%)".

F-040(b): **`variables =` already shipped** in 0.2.0.0 — documented and tested.
MM's complaint was discoverability, not absence. Added `max_panels` (default 9)
and `page`, so a 17-variable fit draws the nine worst and says:

```
Showing page 1 of 2 (17 variables, worst departure first). Draw the rest with
page = 2, or select panels with variables = .
```

F-041: reproduced. `refit_convergence_gate()` (`R/refit_engine.R:401-413`)
reduced with `min()/max()` over `wc$ebfmi_*`, `wc$var_ratio`, `rhats`, `esss`
and `ind[, "Rhat"]`, all of which a degenerate source fit can leave entirely
`NA` — hence the warnings, and an `±Inf` that reads as a *passing or failing*
criterion rather than an absent one. New `finite_reduce()` reduces the finite
entries and returns a stated value otherwise, so an unassessable criterion
abstains exactly as an absent warmup check already does. The R-hat criteria
still decide the gate.

F-041, second part (**report only, no change**): the refit-iteration floor is
**intended and documented**. `refit_run_length()` (`R/refit_engine.R`) uses
warmup 500 / iter 1000 for warm NUTS refits regardless of the source fit's
schedule — 1500 iterations, matching MM's observation from an `iter = 10` fit.
`?prior_sensitivity_check` `@param iter,warmup` says so: "`NULL` (default) to
use the validated short schedule for warm NUTS refits and inherit the original
fit's schedule otherwise."

F-047:

```
Before: Arguments 'variable1' and 'variable2' must name two different variables.
After:  Arguments 'variable1' and 'variable2' both resolve to 'intrusion', but an
        edge joins two different variables. Name the other end of the edge, for
        example plot_edge_posterior(fit, 'intrusion', 'dreams').
```

---

## 4. New findings

### 4.1 NEW (major) — `bgm(variable_type = "continuous", precision_graph_prior = "hierarchical")` crashes at 2 or 3 variables

Found while writing the F-036 tests. A documented, reachable configuration
errors outright:

```
q=2 -> ERROR: replacement has 1 row, data has 0
q=3 -> ERROR: replacement has 1 row, data has 0
q=4 -> OK
q=5 -> OK
q=6 -> OK
```

Pre-existing on `develop` and untouched by this branch (`R/zratio_surfaces.R`,
which I did not modify except for the core guard).

**Root cause, pinpointed.** `zratio_anchor_grids(cap)` filters the bipartite
anchor grid with `bip[bip$n <= cap, ]`; its smallest size is 4, so at `cap ≤ 3`
the result has **zero rows**. `zratio_cap_tier()` returns a 0-row frame
unchanged, and `zratio_build_surfaces()` then does
`bip_jobs$fam = "bip"` (`R/zratio_surfaces.R:411`) — assigning a length-1 value
into a 0-row data frame, which is the error. The bipartite family is genuinely
empty at q ≤ 3 (a bipartite bridge needs 2+2 nodes), so the honest behaviour is
a NULL bip surface, not a crash.

**Not fixed here** — it is in the PR #193/#194 hierarchical machinery, outside
this brief's authorized list, and the right shape of the guard (NULL family vs.
fall back to the additive path) is a call for whoever owns that surface. It
looks like a few lines. Flagging it as a release blocker candidate: `bgm()`
accepts the argument and then dies.

### 4.2 NEW (minor) — the correction-table disk cache key omits the package version

`ggm_ctable_v1_q17_delta1.4166067_eta1_normal_shape1_g120_ns2000_nw500_sd3_gibbs.rds`
carries a **schema** version (`v1`) and the model cell, but not the bgms
version. The Z-ratio surface cache next to it does carry it
(`zratio_surf_v2_0.2.0.0_...`), so the convention exists and this one departs
from it.

Observed concretely: MM's q = 17 table (built by the rc1 install) and mine
(built by `develop`) share a key and differ —
`max|edens_raw diff| = 0.314`, `max|edens diff| = 0.042`, `num_repaired` 1 vs 0.
Within one build the sweep is exactly deterministic (two clean runs agreed to
0), so this is a genuine code difference between the two versions being served
from one cache slot. A user upgrading bgms silently reuses the older version's
table.

### 4.3 NEW (minor) — a bgms computation cannot be interrupted

`Rcpp::checkUserInterrupt()` is called only in `src/mrf_simulation.cpp`. No
sampler translation unit calls it, so Ctrl+C during any chain — a fit, a
correction sweep, a Z-ratio surface build — is queued and not acted on until
the C++ call returns. Demonstrated in §2.3 on a 69-second sweep. Benign for
correctness (§2.3), but it is why MM's interrupt "did nothing".

### 4.4 NEW (note) — the fragility flag was wrong on mixed fits

Consequence of F-037, recorded separately because it is not obvious from that
finding's statement: `mcse` and `draws` reach `build_verdicts()` in the raw
block layout while `log_bf` arrived in row-major order, so on a mixed fit both
standard errors, both boundary distances and the `fragile` column belonged to a
different edge. Repaired by the same fix; no separate change needed.

---

## 5. Consistency flag (F-038 item 4) — `prior_sensitivity_check()` left in log10

As instructed, **not changed**. Its figures ship in the tutorial manuscript.
All log10 usages, for the terminology sync:

*Code, `R/prior_sensitivity.R`:* `371` (`lthr = log10(...)`), `597` (section
comment), `605`, `630` (curve construction), `701`
(`chosen_scale_log10_bf`), `767-768` (`log10_bf`, `log10_bf_mcse` fields),
`820`, `825-826` (noise-band comments), `1152` (the "…differed by up to %.2g
log10 BF" noise line MM saw), `1226` (`thr = log10(...)` in the plot),
`1240`, `1256`, `1306`, `1314`, `1321-1322`, `1343` (plot reads of
`x$log10_bf`).

*Roxygen/Rd:* `R/prior_sensitivity.R:66, 87, 94, 109, 155, 183, 1205`;
`man/prior_sensitivity_check.Rd:67, 102, 138, 159, 166, 181`;
`man/plot.bgms_prior_sensitivity.Rd:22`.

*Vignette:* `vignettes/prior-sensitivity.Rmd:101, 143, 148, 166`.

Note the resulting inconsistency is now visible to a user in one session:
`verdicts(fit)` prints natural log, `prior_sensitivity_check(fit)` prints
log10, and the two speak about the same edges. Whichever way MM settles it,
`tolerance = 0.5` and the `|log10 BF| <= 3` threshold-relevance window are
tuned constants that would need re-expressing (0.5 log10 = 1.15 nats,
3 log10 = 6.91 nats), not just relabelling.

---

## 6. Verification gate

All four gates pass, from a clean install of the branch.

**1. New/changed test files pass.** See gate 3.

**2. Full CRAN-mode suite** (`testthat::test_dir` without `NOT_CRAN`):

```
=== TOTALS ( cran ) ===
 failed   error warning skipped  passed
      0       0       0     234    7331
elapsed: 174.2 s
```

(174 s rather than report 01's 67 s: this run is `test_dir` against the source
tree on a machine also running the `R CMD check` build, not the in-check timing.)

**3. `NOT_CRAN=true` for the touched areas** (`test-verdicts.R`,
`test-calibration-check.R`, `test-plot-methods.R`, `test-prior-sensitivity.R`,
`test-prior-inclusion-probabilities.R`, `test-extractor-functions.R`,
`test-validate-sampler.R`):

```
=== TOTALS ===
 failed   error warning skipped  passed
      0       0       0       5    1045
```

The 5 skips are the `BGMS_RUN_SLOW_TESTS` prior-sensitivity certifications,
skipped on `develop` too.

**4. `R CMD check --as-cran`**, from a `git archive` export of the branch head,
report 01's recipe verbatim:

```
* checking examples ... OK
* checking examples with --run-donttest ... [468s/240s] OK
* checking tests ...
  Running 'testthat.R' [79s/63s] OK
* checking re-building of vignette outputs ... [77s/41s] OK
* checking CRAN incoming feasibility ... NOTE     (Date field over a month old)
* checking HTML version of manual ... NOTE        (HTML Tidy too old)
Status: 2 NOTEs
```

**No new NOTEs, WARNINGs or ERRORs vs report 01's baseline** — the same two
known NOTEs, both environment/metadata, neither a code defect.

One method note: checking a tarball built *directly from a git worktree* adds a
spurious third NOTE ("Found the following hidden files and directories: .git").
A worktree's `.git` is a 98-byte file rather than a directory, so `R CMD build`
does not drop it automatically, and `.Rbuildignore` has no `^\.git$` line. This
is an artifact of the build method, not of the branch — the archive-based check
above is the apples-to-apples comparison. A `^\.git$` line in `.Rbuildignore`
would harden against it; recorded, not changed.

**Pre-commit checks** (per `.github/copilot-instructions.md`):
`styler::style_pkg(style = bgms_style)` — run; changes to files this branch does
not touch were reverted, so the diff stays scoped.
`lintr::lint_package()` — **no lints found**.
`roxygen2::roxygenise()` — run; 4 Rd files regenerated and staged with their
commits.

---

## 7. Open questions for MM

1. **F-036 route (gating question).** Commit `ff29d177` changes the default
   route for hierarchical fits to the closed form. The claim is in §2.2. It is
   the specification `bgm()` documents and PR #193 tests for, so I judge it
   settled — but it is your call, and dropping that one commit reverts it
   cleanly. **The number is wrong today either way**; if the claim were
   rejected, the fallback is not the status quo but an explicit fence.
2. **F-036 has a release consequence.** Any hierarchical-spec result computed
   with 0.2.0.0 as it stands carries the +0.315-nat shift (at 17 variables,
   Bernoulli(0.5); the size depends on the cell). If any manuscript or tutorial
   number came from a hierarchical fit's Bayes factors, it needs re-running.
3. **F-037 has the same consequence, larger.** Every published mixed-fit
   `verdicts()` table from 0.2.0.0 is mislabelled. If the ADHD/mixed verdicts
   are being re-run anyway for the association-scale fix, this folds in.
4. **§4.1 — the q ≤ 3 hierarchical crash.** Do you want it fixed on this branch
   (it looks like a few lines in `zratio_build_surfaces()`), or does it belong
   to whoever owns the Z-ratio surface?
5. **§4.2 — the correction-table cache key.** Adding the package version to the
   key costs one rebuild per cell per release. Worth it before CRAN?
6. **F-038's cap value.** I used 10,000 nats, per your example. It only ever
   fires on `±Inf` in practice — a 17-variable Wenchuan fit's largest is ~300 —
   so nothing real is being truncated. Say if you want it tighter.
7. **Commit granularity.** One commit per finding, as asked. F-038 and F-039
   share `20b56996`: they are the same rename through the same functions and
   splitting them would have left an intermediate commit that does not build.

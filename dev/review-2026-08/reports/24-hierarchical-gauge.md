# Report 24 — hierarchical default, gauge harm wiring, the sweeps program, the mixed degenerate guard

Brief 24. Code changes authorized. Branch `fix/hierarchical-gauge`, based on
`origin/develop` = `5bbec810` (contains `1aefd2ff`), 4 commits — **not pushed**.
Worktree `~/bgms-review/wt-fix15`; the Dropbox tree's branch was never switched
and nothing was built in it.

`origin/develop` advanced to `1cc4857d` while this batch ran (brief 21's
summary-zero-marking merge plus the lead's `PLAN.md` commits). None of those
files are mine — the compare surface was never opened — so the branch rebases
cleanly onto them as far as file overlap goes.

No off-limits file was touched. The branch changes `R/bgm.R`, `R/bgm_spec.R`,
`R/build_spec.R`, `R/build_output_mixed_mrf.R`, `R/zratio_gauge.R`,
`R/zratio_surfaces.R`, `R/sample_ggm_prior.R`,
`src/models/mixed/mixed_mrf_{model,gradient}.cpp`, three regenerated `man/*.Rd`,
seven `tests/testthat/test-*.R` files, and this report. `NEWS.md`,
`vignettes/`, the plot files and every compare-surface file are untouched, and
`tests/testthat/_snaps/` is **byte-unchanged** — no snapshot was re-recorded.

**Headline.** Tasks 1, 2 and 4 landed. Task 3 ran the full measurement and then
took the brief's **STOP branch**: raising the gauge's audit precision does not
move the harm alarm off the threshold, because the dispersion report 13 found
is not the gauge's Monte-Carlo error. At 16 sweeps — eight times the audit
sample — the 12-seed max is 0.0108, still *above* the 0.01 threshold rather
than at or below the 0.005 target, and the between-seed sd is flat (0.0041 →
0.0039). No default was changed and the F-103 block stays in T2. The full table
is §4; the threshold-versus-margin question returns to the maintainer.

Two things beyond the brief's letter are worth the lead's eye. First, flipping
the default made the vacuous-spec advisory fire on **every default ordinal
fit** — the flagship path — because the argument now arrives without being
asked for; §2 explains the fix. Second, the degenerate-block guard closes the
`p_ == 0` route, but **eight more unguarded `p_ - 1` upper-triangle bounds**
exist beyond the seven the brief named; I left them alone rather than widen
scope, and they are listed in §7.

---

## 1. What was done

Numbered as the brief numbers them; the commits are in the order 1, 2, 4, 3,
because the task-3 measurement is the heavy tail and ran last.

| brief task | What | Finding | Commit | Type | § |
|---|------|---------|--------|------|---|
| 1 | Hierarchical default flip | F-010 | `2576f457` | `fix(ggm):` | §2 |
| 2 | Mixed gauge harm wiring | F-022 | `8949ba6b` | `fix(mixed):` | §3 |
| 3 | Sweeps program — un-hardwire, measure, **STOP** | F-103 | `d3cc69f5` | `fix(ggm):` | §4 |
| 4 | Degenerate guard + loop bounds | F-123 (b) | `4aa99821` | `fix(mixed):` | §5 |

---

## 2. Task 1 — F-010, the default flip (`2576f457`)

### The change

`bgm()`'s `precision_graph_prior = c("joint", "hierarchical")` becomes
`c("hierarchical", "joint")`, so `match.arg` resolves an unnamed argument to
`"hierarchical"`. `"joint"` is unchanged and one argument away.

Two structural points came with it.

**(a) `bgm()` now resolves its own default and passes a scalar.** `bgm_spec()`
keeps `c("joint", "hierarchical")`. This is deliberate, not an oversight:
`bgm_spec()`'s other caller is `bgmCompare()` (an off-limits file), which never
names the argument and has no continuous precision block for it to refer to.
Flipping `bgm_spec()`'s own default would have changed what a compare fit
records and would have fired the vacuous-spec advisory on every `bgmCompare()`
call. Passing a two-element vector down instead would have hit
`match.arg`'s `"'arg' must be of length 1"` error, since the vector no longer
matches the callee's own choices. The asymmetry is commented at both ends and
is listed as a follow-up in §7.

**(b) The vacuous-spec advisory now reports a request, not a default.** This is
the one substantive behaviour change beyond the flip itself, and it is not
optional. `zratio_vacuous_spec_notice()` fires when
`precision_graph_prior == "hierarchical"` and the model has no continuous
precision block, telling the user the argument has no effect. With
`"hierarchical"` as the *default*, that condition is met on **every ordinal
fit** — `bgm(x)` at its defaults, the package's flagship call — so the flip as
briefed would have printed

> `precision_graph_prior has no effect for this model: it normalizes the`
> `continuous precision prior across graphs, and this model has no continuous`
> `precision block. The fit is the same under either value.`

on every default ordinal fit with `verbose = TRUE`. The advisory is about a
user's choice; a default is not a choice. `bgm()` therefore records
`hasArg(precision_graph_prior)` before resolving, and `bgm_spec()` passes it to
the notice as `explicit`. An explicit
`bgm(x, precision_graph_prior = "hierarchical")` on ordinal data still fires
it, exactly as before.

### Documentation swept

- `R/bgm.R` `@param precision_graph_prior`: the `(default)` tag moved from the
  `"joint"` item to the `"hierarchical"` item; the trailing `Default: "joint".`
  sentence rewritten to name hierarchical and say what it buys and costs; the
  vacuous-case bullet now states the named-only rule.
- `R/bgm.R` `@param edge_prior`: the normalizing-constant-correction paragraph
  is now scoped to `precision_graph_prior = "joint"` and says that the default
  hierarchical path needs no correction table and builds none. Left unscoped it
  would have told a default-fit reader to expect a minutes-long one-time table
  build that the default never performs.
- `man/bgm.Rd` regenerated. `devtools::document()` clean, **NAMESPACE
  unchanged**.
- Nothing else in `R/` or `man/` names either specification as the default.
  `sample_graph_prior()`/`sample_ggm_prior()` have their own `spec` argument,
  which already defaulted to `"hierarchical"`; `extract_prior_inclusion_-`
  `probabilities()` describes both routes without naming a default; `README.md`
  mentions neither.

### `extract_arguments()` check (gate item 6)

A default GGM fit — `bgm(Y, variable_type = "continuous", …)` with no
`precision_graph_prior` — reports:

```
extract_arguments(fit)$precision_graph_prior  ->  "hierarchical"
is.null(fit$zratio_diag)                      ->  FALSE   (the gauge ran)
eval(formals(bgm)$precision_graph_prior)       ->  c("hierarchical", "joint")
```

### Touched expectations, with derivation and (a)/(b) classification

Classification per the brief: **(a)** = the test's intent is the joint path, so
it now names `precision_graph_prior = "joint"` and its expectations are
untouched; **(b)** = the test's intent is default behaviour, so the expectation
is re-derived under hierarchical.

| # | File / test | Class | Derivation |
|---|-------------|-------|-----------|
| 1 | `test-parameter-prior-cpp.R` — "GGM edge selection works with beta_prime_prior interaction" | **(a)** | *Was failing.* The Z-ratio constants are derived for a Normal or Cauchy slab, so `bgm_spec()` rejects a `beta_prime_prior()` slab under the hierarchical spec. The default flip turned this fit into that error. The test is about the beta-prime slab, which lives on the joint path only. Expectations unchanged; the route is now named. |
| 2 | `test-parameter-prior-cpp.R` — "GGM beta_prime_prior interaction rejects an eta-frame diagonal prior" | **(a)** | *Was failing.* Same cause, and worse: the error under test is `"scale parameter"` (the eta-frame diagonal prior), and the hierarchical slab rejection was pre-empting it with a different message, so `expect_error()` matched nothing. Route named; expectation unchanged. |
| 3 | `test-joint-realized-prior-notice.R` — "bgm wires the realized-prior notice into the joint spec" | **(a)** | *Was failing* under `NOT_CRAN=true` (it is `skip_on_cran()`). `expect_true(fires(fit_messages()))` relied on the bare default being joint. The firing case now names `"joint"`. A third assertion was **added**: the default fit does not fire it — which is the flip's user-visible payoff and was previously untestable. |
| 4 | `test-mixed-correction.R` — all five `bgm()` end-to-end fits | **(a)** | *Two were failing:* the `"correction table"` build announcement no longer appeared (the default path builds no table), and the sbm `"slope curve is not resolvable"` warning no longer fired. The other three (`sbm` fit, one-continuous skip, two-continuous beta-bernoulli) still passed but had silently stopped exercising the route they name — the correction is a joint-specification artifact and the hierarchical path tracks `Z(Gamma)` in the edge moves instead. All five now name `"joint"`; a section comment records why. |
| 5 | `test-prior-inclusion-probabilities.R` — "GGM beta-bernoulli prior PIPs match the prior-only chain" | **(a)** | *Was failing.* Under hierarchical the extractor returns the edge prior's own marginal exactly, so `offdiag[1] < 0.5` is false (it is 0.5) and the prior-only-chain agreement is off by 0.15. Both assertions describe the tilted joint route. Route named. |
| 6 | `test-prior-inclusion-probabilities.R` — "the joint block reweights a fixed Bernoulli prior at delta = 0" | **(a)** | *Was failing.* `expect_lt(pip, 0.35)` got 0.50 — under hierarchical there is no reweighting to observe, which is the point of the specification. The test's name states its route. |
| 7 | `test-prior-inclusion-probabilities.R` — "mixed beta-bernoulli prior PIPs split by edge class" | **(a)** | *Was failing.* `cc[1] < 0.5` got 0.5: the three edge classes stop splitting under hierarchical, because the continuous-continuous class no longer carries the joint-block density. Route named. |
| 8 | `test-prior-inclusion-probabilities.R` — "GGM SBM prior PIPs come from a cached deterministic chain" | **(a)** | *Was passing* — its assertions (`0 < pip < 1`, identity across `recompute`) hold on either route. But under hierarchical there is no prior-only chain to cache, so it had stopped testing the cached deterministic chain it names. Route named to restore the test's intent. |
| 9 | `test-bgm-hier-spec.R` — "the joint default is unchanged" | **(b)** | *Was failing,* and had to: it is the assertion the maintainer decision inverts. Re-derived — an unnamed `precision_graph_prior` on continuous data now resolves to `"hierarchical"`, so `fit@arguments$precision_graph_prior == "hierarchical"` (was `"joint"`) and `fit@zratio_diag` is non-`NULL` (was `NULL`), because `zratio_active` is now true and the gauge runs. Renamed to "a continuous fit defaults to the hierarchical spec (F-010)". The joint behaviour it used to assert is **kept** in the same block as an explicit-argument control, so the flip is a change of default rather than a loss of coverage. |

Two tests were **added** for the advisory change in (b) above:

- `test-bgm-hier-spec.R` — "the vacuous-spec notice reports a request, not the
  default": the three cells of `zratio_vacuous_spec_notice()` (no block +
  explicit → message; no block + default → silent; block present → silent).
- `test-bgm-hier-spec.R` — "bgm keeps the ordinal default fit silent about the
  spec": end-to-end, a default ordinal `bgm()` fit at `verbose = TRUE` prints
  no vacuity message, and an explicit `precision_graph_prior = "hierarchical"`
  on the same data does.

**No snapshot was re-recorded**; `tests/testthat/_snaps/` is byte-unchanged on
the whole branch. The two snapshot files are `plot-methods.md` and
`verdicts.md`, and neither snapshots a default GGM fit's printed output, so no
snapshot could have gained a gauge line. Nothing here needed the brief's
"re-record with a stated reason" clause.

### Vignette / docs-site re-bake list (report-only — nothing edited)

`vignettes/` is off-limits, so this is the list the Phase-3 docs re-bake needs.
There is **no `_pkgdown.yml` in the tree**, so the docs site is built from
`man/` + `vignettes/` + `README.md`; `README.md` names neither specification,
so it needs nothing.

| # | Location | What the flip does to it |
|---|----------|--------------------------|
| 1 | `vignettes/intro.Rmd:107–121`, "Continuous data (GGM)" | The example `fit_ggm = bgm(continuous_data, variable_type = "continuous", seed = 1234)` (line 119, `eval=FALSE`) is now a **hierarchical** fit: it builds the Z-ratio surface, attaches `fit$zratio_diag`, and `summary(fit_ggm)` can gain gauge lines. The chunk is `eval=FALSE` so nothing re-renders, but the section never names a specification and now silently demonstrates one. Needs a sentence saying which. |
| 2 | `vignettes/diagnostics.Rmd:235` | Opens the trust-gauge section with "When you use the hierarchical graph prior (`precision_graph_prior = "hierarchical"`)…", framing the gauge as opt-in. It is now what a continuous fit does by default. |
| 3 | `vignettes/diagnostics.Rmd:239` | "The gauge runs by default." Still true, and now true of the *default* continuous fit rather than only of the opt-in path — the sentence's scope widened without the words changing. Worth restating so a reader does not read "by default" as "by default once you've opted in". |
| 4 | `vignettes/diagnostics.Rmd:273` | "Raise `options(bgms.zratio_gauge_sweeps)` above its default of 2 and re-fit." The number **2 is still correct** after task 3 (see §4: no default changed). Flagged here only so the re-bake re-checks it against whatever the maintainer decides on F-103. |
| 5 | `vignettes/diagnostics.Rmd:275` | "Third, and only if a flag persists, consider the joint prior specification…". The advice is unchanged and still right, but it is now a step *away from* the default rather than back to it. "`bgms` will never make that switch for you" still holds. |
| 6 | `vignettes/checking-your-model.Rmd:237` | Next-steps bullet naming "the hierarchical prior trust gauge" as one more thing the diagnostics vignette covers. It reads as an optional feature and is now part of the default continuous workflow. |

`vignettes/comparison.Rmd:56–62` ("One default differs between the two entry
points…") was checked and needs nothing: `bgmCompare()` has no
`precision_graph_prior` argument at all, so this is an absent argument, not a
differing default.

### Proposed NEWS wording — VERBATIM, for the lead to land

The GGM never shipped, so all of this is **feature** wording for the 0.1.6.3
reader; nothing is phrased as a change from joint. Four edits, all inside the
existing 0.2.0.0 "New features" section.

**(N1)** `NEWS.md:615`, first bullet of the realized-prior group. Replace
`Under the joint precision-graph specification (the default), a fit with edge`
with:

```
* Under `precision_graph_prior = "joint"`, a fit with edge
```

(the rest of that bullet is unchanged).

**(N2)** `NEWS.md:631–635`, the opening of the hierarchical bullet. Replace

```
* `bgm(precision_graph_prior = "hierarchical")` composes the edge prior and the
  precision prior as `p(Gamma) p(K | Gamma)` with `p(K | Gamma)` normalized per
  graph, so the graph marginal is exactly the edge prior; under the default
  `"joint"` specification it is that prior reweighted by the per-graph
  normalizer.
```

with

```
* `precision_graph_prior = "hierarchical"` is the default on continuous and
  mixed models. It composes the edge prior and the precision prior as
  `p(Gamma) p(K | Gamma)` with `p(K | Gamma)` normalized per graph, so the
  graph marginal is exactly the edge prior; under
  `precision_graph_prior = "joint"` it is that prior reweighted by the
  per-graph normalizer. Both specifications are fully supported; the default
  is the one whose graph marginal is the prior you wrote down.
```

**(N3)** `NEWS.md:657–659`, inside the "argument is accepted where the choice is
vacuous" bullet. Replace

```
  refer to; the fit is accepted and a message reports this when
  `verbose = TRUE`.
```

with

```
  refer to; the fit is accepted, and a message reports this when the argument
  was named and `verbose = TRUE`. Inheriting the default is silent: since
  `"hierarchical"` is the default, the value reaches every fit, and the
  message reports a user's choice rather than a default nobody made.
```

**(N4)** `NEWS.md:752`, the gauge bullet's opening sentence. Replace

```
* A trust gauge runs by default on the deployed hierarchical path, with
```

with

```
* A trust gauge runs by default on the deployed hierarchical path — which is
  to say on a default continuous or mixed fit with edge selection — with
```

No other NEWS line names or implies the joint default. `NEWS.md:64–66` ("Under
the joint precision-graph specification the normalizer correction table is
keyed on the interaction prior…") was checked and needs nothing: it names the
specification rather than calling it the default, and it sits under *Changed
defaults* for the **interaction prior**, which is a separate change.

---

## 3. Task 2 — F-022, the mixed harm channel (`8949ba6b`)

### Before / after

`R/build_output_mixed_mrf.R:305` called
`summarize_zratio_gauge(zratio_chains, verbose = TRUE)` with no `harm_inputs`,
so the whole harm block of the summary was dead on the mixed path. Same
fixture, same seed (2 ordinal + 8 continuous, `n = 50`, `beta_bernoulli_prior(9, 1)`,
`precision_graph_prior = "hierarchical"`, 2 gauge sweeps):

| field | before (`origin/develop` build) | after (branch build) |
|-------|--------|-------|
| `harm_pred` | `NA` | **0.001343** |
| `amplification` | `NA` | 2.912 |
| `kappa` | `NA` | 0.6391 |
| `harm_flag` | `FALSE` (vacuously — `NA > threshold` is never `TRUE`) | `FALSE` (measured, 0.00134 against 0.01) |
| `flip_rate` | 0.0003178 | 0.0003178 (unaffected) |
| `n_ent` / `n_ref` | 42 / 42 | 42 / 42 (unaffected) |

Both rows are the same fixture and seed run against a baseline package built
from `git archive origin/develop` and against the branch build. The
`harm_flag = FALSE` in the "before" column is the point: the channel could
never fire on a mixed fit, in either direction, so a mixed fit that genuinely
needed the alarm would have been silent.

After task 1 the gauge is part of the default mixed experience, so this was a
channel that would have shipped permanently blank on the default path.

### Why it is not a one-line mirror of the GGM builder

The GGM builder can pass one `pip` vector because on a GGM the audited block
*is* the whole graph. On a mixed model the two come apart:

- **The audited block** is the continuous-continuous edges only. The Z-ratio
  enters those moves alone, so the gauge's `pair_i`/`pair_j` stream indexes the
  *continuous subgraph*, and `summarize_zratio_gauge()` recovers the subgraph's
  size from `length(pip)` via `q = round((1 + sqrt(1 + 8 n_edges)) / 2)`.
  Passing the full mixed `pip` would have made `q` wrong and mapped every audit
  record onto the wrong edge — the same class of bug the existing
  "maps audit records onto the correct edge" regression pins. The
  `[Gxx_ut | Gyy_ut | Gxy]` indicator layout puts that block at the
  `choose(q, 2)` columns after the discrete ones, in the same upper-triangle
  order, which is what the builder now slices.
- **The feedback pool** is every edge. `BetaBernoulliEdgePrior::update()`
  (`src/priors/edge_prior.h:108`) counts included edges over the whole
  `(p+q) x (p+q)` indicator matrix and draws one shared theta, so a
  perturbation of the continuous-continuous decisions feeds back through *all*
  the edges.

Feeding the audited block as both would have computed the gain
`n m_bar / (theta(1-theta)(a + b + n))` from a pool smaller than the real one.
That is not a rounding difference. On a homogeneous graph the gain collapses to
`n / (a + b + n)` and the amplification to `1 / (1 - gain)`, so on a 5-discrete
+ 5-continuous fixture with `beta_bernoulli_prior(1, 1)` the audited block gives
`n = 10`, gain 0.833, amplification **6.0**, while the real pool gives `n = 45`,
gain 0.957, amplification **23.5** — the wiring would have under-reported the
harm it exists to report by a factor of four.

So `summarize_zratio_gauge()` and `zratio_harm_inputs()` take an optional
`pool_pip` alongside `pip`. The numerator still weights each audit record by
its own audited edge's sensitivity; only the gain reads the pool. **Omitting
`pool_pip` means pool == audited block, so the GGM path is untouched** — the
`test-zratio-gauge.R` harm tests, which pass `pip` alone, are unchanged and
still pass. I did not take the brief's STOP branch, because no ingredient was
missing: both blocks are available in the builder, and the exported contract
grew one optional field rather than being forced.

Tests added: a seconds-scale mixed hierarchical fit whose gauge summary carries
finite `harm_pred`/`amplification`/`kappa` with `n_ref > 0`; and a unit test
that the gain follows the pool while the numerator follows the audited block,
including an explicit assertion that `pool_pip = NULL` reproduces the GGM
numbers exactly.

---

## 4. Task 3 — F-103, the sweeps program: measured, then STOPPED (`d3cc69f5`)

### (a) Un-hardwired — landed

`sample_ggm_prior()`'s `gauge_sweeps = if (isTRUE(zratio_diagnostics)) 2L else 0L`
becomes
`if (isTRUE(zratio_diagnostics)) zratio_gauge_sweeps() else 0L`, the same
resolution `run_sampler_ggm()` and `run_sampler_mixed_mrf()` use.
`zratio_diagnostics` remains the on/off switch and is *not* the precision:
`FALSE` is still 0 sweeps whatever the option says. At the option's default
value of `2L` this is a **no-op**, so no existing fit or test changes.

Pinned by a new test: the audit references a fixed cap of pairs per sweep, so a
saturated audit's `n_ref` grows exactly one cap per sweep — `n_ref(2) = 2 n_ref(1)`
and `n_ref(6) = 6 n_ref(1)` — and `zratio_diagnostics = FALSE` still yields no
gauge at `bgms.zratio_gauge_sweeps = 6L`.

### (b) Measured — the table

`harm_pred` for the negative control `biased_evidence_free_fit(2)`
(`p = 16`, shipped `n_samples = 1200`, `n_warmup = 500`, Gibbs,
`normal_prior(0.5)` slab, `gamma_prior(shape = 2, rate = 6)` diagonal,
`beta_bernoulli_prior(9, 1)`), seeds 1–12, one fit at a time, 48 fits,
326 s of wall time in total. `harm_threshold = 0.01` untouched.

| `gauge_sweeps` | `n_ref` | mean | sd | max | max ÷ threshold | > threshold | s / fit |
|---|---|---|---|---|---|---|---|
| **2** (shipped) | 50 | 0.004832 | 0.004073 | **0.012533** | 1.253 | **1 / 12** | 5.64 |
| 4 | 100 | 0.004762 | 0.004286 | **0.012420** | 1.242 | **1 / 12** | 6.08 |
| 8 | 200 | 0.004486 | 0.003731 | **0.010554** | 1.055 | **1 / 12** | 6.90 |
| 16 | 400 | 0.004534 | 0.003860 | **0.010757** | 1.076 | **1 / 12** | 8.54 |

**Sanity anchor: the sweeps-2 row reproduces report 13 exactly.** Sorted, the
twelve values are 0.00036, 0.00076, 0.00095, 0.00111, 0.00203, 0.00336,
0.00435, 0.00662, **0.00794 (seed 7, the fixture's own seed)**, 0.00876,
0.00921, 0.01253 — report 13 §9b's list to five decimals, in the same order,
with mean 0.00483 and sd 0.00407 matching its table. The seed set is 1–12.

Per seed:

| seed | 2 | 4 | 8 | 16 |
|---|---|---|---|---|
| 1 | 0.00435 | 0.00391 | 0.00318 | 0.00269 |
| 2 | 0.00876 | 0.00968 | 0.00874 | 0.00927 |
| 3 | 0.00036 | 0.00017 | 0.00113 | 0.00112 |
| 4 | 0.00111 | 0.00033 | 0.00064 | 0.00108 |
| 5 | 0.00203 | 0.00123 | 0.00075 | 0.00067 |
| 6 | 0.00095 | 0.00067 | 0.00069 | 0.00073 |
| 7 | 0.00794 | 0.00808 | 0.00932 | 0.00945 |
| 8 | 0.00921 | 0.00826 | 0.00657 | 0.00601 |
| 9 | 0.00336 | 0.00337 | 0.00410 | 0.00380 |
| 10 | 0.00662 | 0.00793 | 0.00672 | 0.00753 |
| 11 | 0.00076 | 0.00110 | 0.00143 | 0.00130 |
| **12** | **0.01253** | **0.01242** | **0.01055** | **0.01076** |

### (c) The choice — **STOP**

The target was a 12-seed max at or below **0.005**, half the threshold. The
achieved max at the largest count tested is **0.01076**, which is not below the
threshold at all, let alone with margin: the margin is **−0.00076**, i.e. the
max still *exceeds* `harm_threshold` by 7.6 %. Five of the twelve seeds sit
above 0.005 at **every** sweep count. So the brief's STOP condition is met —
16 sweeps leaves the max straddling the threshold — and **no default was
changed**. `harm_threshold = 0.01` was not touched. Steps (d) and (e) are
consequently not applicable; §4(e) below records what that means for the tier.

**Why the lever does not work here, stated as evidence rather than as a
verdict.** The lever *does* do what it was aimed at — it sharpens each seed's
estimate:

- Mean absolute change between adjacent columns falls from 0.00052 (2→4) to
  0.00031 (8→16): the per-seed estimates have converged.

But the dispersion report 13 found is not the gauge's Monte-Carlo error:

- The **between-seed sd is flat**: 0.00407, 0.00429, 0.00373, 0.00386 while the
  audit sample grows 8×. Gauge Monte-Carlo error would have shrunk by about
  √8 ≈ 2.8×.
- The per-seed values are a **stable property of the seed**: Pearson r = 0.94
  (Spearman 0.82) between the 2-sweep and 16-sweep columns. Seed 12 stays
  highest at every count and seed 6 stays near the bottom at every count.
- Seed 12's estimate **converges to a value above the threshold**
  (0.01253 → 0.01242 → 0.01055 → 0.01076), so its alarm is not an artifact of a
  small audit. At this cell, on that seed, the harm statistic genuinely exceeds
  the product constant on a healthy negative-control fit.

That is the maintainer's question in its sharpest form: either
`harm_threshold = 0.01` sits below what the shipped gauge genuinely produces on
healthy fits at a cell inside the surface's validated range, or the harm
estimator overstates the distortion. Neither is a sweep-count question, and
neither is a test-constant question. Runtime for the counts, should the answer
turn out to want one: +0.45 s, +1.26 s and +2.91 s per fit at 4, 8 and 16
sweeps against 2, on this `p = 16`, 1200-sample fixture.

### (e) The F-103 block — **not restored**

`tests/testthat/test-zratio-gauge.R:235` keeps its
`skip_unless_certification()`. Restoring it to T1 at the unchanged default
would put back exactly the flake report 13 declined to add: 1 of 12 seeds
exceeds the threshold at every sweep count, and report 13 established that the
Linux CI value (0.01071) is an ordinary draw from that spread rather than
platform drift — so the fixture's fixed seed 7 does not protect the nightly.
The block's header comment now records the second measurement, its conclusion
and the pointer to this report, so the next reader does not re-run the same
program.

The sibling block at `:215` ("a known-biased evidence-free fit fires the harm
channel", the *positive* control at shape 0.4) was parked in T2 in the same
action but is **not** F-103's block and was not in scope. It is listed in §7 as
an open item: it asserts that a flag *fires*, so its restoration needs a
12-seed pass of its own, which this program did not run.

---

## 5. Task 4 — F-123 (b), the degenerate guard and the loop bounds (`4aa99821`)

### The guard

`build_spec_mixed_mrf()` stops when either block is empty, before it touches
anything. Option (a) — making degenerate blocks work — was rejected by the
maintainer, and the guard says why in the message: pure-ordinal data is fitted
by the OMRF and pure-continuous data by the GGM, and `bgm_spec()` routes both
before the mixed builder is reached. Reaching the guard means the model was
chosen internally rather than from the data, so the message says that rather
than blaming the user's `variable_type`.

Unit test on the internal, both degenerate directions
(`tests/testthat/test-bgm-spec.R`): 3 continuous / 0 discrete and 3 discrete /
0 continuous both error, and both messages name the counts. The test drives
`bgms:::build_spec_mixed_mrf()` directly — R's lazy evaluation means only the
arguments the guard reads have to be supplied, since it fires before the block
machinery is touched.

### The loop bounds

The seven bounds the brief named now take the `i + 1 < n` form the RB mirrors
in the same files already use:

| file | line | counter | context |
|---|---|---|---|
| `mixed_mrf_model.cpp` | 851 | `p_` | `get_vectorized_indicator_parameters()`, Gxx upper triangle |
| `mixed_mrf_model.cpp` | 858 | `q_` | same, Gyy upper triangle |
| `mixed_mrf_gradient.cpp` | 31 | `p_` | discrete index cache |
| `mixed_mrf_gradient.cpp` | 87 | `p_` | observed statistics for discrete edges |
| `mixed_mrf_gradient.cpp` | 146 | `p_` | discrete pairwise parameter unpacking |
| `mixed_mrf_gradient.cpp` | 498 | `p_` | discrete pairwise priors |
| `mixed_mrf_gradient.cpp` | 570 | `q_` | continuous pairwise priors |

### Same-seed identity proof

Baseline package built from `git archive origin/develop` into its own library
and run against the branch build. One mixed fit, 3 ordinal + 4 continuous,
`n = 60`, 2 chains, 200 iterations after 300 warmup, seed 314, run on **both**
the joint path (gauge off) and the hierarchical path (2 gauge sweeps):

```
identical(before, after)            TRUE
  joint draws                       TRUE
  hierarchical draws                TRUE
numeric cells compared              62,960
max |diff| main-effect draws        0
max |diff| indicator draws          0
sd of the main-effect draws         0.381   (non-degenerate)
```

The comparison covers every draw matrix returned in `fit$raw_samples` — main,
pairwise, indicator, Rao-Blackwellized inclusion, RB counts — plus the sampled
inclusion parameter and all four posterior-mean summaries. Bit-identical, as
required for all `p, q >= 1`.

No NEWS: the mixed model never shipped.

---

## 6. Verification gate

| # | Gate | Result |
|---|------|--------|
| 1 | Both suite tiers, 0 failures / 0 warnings | **PASS** — T0: `[ FAIL 0 \| WARN 0 \| SKIP 97 \| PASS 8475 ]`; T1 (`BGMS_RUN_SLOW_TESTS=true`): `[ FAIL 0 \| WARN 0 \| SKIP 61 \| PASS 8679 ]`. Both run with `NOT_CRAN=true`, so the `skip_on_cran()` blocks — which is where three of the nine failing expectations lived — actually execute. |
| 2 | `R CMD check --as-cran` on a `git archive` tarball, 2 baseline NOTEs | **PASS** — `Status: 2 NOTEs` on `bgms_0.2.0.0.tar.gz` (1.67 MB, built from `git archive HEAD` with pandoc on `PATH`, vignettes included). Both are the known baselines: "The Date field is over a month old" (F-003) and "Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy" (report 01 finding 10). No ERROR, no WARNING. Examples, tests (`[129s/102s] OK`) and vignette rebuilds (`[93s/50s] OK`) all pass. |
| 3 | `devtools::document()` clean; NAMESPACE unchanged | **PASS** — `roxygen2::roxygenise()` regenerates `man/bgm.Rd`, `man/summarize_zratio_gauge.Rd` and `man/sample_ggm_prior.Rd` and nothing else; `git diff NAMESPACE` is empty |
| 4 | Task-3 table complete with runtimes; F-103 block restored and proven 12/12, **or** the STOP branch taken and documented | **STOP branch taken** — table in §4(b) with per-fit runtimes; no default changed; block left in T2 with the reason recorded in the file and in §4(e) |
| 5 | Task-4 same-seed identity check | **PASS** — §5, 62,960 cells bit-identical on both specification paths |
| 6 | `extract_arguments()` default-fit check | **PASS** — §2, reports `"hierarchical"` |

Every re-tuned expectation is listed with its derivation and its (a)/(b)
classification in §2.

---

## 7. Findings beyond the brief

**F-24-1 (major) — the default flip silences nothing by itself; the vacuous
advisory had to be taught the difference between a request and a default.**
Covered in §2(b). Flagged as major because the flip as briefed would have added
a spurious message to `bgm(x)` — the package's most common call — and the
suite would not have caught it: no existing test asserts that a default ordinal
fit is silent about the specification. One now does.

**F-24-2 (minor) — eight more unguarded `p_ - 1` upper-triangle bounds in
`mixed_mrf_model.cpp`.** The brief named seven sites; the same
`for (size_t i = 0; i < p_ - 1; ++i)` pattern also appears at lines **559, 606,
656, 703, 748, 811, 1121 and 1195**, all over the discrete block and all
unguarded by an enclosing `p_ >= 2`. (The two `q_` analogues at 1130 and 1202
*are* guarded, by `if (q_ >= 2)`.) They are unreachable now that the task-4
guard closes the `p_ == 0` route, exactly as the seven were, so this is
hygiene, not a defect — but leaving the file half-converted invites the next
reader to assume the remaining ones are deliberate. I did not touch them:
widening a maintainer-scoped mechanical change is the lead's call, not mine.
The same-seed identity proof in §5 would cover them unchanged.

**F-24-3 (minor) — `bgm_spec()`'s `precision_graph_prior` default is now
deliberately out of step with `bgm()`'s.** `bgm()` defaults to
`"hierarchical"`, `bgm_spec()` to `"joint"`, because `bgm_spec()`'s other
caller is the off-limits `bgmCompare()`, for which the argument is inert and
for which flipping it would have fired the vacuous advisory on every compare
fit. Both ends carry a comment. Once the two compare batches land, the lead may
want to align them and drop the asymmetry — which needs a decision about what
a compare fit should record in `spec$prior$precision_graph_prior`, if anything.

**F-24-4 (minor) — the F-103 positive control at
`test-zratio-gauge.R:215` is still parked in T2 with no owner.** It was parked
in the same action as the negative control but is a different assertion (a flag
must *fire*), so the sweeps program says nothing about it. Restoring it needs
its own 12-seed pass at shape 0.4. Listed so it does not stay parked by
inheritance.

**F-24-5 (informational) — `harm_pred`'s dispersion is a property of the fit,
not of the audit.** The evidence is in §4(c): r = 0.94 across an 8× change in
audit size, and a flat between-seed sd. Whoever owns the threshold question
should know that the statistic is reproducible per seed, which makes it
tractable: a 12-seed sweep at a candidate threshold is a cheap, stable
measurement (326 s for 48 fits here), not a noise-limited one.

---

## 8. Open questions for the maintainer

1. **F-103, the actual question.** At `harm_threshold = 0.01`, the shipped
   gauge flags 1 of 12 healthy negative-control fits at a cell inside the
   surface's validated range, and the flagged seed's statistic converges to
   0.0108 rather than decaying with the audit size. Two dispositions are
   consistent with the data: raise the threshold to clear the observed
   distribution, or accept that the harm estimator overstates at this cell and
   revisit the projection. If the first, note what the numbers do and do not
   support — a threshold of 0.015 clears all 48 fits measured here (the overall
   max is 0.01253), but it does **not** buy the half-the-threshold margin the
   brief asked for: three to five of the twelve seeds sit above 0.0075 at every
   sweep count (4, 5, 3, 4 at 2, 4, 8, 16 sweeps).
   The margin the fixture can support and the margin the product wants are the
   two numbers to reconcile. Either way the sweeps lever is spent.
2. **Should the F-103 blocks return to T1 at all,** or should the gauge
   detector's tier contract be rewritten to say that a Monte-Carlo negative
   control belongs in the weekly tier? Report 13 called the T2 parking a
   contract violation; a third option is that the contract is wrong.
3. **The `bgm_spec()` / `bgm()` default asymmetry** (F-24-3) — align after the
   compare batches land, or leave it as the deliberate seam it is now?
4. **The eight remaining loop bounds** (F-24-2) — finish the file, or leave the
   named seven as the maintainer scoped them?

# Report 17 — docs + small-defect fix batch

Brief 17 (+ the lead's 2026-08-02 addendum: F-105, F-106). Code changes
authorized. Branch `fix/docs-batch`, 11 commits; archive branch
`archive/zratio-analytic-law`, 1 commit — **neither pushed**.

**Headline.** Every briefed item landed. Both majors are fixed and both are now
drift-proof rather than merely corrected: the inverted Bayes-factor reading
(F-087) reads its numbers through inline R expressions, and the mock
sensitivity report (F-086) is a live seeded chunk. All four gates are green:
both test tiers at 0 failures / 0 warnings / 0 errors, five vignettes rendering
with 0 warnings and 0 messages, and `R CMD check --as-cran` at the 2 baseline
NOTEs. The one judgement call worth the lead's eye is in §2, finding 1: two
Bug-fixes entries the brief listed as "verify at the tag" turned out **not** to
be reachable at `cran-0.1.6.3` in the form they describe, and I dropped them.

---

## 0. Base — brief §Setup, resolved

The brief's stop condition tripped and then resolved. `f24ad3c8` is **not** an
ancestor of the repo's LOCAL `develop`, which is 24 commits stale. It **is**
contained in `origin/develop`, and local `develop` is a strict fast-forward
ancestor of it (`git merge-base --is-ancestor develop origin/develop` → true),
so nothing is lost by basing on the remote ref. I based on
`origin/develop` = `7a662242`, which contains `f24ad3c8` and also carries the
brief-17 amendments (task 13 verify-only, task 16 + `fast-checks.yaml`). The
lead's addendum item C confirms this was pre-authorized. The Dropbox tree's
checked-out branch was never switched and nothing was built in it.

---

## 1. What was done

| # | Item | Commit | Note |
|---|---|---|---|
| 1 | F-087 inverted BF reading | `3d0ee1fa` | numbers now inline R |
| 4 | F-089 transition-count naming | `3d0ee1fa` | |
| 5 | F-092 second NA route | `3d0ee1fa` | uses report 15's observed example |
| 6 | F-094 full-chain ESS sentence | `3d0ee1fa` | |
| 2 | F-086 live sensitivity chunk | `07f72e94` | 56.0 s rendered |
| 2 | F-090 phantom report elements | `07f72e94` | |
| 3 | F-088 `\log_{10}` → natural log | `07f72e94` | |
| 7 | F-002 comparison-vignette sentence | `4146bad7` | |
| 8 | F-026 intro minimum (a)(b)(c) | `4146bad7` | |
| 9 | F-091 verdicts-print guard + 2 ride-alongs | `b207c990` | |
| 10 | F-093 order-of-magnitude messages | `64db6164` | + 2 comment sites |
| 11 | F-095 chain-boundary transitions | `0bd1cd50` | |
| 12 | F-031 two `\examples{}` | `c12e6bbe` | |
| 13 | F-049 Rd/NEWS harmonization | — | **verify only: they agree, no edit** |
| 14 | F-096 NEWS amendments (a)–(g) | `867db989` | per-entry verdicts in §2 |
| 15 | F-099/F-097 archive + remove | `866df819` (archive), `f3c1388c` (removal) | |
| 16 | F-102 CI `paths-ignore` | `311eb490` | |
| A | F-106 NEWS node encoding | `d847ca0d` | addendum |
| B | F-105 `format_log_bf()` negative zero | `d847ca0d` | addendum |

### Task 13 — the F-049 harmonization check (one line, as asked)

**They agree; nothing edited.** The Rd
(`R/prior_sensitivity.R:70-78` → `man/prior_sensitivity_check.Rd:143-151`) and
the NEWS sentence brief 14 merged (`NEWS.md:375-378`) both place exactness on
the anchor fits' own Rao-Blackwellized statistics — Rd: "read straight from each
fit's own Rao-Blackwellized statistics"; NEWS: "come from those fits' own
Rao-Blackwellized statistics, not from a pooled row" — with the between-anchor
curve importance-reweighted in both. That matches `R/anchor_curve.R:207-211`
("Exactness at anchors lives on the anchor fits' own RB statistics … not on
these pooled rows").

### Task 15 — archive-and-remove, in the briefed order

**(a) Archive first.** `archive/zratio-analytic-law` was branched from
`fix/docs-batch` at `867db989`, **before** any removal, and carries one commit
(`866df819`). On it: an archival note at the top of `zratio_law.h` (what was
removed with it, the maintainer's decision and date, the re-wiring condition,
and that regenerating `zratio_law_reference.rds` needs the private companion
implementation), the DORMANT banner kept, and every companion-referencing
comment rewritten self-contained. The brief named three sites; I did **seven**,
because four more `port of mu_law_solve_full / eval_mu_law / re_param2`
referents name private functions the same way and an archive branch that has to
stand alone should not keep them. Not pushed — the lead pushes at integration.

**(b) Removed from the package** (`f3c1388c`): `src/models/ggm/zratio_law.h`,
the `zratio_law_moments` block in `src/zratio_test_interface.cpp` (include line,
comment and function), `tests/testthat/test-zratio-law.R`, the fixture and its
generator. `Rcpp::compileAttributes()` re-run. **Tree-wide grep for
`zratio_law` over `src/`, `R/` and `tests/` returns nothing.**

**(c)** `src/models/ggm/zratio_engine.h:229` no longer says "matching the
companion's gold"; it states that the reference IS the per-component
Monte-Carlo oracle evaluation of the same decomposition the surface scores.

**(d)** One line added to `dev/review-2026-08/MAINTAINERS.md`, architecture
section, naming the branch and the re-wiring condition.

**(e)** No NEWS entry.

**(f) Deltas.** Source tree (`git archive`, `dev/` excluded, gzipped):
**4,333,508 → 4,323,565 bytes, −9,943 (−0.23 %)**; 721 lines of C++/R removed
(58,927 → 58,206 across `src/` + `tests/`). Built tarball is 1,625,738 bytes.
**Test time: 0 s at T0** — the whole file was T1/T2-gated, so the every-run tier
never ran it. At T1 report 12's block table costs it 50 s (fidelity) + 6 s (gate
guard) + 2 s (determinism) = **58 s**, plus 14 s at T2 (the all-MC oracle cell).
`test-zratio-law.R` leaving the tier classification is **authorized by F-099**,
not a silent drop.

---

## 2. Findings

### Finding 1 (needs the lead's eye) — two Bug-fixes entries the brief flagged as "verify" are dev-only, and I dropped them

The brief listed five ambiguous Bug-fixes entries to verify at the tag before
deciding. Three verified as genuine 0.1.6.3 defects and were kept. **Two did
not**, and dropping them is the only place I went beyond the brief's expected
outcome. Both verdicts are cited below; reverse either if the lead reads the tag
differently.

* **The category-scale recode fix — DROPPED.** The entry describes
  `simulate()` returning internal-scale data while `predict()` expects the
  original scale, producing a *"category values not observed in the training
  data"* warning, plus a mixed-model `predict()` path that miscodes silently. At
  the tag: the mixed path does not exist (`variable_type` takes only `"ordinal"`
  and `"blume-capel"`), that warning string does not exist anywhere
  (`git grep` over `cran-0.1.6.3 -- R/` is empty; it lives at
  `R/simulate_predict.R:1076` today), and the round trip is **self-consistent**.
  The tag's `recode_data_for_prediction()` recodes by *min-shift*, not by a
  stored map: Wenchuan is coded 1..5, the fit recodes to 0..4, `simulate()`
  returns 0..4, and `predict()` on that leaves it alone because its min is
  already 0 — while original-scale data shifts by 1 to the same place. Both
  paths land on the correct internal scale.
  **Caveat, stated so the drop is honest:** a *sparse* original coding (values
  with gaps, e.g. 1,2,4,5) would break the tag's min-shift, since the fit's
  recode collapses to 0..3 while the min-shift gives 0,1,3,4. That is a real
  tag-era defect — but it is not the defect the entry describes, and writing it
  up would be a new claim, not a rescue of this one.

* **The two imputation-cache fixes — DROPPED.** The entry is about a stale
  gradient cache and a stale observation transpose after imputation. Neither
  cache exists at the tag: `src/bgm/bgm_logp_and_grad.cpp:322,458` computes
  `const arma::mat obs_double_t = … .t()` fresh inside each call, and there is
  no gradient cache to invalidate. Both are state introduced by the 0.2.0 model-
  class refactor (`observations_double_t_`, `invalidate_gradient_cache()` in
  `src/models/omrf/omrf_model.cpp:1190-1198`) and fixed by `0db388db`, dated
  March 2026 — well after the tag.

### Finding 2 — F-095 changed no snapshot, but it did falsify two existing tests

The repository has exactly two snapshot files (`_snaps/plot-methods.md`,
`_snaps/verdicts.md`); **neither carries transition counts, so no snapshot was
re-recorded.** The re-record list the brief asked for is therefore empty, and
that is a verified emptiness, not an unchecked one. What the fix *did* break is
two live tests that had encoded the buggy pooled scan as their reference:

| test | before | after |
|---|---|---|
| `test-mcmc-diagnostics.R:447` "indicator ESS matches R reference" | R reference paired across the whole pooled vector | pairs within each chain and sums |
| `test-mcmc-diagnostics.R:531` "indicator ESS scales with multiple parameters" | counts sum to `n_total - 1` (999) | counts sum to `(niter-1) * nchains` (998) |

Both were corrected with the semantics, in the same commit. A reviewer should
read them as part of the fix, not as collateral.

### Finding 3 — the F-031 examples inherit uncapped cores, and I left them that way

The brief asked for the sibling extractors' pattern "exactly" and also for "no
uncapped cores" (F-018). These conflict: all 21 sibling examples call
`bgm(x = Wenchuan[, 1:3])` with no `cores` argument, inheriting
`cores = parallel::detectCores()`. I followed the siblings and did **not** add a
cap to only these two, because a lone capped pair would be the odd entry in a
23-example family and F-018 is a package-wide finding with its own owner. The
two examples measure 6.3 s together, in line with `extract_main_effects` (6.9 s)
and `extract_log_odds` (7.3 s). Flagged rather than decided.

### Finding 4 — F-093's wrong arithmetic also sat in two code comments

The brief scoped F-093 to the two runtime messages. The same "two orders below"
claim appeared in the two comments the messages were written from
(`R/zratio_surfaces.R:257`, `src/models/ggm/zratio_engine.h:199`). I corrected
those too — leaving the wrong arithmetic next to the fixed message would
reintroduce it at the next edit. No numbers changed, only the ratio's name.

### Finding 5 (addendum F-106) — the NEWS clause is wrong on this base *and* on the plot branch, but differently

Report 11 proposes replacement text describing a **probability wheel**. That is
the plot branch's behaviour, not this base's. Verified in my tree:
`main_difference_nodes()` (`R/plot_bgms.R:233-259`) fills **qgraph's pie
channel** to the difference indicator's posterior inclusion probability,
coloured by verdict, and returns no ring at all under
`main_difference_selection = FALSE`. I therefore wrote the clause against the
**ring**, which is what a reader of *this* branch would see, rather than adopting
the wheel text verbatim. **Coordination note for the lead:** when the plot branch
merges, this clause needs one more pass (ring → wheel). The Rd
(`R/plot_bgms.R:297-307`) already described the ring correctly — only NEWS was
ever wrong.

### Finding 6 — F-086's live chunk needed longer refits to render clean

The brief cited report 15's measured ~33 s for
`bgm(Wenchuan[, 1:6], chains = 2)` + `prior_sensitivity_check(fit)`. At those
settings the check runs in 25 s but the **1.6x anchor fails its convergence
gate**: the report emits a non-convergence notice, the mover table shows `<NA>`,
and the per-scale verdict counts read `0 0 0` at that anchor. Honest output, but
a poor headline demo, and the surrounding prose about the run-to-run noise band
would not match — the footer falls back to "no run-to-run yardstick".
Lengthening only the refits (`iter = 3000, warmup = 1000`) clears the gate at
**49.5 s measured / 56.0 s rendered**, gives all five anchors, two movers, and a
real noise-band footer. `cores = 2` is passed explicitly on both calls (the
default would be `spec$sampler$chains`, also 2 — pinned rather than inherited).

---

## 3. Evidence

### 3.1 Gate 1 — test suite, both tiers, on the final tree

Run with `devtools::test()`, which is what CI runs (`fast-checks.yaml`) and what
puts internals in scope; `NOT_CRAN=true`, `RCPP_PARALLEL_NUM_THREADS=4`,
sequential, one process at a time.

| tier | failed | warnings | errors | skipped | passed | wall |
|---|---|---|---|---|---|---|
| default (no slow env vars) | **0** | **0** | **0** | 101 | 8245 | 3.6 min |
| slow (`BGMS_RUN_SLOW_TESTS=true`) | **0** | **0** | **0** | 66 | 8437 | 5.3 min |

The slow tier is reported because `src/` changed (F-095). Both tiers were also
run before the addendum landed (8241 / 8433 passed, same zeros); the numbers
above are the final tree.

### 3.2 Gate 2 — vignette renders, sequential, one R process at a time

Every warning and message trapped at the render call with
`withCallingHandlers`, plus a post-hoc scan of the emitted HTML.

| vignette | seconds | warnings | messages | "Warning"/"Error" in HTML |
|---|---|---|---|---|
| intro | 10.1 | 0 | 0 | 0 |
| comparison | 6.2 | 0 | 0 | 0 |
| diagnostics | 10.2 | 0 | 0 | 0 |
| checking-your-model | 16.6 | 0 | 0 | 0 |
| **prior-sensitivity** | **56.0** | 0 | 0 | 0 |
| **total** | **99.1** | | | |

**Budget statement, as asked.** prior-sensitivity went 0.4 s → 56.0 s, and the
five-vignette total went ~46 s → 99.1 s. Inside `R CMD check`, vignette
re-building measures `[195s/102s]` — the elapsed figure is the one to read, and
it is comfortably inside the check. The cost buys a report that cannot go stale
again.

### 3.3 Gate 3 — `R CMD check --as-cran` on a `git archive` tarball

Built and checked outside both the Dropbox tree and the worktree, with RStudio's
bundled pandoc on `PATH` (F-045), `_R_CHECK_LIMIT_CORES_=TRUE`. Run twice — once
before the addendum and once on the final tree, since the addendum touched
`R/verdicts.R` and `NEWS.md`. Both:

```
Status: 2 NOTEs
```

Both are the recorded baseline pair, verbatim:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Maarten Marsman <m.marsman@uva.nl>'
The Date field is over a month old.

* checking HTML version of manual ... NOTE
Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy.
```

No ERROR, no WARNING. Supporting lines from the final-tree run:
`checking examples ... OK`,
`checking examples with --run-donttest ... [579s/295s] OK`,
`checking tests ... [112s/84s] OK`,
`checking re-building of vignette outputs ... [195s/102s] OK`. Tarball
1,626,038 bytes.

### 3.4 Gate 4a — F-087, the printed numbers and the prose, from the rendered HTML

**Before** (`origin/develop:vignettes/diagnostics.Rmd:146,152-153`):

> Here the Bayes factor in favor of inclusion (H1) is **small, meaning that
> there is little evidence for inclusion**. […] This Bayes factor shows that
> there is **strong evidence for the absence** of a network relation between the
> variables `intrusion` and `physior`.

— against a chunk printing `BF_10 = 24.59443`. Both readings were the reciprocal
of the truth.

**After**, quoting the rendered `diagnostics.html` — chunk output and prose
together:

```
BF_10
#> [1] 24.59443

Here the Bayes factor in favor of inclusion (H1) is 24.59: the data are about
25 times more likely under a model that includes the intrusion-physior edge
than under one that excludes it, which is substantial evidence for inclusion.
On the natural log scale that the package uses elsewhere this is a log Bayes
factor of 3.2, and verdicts() — see the model-checking vignette — labels this
edge a presence.

1 / BF_10
#> [1] 0.04065962

This is the Bayes factor for absence, BF_01, and at 0.041 it sits well below 1:
it is evidence against absence rather than for it — the mirror image of the
reading above. Absence is what a large BF_01 would show: a value of, say, 10 or
more, which is the same as a BF_10 of 0.1 or less, would be strong evidence
that there is no network relation between two variables.
```

`24.59` / `25` / `3.2` / `0.041` are `r round(BF_10, 2)`,
`r round(BF_10)`, `r round(log(BF_10), 2)` and `r signif(1 / BF_10, 2)` — the
prose cannot disagree with the chunk again. `3.2` and `presence` agree with
report 15's `verdicts()` reading (`pip 0.961`, `log_bf 3.20`, `presence`). The
1/BF sentence is kept, as briefed, as the correct illustration of what a *large*
`BF_01` would mean.

### 3.5 Gate 4b — F-086, live output, no pasted block

`grep` for the old mock's marker strings over `vignettes/` returns nothing: no
`#> Prior sensitivity check:` fence, no "Details:", no "not certifiable (chains
disagree)". What the rendered `prior-sensitivity.html` now carries is the shipped
print method's own output, produced by the visible chunk above it:

```
fit = bgm(data[, 1:6], seed = 1234, chains = 2, cores = 2,
          display_progress = "none", verbose = FALSE)
ps = prior_sensitivity_check(fit, seed = 1234, iter = 3000, warmup = 1000,
                             cores = 2)
ps
#> Prior sensitivity check: are the edge verdicts robust to the slab scale?
#> …
#>   robust (same verdict at every scale)     13
#>   changed, within run-to-run noise          0
#>   changed, beyond run-to-run noise          2
#> …
#> Method:  43-point curve from 5 anchor fits (0.4x to 2.5x the chosen scale),
#>          joined by importance reweighting; the 1x anchor is the original fit.
#> Refits:  5 nuts refits, warm-started from the original fit, 44 s total.
#> Noise:   two identical refits at 1.6x differed by up to …
```

The four divergences report 15 tabulated are gone by construction. F-090's two
sentences were rewritten against this output: the footer is named as the
`Noise:` line (there is no "Details line"), and the uncertifiable wording is
quoted as `"not certifiable (too noisy to assess)"`
(`R/prior_sensitivity.R:1005`).

### 3.6 F-092 — the second NA route, in the vignette's own rendered table

```
#>                         mean        mcse           sd    n_eff      Rhat n0->1
#> intrusion-dreams  1.00000000          NA 0.000000e+00       NA        NA     0
#> intrusion-flash   1.00000000          NA 2.206767e-11       NA 0.9999966     0
```

`intrusion-flash`: non-zero `sd`, `NA` `mcse`/`n_eff`, **finite** `Rhat` — the
route `R/mcmc_summary.R:268-270` names ("variance falls below the autocovariance
kernel's numerical floor"), distinct from `intrusion-dreams`'s exactly-constant
chain. The vignette now names both and points at this row.

### 3.7 F-091 — before and after

```r
v = verdicts(fit)
subset(v, fragile)
# before: Error in log(threshold) : non-numeric argument to mathematical function
# after:  prints the three fragile rows as a plain data frame
```

Ride-alongs, from the same print: `(1 indicators)` → `(1 indicator)` (verified
on a one-row table), and `Bayes factor of 30 (and 0.0333333 for absence)` →
`(and 0.0333 for absence)`. The existing `_snaps/verdicts.md` fixtures use
threshold 10 and 3/6 rows, so neither snapshot moved.

### 3.8 F-095 — the pinned case

```
2 chains x 50 draws; chain 1 all zeros, chain 2 all ones
before: n0->1 = 1   (the chain boundary counted as a transition)
after:  n0->1 = 0,  n1->0 = 0,  n00 = 49, n11 = 49, n_eff_mixt = NA, mean = 0.5
```

Plus a 3-chain alternating case pinning that the counts are within-chain totals
summed (`n01 = 6`, `n10 = 3` over three chains of `0,1,0,1`).

### 3.9 F-096 — per-entry tag verdict table

All verdicts established at `cran-0.1.6.3` (= `18e660a2`).

**Lead-supplied facts, independently re-verified before use:** both released
samplers use the old scale with no factor 2 —
`src/bgmCompare/bgmCompare_sampler.cpp:118-120` builds `rest_score` off plain
`observations.row(person) * group_pairwise_effects.col(variable)` and adds
`category * rest_score`; `src/bgm/bgm_sampler.cpp:593` builds
`residual_matrix = obs_double * pairwise_effects`;
`src/mrf_simulation.cpp:98` accumulates `(obs - ref) * pairwise_safe(...)`. The
released pair was **mutually consistent**. The tag's `NAMESPACE` exports
`S3method(predict, …)` / `S3method(simulate, …)` for both classes plus
`export(simulate_mrf)` and `export(mrfSampler)`.

| Bug-fixes entry | evidence at the tag | verdict |
|---|---|---|
| category-scale `simulate()`/`predict()` recode | mixed path absent; warning string absent (`git grep` empty); tag recodes by min-shift (`R/simulate_predict.R:1238-1252`) so the round trip is self-consistent | **DROP** (see §2 finding 1) |
| mixed-MRF PIP block order | mixed models do not exist: tag `variable_type` allows `"ordinal"` / `"blume-capel"` only (`R/bgm.R`, `@param variable_type`) | **DROP** |
| SBM number-of-blocks summary | sampler uses shifted Poisson `R::dpois((k+1)-1, lambda, true)` (`src/sbm_edge_prior_interface.cpp`, `compute_Vn_mfm_sbm`); summary uses zero-truncated `dpois(K, lambda)/(1 - dpois(0, lambda))` (`R/mcmc_summary.R:409-410`) — genuine mismatch | **KEEP**, reframed |
| cross-indicator asymmetry (mixed + SBM) | mixed absent | **DROP** |
| `delta = NULL` mixed default | mixed absent; determinant tilt is 0.2.0 (`cb68f0ea`) | **DROP** |
| Cholesky downdate (GGM/mixed) | continuous models absent | **DROP** |
| Alpine/musl include | `tbb::global_control` used at `src/mrf_simulation.cpp:454` with no `<tbb/global_control.h>` include | **KEEP**, 0.1.6.3 framing |
| stale gradient cache after imputation | no cached transpose, no gradient cache at the tag (`bgm_logp_and_grad.cpp:322,458` recompute per call); fixed by `0db388db`, Mar 2026 | **DROP** |
| stale observation transpose after imputation | same | **DROP** |
| NUTS `target_accept` pass-through | user's value reaches dual averaging (`bgm_sampler.cpp:1361-1364`) but the post-mass-update step-size heuristic hard-codes `0.625` (`:607`, `:755`); current code passes `target_acceptance_` (`src/mcmc/samplers/nuts_sampler.h:86`) | **KEEP**, narrowed to what is true |
| NUTS acceptance accumulation | tag overwrites: `alpha = result.alpha` (`src/mcmc/mcmc_nuts.cpp:367-368`); current sums: `sum_metro_prob += result.alpha` (`src/mcmc/algorithms/nuts.cpp:353`) | **KEEP**, 0.1.6.3 framing |

Amendments (a)–(f) applied as briefed. (c) confirmed by the tag's `NAMESPACE`
and left as merged. (f) mirrored in the Rd at
`R/extractor_functions.R` `@return` for `extract_inclusion_bf` — a natural slot
existed ("Entries are `NA` for indicators that were never updated"), so the
"never updated means never *proposed*" clarification sits directly under it.

**Parse gate** (report 14's tool, `tools:::.build_news_db_from_package_NEWS_md`),
re-run after all NEWS edits including the addendum:

```
rows: 33   versions: 12   NA version/category rows: 0
0.2.0.0 sections:  (preamble) 1 | Breaking changes 1 | New features 1 |
                   Other changes 1 | Bug fixes 1 | Deprecated 1
```

Identical to report 14's recorded shape.

### 3.10 F-105 — the negative-zero guard

```
format_log_bf(-0.004)   before: "= -0.0"    after: "= 0.0"
format_log_bf(-0.04)    before: "= -0.0"    after: "= 0.0"
format_log_bf(-0.06)    before: "= -0.1"    after: "= -0.1"   (unchanged)
```

The guard is inside `format_log_bf()` (`R/verdicts.R:113-119`), so
`print.bgms_verdicts()` and the summary tables — which reach it directly — are
covered. **No plot file was touched.** The rule and its comment mirror
`estimate_lines()` (`R/plot_bgms.R:989-992`), which had already decided this
question for a weight.

### 3.11 F-102 — the four workflows

`paths-ignore: ['dev/**']` added to the `on: push` trigger of `lint.yaml`,
`R-CMD-check.yaml`, `test-coverage.yaml` and `fast-checks.yaml`.
`fast-checks.yaml`'s `pull_request:` trigger is byte-identical to before.
`nightly-validation.yaml` and `weekly-certification.yaml` untouched
(schedule-triggered; no push trigger). All nine workflow files re-parse under
`yaml::read_yaml`.

---

## 4. Open questions

1. **The two dropped Bug-fixes entries (§2 finding 1).** These are the batch's
   only judgement calls that go past the brief's expected outcome. If the lead
   or MM reads the tag differently on either, both are one revert away. The
   sparse-coding caveat under the category-scale entry may deserve its own row
   if the maintainer wants tag-era `predict()` behaviour on record.

2. **F-106 will need a second pass.** The clause I wrote describes the pie ring,
   which is what this base draws. When the plot branch (wheel) merges, the same
   sentence needs ring → wheel. Whoever integrates that branch should own it;
   report 11's proposed text is already correct *for that branch*.

3. **F-031 and F-018 (§2 finding 3).** The two new examples inherit uncapped
   cores like their 21 siblings. If F-018's fix is expected to reach examples
   rather than only the sampler path, these two should ride along with the other
   23 rather than be capped alone.

4. **The vignette budget.** prior-sensitivity is now 56 s of the five-vignette
   99 s. It is inside the check, but it makes that one vignette the dominant
   cost of `R CMD build`. If the maintainer wants it cheaper, the lever is the
   refit length (`iter = 3000, warmup = 1000`): dropping back to the defaults
   costs 25 s instead of 50 s but reintroduces the failed 1.6x anchor described
   in §2 finding 6.

5. **`plot(ps)` and `head(ps$edges)` in prior-sensitivity.Rmd remain
   `eval = FALSE`,** even though `ps` is now a live object and the plot would be
   nearly free. Left alone deliberately: the brief scoped F-086 to the report
   block, and F-067 is going to change that plot's title. Cheap to switch on
   later.

---

## 5. Branch state

* `fix/docs-batch` — 11 commits on top of `origin/develop` (`7a662242`). Not
  pushed.
* `archive/zratio-analytic-law` — branched from `fix/docs-batch` at `867db989`
  before the removal; 1 commit (`866df819`). **Not pushed**, per the brief; the
  lead pushes it at integration.
* Nothing was built in the Dropbox tree, and its checked-out branch was not
  switched.

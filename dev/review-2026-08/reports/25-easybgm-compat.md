# Report 25 — easybgm compatibility: what bgms 0.2.0.0 changes downstream (Opus agent)

Report only. Nothing in `R/`, `src/`, `tests/`, `man/`, `vignettes/`, or
`NEWS.md` was touched, and **easybgm was not modified** — the CRAN tarball was
unpacked read-only for source citation and the CRAN binary installed unchanged
into a private library.

**Headline (the ordering answer, one sentence):** bgms 0.2.0.0 does **not**
break the CRAN easybgm — easybgm 0.4.0's test suite passes 0-failure /
0-error against the new build, every workflow that runs on 0.1.6.3 also runs
on 0.2.0.0, and `R CMD check --as-cran` returns the **identical** status
(1 WARNING, 2 NOTEs, byte-identical text, all three artifacts of re-checking a
published tarball) against both builds — so the submission order is
**bgms first, with a notification to the easybgm maintainer**.

**Second headline:** the one thing CRAN's reverse-dependency check *will*
notice is time, not correctness. easybgm's `R CMD check` goes from **83 s to
196 s (2.36×)**, its `--run-donttest` examples pass from **[78 s CPU / 25 s
elapsed] to [252 s / 128 s]**, and its two heaviest examples from **9.1 s to
58 s elapsed each**. Measured, not inferred; the BGGM-backed examples in the
same check are unchanged to within 1%, which is what makes the rest
attributable to bgms.

**Third headline:** the shim's gate is evaluated when a fit is *built*, not
when it is *read*, so a fit made before `library(easybgm)` is an S7 object
that easybgm then mishandles (**F-025-1**). The live run settles the question
the severity turns on: it fails **LOUDLY** — every downstream call errors with
`this object class is not subsettable`, and no wrong number can escape. That
drops F-025-1 to LOW.

---

## 0. Environment

```sh
git -C <dropbox>/bgms worktree add ~/bgms-review/wt-val25 \
    -b review/easybgm-compat origin/develop
R CMD INSTALL --library=~/bgms-review/lib-val25   ~/bgms-review/wt-val25
R CMD INSTALL --library=~/bgms-review/lib-anchor25 ~/bgms-review/val25/src/anchor/bgms
```

| item | value |
|---|---|
| bgms commit | **`8f28d49e`** (`origin/develop`, contains `a37c045d`) |
| bgms version under test | **0.2.0.0** — library `~/bgms-review/lib-val25` |
| bgms anchor | **0.1.6.3** (the version CRAN easybgm was built against) — library `~/bgms-review/lib-anchor25`, built from `~/bgms-review/bgms_0.1.6.3.tar.gz` |
| easybgm version | **0.4.0** (CRAN, published 2026-04-02; binary `easybgm_0.4.0.tgz` + source `easybgm_0.4.0.tar.gz`) |
| easybgm maintainer | Karoline Huth `<k.huth@uva.nl>` |
| R | 4.6.0 (2026-04-24), darwin/aarch64, 15 cores |
| scratch | `~/bgms-review/val25/` |

Both libraries were verified to resolve exactly one bgms each
(`find.package("bgms")`), sharing the one easybgm 0.4.0 installation.

**Machine budget, honoured.** A lead measurement job (`20_phaseA_normal.R`,
PID 3345, ~400% CPU) ran from 14:09 to **18:23**. Tasks 1 and 2, the fit-free
gate demonstration, both library builds, and all script staging were done while
it ran; **no fit was started until 18:23:19**. All twelve runs then executed
strictly sequentially (`run-all.sh`), one job at a time, finishing 18:29:46.
The per-fit timing runs pin `cores = 4`; the two `R CMD check` runs are
unpinned, because that is what CRAN measures.

**Scripts and logs** (all in `~/bgms-review/val25/`): `01-shim.R` (gate, three
states), `02-workflows.R` (14 workflows × 2 builds), `03-testsuite.R` (easybgm
test suite × 2 builds), `04-trap.R` (F-025-1), `05-timing.R` (13 timed fits ×
2 builds), `06-check.sh` (`R CMD check --as-cran` × 2 builds), `run-all.sh`,
plus the matching `.log`, `.rds`, and `timing-merged.csv`.

easybgm file:line citations are to the **CRAN 0.4.0 source tarball** at
`~/bgms-review/val25/src/easybgm/`; anchor citations to
`~/bgms-review/val25/src/anchor/bgms/`; bgms citations to the worktree at
`8f28d49e`.

---

## 1. The compatibility shim, mapped then exercised

### 1.1 What the gate is

```r
needs_easybgm_s3_compat = function() {
  if(!"easybgm" %in% loadedNamespaces()) {
    return(FALSE)
  }
  ebgm_version = utils::packageVersion("easybgm")
  if(ebgm_version < "0.5.0") {
    warning(
      "easybgm ", ebgm_version, " is not compatible with S7-based bgms objects. ",
      "Running in S3 compatibility mode. ",
      "Please update easybgm to version 0.5.0 or later.",
      call. = FALSE
    )
    return(TRUE)
  }
  FALSE
}
```

[build_output.R:317-332](R/build_output.R#L317-L332). Called from exactly
three places, each the last statement of an output builder:

* [build_output_bgm.R:383-387](R/build_output_bgm.R#L383-L387) (GGM + OMRF)
* [build_output_mixed_mrf.R:355](R/build_output_mixed_mrf.R#L355) (mixed MRF)
* [build_output_compare.R:240-244](R/build_output_compare.R#L240-L244) (`bgmCompare()`)

each reading `if(needs_easybgm_s3_compat()) results else s3_list_to_bgms(results)`.

### 1.2 The three states and the object shape each hands easybgm

The gate has no third branch in code — `>= 0.5.0` and "not loaded" both return
`FALSE` — but they differ at the point of the warning, so all three were
exercised. State C needed an easybgm at a version CRAN does not carry: a
**stub package** named `easybgm` version 0.5.0, containing no easybgm code,
was built into a *separate* library (`~/bgms-review/lib-stub`, source at
`~/bgms-review/val25/stub/easybgm/`) and put first on `.libPaths()` for that
one session. The CRAN easybgm installation was never touched.

| state | gate | warning | object returned |
|---|---|---|---|
| **A** easybgm not loaded | `FALSE` | none | **S7**, `class = c("bgms", "S7_object")`, `typeof = "object"` |
| **B** easybgm **0.4.0** loaded | `TRUE` | the four-sentence warning, once per fit, `call. = FALSE` | **plain S3 list**, `class = "bgms"`, `typeof = "list"` |
| **C** easybgm **0.5.0** (stub) | `FALSE` | none | **S7**, identical to A |

In state B the `posterior_summary_*` fields are literal `NULL` placeholders in
the list ([build_output_bgm.R:322-326](R/build_output_bgm.R#L322-L326),
[build_output_compare.R:241-247](R/build_output_compare.R#L241-L247)) so that
the *names* appear in `names(fit)` — that is what the comment at
[build_output_bgm.R:319](R/build_output_bgm.R#L319) means by "for easybgm
compat" — while `` `$.bgms` `` intercepts any name starting
`posterior_summary_`, calls `ensure_summaries()`, and returns the cached value
instead of the `NULL` ([methods_bgms.R:303-319](R/methods_bgms.R#L303-L319));
`` `[[.bgms` `` does the same ([methods_bgms.R:325-345](R/methods_bgms.R#L325-L345)).
This is the mechanism easybgm depends on at `functions.bgms.R:250` and `:255`.

### 1.3 The deprecated easybgm-compatibility section of the S7 class

[class_s7.R:92-95](R/class_s7.R#L92-L95):

```r
    # --- easybgm compatibility (deprecated) ---
    indicator = new_property(class_any, default = NULL),
    interactions = new_property(class_any, default = NULL),
    thresholds = new_property(class_any, default = NULL),
```

Three properties on `bgms_class` only (`bgmCompare_class` has no such
section), populated by `s3_list_to_bgms()`
([class_s7.R:134-136](R/class_s7.R#L134-L136)). They exist so a *legacy* fit
from bgms 0.1.4–0.1.5 — where these were the top-level sample fields —
survives S7 conversion and stays readable by the deprecation branches of the
extractors ([extractor_functions.R:300-306](R/extractor_functions.R#L300-L306),
[:924-930](R/extractor_functions.R#L924-L930),
[:1002-1006](R/extractor_functions.R#L1002-L1006)).

**They are not what easybgm 0.4.0 reads.** A grep of the whole tarball finds no
read of `$indicator`, `$interactions`, or `$thresholds` on a bgms object;
easybgm's own result list has similarly named fields (`bgms_res$thresholds`,
`functions.bgms.R:149`) but never reads them off a fit. A 0.2.0.0 fit sets all
three to `NULL`. With respect to CRAN easybgm the section is dead weight
(**F-025-4**).

### 1.4 The gate exercised live — verbatim

Fit-free gate call, all three states (`01-shim.R`):

```
############ STATE: none ############
bgms          : 0.2.0.0
easybgm loaded: FALSE
gate returns  : FALSE
############ STATE: cran ############
bgms          : 0.2.0.0
easybgm loaded: TRUE
easybgm ver   : 0.4.0
WARNING VERBATIM: easybgm 0.4.0 is not compatible with S7-based bgms objects. Running in S3 compatibility mode. Please update easybgm to version 0.5.0 or later.
gate returns  : TRUE
############ STATE: stub ############
bgms          : 0.2.0.0
easybgm loaded: TRUE
easybgm ver   : 0.5.0
gate returns  : FALSE
```

With a live `bgm()` fit behind it — **state B, the CRAN easybgm**
(`01-shim-cran.log`):

```
--- gate value ---------------------------------------------------
GATE WARNING: easybgm 0.4.0 is not compatible with S7-based bgms objects. Running in S3 compatibility mode. Please update easybgm to version 0.5.0 or later.
[1] TRUE

--- bgm() fit ----------------------------------------------------
FIT WARNING: warmup = 200: limited proposal SD tuning. Consider >= 300.
FIT WARNING: easybgm 0.4.0 is not compatible with S7-based bgms objects. Running in S3 compatibility mode. Please update easybgm to version 0.5.0 or later.

class(fit)        : bgms
typeof(fit)       : list
is S7             : FALSE
length(names(fit)): 13

--- extractor contract on this object -----------------------------
extract_arguments(fit)                         OK   list (length 32)
extract_pairwise_interactions(fit)             OK   matrix (200 x 6)
  [warn] `extract_category_thresholds()` was deprecated in bgms 0.2.0.
ℹ Please use `extract_main_effects()` instead.
extract_category_thresholds(fit)               OK   matrix (4 x 4)
extract_posterior_inclusion_probabilities(fit) OK   matrix (4 x 4)
extract_indicators(fit)                        OK   matrix (200 x 6)
fit$posterior_summary_pairwise                 OK   data.frame (6 x 6)
fit$posterior_summary_indicator                OK   data.frame (6 x 7)
fit$arguments                                  OK   list (length 32)

--- what easybgm 0.4.0 does to the object -------------------------
after class(fit) <- "bgms":  class = bgms
  extract_arguments()$no_variables present : TRUE
  extract_pairwise_interactions() dim      : 200 x 6
```

**State A / C, the S7 object, same script, same statement** (`01-shim-none.log`,
`01-shim-stub.log`):

```
class(fit)        : bgms, S7_object
typeof(fit)       : object
is S7             : TRUE
...
--- what easybgm 0.4.0 does to the object -------------------------
after class(fit) <- "bgms":  ERROR: this object class is not subsettable
```

That contrast is the whole shim in two lines: the same easybgm statement is a
no-op on the S3 list and fatal on the S7 object.

### 1.5 The designed behaviour, stated plainly

The brief anticipated that a sub-0.5.0 CRAN easybgm might meet an "instructive
error". **It does not, and the report's map is authoritative here.** The
easybgm < 0.5.0 path is a `warning()` followed by S3 compatibility mode: bgms
deliberately keeps old easybgm *working* rather than refusing to serve it. The
designed behaviour for easybgm 0.4.0 is one warning per fit, naming the
version and the required version, plus a plain S3 list that every easybgm
accessor consumes exactly as it did under 0.1.6.3 — demonstrated above and
confirmed end-to-end in §3. That is the correct design: a hard error would
take the CRAN easybgm offline the day bgms 0.2.0.0 landed.

### 1.6 F-025-1 — the gate fires at fit-construction time (LOUD)

`needs_easybgm_s3_compat()` runs inside the output builders, so the object
shape is fixed when `bgm()` returns. A session that fits *first* and loads
easybgm *second* holds an S7 object that easybgm then destroys. This is a
supported user order — bgms's own print methods send users to easybgm
([methods_bgms.R:74](R/methods_bgms.R#L74),
[:140](R/methods_bgms.R#L140), [:231](R/methods_bgms.R#L231);
[methods_bgmcompare.R:194](R/methods_bgmcompare.R#L194),
[:271](R/methods_bgmcompare.R#L271),
[:421](R/methods_bgmcompare.R#L421)).

`04-trap.R` runs that order in one process. Verbatim:

```
--- step 1: fit with easybgm NOT loaded --------------------------------
class(fit)     : bgms, S7_object
typeof(fit)    : object
is S7          : TRUE
reference: mean |pairwise| = 0.184612;  pip[1,2] = 1

--- step 2: NOW load easybgm 0.4.0 -------------------------------------
easybgm version: 0.4.0
gate would now return: TRUE  -- but `fit` was built before this point.
class(fit) is still: bgms, S7_object

--- step 3: the exact statement easybgm executes -----------------------
functions.bgms.R:54 / :61  ==>  class(fit) <- "bgms"

>>> class(trapped) <- "bgms"
    OUTCOME: NO ERROR (returned bgms)
    class(trapped) after: bgms
    typeof(trapped)     : object
    length(unclass-len) : 1

--- step 4: what the corrupted object then does ------------------------
>>> extract_arguments(trapped)                        ERROR: this object class is not subsettable
>>> trapped$arguments                                 ERROR: this object class is not subsettable
>>> extract_pairwise_interactions(trapped)            ERROR: this object class is not subsettable
>>> extract_posterior_inclusion_probabilities(trapped) ERROR: this object class is not subsettable
>>> trapped$posterior_summary_pairwise                ERROR: this object class is not subsettable

--- step 5: the whole easybgm entry point, end to end ------------------
>>> easybgm:::bgm_extract.package_bgms(fit = <S7 fit>, ...) ERROR: this object class is not subsettable
>>> plot_network(<S7 fit>)   [plot_network.bgms]            ERROR: this object class is not subsettable
>>> plot_edgeevidence(<S7 fit>)                             ERROR: this object class is not subsettable

--- step 6: VERDICT ----------------------------------------------------
VERDICT: LOUD -- at least one easybgm call errors; the user cannot get a
         silently wrong number out of this path.
```

**The class assignment itself succeeds silently** — `class(trapped) <- "bgms"`
raises nothing — but the object's base type is `object` (R's OBJSXP, which S7
uses), not `list`, so the very next `.subset2()` in `` `$.bgms` `` fails. Every
one of the eight downstream calls errors with the same message. **No corrupted
number can reach a user.** The failure is annoying and its message is unhelpful
(it names neither easybgm nor the load order), but it is not silent, which is
what the severity turned on. **F-025-1: LOW.**

---

## 2. The delta enumeration — every 0.2.0.0 change that reaches easybgm

### 2.1 The complete bgms API surface easybgm 0.4.0 touches

Exhaustive, by grep over the tarball. Anything not listed is **NOT REACHED**
by construction.

| what | easybgm call sites | count |
|---|---|---|
| `bgm()` | `functions.bgms.R:14-19`, `:21-26` | 2 |
| `bgmCompare()` | `functions.bgmscompare.R:20-25`, `:37-42`, `:53-59` | 3 |
| `extract_arguments()` | 21 sites (17 as `bgms::extract_arguments`) | 21 |
| `extract_indicators()` | `functions.bgms.R:185`, `:230`, `:252`; `functions.bgmscompare.R:212`, `:216`, `:218`, `:290`, `:295`, `:297` | 9 |
| `extract_pairwise_interactions()` | `functions.bgms.R:146`, `:148`, `:194`; `functions.bgmscompare.R:233`, `:300` | 5 |
| `extract_group_params()` | `functions.bgmscompare.R:229`, `:230`, `:231`, `:309` | 4 |
| `extract_posterior_inclusion_probabilities()` | `functions.bgms.R:153`, `:201` | 2 |
| `extract_category_thresholds()` | `functions.bgms.R:149`, `:196` | 2 |
| `extract_sbm()` | `functions.bgms.R:77` | 1 |
| `$` on a fit | `$arguments` `functions.bgms.R:56`; `$posterior_summary_pairwise` `:250`; `$posterior_summary_indicator` `:255`; `$posterior_summary_pairwise_differences` `functions.bgmscompare.R:206`, `:235`, `:272`; `$posterior_summary_pairwise_baseline` `:284`, `:308` | 8 |
| datasets | `Wenchuan`, `ADHD` | 2 |

**Not touched at all:** `extract_ess()`, `extract_inclusion_bf()`,
`extract_rhat()`, `extract_main_effects()`,
`extract_prior_inclusion_probabilities()`, `extract_precision()`,
`extract_partial_correlations()`, `extract_log_odds()`,
`extract_edge_indicators()`, `extract_centrality()`, `verdicts()`,
`calibration_check()`, `prior_sensitivity_check()`, and
`plot()`/`summary()`/`coef()`/`predict()`/`simulate()` on a bgms fit,
`simulate_mrf()`, `mrfSampler()`, `sample_*_prior()`,
`plot_edge_posterior()`, `summarize_zratio_gauge()`, every prior constructor.

**Version gates in easybgm:** `< 0.1.3`, `< 0.1.4`, `> 0.1.4.2`, `< 0.1.6`,
`< 0.1.6.0`, `> 0.1.6` — and **nothing above `0.1.6`**. Every 0.2.0.0 fit
therefore falls into the same branch a 0.1.6.3 fit did. easybgm has no place
to notice that the release changed.

### 2.2 The classified delta table

Anchor: bgms 0.1.6.3. **BREAKS** = easybgm errors or produces structurally
wrong output; **CHANGES NUMBERS** = same code runs, different numbers a user
or test would see; **COSMETIC** = warning/label only; **NOT REACHED** = no
easybgm code path arrives.

#### BREAKS — none

None on any path easybgm executes. The shim removes the only candidate (the
S7 class) before easybgm sees it, and every extractor easybgm calls still
exists with the same S3 dispatch and the same `$arguments` field names.
Executed proof in §3.1–§3.4; the one reachable configuration where the shim
does not fire is F-025-1 (§1.6), and it fails loudly.

#### CHANGES NUMBERS

| # | 0.2.0.0 change | easybgm site | what an easybgm user sees |
|---|---|---|---|
| N1 | Pairwise parameters on the **association scale** — `extract_pairwise_interactions()` and the raw pairwise draws are ~half their 0.1.6.3 values | `functions.bgms.R:146-148`, `:194-195` | `res$parameters`, `res$samples_posterior`, every edge weight in `plot_network`/`plot_parameterHDI`, and `res$centrality` (strength = row sums of `abs(samples_posterior)`, `AuxiliaryFunctions.R:100-112`) all halve |
| N2 | Same for `bgmCompare()`, incl. `extract_group_params()` | `functions.bgmscompare.R:207`, `:229-231`, `:233`, `:282`, `:284`, `:300`, `:309` | `res$parameters`, `parameters_g1`, `parameters_g2`, `group_estimates`, `overall_estimate`, `samples_posterior` all halve |
| N3 | `bgm()` interaction prior **Cauchy(2.5) → `normal_prior(1)`** | no site — easybgm never passes `interaction_scale`, so it inherits the default | a default `easybgm(type="binary", package="bgms")` is now a **different model**. easybgm's documentation (`easybgm.R:86-87`, "Cauchy … default is 2.5") is now factually wrong (**F-025-2**) |
| N4 | `bgmCompare()` baseline `interaction_prior = normal_prior(1)` **and** `difference_family = "Normal"` (both Cauchy before) | inherited | every `easybgm_compare(..., package="bgms")` result changes; difference verdicts on categories one group never observed move most |
| N5 | `iter`/`warmup` **1e3 → 2e3** | easybgm passes `iter` (`functions.bgms.R:15`, `functions.bgmscompare.R:38`/`:54`) but **never `warmup`** | a different chain, and the dominant term in the runtime result of §3.5 |
| N6 | NUTS `target_accept` **0.60 → 0.80** (`bgm`), **0.65 → 0.80** (`bgmCompare`) | inherited | smaller steps, more leapfrogs; compounds N5 |
| N7 | **Rao-Blackwellized inclusion probability is canonical** — `extract_posterior_inclusion_probabilities()` defaults to `estimator = "rb"` ([extractor_functions.R:268-296](R/extractor_functions.R#L268-L296)) | `functions.bgms.R:153`, `:201` | `res$inc_probs` changes; `res$inc_BF` (hand-computed `:155-178`, `:203-226`) with it; `res$structure = 1*(inc_probs > 0.5)` (`:184`) can flip an edge. Most visible effect: raw indicator averages saturate at 0/1 on a short chain and gave `inc_BF` of `0`/`Inf`; RB averages do not, so easybgm's Bayes factors are now **finite where they used to be infinite** |
| N8 | Inclusion table's `n_eff` is now the **RB** ESS; 0.1.6.3's was transition-based `T/tau_int` (anchor `mcmc_summary.R:73`, `:95`) | `functions.bgms.R:255` → `BF_MCSE(ess = ...)` | `res$MCSE_BF` changes — and worse than changes: **F-025-3**. `BF_MCSE` now mixes an RB `ess` with a **raw** `p_hat` recomputed from `gamma_mat` (`AuxiliaryFunctions.R:286`, `:322-326`) and an RB-derived `BF_vec`, where under 0.1.6.3 all three were one estimator. RB columns are `NA` where the RB draws are constant; `BF_MCSE` degrades those to `NA` CI rows (`AuxiliaryFunctions.R:357-358`) rather than erroring — confirmed live, `MCSE_BF` came back `6 x 2` with `anyNA=TRUE` (§3.1) |
| N9 | **R-hat is classic split-R-hat**; the Brooks-Gelman adjustment that pinned near-saturated indicators at `sqrt(5/3) ≈ 1.29` is gone | `functions.bgms.R:250`; `functions.bgmscompare.R:235`, `:308` | `res$convergence_parameter` moves slightly everywhere, substantially on decisive edges; `NA` when chains are identical, `+Inf` for chains stuck at different constants |
| N10 | **`bgmCompare()` ordinal-category union fix** — 0.1.6.3 kept only categories observed in *every* group and merged the rest | every `easybgm_compare(..., package="bgms")` call | (i) `res$parameters`, `parameters_g1/g2` change, often a lot, since the old rule overestimated the affected variable's pairwise parameters and inflated its neighbours'; (ii) `num_categories` changes ⇒ the `(main)`/`(diff)` column count of `extract_indicators()` changes ⇒ the width of the `structures` strings pasted at `functions.bgmscompare.R:216`/`:295` ⇒ `structure_probabilities`, `graph_weights`, `sample_graph` change; (iii) a `message()` about renumbering and a classed warning **`bgms_group_support_warning`** ([validate_data.R:519-537](R/validate_data.R#L519-L537)) reach the user — `easybgm_compare`'s `tryCatch` (`easybgm_compare.R:184-193`) catches errors only. *Does not fire in easybgm's own examples: `ADHD[1:10,1:3]` and `[11:20,1:3]` observe both categories in every column in both groups (verified).* |
| N11 | **SBM `posterior_num_blocks` prior-convention fix** (zero-truncated → shifted Poisson) | `functions.bgms.R:77` → `extract_sbm()`; `summary.easybgm.R:234` | `res$sbm$posterior_num_blocks` changes |
| N12 | `target_accept` heuristic fix; NUTS acceptance accumulation fix; multinomial candidate weighting; Stage-2 warmup windowing | inherited by every fit | the chain differs from 0.1.6.3 at matched settings; subsumed into N1–N11 for any single number, listed for completeness |

#### COSMETIC

| # | 0.2.0.0 change | easybgm site | effect |
|---|---|---|---|
| C1 | **S7 fit objects** + the shim | the four `class(fit) <- ...` lines S7 cannot survive: `functions.bgms.R:54`, `:61`; `functions.bgmscompare.R:169`, `:242` | one `warning()` per fit; otherwise easybgm receives exactly the 0.1.6.3 object shape. `.subset2` — the shim comment's second justification ([build_output.R:310-311](R/build_output.R#L310-L311)) — appears **nowhere** in easybgm 0.4.0 (**F-025-5**) |
| C2 | Scalar prior arguments deprecated but translated: `pairwise_scale`, `main_alpha`, `main_beta`, `inclusion_probability`, `beta_bernoulli_*`, `dirichlet_alpha`, `lambda`; character `edge_prior` | **exactly** the arguments easybgm documents (`easybgm.R:86-108`) and its tests pass (`edge_prior = "Stochastic-Block"`, `test-easybgm.R:66`) | all still work ([bgm.R:627-670](R/bgm.R#L627-L670)); one `lifecycle` warning each per session. Caveat with teeth: `pairwise_scale = s` → `cauchy_prior(scale = s)`, so a user following easybgm's docs and passing `interaction_scale = 2.5` gets the 0.1.6.3 *family* at the *new* coordinate — a third model, neither old default nor new |
| C3 | `extract_category_thresholds()` deprecated, forwards to `extract_main_effects()` ([extractor_functions.R:1080-1087](R/extractor_functions.R#L1080-L1087)) | `functions.bgms.R:149`, `:196` | one deprecation warning per session (captured verbatim in §1.4 and §3.2). Return shape unchanged for OMRF (`4 x 4` = `p × max_categories`, verified); the value lands in `bgms_res$thresholds`, which **nothing else in easybgm reads** |
| C4 | `extract_arguments()` gains `main_effect_indices` for compare fits ([extractor_functions.R:86-99](R/extractor_functions.R#L86-L99)) | additive | easybgm reads named fields only |

#### NOT REACHED

| # | 0.2.0.0 change | why it does not arrive |
|---|---|---|
| X1 | `update_method = "gibbs"`, `delta`, `precision_scale_prior` eta frame, `precision_graph_prior = "hierarchical"` default, the normalizer-correction tables, the z-ratio gauge and `fit$zratio_diag` | `easybgm()` reroutes `type = "continuous"`/`"mixed"` away from bgms with a warning (`easybgm.R:226-232`); `easybgm_compare()` likewise (`:167-173`) — confirmed live in §3.1 (workflow 7 came back `package_bdgraph`). `zratio_active` is set only on the hierarchical continuous path ([build_output_bgm.R:29](R/build_output_bgm.R#L29)), so an ordinal fit never sees the gauge. **Scoped exception:** a raw GGM/mixed fit handed to a `plot_*.bgms` method does carry these — see H1 and F-025-7 |
| X2 | `update_method = "hamiltonian-mc"` and `hmc_num_leapfrogs` removed | no easybgm call site passes either. A *user* passing them through `easybgm(...)`'s dots gets a hard error — a user-level break, not an easybgm-code break |
| X3 | `standardize` deprecated (`FALSE` warns, `TRUE` errors) | no easybgm call site; not in easybgm's documented prior list |
| X4 | `extract_ess()` RB semantics; `estimator = "mixt"` deprecation; `n_eff_mixt` column removed | easybgm never calls `extract_ess()` and reads `n_eff`, not `n_eff_mixt` |
| X5 | `extract_inclusion_bf()` and its `log=`; natural-log BF display conventions | easybgm hand-computes `inc_BF` (`functions.bgms.R:155-178`, `:203-226`). bgms's log-BF display lives in `summary()`/`print()`/`verdicts()`, none of which easybgm calls |
| X6 | `bgmCompare()` RB inclusion probabilities are `NA` for never-updated main-difference indicators | easybgm's compare path uses `colMeans()` of the **raw** indicator draws restricted to `grep("\\(pairwise\\)")` (`functions.bgmscompare.R:213`, `:291`); the RB matrix and its `NA`s never enter |
| X7 | `category_support`, the `*` marking in `summary.bgmCompare` | easybgm never reads `category_support` and never calls bgms's `summary()`. (The *warning* accompanying the fix does reach the user — that is N10(iii)) |
| X8 | `summary()` gains `quadratic` + `main_label`; and every extractor / method in the "not touched" list of §2.1 | no call site anywhere in the tarball |
| X9 | `coda` dropped from bgms Imports | easybgm imports `coda` itself and calls `coda::effectiveSize` only when `ess = NULL` (`AuxiliaryFunctions.R:304-312`), which the bgms path never triggers |
| X10 | `progress_callback`; C++ backend refactor; Windows/RcppParallel bitwise-reproducibility note; Alpine/musl fix; `predict()` category-map fix | additive or platform-specific |
| X11 | `bgmCompare(difference_probability =)` deprecated with no default | easybgm reads `difference_probability` from **its own** dots (`functions.bgmscompare.R:182-184`, `:247-249`) and never forwards it |
| X12 | `plot_edge_posterior(binwidth =)` | no call site |

#### The one genuinely new reach

| # | what | verdict |
|---|---|---|
| H1 | `bgm()` can now return **GGM** and **mixed-MRF** fits, classes that did not exist in 0.1.6.3. `easybgm()` cannot produce them (X1), but easybgm registers S3 methods **on class `bgms`** — `plot_network.bgms`, `plot_edgeevidence.bgms`, `plot_structure.bgms`, `plot_structure_probabilities.bgms`, `plot_complexity_probabilities.bgms`, `plot_parameterHDI.bgms`, `plot_centrality.bgms`, `plot_centrality.list` — each funnelling a **raw bgms fit** into `bgm_extract.package_bgms` (`plottingfunctions.bgms.R:21-25`, `plottingfunctions.easybgm.R:667-672`, `:778-783`) | **runs clean and produces real plots** (§3.3) — which is the problem. **F-025-7**, and the anchor control proves it is new: on 0.1.6.3 the same two workflows error with *"The bgm function supports variables of type ordinal and blume-capel"* |

---

## 3. Task 3 — the CRAN easybgm run against this bgms

Everything below was run twice, once per bgms build, adjacently, one job at a
time. Nothing was dropped.

### 3.1 Workflows — 14 per build (`02-workflows-{anchor,new}.log`)

| # | workflow | 0.1.6.3 | 0.2.0.0 |
|---|---|---|---|
| 1 | `easybgm(type="binary", save=TRUE, centrality=TRUE)` | OK | OK |
| 2 | `easybgm(type="ordinal", save=FALSE)` | OK | OK |
| 3 | `easybgm(type="blume-capel", baseline_category=2)` | OK | OK |
| 4 | `easybgm(edge_prior="Stochastic-Block")` | OK | OK |
| 5 | `easybgm(edge_prior="Beta-Bernoulli")` | OK | OK |
| 6 | `easybgm(edge_selection=FALSE)` | OK | OK |
| 7 | `easybgm(type="continuous", package="bgms")` | rerouted → `package_bdgraph` | rerouted → `package_bdgraph` |
| 8 | 7 plot methods on an easybgm result | all OK | all OK |
| 9 | raw `bgm()` fit + 6 `plot_*.bgms` methods | all OK | all OK |
| 10 | `easybgm_compare(list(g1,g2))` | OK | OK |
| 11 | `easybgm_compare(group_indicator=)` | OK | OK |
| 12 | compare `summary()` + 3 plots | all OK | all OK |
| 13 | raw GGM fit + 3 plots; raw mixed fit + 2 plots | **ERROR** ×2 (model class does not exist) | **all OK** |
| 14 | — (F-025-1 moved to `04-trap.R`) | — | — |

**Errors: 2 on the anchor, 0 on the new build.** Both anchor errors are
workflow 13, i.e. model classes 0.1.6.3 does not have. Every workflow that
runs on the version CRAN carries also runs on 0.2.0.0.

Conditions raised on 0.2.0.0 and not on 0.1.6.3, verbatim:

```
  * easybgm 0.4.0 is not compatible with S7-based bgms objects. Running in S3
    compatibility mode. Please update easybgm to version 0.5.0 or later.
  * `extract_category_thresholds()` was deprecated in bgms 0.2.0.
    ℹ Please use `extract_main_effects()` instead.
    ℹ The deprecated feature was likely used in the easybgm package.
      Please report the issue at <https://github.com/KarolineHuth/easybgm/issues>.
```

Result shape on 0.2.0.0 (workflow 1) — unchanged field set, and note `MCSE_BF`
carrying `NA` (N8 / F-025-3):

```
RESULT names : edge.prior, parameters, samples_posterior, thresholds, structure,
               inc_probs, inc_BF, structure_probabilities, graph_weights,
               sample_graph, centrality, convergence_parameter, MCSE_BF, model,
               fit_arguments
   parameters               matrix     4 x 4     anyNA=FALSE
   inc_probs                matrix     4 x 4     anyNA=FALSE
   inc_BF                   matrix     4 x 4     anyNA=FALSE
   samples_posterior        matrix     200 x 6   anyNA=FALSE
   centrality               matrix     200 x 4   anyNA=FALSE
   convergence_parameter    numeric    len 6     anyNA=FALSE
   MCSE_BF                  data.frame 6 x 2     anyNA=TRUE
```

### 3.2 easybgm's own test suite (`03-testsuite-{anchor,new}.log`)

| | bgms 0.1.6.3 | bgms 0.2.0.0 |
|---|---|---|
| failed | **0** | **0** |
| errors | **0** | **0** |
| skipped | 0 | 0 |
| passed | 182 | 185 |
| warnings | 0 | 3 |
| suite wall | 4.6 s | 8.6 s (**1.87×**) |

All three test blocks — `easybgm returns expected structure…` (113 expectations),
`plotting functions work…`, `easybgm_compare returns expected structure…` —
pass on both. Which tests reach bgms: all three blocks do; of the 10
`easybgm()` combos, 5 are bgms and 5 are BGGM/BDgraph, and of the 6
`easybgm_compare` combos, 3 are bgms.

The three new warnings are the shim warning and the
`extract_category_thresholds()` deprecation, surfaced because the tests wrap
those blocks in `suppressMessages` rather than `suppressWarnings`
(`test-easybgm.R:145-169`, `:235-245`). They are warnings, not failures, and
do not fail `R CMD check`.

### 3.3 The new-model-class path (F-025-7)

A raw GGM fit handed to easybgm's `plot_*.bgms` methods returns real plot
objects:

```
model_type: ggm  is_continuous: TRUE
plot_network                   -> qgraph
plot_edgeevidence              -> qgraph
plot_structure_probabilities   -> ggplot2::ggplot,ggplot,...

easybgm-extracted GGM result:
  thresholds (from extract_category_thresholds): NULL
  parameters (assoc scale):
            intrusion dreams   flash   upset
  intrusion    0.0000 0.8434  0.5891  0.0153
  dreams       0.8434 0.0000  0.1015  0.4674
  flash        0.5891 0.1015  0.0000  0.1927
  upset        0.0153 0.4674  0.1927  0.0000
  bgms partial correlations (posterior mean):
  intrusion    dreams     flash     upset
     0.5065    0.5093    0.4328    0.4020
```

No error, no warning. But easybgm documents `parameters` as "a p × p matrix
containing **partial associations**" (`easybgm.R:23`), and for a GGM what it
receives is the association-scale coupling (`-0.5 ×` precision off-diagonal),
not the partial correlation a user reading a Gaussian network would expect.
The two differ, and nothing tells the user which one is on the plot. The mixed
fit behaves the same and additionally carries `zratio_diag` and emits the
z-ratio message. This is not a regression against 0.1.6.3 — the path did not
exist — but it is new, unguarded, untested surface.

### 3.4 `R CMD check --as-cran` — the reverse-dependency dress rehearsal

Both runs on `easybgm_0.4.0.tar.gz`, `--as-cran --no-manual`, differing only in
which bgms is on `R_LIBS`.

| | bgms 0.1.6.3 | bgms 0.2.0.0 |
|---|---|---|
| **Status** | **1 WARNING, 2 NOTEs** | **1 WARNING, 2 NOTEs** |
| `checking examples` | OK | OK |
| `checking examples with --run-donttest` | OK **[78 s / 25 s]** | OK **[252 s / 128 s]** |
| `checking tests` | OK (below timing threshold) | OK **[31 s / 17 s]** |
| total wall | **83 s** | **196 s** |

The WARNING and both NOTEs are **byte-identical** across the two runs and none
is caused by bgms:

```
* checking CRAN incoming feasibility ... WARNING
Maintainer: ‘Karoline Huth <k.huth@uva.nl>’
Insufficient package version (submitted: 0.4.0, existing: 0.4.0)
This build time stamp is over a month old.

* checking top-level files ... NOTE
Files ‘README.md’ or ‘NEWS.md’ cannot be checked without ‘pandoc’ being installed.

* checking dependencies in R code ... NOTE
Namespace in Imports field not imported from: ‘igraph’
```

— i.e. re-checking an already-published tarball locally (version not
incremented, stale timestamp), no `pandoc` on this machine, and a pre-existing
unused `igraph` import. A `--no-examples --no-tests` baseline run before the
measurements produced the `igraph` NOTE alone, which is how these were
separated from anything bgms could cause.

**There is no ERROR, no new WARNING, and no new NOTE for CRAN's
reverse-dependency check to flag.** Note also that `--as-cran` runs
`\donttest` as a *second* examples pass on top of the normal one
(`tools:::.check_packages`: `test_donttest <- !run_donttest && as_cran`), so
easybgm's seven bgms-reaching example fits are paid for twice per check.

### 3.5 The runtime result — measured, not inferred

**Per-fit, `cores = 4` pinned on both sides, `warmup` deliberately not pinned**
(`timing-merged.csv`). Every workload is a verbatim easybgm call from the test
suite or the `\donttest` examples.

| workload (source) | 0.1.6.3 | 0.2.0.0 | ratio |
|---|---|---|---|
| T1 `easybgm` binary sv=F cnt=F (`test:27`) | 0.27 | 0.43 | 1.59 |
| T2 `easybgm` binary sv=T cnt=T (`test:28`) | 0.23 | 0.41 | 1.78 |
| T3 `easybgm` binary sv=F cnt=T (`test:29`) | 0.24 | 0.41 | 1.71 |
| T4 `easybgm` blume-capel (`test:30`) | 0.17 | 0.26 | 1.53 |
| T5 `easybgm` binary SBM (`test:31`) | 0.26 | 0.44 | 1.69 |
| T6 `easybgm` binary, plot fixture (`test:135`) | 0.24 | 0.40 | 1.67 |
| T7 `easybgm_compare` 2-group sv=F (`test:222`) | 0.81 | 1.65 | 2.04 |
| T8 `easybgm_compare` 2-group sv=T (`test:223`) | 0.70 | 1.62 | 2.31 |
| T9 `easybgm_compare` multi-group (`test:224`) | 0.85 | 2.04 | 2.40 |
| E1 `easybgm` ordinal Wenchuan[1:50,1:5] (`HDI`/`structure`/`centrality.Rd`) | 0.47 | 0.96 | 2.04 |
| **E2 `easybgm` ordinal FULL Wenchuan** (`complexity_probs`/`structure_probs.Rd`) | **9.08** | **33.41** | **3.68** |
| E3 `easybgm_compare` ADHD 2-group (`easybgm_compare.Rd`) | 0.38 | 0.06 | **0.16** |
| E4 `easybgm_compare` ADHD 4-group (`easybgm_compare.Rd`) | 0.61 | 1.29 | 2.11 |
| **test-suite bgms fits, total** | **3.8 s** | **7.7 s** | **2.03×** |
| **example bgms fits, as written** | **10.5 s** | **35.7 s** | **3.40×** |
| **example fits as the check pays them** (E1 ×3, E2 ×2) | **20.6 s** | **71.0 s** | **3.45×** |

**Per-example, from the checks' own `easybgm-Ex.timings` (elapsed s):**

| example | backend | 0.1.6.3 | 0.2.0.0 | ratio |
|---|---|---|---|---|
| `complexity_probs` | bgms | 9.065 | 58.129 | **6.41** |
| `structure_probs` | bgms | 8.924 | 58.013 | **6.50** |
| `structure` | bgms | 0.485 | 1.894 | 3.91 |
| `centrality` | bgms | 0.557 | 1.896 | 3.40 |
| `HDI` | bgms | 0.613 | 1.984 | 3.24 |
| `easybgm_compare` | bgms | 1.090 | 2.472 | 2.27 |
| `edgeevidence` | **BGGM** | 0.597 | 0.486 | 0.81 |
| `network` | **BGGM** | 0.237 | 0.236 | 1.00 |
| `easybgm` | **BGGM** | 0.026 | 0.026 | 1.00 |
| `prior_sensitivity` | — (body commented out) | 0.004 | 0.004 | — |

The three BGGM-backed examples are unchanged to within 1%. That is the control
that makes every other row attributable to bgms rather than to machine drift.

Three things the numbers say that the default alone did not:

1. **"Roughly 2×" is right only for the test suite** (2.03× on the fits, 1.87×
   on the suite wall). On the examples it is **3.4×**, and on the two heavy
   ones **~6.4× elapsed**. `warmup` doubling is the floor, not the whole
   effect; `target_accept` 0.60 → 0.80 buys more leapfrog steps per iteration
   on top of twice as many iterations.
2. **The elapsed ratio exceeds the CPU ratio on the heavy examples.** For
   `complexity_probs`, CPU (user) goes 32.588 → 115.662 s (**3.55×**, matching
   the pinned-core measurement of 3.68×) while elapsed goes 9.065 → 58.129 s
   (**6.41×**). The parallel speedup over 4 chains therefore degraded from
   ~3.6× to ~2.0× in this configuration. I did not establish the cause and do
   not speculate; it is recorded as an observation worth a look
   (**F-025-6b**). Both check runs were unpinned and used all 15 cores, as
   CRAN's would; CRAN typically limits to 2, so its absolute numbers will
   differ from these while the ratio should not.
3. **The cost is size-dependent, and at the smallest sizes 0.2.0.0 is
   *faster*.** E3 (3 binary variables, 10+10 observations) goes 0.38 → 0.06 s.
   Re-probed in isolation to confirm it was not a swallowed error — it is
   real: 0.194 s vs 0.066 s, `iter/warmup` correctly reported as 50/1000 and
   50/2000, both returning a well-formed 3×3 result. At that size the 0.2.0.0
   backend's per-iteration savings more than pay for twice the warmup; the
   ratio rises monotonically with problem size (0.16 at 3 variables, ~1.7 at
   5, ~2.0–2.4 at 5 with groups, 3.68 at 17). The honest summary is **"up to
   ~3.5× CPU on realistic sizes, and it grows with the model"**, not a flat
   multiplier.

---

## 4. Recommendations

### 4.1 Changes in easybgm — per file, per line

Line numbers are the **CRAN 0.4.0** source. Written so an agent in the easybgm
repo can execute them without reading bgms first.

#### `R/functions.bgms.R`

**E1 — REQUIRED for 0.5.0. Remove the two `class(fit) <- "bgms"` assignments.**
Lines 50-62 currently read:

```r
  if (!inherits(fit, "bgms")) {
    varnames <- fit$var_names
    fit <- fit$packagefit
    class(fit) <- "bgms"                                   # <- line 54
  } else {
    varnames <- fit$arguments$data_columnnames             # <- line 56
    if (is.null(varnames)) {
      varnames <- paste0("V", 1:fit$arguments$no_variables)}
  }
  if(packageVersion("bgms") > "0.1.4.2"){
    class(fit) <- "bgms"                                   # <- line 61
  }
```

Replace with:

```r
  if (!inherits(fit, "bgms")) {
    varnames <- fit$var_names
    fit <- fit$packagefit
  } else {
    fit_args <- bgms::extract_arguments(fit)
    varnames <- fit_args$data_columnnames
    if (is.null(varnames)) {
      varnames <- paste0("V", seq_len(fit_args$no_variables))
    }
  }
```

Why, with the evidence: a 0.2.0.0 fit has `class = c("bgms", "S7_object")` and
`typeof = "object"`. `class(fit) <- "bgms"` strips `"S7_object"` **without
raising anything**, and the next `.subset2()` in `` `$.bgms` `` then fails with
`this object class is not subsettable` — every downstream call, verbatim in
§1.6. `fit$packagefit` already carries the correct class, so the assignment was
never needed. Until this lands, bgms's shim hands easybgm a plain S3 list and
everything works, at one warning per fit.

**E2 — REQUIRED for 0.5.0. Bump `DESCRIPTION: Version` to 0.5.0** in the *same*
release as E1 and E9. The shim keys on `packageVersion("easybgm") < "0.5.0"`;
bumping without E1/E9 converts a warning into a hard failure.

**E3 — RECOMMENDED. Bump the bgms dependency.** `DESCRIPTION:26` has
`bgms (>= 0.1.4)`; easybgm's highest gate anywhere is `> "0.1.6"`, so a
0.2.0.0 fit is indistinguishable from a 0.1.6.3 one inside easybgm. Set
`bgms (>= 0.2.0.0)` once E1/E9 land and drop the dead `< 0.1.4`, `< 0.1.6`,
`< 0.1.6.0` branches (`functions.bgms.R:20-26`, `:47-49`, `:162`;
`functions.bgmscompare.R:7-10`, `:19-34`, `:103-159`).

**E4 — RECOMMENDED. Replace the deprecated extractor.** Lines 149 and 196:
`extract_category_thresholds(fit)` → `bgms::extract_main_effects(fit)`. Return
shape for OMRF is unchanged (verified: `4 × 4` = `p × max_categories`) and
nothing else in easybgm reads `bgms_res$thresholds`, so this is a one-line swap
that removes a user-visible warning. Note it returns `NULL` for a GGM fit —
harmless here, relevant to E14.

**E5 — DECISION NEEDED. The inclusion-probability estimator.** Lines 153, 201.
Either adopt the new RB default (recommended — better estimates, and `inc_BF`
stops going to `Inf`/`0` on short chains) and fix E6 so the block is
RB-consistent, or pin `estimator = "raw"` to reproduce 0.4.0's numbers exactly.
Do not leave it half-and-half.

**E6 — REQUIRED if E5 adopts the new default. `MCSE_BF` mixes two estimators.**
Lines 250-256. Under 0.1.6.3 all three inputs were raw: `gamma_mat` gave
`p_hat`, `BF_vec` came from the raw inclusion probability, `n_eff` was the
transition ESS of the same chain. Under 0.2.0.0 `BF_vec` is RB-derived and
`n_eff` is the RB ESS, while `BF_MCSE` still recomputes `p_hat` from
`gamma_mat` on the raw scale (`AuxiliaryFunctions.R:286`, `:322-326`) and uses
it for the variance. The interval is a valid MCSE for neither. Simplest fix —
give `BF_MCSE()` an optional `p_hat` argument that overrides
`p_raw`/`p_for_variance`, and pass the RB probability:

```r
    rb_p <- extract_posterior_inclusion_probabilities(fit)          # RB
    bgms_res$MCSE_BF <- BF_MCSE(
      gamma_mat = extract_indicators(fit),
      BF_vec    = bgms_res$inc_BF[lower.tri(bgms_res$inc_BF)],
      p_hat     = rb_p[lower.tri(rb_p)],
      ess       = fit$posterior_summary_indicator$n_eff,
      return    = "ci", smooth_bf = FALSE)
```

Alternatively drop easybgm's interval and read
`bgms::extract_inclusion_bf(fit, log = TRUE)` with the `mcse` column of
`fit$posterior_summary_indicator`, which bgms computes on the RB scale
throughout. Either way, document the `NA` rows: bgms reports `NA` where the RB
draws are constant, and `BF_MCSE` already degrades those to `NA` CIs rather
than erroring (confirmed live, §3.1).

**E7 — OPTIONAL but strongly recommended. Stop hand-rolling the inclusion Bayes
factor.** Lines 154-181 and 202-229 reimplement the prior-odds division,
including the SBM within/between case. bgms 0.2.0.0 exports
`extract_inclusion_bf(fit, log = FALSE)`, which divides out the exact prior
inclusion odds on the accumulator's log scale, stays finite down to
log-acceptances of about −745, handles all three edge priors uniformly, and
returns a `p × p` matrix with the same dimnames easybgm expects (verified at
[extractor_functions.R:469-509](R/extractor_functions.R#L469-L509)). Both
blocks collapse to `bgms_res$inc_BF <- bgms::extract_inclusion_bf(fit)`, which
also removes the `stop("Unknown edge prior type.")` dead ends at lines 180 and
228.

**E8 — OPTIONAL. Update the documented prior defaults** (`easybgm.R:86-108`).
`interaction_scale … Cauchy … default is 2.5` is now wrong: the default is
`normal_prior(scale = 1)` on the association scale, and `interaction_scale`
still works but is deprecated and maps to `cauchy_prior()` at the *new*
coordinate. `threshold_alpha`/`threshold_beta` default to 1 is also wrong
(bgms uses `beta_prime_prior(0.5, 0.5)`; pre-existing drift). Replacement
names: `interaction_prior`, `threshold_prior`, `edge_prior` (a prior object).

#### `R/functions.bgmscompare.R`

**E9 — REQUIRED for 0.5.0. Remove `class(fit) <- c("bgmCompare")`** at lines
169 and 242. Same reason as E1; nothing replaces them, since
`fit <- fit$packagefit` (line 89) already yields class
`c("bgmCompare", "S7_object")` and `extract_arguments(fit)` on the next line
dispatches on it.

**E10 — BUG, pre-existing, still wrong. The two-group Beta-Bernoulli branch
reads argument names that have never existed.** Lines 189-194 use
`args$beta_bernoulli_alpha`/`args$beta_bernoulli_beta`. A `bgmCompare` fit's
`$arguments` has never carried those — not in 0.2.0.0
([build_arguments.R:150-190](R/build_arguments.R#L150-L190)) and not in 0.1.6.3
(anchor `output_utils.R:293-321`), where they are `difference_selection_alpha`
and `difference_selection_beta`, which the **multi-group** branch at lines
254-256 already uses correctly. Fix by copying those names. Without it
`edge.prior` is `numeric(0)` and `inc_BF` (line 214) is a zero-length matrix.
Not a 0.2.0.0 regression — flagged because the maintainer will be in this file.

**E11 — OPTIONAL. Adopt the RB estimator on the compare side** (lines 213,
291). Note the RB matrix carries `NA` for main-effect difference indicators the
sampler never updated unless `main_difference_selection = TRUE`, so keep the
existing `(pairwise)` restriction or filter on `is.na()`.

#### `R/plottingfunctions.easybgm.R`

**E12 — BUG, pre-existing. `plot_centrality.list` reads a field that does not
exist.** Line 663: `if(!fit_args$save)`. `save` has never been a field of a
bgms fit's `$arguments` in either version, so this is `if(logical(0))` →
`argument is of length zero`. Every sibling method guards it first (e.g.
`plottingfunctions.bgms.R:10-12`). Add the same guard here and at `:776-783`.

#### `tests/testthat/test-easybgm.R`

**E13 — RECOMMENDED. Pin `warmup`.** Every bgms fit in the suite passes
`iter = 10` and no `warmup`, and bgms's default moved 1e3 → 2e3. Add
`warmup = 10` beside each `iter = itr` for the bgms combos (lines 45-82,
236-260). Measured effect of not doing so: suite wall 4.6 s → 8.6 s. The same
applies with much more force to the `\donttest` examples in `man/`
(`complexity_probs.Rd`, `structure_probs.Rd`: 9.1 s → 58 s each) — adding
`warmup = 100` there would cut the check's dominant cost.

#### New

**E14 — NEW, worth a decision. Guard or support the new bgms model classes.**
A raw GGM or mixed-MRF fit now passes through every `plot_*.bgms` method
without error and produces a plot (§3.3), but `res$parameters` is the
association-scale coupling while easybgm documents that slot as "partial
associations" — for a Gaussian network the natural reading is the partial
correlation, and the two differ. Either (a) detect
`extract_arguments(fit)$model_type %in% c("ggm", "mixed_mrf")` at the top of
`bgm_extract.package_bgms` and `stop()` with a pointer, or (b) support it
properly and use `bgms::extract_partial_correlations()` for the GGM case.
Doing nothing leaves a plausible-looking but differently-scaled network.

### 4.2 Changes in bgms

**B1 — F-025-1, document rather than code.** The gate fires at fit-construction
time, so fit-then-load-easybgm yields an S7 object easybgm cannot read. Because
the failure is **loud** (§1.6), the cheap fix is the right one: one sentence in
the shim's comment block and in the release notes — *load `easybgm` before
fitting if you intend to use easybgm on the result*. Making `` `$.bgms` ``
tolerate the stripped class would work too but entrenches the idiom the shim
exists to retire. Recommendation: **document only**.

**B2 — F-025-5, trim the shim comment.**
[build_output.R:310-311](R/build_output.R#L310-L311) says easybgm "overwrites
`class(fit)` **and uses `.subset2` directly**". `.subset2` appears nowhere in
easybgm 0.4.0. Drop the second clause so a future reader does not hunt for an
incompatibility that is not there.

**B3 — F-025-4, leave but do not treat as contract.** `indicator`,
`interactions`, `thresholds` on `bgms_class`
([class_s7.R:92-95](R/class_s7.R#L92-L95)) serve the *legacy-fit* deprecation
branches, not easybgm — easybgm 0.4.0 reads none of them. Keep them for now;
retire them on the same release that removes the shim.

**B4 — F-025-6b, worth a look before release.** On easybgm's two heaviest
examples the parallel speedup over 4 chains fell from ~3.6× to ~2.0× between
builds while CPU time rose 3.55×. Cause not established here. It is not a
correctness issue and not a CRAN gate, but it is the difference between a 3.5×
and a 6.4× wall-clock regression on the one workload CRAN actually times.

### 4.3 Draft notification memo to the easybgm maintainer

> **Subject: bgms 0.2.0.0 — easybgm 0.4.0 keeps working, but its numbers move**
>
> bgms 0.2.0.0 is going to CRAN. easybgm 0.4.0 does **not** break: your test
> suite passes 0-failure/0-error against it, and `R CMD check --as-cran`
> returns the identical status against old and new bgms. A shim in bgms detects
> easybgm < 0.5.0 and returns the old S3 fit object instead of the new S7 one —
> you will see one warning per fit saying so.
>
> What changes is numbers. Pairwise parameters are now on the association
> scale, so every edge weight easybgm reports — `parameters`,
> `samples_posterior`, `centrality`, the HDI panels — is about **half** its old
> value. The default interaction prior moved from Cauchy(2.5) to Normal(1),
> `warmup` doubled to 2000, and inclusion probabilities are now
> Rao-Blackwellized, so `inc_probs`, `inc_BF` and `structure` shift; `inc_BF`
> in particular stops hitting `Inf`/`0` on short chains. R-hat is now classic
> split-R-hat. On the comparison side a bug that silently merged ordinal
> categories groups did not share is fixed, which moves `easybgm_compare()`
> results substantially whenever groups observed different categories.
>
> Four things worth doing at your convenience:
> **(1)** For an S7-native 0.5.0, delete the four `class(fit) <- ...`
> assignments (`functions.bgms.R:54`, `:61`; `functions.bgmscompare.R:169`,
> `:242`) and bump the version to 0.5.0 — that is all the shim keys on.
> **(2)** `MCSE_BF` now mixes a Rao-Blackwellized `n_eff` with a raw `p_hat`;
> it needs one or the other (`functions.bgms.R:250-256`).
> **(3)** Your bgms fits pass `iter` but not `warmup`. Your `R CMD check` goes
> from 83 s to 196 s, and `complexity_probs.Rd` / `structure_probs.Rd` from
> 9 s to 58 s each. Pinning `warmup` in those two examples and in
> `tests/testthat/test-easybgm.R` recovers most of it.
> **(4)** GGM and mixed bgms fits now reach your `plot_*.bgms` methods and plot
> without complaint, but `parameters` there is not a partial correlation.
>
> Two pre-existing bugs, unrelated to this release, are E10 and E12 in the
> attached report. Happy to send a PR for any of it.

### 4.4 The ordering question, answered

**Does bgms 0.2.0.0 break the CRAN easybgm outright? No.** Evidence, all
executed:

1. easybgm 0.4.0's test suite: **0 failed, 0 errors** against 0.2.0.0
   (185 expectations), same as against 0.1.6.3 (§3.2).
2. **14 workflows × 2 builds**: 0 errors on 0.2.0.0; the only 2 errors anywhere
   are on the *anchor*, for model classes 0.1.6.3 does not have (§3.1).
3. `R CMD check --as-cran`: **identical status** — 1 WARNING, 2 NOTEs,
   byte-identical text, all three artifacts of re-checking a published tarball
   locally, none bgms-caused; verified against a `--no-examples --no-tests`
   baseline (§3.4).
4. The compatibility shim is what buys this, and it was exercised in all three
   gate states with the live object shape recorded each time (§1.4).

**Therefore: submit bgms first, with notification.**

1. **bgms 0.2.0.0 → CRAN now.** The reverse-dependency check has no error to
   flag.
2. **Notify Karoline Huth at submission time**, not after — the memo in §4.3.
   Numbers move for every easybgm user the day bgms lands, and the maintainer
   should be able to answer the first email about it.
3. **easybgm 0.5.0 follows at its own pace.** Nothing forces simultaneity; the
   shim exists precisely so it does not have to be. bgms should not remove the
   shim until easybgm 0.5.0 has been on CRAN long enough that
   `Imports: bgms (>= 0.2.0.0)` is the norm — one bgms minor release later at
   the earliest.

Alternatives rejected: *easybgm first* is impossible, since an S7-native
easybgm 0.5.0 would be broken against the bgms 0.1.6.3 that CRAN carries today;
*simultaneous* couples two release schedules for no gain, and the only thing
that would force it is a hard break, which there is not.

**The one caveat, and it is a budget not a verdict.** easybgm's check time
goes **83 s → 196 s (2.36×)**, driven by the `--run-donttest` pass
(**[78 s/25 s] → [252 s/128 s]**) and within it by two examples that fit full
17-variable Wenchuan (**9.1 s → 58 s elapsed each**). That is well inside
CRAN's limits at this size, so it does not change the order. It does change
the memo: item (3) should reach the maintainer *before* submission, because the
fix is one argument in two `.Rd` files and it is theirs to make.

---

## 5. Findings

| id | severity | where | finding |
|---|---|---|---|
| **F-025-1** | **LOW** | bgms | The shim gate is evaluated at fit-construction time, so `bgm()` → `library(easybgm)` → any easybgm call yields an S7 object easybgm mishandles. Reachable, and the order bgms's own print methods invite. **Fails loudly**: `class(fit) <- "bgms"` succeeds silently but every one of the eight downstream calls errors with `this object class is not subsettable` (§1.6). No wrong number escapes; the message names neither easybgm nor the load order. → B1, document. |
| **F-025-2** | **LOW** | easybgm | easybgm's documented bgms prior arguments are now wrong: `easybgm.R:86-87` states the interaction prior is Cauchy with default 2.5; it is `normal_prior(1)` on a different coordinate. `threshold_alpha`/`threshold_beta` documented as 1 vs bgms's `beta_prime_prior(0.5, 0.5)` is pre-existing drift. → E8. |
| **F-025-3** | **MEDIUM** | easybgm | `MCSE_BF` silently mixes estimators under 0.2.0.0: RB `ess` + RB-derived `BF_vec` + a raw `p_hat` recomputed inside `BF_MCSE` (`functions.bgms.R:250-256`, `AuxiliaryFunctions.R:286`, `:322-326`). Coherent under 0.1.6.3, not now. Produces a plausible interval that is a valid MCSE for neither estimator; `NA` rows appear where RB draws are constant (observed live). → E5/E6. |
| **F-025-4** | **INFO** | bgms | The "easybgm compatibility (deprecated)" properties `indicator`/`interactions`/`thresholds` ([class_s7.R:92-95](R/class_s7.R#L92-L95)) are read by *legacy-fit* deprecation branches, not by easybgm — easybgm 0.4.0 reads none of the three off a bgms fit. → B3. |
| **F-025-5** | **INFO** | bgms | The shim's comment ([build_output.R:310-311](R/build_output.R#L310-L311)) says easybgm "uses `.subset2` directly". It does not; `.subset2` appears nowhere in easybgm 0.4.0. Only the class overwrite is real. → B2. |
| **F-025-6** | **MEDIUM** | both | easybgm's `R CMD check` cost rises **2.36×** (83 s → 196 s); the `--run-donttest` pass **3.2× CPU / 5.1× elapsed**; the two full-Wenchuan examples **~6.4× elapsed** each. Cause: `warmup` 1e3 → 2e3 (easybgm passes `iter` but never `warmup`) plus `target_accept` 0.60 → 0.80. Not a CRAN gate at this size, but it is the reverse-dependency check's only visible change. BGGM-backed examples unchanged to within 1%, which makes it attributable. → E13, memo item (3). |
| **F-025-6b** | **LOW** | bgms | On the two heavy examples the 4-chain parallel speedup fell from ~3.6× to ~2.0× between builds (CPU 3.55× vs elapsed 6.41×). Cause not established; recorded, not diagnosed. → B4. |
| **F-025-7** | **MEDIUM** | both | GGM and mixed-MRF fits — model classes new in 0.2.0.0 — pass through easybgm's `plot_*.bgms` methods without error and produce real `qgraph`/`ggplot` output, but `res$parameters` is the association-scale coupling while easybgm documents that slot as "partial associations". For a GGM the expected quantity is the partial correlation, and they differ (§3.3). New, unguarded, untested surface; the anchor control confirms it did not exist under 0.1.6.3. → E14. |
| **F-025-8** | **LOW** | easybgm, pre-existing | `functions.bgmscompare.R:189-194` reads `args$beta_bernoulli_alpha`/`beta` from a `bgmCompare` fit, which has never carried them in either bgms version; the correct names are `difference_selection_alpha`/`beta`, used correctly by the multi-group branch at `:254-256`. → E10. |
| **F-025-9** | **LOW** | easybgm, pre-existing | `plottingfunctions.easybgm.R:663` does `if(!fit_args$save)` where `save` is not a field of `$arguments` in either version → `argument is of length zero`. Sibling methods guard it first. → E12. |

Nothing in this report requires a change to bgms before submission. B1, B2 and
B3 are documentation and comment hygiene; B4 is worth a look but is not a gate.

---

## 6. Open questions

1. **F-025-6b's cause.** Why did the 4-chain parallel speedup degrade on the
   full-Wenchuan example between builds? Both runs used the same core count and
   the same `chains = 4`. Candidates not tested: lazy summary computation
   moving work out of the parallel region and into single-threaded `$` access;
   oversubscription (`cores = detectCores() = 15` against 4 chains); the
   correction-table cache. One profiled run would settle it.
2. **`cores = parallel::detectCores()` under CRAN's 2-core limit.** Neither
   check NOTEd it here, but bgms's default is what easybgm inherits, and CRAN's
   reverse-dependency machines are core-limited. Whether the measured ratios
   hold at 2 cores is untested — the CPU-time ratio (3.55×) should, the elapsed
   ratio (6.41×) probably will not.
3. **Whether easybgm 0.5.0 should adopt `extract_inclusion_bf()` (E7).** This
   changes easybgm's Bayes factors a second time, on top of N7. Doing both in
   one release is cleaner for users than doing them one release apart, but it
   is the maintainer's call and depends on whether they want to keep easybgm's
   BF definition independent of bgms's.
4. **When the shim comes out.** The report recommends "one bgms minor release
   after easybgm 0.5.0 is on CRAN"; the actual trigger should probably be a
   reverse-dependency check showing no CRAN package still pinned below
   easybgm 0.5.0, which nobody has run yet.
5. **N10's reach in the wild.** The category-union fix does not fire in
   easybgm's own examples (verified), so nothing in this report exercises the
   largest compare-side number-mover end-to-end through easybgm. A real
   two-group dataset with unequal category support would be the honest test,
   and none was available here.

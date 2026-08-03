# Report 15 — vignette + man-page accuracy sweep

Brief 15. Analysis only; nothing in the repo was modified. Produced 2026-08-01
against a clean export of `develop` at `2aa6d6ef` (≥ `5b850410`, as required).

---

## 1. What was done

**Oracle build.** `git archive develop` → `~/bgms-review/val15/bgms-val`;
`R CMD build` + `R CMD INSTALL --library=~/bgms-review/lib-val15`. Nothing was
built in the Dropbox tree (F-009).

* `R CMD build` **fails immediately without pandoc on PATH** — confirms F-045.
  Worked around with RStudio's bundled pandoc 3.8.3
  (`/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64`).
  All five vignettes knit live during `build`; `creating vignettes ... OK`.
* Compilation: 0 warnings at `-Wall` across the TU set (`grep -ci "warning:"
  install.log` → 0), consistent with report 01 f9.

**Vignette renders.** All five rendered **sequentially**, one R process at a
time, against the installed build, with `withCallingHandlers` trapping every
warning and message at the render call and a post-hoc scan of the emitted HTML.
Harness: `~/bgms-review/val15/render_one.R`; logs
`~/bgms-review/val15/vig/render-*.log`.

**Man-page examples.** Examples for the 23 new exports extracted with
`tools::Rd2ex(commentDontrun = TRUE, commentDonttest = FALSE)` and run one
topic per R process with `source(echo = TRUE, print.eval = TRUE)` — printing
deliberately enabled, because `print.eval = FALSE` silently skips every S3
print method and would have hidden D-6 below. `_R_CHECK_LIMIT_CORES_=TRUE` set
throughout, as a budget measure — not as an F-018 re-measurement. Logs
`~/bgms-review/val15/exlogs2/*.log`.

**Machine budget.** Renders sequential; example topics sequential; the only
parallelism is the sampler's own (≤ 2 cores under `_R_CHECK_LIMIT_CORES_`).
Peak footprint ≈ 4 threads. Total wall clock ≈ 12 min of sampling.

**Known rows.** F-005, F-006, F-026, F-031 were read first and are not
re-discovered. F-006's substance was verified independently against
`src/mcmc_diagnostics.cpp` (§5). The two Rd surfaces landing in a parallel
brief — F-049's `?prior_sensitivity_check` extrapolation-bound line and
F-057's `?bgmCompare` group-numbering line — were **skipped entirely** and are
excluded from the coverage claim in §7.

---

## 2. Findings

Severity uses the master list's scale. Cross-references to existing rows are
marked; everything tagged **D-n** is new in this report.

| ID | Sev | Where | Finding |
|---|---|---|---|
| **D-1** | major | `vignettes/prior-sensitivity.Rmd:62-103` | The pasted `ps` report does not match what `print.bgms_prior_sensitivity()` produces, in four places — and never did: the block and the print method were introduced by the *same* commit (`4a7212fd`, #183). Full diff in §4.1. |
| **D-2** | major | `vignettes/diagnostics.Rmd:146,152-153` | The Bayes-factor reading is **inverted**. The chunk prints `BF_10 = 24.59` (evidence *for* inclusion); the prose calls it "small, meaning that there is little evidence for inclusion" and reads `1/BF_10 = 0.041` as "strong evidence for the absence of a network relation between `intrusion` and `physior`". |
| **D-3** | minor | `vignettes/prior-sensitivity.Rmd:143` | Lingering `$\|\Delta \log_{10}\mathrm{BF}\|$` — the one site the F-039 nats conversion (`4cc8ce41`) missed; the two sentences on either side of it were converted. |
| **D-4** | minor | `vignettes/diagnostics.Rmd:58` | Names four transition counts "`n0->0`, `n0->1`, `n1->0`, `n1->1`, reported in the indicator summary". `summary(fit)$indicator` carries only `n0->1` and `n1->0`. |
| **D-5** | minor | `vignettes/prior-sensitivity.Rmd:144-145, 169-170` | Prose names report elements that do not exist: "the report's Details line" (there is no Details line) and 'the report prints these as "not certifiable"' (it prints "are too noisy to assess"). |
| **D-6** | minor | `R/verdicts.R:451-470` | `print.bgms_verdicts()` **errors** on a subset that keeps all five display columns: `subset(v, fragile)` → `Error in log(threshold) : non-numeric argument to mathematical function`. The guard at `:454` was written for exactly this case but tests for a missing *column*, not the dropped `evidence_threshold` *attribute*. |
| **D-7** | minor | `vignettes/diagnostics.Rmd:84-86` | "They are `NA` **only** for an indicator whose Rao-Blackwellized draws are constant to double precision" — falsified by the vignette's own rendered table: `intrusion-flash` has `sd = 2.21e-11` (non-zero) with `mcse`/`n_eff` `NA`. A second NA route exists and the source names it. |
| **D-8** | minor | `R/zratio_surfaces.R:545-546`, `R/zratio_gauge.R:394` | Two shipped runtime messages say 0.00028 nats is "two orders below" 0.003 nats. It is ~10.7×, i.e. one order. The diagnostics vignette (`:237`) says "roughly a tenth" — correct, and it disagrees with the runtime message. |
| **D-9** | note | `vignettes/diagnostics.Rmd:64-66` | "Because the parameter has a well-defined value at every iteration, the full chain — including zeros — is a valid sequence for computing ESS." The reported `n_eff` is explicitly *not* a full-chain ESS; `R/mcmc_summary.R:334-336` rejects that route ("a raw-chain ESS on the effect inflates in inclusion-dominated cells"). `mean`/`sd`/`Rhat` do use the full chain, so the sentence is half-right and reads as a justification for something the package deliberately does not do. |
| **D-10** | note | `src/mcmc_diagnostics.cpp:352-371` | The pooled transition scan carries `prev` across chain boundaries (`int start = (c == 0) ? 1 : 0;`), so `nchains - 1` spurious transitions are counted per parameter. Inert at scale (1 of 3999 for the vignette's 2×2000 fit) but it inflates the very `n0->1`/`n1->0` counts the vignette teaches users to read, and biases `n_eff_mixt` upward for a barely-flipping indicator. |
| **D-11** | note | `vignettes/comparison.Rmd` | The F-002-mandated sentence on the `bgm()` = `normal_prior(1)` vs `bgmCompare()` = `cauchy_prior(1)` asymmetry is **still absent**. Verified against formals: `R/bgm.R:481` vs `R/bgmCompare.R:197`. Cross-ref F-002 (open) — reported here as unlanded, not as new. |
| **D-12** | note | `NEWS.md:93` | Still states the fragility operating point as "0.25 of a threshold on the log10 Bayes factor scale" and that "distances are measured on the log10 Bayes factor scale". The shipped `boundary_distance()` (`R/verdicts.R:70-88`) is natural log, and `?verdicts` + `checking-your-model.Rmd:67` both say 0.58 nats. Numerically the same operating point; the *unit claim* is stale. NEWS is brief 14's lane — flagged here only as a cross-consistency disagreement. |

### Known rows, re-confirmed (no new work needed)

* **F-006** — verified independently in §5. The wrong mixture-ESS gloss is
  indeed gone; the soft survivor at `diagnostics.Rmd:88` is what remains, and
  my independent read of `src/mcmc_diagnostics.cpp:375-387` agrees with the
  lead's corrected reading.
* **F-026** — `intro.Rmd` confirmed pre-0.2: no prior constructors, no
  `verdicts()`, no `precision_graph_prior`, no mixed-MRF mention, no
  bibliography header. It renders and is *not wrong*; it is incomplete.
* **F-031** — 2 of the 23 new exports have no `\examples{}`:
  `extract_inclusion_bf`, `extract_prior_inclusion_probabilities`. All 23 have
  `\value{}`.
* **F-045** — reproduced: `R CMD build` fails without pandoc.

### Clean

* **All five vignettes render with zero errors, zero warnings, zero
  messages-to-stderr.** Runtimes: intro 9.3 s, comparison 5.3 s, diagnostics
  8.3 s, checking-your-model 22.8 s, prior-sensitivity 0.4 s (all chunks
  `eval = FALSE`). The rendered HTML contains no "Warning"/"Error" strings
  other than the word "warnings" in prose (`diagnostics.html:698`).
* **21 of 23 new exports' examples run clean** (0 errors, 0 warnings each); the
  other 2 have no examples. Every documented default equals the actual formal
  (§6).
* The natural-log convention holds everywhere except D-3: `grep` for
  `log10|log-10|base-10|log_{10}` over `R/`, `man/`, `vignettes/` returns
  exactly three hits — the D-3 site, a deliberate historical note in
  `?prior_sensitivity_check` ("former 0.5 carried in log10 units"), and a code
  comment. No stale `1.15/6.91`-style caps survive; `verdicts()` prints
  `presence: log BF > 2.30` (= ln 10) and the sensitivity plot's y-axis reads
  "natural log Bayes factor".

---

## 3. Evidence — claim tables

Every row states the check. "run" = executed against the installed build;
`file:line` = source read. Verbatim outputs are in
`~/bgms-review/val15/vig/out/*.html` and `~/bgms-review/val15/exlogs2/*.log`.

### 3.1 `intro.Rmd`

| Claim | How checked | Verdict |
|---|---|---|
| Supports ordinal / Blume–Capel / continuous (`:22-26`) | `R/bgm.R` `variable_type`; `R/build_spec.R:216` mixed builder | **incomplete** — a mixed (discrete+continuous) MRF is a fourth, unlisted case. F-026 |
| `bgm()` one-sample, `bgmCompare()` groups (`:32-33`) | run | ✅ |
| `Wenchuan` = 362 rows, 17 PTSD items; cols 1–5 are intrusion…physior | run (`head(data)`), `R/datasets.R:10-40` | ✅ |
| `summary(fit)` / `coef(fit)` shapes (`:72,78`) | run | ✅ |
| Median-probability network = threshold PIP at 0.5 (`:83`) | run | ✅ |
| GGM pairwise are unstandardized, "stored as `-0.5 * K_ij`" (`:112`) | `R/extractor_functions.R:1711-1713` (`precision = -2 * associations`) | ✅ exact |
| `extract_partial_correlations()` / `extract_precision()` convert (`:113-115`) | run (example logs) | ✅ |
| `na_action = "impute"` (`:116`) | run — the fit emits "To impute missing values instead, use `na_action = "impute"`" | ✅ |
| GGM chunk `:106-109` | `eval = FALSE`; `continuous_data` is undefined | never executes — no error, but the vignette's only GGM demo is unrunnable as written |

### 3.2 `comparison.Rmd`

| Claim | How checked | Verdict |
|---|---|---|
| "edge weights and category thresholds differ across groups in an **ordinal** MRF" (`:22-23`) | `R/bgmCompare.R:41-42` — `variable_type` is `"ordinal"` or `"blume-capel"` only | ✅ (correctly does *not* claim continuous) |
| `ADHD` subset by `group == 1` / `0` (`:37-40`) | run; `R/datasets.R:63` documents `group` 1 = diagnosed | ✅ |
| `summary(fit)` "shows both baseline effects and group differences" (`:58`) | run — four blocks: thresholds, pairwise, group differences (main), group differences (pairwise) | ✅ |
| `coef(fit)$pairwise_effects_groups[, 1]` is the ADHD network (`:83`) | `R/build_spec.R:447-449`: with `y` supplied, `group = 1` for `x`; `?bgmCompare` documents `x` = Group 1 | ✅ — but the mapping is stated nowhere in the vignette (see F-072) |
| `lower.tri()` fill reproduces the package edge order | run + hand-check: column-major lower-triangle order == row-major upper-triangle order == the printed `pairwise` order | ✅ correct |
| "Difference verdicts are scale-contingent … `difference_scale` … calibration under study" (`:70-73`) | `R/verdicts.R` compare-print caveat (F-060b) | ✅ consistent |
| interaction-prior asymmetry vs `bgm()` | `R/bgm.R:481` `normal_prior(1)` vs `R/bgmCompare.R:197` `cauchy_prior(1)` | **absent** — D-11 / F-002 |

### 3.3 `diagnostics.Rmd`

| Claim | How checked | Verdict |
|---|---|---|
| R-hat is split-R̂ of Vehtari et al. (`:54`) | `.compute_rhat_cpp(split_chains(...))`, `R/mcmc_summary.R:398` | ✅ |
| 1.01 guideline for continuous params (`:54`) | prose guidance | ✅ |
| Indicator R-hat can be `NA`/large; read with transition counts (`:58`) | run — `intrusion-dreams` Rhat `NA` | ✅ on substance |
| "…the transition counts (`n0->0`, `n0->1`, `n1->0`, `n1->1`, reported in the indicator summary)" (`:58`) | run — `colnames(summary(fit)$indicator)` = mean, mcse, sd, n_eff, Rhat, `n0->1`, `n1->0`; source `R/mcmc_summary.R:279-284` keeps only `n01`,`n10` | ❌ **D-4** |
| "the full chain — including zeros — is a valid sequence for computing ESS" (`:64-66`) | `R/mcmc_summary.R:333-336` | ⚠️ **D-9** |
| `n_eff` = composite ESS combining conditional-weight MC error and inclusion MC error (`:68-70`) | `derived_weight_ess()`, `R/mcmc_summary.R:342-354` | ✅ exact |
| `share_incl` = inclusion part's share of composite MC variance (`:71-73`) | `R/mcmc_summary.R:353` | ✅ exact |
| RB columns `NA` **only** when RB draws constant to double precision (`:84-86`) | run — `intrusion-flash` `sd = 2.21e-11`, `mcse`/`n_eff` `NA`; `R/mcmc_summary.R:263-272` names a second route | ❌ **D-7** |
| `n0->1`/`n1->0` are directional flip counts (`:87-88`) | `src/mcmc_diagnostics.cpp:365-370` | ✅ |
| "They record how much the chain explored the two states" (`:88-89`) | §5 | ⚠️ known **F-006(a)** |
| `fit$raw_samples$pairwise` traceplot path (`:101`) | run | ✅ |
| PIP near 1/0/0.5 readings (`:127-129`) | run | ✅ |
| BF from PIP valid when prior inclusion = 0.5, incl. symmetric Beta-Bernoulli (`:133-137`) | `R/extractor_functions.R:385-388` — prior odds are divided out; the two coincide at 1/2 | ✅ |
| `BF_10 = p/(1-p)` for `[1,5]` reproduces `extract_inclusion_bf()` | run: manual 24.59443, `extract_inclusion_bf(fit)[1,5]` = 24.59443, `log = TRUE` = 3.20252 | ✅ identical |
| "the Bayes factor … is small … little evidence for inclusion" (`:146`) | run: `BF_10 = 24.59` | ❌ **D-2** |
| "strong evidence for the **absence** of a network relation between `intrusion` and `physior`" (`:152-153`) | run: `1/BF_10 = 0.0407`; `verdicts(fit)` gives `intrusion-physior  pip 0.961  log_bf 3.20  presence` | ❌ **D-2** |
| `update_method = "nuts"` is the default (`:157`) | `R/bgm.R:490` | ✅ |
| E-BFMI < 0.3 is the concern level (`:167`) | `R/diagnostics_nuts.R:90` | ✅ |
| Divergences "fewer than 0.1% of samples" acceptable (`:173`) | `R/diagnostics_nuts.R:191` `divergence_rate > 0.001` | ✅ exact |
| `nuts_max_depth` default 10 (`:179`) | `R/bgm.R:492` | ✅ |
| Non-reversible steps: projection + round-trip tolerance ∝ step² (`:183`) | `R/bgm.R` docs + NUTS diag block | ✅ |
| `slope_significant` at p < 0.01 (`:201`) | `R/diagnostics_nuts.R:80` — `|t| > 2.58` | ✅ (two-sided; prose singles out negative slopes, harmless) |
| `ebfmi_first_half` "below 0.3" and `var_ratio` "above 2" trip it (`:202-203`) | `R/diagnostics_nuts.R:90` | ✅ exact |
| `var_ratio` = var(first half)/var(second half) (`:203`) | `R/diagnostics_nuts.R:87` | ✅ |
| Gauge runs by default; `options(bgms.zratio_gauge_sweeps = 0L)` disables; default 2 (`:213,247`) | `R/zratio_surfaces.R:641-642`; `?bgm` `R/bgm.R:236-238` | ✅ |
| Summary stored in `fit$zratio_diag`, one row per chain (`:213`) | `R/build_output_bgm.R:355`, `R/class_s7.R:90` | ✅ |
| Only blocks with ≥ 2 variables audited (`:215`) | run — `summarize_zratio_gauge` example: `block_lo 2` | ✅ |
| `flip_rate` flagged above 1% tolerance (`:217`) | `R/zratio_gauge.R:154` `threshold = 0.01` | ✅ |
| `harm_pred` flagged above 0.01 (`:219`) | `R/zratio_gauge.R:155` `harm_threshold = 0.01` | ✅ |
| harm channel only for Bernoulli / Beta-Bernoulli, else `NA` (`:219`) | `R/zratio_gauge.R:67-72` | ✅ |
| All 16 `per_chain` fields (`:223-231`) | run — printed frame carries exactly `flip_rate, flag, se_mean, se_sd, se_mcse, se_se, noise_floor, n_ent, n_ref, n_capped, block_lo, block_hi, amplification, kappa, harm_pred, harm_flag` | ✅ 16/16 |
| `gamma_prior()` shape above 10 → isolated-edge route (`:235-237`) | `R/zratio_surfaces.R:226` `.zratio_surface_shape_hi = 10` | ✅ |
| "at most 0.00028 nats … roughly a tenth of the 0.003 nats" (`:237`) | `R/zratio_surfaces.R:544-546` | ✅ vignette correct; the **runtime message** says "two orders" — **D-8** |
| `n_isolated` in per-chain counters (`:237`) | `R/zratio_gauge.R:376` | ✅ |
| Bound measured only at rates up to `eta = 2` (`:239`) | `R/zratio_surfaces.R:234` `.zratio_mediation_off_eta_hi = 2` | ✅ |
| Surface anchored to components up to 80 variables (`:245`) | `R/zratio_surfaces.R:364` `.zratio_surface_size_cap = 80L` | ✅ |
| Joint spec's graph marginal is the edge prior reweighted by the normalizer (`:249`) | run — the GGM examples emit exactly that notice | ✅ |
| Prior-only mean PIP must equal `a/(a+b)` (`:253`) | run — `sample_graph_prior` example: `bernoulli_prior(0.3)` → 0.2953; `beta_bernoulli_prior(2,4)` → `mean(theta)` 0.3311 vs 1/3 | ✅ |

### 3.4 `checking-your-model.Rmd`

| Claim | How checked | Verdict |
|---|---|---|
| Three verdicts at a threshold (`:52-55`) | run — print header `presence: log BF > 2.30; absence: log BF < -2.30` | ✅ |
| Fragile = a threshold within two standard errors (`:63-65`) | `R/verdicts.R:272-274`, `boundary_distance()` `:83-88` | ✅ |
| "37,010 edge-fits … every verdict error sat within **0.58** of a threshold on the **log** Bayes factor scale" (`:66-68`) | `man/verdicts.Rd:55-57`; 0.25·ln10 = 0.5756 | ✅ ln-based, agrees with `?verdicts` (NEWS disagrees — D-12) |
| Network plot: solid = presence (width ∝ model-averaged weight), dotted grey = undecided, absence not drawn (`:84-86`) | `R/plot_bgms.R:59-80` (`drawn = verdict != "absence"`, `lty = ifelse(presence, 1L, 3L)`) | ✅ exact |
| `plot_edge_posterior(fit, "intrusion", "upset")` (`:90`) | run | ✅ |
| `extract_centrality()` returns the centrality posterior over draws (`:93-95`) | run — `summary(strength)` gives node/mean/lower/upper/p_most_central | ✅ |
| `calibration_check()`: isotonic fit for discrete, PIT for continuous, `kind` column (`:105-141`) | run — `kind` column present, all `pav` on this all-ordinal fit | ✅ |
| Band is built by resampling each case's **category**, not threshold events (`:125-130`) | `R/calibration_check.R` band construction | ✅ |
| `simulate()` with `method = "posterior-sample"`, one replicate per draw (`:162-171`) | run | ✅ |
| Sum-score scale "has 46 points" (`:213`) | run — `scores = 0:(9 * 5)` → 46 | ✅ |
| "about a tenth of the scale sits outside a 95% pointwise band" (`:211`) | run — printed `0.1086957` | ✅ |
| "bgms ships no joint-level check, and the omission is deliberate" (`:150`) | NAMESPACE — no such export | ✅ |

### 3.5 `prior-sensitivity.Rmd`

| Claim | How checked | Verdict |
|---|---|---|
| "A pairwise interaction in `bgm()` is half the log odds ratio between adjacent response categories" (`:30-31`) | `R/extractor_functions.R:1832` — "log adjacent-category odds ratio = 2 * association" | ✅ exact |
| `normal_prior(scale = 1)` puts ~2/3 of prior mass on OR between ~1/7 and 7 (`:32-34`) | arithmetic: ±1 on the association → ±2 on log OR → OR ∈ [0.135, 7.39] | ✅ |
| Works on any `bgm()` fit made with edge selection (`:53-54`) | `R/prior_sensitivity.R` guard | ✅ |
| Anchors default `c(0.4, 0.63, 1, 1.6, 2.5)` (`:117-118`) | `R/prior_sensitivity.R:205` | ✅ |
| `ess_floor` default 400, points below are masked `NA` (`:132-133`) | `R/prior_sensitivity.R` formals; run — "Points with reweighting effective sample size below 400 are not shown." | ✅ |
| 1x anchor is the original fit, never refit (`:123-126`) | run — "the 1x anchor is the original fit" | ✅ |
| Mover categories stored as `stable`, `indistinguishable-from-wobble`, `moved-beyond-wobble` (`:152-154`) | `R/prior_sensitivity.R:993-995` | ✅ |
| Plot: title states the answer; anchor dots; shaded undecided band; verdict zones labeled; movers colored + named; low-ESS gaps (`:160-164`) | `R/prior_sensitivity.R:1286-1330` — `main = title`, `axis(1, at = anchors)`, `rect(...)`, zone text "presence"/"undecided"/"absence", `mover_palette()` | ✅ (F-067 will remove the title later; accurate *today*) |
| Table carries chosen-scale verdict + **natural log** BF + MCSE, per-anchor verdicts, stability range, mover, `insufficient` (`:166-171`) | `R/prior_sensitivity.R` `$edges` build | ✅ |
| "the report prints these as **not certifiable**" (`:170`) | run — the report says "N edges are too noisy to assess" | ❌ **D-5** |
| Noise band = spread of $\|\Delta\log_{10}\mathrm{BF}\|$ (`:143`) | `4cc8ce41` converted the neighbours, missed this | ❌ **D-3** |
| "…in the report's **Details line**" (`:144-145`) | run — the footer is `Method:` / `Refits:` / `Noise:` | ❌ **D-5** |
| Refit cost: chosen-scale free, others warm-started, short warmup; AM/Gibbs refit with the same sampler (`:188-200`) | run — "5 nuts refits, warm-started from the original fit" | ✅ |
| Warm starts weaken split-R̂; gate also leans on ESS + per-chain verdict agreement; failing anchor dropped (`:202-206`) | run — the example emitted "The repeated 1.6x-scale refit did not converge … the mover rule falls back to the tolerance and Monte Carlo error floors" | ✅ (and the code path fires) |
| The pasted `ps` output block (`:62-103`) | run | ❌ **D-1**, §4.1 |

---

## 4. Reproductions

### 4.1 D-1 — the sensitivity report block

Reproduction: run the `?prior_sensitivity_check` example
(`~/bgms-review/val15/exlogs2/prior_sensitivity_check.log`), or

```r
fit = bgm(Wenchuan[, 1:6], chains = 2)
prior_sensitivity_check(fit)
```

Real output (verbatim, trimmed to the divergent parts):

```
  robust (same verdict at every scale)     14
  changed, within run-to-run noise          0
  changed, beyond run-to-run noise          1
...
Method:  43-point curve from 5 anchor fits (0.4x to 2.5x the chosen scale),
         joined by importance reweighting; the 1x anchor is the original fit.
         Points with reweighting effective sample size below 400 are not shown.
Refits:  5 nuts refits, warm-started from the original fit, 20 s total.
Noise:   no threshold-relevant edge had a measurable spread between two identical
         refits at 1.6x, so there is no run-to-run yardstick; verdict moves are
         judged against the tolerance and Monte Carlo error alone.
See ?prior_sensitivity_check for the full construction.
```

| Vignette (`:62-103`) | Shipped `print.bgms_prior_sensitivity` | Source |
|---|---|---|
| `not certifiable (chains disagree)` (`:73`) | `not certifiable (too noisy to assess)` | `R/prior_sensitivity.R:1000` |
| `12 edges cannot be certified from this run: …` (`:83-84`) | `12 edges are too noisy to assess: …` | `:1039-1044` |
| — (absent) | a follow-on sentence naming the cause and the remedy ("Their Bayes factor sits within Monte Carlo error of an evidence threshold, or their chains disagree … Run more iterations to settle these verdicts…") | `:1046-1068` |
| single wrapped `Details: …` paragraph ending `?prior_sensitivity_check for how to read this.` (`:98-102`) | three labeled blocks `Method:` / `Refits:` / `Noise:` + `See ?prior_sensitivity_check for the full construction.` | `:1143-1169` |

Row 4 (the footer) is **observed** in the run above. Rows 1–3 all belong to the
uncertifiable-edge block, which the small demo fit did not trigger (0
uncertifiable edges, so the fourth count label is suppressed too); they are read
from source at the cited lines rather than seen in output. The mover table
(`:76-81`) and the per-scale verdict counts (`:86-93`) *were* observed and their
formats match the vignette.

Provenance: `git show 4a7212fd:vignettes/prior-sensitivity.Rmd` carries the same
block, and `git log -S 'Method:  %d-point curve'` and
`git log -S 'too noisy to assess'` both return only `4a7212fd`. The block is a
hand-written mock that has never corresponded to any shipped print method. The
`0.31 log10 → 0.71 log BF` edit in `4cc8ce41` touched a number inside a
paragraph that does not exist in the output.

### 4.2 D-2 — inverted Bayes-factor reading

```r
fit = bgm(Wenchuan[, 1:5], seed = 1234, chains = 2,
          display_progress = "none", verbose = FALSE)
p = coef(fit)$indicator[1, 5]          # intrusion-physior
p / (1 - p)                            # 24.59443
1 / (p / (1 - p))                      # 0.04065962
extract_inclusion_bf(fit)[1, 5]        # 24.59443   (identical)
verdicts(fit)                          # intrusion-physior  pip 0.961  log_bf 3.203  presence
```

The vignette reads 24.59 as "little evidence for inclusion" and 0.041 as
"strong evidence for the absence". Both readings are the reciprocal of the
truth, and `verdicts()` calls the same edge `presence`. Most likely a casualty
of the 0.2.0.0 slab-default change (Cauchy → `normal_prior(1)`, F-002/F-029):
the narrower slab strengthens inclusion evidence on this edge, and the prose
was not re-read. **This is the vignette's only worked Bayes-factor example.**

### 4.3 D-6 — `print.bgms_verdicts` errors on a column subset

```r
fit = bgm(Wenchuan[, 1:5], seed = 1234, chains = 2,
          display_progress = "none", verbose = FALSE)
v = verdicts(fit)
subset(v, fragile)
#> Error in log(threshold) : non-numeric argument to mathematical function
```

Behaviour map (`~/bgms-review/val15/verd_bug.R`):

| expression | result |
|---|---|
| `v` | OK |
| `v[v$fragile, ]` | OK — row-only subsetting keeps the `evidence_threshold` attribute |
| `v[1:3, ]`, `head(v)`, `v[order(v$log_bf), ]` | OK |
| `v[, c("parameter","log_bf","verdict")]` | OK — falls into the `:454` plain-data-frame guard |
| `subset(v, fragile)` | **ERROR** |
| `v[i, c("parameter","pip","log_bf","verdict","fragile")]` | **ERROR** |

`[.data.frame` drops non-standard attributes as soon as `j` is supplied, so
`attr(x, "evidence_threshold")` is `NULL` while all five required columns
survive — the one combination the guard does not cover. `?verdicts`'s own
example (`v[v$fragile, ]`) happens to take the safe path, which is why this
survived. Fix shape: test the attribute, not just the columns.

### 4.4 D-4 / D-7 — indicator summary columns

```r
summary(fit)$indicator
#>                         mean        mcse           sd    n_eff      Rhat n0->1 n1->0
#> intrusion-dreams  1.00000000          NA 0.000000e+00       NA        NA     0     0
#> intrusion-flash   1.00000000          NA 2.206767e-11       NA 0.9999966     0     0
#> ...
```

Seven columns; `n0->0` and `n1->1` are absent (D-4). `intrusion-flash` has
non-zero `sd` yet `NA` `mcse`/`n_eff` and a *finite* `Rhat` — a different NA
route from `intrusion-dreams`'s exactly-constant chain, so the vignette's
"only … constant to double precision" is too narrow (D-7). The source names
the second route at `R/mcmc_summary.R:263-268` ("Draws whose variance falls
below the autocovariance kernel's numerical floor pick up an NA ESS there
rather than here").

### 4.5 D-8 — order-of-magnitude claim

`R/zratio_surfaces.R:544-546` and `R/zratio_gauge.R:392-395` both print
"0.00028 nats … two orders below the 0.003 nats". 0.003 / 0.00028 = 10.7 — one
order. `diagnostics.Rmd:237` says "roughly a tenth", which is right. Fix the
message, not the vignette.

---

## 5. F-006 — independent read of the mixture ESS

`src/mcmc_diagnostics.cpp:375-387` (inside `IndicatorESSWorker::operator()`):

```cpp
double p_hat = sum_x / n_total;
double sd = std::sqrt(p_hat * (1.0 - p_hat));
...
double a = (double)c01 / (c00 + c01);   // empirical P(0 -> 1)
double b = (double)c10 / (c10 + c11);   // empirical P(1 -> 0)
double tau_int = (2.0 - a - b) / (a + b);
n_eff_mixt = n_total / tau_int;
mcse = sd / std::sqrt(n_eff_mixt);
```

**What the code computes.** It fits a first-order two-state Markov chain to the
pooled binary indicator draws, reads off the two directional switching
probabilities `a` and `b`, forms that chain's integrated autocorrelation time
`τ = (2 − a − b)/(a + b)`, and returns `n_eff_mixt = n_total / τ` with the
matching Bernoulli MCSE. It is a **mixing-rate ESS for the binary indicator
chain**, not a precision statement about the inclusion probability. Its
degenerate case is explicit: `c01 + c10 == 0` ⇒ `n_eff_mixt = NA`,
`mcse = NA` — a never-flipping indicator gets no transition ESS at all. The
four raw counts `n00, n01, n10, n11` are returned alongside (columns 4–7).

**What the vignette says.** The rewritten `diagnostics.Rmd` no longer describes
`n_eff_mixt` at all — it is not in the shipped `summary(fit)$indicator`, and the
vignette correctly attributes `mean/mcse/sd/n_eff/Rhat` to the
**Rao-Blackwellized** inclusion probability (`:81-86`), keeping `n0->1`/`n1->0`
as separate exploration counters (`:87-92`). **That matches the code.** My read
agrees with the lead's corrected reading; F-006's substance is fixed.

**The two residues I can confirm.**

1. `diagnostics.Rmd:88-89` — "They record how much the chain explored the two
   states". A raw count of transitions is a *count*, not a measure of
   exploration; the code's own exploration measure is `n_eff_mixt`, which is
   not shown. This is F-006(a), unchanged.
2. New, code-side (**D-10**): the pooling loop carries `prev` across chain
   boundaries. At `c = 1, i = 0`, `prev` is still the last draw of chain 0, so
   a transition is counted between two independent chains. That is
   `nchains − 1` fabricated transitions per parameter. For a decisive edge with
   `n0->1 = n1->0 = 0` it can turn an exactly-zero count into 1 and flip
   `n_eff_mixt` from `NA` to a finite (and meaningless) value — the precise
   corner the vignette tells readers to interpret. Negligible in magnitude,
   cheap to fix (`int start = 1;` with `prev` reset per chain), and worth a
   backlog line rather than a release gate.

I did not re-check the docs-site copy or `doc/diagnostics.Rmd` — F-006(b)
covers the stale pre-build artifact.

---

## 6. Man-page spot-checks — the 23 new exports

`E`/`W` = errors / warnings raised while running the example.
Defaults verified formal-by-formal against the `\usage{}` block.

| Export | Example runs? | E/W | Claims verdict |
|---|---|---|---|
| `cauchy_prior` | yes | 0/0 | ✅ `scale = 1` (`R/priors.R:49`); prints `Cauchy(0, 1)` |
| `normal_prior` | yes | 0/0 | ✅ `scale = 1` (`:93`) |
| `beta_prime_prior` | yes | 0/0 | ✅ `alpha = 0.5, beta = 0.5` (`:138`) |
| `gamma_prior` | yes | 0/0 | ✅ `shape = 1, rate = NULL, eta = NULL` (`:208`); rate-XOR-eta honoured — `gamma_prior()` prints "eta = 1, standardized frame", `gamma_prior(shape = 2, rate = 0.5)` prints the raw frame |
| `exponential_prior` | yes | 0/0 | ✅ `rate = NULL, eta = NULL` (`:297`); `exponential_prior()` → `Exponential(eta = 1, standardized frame)` |
| `bernoulli_prior` | yes | 0/0 | ✅ `inclusion_probability = 0.5` (`:338`) |
| `beta_bernoulli_prior` | yes | 0/0 | ✅ `alpha = 1, beta = 1` (`:375`) |
| `sbm_prior` | yes | 0/0 | ✅ all six defaults = 1 (`:436-438`); print block matches |
| `prior_sensitivity_check` | yes (32.5 s) | 0/0 | ✅ all 15 formals match `\usage{}`; example exercises the non-converged-replicate branch and reports it honestly. **F-049's Rd bound line skipped per brief.** |
| `verdicts` | yes (18.4 s) | 0/0 | ✅ `evidence_threshold = 10`; the `evidence_threshold = 30` variant re-reads the same fit correctly (boundaries move to ±3.40 = ln 30). Details' 37,010 / 0.58 / 74% / 94% / 3.0% figures agree with `checking-your-model.Rmd`. Caveat: **D-6** on subsetting. |
| `calibration_check` | yes (19.9 s) | 0/0 | ✅ `nrep = 200, probs = c(0.025, 0.975), grid_size = 101, seed = NULL, ndraws = 500` (`R/calibration_check.R:463-469`); the F-040 units footer prints ("0.3 is 30% of the curve, not 0.3%") |
| `extract_centrality` | yes (18.3 s) | 0/0 | ✅ `measure = "strength", group = 1` (`R/centrality.R:74`); returns node/mean/lower/upper/p_most_central as documented |
| `plot_edge_posterior` | yes (18.6 s) | 0/0 | ✅ `evidence_threshold = 10`, `binwidth = lifecycle::deprecated()` — Rd `\usage{}` shows the deprecated form, i.e. **the Rd is current with brief 10**, and report 00c's `binwidth = 0.01` is the stale entry |
| `sample_graph_prior` | yes | 0/0 | ✅ every commented expectation holds: `bernoulli_prior(0.3)` → 0.2953; `beta_bernoulli_prior(2,4)` `mean(theta)` → 0.3311 vs 1/3; joint spec → 0.2153, "shifted away from 0.3 by the Z(Gamma) tilt" |
| `sample_ggm_prior` | yes | 0/0 | ✅ `dim(K_offdiag)` = 200 × 6 as the comment says; the withheld-edge check returns `TRUE` |
| `sample_sbm_prior` | yes | 0/0 | ✅ `edge_prior = sbm_prior(), seed = 1L` |
| `summarize_zratio_gauge` | yes (2.7 s) | 0/0 | ✅ `threshold = 0.01, verbose = TRUE, harm_inputs = NULL, harm_threshold = 0.01` (`R/zratio_gauge.R:154-155`); the printed `per_chain` frame contains all 16 documented fields |
| `extract_inclusion_bf` | **no example** | — | Rd claims verified by computation instead: the RB identity reproduces `p/(1-p)` exactly on a prior-0.5 fit (24.59443 both ways), and `log = TRUE` returns natural log (3.20252 = ln 24.59). F-031 |
| `extract_main_effects` | yes (6.9 s) | 0/0 | ✅ no formals beyond `bgms_object`; returns the threshold matrix |
| `extract_precision` | yes | 0/0 | ✅ returns the full precision matrix for a GGM; consistent with `intro.Rmd`'s `-0.5 * K_ij` statement |
| `extract_partial_correlations` | yes | 0/0 | ✅ `rho_ij = -Theta_ij / sqrt(Theta_ii Theta_jj)`; diagonal 1, values in (−1, 1) |
| `extract_log_odds` | yes (7.3 s) | 0/0 | ✅ returns `2 * associations` (`R/extractor_functions.R:1832`), matching the vignette's half-log-OR framing |
| `extract_prior_inclusion_probabilities` | **no example** | — | Formals `iter = 4000L, warmup = 1000L, recompute = FALSE` match `\usage{}`. F-031 |

Notes from the run that are not findings but belong on the record:

* `extract_centrality`, `calibration_check` and `prior_sensitivity_check`
  examples emit a genuine `NUTS issues: Warmup may be incomplete …` notice on
  the short default fits. Correct behaviour, but a CRAN reader sees a
  convergence complaint in three of the package's flagship examples.
* Every example was run with `_R_CHECK_LIMIT_CORES_=TRUE` and none failed or
  stalled under it. I did **not** instrument concurrency, so this is not an
  independent re-measurement of F-018 — report 05's measured 3.93 → 1.98 stands
  as that evidence.
* `verdicts()`'s print emits `(1 indicators)` and
  `Bayes factor of 30 (and 0.0333333 for absence)` — grammar and `%g`
  formatting, note-level cosmetics only.

---

## 7. Cross-consistency (vignette ↔ man page ↔ runtime)

| Subject | Vignette | Man page | Runtime | Verdict |
|---|---|---|---|---|
| Verdict categories | presence / undecided / absence (`checking-your-model.Rmd:52-55`) | `?verdicts` `\value` factor levels | `presence 11 \| undecided 15 \| absence 10` | ✅ agree |
| Evidence threshold unit | "log Bayes factor scale", 0.58 (`:67`) | `man/verdicts.Rd:55-57`, same | `presence: log BF > 2.30` (= ln 10) | ✅ agree — but **NEWS.md:93** says log10 / 0.25 → **D-12** |
| Sensitivity noise-band unit | `\log_{10}` at `:143`, "natural log" at `:147` | `?prior_sensitivity_check` (nats) | `Noise: … log BF` (nats) | ❌ vignette internally inconsistent → **D-3** |
| "not certifiable" wording | `:73`, `:170` | — | "too noisy to assess" | ❌ **D-1 / D-5** |
| Isolated-route error bound | "roughly a tenth" (`:237`) | `?summarize_zratio_gauge` (states the numbers, no ratio) | message says "two orders below" | ❌ **D-8** (runtime is the wrong one) |
| Gauge per-chain fields | 16 named (`:223-231`) | `?summarize_zratio_gauge` `\value` | printed frame | ✅ all three agree, 16/16 |
| Fragility caveat on compare fits | comparison vignette silent | `?verdicts` Details: operating point is single-network only | compare print carries the disclaimer (F-060) | ✅ no overclaim; the vignette simply does not raise it |
| Group numbering | `pairwise_effects_groups[, 1]` = ADHD, unstated (`comparison.Rmd:83`) | `?bgmCompare`: `x` = Group 1, `y` = Group 2 | `group1`/`group2` colnames | ✅ consistent for the x/y path. The `group_indicator` path is the opaque one — F-072, and **F-057's Rd line was skipped per brief** |
| `bgm()` vs `bgmCompare()` interaction prior | not mentioned | `?bgm` / `?bgmCompare` each state their own default; neither flags the asymmetry | `normal_prior(1)` vs `cauchy_prior(1)` | ❌ **D-11 / F-002** — a user reading either page alone cannot see it |
| Indicator transition counts | four named (`diagnostics.Rmd:58`) | — | two present | ❌ **D-4** |

---

## 8. Coverage statement

**Vignettes: 5 of 5 covered** (`intro`, `comparison`, `diagnostics`,
`checking-your-model`, `prior-sensitivity`), rendered and claim-swept.

**Man pages: 23 of 69 covered** — the 23 new exports from report 00c §1, as
briefed. Examples were run for the 21 that have them; the 2 without
(`extract_inclusion_bf`, `extract_prior_inclusion_probabilities`) had their
Details/Value claims checked by computation and source read instead.

**Deliberately skipped (2 line-level surfaces, per brief):**
`?prior_sensitivity_check`'s F-049 extrapolation-bound line and
`?bgmCompare`'s F-057 first-appearance group-numbering line. Their *containing*
Rd files were otherwise checked.

**The 46 man pages NOT checked**, listed in full:

`ADHD`, `Boredom`, `Wenchuan`, `bgm`, `bgmCompare`, `bgms-package`,
`cash-.bgmCompare`, `cash-.bgms`, `coef.bgmCompare`, `coef.bgms`,
`extract_arguments`, `extract_category_thresholds`, `extract_edge_indicators`,
`extract_ess`, `extract_group_params`, `extract_indicator_priors`,
`extract_indicators`, `extract_pairwise_interactions`,
`extract_pairwise_thresholds`, `extract_posterior_inclusion_probabilities`,
`extract_rhat`, `extract_sbm`, `mrfSampler`, `plot.bgmCompare`, `plot.bgms`,
`plot.bgms_calibration`, `plot.bgms_centrality`, `plot.bgms_prior_sensitivity`,
`predict.bgmCompare`, `predict.bgms`, `print.bgmCompare`, `print.bgms`,
`print.bgms_calibration`, `print.bgms_prior_sensitivity`,
`print.bgms_verdicts`, `simulate.bgmCompare`, `simulate.bgms`, `simulate_mrf`,
`summary.bgmCompare`, `summary.bgms`, `summary.bgms_centrality`,
`unpack_indicator_prior`, `unpack_interaction_prior`,
`unpack_parameter_prior`, `unpack_scale_prior`, `unpack_threshold_prior`.

Partial exception, stated so the claim is honest: `?bgm`, `?bgmCompare`,
`?verdicts` and `?summarize_zratio_gauge` were read *for the specific claims in
the tables above* (defaults, `target_accept`, `zratio_diag`, group numbering,
the 37,010 study) without a full page sweep. `print.bgms_verdicts.Rd` was not
swept, but its implementation carries D-6.

Also not covered by this brief: NEWS.md (brief 14), the docs site, and
`doc/diagnostics.Rmd` (F-006(b)).

---

## 9. Open questions

1. **D-1 rewrite policy.** The pasted block is a mock, and a mock will go stale
   again. Options: (a) make the chunk live with a small fit (`Wenchuan[, 1:6]`,
   `chains = 2` runs the whole check in ~33 s — measured, and inside the
   existing show-uncapped/run-capped pattern F-045 asks us to preserve); or
   (b) regenerate the block from a real run and add a note that it is pasted.
   I lean (a): it is the only way the vignette stops drifting, and the cost is
   already inside the vignettes' 78 s budget. MM's call.
2. **D-2 scope.** The inverted reading is almost certainly a default-change
   casualty. Worth grepping the tutorial and docs-site copies for the same
   paragraph before the fix batch — it may have been copied.
3. **D-6 fix shape.** Guard on `is.null(threshold)` in addition to the missing
   column, or re-attach the attributes in a `[.bgms_verdicts` method? The
   second is more work but also fixes `subset()` keeping a table that prints as
   a verdicts table. Lead's call; I did not implement either.
4. **D-10.** Backlog or 0.2.0.0? It is a two-line C++ change with a visible
   (if tiny) effect on a user-facing column. My read: backlog, unless the
   0.2.0.0 batch is already touching that file.
5. **F-002's comparison-vignette sentence (D-11)** is still unwritten. It is
   one sentence and this brief is analysis-only; it should ride along with
   whichever batch fixes D-1/D-2 in the same file set.
</content>
</invoke>

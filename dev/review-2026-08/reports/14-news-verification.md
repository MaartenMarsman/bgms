# Report 14 — NEWS.md verification and reconstruction for 0.2.0.0

Brief 14. Branch `fix/news-0.2.0.0`, worktree `~/bgms-review/wt-fix6`, based
on `develop` at `2aa6d6ef`. One commit, `80a28dfb`, touching `NEWS.md` only
(+915 / −145). Not pushed. No code changed, so no test suite and no
`R CMD check` was run.

Mandate: F-001 (blocker), F-023 (minor), F-029 (major).

---

## What was done

1. Inventoried the release: 290 commits `cran-0.1.6.3..HEAD`, 118 of them
   squash-merged PRs; the 23 new exports from report 00c; the F-029 list of
   silent behaviour changes.
2. Verified every pre-existing `# bgms 0.2.0.0` NEWS entry (98 bullets across
   five sections) against the code at HEAD. Claims table below.
3. Rewrote the section as one coherent `# bgms 0.2.0.0` addressed to a
   0.1.6.3 user, ordered Breaking changes → New features → Other changes →
   Bug fixes → Deprecated. Sections `# bgms 0.1.6.3` and older are untouched
   byte-for-byte.
4. Ran the parse gate.

Severity legend for findings: **blocker** = would ship a false statement in a
release artefact; **major** = a silent behaviour change with no user-facing
record; **minor** = misleading but not false.

---

## Findings

### F-001 — already fixed on `develop`; entry further sharpened. Status: resolved.

The text F-001 quotes ("read from the anchor with the highest importance
ESS") is **not present at HEAD**. It was corrected on 2026-07-30 by
`45a4fae7` ("Release meta for 0.2.0.0", PR #187), whose first commit message
is *"docs: describe the sensitivity curve as pooled over anchors in NEWS"*.
`git log -S "highest importance ESS" -- NEWS.md` returns nothing at HEAD;
`git show 45a4fae7 -- NEWS.md` shows exactly the one-line swap.

I re-verified the surviving text against the code rather than taking the fix
on trust:

| NEWS claim at HEAD | Code | Verdict |
|---|---|---|
| "pools every anchor that clears `ess_floor` (default 400) there" | `R/anchor_curve.R:233` `contrib = which(ess_mat[p, ] >= ess_floor)`; default at `R/prior_sensitivity.R` signature (`ess_floor = 400`, report 00c §2) | accurate |
| "weighting each anchor's inclusion-probability estimate by its inverse variance" | `R/anchor_curve.R:244` `wa = ess_mat[p, a] / pmax(pa * (1 - pa), 1e-6)` — ESS/(p(1−p)), the reciprocal of p(1−p)/ESS | accurate |
| "transforming the pooled probability to the Bayes-factor scale" | `R/anchor_curve.R:251` pools on the PIP scale; transform happens downstream | accurate |
| "points where no anchor clears the floor are reported `NA`" | `R/anchor_curve.R:234` `if(length(contrib) == 0L) next` over an `NA_real_`-initialised matrix | accurate |
| man page agreement | `man/prior_sensitivity_check.Rd:135-144` states the same pooling and the same weight | accurate |

One thing the entry did **not** say, which `assemble_curve()`'s own header
comment (`R/anchor_curve.R:207-209`) and the man page
(`man/prior_sensitivity_check.Rd:143-144`) both make a point of: the pooled
curve is pooled *everywhere*, including at an anchor's own grid point.
Exactness lives on the anchor fits' Rao-Blackwellized statistics (the
per-anchor verdict columns and the chosen-scale quantities,
`R/prior_sensitivity.R:706-709`), not on a pooled row. The rewritten entry
now says so, which is the "anchors exact, between-anchor curve reweighted"
distinction the brief asks for, stated precisely.

**F-001 can be closed.** It was fixed before this brief; nothing in it
survived to HEAD.

### F-100 (new, minor) — a stale NEWS claim: `pairwise_scale` has no default at 0.2.0.0

The pre-existing Breaking-changes bullet read:

> Default `pairwise_scale` changed from 2.5 to 1 to match the
> association-scale reparameterization.

This is false as written. `pairwise_scale` is a **bare deprecated formal**
with no default (`R/bgm.R:501`); passing it warns and is translated to
`cauchy_prior(scale = pairwise_scale)` (`R/bgm.R:561-569`) — i.e. a
deprecated call keeps the *Cauchy* family, so it does not land on the new
default at all. What changed is the default of its replacement:
`interaction_prior = normal_prior(scale = 1)` (`R/bgm.R:481`) where
0.1.6.3 had `pairwise_scale = 2.5` on a Cauchy (`cran-0.1.6.3:R/bgm.R:389`).
A user following the old bullet would conclude that
`bgm(pairwise_scale = 1)` reproduces the 0.2.0.0 default; it does not.

Replaced by an accurate default-change entry, and the deprecation path is
described in its own bullet.

### F-101 (new, minor) — a superseded NEWS claim: the old `plot_edge_posterior()` description

The New-features section still carried the pre-redesign description:

> `plot_edge_posterior()` draws one edge's spike-and-slab posterior: a stem
> at zero whose height is the posterior probability that the edge is absent
> … drawn as probability per weight bin.

At HEAD the panel draws a density against the closed-form prior with the
mass at zero on a probability wheel (`R/plot_bgms.R:534+`, roxygen at
`:529-560`), and `binwidth` is `lifecycle::deprecated()`
(`R/plot_bgms.R:628-632`). The old bullet and the redesign bullet sat in the
same released file describing the same function two different ways; the old
one is now false. Removed; the shipped panel is described once.

### F-102 (new, minor) — missing-data imputation was not new in 0.2.0.0

> Missing data imputation: `na_action = "impute"` integrates over missing
> values during MCMC sampling for ordinal, continuous, and mixed models.

`na_action = c("listwise", "impute")` is already a 0.1.6.3 formal
(`cran-0.1.6.3:R/bgm.R:402`) and imputation is announced in the 0.1.1
section of this same file. What is new is the coverage of the new model
classes. Restated as "now covers the new model classes … for continuous and
mixed models as well as ordinal ones."

### F-103 (new, minor) — a log10 figure in a natural-log package

The `verdicts()` entry quoted the calibration study as "within 0.25 of a
threshold on the log10 Bayes factor scale". Per F-039 the package converted
to natural log everywhere, and neither `R/verdicts.R` nor `man/verdicts.Rd`
contains the string `log10`; both state 0.58 nats
(`R/verdicts.R:286`, `man/verdicts.Rd:57`). NEWS was the last log10 residue.
Restated as 0.58 on the log Bayes factor scale, matching the shipped Rd
exactly.

### F-023 — closed; sentence added, mechanism verified in source

Verified rather than assumed. The discrete Gibbs conditional in
`simulate_mrf()`:

- HEAD, `src/mrf_simulation.cpp:97`:
  `rest_score += 2.0 * (obs - ref) * pairwise_safe(vertex, variable);`
- `cran-0.1.6.3:src/mrf_simulation.cpp:98`:
  `rest_score += (obs - ref) * pairwise_safe(vertex, variable);`

So a user-supplied `pairwise` matrix produces twice the coupling it did at
0.1.6.3. The new Breaking-changes bullet states the change, the user action
(halve a hand-specified matrix to reproduce 0.1.6.3 output), and the
exemption (a matrix taken from a 0.2.0.0 `bgm()` fit needs no adjustment,
because fit and simulator are now on one scale; the `variable_type =
"continuous"` branch takes a precision matrix and is unaffected).

The same factor-2 convention was confirmed for the ordinal sampler itself —
`src/models/omrf/omrf_model.cpp:618` `+ 2.0 * obs_other * delta` against
`cran-0.1.6.3:src/bgm/bgm_logp_and_grad.cpp:178` `+ obs_other * delta` —
which is the evidence behind the "stored on the association scale, about
half the 0.1.6.3 values" claim that was already in NEWS.

### F-029 — closed; seven of eleven items had no NEWS entry

Checklist with the citation for each, and where it now lives.

| F-029 item | Verified at | Pre-existing NEWS? | Now |
|---|---|---|---|
| `iter`/`warmup` 1e3 → 2e3 | `R/bgm.R:479-480` = `2e3`; `R/bgmCompare.R:199-200` = `2e3`; `cran-0.1.6.3:R/bgm.R:386-387` and `:R/bgmCompare.R` = `1e3` | **no** | Breaking changes, "Changed defaults" |
| NUTS `target_accept` 0.60/0.65 → 0.80 | `R/validate_sampler.R:147-155` (`nuts` = 0.80, `adaptive-metropolis` = 0.44, `gibbs` = `NA_real_`); `cran-0.1.6.3:R/bgm.R:463-467` (nuts 0.60) and `:R/bgmCompare.R:346-350` (nuts 0.65) | **no** | Breaking changes, with the runtime consequence stated |
| `hamiltonian-mc` + `hmc_num_leapfrogs` removed | zero hits for either string in `R/` at HEAD; `cran-0.1.6.3:R/bgm.R:402,404` had both | yes | Breaking changes (kept, merged into one bullet) |
| new formal `means_prior` | `R/bgm.R:483` (`normal_prior(scale = 1)`); introduced `586df5f9` (#91) | **no** (0 hits in NEWS) | New features, "Prior specification" |
| new formal `threshold_prior` | `R/bgm.R:482`, `R/bgmCompare.R:197` (`beta_prime_prior(0.5, 0.5)`); supersedes `main_alpha`/`main_beta` (`R/bgm.R:571+`) | **no** (0 hits) | New features, "Prior specification"; deprecation of the scalars in Breaking changes + Deprecated |
| new formal `delta` | `R/bgm.R:485` (`NULL`, resolved to `0.5*log(p)`); introduced `9480d46c` (#107) | mentioned only in passing | New features, "Prior specification" |
| new formal `progress_callback` | `R/bgm.R:498`, `R/bgmCompare.R:211`, threaded at `:474`; introduced `e8d078d7` (#88) | **no** (0 hits) | Other changes |
| new formal `difference_family` | `R/bgmCompare.R:194` `c("Cauchy","Normal")`, resolved `:424-425`; introduced `53441521` (#154) | **no** (0 hits) | New features, "Prior specification" |
| `extract_ess` default now RB, no warning on new fits | `R/extractor_functions.R:1539` `estimator = c("rb","mixt")`; deprecation fires only on an explicit `"mixt"` (`:1540`) | yes | Breaking changes, "Changed output and semantics" (kept, plus the legacy-fit fallback the Rd documents) |
| prediction/simulation on the ×2 association-scale convention | `src/mrf_simulation.cpp:97`; `src/models/omrf/omrf_model.cpp:203,618`; compare fix `d962d778` (#189) | partly (compare only) | Breaking changes; the `simulate_mrf()` gap is F-023 above |
| bgmCompare RB PIPs `NA` for unselected indicators | `R/extractor_functions.R:342` `out[(n0_visits + n1_visits) == 0] = NA_real_`; `:660-666` rb branch with the comment "difference indicators have no draws and average to NA"; `R/mcmc_summary.R:236-240` same rationale | **no** | Breaking changes, flagged as breaking an all-numeric assumption, with `estimator = "raw"` as the escape |

All eleven now have a line. The last one is also the item the brief routes to
the F-002 downstream memo; the NEWS bullet states it in the terms a
downstream consumer needs (previously all-numeric matrix, `is.na()` filter,
`"raw"` fallback).

---

## Evidence

### A. Claims table — 100% of pre-existing `# bgms 0.2.0.0` entries

98 bullets. IDs: `D*` = the "Changes since the 0.2.0.0 development build"
section, `B*` = Breaking changes, `N*` = New features, `O*` = Other changes,
`X*` = Bug fixes. Verdicts: **accurate** (carried, possibly reworded or
merged), **corrected**, **removed**.

Verification depth is stated per row. Where an entry cites a measured
simulation result (nats-to-gold, recall rates, timings), I verified that the
mechanism, argument, option, field or counter the sentence is *about* exists
and behaves as described, and did not re-run the study; those rows are marked
"mechanism verified, measurement not re-run". No such measurement was
altered.

#### Development-build section (D1–D37)

| ID | Entry | Verified against | Verdict |
|---|---|---|---|
| D1 | `plot_edge_posterior()` redrawn to the JASP panel | `R/plot_bgms.R:529-560` roxygen; `:534` signature | accurate — carried |
| D2 | Estimator-dependent panel (no Savage-Dickey under selection) | `R/plot_bgms.R` Details block, `:556-600` | accurate — carried, merged into D1 |
| D3 | Ruled-out edge gets a figure, not an error | `R/plot_bgms.R` (caption branch) | accurate — merged into D1 |
| D4 | Fixed: mixed cross-edge named by its discrete end, both orientations | fix to a 0.2.0.0-only function | accurate — removed (never shipped to a 0.1.6.3 user) |
| D5 | Refuses a `bgmCompare()` fit | `R/plot_bgms.R` guard | accurate — merged into D1 |
| D6 | `binwidth` deprecated and ignored | `R/plot_bgms.R:628-632` `lifecycle::is_present(binwidth)` → `deprecate_warn("0.2.0.0", ...)` | accurate — carried to Deprecated (live formal, warns) |
| D7 | Classic split-R-hat, df adjustment dropped | `src/mcmc_diagnostics.cpp:155-160,227-230` — `var_plus = (n-1)/n*W + B/n`, `Rhat = sqrt(var_plus/W)`, with the Brooks-Gelman/`coda::gelman.diag` note in the header comment | accurate — carried to Breaking changes (0.1.6.3 imported `coda`) |
| D8 | Windows + RcppParallel ≥ 6.0.0 not bitwise across core counts | `fa5c0463` (#176); DESCRIPTION `Suggests: coda`, no version pin | accurate — carried to Other changes |
| D9 | `interaction_prior` default `normal_prior(1)` in `bgm()`/`sample_ggm_prior()`; compare keeps Cauchy | `R/bgm.R:481`; `R/bgmCompare.R:197` `cauchy_prior(scale = 1)`; `sample_ggm_prior` default per 00c §2 | accurate — carried to Breaking changes, extended with the F-002/report-09 Q5 asymmetry warning |
| D10 | Trust gauge on by default; `bgms.zratio_gauge_sweeps` off switch; cost; ladder; surface cache/announce | `R/zratio_gauge.R`, `R/zratio_surfaces.R`, `R/bgm.R` all reference `bgms.zratio_gauge_sweeps`; `bgms.zratio_surface_cores` at `R/zratio_surfaces.R`; `bgms.correction_table_cache` at `R/correction_tables.R`, `R/zratio_surfaces.R`; `block_lo`/`block_hi` at `R/zratio_gauge.R` | accurate (mechanism verified, timings not re-run) — carried, split across three consolidated hierarchical bullets |
| D11 | Absolute-moment surface replaces the online calibrator; `calibration_window` removed; counters dropped | zero hits for `calibration_window`, `n_clamp`, `n_oracle`, `n_anchors` anywhere in `R/` or `src/` | accurate — carried; the "replaces the development build's calibrator" framing dropped (the calibrator never shipped) |
| D12 | `harm_pred` second alarm channel | `harm_pred`, `harm_flag`, `amplification`, `se_mcse`, `se_se` all in `R/zratio_gauge.R` | accurate (mechanism verified) — carried into the gauge bullet |
| D13 | Surface validated at gamma shapes 0.5, 1, 2 | superseded by D16 (range to 10) | accurate at the time — removed as superseded |
| D14 | Fixed: anchor helpers did not pass the cell's Gamma shape | subsystem absent from `cran-0.1.6.3` (`git ls-tree cran-0.1.6.3 \| grep zratio` = empty) | accurate — removed (never shipped) |
| D15 | Fixed: non-unit-shape surface built, cached, then ignored | same | accurate — removed (never shipped) |
| D16 | Shape range to 10; scored at 0.5/1/2/3/5; interior interpolated | `R/zratio_surfaces.R` deployment gate | accurate (mechanism verified) — carried, consolidated |
| D17 | Above shape 10 the isolated-edge route; bound 0.00028/0.00022 nats; `n_isolated` | `n_isolated` present in `R/zratio_gauge.R`, `src/models/ggm/zratio_engine.{h,cpp}`, `ggm_model.cpp`, `mixed_mrf_model.cpp`, `chain_runner.cpp` | accurate (mechanism verified) — carried, consolidated |
| D18 | Additive kernel zero-collapse; `n_collapsed`, `max_collapse_size`; 12–32-variable boundary | both counters present in `R/zratio_gauge.R` and five `src/` files; `?summarize_zratio_gauge` documents it | accurate (mechanism verified) — carried, consolidated (relates to F-025(a)) |
| D19 | `"hierarchical"` accepted where the choice is vacuous | `R/bgm_spec.R:308` `zratio_vacuous_spec_notice()` (report 00c §4) | accurate — carried |
| D20 | Certified constants range [0.5, 20]; warns above 20 | `R/zratio_tables.R` | accurate (mechanism verified) — carried |
| D21 | Anchor simulator now slice-samples the diagonal | internal to a subsystem absent at 0.1.6.3; no user-visible surface | accurate — removed (never shipped); the build timings it quotes are carried on the cache bullet |
| D22 | Fixed: quadrature grid accuracy at large shape | same subsystem | accurate — removed (never shipped) |
| D23 | Fixed: tabulated pair integrals ended before the range read | same subsystem. Note its clause "The plateau is present in released `bgms`" is **not true of 0.1.6.3** — the tables do not exist there | corrected by removal; see Open questions Q1 |
| D24 | Boundary-slope extension past the anchored range; `n_slope_floor` | `n_slope_floor` in `src/models/ggm/zratio_engine.{h,cpp}` + 3 more | accurate (mechanism verified) — carried, consolidated |
| D25 | Anchored range to 80 variables | `R/zratio_surfaces.R` anchor grid | accurate (mechanism verified) — carried, consolidated |
| D26 | Extrapolation notice separates warmup from retained; three counters | `n_pred_retained`, `n_extrap_retained`, `max_extrap_size_retained` in `R/zratio_gauge.R` + 4 `src/` files | accurate — carried |
| D27 | Surface anchors at the size cap | `R/zratio_surfaces.R` | accurate — carried, folded into the anchored-range bullet |
| D28 | Joint spec reports the realized edge prior `pi(Gamma)*Z(Gamma)` | `R/bgm_spec.R:256` `zratio_joint_realized_prior_notice()` | accurate — carried |
| D29 | `extract_inclusion_bf(log =)`, default `FALSE` | `R/extractor_functions.R:419` `function(bgms_object, log = FALSE)` | accurate — carried |
| D30 | Dev-build scripts silently change meaning on 0.2.0.0 | true, but addressed to develop-trackers, not 0.1.6.3 users | accurate — removed (see Open questions Q2) |
| D31 | Inclusion table RB `mcse`/`n_eff`/`Rhat` masked on constant RB draws, not zero flips | `R/mcmc_summary.R:273-274` `constant = is.finite(mat[,"sd"]) & mat[,"sd"] == 0` then mask | accurate — carried, restated for a 0.1.6.3 reader (who never saw the flip-count masking) |
| D32 | `n_eff_mixt` removed from the inclusion table; `estimator = "mixt"` deprecated; three grounds | inclusion table for a 0.2.0.0 fit is built by `summarize_rb_inclusion()` (`R/mcmc_summary.R:55,121`) whose columns are `mean, mcse, sd, n_eff, Rhat, n0->1, n1->0` (`:276-283`) — no `n_eff_mixt`. The surviving `n_eff_mixt` at `R/mcmc_summary.R:217,494` is the **legacy** path taken only when `raw[[1]]$rb_inclusion_samples` is `NULL`, i.e. a pre-0.2.0.0 fit | accurate — carried to Deprecated |
| D33 | `plot()` on a `bgmCompare()` fit | `R/plot_bgms.R:500-527` (`type = "groups"` panels, `"difference"`, `"centrality"`) | accurate — carried |
| D34 | `extract_centrality()` on `bgmCompare()` fits | `R/centrality.R:68` `group` formal | accurate — carried |
| D35 | `prior_sensitivity_check()` on `bgmCompare()` fits | `R/prior_sensitivity.R:329` compare unit with `difference_scale` | accurate (mechanism verified) — carried |
| D36 | `vary` argument | `R/prior_sensitivity.R:203` `vary = c("auto","slab","slab-and-diagonal")` (00c §2) | accurate — carried, folded into the sensitivity feature |
| D37 | Noise band no longer `Inf` on a saturated edge | `R/prior_sensitivity.R` band computation | accurate — carried as a clause in the sensitivity feature (the `Inf` bug is dev-build-only) |

#### Breaking changes (B1–B10)

| ID | Entry | Verified against | Verdict |
|---|---|---|---|
| B1 | `update_method = "hamiltonian-mc"` removed | 0 hits in `R/`; `cran-0.1.6.3:R/bgm.R:402` had it; `8bbf6b9f` (#93) | accurate — carried |
| B2 | `hmc_num_leapfrogs` removed | 0 hits in `R/`; `cran-0.1.6.3:R/bgm.R:404` `= 100` | accurate — carried, merged with B1 |
| B3 | `standardize`: `FALSE` warns, `TRUE` errors | `R/bgm.R:589-609` — `deprecate_warn` on `isFALSE`, `deprecate_stop` otherwise, with the `interaction_prior`/`difference_scale` pointer in the message | accurate — carried |
| B4 | Ordinal pairwise on the association scale, half the 0.1.6.3 sigma | `src/models/omrf/omrf_model.cpp:203,618` (`2.0 *`) vs `cran-0.1.6.3:src/bgm/bgm_logp_and_grad.cpp:178` (no factor) | accurate — carried, with the "no automatic conversion" user action added |
| B5 | "Default `pairwise_scale` changed from 2.5 to 1" | `R/bgm.R:501` (bare deprecated formal, no default), `:561-569` (translates to `cauchy_prior`) | **FALSE — corrected**; see F-100 |
| B6 | `bgmCompare()` stores pairwise on the association scale; `omega*x` vs `2*omega*x` | `d962d778` (#189); F-015 records the landing | accurate — carried, promoted to lead the section, with "refit rather than rescale" stated as the user action |
| B7 | compare `predict()`/`simulate()` were wrong as a consequence; Wenchuan numbers | `d962d778`; F-015 | accurate (measurement not re-run) — carried |
| B8 | `interaction_prior`/`difference_scale` keep numeric defaults; same number is now a tighter prior; verdicts scale-contingent | `R/bgmCompare.R:193,197`; `cran-0.1.6.3:R/bgmCompare.R` `difference_scale = 1` unchanged; caveat shipped in `vignettes/comparison.Rmd:70-75` (F-015) | accurate — carried, merged into B6 |
| B9 | `extract_category_thresholds()` deprecated | `R/extractor_functions.R:1037-1044` `deprecate_warn` → `extract_main_effects()` | accurate — carried |
| B10 | `extract_ess()` reports the RB ESS | `R/extractor_functions.R:1539-1556` | accurate — carried, extended with the legacy-fit fallback |

#### New features (N1–N34)

| ID | Entry | Verified against | Verdict |
|---|---|---|---|
| N1 | `verdicts()` + fragility flag + 37,010-fit study | `R/verdicts.R:263`; distance figure at `R/verdicts.R:286` / `man/verdicts.Rd:57` = **0.58 nats**, not 0.25 log10 | accurate except the units — **corrected** (F-103) |
| N2 | `plot()` on a `bgm()` fit | `R/plot_bgms.R` | accurate — carried |
| N3 | `plot_edge_posterior()` "stem at zero … probability per weight bin" | superseded by D1; `binwidth` deprecated at `R/plot_bgms.R:628` | **FALSE at HEAD — removed** (F-101) |
| N4 | `extract_centrality()` | `R/centrality.R:68`, `measure = "strength"` only | accurate — carried |
| N5 | `calibration_check()` isotonic reliability curve | `R/calibration_check.R:277` | accurate — carried |
| N6 | Continuous panels via the PIT; `curves$pav` renamed `curve` | `R/calibration_check.R` `kind` column | accurate — carried, merged into N5; **extended**: `S3method(calibration_check, bgmCompare)` is in NAMESPACE (F-021 shipped it as `82edaa4c`), which no NEWS entry mentioned |
| N7 | New vignette "Checking your fitted model" | `vignettes/checking-your-model.Rmd` (00c §6) | accurate — carried |
| N8 | `prior_sensitivity_check()` anchored curve | see F-001 table above | accurate — carried |
| N9 | `sample_graph_prior()` | `R/sample_graph_prior.R:230` | accurate — carried |
| N10 | `sample_sbm_prior()` | `R/sample_graph_prior.R:379` | accurate — carried |
| N11 | Correction table announces only on build; progress bar follows `display_progress` | `R/correction_tables.R` | accurate — carried |
| N12 | GGM via `variable_type = "continuous"`, theta-space NUTS | `d1fe7522` (#111), `9fcb5c11` (#78) | accurate — carried |
| N13 | GGM Gibbs `update_method = "gibbs"` | `R/validate_sampler.R:116` GGM-only gate; `3d18910f` (#146) | accurate — carried |
| N14 | Gibbs warmup staging 15% / 85% | `src/mcmc/execution/warmup_schedule.h:173` `gibbs_settle_fraction = 0.15`, used at `:81` | accurate — carried |
| N15 | Standardized-frame diagonal prior; default `exponential_prior(eta = 1)`; `beta_prime_prior()` needs `rate` | `R/priors.R:208,297` (`rate` XOR `eta`, `validate_rate_eta()` `:240`); `R/bgm.R:484` | accurate — carried, and the **default change** promoted to Breaking changes (it changes numbers for a 0.1.6.3 user at a non-unit slab scale) |
| N16 | Mixed MRF models | `R/build_spec.R:216` `build_spec_mixed_mrf()` | accurate — carried |
| N17 | "Missing data imputation: `na_action = "impute"` …" | `cran-0.1.6.3:R/bgm.R:402` already has `na_action = c("listwise","impute")`; 0.1.1 section of this file announces it | misleading — **corrected** (F-102) |
| N18 | `extract_precision()` | `R/extractor_functions.R:1682` | accurate — carried |
| N19 | `extract_partial_correlations()` | `R/extractor_functions.R:1753` | accurate — carried |
| N20 | `extract_log_odds()` | `R/extractor_functions.R:1806` | accurate — carried |
| N21 | `extract_main_effects()` | `R/extractor_functions.R:932` | accurate — carried |
| N22 | RB inclusion machinery, `J = gamma + (1-2 gamma) alpha`, −745 floor, prior odds divided out, RB is canonical | `R/extractor_functions.R:329-343` (`rb_log_odds_from_counts`); `R/mcmc_summary.R:241-283`; `5fc618d3` (#182) | accurate — carried |
| N23 | NUTS `accept_prob` + `mean_accept_prob` | `R/diagnostics_nuts.R` | accurate — carried to Other changes |
| N24 | `sample_ggm_prior()` gains `"gibbs"` and `"beta-bernoulli"` | `R/sample_ggm_prior.R:190` signature (00c §2) | accurate — carried |
| N25 | Corrected BB inclusion updates on GGM (0.37 vs 0.5 at p = 5) | `R/correction_tables.R`; `93690fc8` (#151) | accurate (mechanism verified) — carried |
| N26 | Corrected SBM updates (1.0 vs 1.87 blocks at p = 20) | same table; `6ebcfc9a` (#99), `93690fc8` | accurate (mechanism verified) — carried |
| N27 | Correction extends to mixed models; <2 continuous → plain updates; ==2 → `sbm_prior()` warns | `R/build_spec.R:216+`, `R/correction_tables.R` | accurate — carried |
| N28 | Hierarchical spec: `p(Gamma) p(K\|Gamma)`, per-graph normalization, requirements, mixed handling, +0.002 check | `cfa17068` (#153), `cb8e81c1` (#193); `R/bgm.R:489` `precision_graph_prior = c("joint","hierarchical")` | accurate (mechanism verified) — carried as the lead hierarchical bullet |
| N29 | Trust gauge `flip_rate` > 1% flag; `fit$zratio_diag` | `R/zratio_gauge.R:154` `summarize_zratio_gauge(chains, threshold = 0.01, ...)`; `zratio_diag` in `R/build_output_bgm.R`, `R/class_s7.R` | accurate — carried, merged with D10/D12 |
| N30 | Constants built at unit slab scale; the 0.292-for-0.30 bias it fixed | `42c277cd` (#162) | accurate — carried; the "previously" clause dropped (never shipped) |
| N31 | Gamma diagonal of any shape; constants cell `(delta, eta, shape, slab)` | `5d0a2120` (#168), `cb8e81c1` (#193) | accurate — carried, consolidated with D16/D17 |
| N32 | Cauchy slab gets its own constants; 0.247-for-0.30 bias fixed | `cef30b88` (#159), `42c277cd` (#162) | accurate — carried; "previously" clause dropped |
| N33 | `warmup_incomplete` flag now printed | `R/diagnostics_nuts.R` | accurate — carried to Other changes |
| N34 | `extract_prior_inclusion_probabilities()`; joint-spec 0.37 / 0.27 figures | `R/extract_prior_inclusion_probabilities.R:443`; F-036 confirms the hierarchical branch now returns the nominal prior (`ff29d177`) | accurate — carried |

#### Other changes (O1–O6) and Bug fixes (X1–X11)

| ID | Entry | Verified against | Verdict |
|---|---|---|---|
| O1 | S7 objects; `$`/`[[`/`names()` still work; easybgm shim | `R/class_s7.R:29,149`; shims at `R/methods_bgms.R:303,326,351` | accurate — carried to Breaking changes (a class change a 0.1.6.3 user sees) |
| O2 | C++ backend refactor, unified hierarchy | `src/models/` layout; `9fcb5c11` (#78) | accurate — carried; "NUTS/HMC infrastructure" reworded to "NUTS" since HMC is gone |
| O3 | NUTS multinomial candidate weighting | `53c045c6` (#96) | accurate — carried |
| O4 | Stage-2 warmup windowing matches Stan | `src/mcmc/execution/warmup_schedule.h`; `53c045c6` | accurate — carried |
| O5 | `coda` dropped from Imports | DESCRIPTION: `coda` now in Suggests only; `4f70155a` (#81) | accurate — carried |
| O6 | `$`/`[[` trigger lazy diagnostics | `R/class_s7.R` property getters calling `ensure_summaries()` | accurate — carried, merged with O5 |
| X1 | simulate/predict category-scale mismatch | `54cea52b` (#152), `965325e7` (#116), `10ac178b` (#114) | accurate — carried |
| X2 | Mixed PIP extractor block order | `R/extractor_functions.R` block filler | accurate — carried |
| X3 | SBM number-of-blocks summary (shifted vs truncated Poisson) | `R/mcmc_summary_sbm.R` | accurate — carried |
| X4 | Mixed cross-indicator upper-triangle asymmetry | `4b55d068` (#150) | accurate — carried |
| X5 | `delta = NULL` default counted blume-capel variables | `4b55d068` (#150) | accurate — carried |
| X6 | Cholesky downdate guard | `a7130ebf` (#134) | accurate — carried |
| X7 | Alpine/musl `<tbb/global_control.h>` | `src/mrf_simulation.cpp` includes | accurate — carried |
| X8 | Stale gradient cache after imputation | `53441521` (#154) | accurate — carried, merged with X9 |
| X9 | Stale observation transpose after imputation | `53441521` (#154) | accurate — carried, merged with X8 |
| X10 | `target_accept` not passed to lower-level NUTS | `53c045c6` (#96) | accurate — carried |
| X11 | Acceptance-probability accumulation overwritten by the last subtree | `53c045c6` (#96) | accurate — carried |

**Totals.** 98 entries: 89 accurate (74 carried as-is or reworded, 15 merged
into a consolidated entry), 4 corrected (B5, N1, N3, N17 — of which B5 and N3
were false at HEAD), 8 removed as never-shipped development-build material
(D4, D14, D15, D21, D22, D23, D30, plus D13 as superseded). Some rows are
counted in two of those buckets where an entry was both merged and corrected.

### B. Mapping table — entries added by this brief

Every new entry maps to a commit, PR, or FINDINGS row.

| New entry | Maps to |
|---|---|
| `simulate_mrf()` factor-2 input scale | F-023; `src/mrf_simulation.cpp:97` vs `cran-0.1.6.3:…:98`; `243c8def` (#84) |
| `iter`/`warmup` 1e3 → 2e3 | F-029; `R/bgm.R:479-480`, `R/bgmCompare.R:199-200` |
| NUTS `target_accept` → 0.80 | F-029; `R/validate_sampler.R:147-155` |
| Prior-object API as a unit (8 constructors, 5 replaced formals) | F-029; `586df5f9` (#91); report 00c §1 (the 8 constructors are 8 of the 23 new exports) |
| `means_prior` | F-029; `586df5f9` (#91); `R/bgm.R:483` |
| `threshold_prior` | F-029; `586df5f9` (#91); `R/bgm.R:482` |
| `delta` | F-029; `9480d46c` (#107); `R/bgm.R:485` |
| `difference_family` | F-029; `53441521` (#154); `R/bgmCompare.R:194` |
| `progress_callback` | F-029; `e8d078d7` (#88); `R/bgm.R:498`, `R/bgmCompare.R:211` |
| bgmCompare RB PIPs return `NA` | F-029; `R/extractor_functions.R:342,660-666`; `5fc618d3` (#182) |
| `posterior_mean_indicator` is now the RB estimate | F-029 (extract_ess sibling); `5fc618d3` (#182); `R/mcmc_summary.R:55,121` |
| `summary()` gains a `quadratic` element | report 00c §3; `R/methods_bgms.R:104,117,158-161` |
| `precision_scale_prior = exponential_prior(eta = 1)` promoted to a default change | F-029 (defaults freeze) / F-002; `5d0a2120` (#168), `8b60ca42` (#148); `R/bgm.R:484` |
| `calibration_check()` has a `bgmCompare` method | F-021; `82edaa4c`, merged `4969f843`; `NAMESPACE:7` |
| `bgmCompare(difference_probability =)` deprecated | report 00c §3; `R/bgmCompare.R:196` (bare formal), roxygen `:61` |
| `bgm()`/`bgmCompare()` interaction-prior default asymmetry called out | F-002 / report 09 Q5; `R/bgm.R:481` vs `R/bgmCompare.R:197` |
| Natural-log convention stated once at the top | F-039 |

### C. Verification gate

1. **Claims table covers 100% of pre-existing NEWS entries** — 98 of 98,
   table A.
2. **F-029 checklist** — all eleven items have an entry, each cited; table in
   the F-029 section.
3. **Parse.** `tools:::news2Rd()` is the *plain-text* NEWS converter and does
   not read a markdown `NEWS.md`; it fails identically on the pre-change file
   and on mine (`ERR: No news found in given file using package default
   format`, reproduced on `HEAD:NEWS.md` in a package-shaped temp directory),
   so it is not the applicable tool and not a regression. The reader R
   actually uses for a markdown `NEWS.md` — `tools:::.build_news_db_from_package_NEWS_md()`,
   the same one behind `utils::news()` on an installed package — parses the
   file cleanly, run against the file alone with no package build:

   ```
   rows: 33   versions: 12   NA version/category rows: 0
   0.2.0.0 sections:  (preamble) 1 | Breaking changes 1 | New features 1 |
                      Other changes 1 | Bug fixes 1 | Deprecated 1
   ```

   (Baseline for comparison: 32 rows. The extra row is the 0.2.0.0 preamble
   paragraph, which parses as an uncategorised entry — the same shape the
   existing `# bgms 0.1.4.1` section already has.) R 4.6.0.

---

## Open questions

**Q1 — D23's claim that the tabulated-pair-integral plateau "is present in
released `bgms`".** The whole z-ratio subsystem is absent from
`cran-0.1.6.3` (`git ls-tree -r cran-0.1.6.3 | grep zratio` is empty), so
that sentence cannot be about the released version. It reads as though it
were written against an internal reference build or the companion GGM-paper
implementation. I removed the entry with the other never-shipped fixes. If
"released" means something the maintainer intends users to be able to
identify, the sentence needs a subject and should come back.

**Q2 — how much to say to development-build users.** The old D29/D30 pair
warned that a script written against a 0.2.0.0 development build silently
changes meaning (`extract_inclusion_bf()` returning the Bayes factor rather
than its log). A 0.1.6.3 user has no such script, so I dropped the warning
paragraph and kept only the `log =` description. But JASP and easybgm track
`develop`, and the F-002 memo is the natural place for it. Confirm the memo
carries it, or say the word and I will restore a short "if you tracked the
development build" note under Deprecated.

**Q3 — how much hierarchical-prior detail belongs in a released NEWS.** The
pre-existing text spent roughly 15 dense paragraphs (about 40% of the whole
0.2.0.0 section) on a non-default specification. I consolidated it to nine
bullets under New features, keeping every fence, every validated range, every
measured envelope, and every diagnostic counter, and dropping only the
"was X, now Y relative to the development build" framing. That is still the
longest block in the file. If you would rather it lived in
`?bgm` / the diagnostics vignette with a three-sentence NEWS pointer, that is
a straightforward cut — but it removes the measured accuracy claims from the
artefact a reviewer reads first, so I did not make that call.

**Q4 — the `pairwise_scale` deprecation path.** `pairwise_scale = s`
translates to `cauchy_prior(scale = s)` (`R/bgm.R:567`), which preserves the
0.1.6.3 *family*, so a deprecated call and a call at the new default are
different models. I stated that explicitly. Confirm it is the intent (rather
than, say, translating to the new default family), because it is the kind of
thing a user will only discover from this sentence.

**Q5 — `simulate_mrf()` user action.** I wrote "halve a hand-specified
`pairwise` matrix to reproduce 0.1.6.3 output". That is the arithmetic, and I
verified the two source lines, but I did not run a paired simulation to
confirm the distributions match. If you want that asserted rather than
derived, it is a seconds-long check at p = 5 and I will run it.

**Q6 — DESCRIPTION `Date:` (F-003) is untouched here.** It stays stale at
`2026-03-26`; per F-003 the bump belongs in the submission commit, not this
one. Flagging so it is not assumed covered by "the NEWS brief".

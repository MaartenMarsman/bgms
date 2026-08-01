# Report 06 — bgmCompare defect + units batch (Opus agent)

Branch `fix/bgmcompare-batch`, based on `origin/develop` at `fe006909`
(contains merge `b04dbd06`; the local Dropbox checkout was behind, so the
branch was cut from the fetched remote). Worktree `~/bgms-review/wt-fix3`;
installed builds `~/bgms-review/lib-fix3` (branch) and `lib-fix3-base`
(base, for before/after). Ten commits, one per item plus two follow-ups;
nothing pushed.

---

## Item 1 first: F-061 diagnosis and blast radius

**Root cause (confirmed on both Boredom fits).** The run-to-run noise
estimator `own_lbf` inside `prior_sensitivity_engine`
(`R/prior_sensitivity.R`) averaged `anchor_draws(fit)$gamma` — the
**per-gated-parameter** indicator expansion — and divided by the
**per-indicator** prior odds. The actual layouts, established empirically:

| object | default fit | main-sel fit |
|---|---|---|
| `parameter_names$indicator` | 36 | 36 |
| `raw$indicator` draw columns | 36 | 36 |
| `raw$rb_inclusion` columns | 36 (mains `NaN`) | 36 |
| `anchor_draws()$gamma` / `$theta` | **76** | **76** |
| `anchor_draws()$indicator` | 36 | 36 |

76 = 28 pairwise differences + 48 main-threshold difference parameters
(8 variables x 6 thresholds); each main indicator is repeated once per
threshold it gates. 76 against 36 is not a multiple, hence the two
recycling warnings — one per fit of the replicate pair. The engine's other
vectors (`refit_edge_stats`, the curve via `$indicator`, `lbf_of`) were all
36/36-aligned; the trace in the brief was right about the mechanism but the
misalignment lives in the noise estimator, not in `refit_edge_stats`.

**What was numerically wrong.** Because the difference prior is
exchangeable, `prior_odds` is a constant vector, so recycling corrupted no
individual value — it corrupted the *set*: the wobble quantile `q95` pooled
one entry per gated parameter instead of one per indicator (main-effect
differences overweighted sixfold when selected), `$wobble$per_edge` was a
76-vector misaligned with the 36-row edges table, and the printed noise
band and the mover rule's threshold (`max(tolerance, 2*MCSE, q95)`) fed on
it. The curve, the chosen-scale columns, every verdict, and the anchor
tables were never touched. On the Boredom fits the movers happened not to
reclassify after the fix, but that is luck of this dataset, not of the code.

**Silent variant already in the suite.** The existing slow-tier compare
test uses 4 variables: 22 gated parameters against 10 indicators on
Wenchuan (30 vs 10 on Boredom-shaped 7-category data) — when the counts
divide evenly R recycles *without* a warning, so the shipped test never
saw the defect. The new regression test pins both the warning and the
yardstick length.

**Fix** (`7e31a491`): `own_lbf` averages `d$indicator`, the unit every
other quantity in the check is keyed to.

**Blast radius on single-group fits: none.** `anchor_draws()` for a `bgms`
fit returns `indicator = gamma` — the same object
(`R/anchor_curve.R:38-43`) — so the change is an exact no-op for `bgm()`
sensitivity results; no released single-group number ever passed through
the defective path. Confirmed empirically: the item-8 invariance run
reproduced the single-group check exactly.

**Verification (all four parts + amendment):**

1. Both maintainer fits (`seed = 1`, defaults / `main_difference_selection
   = TRUE`) run `prior_sensitivity_check()` with **zero warnings**;
   `length($wobble$per_edge)` = 36 = `nrow($edges)` on both.
2. Ground-truth cross-check: fresh full-defaults refits at the off-anchor
   multipliers 0.8 and 1.3 (plus an independent repeat at each, for the
   noise scale):

   ```
   multiplier 0.80 (grid 0.795): max|curve - refit| = 0.0451, refit-vs-refit noise = 0.0503, ratio = 0.90
   multiplier 1.30 (grid 1.316): max|curve - refit| = 0.0798, refit-vs-refit noise = 0.0680, ratio = 1.17
   ```

   The reweighted curve agrees with brute-force refits within their own
   run-to-run spread — comfortably inside the ~4x band the single-group
   path is held to.
3. The maintainer's anomaly **survives the fix and is real**:
   `sit_around-half_dead_dull` rises monotonically from 1.88 to 3.79 nats
   across 0.4x-2.5x (undecided -> presence). The per-anchor verdict columns
   are each anchor fit's own exact Rao-Blackwellized analysis — no
   importance reweighting involved — and they show the same trend
   (undecided at 0.4x, presence from 0.63x on), and part 2's refits bracket
   the curve. Evidence for this difference genuinely increases with the
   difference scale; with a wider slab the other five movers drain toward
   absence while this one strengthens. After plot:
   `assets/sensitivity_compare_fixed.pdf` (before: `sensitivity_compare.pdf`).
4. Blast radius: above.

**F-049 gate re-run (amendment 2) — still fails, and the 20-seed
diagnostic says why.** On the branch build, the shipped construction gives

```
identity gap (must be < 0.02): 0.00323
ESS at 2x (must be > 400):     3264.5
noise (refit vs refit at 2x):  0.00505
gap (reweight vs refit at 2x): 0.0255
gap/noise ratio (must be < 4): 5.052   -> FAIL (unchanged from the shipped 5.05)
```

F-061 was never in this code path, and the lead's layout suspicion is
refuted: names, RB columns, and indicator draws all share one 10-column
layout on this fit (verified above at 36/36 for Boredom, 10/10 for
Wenchuan). The pre-approved diagnostic — 20 fresh seed triples, each one
reweighting prediction plus two independent refits — gives:

```
ratios sorted: 1.02, 1.43, 1.52, 1.69, 2.05, 2.06, 2.18, 2.41, 2.43, 2.70,
               2.85, 2.90, 2.93, 3.11, 3.26, 3.58, 4.10, 4.35, 5.17, 6.72
median 2.78, mean 2.92, share above the x4 gate: 4/20
shipped seed's 5.05 = 90th percentile of this distribution
```

Reading: **both effects are present**. The gap is systematic — every one
of the 20 seeds has gap > noise (gaps 0.006-0.017 pip at the worst edge
against noises 0.002-0.006), median ratio ~2.8 — so a 2x extrapolation
carries a real, small importance-reweighting bias of roughly 0.01 pip at
the worst pairwise difference. On top of that the pass/fail is seed luck:
the gate trips on ~20% of seeds, and the shipped seed sits in that tail
(q90). The two lowest-ESS triples produced the two worst ratios, which
fits ordinary importance-sampling error growing with extrapolation. Per
the amendment the test is untouched; the call on the gate (widen the
tolerance, pool the noise estimate over more refit pairs, or accept and
document the ~0.01-pip bound at 2x) stays with the maintainer.

---

## What was done, per item

| Item | Finding | Commit(s) | One line |
|---|---|---|---|
| 1 | F-061 | `7e31a491` | noise yardstick reads the indicator draws |
| 2 | F-057 | `fd191193` | `match()` recode at both sites; fixtures pass `language` directly |
| 3 | F-056 | `1d53bf39` | `data_check(x)` before `ncol(x)` in `bgmCompare()` |
| 4 | F-060 | `ded8418d` | print excludes unselected mains from counts/table; scale caveat; "Consider a longer run." |
| 5 | F-021 | `82edaa4c` (+ `2934f3b4`) | **full per-group `calibration_check.bgmCompare` landed** (not the stub) |
| 6 | F-063 | `b0398b11` | difference-centrality plot refused with a reasoned message; per-group kept |
| 7 | F-062 | `1b99e4cd` | legend strip below the panel; pie rings for main-difference pips (d implemented); wording |
| 8 | units | `4cc8ce41` | sensitivity surface in natural log; constants re-expressed exactly |
| 9 | F-065 | `8019134e` | progress bar renders at iteration 1; the silence was the first-50-iteration render gate |
| — | follow-up | `1508a888` | arguments field whitelist admits `baseline_category` |

Details that go beyond the brief's text:

**Item 2.** `unique()`/`match()` keeps first-appearance numbering exactly,
so integer indicators keep their group numbers. One consequence worth
knowing: Boredom is stored fr-first (rows 1-490 fr), so the fixtures that
previously hand-coerced via `as.integer(as.factor(...))` (alphabetical:
en = 1) now number fr = 1. The four fixture consumers are purely
structural (dims, ranges, symmetry), so nothing else moved. The natural
shipped-data call `bgmCompare(Boredom[,-1], group_indicator =
Boredom$language)` is the new slow-tier regression test; character,
factor, integer and 0/1 codings are asserted to agree on membership at the
spec level, the post-listwise re-recode is covered, and a vector `y` gets
`data_check`'s message ("y must be a matrix or data.frame.").

**Item 4.** The unselected main rows stay in the *returned* table (the
plot machinery and any user code indexing 36 rows keep working); only the
printed counts and body exclude them, via a `unselected_main` attribute
the compare method attaches. Genuine never-updated-but-selected indicators
keep the old "never updated" line. The scale-contingency caveat prints on
every compare-fit verdict table, worded from the comparison vignette. New
attributes are documented in `?verdicts`; snapshot tests cover both prints
and the footer on synthetic (platform-stable) tables.

**Item 5 (amended) — the full method shipped.** Reconciliation onto the
post-`0131bcdc` internals, as instructed: no second decode path —
`fitted_observed_data()` gets a compare branch that adds the Blume-Capel
baseline back (compare stores BC columns as `original - shift - baseline`,
which the shipped decode did not know) and then reuses
`recode_simulated_to_original()`; `discrete_category_index()` learns
bgmCompare's named many-to-one recode lookup; the per-panel builders are
extracted (`calibration_panel`, `pav_panels`, `calibration_result`) and
shared with the `bgms` method; the plot keeps the F-040 pagination and
pages over variable-by-group panels. User-supplied `newdata` arrives in
input order and is aligned to the fit's group-sorted internal order by the
same stable permutation the spec build used.

While wiring this I found and fixed a **new wrong-numbers defect**
(NEW-1 below): `build_arguments_compare()` never stored
`baseline_category`, so `predict.bgmCompare()` centered every Blume-Capel
term at 0 — up to **0.23 probability error** on a Boredom BC fit
(`2934f3b4`, with a manual-reference test that pins the sampler's
convention).

Verification at the F-035 standard: round-trip tests for non-contiguous
ordinal scores (2/5/9) and BC baselines; a layout-divergent regression test
on a many-to-one collapsed lookup; per-group curves on the Boredom fit run
clean, and per group the mean predicted category probabilities match the
observed margins to **0.0070 (group 1) / 0.0098 (group 2)** in-sample.
Print names the groups; `plot()` output in
`assets/calibration_compare_fixed.pdf`.

**Item 6.** What the panels drew: `plot(fit, type = "centrality")` (default
`group = 1`) draws group 1's posterior strength centrality — per-node sum
of absolute model-averaged edge weights, mean and 95% interval, rebuilt
per draw as `baseline + P %*% differences`; `group = c(1, 2)` drew the
per-draw *difference* in that quantity against zero. The difference
display is what the maintainer rejected, and it is now refused (in
`plot.bgmCompare` and on a plotted difference-centrality object) with the
reason: strength sums absolute weights, so a between-group difference
conflates which edges differ with how signs cancel — a zero can be
identical networks or compensating differences. The extractor keeps the
difference path; `summary()` still reports it; `?extract_centrality`
carries the interpretation caveat.

**Item 7.** (a) Established: the difference network **already encodes
verdicts**, the bgm convention — presence solid and sign-colored,
undecided dotted grey, absence undrawn, all at the same
`evidence_threshold` as `verdicts()`; edge *width* carries the
posterior-mean difference magnitude, exactly as bgm's width carries the
model-averaged weight. So the answer to "differences or Bayes factors?"
is: the drawn/dotted/absent channel is the Bayes-factor verdict, width is
the magnitude — one visual language across both functions already; no
encoding change was needed or made. (b) The legend overlapped the nodes;
the panels now reserve a bottom strip via qgraph's `mar` (which extends
the coordinate range) and the subtitle + legend draw there — verified
renders at default size, `assets/compare_plots_fixed.pdf`. (c) Banner is
now "main-effect differences not under selection (main_difference_selection
= FALSE)", matching item 4's line. (d) **Implemented** — it stayed small:
qgraph's `pie` channel draws each node's ring filled to its
main-difference indicator's posterior inclusion probability (full ring =
1, half = 0.5, as proposed), colored by verdict, replacing the square-node
encoding; the group panels carry the same rings, and the legend entry
reads "node ring: P(main-effect difference); full ring = 1".

**Item 8.** All sites from the report-03 §5 inventory converted:
`chosen_scale_log10_bf -> chosen_scale_log_bf`, `$log10_bf -> $log_bf`,
`$log10_bf_mcse -> $log_bf_mcse`; tolerance default `0.5 -> 0.5 * log(10)`,
threshold-relevant window `3 -> 3 * log(10)`, plot cap now reads
ln(1e6) ≈ 13.8 (it falls out of the same 1e-6 pip clamp); axis, zone
labels, print lines, Rd text, the F-051 details block, and the vignette
all in nats. No deprecation duplicate: 0.2.0.0 has not shipped, matching
the `log10_bf -> log_bf` precedent in `20b56996`. **Invariance check
(mandatory) — exact.** On a fixed single-group fit and a fixed compare
fit, before vs after: every verdict, every mover/insufficient
classification, every stability interval identical; every curve equals the
old curve times ln(10) to 1.8e-15; the wobble q95 rescales exactly.

**Item 9.** The silence is not initialization. Timing the phases: R-side
spec build is negligible (an adaptive-metropolis run at `iter = 10`
completes in 0.1 s total), the NUTS step-size heuristic is 4 ms, and the
whole cost sits in the first sampling iterations: the progress manager
polled only every 50 main-thread iterations, and the earliest warmup
iterations are the slowest NUTS takes (step size not yet adapted, deepest
trees). Measured on the maintainer's fit: **6.8 s** from call to first
bar render on the base build; immediate on the branch. Fix: poll on the
first main-thread update (print throttle backdated so it draws); all
other cadences — 50-iteration polling, half-second print throttle,
interrupt checks — unchanged. All samplers share the manager, so `bgm()`
gains the same immediate feedback.

---

## Findings (new, severity-tagged)

- **NEW-1 (major, wrong numbers, FIXED `2934f3b4`):**
  `predict.bgmCompare()` centered Blume-Capel variables at baseline 0
  because `build_arguments_compare()` never stored the baseline. Both the
  target variable's own quadratic term and the *rest scores of every
  neighbouring variable* were mis-centered, so predictions were wrong
  whenever any BC variable was in a compare fit — up to 0.23 in category
  probability on a Boredom BC fit. Only new fits benefit (arguments are
  stored at build time); nothing shipped, since 0.2.0.0 is unreleased.
- **NEW-2 (minor, analysis only — decision with the maintainer):** the
  F-049 gate fails for two compounding reasons: a genuine small (~0.01
  pip) systematic reweighting bias at a 2x extrapolation (median gap/noise
  2.78 over 20 seeds), and a fragile threshold (the x4-noise gate trips on
  ~20% of seeds; the shipped seed is at q90). Full numbers above.
- **NEW-3 (observation, no change):** in-sample calibration of the Boredom
  compare fit leaves several variables outside the 95% consistency band
  for up to ~30% of the curve (worst: `entertain`, `loose_ends` in group 1,
  max_dev ≈ 0.17-0.19). The margins match to <0.01, so this is
  conditional-shape miscalibration, not marginal bias. Worth the
  maintainer's eye now that the check exists for compare fits.
- **NEW-4 (cosmetic, not addressed):** the sensitivity plot's right-margin
  edge labels can clip at the device edge at default width (pre-existing;
  visible in both the before and after PDFs).

## Evidence

- **Before/after consoles.** Before (maintainer's session, report 07): two
  `longer object length is not a multiple of shorter object length`
  warnings from `prior_sensitivity_check()` on both fits; verdicts print
  interleaving `NaN` mains counted in "36 indicators". After (branch
  build, same calls): zero warnings, `presence 2 | undecided 9 | absence
  17 (28 indicators)` plus "Main-effect differences are not under
  selection (main_difference_selection = FALSE).", the vignette's
  scale-contingency caveat, and "Consider a longer run." — full transcript
  in the verification log (`verify-f061.log`, scratchpad).
- **Plots.** `assets/sensitivity_compare.pdf` (maintainer, before) vs
  `assets/sensitivity_compare_fixed.pdf` (after, natural-log axis);
  `assets/compare_plots_fixed.pdf` (difference network + rings, legend in
  its strip); `assets/calibration_compare_fixed.pdf`.
- **Gate 1 (regression tests fail on base, pass on branch), spot-checked
  one per item:** F-057: base errors `'bin' must be numeric or a factor`
  (character and factor alike), branch fits; F-056: base `non-numeric
  matrix extent`, branch gives the validator's message; F-061: base emits
  2 recycling warnings and a 22-long yardstick for 10 edges, branch 0 and
  10/10; F-060: the two updated print expectations fail on the old
  wording; item 5: base has no method (`UseMethod` error) and
  `arguments$baseline_category` is NULL on a base-build BC fit.
- **Gate 2 (full slow tier):** `TOTALS: failed 1 warnings 0 skipped 7
  passed 8823`. The one failure is the F-049 test, the only failure the
  amended gate permits ("the difference-scale reweighting reproduces a
  refit at that scale", now at `test-prior-sensitivity.R:460` after this
  batch's test insertions shifted it from line 432): `Expected
  max(abs(rw$pip[2, pairwise] - pip_of(f2)[pairwise])) < 4 * noise.
  Actual comparison: 0.0255 >= 0.0202`. The 7 skips are the pre-existing
  "golden fixtures not found" skips, unchanged by this batch.
- **Gate 3 (CRAN-mode suite):** `TOTALS: failed 0 warnings 0 skipped 240
  passed 7391` (after one legitimate follow-up: the arguments field
  whitelist gained `baseline_category`, `1508a888`).
- **Gate 4 (R CMD check --as-cran on a git-archive tarball):** `Status: 2
  NOTEs` — the two baseline NOTEs only ("The Date field is over a month
  old"; "'tidy' doesn't look like recent enough HTML Tidy"), nothing new.
  Examples (incl. --run-donttest), tests, and vignette rebuilds all OK.
  (A first run reported a third NOTE — "Files 'README.md' or 'NEWS.md'
  cannot be checked without 'pandoc'" — which was an environment
  artifact: the top-level-files check reads pandoc from PATH, not
  RSTUDIO_PANDOC; re-run with pandoc on PATH reproduced the baseline
  exactly.)
- **Gate 5 (targeted):** natural Boredom call runs (slow-tier test +
  spec-level 4-coding agreement); both sensitivity checks warning-free
  with the refit cross-check table (above); item-8 invariance exact.

## Open questions

1. **F-049 / the reweighting gate** — the 20-seed numbers are in; the
   choice between widening the gate, strengthening the noise estimator,
   and documenting the ~0.01-pip 2x-extrapolation bound is the
   maintainer's.
2. **Boredom in-sample calibration** (NEW-3): acceptable for 0.2.0.0, or
   worth a look before release?
3. Fixture group numbering flipped to first-appearance order (fr = 1) —
   intended per the brief, but flagged in case any downstream reading of
   the cached fixture fits assumed en = 1.
4. `prior_sensitivity_check()` on compare fits runs cold refits (compare
   is not omrf, so no warm start): 250 s on the Boredom fit. Fine for now;
   a warm-start port is a possible 0.2.1 item.

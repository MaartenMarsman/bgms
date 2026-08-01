# Brief 06 — bgmCompare defect + units batch (Opus agent)

You are the execution agent in the pre-CRAN review of the R package **bgms**.
You share no context with the review lead; everything you need is here. CODE
CHANGES authorized on a fix branch. This batch fixes defects found in the
maintainer's hands-on bgmCompare pass (report 07) plus two input-handling bugs
the lead verified, and executes a decided package-wide unit conversion. Item 1
is a suspected wrong-number defect — do it first and report its diagnosis
prominently.

## Setup

- Repo (Dropbox; do NOT switch its checked-out branch or build in it):
  `/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms`
- Base your branch on `develop` AT OR AFTER merge `b04dbd06` (contains both
  merged fix batches):
  ```sh
  git -C /Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms \
      worktree add ~/bgms-review/wt-fix3 -b fix/bgmcompare-batch develop
  cd ~/bgms-review/wt-fix3
  ```
- One commit per item, `fix:`/`feat:`/`docs:` prefixes, reference the F-number.
  No pushing; no attribution trailers. Every behavioral change gets a test that
  fails before, passes after.
- Context: reports 05/07 and FINDINGS.md are on this worktree under
  `dev/review-2026-08/`. The maintainer's console evidence is in report 07;
  his sensitivity plot is `reports/assets/sensitivity_compare.pdf`.
- Reproduction data ships with the package: `data("Boredom")` — column 1
  `language` (character, "en"/"fr", 496/490), columns 2–9 ordinal 1–7.
- Do NOT touch `tests/testthat/test-prior-sensitivity.R:432` (F-049, pending
  a maintainer decision; it is the slow tier's one expected failure).

## Items, in order

### 1. F-061 — sensitivity check on bgmCompare fits: recycled division (SUSPECTED WRONG NUMBERS)

Reproduction (maintainer, on develop ≥ `adf87013`):

```r
fit  = bgmCompare(Boredom[,-1], group_indicator = as.integer(factor(Boredom[,1])), seed = 1)
ps   = prior_sensitivity_check(fit)
# Warning messages:
# 1: In (p/(1 - p))/prior_odds :
#   longer object length is not a multiple of shorter object length
# 2: ...same...
```

Fires identically with `main_difference_selection = TRUE`. A recycling warning
inside a Bayes-factor computation means misaligned vectors, which means wrong
log-BFs — treat as a correctness defect until proven display-only.

Where the lead's trace points (verify, don't trust): `R/refit_engine.R:315-336`
builds `enm` from `raw$parameter_names$indicator %||% ...$pairwise`, builds
`pc`/`pbar` from the RB draw layout (`raw$rb_inclusion` column count), and then
sizes `prior_odds = rep(difference_prior_inclusion(fit), length(enm))` — names
and draws can disagree on a compare fit (36 indicator names incl. mains vs 28
pairwise RB columns on the Boredom fit, or whatever the true layouts are —
establish them). `R/prior_sensitivity.R:605` and `:630` divide by the same
`prior_odds`. Diagnose the actual layouts on both fits (default and
mains-selected), then fix so every vector is built off the SAME layout, keyed
to the fit's own parameter table — the F-037 fix (`indicator_pair_index()`
reading positions off the fit's layout) is the in-repo precedent for how to do
this right.

Verification, three parts:
1. Both Boredom fits run `prior_sensitivity_check()` with ZERO warnings.
2. Cross-check the curve against ground truth: brute-force refits at 2–3
   off-anchor `difference_scale` values; reweighted curve values must agree
   with the refit pips within the reweighting tolerance the single-group path
   already meets.
3. Re-run the maintainer's plot. His anomaly: `sit_around-half_dead_dull`
   shows evidence INCREASING with difference scale
   (`reports/assets/sensitivity_compare.pdf`). Report whether that survives
   the fix (it may be real; the point is to know), and save the after plot to
   `reports/assets/sensitivity_compare_fixed.pdf`.
4. State the blast radius in your report: were single-group (`bgm`) sensitivity
   results ever touched by this path, yes or no, with the code evidence.

### 2. F-057 — character/factor `group_indicator` crashes (major)

`bgmCompare(x = Boredom[,-1], group_indicator = Boredom$language)` dies with
`'bin' must be numeric or a factor`. Cause: `R/build_spec.R:439-442` — `group =
group_indicator` inherits character storage, the recode loop's integers coerce
back to character, `tabulate(group)` rejects. A factor indicator fails the same
way (`as.vector(factor)` is character). The repair `as.integer(group)` exists
only far downstream (`:652`). The post-listwise re-recode (`~:505-507`) repeats
the pattern — fix BOTH sites.

Fix: coerce once — `group = match(group_indicator, unique_g)` — which keeps the
existing first-appearance-order group numbering exactly (do not switch to
factor-level order; that would silently renumber groups for existing users) and
deletes the loop. Tests: character, factor, integer, 0/1-coded indicators all
fit and agree on group sizes; PLUS the natural shipped-data call
`bgmCompare(Boredom[,-1], group_indicator = Boredom$language)` runs — that
exact call is the regression test this defect earns. Also update the
`helper-fixtures.R` hand-coercions (`as.integer(as.factor(...))` at :216, :328,
:351 and the comment) to pass `language` directly — the workaround must not
outlive the defect, and switching the fixtures IS the coverage.
Check `y`-provided-as-vector too while in there (the x/y path).

### 3. F-056 — vector `x` gives an internal error instead of the validator's message (minor)

`bgmCompare(Boredom[,1], group_indicator = ...)` → `non-numeric matrix extent`,
because `R/bgmCompare.R:428-429` computes `num_variables = ncol(x)` and unpacks
the difference prior BEFORE `bgm_spec()` reaches `data_check(x, "x")`
(`R/bgm_spec.R:409`). `bgm()` with the same input answers
`x must be a matrix or data.frame.` — the message bgmCompare owes.

Fix: validate x before using its dimensions — run `x = data_check(x, "x")` at
the top of `bgmCompare()` (it is idempotent; `bgm_spec` runs it again
harmlessly) or hoist nothing and unpack inside the spec as `bgm()` does; your
call, smallest honest diff. Test: vector x → the validator's message, on both
`bgmCompare` and (as a guard) `bgm`.

### 4. F-060 — bgmCompare `verdicts()` print (major)

Report 07 shows the console verbatim. Three parts:

(a) On the DEFAULT fit (`main_difference_selection = FALSE`) the table
interleaves 8 main-effect rows as `pip NaN / log_bf NA / verdict <NA>`, counts
them in "36 indicators", and adds "8 indicator(s) were never updated and carry
no verdict." Those mains are not under selection — no indicator exists. Exclude
them from the table and the counts, and print one honest line instead, e.g.
"main-effect differences not under selection (main_difference_selection =
FALSE)". Keep genuine never-updated-but-selected indicators in the "never
updated" mechanism. With selection ON the current table is correct — do not
change it.

(b) The comparison vignette (`vignettes/comparison.Rmd:70-75`) promises a
scale-contingency caveat for difference verdicts (near-threshold difference
verdicts are scale-contingent; calibration under study). The print carries no
such line. Add it to compare-fit verdict prints — one or two sentences, matched
to the vignette's wording.

(c) Footer tone: "Run longer." → something like "Consider a longer run." Same
footer serves single-group fits; change it once, everywhere.

Snapshot tests for the default-fit print, the selection-on print, and the
footer.

### 5. F-021 — `calibration_check()` on a compare fit: graceful refusal (decided)

Today: `Error in UseMethod("calibration_check"): no applicable method ... class
"c('bgmCompare', 'S7_object')"`. The full method is DELIBERATELY deferred
(maintainer decision; a prepared patch exists and was not applied — do NOT
implement the method). Add a stub that fails informatively: calibration checks
are not yet available for bgmCompare fits, planned for a future release. Match
however the package words its other not-supported errors. Test asserts the
message, not the raw dispatch error.

### 6. F-063 — difference centrality: remove from the plot surface (decided)

Maintainer on `plot(fit, type = "centrality")` for compare fits: "No idea what
these mean. Lets not do centrality for differences." Establish what those
panels currently draw (report it — one paragraph), then remove the
difference-scale centrality display for compare fits. If per-group centrality
panels exist and are coherent, keep them; requesting a difference-centrality
plot should say why it is not offered rather than draw one. Update the Rd for
`plot.bgmCompare`/`extract_centrality` accordingly (extractor itself: leave
per-group; if a difference-centrality extractor path exists, leave it but
document the interpretation caveat — the DECIDED removal is the plot).

### 7. F-062 — compare-fit plots (a–c decided, d if cheap)

(a) DECIDED: the difference network must encode evidence the way `bgm`'s
network plot does (verdict/BF-encoded edges). First establish what it encodes
today (posterior-mean differences? weights?) and say so in the report; then
converge on the bgm convention so one visual language covers both functions.
(b) The legend overlaps the plot at default size — fix layout.
(c) Reword: the "main-effect differences not selected" banner and the
"node: main-effect difference" legend label (maintainer: both read wrong).
Align wording with item 4(a)'s line.
(d) Maintainer proposal, implement IF it stays small with qgraph's `pie`
argument (qgraph is already the Suggests renderer): when main selection is ON,
draw each node's main-difference evidence as a pie ring keyed to pip (1 = full,
0.5 = half). If it is not a small change, deliver a feasibility note in the
report instead — do not sink the batch into it.
Applies to the default plot and `type = "groups"`.

### 8. Natural log everywhere (DECIDED; closes the F-039 residue, refreshes F-051)

Maintainer, report 07: the verdicts-in-nats / sensitivity-in-log10 mix is
jarring — "Lets make sure we have only natural log everywhere." Convert
`prior_sensitivity_check()` and its print/plot/summary surface from log10 to
natural log. Report 03 §5 has the complete site inventory. Mechanics:

- Column rename `chosen_scale_log10_bf` → `chosen_scale_log_bf` (values in
  nats); grep every reader, including the F-051 docs block in
  `R/prior_sensitivity.R` (`?plot.bgms_prior_sensitivity` Details) merged in
  brief 05 — rewrite it in nats (cap becomes ln(1e6) ≈ 13.82).
- Tuned constants re-express EXACTLY, not renamed: tolerance 0.5 log10 =
  0.5·ln(10) ≈ 1.1513 nats; the |log10 BF| ≤ 3 window → 3·ln(10) ≈ 6.9078
  nats; the plot cap 6 log10 → 13.8155 nats. Same decisions, new unit.
- Plot axis, zone labels, print headers, Rd text; boundaries print as
  ±log(threshold) like `verdicts()` does (±2.30 at 10).
- Deprecation courtesy: if the old column name is part of the documented
  return, keep it one release as a duplicate with a deprecation note in the
  Rd — your judgment given how the package handled `log10_bf` → `log_bf` in
  brief 03 (look at that precedent, `20b56996`).

INVARIANCE CHECK, mandatory: on one fixed single-group fit and one compare fit
(after item 1's fix), every verdict, every scale-dependence classification, and
the curve SHAPES are identical before/after — only labels and units move.
Report the check.

### 9. F-065 — the 5–10 s silent startup (diagnose; fix only if small)

`bgmCompare(...)` at Boredom defaults sits silent 5–10 s before the progress
bar. Find what occupies it (timestamps around the spec build / validation /
warmup init are enough). If it is one obvious thing, either narrate it (a
verbose-gated message like the zratio fence message) or shorten it; if it is
structural, report the breakdown and stop — no redesign in this batch.

## Verification gate

From a clean install of the branch:
1. Changed-file tests pass; new regression tests fail on `b04dbd06` (spot-check
   one per item, state which).
2. Full slow tier: exactly ONE failure (`test-prior-sensitivity.R:432`, F-049),
   ZERO warnings.
3. Full CRAN-mode suite: 0 failures, 0 warnings.
4. `R CMD check --as-cran` on a `git archive` tarball: the same 2 baseline
   NOTEs (stale Date, old HTML Tidy), nothing new.
5. Targeted: the natural Boredom call runs; both sensitivity checks
   warning-free with the refit cross-check table; the invariance check from
   item 8.

## Deliverable

Write
`/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms/dev/review-2026-08/reports/06-bgmcompare-defect-batch.md`
(write it even if the checkout is not on your branch — the lead collects it):
What was done per item (commit SHAs) / Findings (new ones severity-tagged;
item 1's diagnosis and blast radius FIRST) / Evidence (gate outputs verbatim,
before/after console excerpts, the two sensitivity plots) / Open questions.

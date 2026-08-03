# Report 21 — mark prior-only rows in bgmCompare summaries (F-117)

Brief 21 + the maintainer rulings on the first version of this report. Code
changes authorized. Branch `fix/summary-zero-marking`, based on
`origin/develop` = `afe11510` (contains `acc19428`), 4 commits — **not pushed**.
Worktree `~/bgms-review/wt-fix12`; the Dropbox tree's branch was never switched
and nothing was built in it. Files touched: `R/methods_bgmcompare.R`,
`R/bgmCompare.R`, `R/validate_data.R`, `man/summary.bgmCompare.Rd`,
`man/bgmCompare.Rd`, `tests/testthat/test-collapse-categories.R`, the new
`tests/testthat/test-summary-marking.R` + its snapshot, and this report. No
off-limits file was touched (`R/validate_data.R` was released by the ruling).

> **Revision note.** §§1–4 describe the shipped behaviour, which is the
> **widened** rule the maintainer ruled for: a variable's rows are all marked
> when any group's *reference*-category cell is empty, in addition to each
> row's own cell, with the footnote text the ruling specifies. §5 records what
> the narrow rule looked like and what the ruling changed. The F-121
> diagnosis is a separate addendum at the end.

**Headline.** The mark ships: a leading `*` on the affected rows plus one
footnote line, printed only when a rendered row carries the mark. Rows are
mapped to support cells by index and both difference-row layouts are handled.
The fit-time warning is widened to match. Gate is clean: 0 failures, 0
warnings, no existing snapshot changed, `document()` idempotent, NAMESPACE
unchanged.

**One thing still wants the maintainer's eye: §5.2 — the footnote and the
warning both speak of "a group", but the rows are contrasts, not groups.** A
row is `A (diff1; 3)`, an orthogonal contrast over all groups. The ruled
wording ("a group lacks observations…") is a real improvement on the original
("no observations in *this* group…") because it no longer promises that the
row names one — but it is worth knowing that no row ever will.

**F-121 (the build-to-build divergence reported in §3.4) is diagnosed in the
addendum. Verdict: legal floating-point compilation variance, not undefined
behaviour — and narrower than it looked.** One of the two builds was
unoptimized; two properly optimized builds at different paths agree to the
last digit. Sanitizers and the object bisection are in the addendum.

---

## 1. What was done

| # | Item | Commit | Type |
|---|---|---|---|
| 1 | Mark + footnote + Rd | `20454889` | `fix(compare):` |
| 2 | New test file + snapshot | `90cee65d` | `fix(compare):` |
| 3 | Widen to the reference category; widen the fit-time warning; Blume–Capel test | `97605121` | `fix(compare):` |
| 4 | Make the snapshot invariant to the build (§2) | `4ebabd24` | `fix(compare):` |

```
 R/bgmCompare.R                            |  10 +-
 R/methods_bgmcompare.R                    | 132 +++++++++++++++
 R/validate_data.R                         |  39 ++++-
 man/bgmCompare.Rd                         |  10 +-
 man/summary.bgmCompare.Rd                 |  23 +++
 tests/testthat/_snaps/summary-marking.md  |  63 +++++++
 tests/testthat/test-collapse-categories.R |  26 +++
 tests/testthat/test-summary-marking.R     | 271 ++++++++++++++++++++++++++++++
```

### 1.1 Where the mark goes, and why

**A leading `*` in front of the parameter label, not a trailing marker.** Two
reasons. The mark qualifies *which parameter the row is* — it belongs next to
the row's identity, not after its numbers. And the difference table's last
column (`Rhat`, and `share_incl` before it) is routinely blank, so a trailing
`*` would float in whitespace next to `mean` on the *next* line as often as it
would read as a column of its own. `print.data.frame` right-aligns the label
column, so prefixing `"* "` to the marked rows and `"  "` to the rest produces
a clean two-character gutter at no cost to the existing alignment.

Only the printed table changes. `summary(fit)$main_diff` is untouched — the
labels there stay exactly as they were, which the test asserts.

### 1.2 How rows are mapped to support cells

By index, never by parsing the label. `compare_prior_only_main_diff()`
(`R/methods_bgmcompare.R:74`) rebuilds the main-effects row layout from
`num_categories` and `is_ordinal_variable` — `num_categories[v]` threshold rows
for an ordinal variable, two (`linear`, `quadratic`) for a Blume–Capel one —
and threshold `k` of variable `v` reads `category_support[[v]][k + 1, ]`
(row 1 is the reference category 0).

Four things the mapping has to get right:

- **Two ways a row can be prior-driven.** Its own cell can be empty: that group
  observes no one in that category. Or the *reference* cell can be empty: every
  threshold is identified relative to category 0, so a group that never used
  category 0 has no data fixing the level of its threshold vector, and *all* of
  that variable's thresholds go with it. The second condition marks the whole
  variable. In code that is one clause —
  `by_row[[v]] = own | any(cells[1L, ] == 0L)`.
- **Two row layouts.** The two summarizers order their rows differently. With
  selection, `summarize_main_diff_compare()` walks variable → threshold →
  contrast, so the contrasts of one threshold are adjacent. Without selection,
  `summarize_manual_compare()` runs over the raw difference columns, which are
  contrast-major: the whole main-effects matrix once per contrast. Both are
  handled, keyed on `arguments$difference_selection`. Two groups cannot tell
  the two orderings apart (one contrast), so the test uses three — see §2.
- **Rows are contrasts, not groups.** A group's effect is
  `baseline + projection[g, ] %*% differences`, and the projection is a dense
  eigen contrast basis, so an empty cell in *any* group reaches *every*
  contrast of that variable-by-category pair. A pair is therefore marked in all
  of its contrasts or in none.
- **Blume–Capel variables are exempt** from the union recode and carry no
  support matrix, so their rows are never marked.

The function returns `NULL` — meaning "print exactly as before" — when
`category_support` is absent, when `num_groups < 2`, and on any layout it
cannot verify (wrong list length, wrong matrix shape, row count that does not
match the summary). That is the graceful-degradation path for old fits, and it
is a positive check rather than an error handler.

### 1.3 The fit-time warning, widened to match

`collapse_categories_across_groups()` listed every empty cell the same way, so
`variable 'A', category 0, group 2` read as one cell among others when it is in
fact the whole threshold vector of `A` for group 2. Three changes:

- **Reference cells say what they are.** The line becomes
  `variable 'A', category 0, group 2 -- the reference category; every threshold
  of this variable is affected for that group`.
- **Reference cells are listed first.** The list is capped at ten with a
  `... and N more`; emitting reference cells ahead of ordinary ones means the
  cap can no longer bury the consequential ones behind routine ones.
- **The prose gained the reference case** — "An empty reference category is
  worse: every threshold is measured relative to category 0, so all of that
  variable's threshold differences for that group rest on the prior, not just
  one."

**One thing this turned up.** The widened message overran
`getOption("warning.length")` (1000 by default), and R truncates there — so the
tail, which is the `see ?summary.bgmCompare` pointer, was being cut. The prose
is tightened to land the two-cell case at 906 characters, and
`test-collapse-categories.R` now pins `nchar(w) < getOption("warning.length")`.
Long cell lists can still truncate, as they could before this change; the
reference-first ordering is what makes that survivable.

### 1.4 Documentation

`man/summary.bgmCompare.Rd` gains a `\details` block naming the mark, giving
the footnote text, distinguishing the two conditions, and pointing at
`extract_arguments(fit)$category_support`. `man/bgmCompare.Rd`'s
union-semantics paragraph gains the cross-reference sentence and one sentence
on the reference category — needed because that paragraph is where the union
semantics are explained, and it otherwise described only the narrow case that
the ruling widened. That is one sentence more than brief 21 allowed in that
file; revert it if you would rather keep it minimal.

---

## 2. Tests — `tests/testthat/test-summary-marking.R`

A new file, per the brief, to stay clear of the files the other batch edits.
Five blocks, six fits, **~5 s of sampling in total** (0.35 s for each two-group
fit, 1.7 s for each three-group one; `n = 160`, `180` or `210`, three
variables, `iter = warmup = 200`, `chains = cores = 1`). The widened-warning
assertions live in `test-collapse-categories.R`, which already owns that
warning's tests.

| Block | What it fixes |
|---|---|
| (a) adverse fit | Two groups, four-level variable A on {0,1,2} vs {1,2,3}. Asserts the fixture is what it claims (`support[[1]][4,1] == 0`, `support[[1]][1,2] == 0`, B and C complete), that the rendered mark pattern is exactly `c(T,T,T,F,F,F)` and the three marked rows are A's thresholds 1–3, that the footnote line is present verbatim, that `summary(fit)$main_diff$parameter` carries no `*`, and a snapshot of the whole printed summary. |
| (b) layout | Three groups, run twice with `difference_selection` TRUE and FALSE. Checks `compare_prior_only_main_diff()` against `grepl("^A \\(diff\\d+; \\d+\\)$", labels)` — the labels the marking code never reads — and then that the rendered gutter agrees with it. **Two groups cannot catch a swapped layout; this is the block that can.** |
| (c) shared support | Both groups use every category. No mark, no footnote, and the labels print with no gutter at all. |
| (d) **Blume–Capel** | A mixed fit: ordinal `A` and Blume–Capel `D` drawn from *the same* split-support data, plus a complete ordinal `B`. Asserts `is_ordinal_variable[2]` is `FALSE`, `category_support[[2]]` is `NULL`, that `D` contributes a `linear` and a `quadratic` row per contrast, and that the flag vector is exactly `grepl("^A \\(", labels)` — A marked throughout, D never, on identical data. Pins the exemption. |
| (e) old fit | The adverse fit with `category_support` stripped from its arguments. No error, no warning, no mark, no footnote, labels unpadded. |

The classed condition is **named**, not muffled wholesale:
`expect_warning(fit <- ..., class = "bgms_group_support_warning")`, so an
unrelated warning from the fit still fails the block.

**The snapshot drops the MCMC numbers.** The `transform` strips the trailing
run of numeric cells from each table row and then collapses runs of spaces,
keeping the row name, the mark gutter, the label and the prose. What is under
test is which rows carry the mark and the exact footnote text; the F-121
addendum is the reason posterior means have no business in a snapshot.

**The first version of this transform was wrong, and F-121 caught it.** Dropping
the numbers is not sufficient: `print.data.frame` sizes each column to its
widest value, so the *header* line carries the width of the widest `n_eff` in
the run — `n_eff` is padded by `262.605` and not by `46.143`. The snapshot
recorded against one build and failed against another with
`- parameter mean mcse sd  n_eff …` / `+ parameter mean mcse sd n_eff …`. That
is one space, and it is MCMC output by proxy. Fixed in `4ebabd24` by collapsing
whitespace, and verified the only way that means anything: the same fixture
rendered under two builds that genuinely disagree on the numbers (`-O2` giving
2.865 and `-O2 -ffp-contract=off` giving 3.035) produces **byte-identical**
transformed output.

---

## 3. Evidence

### 3.1 Full local default-tier suite — 0 failures, 0 warnings

```
$ Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_local(reporter="summary")'
...
collapse-categories: ..........................................................................................S
summary-marking: ......................................
verdicts: ..................................................................................................
...
══ Skipped ═════════════════════════════════════════════════════════════════════
(102 skips)
══ DONE ════════════════════════════════════════════════════════════════════════
exit=0
```

Roughly 10 minutes wall. No `Failed` section, no `Warnings` section, working
tree clean afterwards. All 102 skips are tier-gated (`BGMS_RUN_SLOW_TESTS` /
`BGMS_RUN_CERTIFICATION`) or the pre-existing `golden fixtures not found` in
`test-simulate-predict-regression.R` (those fixtures are generated by a
separate script and are not in the repo) — none of them mine and none of them
new.

This is the fourth full run. The third failed, on the snapshot, for the reason
in §2 — worth stating plainly rather than only reporting the green run: the
snapshot as first recorded was build-dependent through one space of column
padding, and it took a run against a differently-compiled `bgms.so` to expose
it. Fixed and re-verified before this run.

### 3.2 No existing snapshot changed

`_snaps/plot-methods.md` and `_snaps/verdicts.md` are untouched and neither
produced a `.new.md`; the only snapshot in the diff is the new
`_snaps/summary-marking.md`. That is the expected result: both existing
snapshots build their input from synthetic deterministic objects, and every
existing compare fit in the suite has shared support, so no fit in the suite
reaches the marking path.

### 3.3 `document()` clean, NAMESPACE unchanged

```
$ Rscript -e 'suppressMessages(devtools::document())'
$ git status --porcelain
 M R/bgmCompare.R
 M R/methods_bgmcompare.R
 M R/validate_data.R
 M man/bgmCompare.Rd
 M man/summary.bgmCompare.Rd
 M tests/testthat/_snaps/summary-marking.md
 M tests/testthat/test-collapse-categories.R
 M tests/testthat/test-summary-marking.R
```

Idempotent (second run produces no further change), `NAMESPACE` absent from the
list, no new `man/` file.

### 3.4 The build-to-build divergence — F-121

Two builds of the identical C++ source gave different posterior means from the
same seed. That observation stands, and it is why the BEFORE/AFTER in §4 is
rendered from **one** build with the R code swapped rather than from two.
Establishing that the divergence was not mine took three runs:

| Build | R code | `A (diff1; 1)` mean |
|---|---|---|
| `wt-fix12-base` @ `afe11510` | base | 2.865 |
| `wt-fix12` | **base** (checked out over mine) | 2.431 |
| `wt-fix12` | mine | 2.431 |

Base R code on my build reproduces my numbers exactly, so the change is
print-only, as designed. Copying the base `bgms.so` into the third library and
changing nothing else flips it straight back to 2.865.

**The diagnosis is in the addendum**, and it narrows the finding considerably:
one of the two builds was compiled *unoptimized*, and two properly optimized
builds at different paths agree to the last digit.

---

## 4. BEFORE and AFTER — verbatim console output

Both rendered by the same script, same seed (117), same data, same build; the
only difference is the R code. The two-group adverse fixture is the one the
test file uses: variable A four-level, group 1 on {0,1,2}, group 2 on {1,2,3};
B and C on {0,1,2} in both.

```
### category_support[[1]] (variable A)

           group 1 group 2
category 0      26       0
category 1      23      27
category 2      31      26
category 3       0      27
```

### BEFORE — `summary(fit)`

```
Posterior summaries from Bayesian grouped MRF estimation (bgmCompare):

groups: 1 = 1 (n = 80), 2 = 2 (n = 80)

Category thresholds:
  parameter   mean  mcse    sd   n_eff  Rhat
1     A (1)  1.381 0.039 0.462 138.151 0.996
2     A (2)  1.439 0.046 0.496 118.300 1.005
3     A (3) -1.431 0.189 1.285  46.090 1.083
4     B (1) -0.081 0.027 0.242  80.302 1.039
5     B (2) -0.309 0.042 0.378  80.619 1.013
6     C (1)  0.014 0.026 0.247  91.185 1.016
... (use `summary(fit)$main` to see full output)

Pairwise interactions:
  parameter   mean  mcse    sd   n_eff  Rhat
1       A-B  0.037 0.004 0.043 107.518 1.029
2       A-C -0.010 0.004 0.045 105.482 1.032
3       B-C -0.058 0.004 0.044 135.890 0.995

Inclusion probabilities:
      parameter  mean  mcse    sd   n_eff  Rhat n0->1 n1->0
       A (main)                                     0     0
 A-B (pairwise) 0.359 0.057 0.414  52.538 0.996    15    15
 A-C (pairwise) 0.713 0.097 0.415  18.289  1.03     7     7
       B (main)                                     0     0
 B-C (pairwise) 0.086 0.018 0.197 124.354 1.028    15    16
       C (main)                                     0     0
Note: NA values are suppressed in the print table; they occur for indicators
that were not updated or whose draws are constant, so ESS/Rhat are undefined.
`summary(fit)$indicator` still contains all computed values.

Group differences (main effects):
    parameter   mean  mcse    sd   n_eff share_incl Rhat
 A (diff1; 1)  2.431 0.101 1.043 106.442          0     
 A (diff1; 2)  1.738 0.134 1.093  66.864          0     
 A (diff1; 3)  7.067 0.435 2.953  46.143          0     
 B (diff1; 1) -0.106 0.035 0.393 125.242          0     
 B (diff1; 2)  0.019 0.062 0.534  73.770          0     
 C (diff1; 1) -0.130 0.074 0.526  50.065          0     
... (use `summary(fit)$main_diff` to see full output)
Note: NA values are suppressed in the print table. They occur for differences
that were never selected, so the composite ESS and share are undefined;
`summary(fit)$main_diff` still contains the NA values.

Group differences (pairwise effects):
   parameter  mean  mcse    sd   n_eff share_incl Rhat
 A-B (diff1) 0.054 0.010 0.090  84.359      0.923     
 A-C (diff1) 0.153 0.022 0.119  29.872      0.893     
 B-C (diff1) 0.004 0.002 0.026 262.605      0.203     
Note: NA values are suppressed in the print table. They occur for differences
that were never selected, so the composite ESS and share are undefined;
`summary(fit)$pairwise_diff` still contains the NA values.

Use `summary(fit)$<component>` to access full results.
See the `easybgm` package for other summary and plotting tools.
```

### AFTER — `summary(fit)`

```
Posterior summaries from Bayesian grouped MRF estimation (bgmCompare):

groups: 1 = 1 (n = 80), 2 = 2 (n = 80)

Category thresholds:
  parameter   mean  mcse    sd   n_eff  Rhat
1     A (1)  1.381 0.039 0.462 138.151 0.996
2     A (2)  1.439 0.046 0.496 118.300 1.005
3     A (3) -1.431 0.189 1.285  46.090 1.083
4     B (1) -0.081 0.027 0.242  80.302 1.039
5     B (2) -0.309 0.042 0.378  80.619 1.013
6     C (1)  0.014 0.026 0.247  91.185 1.016
... (use `summary(fit)$main` to see full output)

Pairwise interactions:
  parameter   mean  mcse    sd   n_eff  Rhat
1       A-B  0.037 0.004 0.043 107.518 1.029
2       A-C -0.010 0.004 0.045 105.482 1.032
3       B-C -0.058 0.004 0.044 135.890 0.995

Inclusion probabilities:
      parameter  mean  mcse    sd   n_eff  Rhat n0->1 n1->0
       A (main)                                     0     0
 A-B (pairwise) 0.359 0.057 0.414  52.538 0.996    15    15
 A-C (pairwise) 0.713 0.097 0.415  18.289  1.03     7     7
       B (main)                                     0     0
 B-C (pairwise) 0.086 0.018 0.197 124.354 1.028    15    16
       C (main)                                     0     0
Note: NA values are suppressed in the print table; they occur for indicators
that were not updated or whose draws are constant, so ESS/Rhat are undefined.
`summary(fit)$indicator` still contains all computed values.

Group differences (main effects):
      parameter   mean  mcse    sd   n_eff share_incl Rhat
 * A (diff1; 1)  2.431 0.101 1.043 106.442          0     
 * A (diff1; 2)  1.738 0.134 1.093  66.864          0     
 * A (diff1; 3)  7.067 0.435 2.953  46.143          0     
   B (diff1; 1) -0.106 0.035 0.393 125.242          0     
   B (diff1; 2)  0.019 0.062 0.534  73.770          0     
   C (diff1; 1) -0.130 0.074 0.526  50.065          0     
... (use `summary(fit)$main_diff` to see full output)
* a group lacks observations in this category or in the reference category; the estimate reflects the prior, not the data
Note: NA values are suppressed in the print table. They occur for differences
that were never selected, so the composite ESS and share are undefined;
`summary(fit)$main_diff` still contains the NA values.

Group differences (pairwise effects):
   parameter  mean  mcse    sd   n_eff share_incl Rhat
 A-B (diff1) 0.054 0.010 0.090  84.359      0.923     
 A-C (diff1) 0.153 0.022 0.119  29.872      0.893     
 B-C (diff1) 0.004 0.002 0.026 262.605      0.203     
Note: NA values are suppressed in the print table. They occur for differences
that were never selected, so the composite ESS and share are undefined;
`summary(fit)$pairwise_diff` still contains the NA values.

Use `summary(fit)$<component>` to access full results.
See the `easybgm` package for other summary and plotting tools.
```

All three of A's thresholds carry the mark. Threshold 3 has its own empty cell
(group 1 never uses category 3); thresholds 1 and 2 are marked because group 2
never uses the reference category. Under the narrow rule only the 7.067 row was
marked and the 2.431 and 1.738 rows went out clean — those two are what the
ruling recovers.

### The fit-time warning on the same fit

```
Warning message:
Some categories were not used by every group:
  variable 'A', category 0, group 2 -- the reference category; every threshold of this variable is affected for that group
  variable 'A', category 3, group 1
These categories are kept, because the other groups do use them. But a group with no observations in a category has nothing to say about where its threshold for that category lies, so the reported difference for that group and that category is set by the prior, not by the data. An empty reference category is worse: every threshold is measured relative to category 0, so all of that variable's threshold differences for that group rest on the prior, not just one. Expect large, very uncertain numbers, and do not read them as evidence of a group difference. Only the category thresholds are affected, not the pairwise (edge) differences. The printed summary marks the rows; see ?summary.bgmCompare.
```

906 characters — inside R's 1000-character warning cap, so the pointer at the
end survives. Before the prose was tightened it did not (§1.3).

### The whole diff

`diff BEFORE AFTER` is three things and nothing else — every number is
byte-identical:

```
44,50c44,50
<     parameter   mean  mcse    sd   n_eff share_incl Rhat
<  A (diff1; 1)  2.431 0.101 1.043 106.442          0     
<  A (diff1; 2)  1.738 0.134 1.093  66.864          0     
<  A (diff1; 3)  7.067 0.435 2.953  46.143          0     
<  B (diff1; 1) -0.106 0.035 0.393 125.242          0     
<  B (diff1; 2)  0.019 0.062 0.534  73.770          0     
<  C (diff1; 1) -0.130 0.074 0.526  50.065          0     
---
>       parameter   mean  mcse    sd   n_eff share_incl Rhat
>  * A (diff1; 1)  2.431 0.101 1.043 106.442          0     
>  * A (diff1; 2)  1.738 0.134 1.093  66.864          0     
>  * A (diff1; 3)  7.067 0.435 2.953  46.143          0     
>    B (diff1; 1) -0.106 0.035 0.393 125.242          0     
>    B (diff1; 2)  0.019 0.062 0.534  73.770          0     
>    C (diff1; 1) -0.130 0.074 0.526  50.065          0     
51a52
> * a group lacks observations in this category or in the reference category; the estimate reflects the prior, not the data
```

### Both row layouts, three groups

Group 1 on {0,1,2}, group 2 on {1,2,3}, group 3 on {0,1,2,3}; two contrasts.
This is the check that the mark tracks the *layout* and not the label order.

Group 2 never uses A's reference category, so under the ruled rule every row of
A is marked in both layouts — which makes the *ordering* the thing to read
here: which six rows land above the fold differs between the two summarizers,
and the gutter tracks that, not the labels.

**`difference_selection = TRUE`** — variable-major, contrasts adjacent, so all
six visible rows belong to A:

```
Group differences (main effects):
      parameter   mean  mcse    sd   n_eff share_incl Rhat
 * A (diff1; 1) -1.667 0.053 0.568 113.106          0     
 * A (diff2; 1) -1.111 0.029 0.369 159.675          0     
 * A (diff1; 2) -1.382 0.059 0.583  96.305          0     
 * A (diff2; 2) -0.842 0.029 0.366 162.717          0     
 * A (diff1; 3) -1.388 0.058 0.590 104.848          0     
 * A (diff2; 3) -6.397 0.313 2.799  79.969          0     
... (use `summary(fit)$main_diff` to see full output)
* a group lacks observations in this category or in the reference category; the estimate reflects the prior, not the data
```

**`difference_selection = FALSE`** — contrast-major, the whole matrix once per
contrast, so A's three rows are followed by B and C, and contrast 2 is below
the fold:

```
Group differences (main effects):
      parameter   mean  mcse    sd   n_eff  Rhat
 * A (diff1; 1) -1.574 0.094 0.651  47.721 1.030
 * A (diff1; 2) -1.143 0.090 0.690  58.617 1.028
 * A (diff1; 3) -0.870 0.093 0.744  64.141 1.020
   B (diff1; 1) -0.265 0.024 0.341 200.000 0.995
   B (diff1; 2) -0.183 0.040 0.466 138.209 1.001
   C (diff1; 1) -0.251 0.023 0.331 200.000 1.006
... (use `summary(fit)$main_diff` to see full output)
* a group lacks observations in this category or in the reference category; the estimate reflects the prior, not the data
```

---

## 5. Open questions

### 5.1 The reference category — RULED AND SHIPPED

The first version of this report shipped the narrow rule from F-117 (mark a row
only when its *own* cell is empty) and flagged that it left the worse case
unmarked. In the adverse fixture group 2 never uses category 0, so all three of
A's threshold differences are prior-driven, but only threshold 3 carried a
mark:

| row | mean | sd | narrow rule | ruled rule |
|---|---|---|---|---|
| `A (diff1; 1)` | **2.431** | 1.043 | unmarked | `*` |
| `A (diff1; 2)` | **1.738** | 1.093 | unmarked | `*` |
| `A (diff1; 3)` | **7.067** | 2.953 | `*` | `*` |
| `B (diff1; 1)` | −0.106 | 0.393 | — | — |
| `C (diff1; 1)` | −0.130 | 0.526 | — | — |

The maintainer ruled to widen. Shipped in `97605121`: one clause in
`compare_prior_only_main_diff()` (`by_row[[v]] = own | any(cells[1L, ] == 0L)`),
the ruled footnote text, and the same widening in the fit-time warning (§1.3).
The rendering is in §4. **Closed.**

### 5.2 The rows are contrasts, and no wording will make one name a group

The ruled footnote — "*a group* lacks observations in this category or in the
reference category" — is a real improvement on the original "no observations in
*this* group", because it no longer promises that the row identifies one. It is
worth being explicit that none ever will: `A (diff1; 3)` is an orthogonal
contrast whose group loadings here are `projection[, 1] = (-0.5, 0.5)`, so
every group appears in every difference row. The same is true of the fit-time
warning, which *does* name the group (`group 2`) — correctly, because it is
talking about the data, not about a row.

So the two messages are consistent and both true; they just refer to different
objects, and a reader moving from the warning to the table cannot map "group 2"
onto any particular marked row. If you want that mapping to exist, the honest
route is not wording but data: expose the offending groups per marked row, e.g.
a `prior_only` attribute on the summary carrying the group indices. **Nothing
shipped for this — flagging it as the residual imprecision, not proposing more
text.**

### 5.3 Smaller things

- **The footnote prints between the truncation line and the NA note.** With
  both present the block ends with three consecutive advisory lines. It reads
  fine (see §4) but if the ordering matters to you, say which.
- **Only `main_diff` is marked.** Correct per the warning text — pairwise
  differences use all the data — and `pairwise_diff` is untouched.
- **`summary(fit)$main_diff` deliberately carries no mark.** The marking is a
  display device; the data frame is a downstream contract (easybgm, JASP). If
  you would rather expose it programmatically, a `prior_only` logical on the
  summary object is the cheap version — say so and it lands.
- **Blume–Capel variables are never marked — now pinned.** They are exempt from
  the union recode and carry no support matrix, so the `NULL` entry in
  `category_support[[v]]` maps to two unmarked rows. Test block (d) builds an
  ordinal and a Blume–Capel variable from *the same* split-support data and
  asserts the ordinal one is marked throughout while the Blume–Capel one never
  is. **Closed.**
- **The warning can still truncate on a long cell list.** `head(…, 10L)` plus
  the prose can exceed 1000 characters when many cells are empty — true before
  this change too. What is new is that reference cells are emitted first, so
  the cells that matter most survive the cut. Raising
  `getOption("warning.length")` is the caller's decision, not ours.

---

# Addendum — F-121 diagnosis

**Verdict: legal floating-point compilation variance, not undefined behaviour.**
No code changes made, per the brief.

And the finding is narrower than it looked. **Two independent optimized builds
of the same commit, at different filesystem paths, agree to every reported
digit.** bgmCompare *is* reproducible across builds at the package's own
settings. The divergence in §3.4 came from one build that was compiled
**unoptimized**, which is not a configuration the package ever asks for and
which I could not reproduce afterwards from the same command in the same
directory.

## A.1 What the two builds actually were

Object-file bisection (§A.3) landed on `src/mcmc/algorithms/hamiltonian_utils.o`
and, once I looked at it, the two objects were not subtly different — they were
different by 3×:

| | anomalous build (`wt-fix12`, 02:09) | good build (`wt-fix12-base`, 02:41) |
|---|---|---|
| `__TEXT,__text` | 33 488 B | 11 548 B |
| `__LD,__compact_unwind` | 8 032 B | 1 504 B |
| `__DWARF,__debug_loclists` | absent | 10 260 B |
| `__TEXT,__literal8/16` | absent | present |
| symbols | 395 | 181 |

Same source (`shasum` identical), same compiler
(`Apple clang version 21.0.0 (clang-2100.1.1.101)` in both `DW_AT_producer`),
byte-identical `src/Makevars`, no `~/.R/Makevars` on this machine. Every one of
the 44 objects in the anomalous build is 2–4× its counterpart.

The symbols only the big object carries say what happened: `arma::arma_check<…>`,
`rnorm(SafeRNG&, double, double)`, `Rcpp::Rstreambuf<…>::Rstreambuf()`,
`LeapfrogJointResult::~LeapfrogJointResult()` — helpers that an optimizing
compiler inlines away and an unoptimizing one emits out of line, plus 18 extra
`GCC_except_table` entries.

**Reproduced exactly.** Recompiling the same file by hand:

```
-O2 -DNDEBUG -DARMA_NO_DEBUG   -> text=11548  sha=348254cba7c5faa7   (= good build, modulo comp_dir)
-O0 -DNDEBUG -DARMA_NO_DEBUG   -> text=33360  sha=7c5493d6a24cb994
-O0          -DARMA_NO_DEBUG   -> text=33488  sha=291d99639266c383   <-- byte-identical to the anomalous object
```

So the anomalous compile ran at **`-O0` with `-DNDEBUG` absent**. R's own
`Makeconf` hardcodes both (`CXX20FLAGS = … -O2`, `-DNDEBUG` in `ALL_CPPFLAGS`),
so something overrode them for that invocation.

**I could not reproduce it.** Deleting the object and re-running
`R CMD INSTALL` in the same directory now produces the `-O2` object
byte-for-byte, and echoes `-Wall -g -O2` — which is also what the original
build's log echoed. I cannot reconcile the echoed flags with the object; the
overriding mechanism (a transient `~/.R/Makevars` or `R_MAKEVARS_USER` from one
of the other agents on this machine is the obvious candidate — none is present
now) is **unidentified**. That is the one loose end here.

## A.2 Sanitizers

**UBSan: clean.** Built from `afe11510` at `-O1 -fno-omit-frame-pointer
-fsanitize=undefined -fno-sanitize=vptr` (`R_MAKEVARS_USER` overriding
`CXX20FLAGS` / `SHLIB_CXX20LDFLAGS`; the package is `CXX_STD = CXX20`, so the
plain `CXXFLAGS` names are ignored — worth knowing for the next agent).

```
$ Rscript --vanilla probe.R ~/bgms-review/f121/lib-ubsan     # the exact seeded fit from §3.4
2.431
```

No diagnostics — no output at all beyond the value. Then the compare test
files:

```
$ UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=0:report_error_type=1 \
  R_LIBS=…/lib-ubsan Rscript --vanilla -e 'testthat::test_dir("tests/testthat", package="bgms",
      filter="bgmCompare|collapse-categories|methods|mcmc-diagnostics",
      load_package="installed", reporter="summary")'
bgmCompare-gradient: ..............
bgmCompare: ............................S..............S......................
collapse-categories: ....................................................................................S
mcmc-diagnostics: ......................................................S...........................…
methods: ………(all pass)
plot-methods: ......SSS.....SS.........................SSSS...S.S..S.....S......SSSSSSS....S
```

`grep -c "runtime error"` over the whole log: **0**. All tests pass.

The instrumentation is verifiably present, so the zero means something:

```
$ nm -u lib-ubsan/bgms/libs/bgms.so | grep -c __ubsan_handle
13
  ___ubsan_handle_add_overflow          ___ubsan_handle_out_of_bounds
  ___ubsan_handle_alignment_assumption  ___ubsan_handle_pointer_overflow
  ___ubsan_handle_builtin_unreachable   ___ubsan_handle_shift_out_of_bounds
  ___ubsan_handle_divrem_overflow       ___ubsan_handle_sub_overflow
  ___ubsan_handle_float_cast_overflow   ___ubsan_handle_load_invalid_value
  ___ubsan_handle_mul_overflow          ___ubsan_handle_negate_overflow
$ nm -u lib12-base/bgms/libs/bgms.so | grep -c __ubsan_handle     # control: the ordinary build
0
```

**ASan: BLOCKED, not clean — this needs saying plainly.** The package builds and
links against `libclang_rt.asan_osx_dynamic.dylib`, but the runtime has to be
loaded before R, and it cannot be on this machine:

```
==28587==ERROR: Interceptors are not working. This may be because AddressSanitizer is loaded too late (e.g. via dlopen). Please launch the executable with:
DYLD_INSERT_LIBRARIES=/Library/Developer/CommandLineTools/usr/lib/clang/21/lib/darwin/libclang_rt.asan_osx_dynamic.dylib
"interceptors not installed" && 0
ERROR: loading failed
```

Setting `DYLD_INSERT_LIBRARIES` does not help, because this R is signed with the
hardened runtime and dyld drops the variable:

```
$ codesign -dv /Library/Frameworks/R.framework/Resources/bin/Rscript
CodeDirectory v=20500 size=755 flags=0x10000(runtime) …
$ DYLD_INSERT_LIBRARIES=$ASAN_LIB Rscript -e 'cat(Sys.getenv("DYLD_INSERT_LIBRARIES"))'
                                    # empty
```

Re-signing the user's R installation is not something I will do for a
diagnostic. **So: no ASan evidence either way.** The way to get it is a
container with an ASan-instrumented R — CRAN's `clang-ASAN`/`gcc-ASAN` check
flavours, or `rhub`'s `clang-asan` image. Recommend running that once before
release regardless of this finding; it is cheap and it is the one class of bug
UBSan does not cover.

## A.3 Object-file bisection

Harness: relink `bgms.so` from a mixture of the two builds' objects (identical
link line, only the object paths vary), drop it into a fixed library, run the
seeded fit. Both endpoints reproduce, so the harness is sound:

```
all 44 objects from the fix build   -> 2.431
all 44 objects from the base build  -> 2.865
```

Binary search over the 44 objects, sorted:

```
base objects [1..22] (22 files) -> 2.865
base objects [1..11] (11 files) -> 2.865
base objects [1..6]  (6 files)  -> 2.431
base objects [7..9]  (3 files)  -> 2.431
base objects [10..10] (1 file)  -> 2.431
CULPRIT: mcmc/algorithms/hamiltonian_utils.o
```

Necessary *and* sufficient, confirmed both ways:

```
ONLY hamiltonian_utils from base     -> 2.865
ALL BUT hamiltonian_utils from base  -> 2.431
```

One translation unit, swapped alone, flips the whole result in either
direction. That is expected of a chaotic sampler: `hamiltonian_utils` is on the
innermost NUTS path, so a one-ULP difference there reaches every subsequent
draw.

## A.4 Disassembly

Whole-unit instruction mix for `hamiltonian_utils.o`:

| | anomalous | good (`-O2`) |
|---|---|---|
| instructions | 8 427 | 2 980 |
| `ldr` / `str` | 1 447 / 1 088 | 269 / 125 |
| `fmul` / `fadd` | 27 / 8 | 63 / 25 |
| `fmadd` | 9 | 9 |

The load/store ratio is the tell — the anomalous build spills every value to
the stack and reloads it. Its entry sequence is textbook `-O0`, including a
branch to the next instruction:

```
0000000000000500 <__Z13__ieee754_logd>:      # anomalous
     500:  sub  sp, sp, #0x90
     504:  str  d0, [sp, #0x80]
     508:  b    0x50c                        # <- branch to the following instruction
     50c:  ldr  d0, [sp, #0x80]
     510:  str  d0, [sp, #0x18]
     514:  ldr  w8, [sp, #0x1c]
```

```
00000000000002c0 <__ieee754_log(double)>:    # -O2
     2c0:  fmov x8, d0
     2c4:  lsr  x10, x8, #32
     2c8:  cmp  w10, #0x100, lsl #12
     2cc:  b.ge 0x2fc
```

**No sign of uninitialized reads.** No reads of stack slots that are never
written, no `mov` of an undefined register into an FP operand — the anomalous
version's loads all pair with a preceding store of the same slot, which is what
a spill looks like and not what garbage looks like. UBSan agrees (`
__ubsan_handle_load_invalid_value` never fires).

**FMA contraction is where the arithmetic actually diverges.** Tail of
`__ieee754_log` in `math/custom_explog.o`, same function, both builds:

```
-O2:            unoptimized:
  fmul  d4, d1, d4        fsub  d1, d1, d3
  fmadd d2, d2, d3, d4    fmsub d0, d0, d1, d2
  fsub  d0, d2, d0        fsub  d2, d2, d3
  fmsub d1, d2, d4, d3    fmul  d3, d3, d4     <- fmul + fsub, two roundings
  fsub  d0, d0, d1        fsub  d2, d1, d2
  fmsub d0, d2, d3, d0    <- fused, one rounding
```

## A.5 Confirming the mechanism

If contraction is the mechanism, turning it off should move the answer again —
and it does. Four legal compilations of `afe11510`, same seed, same data:

| build | `A (diff1; 1)` | `A-B` (pairwise) | `fmadd` in culprit TU |
|---|---|---|---|
| `-O2` (`wt-fix12-base`) | 2.8655 | 0.0418 | 9 |
| `-O2` (`f121/wt-o2-third`, **different path**) | **2.8655** | **0.0418** | 9 |
| `-O1 -fsanitize=undefined` | 2.4306 | 0.0374 | — |
| `-O2 -ffp-contract=off` | 3.0355 | 0.0396 | **0** |

Three things to read off this table.

1. **The two `-O2` builds at different paths are identical to the last reported
   digit** — every parameter, every `n_eff`, not just this one. The package is
   reproducible across builds at its own settings. The §3.4 finding as
   originally stated ("bgmCompare draws are not reproducible across builds") is
   **too strong and should be retired** in favour of this.
2. **`-ffp-contract=off` zeroes the `fmadd` count and moves the answer.** That
   is the mechanism, demonstrated rather than asserted: fused multiply-add
   rounds once where the unfused pair rounds twice, the sampler's trajectory is
   chaotic, and a one-ULP divergence at iteration 1 gives an effectively
   independent chain by iteration 20.
3. **The spread is confined to the parameters that have no data behind them.**
   The well-identified pairwise effects move by at most 0.006 across all four
   builds (`A-C`: −0.0040 to −0.0100; `A-B`: 0.0374 to 0.0418) — around two
   Monte-Carlo standard errors, which is what two short chains from the same
   posterior look like. The threshold *differences*
   swing 2.43 → 3.04, 25 % — and those are exactly the prior-driven rows that
   F-117 now marks with a `*`. Their `sd` is 1.0–4.8 on `n_eff` 20–106, so even
   the 25 % swing is roughly one posterior SD: two short chains from the same
   posterior, not two different posteriors.

## A.6 Verdict and what I would do about it

**Legal floating-point compilation variance.** Identical source, identical
semantics, different rounding from a legal optimization (FMA contraction, and
whatever else `-O2` reorders), amplified into a visibly different short chain by
a chaotic sampler. Not undefined behaviour: UBSan is clean across the compare
suite with the instrumentation verified present, and the disassembly shows
spills, not garbage. **No fix is called for in the sampler.** The one honest gap
is ASan, which this machine cannot run (§A.2).

Three suggestions, none of them code, all for a separate batch:

- **Retire the §3.4 claim as originally worded.** Two `-O2` builds agree
  exactly. What is true is the weaker and more useful statement: *bgmCompare
  results are bit-reproducible for a given compilation, and are not portable
  across different optimization settings or compilers.* That belongs in the
  docs next to `seed`, because users will compare a Mac result to a Linux one
  and find they differ.
- **Never snapshot a posterior mean.** Already the practice in this file
  (`labels_only()`), and the reason is now measured rather than suspected.
- **Run one ASan check before release**, in a container with an instrumented R.
  It is the only class of defect this diagnosis leaves unexamined.

## A.7 Artifacts

Everything is under `~/bgms-review/f121/` (~100 MB), nothing inside the repo:

- **`obj-fix12/`** — the anomalous build's 44 objects, preserved before any
  rebuild. **The only surviving copy**, and not reconstructible: the compile
  that produced them is the one thing here I could not reproduce (§A.1). If any
  of this is to be kept, keep this.
- `relink.sh` + `bisect.sh` + `all.txt` — the bisection harness. It reads the
  `-O2` reference objects from the worktree `~/bgms-review/wt-fix12-base`, so
  that worktree has to stay for the bisection to be re-runnable.
- `probe.R` / `probe2.R` / `snapcheck.R` — the seeded fit reduced to its
  numbers, and the build-invariance check on the snapshot transform.
- `dis-*.txt`, `sym-*.txt`, `s-*.txt` — disassembly and symbol diffs.
- `lib-ubsan/`, `lib-san/`, `lib-nofma/`, `lib-o2-third/` — the installed
  sanitizer and comparison builds, still runnable via `probe.R`.
- `makevars-{san,ubsan,nofma}` — the flag overrides. Note these set
  `CXX20FLAGS` / `SHLIB_CXX20LDFLAGS`, not `CXXFLAGS`, because the package is
  `CXX_STD = CXX20`; the plain names are silently ignored, which cost me a
  build.

The three throwaway worktrees the diagnosis used (`wt-san`, `wt-nofma`,
`wt-o2-third`) are already removed — they are trivially recreated with
`git worktree add --detach <path> afe11510` plus the matching `makevars-*`.
`~/bgms-review/wt-fix12` (the branch) and `~/bgms-review/wt-fix12-base` (the
`-O2` reference) are still registered. Delete the lot when the finding closes.

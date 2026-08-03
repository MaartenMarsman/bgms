# Report 10 — the edge-posterior panel, redrawn to the JASP standard

**Brief:** 10 — edge-posterior panel redesign (F-066)
**Branch:** `feat/edge-posterior-panel`
**Base:** `fix/bgmcompare-batch` @ `1508a888` (see *Base of the branch*, below)
**Commits:** `a682b57e`, `7cb2b837`, `f9ac7526`
**Worktree:** `~/bgms-review/wt-fix4`

---

## Summary

`plot_edge_posterior()` is now the R Graph Compendium's JASP prior-and-posterior
figure, and the style it is drawn in has been factored out as
`R/plot_style.R`, written as the package's plotting law rather than as this
panel's local choices — brief 11 restyles everything else onto it.

The panel gained a prior, a numbered density axis, a credible interval, a
median, a probability wheel, and a caption; it lost its headline title and its
`binwidth` argument. Which estimator it displays follows the model the fit
actually used: a spike-and-slab fit gets the Rao-Blackwellized indicator Bayes
factor and **no** Savage-Dickey dots; a fit without edge selection gets the JASP
figure exactly, ordinate dots and all.

Three defects were found and fixed while rewriting: a decisively-absent edge
raised an error instead of drawing, a mixed fit's cross edge could not be named
in variable order, and a `bgmCompare()` fit was silently drawn as though its
baseline pairwise effect were one network's edge.

---

## What was done

### 1. `R/plot_style.R` — the style as a module (new file, 353 lines)

Roxygen-commented internal helpers, each documented as a rule rather than a
utility. The file header states the six conventions the package now follows,
in the order they matter:

1. **No big titles.** A JASP plot carries annotations, not a headline.
2. **Offset axes.** Each axis spans only its own tick range (a Tufte range
   frame; `geom_rangeframe()` in `jaspGraphs`), with the data region grown by
   an epsilon so the frame cannot close at the corner.
3. **No box** (`bty = "n"`).
4. **Large fonts.** `jaspGraphs` works at a 17 pt base with axis titles at
   1.2× and legends at 1.25×; the compendium's base-R figures run
   `cex.axis = 1.2`, `cex.lab = 1.5`, `lwd = 2`. Both are distilled into
   `bgms_style()`.
5. **Needless ink omitted.** Grey ink, no grid, no legend box, low-alpha fills,
   one accent from `mover_palette()`.
6. **Probabilities are wheels.** Any probability a panel reports is also a
   filled wheel; printed text carries log Bayes factors, which are not
   probabilities and get no wheel. This is the same encoding the compare plots
   put on their main-difference node rings (F-062d).

| helper | what it owns |
| --- | --- |
| `bgms_style()` | the constants: colours, `cex` multipliers, line weights, `mar`, the axis epsilon, the wheel radius |
| `bgms_panel_par()` | opens a panel in the style; returns the constants plus `old_par` for the caller's `on.exit()` |
| `bgms_axis_range()` | `pretty()` ticks plus the wider data range that offsets the axis from them |
| `bgms_axis()` | draws one offset axis |
| `probability_wheel()` | the `cos()`/`sin()` + `polygon()` wheel |
| `annotation_block()` | a corner stack of annotation lines, spaced by their own height |
| `panel_x()` / `panel_y()` | figure-relative placement in user coordinates (JASP's `grconvertY(..., "ndc", "user")` device) |
| `format_probability()` | a probability printed with its relation, mirroring `format_log_bf()` |

Two implementation notes worth the maintainer's eye:

- **The wheel radius is in inches**, converted per axis through
  `par("usr")`/`par("pin")`, so a wheel is round on any device and at any axis
  scale, and can be drawn in a margin under `xpd = NA`.
- **The minority share is always the wedge.** A wedge covering nearly the whole
  circle closes on itself and leaves a visible seam; drawing the small share on
  a full disc of the majority colour keeps every wedge under half a turn. The
  accented share still sits at the top of the wheel either way, because the
  pale complement of a top-centred wedge is a bottom-centred one.

### 2. The panel

**Layout** (top margin, figure-relative, so it holds while the axis scale
moves under it — JASP's own device):

```
  intrusion-upset                              median = 0.09
  undecided                                    95% CI [0.03, 0.15]

  (wheel)  PIP = .79
           log BF = 1.3
  ------------------------------- plot region -------------------------------
  legend "Posterior / Prior"; prior dashed grey, posterior solid accent + fill
  95% CI bar above the curve
  ------------------------------------------------------------------------
  Edge weight
  Density: the weight given inclusion. Pale share of the wheel: P(absent).
```

**The prior** is computed in closed form from the fit's own
`interaction_prior`, read from `get_fit_spec()$prior` — Normal, Cauchy, or
beta-prime (the last through its logistic Jacobian, with the tails floored at
zero where the Beta density overflows while the Jacobian underflows). Nothing
is hardcoded; the default changed to `normal_prior(scale = 1)` in 0.2.0 and a
`cauchy_prior(scale = 2.5)` fit now draws a visibly different curve.

**No change of variable is needed**, which was not obvious and is worth
recording. The sampler evaluates the slab on the association-scale parameter —
for a continuous block that is `-K_ij/2`, not the precision entry
(`ggm_model.cpp:310`, `mixed_mrf_gradient.cpp:574`) — and
`extract_pairwise_interactions()` reports the draws in exactly that frame for
every model type: it returns, element for element, the `theta` that
`anchor_draws()` reweights the anchored sensitivity curve on, which is *defined*
as the frame the slab applies to. A test pins that identity so a future change
to either side breaks loudly rather than silently misdrawing the prior.

**The dots rule**, as briefed:

| fit | evidence shown | ordinate dots | wheel filled to | wheel labels |
| --- | --- | --- | --- | --- |
| `edge_selection = TRUE` | Rao-Blackwellized indicator log BF | **none** | posterior inclusion probability | none (caption carries the pale share) |
| `edge_selection = FALSE` | Savage-Dickey log BF | **both, grey, at zero** | `BF/(1 + BF)` | `data\|H1` / `data\|H0` |

The posterior ordinate at zero comes from `stats::density(bw = "SJ")`
interpolated at 0; the prior ordinate is exact. The Rd documents that JASP's
own implementations use a logspline fit and says why `bgms` does not adopt it
(a dependency for one number), and where the two estimators disagree.

### 3. Edge cases

| case | behaviour |
| --- | --- |
| PIP ≈ 1 (saturated) | wheel drawn full, `PIP > .99`, `log BF > 10,000`; median/CI of the conditional posterior unaffected |
| PIP ≈ 0 / < 2 included draws | **no longer an error**: prior curve, near-empty wheel, evidence, and a caption saying no retained draw included the edge; no median/CI, because none was sampled |
| no edge selection | the Savage-Dickey branch above; snapshotted at both ends — a posterior that has left zero (dots at 0.399 and ≈ 0.000, log BF = 7.3) and one piled on zero (0.399 and 13.07, log BF = −3.5, wheel 3 % filled) |
| Blume-Capel / mixed | nothing special, as briefed; both are exercised end to end, and the mixed one turned up a real bug (below) |

### 4. API

- `binwidth` is `lifecycle::deprecated()`: warned about, ignored, badged in the
  Rd, with a message saying what replaced it.
- `evidence_threshold` still feeds the verdict, now printed as an annotation
  rather than in a title.
- No other signature change. `edge_panel_title()` is retired with the titles it
  built, along with its unit test and its snapshot.

---

## Findings

Three of these are defects the redesign uncovered rather than introduced, so
they want their own IDs in the ledger; the proposed numbers are suggestions for
the review lead, F-065 being the current high-water mark.

| # | proposed ID | severity | status | finding |
| --- | --- | --- | --- | --- |
| 1 | F-067 | **major** | fixed | a decisively-absent edge raised an error instead of drawing |
| 2 | F-068 | **major** | fixed | a mixed fit's cross edge could not be named in variable order |
| 3 | F-069 | **major** | fixed | a `bgmCompare()` fit was drawn as though its baseline pairwise effect were one network's edge |
| 4 | F-066 | major | fixed | the shipped panel had no prior, no interval, no median and an unnumbered y axis — the maintainer's complaint, and this brief |
| 5 | — | note | documented | on a continuous block the drawn prior is the slab, not the marginal |
| 6 | — | note | decision wanted | the RB inclusion probability on the wheel and the empirical share of nonzero draws behind the density are two different numbers |

### Finding 1 (F-067) — decisive absence was an error (major; fixed)

The shipped function ended with

```r
if(length(slab) < 2L) stop("The edge ... was included in fewer than two retained draws ...")
```

so the one figure a reader most wants for a settled null — "the data rule this
edge out" — was the one the package refused to draw. It now draws the prior,
the wheel (nearly all pale) and the evidence, and its caption says no retained
draw included the edge. Reproduced with four independent binary variables under
`bernoulli_prior(0.001)`, which leaves five of six edges at exactly zero
included draws.

### Finding 2 (F-068) — a mixed fit's cross edge was reported missing (major; fixed)

`plot_edge_posterior()` built the label as `variables[first]-variables[second]`
with `first < second` in *column* order and tested membership in one
orientation only. A mixed fit lays its pairwise draws out by block —
discrete-discrete, continuous-continuous, cross — and names a cross edge by its
**discrete** end, which need not be the earlier column. On the four-variable
mixed fixture the pair `(c1, d2)` is stored as `"d2-c1"`, so

```r
plot_edge_posterior(fit, "c1", "d2")
#> Error: No edge between 'c1' and 'd2'.
```

for an edge that exists and is decisive (log BF ≈ 30). `edge_column_label()`
now tries both orientations, in the pairwise draws and in the `verdicts()`
table. This is a pre-existing defect, not one the redesign introduced; the same
family as F-047 (both ends resolving to one variable) and F-037 (reading a
mixed fit's indicators in its own layout).

### Finding 3 (F-069) — a compare fit was drawn as a single network's edge (major; fixed)

Nothing in the shipped function checked the class. Given a `bgmCompare()` fit
it read `extract_pairwise_interactions()`, which for a compare fit returns the
**baseline** pairwise effects, and drew them with a zero spike (baselines are
never excluded), under a title carrying a *difference* verdict from
`verdicts()`. The picture and the number in it were about different parameters.

The redesign would have made this worse rather than better: a compare fit's
`extract_arguments()` has no `edge_selection` field, so `isTRUE(...)` is
`FALSE` and the fit would have fallen into the Savage-Dickey branch — a
licensed-estimator claim on a model whose selected quantity is the difference
indicator. It is now refused with a message pointing at `plot(fit)` and
`verdicts(fit)`.

This is the one behaviour change beyond the brief's scope; the dots rule forced
it, and it is flagged for the maintainer rather than assumed.

### Finding 4 (F-066) — the panel itself (major; addressed)

Addressed by the redesign; the before/after figures are the evidence.

### Finding 5 — the slab is not the marginal on a continuous block (note)

On a block of continuous variables the slab sits on an entry of a precision
matrix, whose joint prior is the slab times a prior on the diagonal, restricted
to the positive-definite cone. The drawn curve is the slab — the density the
sampler evaluates for that edge, and the one the prior-sensitivity machinery
reweights on — which is the marginal prior of that entry only up to the
restriction.

Measured, so the size is on record rather than asserted: with `p = 4`,
`normal_prior(scale = 1)`, `exponential_prior(eta = 1)`, `spec = "joint"`, a
40,000-draw `sample_ggm_prior()` run gives an empirical SD of 1.34 on the
included `K_offdiag` draws against the unconstrained slab's 2.0 — the
restriction pulls the true marginal in by about a third. On a discrete block
the parameters are unconstrained and the slab is the marginal exactly.
`?plot_edge_posterior` says all of this.

### Finding 6 — two inclusion numbers (note)

The wheel is filled to the Rao-Blackwellized `pip` from `verdicts()`, so that it
agrees with the log Bayes factor printed beside it (both are RB quantities).
The density beside it is conditional on the *empirical* inclusion — it is the
kernel density of the nonzero draws. The two can differ by Monte Carlo error.
This is the right pairing (evidence with evidence), but it is a choice, and it
is listed here so the maintainer can overrule it. If the wheel is instead
filled to `mean(draws != 0)`, the wheel and the curve become exactly
commensurable and the wheel stops agreeing with the log BF.

---

## Evidence

### Before / after

Same fit (`Wenchuan[, 1:6]`, 2 chains, 400/400, seed 1), same three edges.

| edge | before | after |
| --- | --- | --- |
| decisive (`intrusion-dreams`, log BF 277) | `assets/edge-panel-before-decisive.png` | `assets/edge-panel-after-decisive.png` |
| undecided (`intrusion-upset`, log BF 1.3) | `assets/edge-panel-before-undecided.png` | `assets/edge-panel-after-undecided.png` |
| absent (`intrusion-avoidth`, log BF −3.5) | `assets/edge-panel-before-absent.png` | `assets/edge-panel-after-absent.png` |

The undecided pair is the clearest read: before, a bold two-line title, a black
stem labelled `.21`, a slab with no scale, an unnumbered y axis and a caption
apologising for the units. After, a numbered density axis, the prior it was
updated from, `PIP = .79` on a wheel filled to that share, `log BF = 1.3`,
`median = 0.09`, `95% CI [0.03, 0.15]` with the interval drawn as a bar, and a
one-line caption stating the conditioning.

### Gate

All four gate items pass. Everything below was run outside the Dropbox tree
(F-009 standing practice) on a `git archive` export of this branch's HEAD
(`f9ac7526`).

**1. Test suite** — `~/bgms-review/run-fix4-suite.R`, against the tarball
install in `~/bgms-review/lib-fix4`:

| tier | files | tests | failed | errors | warnings | skipped | elapsed |
| --- | --- | --- | --- | --- | --- | --- | --- |
| full (`NOT_CRAN=true`) | 78 | 1187 | **0** | **0** | 0 | 105 | 8.2 min |
| CRAN mode | 78 | 1187 | **0** | **0** | 0 | 251 | 2.1 min |

Logs: `~/bgms-review/suite-fix4-{notcran,cran}.log`.

The 15 new tests in `test-plot-methods.R` and the 9 in `test-plot-style.R` all
run in the full tier; seven snapshots are recorded in
`tests/testthat/_snaps/plot-methods.md` (the six briefed cases plus the
Blume-Capel structural one). The retired `edge_panel_title` snapshot is gone
from that file.

> **A false alarm worth recording, because the next agent will hit it.** The
> first run of the full tier reported 26 failures and 64 errors. All of them
> were an artefact of the runner: `run-fix4-suite.R` prepends `lib-fix4` to
> `.libPaths()`, but the suite was launched before `R CMD INSTALL` had finished
> writing it, so `library(bgms)` silently fell back to the **system** library's
> older 0.2.0.0 build — one without the F-057 group-indicator fix, the F-021
> calibration stub, or the current zratio helpers. The runner now hard-stops
> unless `system.file(package = "bgms")` resolves inside `lib-fix4`; the
> numbers above are from the guarded rerun. Nothing in the package was at
> fault, and no test was changed to make this pass.

**2. `R CMD check --as-cran`** on `bgms_0.2.0.0.tar.gz` (1.62 MB, built with
pandoc on `PATH`, vignettes included):

```
Status: 2 NOTEs
```

Both are the known baseline NOTEs from report 01 and neither is a code defect:

- `checking CRAN incoming feasibility ... NOTE` — "The Date field is over a
  month old" (F-003, to be bumped in the submission commit).
- `checking HTML version of manual ... NOTE` — macOS ships an HTML Tidy too old
  to run, so validation is *skipped* rather than failed (report 01 finding 10).

No ERROR, no WARNING. Notable lines:

```
* checking R code for possible problems ... [9s/20s] OK
* checking Rd files ... OK
* checking for code/documentation mismatches ... OK
* checking examples ... OK
* checking examples with --run-donttest ... [11m/12m] OK
* checking tests ...
  Running 'testthat.R' [127s/136s] OK
* checking re-building of vignette outputs ... OK
```

Log: `~/bgms-review/check-fix4.log`; check dir `~/bgms-review/checkfix4`.

**3. The Rd example** runs against the installed build:

```
example("plot_edge_posterior", package = "bgms", run.donttest = TRUE)
#> RD EXAMPLE OK
```

It also runs inside the check, under `--run-donttest`, above. Rd regenerated
with roxygen2 8.0.0; `man/plot_edge_posterior.Rd` carries the deprecation badge
on `binwidth`.

**4. The tarball is clean** — no `dev/`, no assets, and both new files ship:

```
bgms/R/plot_style.R
bgms/tests/testthat/test-plot-style.R
bgms/tests/testthat/test-plot-methods.R
```


---

## Deviations from the design, and why

1. **Legend placement.** The brief says "top-left". Implemented as the
   compendium's own rule instead: left unless most of the posterior mass sits
   left of centre, in which case right (`JASPttest1.R` chooses by
   `mean(delta > mean(range(xticks)))`). Edge weights can be negative, so a
   fixed top-left legend would sometimes be drawn over the curve. **This is the
   one deviation most worth a decision.**
2. **A credible-interval bar was added.** The brief lists five elements and
   does not mention it; the reference figure has it (its element 3), and the
   panel would otherwise print an interval it does not show. Drawn as JASP
   draws it, `arrows(angle = 90, code = 3)` above the curve.
3. **No title at all.** The brief left this to judgement against the reference
   figures; the reference has none, so the panel has none. The edge and its
   verdict became the first two annotation lines. `edge_panel_title()` and its
   snapshot are gone.
4. **The prior is read from the spec, not from `extract_arguments()`.** The
   brief says to read the family and scale from `extract_arguments()`; the
   family is not in `extract_arguments()` at all, and a GGM fit's arguments
   carry no `pairwise_scale` either. `get_fit_spec()$prior` has both for every
   model type. A fit carrying no spec draws no prior curve rather than a
   guessed one.
5. **`binwidth` uses `lifecycle::deprecated()`/`is_present()`**, not the
   `hasArg()` form `bgm()` uses for its deprecated arguments. `bgm()`'s
   deprecated arguments have no defaults; this one does, and the `deprecated()`
   sentinel keeps the signature self-documenting in the Rd. Same `lifecycle`
   dependency, same warning class.
6. **The `bgmCompare` guard is a behaviour change the brief did not ask for.**
   Finding 3 above.
7. **Eight helpers, not three.** The brief named `bgms_panel_par()`,
   `probability_wheel()` and `annotation_block()`. The offset-axis convention
   needs a range computer and a drawer to be a convention rather than a
   copy-paste, and figure-relative placement needs `panel_x()`/`panel_y()`;
   `bgms_style()` holds the constants all of them read.
8. **The compare node rings still use qgraph's `pie` channel.**
   `probability_wheel()` is written to serve them, but converting them is a
   change to `plot.bgmCompare()`, which is brief 11's territory. The visual
   language is already shared (fill fraction = probability, accent = the thing
   asserted); only the drawing code is not.
9. **Snapshots are of the panel's content, not of a raster.** The six cases are
   snapshotted through a `describe_panel()` helper that prints the label,
   subtitle, evidence lines, estimate lines, wheel value, wheel tags, which
   curves and dots are present, the window, and the caption. Snapshotting
   pixels would pin the device rather than the figure and would fail on every
   platform but this one.

### Base of the branch

The brief says to base on `develop` at or after the brief-06 merge. That merge
has not happened: `develop` is at `289d028d` and `fix/bgmcompare-batch` is at
`1508a888`. The latter is a **strict superset** of `develop`
(`git log fix/bgmcompare-batch..develop` is empty), so it is exactly what
`develop` becomes when brief 06 merges — and it carries the F-062d ring work
this panel's wheel convention is meant to match. The branch is based there. It
will fast-forward cleanly once brief 06 lands; if brief 06 is revised before
merging, this branch must be rebased.

---

## Open questions

1. **Legend placement** — adaptive (as built, and as the compendium does) or
   pinned top-left (as briefed)?
2. **Two or three decimals on the estimate?** The brief's example reads
   `median = 0.31`, so the panel prints two. JASP prints three
   (`formatC(, 3, format = "f")`). Two decimals lose resolution on a weak edge:
   a weight of 0.045 prints as `0.04`, and a decisive-absence interval can
   print as `[-0.06, 0.06]` around `0.00`.
3. **Which inclusion number fills the wheel** — Finding 6.
4. **Should a Savage-Dickey reading ever be offered for a selection fit?** The
   dots rule says no by default, and no argument was added. If the maintainer
   wants it available for comparison (it is what other packages report), it
   would be an argument on this function rather than a default.
5. **Does `plot_edge_posterior()` want a `bgmCompare` sibling** for one
   difference, now that the compare fit is refused rather than mis-drawn? That
   is a new function, not a fix; noting it because the refusal makes the gap
   visible.
6. **`R/plot_style.R` is not exported.** Brief 11 will use it internally.
   If the JASP style is meant to be reusable by JASP itself or by `easybgm`,
   some of it (`probability_wheel()` in particular) may want to be public.

# Report 11 — package-wide plot restyle and compare evidence parity

Branch `fix/plot-sweep`, based on `origin/develop` at `a930b5a5`.

**Base note.** The brief asked for a base at or after `f24ad3c8` (the brief-12
merge). The repository's local `develop` ref is 23 commits stale and does *not*
contain it; `origin/develop` (`a930b5a5`) does, and `f24ad3c8` is its
grandparent. The worktree was cut from `origin/develop`. Nothing was checked
out or built in the Dropbox repository.

---

## What was done

### F-067 — every remaining plot onto the style module (task 1)

| entry point | commit | what moved |
| --- | --- | --- |
| `plot.bgms` (network) | `3f526930` | Whole method rewritten onto a shared renderer; module ink, type, legend, caption, annotation band |
| `plot.bgmCompare` (difference, groups) | `3f526930` | Same renderer; group panels lose their qgraph `title` banner for module labels, type scaled for the three-panel layout |
| `plot.bgms_centrality` | `b6633d71` | Module par/axis/ink; measured left margin; quantity on the axis, marks in a caption |
| `plot.bgms_calibration` | `fe96e192` | Per-panel `main =` → module label + subtitle; offset axes on the unit square; whole style scaled for the grid |
| `plot.bgms_prior_sensitivity` | `c5afad70` | Headline `main =` retired to a module label; module axes and ink; corner notes to the caption |

A follow-up, `0957cd8a`, keeps the calibration figure's outer axis labels and
caption at full size: they describe the whole figure rather than one panel, so
they must not shrink with the grid. At nine panels the caption had come out at
0.44 of the base size.

Supporting module work is in `e9d1d153`: `bgms_style(scale =)` so a small
multiple is the same style seen smaller rather than a second style;
`bgms_caption()` and `bgms_panel_label()` as the two placements a title-less
panel is allowed; `margin_lines_for()` to size a margin from the text it has
to hold; `display_log_bf()` to round away a `-0.0` (see findings).

No `main =` banner is left on any figure the package draws. Grep for it:
`plot.bgms_calibration` and `plot.bgms_prior_sensitivity` were the only two
that still passed one, and both now pass `main = ""`.

### Task 2 — the verdict word off the edge panel (`7b17c96e`)

`edge_panel_selection()` sets `subtitle = NULL` in every case. `verdict_phrase()`
is deleted, the `verdict` argument is gone from `edge_panel_selection()`, and
`edge_selection_evidence()` no longer reads a verdict. `"no edge selection"`
stays on the Savage-Dickey variant: it names the model the panel is drawing,
not a verdict about the edge.

The network legends at `plot_bgms.R:187` / `:445` were **not** changed
unilaterally — see open questions, and the two rendered variants.

### Task 3 — "PIP" out of the rendered text (`e9d1d153`, `7b17c96e`)

`bgms_style()$label_inclusion` is `"P(included)"` and `format_inclusion()`
composes it with `format_probability()`'s relation, so the panel prints
`P(included) = .786` / `P(included) > .999`. No rendered string in `R/` says
"PIP" any more (`grep -rn "PIP" R/` returns only internal comments about the
quantity in `anchor_curve.R` and `prior_sensitivity.R`).

### Task 4 — three decimals via the style constant (`e9d1d153`)

`bgms_style()$prob_digits = 3L`; `format_probability(p, digits = ...)` derives
its cut-offs and its end strings from that constant rather than writing them
out, so `> .999` / `< .001` follow the digit count. Verified at two digits in
`test-plot-style.R`, which pins the old strings as the `digits = 2` case.

`format_log_bf()`'s precision is **unchanged** at one decimal. The renders did
not argue for more: a log Bayes factor at one decimal is already finer than the
Monte Carlo error of the fits these figures are drawn from, and the reporting
cap (`> 10,000`) is what actually governs the extreme cases. What the renders
*did* surface is a `-0.0` at the rounding boundary — reported as a finding
rather than fixed in `R/verdicts.R`, which this brief does not own.

### Task 5 — F-079, compare evidence parity (`3f526930`)

See the parity paragraph below.

### Task 6 — compare node rings become wheels (`3f526930`)

`main_difference_nodes()` returns `prob`/`color` instead of `pie`/`pie_color`;
`qgraph`'s `pie`/`pieColor` arguments are no longer passed. Each node's
`P(main-effect difference)` is drawn by `probability_wheel()` beside the node,
pushed radially outward from the centre of the layout so it lands in open space
and clears the node's own border and label. qgraph's `mar` widens when the
nodes carry wheels, so the wheels have somewhere to sit that is not the
annotation band.

### Task 7 — F-070, sensitivity right-margin clipping (`c5afad70`)

The old code fixed `mar[4] = 7.5` lines and placed each label at
`10^(usr[2] + 0.015 * diff(usr[1:2]))` — a fraction of the *data* range, in a
margin sized by a constant. Whether a name fit therefore depended on how long
the name was and how wide the device happened to be, which is the clip.

Now: the right margin is `margin_lines_for(edges$edge[named], cex = name_cex)`
— the widest name actually being drawn, measured in inches and divided by the
device's own line height `par("csi")` — plus the leader length, also in inches.
The leader stub and the label are placed at `10^(usr[2] + inches / per_log_unit)`
where `per_log_unit = pin[1] / diff(usr[1:2])`, i.e. device lengths converted
into the log-scale coordinate, not fractions of a range. Nothing is
hard-coded to a width, and the fix is visible at the default device size in
`sensitivity-differences-{before,after}.png`.

---

## The parity paragraph (task 5)

Parity ended up meaning **one renderer, not two agreeing ones**. `plot.bgms`
and `plot.bgmCompare` previously reached `qgraph` through separate bodies that
happened to encode edges the same way; nothing structural stopped them from
drifting, and in fact only the compare side had ever been given anything to say
about its evidence (two subtitle notes), while neither printed a Bayes factor
at all. Both now call `draw_verdict_network()`, which builds the qgraph input,
draws the node probability marks, composes the evidence band, places the legend
and the caption. The two methods differ in exactly two arguments: the strings
`network_unit("edge")` versus `network_unit("difference")` hands over, and
whether `nodes` carries a probability per node. The evidence band itself is
`evidence_band()`, one function called by both: it tallies the units by verdict,
prints the threshold as `format_log_bf(log(evidence_threshold))`, and prints the
strongest log Bayes factor for and against a unit through the same
`format_log_bf()` the edge panel uses — so a difference Bayes factor is on the
natural-log scale, capped at the same reporting cap, rendered by the same
formatter, and laid out in the same top-margin annotation block as an edge
Bayes factor, because it is literally the same code path. `test-plot-methods.R`
pins that with an assertion that `evidence_band()`'s right-hand block is
*identical* between the edge and difference units for the same numbers.
`ln 10 ≈ 2.3` appears as the printed threshold on both by construction rather
than by two call sites agreeing. The single-DIFFERENCE-edge panel remains out
of scope (post-release item 42) and was not built.

---

## Proposed NEWS clauses (verbatim; the lead lands these)

Under **New features → Plots and centrality**, replacing the `plot()` on a
`bgmCompare()` fit bullet's main-effect sentences and extending both network
bullets:

> * `plot()` on a `bgm()` fit and `plot()` on a `bgmCompare()` fit now draw
>   their evidence through one routine, so an edge Bayes factor and a
>   difference Bayes factor are displayed the same way rather than similarly.
>   Both carry a band above the network naming the display, tallying the units
>   by verdict, stating the threshold on the natural-log Bayes factor scale,
>   and printing the strongest Bayes factor for and against a unit in the fit.
>   Neither picture stated any of that before, and the compare picture in
>   particular showed group differences with no statement of the evidence
>   behind them.

> * Main-effect difference evidence rides on the nodes as the package's own
>   probability wheel, filled to the difference indicator's posterior inclusion
>   probability and placed beside the node, rather than as a qgraph pie ring.
>   It is the same wheel `plot_edge_posterior()` draws, so a probability looks
>   the same on every figure the package produces. Under the default
>   `main_difference_selection = FALSE` those indicators do not exist, no wheel
>   is drawn, and a note under the tally names the setting.

> * Every figure the package draws now follows one set of conventions: no
>   headline titles, offset axes, no plot box, muted grey ink, larger type, and
>   probabilities shown as filled wheels. `plot()` on a calibration check, on a
>   centrality object and on a prior-sensitivity check moved onto them, along
>   with both network methods. The prior-sensitivity panel's title becomes a
>   compact label above the panel; its answer is still the first thing read.

Under **Bug fixes**:

> * The prior-sensitivity plot's edge names no longer clip at the right-hand
>   edge of the device. The margin is now measured from the widest name being
>   drawn and the labels are offset from the axis by a length in inches, so the
>   names fit at any device width rather than at the one the margin constant
>   was chosen for.

Under **Other changes**:

> * `plot_edge_posterior()` no longer prints a verdict word on the panel. What
>   the panel shows is the evidence — the wheel, the inclusion probability and
>   the log Bayes factor — and the reading those license is the reader's to make
>   at a threshold they choose; `verdicts()` is where the package states
>   verdicts, and it names the threshold it used. `plot_edge_posterior()`'s
>   `evidence_threshold` is still accepted and still validated, but neither
>   number the panel prints depends on it.

> * Posterior inclusion probabilities are printed as `P(included)` rather than
>   `PIP`, to three decimals rather than two, wherever the package draws one.

Two existing 0.2.0.0 clauses are now **factually wrong** and want correcting in
the same pass — note that one of them was already wrong before this branch:

* The `plot()` on a `bgmCompare()` fit bullet says main-effect differences show
  as "a square node with an accented border where the data settle a main-effect
  difference, a circle where they do not". The shipped code drew a qgraph pie
  ring filled to the inclusion probability, not squares and circles; that
  sentence was stale independently of this branch and is now doubly so.
* The `plot_edge_posterior()` bullet says the wheel matches "the ring encoding
  the `bgmCompare()` panels use". It is now a wheel on both sides.

---

## Findings

| id | severity | finding |
| --- | --- | --- |
| A | minor | `format_log_bf()` (`R/verdicts.R:104`) prints `"= -0.0"` for a Bayes factor that rounds to nothing — a sign the run did not establish. It is visible on a compare network whose strongest positive difference is near zero. `estimate_lines()` already rounds exactly this away for a weight ("the sign of a rounded-away quantity is not information the run established"), so the package has already decided the question; `format_log_bf()` was simply not given the same treatment. **Not fixed here**: `R/verdicts.R` is off-limits for this brief. Worked around at the figure call sites by `display_log_bf()` in `R/plot_style.R`, which zeroes a magnitude under 0.05 before formatting. The workaround should be deleted and the guard moved into `format_log_bf()` once its owner can take it, because `print()` and the summary tables reach `format_log_bf()` directly and are still exposed. |
| B | minor | `NEWS.md` describes the `bgmCompare()` main-effect node encoding as square-versus-circle nodes with accented borders. The code has drawn a pie ring filled to the inclusion probability since F-062d, and now draws a wheel. The NEWS sentence has never matched the shipped behaviour. Flagged rather than fixed (`NEWS.md` is off-limits); clause text above. |
| C | minor | `plot_edge_posterior(evidence_threshold =)` no longer changes anything the panel draws, now that the verdict word is gone: it is passed to `verdicts()` to locate the edge's row, and neither the inclusion probability nor the Bayes factor on that row depends on it. It is still accepted and validated, and the Rd now says plainly that it does not change the figure. Recommend retiring or soft-deprecating it at the next breaking change rather than silently keeping an inert argument; not done here because removing a documented argument is the maintainer's call, not a restyle's. |
| D | note | `plot.bgms` and `plot.bgmCompare` disagree about the same situation. When nothing is drawable, `plot.bgms` errors ("No edge reaches evidence of presence or sits undecided at this threshold"); `plot.bgmCompare` draws the nodes alone and says so under the tally, on the stated principle that an empty network is a result rather than a failure. The principle applies equally to a single network at a high threshold. Behaviour deliberately left alone — the brief's parity item is about the evidence display — but the two should be made to agree, and the compare behaviour is the better one. |
| E | note | The prior-sensitivity panel still prints the verdict words `presence` / `undecided` / `absence` as zone labels down its left edge. They were left because that plot's whole subject is whether a *verdict* moves with the scale, so the zones are the axis of the question rather than an assertion about one edge — but it is the one place a verdict word survives on a figure, and it is worth a ruling alongside the legend question below. |
| F | note | `main_difference_nodes()` still colours the node wheel by verdict as well as filling it to the probability — two channels for one quantity. On an undecided or ruled-out node that means a grey wedge on a pale grey wheel, which is the weakest contrast on any figure in the package. Deviation 8 authorised the pie→wheel change only, so the colour channel was kept; the module's own rule ("the fraction, not the colour, is the message") argues for dropping it and drawing every node wheel in the style accent. |
| G | note | qgraph exposes no accessor for the radius it actually drew a node at; `graphAttributes$Nodes$width` is the size parameter it was given. `node_probability_wheels()` therefore places the wheels against a factor calibrated on qgraph's output (`1.9 * width / 100`), chosen on the generous side, because a wheel placed too far out is merely further out while one placed too close is drawn over the node's label. It is documented as such in the source. A qgraph release that changes its node scaling would move the wheels, not break them. |
| H | note | `evidence_band()`'s tally counts the three verdict classes and silently omits an `NA` verdict, so a tally would not sum to the number of units if one ever appeared. It cannot on either method's pairwise family today (only the main-effect indicators go `NA`, and they are not tallied), and `verdict_network_input()` already guards `is.na` on the drawing side. Recorded so the assumption is written down rather than assumed. |

---

## Evidence

### Test suite

`devtools::test()`, local default tier, no slow environment variables set,
sequential.

```
$ Rscript -e 'devtools::test(reporter = "summary")'

... 78 context lines, no failure, warning or error mark in the stream ...

══ DONE ════════════════════════════════════════════════════════════════════════
exit=0
```

Counted from the reporter's own marks: **8253 passing expectations, 105 skips,
0 failures, 0 warnings, 0 errors** across 78 test files. Every skip is
tier-gated or `skip_on_cran()`; none is new to this branch. No slow-tier
environment variable was set (`BGMS_RUN_SLOW_TESTS` and
`BGMS_RUN_CERTIFICATION` both unset), and the run was sequential at
`Ncpus = 2`.

**Snapshots re-recorded** — one file, `tests/testthat/_snaps/plot-methods.md`,
five cases, all for the same two ratified wording changes and nothing else
(diff verified line by line):

| snapshot case | reason |
| --- | --- |
| `the panel of a decisive edge carries the wheel, not a stem` | `subtitle` drops "evidence of presence"; `PIP > .99` → `P(included) = .995` |
| `the panel of an undecided edge splits its wheel` | `subtitle` drops "undecided"; `PIP = .60` → `P(included) = .600` |
| `a saturated edge prints the capped Bayes factor and a full wheel` | `subtitle` drops "evidence of presence"; `PIP > .99` → `P(included) > .999` |
| `a decisive absence with no included draw is a figure, not an error` | `subtitle` drops "evidence of absence"; `PIP < .01` → `P(included) < .001` |
| `the panel reads a Blume-Capel fit like any other` | `subtitle` is now absent, so the snapshotted expression gains a `%||% "(none)"` |

No other snapshot file was touched. `tests/testthat/_snaps/verdicts.md` is
unchanged.

### Documentation

`devtools::document()` runs clean. Rd drift is confined to the four methods
whose behaviour changed:

```
man/plot.bgmCompare.Rd             | 24 ++++++++++++++++--------
man/plot.bgms.Rd                   |  6 ++++++
man/plot.bgms_prior_sensitivity.Rd |  5 +++--
man/plot_edge_posterior.Rd         | 17 ++++++++++++-----
```

`NAMESPACE` is unchanged. `grep -rn "PIP\|verdict printed on the panel\|node
ring\|pie channel" man/plot*.Rd` returns nothing.

### Renders

21 figures, 42 files, in `dev/review-2026-08/renders/11/`, with `INDEX.md`
giving one row per figure (what changed, which task did it) and
`make_renders.R` producing both sides. Coverage: all five plot methods at every
meaningful `type`/variant, plus the single-edge panel in its presence, absence,
undecided and no-selection cases, plus the two legend-wording variants.

**Regeneration.** Each side is one fresh `Rscript` session against its own
source tree; the fits are cached per side and the caches are gitignored:

```sh
Rscript dev/review-2026-08/renders/11/make_renders.R before ~/bgms-review/wt-fix8-base
Rscript dev/review-2026-08/renders/11/make_renders.R after  ~/bgms-review/wt-fix8
```

**Byte-stability.** Re-running the AFTER side in a fresh session reproduces
every file bit for bit:

```
$ md5 *-after.png | md5      # run 1
8197f91ee518905df9d7a861694f0dd6
$ md5 *-after.png | md5      # run 2, fresh session
8197f91ee518905df9d7a861694f0dd6
```

Every figure is seeded immediately before it is drawn, which the renders need
because qgraph's spring layout is randomised; without it a BEFORE/AFTER pair
would differ in layout as well as in style and be unreadable as a comparison.

**The fits reproduce across trees, not just across runs.** The two sides built
their caches independently, cold, in two different source trees. The resulting
figures show the same networks with the same layouts and the same edges — see
`bgms-network-{before,after}.png`, where the six nodes sit in identical
positions and the nine present and two undecided edges are the same eleven
edges. A BEFORE/AFTER pair therefore differs only in drawing code.

**F-070, shown at the default device width.** In
`sensitivity-differences-before.png` the label `intrusion-physior (pairwise)`
runs off the right-hand edge of the device and is cut mid-word. In
`sensitivity-differences-after.png` both labels sit inside a margin measured
for them.

### Files this brief did not touch

`NEWS.md`, `vignettes/`, `R/verdicts.R`, `src/`, `.github/workflows/`,
`tests/testthat/test-mcmc-diagnostics.R`, and the roxygen Details block at the
top of `R/prior_sensitivity.R` — confirmed by `git diff --stat origin/develop`.

---

## Open questions

1. **The legend word (task 2's flagged consistency question).** The panel no
   longer prints "undecided"; both network legends still do, as a line-type
   key. Two readings are defensible: a legend key names an *ink*, not a
   verdict, so "undecided" there is a caption for a dotted grey line rather
   than an assertion about an edge — or the same argument that took the word
   off the panel takes it off the legend. Both variants are rendered on the
   same figure, at the same seed and size, for a decision by eye:
   `bgms-network-{before,after}.png` versus
   `bgms-network-legend-evidence-after.png`, and the same pair on the compare
   side (`compare-difference-legend-evidence-after.png`). The variant keys the
   three lines to the Bayes factors that produce them
   (`log BF > 2.3, positive` / `|log BF| < 2.3`). It is reachable behind
   `getOption("bgms.network_legend")`, which is a review affordance: it should
   be promoted to the chosen default and the option deleted, not shipped as a
   hidden switch.

2. **Finding E** — does the ruling on verdict words extend to the
   prior-sensitivity plot's zone labels, or is that plot exempt because
   verdict movement is its subject?

3. **Finding F** — should the node wheel keep its verdict colour, or go to the
   module accent like every other wheel in the package?

4. **Finding D** — should `plot.bgms` stop erroring on an undrawable network
   and adopt `plot.bgmCompare`'s "this is a result" behaviour?

5. **Finding C** — retire, soft-deprecate, or keep
   `plot_edge_posterior(evidence_threshold =)`?

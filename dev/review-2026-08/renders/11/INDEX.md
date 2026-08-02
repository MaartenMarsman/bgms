# Brief 11 renders — every figure bgms can draw, BEFORE and AFTER

**Round 3.** BEFORE is `develop` at `a930b5a5` and has not moved. AFTER is
`fix/plot-sweep` after the maintainer's round-2 verdict and the lead's K > 2
ruling. Both sides come from the one committed script here, run once per source
tree:

```sh
Rscript dev/review-2026-08/renders/11/make_renders.R before ~/bgms-review/wt-fix8-base
Rscript dev/review-2026-08/renders/11/make_renders.R after  ~/bgms-review/wt-fix8
```

Same seeds, same fits, same device size on both sides, so a pair differs only
where the drawing code differs. Fits are cached per side under `.cache-<side>/`
(untracked); deleting the cache re-fits from the same seeds. Where the BEFORE
tree cannot draw a figure at all — a selection-off fit, a paged group display —
the pane says so in place of the drawing, which is itself the difference.

**The network display is now the documentation's edge evidence plot**: three
panels on one shared layout, node for node — the pairs the data support, the
pairs the data rule out, the pairs the data cannot decide — with each panel
titled `<what it holds>: <count>` and, under that, the rule that put them there
in raw Bayes factors (`BF > 10`, `BF < 0.1`, `0.1 < BF < 10`). **Only the first
panel is weighted**: width is the posterior mean pairwise association (the
posterior mean difference, on the compare side) and colour its sign. The other
two are drawn at uniform width, dashed for ruled out and dotted for undecided,
because there the classification is the result. A fit run **without selection**
has no inclusion Bayes factor to split by and gets one weighted panel titled
`Edge weights` / `Difference weights`.

**What round 3 changed on top of that.** Panel titles are much larger, sized by
measurement so a narrow device shrinks a title rather than clipping it. Every
sign key is gone -- the titles carry the evidence classes and their rules, and
the colour convention is in the Rd -- which also removes the key that sat on top
of the `physior` node on `compare-difference-no-selection`. The `node ring:`
explainer is gone from the main-selection figure. The bottom strip the key used
to occupy went with it, so the networks fill more of each panel. All four edge
panels lose their bottom caption line, which moves into the Rd. And a
`bgmCompare` fit on more than two groups now draws at all: the same three
panels, with the supported panel **unweighted**, because one indicator per pair
is shared across `K - 1` contrasts and no pair has a single magnitude.

**One geometry limit, unfixed and reported.** Three panels in one row want a
wide device. On R's default 7x7, each panel is about 1:3 and qgraph -- which
writes its coordinate range straight from `mar` while drawing nodes at a
physical size -- puts the outer nodes of one panel across its neighbour's
boundary. This predates round 3 and is not the legend overlap, which is fixed.
The figures here are rendered at 13.5 x 5.2, where it does not arise.

Earlier items that survive: no top-right evidence block, no in-figure captions,
no answer-sentence titles, no jargon in figure text, the `group_tag`
vacuous-parenthetical fix, the qgraph-style ring for main-effect difference
marks, connected sensitivity trajectories, and the renamed `1x` marker.

One figure kind still carries a key: the edge panel's `Posterior` / `Prior`
curve labels. Two curves on one pair of axes cannot be told apart by anything
else on the panel, which is the one case the style module's no-keys rule
exempts.

## `plot.bgms`

| figure | what changed |
| --- | --- |
| `bgms-network` | The three-panel edge evidence plot replaces the single verdict-encoded network. |
| `bgms-network-all-absence` | Independent data: nothing is supported and the absence panel is the one that fills. On BEFORE this fit is an error (nothing drawable); it is now a figure, which is the honest picture of it. |
| `bgms-network-sparse` | `evidence_threshold = 1000`: the rule lines read `BF > 1000` / `BF < 0.001`, and the absence panel is the empty one. |
| `bgms-network-no-selection` | `edge_selection = FALSE`: one weighted panel titled `Edge weights`. BEFORE could not draw this fit at all. |
| `bgms-centrality` | `type = "centrality"`; the bottom caption is gone, the axis carries the quantity. |

## `plot.bgmCompare`

| figure | what changed |
| --- | --- |
| `compare-difference` | The same three panels, read for differences: `difference supported` / `difference ruled out` / `undecided`. |
| `compare-difference-main-selection` | `main_difference_selection = TRUE`. The beside-node probability wheel is gone; the ring is back around the node circle, qgraph-style, filled to `P(main-effect difference)`. |
| `compare-difference-empty` | `evidence_threshold = 1e6`: everything lands undecided, and the figure says so by which panel is full. |
| `compare-difference-no-selection` | `difference_selection = FALSE`: one weighted panel titled `Difference weights`. BEFORE could not draw this fit. The sign key that overlapped the `physior` node in round 2 is gone. |
| `compare-difference-k3` | Three groups, so two contrasts under one indicator per pair. The three panels are unchanged and the supported panel is **unweighted** -- uniform width, plain ink, no sign colour -- because no pair has a single magnitude. The fixture is built rather than sliced: groups 1 and 2 share a chain, group 3 trades one link for another, and the figure reports exactly those three pairs as supported. BEFORE refuses the fit outright, so its pane says so. |
| `compare-groups` | Group networks only, on the difference display's shared layout, with a character-labelled fit so `group_tag` shows its informative form. The difference panel no longer rides along; it is the default display. |
| `compare-groups-paged-1`, `-2` | Four groups, `max_panels = 3`: page 1 holds three networks, page 2 the fourth at the same panel width. BEFORE has no pager and errors. |
| `compare-centrality` | `type = "centrality"` for one group; caption gone. |

## `plot.bgms_centrality`

| figure | what changed |
| --- | --- |
| `centrality-panel` | Module type, ink and offset axis; measured left margin; no caption. |

## `plot.bgms_calibration`

| figure | what changed |
| --- | --- |
| `calibration-isotonic` | Module style and scaling; no bottom caption; the `isotonic` kind-mark is gone, and the outer axis labels no longer name the estimators. |
| `calibration-mixed` | Both variable kinds on one page, so each panel says which it is — `discrete` / `continuous`, not `isotonic` / `PIT`. |
| `calibration-groups` | A `bgmCompare` check, paged; same changes. |
| `calibration-single` | The one-panel case. |

## `plot.bgms_prior_sensitivity`

| figure | what changed |
| --- | --- |
| `sensitivity-edges` | Headline sentence gone; module type and offset axes; the `1x` marker now says `your fit` and the axis label agrees; one connected trajectory per edge through the anchor dots; off-range note removed (it is in `print()`). |
| `sensitivity-differences` | The same on a compare check. This is also the F-070 pair: on BEFORE the edge names run off the right of the device, on AFTER a measured margin holds them. |

## `plot_edge_posterior`

| figure | what changed |
| --- | --- |
| `edge-panel-presence` | The lock is lifted by maintainer instruction. The only change from round 2 is the removed bottom caption; every other pixel is unmoved. |
| `edge-panel-undecided` | As above: caption removed, nothing else. |
| `edge-panel-absence` | As above: caption removed, nothing else. |
| `edge-panel-no-selection` | The `data|H1` / `data|H0` tags went in round 2; the two caption lines go now. What each mark means is in the Rd. |

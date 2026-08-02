# Brief 11 renders — every figure bgms can draw, BEFORE and AFTER

**Round 2.** BEFORE is `develop` at `a930b5a5`. AFTER is `fix/plot-sweep` after
the maintainer's render judgment. Both sides come from the one committed script
here, run once per source tree:

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

Round-1 items that survive: no top-right evidence block, no in-figure captions,
no answer-sentence titles, no jargon in figure text, the `group_tag`
vacuous-parenthetical fix, the qgraph-style ring for main-effect difference
marks, connected sensitivity trajectories, and the renamed `1x` marker.

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
| `compare-difference-no-selection` | `difference_selection = FALSE`: one weighted panel titled `Difference weights`. BEFORE could not draw this fit. |
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
| `edge-panel-presence` | **Locked.** Byte-identical to round 1. |
| `edge-panel-undecided` | **Locked.** Byte-identical to round 1. |
| `edge-panel-absence` | **Locked.** Byte-identical to round 1. |
| `edge-panel-no-selection` | The `data|H1` / `data|H0` tags are gone and the wheel corner is de-crowded: the wheel sits where it does on the other three panels and the caption says in words what its two shares are. |

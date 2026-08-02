# Brief 11 renders — every figure bgms can draw, BEFORE and AFTER

BEFORE is `develop` at `a930b5a5` (contains the brief-12 merge `f24ad3c8`).
AFTER is `fix/plot-sweep`. Both sides are produced by the one committed script
in this directory:

```sh
Rscript dev/review-2026-08/renders/11/make_renders.R before ~/bgms-review/wt-fix8-base
Rscript dev/review-2026-08/renders/11/make_renders.R after  ~/bgms-review/wt-fix8
```

Same seeds, same fits, same device size on both sides, so a pair differs only
where the drawing code differs. The fits are cached per side under
`.cache-<side>/` (untracked); deleting the cache re-fits from the same seeds.

Task numbers are the brief's: **1** = F-067 package-wide restyle, **2** =
verdict word off the edge panel, **3** = "PIP" → `P(included)`, **4** = three
decimals, **5** = F-079 compare evidence parity, **6** = qgraph pie → wheel,
**7** = F-070 sensitivity label margin.

## `plot.bgms`

| figure | what changed | tasks |
| --- | --- | --- |
| `bgms-network` | Evidence band added above the network: display name, verdict tally, threshold and the strongest log BF for and against an edge. Legend and caption moved to the module's type, colour and placement; nodes drawn in module ink on white. | 1, 5 |
| `bgms-network-legend-evidence` | The same figure with the legend keyed to the Bayes factors (`log BF > 2.3, positive`) instead of the verdict words, for the legend-wording decision in the report's open questions. BEFORE has no such variant and shows the ordinary figure. | 2 (question) |
| `bgms-network-sparse` | `evidence_threshold = 1000`: most edges ruled out. Shows the tally and the band on a nearly empty network. | 1, 5 |
| `bgms-centrality` | `type = "centrality"`, which delegates to `plot.bgms_centrality`. | 1 |

## `plot.bgmCompare`

| figure | what changed | tasks |
| --- | --- | --- |
| `compare-difference` | Same evidence band as `plot.bgms`, composed by the same code: difference tally, threshold in natural-log Bayes factors, strongest difference Bayes factor for and against. Restyled legend, caption and nodes. | 1, 5 |
| `compare-difference-legend-evidence` | The legend-wording variant on the compare side, so the choice can be judged on both figures. | 2 (question) |
| `compare-difference-main-selection` | `main_difference_selection = TRUE`. The per-node qgraph `pie` ring is gone; each node carries the package's own `probability_wheel()` beside it, filled to `P(main-effect difference)`. | 1, 5, 6 |
| `compare-difference-empty` | `evidence_threshold = 1e6`: nothing reaches presence. The note under the tally says so; the picture is still drawn. | 1, 5 |
| `compare-groups` | `type = "groups"`. Each group panel names itself with the module's compact label instead of a qgraph `title` banner; the difference panel carries the shared evidence band. Type scaled down as one for the three-panel layout. | 1, 5 |
| `compare-centrality` | `type = "centrality"` for one group — the group-labelled variant of the centrality panel. | 1 |

## `plot.bgms_centrality`

| figure | what changed | tasks |
| --- | --- | --- |
| `centrality-panel` | Module type, ink and offset axis. The left margin is measured from the node names actually being drawn. The axis carries the quantity; what the dot and the bar are moved into a muted caption. | 1 |

## `plot.bgms_calibration`

| figure | what changed | tasks |
| --- | --- | --- |
| `calibration-isotonic` | Six isotonic panels. Per-panel `main =` replaced by the module's compact label with the panel kind as its subtitle; offset axes on the unit square; module colours, band alpha and line weights, scaled as one for the grid. Outer caption added. | 1 |
| `calibration-mixed` | A mixed fit: PIT panels for the continuous variables, isotonic for the discrete. Same changes. | 1 |
| `calibration-groups` | A `bgmCompare` check, one panel per variable per group, paged. Same changes. | 1 |
| `calibration-single` | One variable — the single-panel case, where the grid scaling has to degrade gracefully. | 1 |

## `plot.bgms_prior_sensitivity`

| figure | what changed | tasks |
| --- | --- | --- |
| `sensitivity-edges` | The headline `main =` title is retired; the answer ("No edge verdict depends on the scale") is now the module's compact label. Module type, ink, offset axes on both scales. The right margin is measured from the widest edge name, and the leader and label offsets are device lengths rather than fractions of the data range. | 1, 7 |
| `sensitivity-differences` | The same on a `bgmCompare` check, where labelled movers exist — this is the pair that shows F-070: on BEFORE the names run into the device edge at the default width, on AFTER the margin holds them. | 1, 7 |

## `plot_edge_posterior` (already on the module; the ratified changes only)

| figure | what changed | tasks |
| --- | --- | --- |
| `edge-panel-presence` | The verdict word ("evidence of presence") is gone from under the edge name. `PIP` is now `P(included)`, printed to three decimals. | 2, 3, 4 |
| `edge-panel-undecided` | As above; the "undecided" annotation is gone. | 2, 3, 4 |
| `edge-panel-absence` | As above; the "evidence of absence" annotation is gone. | 2, 3, 4 |
| `edge-panel-no-selection` | The Savage-Dickey variant. `"no edge selection"` stays — it names the model, not a verdict — and the wheel keeps its `data|H1` / `data|H0` tags. | 2, 4 |

# Report 26b — reframe the vignette set under the maintainer's terminology rulings

Brief 26b, the remediation round of brief 26. Branch `fix/intro-vignette-2`,
worktree `~/bgms-review/wt-fix17`, based on `origin/develop` at `00ea851b`
(at or after `d9b5a4dd`, as the brief required). Four commits — **not pushed**.
`origin/develop` has since moved to `d05c03b0` (report-25 work, another agent);
nothing of mine touches it and nothing under `dev/review-2026-08/reports/25*`
was opened in any tree. The Dropbox tree's branch was never switched and
nothing was built in it.

**Surface touched:** `vignettes/intro.Rmd`, `vignettes/checking-your-model.Rmd`,
`vignettes/comparison.Rmd`, `vignettes/diagnostics.Rmd`, plus this report and
its assets. Nothing under `R/`, `src/`, `tests/`, `man/`, `NEWS.md`,
`vignettes/refs.bib`, `vignettes/prior-sensitivity.Rmd` or the workflows was
modified. Every citation the rewrites add was already a key in `refs.bib`, so
the shared bibliography is byte-unchanged.

```
 vignettes/checking-your-model.Rmd |  26 ++-
 vignettes/comparison.Rmd          |  53 +++---
 vignettes/diagnostics.Rmd         |  10 +-
 vignettes/intro.Rmd               | 375 +++++++++++++++++++++++++++-----------
 4 files changed, 317 insertions(+), 147 deletions(-)
```

**Headline.** The four rulings are now carried by every sentence in the
vignette set. `intro.Rmd` opens with the maintainer-approved text, byte-for-byte
apart from the file's `--` and `Blume--Capel` typography; the regularization
contrast is gone with no replacement; the models section presents one Markov
random field family and its members rather than a discrete/continuous split; and
"spike-and-slab" appears once, under *Priors*. **The word "network" does not
appear in any of the five vignette sources** — the count is zero, not merely
zero for "network model", so the retained-instance table in §3 is empty of
prose. Two places still put "network" in front of a reader's eye and neither is
mine to fix: three published paper titles in the bibliographies, and one string
the plot code itself draws into the `comparison.Rmd` figure (§7.1 — this one is
release-visible and I recommend the lead act on it before merge).

F-127 is fixed and measured: the narrow-device advisory is **gone** from the
released `checking-your-model.html` and the three panel titles no longer
collide (§4, assets). F-128 is fixed: `comparison.Rmd` now calls
`plot(fit, type = "groups")` and carries the bibliography header.

Gates: `R CMD check --as-cran` at the 2 baseline NOTEs; the default-tier suite
at 8518 pass / 0 fail / 0 error / 0 warning, identical to round 1 — nothing
moved and nothing should have; intro knit **16.5 s** against a ~30 s budget;
full five-vignette rebuild inside the check **`[94s/50s]`**.

---

## 1. What was done

| # | Commit | What | Finding |
|---|--------|------|---------|
| 1 | `ca6c33b2` | Reframe `vignettes/intro.Rmd` | F-026 |
| 2 | `18274dc9` | `checking-your-model.Rmd`: plot chunk width, three-panel prose, model naming | F-127 |
| 3 | `52620f60` | `comparison.Rmd`: `plot(fit, type = "groups")`, bibliography header, terminology | F-128 |
| 4 | `1e5ad7c6` | `diagnostics.Rmd`: terminology sweep | F-026 |
| 5 | (below) | This report + six assets | — |

The intro is one commit because it is one file rewritten under one set of
rulings; splitting the frame from the terminology would produce an intermediate
state that satisfies neither.

---

## 2. Per-file change map

### 2.1 `vignettes/intro.Rmd` (F-026)

The round-1 machinery survives intact: the nine-item Wenchuan worked example,
the show-uncapped / run-capped chunk pair with `cores = 2`, the five-row
defaults table, the `summary()` → `verdicts()` → `plot()` flow, the
`fig.width = 12, fig.height = 4.5, out.width = "100%"` plot chunk, the
GGM/mixed and `bgmCompare()` sections, the bibliography header, and the
website / `easybgm` / `NEWS.md` pointers. The frame and the vocabulary are new.

| Round-1 element | Disposition under the rulings |
|---|---|
| `# What the package answers` header | **Dropped.** The approved opening is lead text under the vignette's own title; inventing a header above it would have put a frame in front of the frame. |
| Para 1: "estimates network models -- Markov random fields for discrete data, Gaussian graphical models for continuous data, and mixed networks that hold both" | **Replaced verbatim** by the approved opening's first paragraph. This was the sentence that carried both ruling-1 and ruling-2 violations at once. |
| Para 2: "A penalized estimate returns one network, and an absent edge in it can mean either of two things…" | **Deleted, not rewritten** (ruling 3). Nothing replaces the regularization contrast; the approved opening's second paragraph states what the analysis returns and stops there. |
| Para 2: "puts a spike-and-slab prior on every edge" | **Deleted** (ruling 4). Spike-and-slab now appears once, opening *Priors*. |
| Scope contract + website pointer | **Kept**, and now carries `@HuthEtAl_2023` / `@SekulovskiEtAl_2023`, which the old paragraph 2 carried. This keeps the approved opening free of any insertion of mine. |
| `# The models`: "fits a single network", "where their networks differ" | → "fits one Markov random field to one sample", "where their graphs differ". |
| `# The models`: the four `variable_type` bullets | **Kept as four bullets, reframed as members.** "Which member `bgm()` fits follows from `variable_type`", "fits the **ordinal Markov random field**", "fits the **Gaussian graphical model**", "mixes the members". The GGM bullet no longer says "the network is the precision matrix". |
| `# The models` closing sentence | Was "Everything downstream … works the same way whichever of these is fitted." Now leads with the ruling: "These are members of one family rather than separate modelling frameworks, and everything downstream … is the same whichever member is fitted." |
| *Which edges the evidence settles*: "the point of running an evidence analysis rather than an estimation one" | → "the point of treating the graph as unknown". The old phrasing was the regularization contrast in miniature. |
| *The edge evidence plot*: "A single network drawing would have to collapse…" | → "A single drawing of the graph would have to collapse…". |
| *The edge evidence plot*: "Three networks side by side want a wide device" | → "Three panels side by side want a wide device". |
| *The edge evidence plot*: "`plot_edge_posterior()` … showing its spike-and-slab posterior against the prior it was updated from" | → "showing how its posterior mass divides between the edge being absent and the edge being present". This was the second spike-and-slab mention. |
| `# Priors` opening | **New sentence**, the single licensed naming: "Edge selection rests on a spike-and-slab prior for each pairwise interaction: a point mass at zero for the edge being absent, and a slab for the values it can take when the edge is present." |
| `# Continuous and mixed data` header | → `# The Gaussian and mixed members`. The old header read as a second kind of model; the new one names them as members of the family the opening established. |
| That section: "for all three model classes" | → "for all members of the family". |
| `# Comparing groups`: unchanged in substance | No violations found; left as round 1 wrote it. |
| `# Where to go next` | Unchanged. |

### 2.2 `vignettes/checking-your-model.Rmd` (F-127)

| Line (pre) | Change |
|---|---|
| 24 | "A fitted **network** raises questions of two kinds." → "A fitted **graphical model** raises questions of two kinds." |
| 80 | Chunk options `fig.height = 5.5, fig.width = 6` → **`fig.width = 12, fig.height = 4.5, out.width = "100%"`** — brief task 2(a). Clears the advisory and the title collision; see §4. |
| 84–87 | **Rewritten.** The old paragraph described the pre-0.2.0.0 single-panel display ("dotted grey edges are undecided; edges with evidence of absence are not drawn at all"), which the figure above it has contradicted since 0.2.0.0. It now describes the three-panel display: one panel per verdict, shared layout, only the first weighted, dashed for ruled-out and dotted for undecided, and why a single drawing cannot say both kinds of blank. |
| 93 | "Centrality is a function of the **network**" → "of the **graph**". |
| 205 | "which is what a **network model** can miss" → "which is what a **graphical model** can miss" — brief task 2(c). |
| 227 | "whether to trust the **network** as a description of the joint distribution" → "whether to trust the **model** as …". |

Reading the whole file for the same two diseases turned up nothing else: no
other stale display description, and no other model naming. The remaining
"spike-and-slab posterior" at line ~86 is left standing — it introduces
`plot_edge_posterior()`, whose subject *is* the split of the mass, and ruling 4
governs what the front door leads with rather than banning the term where
priors are the topic. Flagged here so the maintainer can overrule.

### 2.3 `vignettes/comparison.Rmd` (F-128)

| What | Change |
|---|---|
| YAML | **Added** `bibliography: refs.bib`, `csl: apa.csl`, `link-citations: TRUE`, and a `# References` section at the foot. Three keys cited, all pre-existing: `@MarsmanHaslbeck_2023_ordinal`, `@HuthEtAl_2023`, `@SekulovskiEtAl_2023`. `refs.bib` is byte-unchanged. |
| Intro ¶1 | "estimates whether **edge weights** and category thresholds differ across groups in an ordinal Markov random field (MRF)" → "estimates whether **pairwise interactions** and category thresholds differ across groups in an ordinal Markov random field […], **the member of the family that `bgm()` fits to binary and ordinal data**". The parenthetical "(MRF)" is dropped — the abbreviation is never used again in the file. |
| Intro ¶2 | Was "Posterior inclusion probabilities indicate how plausible it is … These can be converted to Bayes factors for hypothesis testing." Now states the inclusion Bayes factor directly and names `verdicts()` and its three answers, matching the intro's vocabulary. |
| `# Visualizing group networks` | → `# Visualizing the groups`. |
| The `qgraph` block | **Replaced.** The hand-rolled 14-line block (`library(qgraph)`, an `adhd_network` matrix assembled from `coef(fit)$pairwise_effects_groups[, 1]`, `theme = "TeamFortress"`) drew group 1 only, by hand, at `fig.width = 7, fig.height = 7`. It is now `plot(fit, type = "groups")` at `fig.width = 10, fig.height = 5, out.width = "100%"` — both groups, shared layout, the package's own encoding. Rendered figure at `assets/f026b-comparison-groups.png`. |
| Prose around it | New: what the default display is, what `type = "groups"` gives, that qgraph is Suggested rather than required, and pointers to `plot(fit)` and `plot(fit, type = "centrality")` with `?plot.bgmCompare`. |

`plot(fit, type = "groups")` does not call `warn_narrow_device()` —
`R/plot_bgms.R` calls it only from `draw_evidence_panels()`, and
`?plot.bgmCompare` says so explicitly ("The single-panel and `type = "groups"`
displays are content with the default device"). Verified in the render: no
message.

### 2.4 `vignettes/diagnostics.Rmd` (F-026, terminology only)

Five instances, all prose, no chunk or output touched:

| Line (pre) | Change |
|---|---|
| 179 | "strong evidence that there is no **network relation** between two variables" → "…that there is no **edge** between two variables" |
| 235 | "the recovered **network** reflects the shortcut" → "the recovered **graph** reflects the shortcut" |
| 252 | "consequences for the recovered **network**" → "…for the recovered **graph**" |
| 255 | "leaning the recovered **network** in one direction" → "…the recovered **graph**…" |
| 256 | "what moves a recovered **network**" → "what moves a recovered **graph**" |

### 2.5 `vignettes/prior-sensitivity.Rmd`

Swept, zero instances of "network", zero model-naming violations. **Not
touched.**

---

## 3. The retained-"network" table

The brief asks for every retained instance across the five vignettes with its
sentence, marked keep/changed. **There are none.** `grep -in "network"` over
`vignettes/*.Rmd` returns nothing:

```
$ grep -ci "network model" vignettes/*.Rmd
vignettes/checking-your-model.Rmd:0
vignettes/comparison.Rmd:0
vignettes/diagnostics.Rmd:0
vignettes/intro.Rmd:0
vignettes/prior-sensitivity.Rmd:0

$ grep -in "network" vignettes/*.Rmd
(no output)
```

The gate the brief sets — "network model" count 0 — is met with margin: the bare
noun is gone too. No vignette uses `plot.bgms`'s own `type = "network"`
argument, so the exemption the brief allowed was never needed.

Because the maintainer judges the **render**, the honest accounting is what a
reader still sees, not what the source still says. Two channels remain, neither
of them vignette prose:

| Where | Text | Status | Why |
|---|---|---|---|
| `intro.html`, `comparison.html` — reference lists | "Bayesian analysis of cross-sectional **networks**: A tutorial in R and JASP" (Huth et al., 2023) | **keep** | Published title, rendered from `refs.bib`, which is off-limits and correct as it stands. |
| `intro.html`, `comparison.html` — reference lists | "Testing conditional independence in psychometric **networks**" (Sekulovski et al., 2024) | **keep** | Same. |
| `intro.html` — reference list | "A **network** approach to posttraumatic stress disorder" (McNally et al., 2015) | **keep** | Same. |
| `comparison.html` — **inside the figure** | Panel subtitle "posterior mean **network**", drawn twice (once per group) | **cannot be changed from my surface** | `R/plot_bgms.R:974`. See §7.1 — this is the one I would ask the lead to act on before merge. |

Rendered-text scan of all five built HTML files, confirming the count above:

```
checking-your-model.html | network mentions in rendered text: 0
comparison.html          | network mentions in rendered text: 2   (both bibliography titles)
diagnostics.html         | network mentions in rendered text: 0
intro.html               | network mentions in rendered text: 3   (all bibliography titles)
prior-sensitivity.html   | network mentions in rendered text: 0
```

---

## 4. F-127 — the plot chunk output, verbatim

Round 1 quoted what the **shipped** `inst/doc/checking-your-model.html`
contained: the sampler's own advisory printed as vignette output, and a figure
whose three titles ran together as
`evidence of presence: 11evidence of absence: 10undecided: 15`.

Here is the same chunk in the HTML this branch built, taken raw from
`bgms.Rcheck/bgms/doc/checking-your-model.html` (base64 image payload elided):

```html
drift apart.</p>
<div class="sourceCode" id="cb4"><pre class="sourceCode r"><code class="sourceCode r"><span id="cb4-1"><a href="#cb4-1" tabindex="-1"></a><span class="fu">plot</span>(fit)</span></code></pre></div>
<p><img role="img" aria-label src="data:image/png;base64,<BASE64 PNG, elided>
```

The source line `plot(fit)` is followed directly by the `<img>`. There is no
output block, no `#>` prefix, no message — the code element closes and the
figure begins. As rendered text:

```
The default picture encodes the same three verdicts, so what a reader
sees and what the table says cannot drift apart.
plot(fit)
[IMG]
There is one panel per verdict – the pairs the data support, the
pairs the data rule out, and the pairs the data cannot decide – ...
```

Programmatic confirmation on the built file:

```
ADVISORY PRESENT: False        # 'three-panel display is drawn on a device' not in file
```

The figure is at `assets/f026b-cym-plot-clean.png` — the direct counterpart of
round 1's `assets/f026-cym-plot-crowded.png`. Titles read cleanly and
separately: `evidence of presence: 11` / `evidence of absence: 10` /
`undecided: 15`, each with its rule (`BF > 10`, `BF < 0.1`, `0.1 < BF < 10`)
beneath it.

---

## 5. Times

| Measurement | Value | Budget |
|---|---|---|
| `intro.Rmd` knit, standalone, installed `-O2` build | **16.5 s** | ~30 s |
| `checking-your-model.Rmd` knit, standalone | 13.4 s | — |
| `comparison.Rmd` knit, standalone | 4.4 s | — |
| Full five-vignette rebuild inside `R CMD check` | **`[94s/50s]`** | — |
| `R CMD check --as-cran`, tests block | `[110s/88s]` | — |
| `R CMD check --as-cran`, examples with `--run-donttest` | `[484s/248s]` | — |

```
KNIT SECONDS: 16.5
Rscript -e  2>&1  31.28s user 0.32s system 190% cpu 16.602 total
```

Round 1 measured 16.6 s on the same fit; the reframe is prose, so the knit time
is unchanged within noise. The build ran with no other R process on the machine
(`ps aux | grep "[e]xec/R"` empty before each stage), one job at a time,
`MAKEFLAGS="-j4"`.

The three standalone renders are **byte-identical** to the files `R CMD build`
put in `inst/doc/`:

```
intro: IDENTICAL
comparison: IDENTICAL
checking-your-model: IDENTICAL
```

so the assets in §9 are the shipped artifacts, not a separate render.

---

## 6. NEWS proposals — VERBATIM

The lead lands all of these. `NEWS.md` was not opened on this branch.

### 6.1 Tag status, established

```
$ git ls-tree --name-only cran-0.1.6.3 vignettes/
vignettes/.gitignore
vignettes/apa.csl
vignettes/comparison.Rmd
vignettes/diagnostics.Rmd
vignettes/intro.Rmd
vignettes/refs.bib
```

| Vignette | At `cran-0.1.6.3`? | NEWS needed? |
|---|---|---|
| `intro.Rmd` | yes | **yes** — §6.2 |
| `comparison.Rmd` | yes | **yes** — §6.3 |
| `diagnostics.Rmd` | yes | **no** — §6.4 |
| `checking-your-model.Rmd` | **no** | **no** — §6.5 |
| `prior-sensitivity.Rmd` | no | untouched |

### 6.2 `intro.Rmd` — under `## Other changes` (0.2.0.0)

Round 1's proposed line, rewritten under the rulings. It said the vignette
"routes the four model classes through `variable_type`" (ruling 1) and led with
prior machinery (ruling 4); both are corrected, and the note about what was
dropped now names the regularization framing as well.

```markdown
* The "Getting Started with bgms" vignette has been rewritten against this
  release. It opens on what the package is for -- Bayesian analysis of
  graphical models, every one of them a Markov random field, with
  `variable_type` selecting the member -- and on what an analysis returns: a
  posterior over graphs and an inclusion Bayes factor per pair, read as
  evidence of presence, evidence of absence, or undecided. It states the
  current defaults (`normal_prior(scale = 1)` for the interaction slab,
  `bernoulli_prior(0.5)` for the edge prior, `"hierarchical"` for
  `precision_graph_prior`) and works through `summary()`, `verdicts()` on the
  natural log Bayes-factor scale, and `plot()`'s three-panel edge evidence
  display on nine Wenchuan items. It gained the bibliography header the other
  vignettes carry, a `bgmCompare()` section, and pointers to the package
  website and to `easybgm`. The hand-thresholded median-probability graph drawn
  with `qgraph`, the `coef()` walkthrough, and the "What's new in 0.2.0"
  section were dropped.
```

### 6.3 `comparison.Rmd` — under `## Other changes` (0.2.0.0)

Tag-visible, and the change is user-facing: the display it teaches is different
code.

```markdown
* The "Model Comparison with bgmCompare" vignette now draws the groups with
  `plot(fit, type = "groups")` -- each group's own graph on one shared layout --
  in place of the hand-assembled `qgraph()` call it built from
  `coef(fit)$pairwise_effects_groups`, which drew one group only. It gained the
  bibliography header the other vignettes carry, and its introduction now states
  the inclusion Bayes factor and the three verdicts `verdicts()` returns.
```

### 6.4 `diagnostics.Rmd` — no entry proposed

Five words of prose ("network" → "graph"/"edge"), no code, no output, no
behaviour. My recommendation is no NEWS line; if the lead wants the terminology
pass recorded once for the whole set rather than per file, this is the wording
I would use, again under `## Other changes`:

```markdown
* Vignette prose throughout now calls the models graphical models and the
  estimated structure the graph, matching the package documentation and the
  `DESCRIPTION`.
```

Landing that line would make §6.2's and §6.3's mentions of terminology
redundant but not wrong; I have not tried to pre-merge them, since which of the
two shapes the lead wants is his call.

### 6.5 `checking-your-model.Rmd` — no entry proposed

The vignette is new in 0.2.0.0 (absent from `cran-0.1.6.3`) and already has its
entry at `NEWS.md:438`, *"New vignette 'Checking your fitted model'…"*. F-127
corrects a defect in unreleased material; a reader upgrading from 0.1.6.3 sees
the vignette for the first time, correct. Nothing to announce.

---

## 7. Findings

### 7.1 The plot code draws "posterior mean network" into the comparison vignette — release-visible

`R/plot_bgms.R:974`, in the `type = "groups"` panel loop:

```r
    draw_network_panel(
      m, variables, shared,
      title = group_tag(labels, g),
      rule = "posterior mean network",
```

That string is rendered as the subtitle under each group's panel title, so the
figure `comparison.Rmd` now draws reads:

```
group 1 (x)                    group 2 (y)
posterior mean network         posterior mean network
```

See `assets/f026b-comparison-groups.png`. This is the exact vocabulary the
maintainer rejected, in a figure he will look at, and F-128's fix is what puts
it there — the hand-rolled `qgraph()` block it replaces drew no subtitle. The
string is in `R/`, outside my surface, so I did not change it.
`"posterior mean graph"` is a one-word edit with no API consequence.

The full audit of user-visible "network" strings in `R/` (excluding code
comments and internal function names), for the lead:

| Site | String | Reaches a user how |
|---|---|---|
| `R/plot_bgms.R:974` | `rule = "posterior mean network"` | **Drawn into every `type = "groups"` figure**, including `comparison.Rmd`'s |
| `R/plot_bgms.R:403` | `"…three networks side by side want at least "` | The narrow-device advisory message |
| `R/plot_bgms.R:591` | `"Drawing the network needs the qgraph package…"` | Error when qgraph is absent |
| `R/plot_bgms.R:1103–1105` | `"plot_edge_posterior() draws an edge of a single network, and a bgmCompare() fit parameterizes differences between networks rather than one edge weight. Use plot(fit) for the difference network…"` | Error on `plot_edge_posterior()` of a `bgmCompare` fit |
| `R/verdicts.R:534` | `"operating point was established on single-network edge indicators only."` | `verdicts()` message |
| `R/plot_bgms.R:540` | `type = c("network", "centrality")` | `plot.bgms()`'s argument — the brief exempts the API name |

Only the first is inside a rendered vignette. The other four are messages and
errors; whether the rulings reach them is the maintainer's call, not something
I would assume.

### 7.2 The approved opening is verbatim, with one structural decision recorded

The two approved paragraphs are reproduced with no word changed. The file's
typography (`--` for the em dash, `Blume--Capel`) matches the rest of the
vignette set, which the brief licensed.

The decision I made rather than the brief: **the citations that round 1 hung on
the inclusion-Bayes-factor sentence are not inserted into the approved text.**
`[@HuthEtAl_2023; @SekulovskiEtAl_2023]` would have sat naturally in the second
paragraph, but any insertion is a change to text the maintainer approved. They
moved to the third paragraph instead: *"The inclusion Bayes factor and the
three-state reading it supports are developed in @HuthEtAl_2023 and
@SekulovskiEtAl_2023."* If the maintainer would rather have them inline, the
edit is two brackets.

### 7.3 "Spike-and-Slab" survives once more than ruling 4 licenses — as a title

`# Where to go next` points at *"Diagnostics and Spike-and-Slab Summaries"*,
which is that vignette's actual `title:` and `\VignetteIndexEntry`. A pointer
that renames the thing it points at is a worse defect than the repetition, so I
kept it accurate. If the maintainer wants the phrase gone from the front door
entirely, the fix is to rename the diagnostics vignette — a `man/`- and
`NEWS`-visible change, outside this brief.

### 7.4 `comparison.Rmd` still calls `?ADHD` in an evaluated chunk

Line 36, unchanged from before this branch. In a knitted vignette it prints the
whole help page to the build console and contributes nothing to the HTML — the
rendered document shows only the source line. It is neither a terminology nor a
display defect, so it stayed outside my three tasks; recording it because it is
one line to delete and the vignette reads better without it.

### 7.5 The round-1 finding that the intro alone lacked a bibliography was half right

Round 1 noted the brief's claim that intro "is the only one without" a
bibliography header and pointed out that `comparison.Rmd` had none either.
F-128 closes that: **all five** vignettes now carry `bibliography: refs.bib`,
`csl: apa.csl`, `link-citations: TRUE`. What they do with it is not uniform,
and I did not make it so:

| Vignette | `# References` section | Keys cited |
|---|---|---|
| `intro.Rmd` | yes | 5 |
| `comparison.Rmd` | yes (added) | 3 |
| `prior-sensitivity.Rmd` | yes | **0 — the section renders empty** |
| `diagnostics.Rmd` | no | 1 (`@VehtariEtAl_2021`; citeproc appends it at the foot) |
| `checking-your-model.Rmd` | no | 0 |

`prior-sensitivity.Rmd`'s empty *References* heading is pre-existing, cosmetic,
and in a file with no terminology violation, so I left it. One line to delete
if the lead wants the set tidy before CRAN.

---

## 8. Verification gate

| # | Gate | Result |
|---|---|---|
| 1 | `R CMD check --as-cran` on a `git archive` tarball — baseline NOTEs only | **PASS** — `Status: 2 NOTEs`, both baseline. §8.1 |
| 2 | Full default-tier suite — 0 failures / 0 warnings | **PASS** — 8518 pass / 0 fail / 0 error / 0 warning. §8.2 |
| 3 | Terminology gates — "network model" 0; retained table complete; approved opening verbatim | **PASS** — §3, §7.2 |
| 4 | Every touched vignette's HTML read end to end; F-127 advisory gone; times measured | **PASS** — §4, §5, §8.3 |

### 8.1 `R CMD check --as-cran`

Tarball built by `R CMD build` from a `git archive` of `HEAD` (`1e5ad7c6`) into
a clean tree, pandoc on `PATH`, `_R_CHECK_FORCE_SUGGESTS_=true`.

```
Status: 2 NOTEs

* checking CRAN incoming feasibility ... [3s/12s] NOTE
Maintainer: 'Maarten Marsman <m.marsman@uva.nl>'
The Date field is over a month old.

* checking HTML version of manual ... NOTE
Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy.
```

Both are the known baselines — F-003 (stale `Date`) and report 01's finding 10
(local HTML Tidy). No ERROR, no WARNING, no third NOTE. The examples-timing NOTE
the brief warns about did not appear
(`checking examples with --run-donttest ... [484s/248s] OK`). Tests
`[110s/88s] OK`; `checking package vignettes ... OK`; `checking for unstated
dependencies in vignettes ... OK`; `re-building of vignette outputs
... [94s/50s] OK`.

### 8.2 Default-tier suite

`testthat::test_check("bgms")` from `tests/` against the installed branch build
(`~/bgms-review/lib17`, `-O2`), `NOT_CRAN=true`, `BGMS_RUN_SLOW_TESTS` and
`BGMS_RUN_CERTIFICATION` unset.

```
FILES: 79
TESTS: 1225
PASS:  8518
FAIL:  0
ERROR: 0
WARN:  0
SKIP:  97
elapsed: 223 s
```

All 97 skips are tier gates. The numbers are identical to round 1's, which is
the point: the branch changes four `.Rmd` files under `vignettes/` and no code
path the suite exercises. **Nothing moved, and nothing should have.**

(For the record: a first attempt used `test_dir("testthat")` with only
`library(bgms)` attached, which cannot see internal functions and reported 132
errors of the form `could not find function "anchor_log_weights"`. That is a
harness mistake of mine, not a branch state — `test_check()` is what
`tests/testthat.R` and `R CMD check` use, and it is the run recorded above.)

### 8.3 The rendered vignettes

All four touched vignettes read end to end in the built HTML.

- **`intro.html`** — every evaluated chunk ran and produced output; nothing
  warned, errored, or leaked a message. Evaluated: the hidden fit,
  `summary(fit)`, `verdicts(fit)`, `plot(fit)`. The `eval = FALSE` chunks (the
  shown `bgm()` call, the prior example, the GGM example, the `bgmCompare()`
  example) are by design and were each re-checked against current formals. One
  figure, one table, five references rendered. The approved opening renders as
  written.
- **`checking-your-model.html`** — five figures, all chunks clean, **advisory
  absent** (§4), the three panel titles separate.
- **`comparison.html`** — one figure (`type = "groups"`, both groups drawn),
  reference list rendered from the new header, no message from `plot()`, no
  qgraph attach in the source.
- **`diagnostics.html`** — unchanged except for the five words; rendered clean,
  zero "network" in the rendered text.
- **`prior-sensitivity.html`** — untouched, rebuilt clean by the check.

---

## 9. Assets

Copied from the built tarball's `inst/doc/`
(`bgms.Rcheck/bgms/doc/`), not separately rendered:

| File | What |
|---|---|
| `assets/f026b-intro.html` | 155,121 B — the reframed front door |
| `assets/f026b-checking-your-model.html` | 334,735 B — F-127 fixed |
| `assets/f026b-comparison.html` | 138,388 B — F-128 fixed |
| `assets/f026b-diagnostics.html` | 236,974 B — terminology sweep |
| `assets/f026b-cym-plot-clean.png` | the F-127 figure, clean; compare `assets/f026-cym-plot-crowded.png` from round 1 |
| `assets/f026b-comparison-groups.png` | the new `type = "groups"` figure — **and the evidence for §7.1** |

---

## 10. Open questions

1. **§7.1 is the one thing I would fix before merge, and it is in `R/`.**
   `"posterior mean network"` at `R/plot_bgms.R:974` prints into the comparison
   vignette's figure. Every vignette source is clean; this string is not, and it
   is the maintainer's own rejected vocabulary appearing in the render he will
   judge. One word.
2. **How far do the rulings reach into package messages?** §7.1's table lists
   five more user-visible "network" strings in `R/` — three error messages, one
   advisory, one `verdicts()` note — plus `plot.bgms(type = "network")`, which
   the brief exempts. None is in a vignette. Whether they are in scope for CRAN
   is a maintainer call I did not make.
3. **The citations in the approved opening** (§7.2) — inline or in the following
   paragraph. I chose not to touch approved text; two brackets either way.
4. **"Diagnostics and Spike-and-Slab Summaries"** (§7.3) is the only remaining
   spike-and-slab mention on the front door and it is a title. Renaming the
   vignette is the only fix and it is outside this brief.
5. **Nine items or five**, carried forward from round 1 unchanged. The worked
   fit is 15 s of the 16.5 s knit and exists so the three-panel figure has
   content in all three panels (11 / 10 / 15 at nine items). Five items give a
   thin figure and halve the knit.

# Brief 02 — Hands-on evaluation of the new user-facing checks (MM)

Time box: ~2–3 hours at the console. You are evaluating as a *user* — a
methodologically sophisticated one — not reading code. Record findings as you
go, each with a severity (blocker / major / minor / note).

## Orientation — where this layer sits

0.2.0.0 adds a user-facing "checking and reporting" layer that post-processes
fit objects; none of it runs during sampling. The components and their call
paths:

- **`verdicts(fit)`** (`R/verdicts.R`) — three-category edge classification
  (presence / absence / undecided) from inclusion BFs; consumes the extractors
  (`extract_inclusion_bf()`, RB-based posterior inclusion probabilities from
  PR #182). Print method `print.bgms_verdicts`.
- **`extract_centrality(fit)`** (`R/centrality.R`) — posterior centrality
  measures with `summary()` and `plot()` methods; consumes posterior draws of
  the pairwise parameters.
- **`calibration_check(fit)`** (`R/calibration_check.R`) — model-calibration
  check with print/plot methods; consumes the simulate/predict machinery. Its
  scope was deliberately narrowed for continuous data — decision record in
  `dev/audit/2026-07-30-ppc-continuous-scope-decision.md`.
- **Network plots** (`R/plot_bgms.R`) — `plot.bgms`, `plot.bgmCompare`,
  `plot_edge_posterior`; the evidence-split display convention comes from the
  docs-site A6 work.
- **`prior_sensitivity_check(fit)`** (`R/prior_sensitivity.R`, PR #183) — the
  anchored sensitivity curve: `refit_engine.R` refits the model at scaled slab
  values (the 1× anchor is your own fit), `anchor_curve.R` assembles the
  continuous inclusion-BF curve by importance reweighting between anchors and
  inverse-variance pooling of anchors above an ESS floor. Print/plot methods.
  Design record: `dev/audit/2026-07-28-anchored-curve-spec.md`.

API design rationale for the whole layer:
`dev/audit/2026-07-30-user-facing-checks-api-proposal.md`.

## Setup

Use a **clean install of the frozen target** — not the Dropbox working tree
(stale `.o` files make it produce unloadable builds). Easiest: after the Opus
agent finishes brief 01, `.libPaths(c("~/bgms-review/lib", .libPaths()))` gives
you its rc1 install. Otherwise replicate its recipe (git archive
`v0.2.0.0-rc1` to a clean dir, `R CMD build` + `R CMD INSTALL`).

Data: `Wenchuan` (shipped with the package), your usual fit settings. Fit once
at defaults, keep the fit object for all components. Note fit runtime.

## Questions to answer (per component, in order)

Below I show the code that I ran, and I describte the issues that I found. 
Some issues were the same across fits, ill mark that specifically. 
So that I am not repeating the same issue/bugs for the different fits. 
I used R CMD BUILD on main.

---
```r
fit = bgm(Wenchuan)
verdicts(fit)
extract_centrality(fit)
```
The last function, to me, is not for users but for other packages to use. 
But that is okay. I did not check the numbers :-)

```r
calibration_check(fit)
```
Unclear what share_outside_band is; perhaps indicate it as a percentage? Or perhaps it already is, but is 0.3 30% or 0.3%?
Ideally we would like to plot this, right? 
Ah, plot(calibration_check(fit)) works! 
Perhaps there is a way to zoom in, as there are 17 plots here!

```r
plot(fit, type = "network", evidence_threshold = 5)
```
This is a qgraph object. 
Then qgraph should not be suggested I guess. 
Problem is that qgraph inherits a lot of additional dependencies that we dont want in our package.
This is a real challenge!

```r
plot(fit, type = "centrality")
plot_edge_posterior(fit, 1, 2, 3)
```
BUG: The title for this function has a massive string of numbers. 
FEATURE: For the density plot we might want to explore the rwmde method used in Bartos's sensitivity paper.

---
```r
fit = bgm(Wenchuan, variable_type = "blume-capel", baseline_category = 1)
verdicts(fit)

Edge verdicts at an inclusion Bayes factor of 10 (and 0.1 for absence):

  presence 37 | undecided 38 | absence 61   (136 indicators)

          parameter   pip log10_bf   verdict fragile
   intrusion-dreams 1.000  130.900  presence   FALSE
    intrusion-flash 1.000   10.784  presence   FALSE
    intrusion-upset 0.808    0.625 undecided   FALSE
  intrusion-physior 0.264   -0.445 undecided   FALSE
  intrusion-avoidth 0.032   -1.474   absence   FALSE
 intrusion-avoidact 0.034   -1.454   absence   FALSE
  intrusion-amnesia 0.038   -1.401   absence   FALSE
  intrusion-lossint 0.052   -1.261   absence   FALSE
  intrusion-distant 0.030   -1.513   absence   FALSE
     intrusion-numb 0.053   -1.249   absence   FALSE
... (126 more rows)

10 verdicts are Monte-Carlo fragile: a verdict boundary lies within two standard
errors of the evidence, so the verdict could change on a rerun. Run longer.
```
BUG: log_bf = -1.474 is called absent, but log(.1) = -2.3, so it should beundecided.
The same issue reoccurs at other numbers (and analyses/fits).
We also found the ame issue for presence classification.
```r
plot(prior_sensitivity_check(fit))
```
This shows one huge drifter intrusion-anger. See pdf in /assets

```r
extract_centrality(fit)
calibration_check(fit)
Error in x[, j] <- levels_list[[j]][x[, j] + 1L] : 
  number of items to replace is not a multiple of replacement length
```
BUG

```
plot(fit, type = "network", evidence_threshold = 5)
plot(fit, type = "centrality")
plot_edge_posterior(fit, 1, 2, 3)
```

---
```r
fit = bgm(Wenchuan, variable_type = "continuous", precision_graph_prior = "hierarchical")
verdicts(fit)
```
BUG: Unexpectedly, this took a very long time. 
After two minutes I killed the function. 
Then it does show results!

```r
extract_centrality(fit)
calibration_check(fit)
plot(fit, type = "network", evidence_threshold = 5)
plot(fit, type = "centrality")
plot_edge_posterior(fit, 1, 2, 3)
```
BUG explained: presence, BF = a thirty number string....

```r
plot_edge_posterior(fit, 1, 1, 3)
Error in plot_edge_posterior(fit, 1, 1, 3) : 
  Arguments 'variable1' and 'variable2' must name two different variables.
```
This could be more graceful.

---
```r
fit = bgm(Wenchuan, variable_type = "continuous", precision_graph_prior = "joint")
verdicts(fit)
extract_centrality(fit)
calibration_check(fit)
plot(fit, type = "network", evidence_threshold = 5)
plot(fit, type = "centrality")
plot_edge_posterior(fit, 1, 2, 3)
```
Nothing to report what is not said before. 

---
```r
fit = bgm(Wenchuan, variable_type = c(rep("continuous", 6), rep("ordinal", 6),
rep("blume-capel", 5)), precision_graph_prior = "joint", baseline_category = 1)
verdicts(fit)
extract_centrality(fit)
calibration_check(fit)
Error in x[, j] <- levels_list[[j]][x[, j] + 1L] : 
  number of items to replace is not a multiple of replacement length
```
BUG

```r
plot(fit, type = "network", evidence_threshold = 5)
plot(fit, type = "centrality")
plot_edge_posterior(fit, 1, 2, 3)
```
This is weird. 
The BF here is very small, while it was incredibly huge everywhere else. 
Must be a bug!

---
```r
fit = bgm(Wenchuan, variable_type = "continuous", update_method = "gibbs")
verdicts(fit)
extract_centrality(fit)
calibration_check(fit)
plot(fit, type = "network", evidence_threshold = 5)
plot(fit, type = "centrality")
plot_edge_posterior(fit, 1, 2, 3)
```


1. **Defaults**: run each function with no arguments beyond the fit. Is the
   default output something you would put in a paper or show a PhD student
   without caveats? If not, what exactly is off?
2. **Wording vs computation**: does every printed sentence match what is
   actually computed? Check each print method against its man page
   (`?verdicts`, `?calibration_check`, `?prior_sensitivity_check`, ...). We
   already know NEWS.md misdescribes the sensitivity pooling (finding F-001) —
   look for the same class of drift in the print methods themselves.
3. **Terminology**: do category names, evidence thresholds, and prior names
   match the tutorial/guidelines vocabulary (presence / absence / undecided;
   the 3/10/30 discussion)? List every mismatch — this is the
   decide-once-with-NS terminology set.
4. **Sensitivity check specifically**: at defaults on the Wenchuan fit —
   runtime; does the stability headline match what the curve shows; does the
   1×-anchor verdict set equal your reported analysis exactly (the design
   property that makes the check defensible); is the plot the one you would
   ship in the tutorial?
5. **Centrality**: which measures, on what scale, are the uncertainty
   intervals interpretable? Anything statistically misleading?
6. **Calibration**: what does it actually check, is that what the name
   promises, and does the output tell a user what to *do*?
7. **Plots**: publication-ready at defaults (labels, legends, sizing)? Is the
   evidence-split network consistent with the tutorial's three-panel figures?
8. **Failure behavior**: run `verdicts()` and `prior_sensitivity_check()` on a
   deliberately bad fit (e.g. 100 iterations, no warmup) — are the
   errors/warnings helpful or misleading?
   
```r
> fit = bgm(Wenchuan, warmup = 10, iter = 10)
18 rows with missing values excluded (n = 344 remaining).
To impute missing values instead, use na_action = "impute".
Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━⦘ 20/20 (100.0%)             
Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━⦘ 20/20 (100.0%)             
Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━⦘ 20/20 (100.0%)             
Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━⦘ 20/20 (100.0%)             
Total   (Sampling): ⦗━━━━━━━━━━━━━━⦘ 80/80 (100.0%)             
Elapsed: 1s | ETA: 0s                                           
NUTS issues:
  - Tree depth: 24 hits (60.0%) - consider max_depth > 10 
  - E-BFMI: 0.150 in chains 2, 3 
See vignette('diagnostics') for guidance.
Warning message:
In validate_sampler(update_method = update_method, target_accept = target_accept,  :
  warmup = 10 is very short for edge selection. Consider >= 300.

> verdicts(fit)
Edge verdicts at an inclusion Bayes factor of 10 (and 0.1 for absence):

  presence 33 | undecided 92 | absence 11   (136 indicators)

          parameter   pip log10_bf   verdict fragile
   intrusion-dreams 1.000   87.644  presence   FALSE
    intrusion-flash 1.000   10.106  presence   FALSE
    intrusion-upset 0.967    1.468  presence    TRUE
  intrusion-physior 0.338   -0.292 undecided    TRUE
  intrusion-avoidth 0.441   -0.104 undecided    TRUE
 intrusion-avoidact 0.201   -0.599 undecided    TRUE
  intrusion-amnesia 0.175   -0.674 undecided    TRUE
  intrusion-lossint 0.164   -0.708 undecided    TRUE
  intrusion-distant 0.474   -0.046 undecided   FALSE
     intrusion-numb 0.281   -0.407 undecided   FALSE
... (126 more rows)

102 verdicts are Monte-Carlo fragile: a verdict boundary lies within two standard
errors of the evidence, so the verdict could change on a rerun. Run longer.
```
Not super bad, though it still has 34 decided verdicts.. 

```r
> prior_sensitivity_check(fit)
Refitting at 4 scales (0.4x, 0.63x, 1.6x, 2.5x), plus a 1.6x repeat for the noise band.
Each refit runs 4 chains; warmup and sampling are shown per chain.

[1/5] 0.4x scale
Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━⦘ 1500/1500 (100.0%)         
Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━⦘ 1500/1500 (100.0%)         
Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━⦘ 1500/1500 (100.0%)         
Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━⦘ 1500/1500 (100.0%)         
Total   (Sampling): ⦗━━━━━━━━━━━━━━⦘ 6000/6000 (100.0%)         
Elapsed: 3m 38s | ETA: 0s                                       

[2/5] 0.63x scale
Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━⦘ 1500/1500 (100.0%)         
Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━⦘ 1500/1500 (100.0%)         
Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━⦘ 1500/1500 (100.0%)         
Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━⦘ 1500/1500 (100.0%)         
Total   (Sampling): ⦗━━━━━━━━━━━━━━⦘ 6000/6000 (100.0%)         
Elapsed: 3m 41s | ETA: 0s                                       

[3/5] 1.6x scale
Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━⦘ 1500/1500 (100.0%)         
Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━⦘ 1500/1500 (100.0%)         
Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━⦘ 1500/1500 (100.0%)         
Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━⦘ 1500/1500 (100.0%)         
Total   (Sampling): ⦗━━━━━━━━━━━━━━⦘ 6000/6000 (100.0%)         
Elapsed: 3m 40s | ETA: 0s                                       

[4/5] 2.5x scale
Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━⦘ 700/1500 (46.7%)           
Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━⦘ 708/1500 (47.2%)           
Chain 3 (Sampling): ⦗━━━━━━╺━━━━━━━⦘ 693/1500 (46.2%)           
Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━⦘ 698/1500 (46.5%)           
Total   (Sampling): ⦗━━━━━━━━━━━━━━⦘ 2799/6000 (46.7%)          
Elapsed: 1m 41s | ETA: 1m 55s 

....

Prior sensitivity check: are the edge verdicts robust to the slab scale?
Bayes-factor curve from 0.4x to 2.5x the chosen scale (anchors at 0.4x, 0.63x,
1x, 1.6x, 2.5x; the 1x anchor is the original fit); 136 edges.

89 of 136 verdicts hold across the whole 0.4x-2.5x range; the exceptions are
named below.

  robust (same verdict at every scale)     35
  changed, within run-to-run noise          3
  changed, beyond run-to-run noise          3
  not certifiable (too noisy to assess)    95

3 edges' verdicts genuinely depend on the scale:
  edge              0.4x       0.63x      1x         1.6x     2.5x   
  dreams-numb       undecided  undecided  undecided  absence  absence
  upset-concen      undecided  undecided  undecided  absence  absence
  physior-avoidact  undecided  absence    undecided  absence  absence

95 edges are too noisy to assess: intrusion-upset, intrusion-physior,
intrusion-avoidact, intrusion-amnesia, intrusion-lossint, intrusion-numb,
intrusion-sleep, intrusion-anger, intrusion-hyper, intrusion-startle, and 85
more (see $edges).
Their Bayes factor sits within Monte Carlo error of an evidence threshold, or
their chains disagree on the verdict, at the chosen scale itself; a rerun with
a fresh seed could flip them without any prior change. Run more iterations to
settle these verdicts before reading their sensitivity.

Verdict counts by scale (at the anchors):
            0.4x 0.63x 1x 1.6x 2.5x
  presence    37    36 33   35   34
  undecided   64    49 92   28   25
  absence     35    51 11   73   77
More absence at wider scales is expected: a wider slab strengthens evidence
against borderline edges.

Note: the chosen scale (1) is much wider than the estimated interactions (about
0.121 [0.103, 0.143]); absence verdicts in particular depend on this choice.

Method:  43-point curve from 5 anchor fits (0.4x to 2.5x the chosen scale),
         joined by importance reweighting; the 1x anchor is the original fit.
         Points with reweighting effective sample size below 400 are not shown.
Refits:  5 nuts refits, warm-started from the original fit, 1094 s total.
Noise:   two identical refits at 1.6x differed by up to 0.41 log10 BF across
         threshold-relevant edges; verdict moves smaller than that are reported
         as run-to-run noise, not prior sensitivity.
         2 edges saturated in one of the two and were left out of that spread.
See ?prior_sensitivity_check for the full construction.
Warning messages:
1: In min(c(wc$ebfmi_first_half, wc$ebfmi_second_half), na.rm = TRUE) :
  no non-missing arguments to min; returning Inf
2: In max(wc$var_ratio, na.rm = TRUE) :
  no non-missing arguments to max; returning -Inf

```
This one is really slow due to the bad initial fit. Which is expected.
Results are difficult to say!
   
9. **The embarrassment test**: anything here that a sharp reviewer of the
   tutorial paper could use against the package?

## Deliverable

Write your report to `dev/review-2026-08/reports/02-user-facing-checks.md`:

1. **What was done** — fit settings, versions, which functions exercised.
2. **Findings** — numbered, severity-tagged (blocker / major / minor / note).
3. **Evidence** — console excerpts / saved plots (drop images in
   `dev/review-2026-08/reports/assets/` if useful).
4. **Open questions** — anything needing a code-level answer; these become
   your next code-read briefs.

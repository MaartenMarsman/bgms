# Brief 07 — bgmCompare hands-on pass (MM)

Time box: ~1 hour. Same register as brief 02: evaluate as a user, record
findings with severities as you go. Runnable now that the checking-layer fixes
are merged (develop `adf87013`) — install from a clean export of develop, not
rc1 and not the Dropbox tree:

```sh
mkdir -p ~/bgms-review/dev-install && cd ~/bgms-review/dev-install
git -C ~/Library/CloudStorage/Dropbox/Projecten/R/bgms archive develop --prefix=bgms-dev/ | tar -x
R CMD build bgms-dev && R CMD INSTALL --library=~/bgms-review/lib-dev bgms_0.2.0.0.tar.gz
```

## Orientation

The bgmCompare surface you have not yet touched: `verdicts()` on difference
indicators, `plot()` group panels + the verdict-encoded difference network,
per-group and difference `extract_centrality()`, and the difference-scale
`prior_sensitivity_check()` trace. Its sampler was rewritten this release
(+912 changed lines), it carries the breaking association-scale change, and it
has no C++ test interface — the Opus cross-path validation (brief 06) covers
the numbers; this pass covers what only you can judge.

Data: your usual two-group example (the comparison vignette's setup is fine).
Fit once at defaults; note runtime.

## Questions

1. `verdicts(fit)`: are difference verdicts labeled unambiguously (which
   group pair, which scale)? 
> Differences are modeled per factor (the vector of group differences) not per
group pair.    
   
   Does the print carry the scale-contingency
   caveat the vignette promises (difference verdicts near a threshold are
   scale-contingent; calibration under study)?
> I have no idea, see the reports below.
   
2. Fragility: the flag was validated for ordinal single-group indicators only
   (UFC-2) — does the print say so for difference indicators, or does it
   overclaim?
> Yes. We should analyze this and validate also for bgmCompare, so that this flag can be removed.

3. Units: verdicts now print natural-log "log BF"; the sensitivity trace still
   speaks log10 (pending the terminology sync). How jarring is the mix in one
   bgmCompare session? This feeds your unit decision.
> Jarring. Lets make sure we have only natural log everywhere.

4. Plots at defaults: group panels + difference network — publication-ready?
   Consistent with the single-group conventions you approved?
> See below

5. `prior_sensitivity_check` over `difference_scale`: runtime, does the
   headline match the curve, do the 1×-anchor verdicts equal your reported
   analysis?
> See below

6. Known, expected: `calibration_check()` on a bgmCompare fit errors with "no
   applicable method" — the prepared patch was deliberately not applied
   (F-021; lead recommends deferring to 0.2.1). Confirm you agree, or say
   apply.
> Correct. I would suggest adding something for 0.2.0
   
7. The embarrassment test, group-comparison edition: anything the companion
   paper's reviewers could use against the package?

## Deliverable

`dev/review-2026-08/reports/07-bgmcompare-user-pass.md` — What was done /
Findings (severity-tagged) / Evidence (console excerpts, plots to
`reports/assets/`) / Open questions.

## Analyses

```r
fit = bgmCompare(Boredom[,-1], group_indicator = as.integer(factor(Boredom[,1])), seed = 1)
```
This construction should be handled internally: group_indicator = as.integer(factor(Boredom[,1]))
Before the pb starts there is a 5-10s silence.

```r
> verdicts(fit)
Edge verdicts at an inclusion Bayes factor of 10 (and 0.1 for absence):
presence: log BF > 2.30; absence: log BF < -2.30

  presence 2 | undecided 9 | absence 17   (36 indicators)
  8 indicator(s) were never updated and carry no verdict.

                            parameter   pip log_bf   verdict fragile
                    loose_ends (main)   NaN     NA      <NA>   FALSE
      loose_ends-entertain (pairwise) 0.164 -1.626 undecided   FALSE
     loose_ends-repetitive (pairwise) 0.010 -4.639   absence   FALSE
    loose_ends-stimulation (pairwise) 0.529  0.117 undecided   FALSE
      loose_ends-motivated (pairwise) 0.055 -2.837   absence   FALSE
  loose_ends-keep_interest (pairwise) 0.018 -3.988   absence   FALSE
     loose_ends-sit_around (pairwise) 0.142 -1.796 undecided   FALSE
 loose_ends-half_dead_dull (pairwise) 0.371 -0.528 undecided   FALSE
                     entertain (main)   NaN     NA      <NA>   FALSE
      entertain-repetitive (pairwise) 0.014 -4.232   absence   FALSE
... (26 more rows)

3 verdicts are Monte-Carlo fragile: a verdict boundary lies within two standard
errors of the evidence, so the verdict could change on a rerun. Run longer.

The fragility flag is not validated for difference indicators: its
operating point was established on single-network edge indicators only.
Read it as an indication that a verdict sits near a boundary, not as a
calibrated error rate.
```
Observe that main effects leaked into the verdicts. 
I believe this could work once "main_difference_selection = TRUE", but it is
off by default. 
Lets try:
```r
fit2 = bgmCompare(Boredom[,-1], group_indicator = as.integer(factor(Boredom[,1])), 
seed = 1, main_difference_selection = TRUE)
> verdicts(fit2)
Edge verdicts at an inclusion Bayes factor of 10 (and 0.1 for absence):
presence: log BF > 2.30; absence: log BF < -2.30

  presence 2 | undecided 8 | absence 26   (36 indicators)

                            parameter   pip log_bf   verdict fragile
                    loose_ends (main) 0.005 -5.199   absence   FALSE
      loose_ends-entertain (pairwise) 0.087 -2.347   absence    TRUE
     loose_ends-repetitive (pairwise) 0.011 -4.457   absence   FALSE
    loose_ends-stimulation (pairwise) 0.707  0.881 undecided   FALSE
      loose_ends-motivated (pairwise) 0.015 -4.199   absence   FALSE
  loose_ends-keep_interest (pairwise) 0.017 -4.054   absence   FALSE
     loose_ends-sit_around (pairwise) 0.754  1.118 undecided   FALSE
 loose_ends-half_dead_dull (pairwise) 0.246 -1.118 undecided   FALSE
                     entertain (main) 0.001 -6.781   absence   FALSE
      entertain-repetitive (pairwise) 0.009 -4.751   absence   FALSE
... (26 more rows)

3 verdicts are Monte-Carlo fragile: a verdict boundary lies within two standard
errors of the evidence, so the verdict could change on a rerun. Run longer.

The fragility flag is not validated for difference indicators: its
operating point was established on single-network edge indicators only.
Read it as an indication that a verdict sits near a boundary, not as a
calibrated error rate.
```
Yay!

The "Run longer." is quite direct.

```r
plot(fit)
plot(fit2)
```
The "main-effect differences not selected" is a bit weird in the plot.

The "node: main-effect difference" is weird.

Legend overlaps the plot.

For "plot(fit2)", could we plot the main effect tests similar to how node 
predictions are shown in qgraph (with the bands around the nodes), perhaps based 
on pip, with pip = 1 showing the filled full circle, and pip = .5 showing half.

Do these plots show differences or bayes factors. 
For bgm it is BF, we should do that here too.

```r
plot(fit, type = "groups")
plot(fit2, type = "groups")
```
Similar issues as above.

```r
plot(fit, type = "centrality")
plot(fit2, type = "centrality")
```
No idea what these mean. Lets not do centrality for differences. 

```r
ps = prior_sensitivity_check(fit)

> Warning messages:
1: In (p/(1 - p))/prior_odds :
  longer object length is not a multiple of shorter object length
2: In (p/(1 - p))/prior_odds :
  longer object length is not a multiple of shorter object length
```

These warning messages are not good!

```r
plot(ps)
```
Weird result, and effect that has increasing evidence for inclusion with scale 
(see pdf in /assets): sit_around-half_dead_dull


```r
ps2 = prior_sensitivity_check(fit2)

> Warning messages:
1: In (p/(1 - p))/prior_odds :
  longer object length is not a multiple of shorter object length
2: In (p/(1 - p))/prior_odds :
  longer object length is not a multiple of shorter object length
```

```r
> calibration_check(fit)
Error in UseMethod("calibration_check") : 
  no applicable method for 'calibration_check' applied to an object of class "c('bgmCompare', 'S7_object')"
```
I do think we need to do something here: Research item.







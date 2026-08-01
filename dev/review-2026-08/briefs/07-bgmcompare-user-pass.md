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
   group pair, which scale)? Does the print carry the scale-contingency
   caveat the vignette promises (difference verdicts near a threshold are
   scale-contingent; calibration under study)?
2. Fragility: the flag was validated for ordinal single-group indicators only
   (UFC-2) — does the print say so for difference indicators, or does it
   overclaim?
3. Units: verdicts now print natural-log "log BF"; the sensitivity trace still
   speaks log10 (pending the terminology sync). How jarring is the mix in one
   bgmCompare session? This feeds your unit decision.
4. Plots at defaults: group panels + difference network — publication-ready?
   Consistent with the single-group conventions you approved?
5. `prior_sensitivity_check` over `difference_scale`: runtime, does the
   headline match the curve, do the 1×-anchor verdicts equal your reported
   analysis?
6. Known, expected: `calibration_check()` on a bgmCompare fit errors with "no
   applicable method" — the prepared patch was deliberately not applied
   (F-021; lead recommends deferring to 0.2.1). Confirm you agree, or say
   apply.
7. The embarrassment test, group-comparison edition: anything the companion
   paper's reviewers could use against the package?

## Deliverable

`dev/review-2026-08/reports/07-bgmcompare-user-pass.md` — What was done /
Findings (severity-tagged) / Evidence (console excerpts, plots to
`reports/assets/`) / Open questions.

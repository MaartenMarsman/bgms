# Brief 16 — F-075: bgmCompare magnitude inflation — replicate, localize, propose the fix (Opus agent; RELEASE GATE)

You are the execution agent in the pre-CRAN review of the R package **bgms**.
You share no context with the review lead; everything you need is here.
**ANALYSIS ONLY — no package code changes.** Your deliverable is the
mechanism and a concrete fix proposal; the fix itself lands in a separately
authorized batch once the maintainer and lead review your evidence. The
maintainer has ruled (2026-08-01): this is cleared up and fixed FOR the
0.2.0.0 release. It gates submission.

## Machine budget

This machine has 15 cores and is SHARED. This brief is the one authorized
HEAVY lane: run fits **strictly sequenced** (one at a time, `cores = 4`,
`chains = 4`), log per-fit wall times, and prefer overnight-style batches.
Phase 0 costs nothing and comes first. Do not launch anything beyond the
stated grid without flagging it in the report instead.

## The finding you are resolving (F-075, blocker)

On simulated two-group data with planted differences, the compare path's
per-group pairwise estimates are badly inflated relative to what two
separate `bgm()` fits recover from the **byte-identical data** at identical
n (2000/group, one data seed):

| n = 2000/group | 2 × `bgm()` | one `bgmCompare()` |
|---|---|---|
| group 1 slope vs truth | 0.963 [0.921, 1.005] | 1.270 [1.159, 1.382] |
| group 2 slope vs truth | 1.034 [0.931, 1.136] | 1.926 [1.599, 2.253] |
| group 2 rmse | 0.051 | 0.214 |
| noise on true-zero differences | sd 0.041, max 0.110 | sd 0.156, max 0.830 |

Already REFUTED (do not re-litigate; report 09 has the evidence): units/×2
convention (four independent confirmations of health), selection artefact
(selection OFF is worse), `simulate_mrf()` mixing (converged by iter 250),
truth recoverability (`bgm()` recovers it), regime conditioning. Also known:
part of the CONTRAST inflation is shared small-sample behaviour (`bgm()`'s
own contrast slope at this n is 1.500 [1.314, 1.687]) — only the group-level
slopes and the noise inflation are compare-specific. Which group inflates
FLIPPED between the n = 400 and n = 2000 cells in report 09, so run-to-run
variability is large: **nothing is scoped off a single seed** — that is why
Phase A exists. On real Wenchuan splits at n = 172/group the effect did NOT
appear (difference means sd 0.046–0.058, near-nominal intervals); the
related-but-distinct F-074 records a small (6–9%) level offset on real data.

## Setup

- Repo (Dropbox; do NOT switch its checked-out branch or build in it):
  `/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms`
- Clean export of `develop` AT OR AFTER `5b850410`, built into a private
  library (as report 09 did; the compare sampler is unchanged since 09's
  export, so results are comparable).
- REQUIRED READING: `dev/review-2026-08/reports/09-crosspath-validation.md`
  §3.2 and §09-3 (the construction you are replicating), and FINDINGS rows
  F-074/F-075/F-076.
- Prior scripts and fit objects live at `~/bgms-review/val09/out/`
  (`item2.rds`, `control.rds`, `16_apples.R` …). REUSE what is sound; the
  planted-truth construction, restated so this brief stands alone:
  ```r
  # truth: fit real Wenchuan, then plant differences on 4 pairwise entries
  base = bgm(Wenchuan[, 1:10], edge_selection = FALSE,
             interaction_prior = cauchy_prior(1), seed = 31)   # -> OM, main
  # pairs 3, 14, 27, 41 (upper-tri order); delta = (0.05, 0.10, 0.20, 0.40)
  # OM1 = OM - delta/2 on those pairs; OM2 = OM + delta/2
  # data: simulate_mrf(n, 10, num_categories = 4, pairwise = OMg,
  #                    main = main, iter = 1000, seed = <data seed>)
  ```

## Phases

### Phase 0 — mine what already exists (no compute)

Load 09's saved fits (`item2.rds`, apples controls). For the bad cells,
report the convergence picture on the DIFFERENCE and GROUP parameters:
split-Rhat, n_eff, indicator transition counts, and (if stored) per-chain
means — is the inflation carried by all chains or by excursions in some?
If the compare fits show poor mixing where the `bgm()` fits are healthy,
that reframes everything downstream — check this FIRST and say so plainly.

### Phase A — the replication gate (the maintainer's precondition)

Ten data seeds at n = 2000/group, same truth, same planted δ. Per seed:
one `bgmCompare()` at defaults (`difference_selection = TRUE`) + two
matched `bgm()` fits on the same halves (priors matched as report 09 §3.3
did). ~30 sequenced fits; log runtimes. Deliver, across seeds:

- the distribution of compare-vs-truth group slopes against bgm-vs-truth
  (mean, sd, and per-seed pairs);
- the noise ratio on true-zero differences (compare/bgm), per seed;
- which group inflates, per seed — RANDOM direction across seeds says
  instability/variance; a SYSTEMATIC direction says structure;
- planted-edge posterior means vs δ (does the δ = 0.40 edge over-estimate
  in most seeds, or was 09's 1.21 a tail draw?);
- detection: PIPs on planted vs unplanted edges per seed (does the clean
  detection story survive replication?).

### Phase B — localize the mechanism (cheapest discriminator first)

On 2–3 seeds chosen from Phase A (one typical, one worst):

1. **`difference_scale` sweep**: 1 (default) → 0.5 → 0.25 → 0.1, defaults
   otherwise. If inflation shrinks toward the realistic ~0.1 scale, the
   PRIOR DEFAULT is the lever; if it persists, the sampler is.
2. **Selection interplay** at the winning scale: `difference_selection`
   TRUE vs FALSE (09 saw FALSE worse at scale 1 — does that invert?).
3. **Adaptation probe** on one bad cell: double `iter`/`warmup`
   (4000/4000). If the inflation melts with longer adaptation, this is a
   warmup/step-size problem on the difference posterior, not a model
   property — a different fix entirely.
4. Only if 1–3 all fail to move it: one cell with main differences held
   out (if the API allows pinning `main_difference_*`), to test whether
   main-difference freedom leaks into the pairwise block. Flag before
   running anything beyond this.

### Phase C — the fix proposal (deliverable, not implementation)

State, with Phase A/B evidence attached: what the mechanism is, what the
release fix is, and what it costs. The admissible shapes, in the lead's
prior order (overturn with evidence freely):

- a `difference_scale` default change (+ refit guidance in NEWS and the
  comparison vignette, + a `verdicts()` caveat update);
- a compare-path adaptation/default change (`iter`/`warmup`/`target_accept`
  for compare fits) with the measured before/after;
- a sampler-side defect you can point at in `src/models/bgmCompare/` —
  name the site and the failure mode; do NOT patch it;
- and in every case: the operating-characteristics table that goes in the
  documentation, because users deserve the number whatever the fix is.

If Phase A does NOT replicate the inflation (09's seed was a tail event),
say so with the distribution as evidence — that outcome changes the release
answer and is a perfectly good result.

## Verification gate

Seeds, runtimes, and derived tolerances stated for every claim (report 09's
standards). Phase A's distribution reported for ALL ten seeds, including
any that look fine — selective reporting is a FAIL. Every Phase B cell
tabulated even when null.

## Deliverable

Write
`/Users/maartenmarsman/Library/CloudStorage/Dropbox/Projecten/R/bgms/dev/review-2026-08/reports/16-f075-resolution.md`:
What was done / Findings (severity-tagged; the mechanism verdict FIRST) /
Evidence (per-phase tables, figures to `reports/assets/`) / The fix
proposal / Open questions.

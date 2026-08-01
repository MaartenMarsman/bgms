# Brief 08 — MM curriculum step 2: the zratio law and the trust gauge (math read)

Reader: the maintainer. Goal: by the end of this read you own the
hierarchical precision-graph prior's surface approximation — the newest,
most intricate machinery in 0.2.0.0 and the subsystem where late defects
were found (PRs #172/#193/#194). Nothing here needs a computer beyond
opening files; the numerical certification already exists (report 04: all
five certificates PASS against the gold bank in `dev/validation/`). Your
job is the part no agent can do: judge whether the MATH is right and
whether its guard rails guard the right thing. Timebox ~2 hours; depth
over completeness.

## Reading order

1. **`src/models/ggm/zratio_law.h` — the law itself.** What is being
   approximated: the normalizing-constant RATIO for an edge move under the
   hierarchical prior. Read for three things: (a) the central identity —
   is the quantity the surface interpolates actually the quantity the MH
   ratio needs; (b) the parameterization of the surface (what are its
   axes, and are they the sufficient statistics of the ratio or a
   convenient projection); (c) the validity envelope the file CLAIMS —
   write it down, you will compare it with what the gauge PROTECTS.
2. **`R/zratio_surfaces.R` + `R/zratio_tables.R` — where surfaces come
   from.** Anchor hulls, the certified-constants band (shape ∈ [0.5, 20]),
   caching. Read for: what happens at and beyond the hull boundary —
   extrapolation, clamping, or refusal — and whether the certified band is
   enforced or advisory.
3. **`src/models/ggm/zratio_engine.{h,cpp}` — routing.** Which requests go
   to the certified surface, which to the additive kernel, which to exact
   computation (isolated edges). Read for: is the routing predicate
   exactly the certification predicate, or can a request be served by a
   path whose guarantees don't cover it.
4. **`R/zratio_gauge.R` + `src/models/ggm/zratio_gauge.h` — the trust
   gauge.** Two assessment sweeps auditing a capped number of edge moves;
   `harm_pred` projecting coherent error; the advice ladder (raise
   `bgms.zratio_gauge_sweeps` → refit → only then mention
   `precision_graph_prior = "joint"`). Read for: (a) does the audited
   quantity bound the error that actually moves a recovered network;
   (b) the documented blind spot — rare edge-specific failures — is it
   stated honestly everywhere the gauge is described; (c) is the cap
   (fixed cost per chain) a sample the user can reason about.

## The decision this read carries — F-059 (hook or accept)

FINDINGS row F-059, verbatim substance: post-#193, the additive kernel's
only serving band (shape < 0.5) coincides EXACTLY with the
uncertified-constants band (certified range [0.5, 20]) — so no certified
cell can exercise the additive kernel end to end, and the gauge's firing
test necessarily runs on unscored constants. That is acceptable for a
detector test and documented in-test. The question: do you want a
deliberate TEST-ONLY surface-withhold hook so the additive path can be
exercised inside the certified envelope, or do you accept the coverage
asymmetry as a documented property?

Lead's framing, for and against: a hook buys an end-to-end certified
exercise of a path real fits do take, at the price of test-only machinery
inside production routing (a seam that must itself be guarded); accepting
costs nothing now and leaves the additive band's certification implicit in
the law's structure rather than tested. Either answer is defensible;
record WHY, not just which.

## Deliverable

`dev/review-2026-08/reports/08-zratio-math-read.md`, informal is fine:
per component, what convinced you and what did not (one paragraph each);
the F-059 decision with its why; any new findings (mismatches between the
law's claimed envelope, the routing predicate, and the gauge's protection
— that triple agreeing is the whole safety argument). Anything you flag
becomes a finding the lead triages.

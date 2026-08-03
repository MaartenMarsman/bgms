# Brief 08 v2 — MM curriculum step 2: the zratio machinery (math read) — REISSUED

Reader: the maintainer. REISSUED 2026-08-02 after v1 stalled — the stall was
diagnostic, and two of its causes were on the review, not on you:

1. **v1 sent you to the wrong file first.** `src/models/ggm/zratio_law.h` is
   DORMANT insurance code — its own header says it is not wired into the
   deployed surface build, is reachable only through a test-only export, and
   "is not used in the paper". It is also the single most mathematically
   dense file in the subsystem (spectral solves, damped Picard iteration).
   Starting there guaranteed the lost feeling. Lead error, recorded (F-098).
2. **The derivations you need were never assembled.** They live across
   `sv/ggm_paper` (the analytic-correction companion — the thing the code
   comments call "the companion"), `sv/Z`, the related SV projects, and the
   bgms dev notes. Brief 18 assembles them into ONE self-contained
   curriculum document, which the lead verifies before you read it.
3. The four code comments that deferred to "the companion" are F-097 and
   get rewritten self-contained (brief 17).

## What you can do TODAY (10 minutes, no prerequisites)

Read `src/models/ggm/zratio_engine.h` lines 111-131 — the engine's own
statement of the problem and its solution shape. It is genuinely good
documentation, and it is the paragraph everything else hangs off: every
between-graph move needs J = Z(Γ⁻)/Z(Γ⁺); J is evaluated from three integer
counts of the toggled edge's mediating neighbourhood through a two-moment
saddle over tabulated pair integrals; no sampling runs inside the acceptance
step; distinct count tuples are cached. If that paragraph reads clearly, the
rest of the read is that paragraph unfolded.

## The read itself — AFTER brief 18's curriculum lands (lead-verified)

Work through `dev/review-2026-08/orientation/zratio-derivation-curriculum.md`
in its order (model → locality → saddle → surface → routing → gauge). Your
job at each step is the one only you can do: **is the mathematical argument
right?** The curriculum's rules require every step to be derived, sourced,
or flagged — the flagged gaps are where your judgment matters most.

With the curriculum open, then read against the code, and judge:

1. **Routing vs certification** (`zratio_engine.{h,cpp}`,
   `R/zratio_tables.R`): is the predicate that routes a request to a path
   exactly the predicate under which that path is certified — or can a
   request be served by a path whose guarantees don't cover it?
2. **The surface's edges** (`R/zratio_surfaces.R`): hull clamps, the log-
   moment clamp, boundary-slope extension, additive fallback — at every
   boundary, what happens, and is it announced or silent?
3. **The gauge** (`R/zratio_gauge.R`, `zratio_gauge.h`): does the audited
   quantity bound the error that actually moves a recovered network
   (`harm_pred`'s projection); is the documented blind spot (rare
   edge-specific failures) stated everywhere the gauge is described; is the
   fixed per-chain cap a sample a user can reason about?
4. **The triple agreement** — the law's claimed envelope, the routing
   predicate, and the gauge's protection: these three agreeing is the whole
   safety argument. Any mismatch is a finding.

The DORMANT analytic law needs only its header note (you have read the
relevant part already) — no deep read for the release.

## The decision this read carries — F-059 (hook or accept)

Unchanged from v1. Post-#193, the additive kernel's only serving band
(shape < 0.5) coincides exactly with the uncertified-constants band
(certified range [0.5, 20]) — so no certified cell exercises the additive
kernel end to end, and the gauge's firing test necessarily runs on unscored
constants (documented in-test). Do you want a deliberate TEST-ONLY
surface-withhold hook so the additive path can be exercised inside the
certified envelope, or do you accept the coverage asymmetry as a documented
property? A hook buys an end-to-end certified exercise of a path real fits
take, at the price of test-only machinery inside production routing;
accepting costs nothing now and leaves the additive band's certification
implicit in the law's structure. Either is defensible; record WHY.

## How to read: interactively

Read with the lead's session open. Paste ANY line, step, or sentence that
loses you — from the curriculum or the code — and you get the orientation
answer on the spot. Your confusion is data: every place you get lost
despite the curriculum becomes a documentation finding with your name on
the discovery. That is not a fallback mode; it is the intended mode.

## Deliverable

`dev/review-2026-08/reports/08-zratio-math-read.md`, informal is fine: per
curriculum section, what convinced you and what did not; the F-059 decision
with its why; any new findings. The safety context, so the weight sits
right: the numerics are already independently certified (report 04: all
five certificates PASS against the gold bank) and the gauge audits every
production fit at runtime — your read is about OWNING the argument, not
about being its last line of defense.

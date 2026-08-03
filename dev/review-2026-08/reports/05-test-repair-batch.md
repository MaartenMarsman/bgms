# Report 05 — Test-repair and release-hygiene batch

Six items from briefs 03 and 04, landed as six commits on a fix branch.
**Code and tests changed; no sampler or engine behaviour changed except the
q ≤ 3 hierarchical crash guard.**

Agent: brief 05. Date: 2026-08-01.

---

## 0. Environment and branch

| | |
|---|---|
| Branch | `fix/test-repair-batch`, based on `develop` at `289d028d` (which contains `adf87013`) |
| Worktree | `~/bgms-review/wt-fix2/` (the Dropbox checkout was never built in and never had its branch switched) |
| Library | `~/bgms-review/lib-fix2` — clean `R CMD INSTALL --preclean` of the branch head |
| R | 4.6.0 (2026-04-24), aarch64-apple-darwin23 |
| Machine | Darwin 25.4.0, 15 cores, Apple clang 21.0.0 |
| Pre-commit | `styler::style_pkg(style = bgms_style)` run; `lintr::lint_package()` → **no lints found**; `roxygen2::roxygenise()` run |

Styler wanted to reformat six files this branch does not touch
(`R/zratio_gauge.R`, `R/zratio_tables.R`,
`tests/testthat/fixtures/make_zratio_law_reference.R`,
`tests/testthat/test-bgmCompare.R`,
`tests/testthat/test-zratio-isolated-edge-routing.R`,
`tests/testthat/test-zratio-surface.R`) and to reflow ~70 pre-existing lines of
`R/zratio_surfaces.R`. All of that was reverted, so the diff stays scoped —
same convention report 03 used.

### Commits

| | Commit | Item |
|---|---|---|
| 1 | `d29989cd` | `test:` retire the pre-#193 surface fence from the slow tier (F-048) |
| 2 | `afd68b56` | `test:` move the slow tier to the 0.2.0 prior arguments (F-050) |
| 3 | `7d24331e` | `fix(correction):` scope the correction-table cache key to the version (F-054) |
| 4 | `d6d804e2` | `fix(zratio):` serve hierarchical fits at 2 and 3 variables (F-055) |
| 5 | `50281f7a` | `build:` keep the compliance tier and the fixture generator out of the tarball (F-017) |
| 6 | `8e76b87e` | `docs(prior-sensitivity):` state what the plotted Bayes-factor cap means (F-051) |

```
 13 files changed, 330 insertions(+), 71 deletions(-)
```

Nothing was pushed.

---

## 1. What was done, per item

### F-048 — the three stale slow-tier assertions (`d29989cd`, test-only)

**`test-zratio-surface-build.R:147-168`.** The assertion expected
`zratio_build_surfaces()` to return NULL at shapes 2.5 and 5, both inside the
shipped `[0.5, 10]` fence. Rewritten to assert the shipped contract, and to read
the range off its own constants rather than restating literals:

```r
lo = bgms:::.zratio_surface_shape_lo
hi = bgms:::.zratio_surface_shape_hi
expect_equal(c(lo, hi), c(0.5, 10))
for(shape in c(lo, 2.5, 5, hi))  expect_false(is.null(build(shape)))   # inside
for(shape in c(0.25, 12))        expect_null(build(shape))             # outside
```

Both fence endpoints are now covered (0.5 and 10), plus the two interior shapes
the brief named and one cell below the fence. Cost is unchanged: 4 real builds
before, 4 after (the two out-of-range cells return before any anchor runs).
Measured on the branch install, `max_size = 8L, cores = 1L`:

```
shape 0.25  build 0.0s -> NULL       shape 2.5   build 0.4s -> surface
shape 0.5   build 0.7s -> surface    shape 5     build 0.4s -> surface
shape 2     build 0.4s -> surface    shape 10    build 0.4s -> surface
                                     shape 12    build 0.0s -> NULL
```

**`test-zratio-gauge.R:190` — repair (b), the preferred one, plus (a).** I
reproduced report 04's diagnosis on this build before touching the test, and
then looked for a cell where the additive kernel still serves. It exists: the
fence's lower edge. Below `.zratio_surface_shape_lo = 0.5` the engine keeps the
additive-counts saddle, and the same fixture is genuinely harmful there.
Measured (p = 16, 1200 samples, no data, `beta_bernoulli_prior(9, 1)`,
`seed = 7`):

| Gamma shape | kernel served | amplification | `harm_pred` | threshold | `harm_flag` | `flip_rate` |
|---:|---|---:|---:|---:|---|---:|
| 0.25 | additive | 12.391 | 0.05251 | 0.01 | **TRUE** | 0.00099 |
| **0.40** | **additive** | **12.144** | **0.04306** | 0.01 | **TRUE** | 0.00115 |
| 2 | surface | 17.194 | 0.00794 | 0.01 | FALSE | 0.00114 |

The shape-2 row reproduces report 04's shipped-fence measurement to the digit
(17.19 / 0.00794 / FALSE), which confirms the harm channel is intact and the
fixture was the stale part.

So the fixture is factored into `biased_evidence_free_fit(shape)` and the file
now carries both halves:

- **(b)** the firing test moves to shape 0.4 — the additive kernel, which is
  what this channel exists to police — and additionally asserts
  `harm_pred > harm_threshold` and that rung 1 stays quiet (`flip_rate < 0.01`),
  which is the property that makes the harm channel a separate channel;
- **(a)** a new test, *"the deployed surface holds that fixture under the harm
  threshold"*, keeps the shape-2 cell as the channel's negative control on a
  cell that is genuinely at risk: amplification still first-order (> 5),
  `harm_pred < harm_threshold`, flag FALSE.

Documented in the test comments, including that the Gamma-shape *constants* are
themselves unscored below 0.5 (the cell warns; the warning is suppressed and the
reason stated) — the fixture is chosen for the kernel it reaches, not as a
certified cell.

**`test-prior-sensitivity.R:432` — not touched.** It is F-049, MM's pending
statistical decision, and it still fails in the slow tier. See §4.

### F-050 — slow tier off the deprecated argument names (`afd68b56`, test-only)

All 1152 occurrences replaced, in the four files report 04 named, using the same
mapping `bgm()` itself applies to the deprecated arguments (`R/bgm.R:561-582`):

- `pairwise_scale = s` → `interaction_prior = cauchy_prior(scale = s)`
- `main_alpha = a, main_beta = b` → `threshold_prior = beta_prime_prior(alpha = a, beta = b)`

None of the migrated calls passed `interaction_prior` or `threshold_prior`
already, which is the condition under which `bgm()` honours the deprecated
argument, so the mapping is exact. Verified rather than argued — a GGM cell and
a mixed cell run both ways:

```
GGM pairwise draws           identical: TRUE
GGM interaction_prior        identical: TRUE
mixed pairwise draws         identical: TRUE
mixed main draws             identical: TRUE
mixed threshold_prior        identical: TRUE
```

Bitwise. The certified cells did not move.

### F-054 — version the correction-table cache key (`7d24331e`)

Key construction moved into `ggm_correction_table_key()` (mirroring
`zratio_surface_cache_key()`, which is the convention this one departed from)
and the package version inserted after the schema tag:

```
ggm_ctable_v1_0.2.0.0_q17_delta1.4166067_eta1_normal_shape1_g120_ns2000_nw500_sd3_gibbs.rds
```

Old-key files are never read again; nothing is migrated or deleted, per the
brief. Regression test: the key embeds `packageVersion("bgms")`, and a sentinel
table planted under a `0.0.0.0` key is not served — the build runs and writes
under the current key instead.

### F-055 — the q ≤ 3 hierarchical crash (`d6d804e2`)

Reproduced first against the rc1 install:

```
q=2 -> ERROR: replacement has 1 row, data has 0
q=3 -> ERROR: replacement has 1 row, data has 0
q=4 -> OK
```

Fixed on the shape the lead recommended, and I found nothing in the code arguing
otherwise — the opposite, in fact. A new predicate,

```r
zratio_anchor_grids_empty(cap)   # TRUE when either family's grid has no rows
```

guards `zratio_build_surfaces()` immediately after the size cap is resolved, so
a 0-row job table is never built and the surface is NULL. That routes the fit
through the engine's existing no-surface path, the additive saddle, and **the
additive saddle is exact there**: a mediating block excludes the toggled edge's
own two endpoints, so at a cap of 3 or less every component is a single node or
a single bridge — below `SurfaceFamily::size_min`, which is precisely the region
the engine already documents as "additive == exact"
(`src/models/ggm/zratio_engine.h:45-49`). Nothing is approximated away.

One thing the brief did not ask for, which I judged a defect I would otherwise
be introducing: the fence message has three branches, and with a NULL surface at
q ≤ 3 a verbose fit would have taken the fourth-wall `else` and told the user
*"the absolute-moment surface build failed"*, which is false. The message
function gains a `size` argument (default `NA`, so its existing one-argument
callers are untouched) and a fourth branch that states the honest thing. Both
the guard and the message read the same predicate, so they cannot drift.

Verified end to end on the branch:

```
z-ratio: at 2 variables no absolute-moment surface is built, and none is needed. …
q=2 -> OK | surface: NULL | bip: none
q=3 -> OK | surface: NULL | bip: none
q=4 -> OK | surface: built | bip: built
q=5 -> OK | surface: built | bip: built
```

Two regression tests: a unit test in `test-zratio-surface-build.R` (predicate,
NULL builds at caps 2 and 3, and the message wording) and an end-to-end smoke in
`test-bgm-hier-spec.R` running `bgm()` at q = 2, 3, 4, 5 and asserting the bip
surface is absent at 2/3 and present at 4/5.

### F-017 + hardening — `.Rbuildignore` (`50281f7a`)

Added `^tests/compliance$`, `^tests/fixtures$`, `^\.git$`; kept
`^tests/testthat/fixtures$`. Neither ignored directory is reachable from
`tests/testthat.R`: the compliance tier runs from
`.github/workflows/weekly-compliance.yaml` against a git checkout, and
`test-extractor-functions.R:829` skips when the legacy fixture directory is
absent. Tarball listing delta in §3.

### F-051 — the sensitivity plot's cap (`8e76b87e`)

A `@details` block on `plot.bgms_prior_sensitivity`, four sentences, no
computation and no unit touched (the log10 → nats question stays open):

> The drawn curves clamp the pooled inclusion probability at \eqn{1 - 10^{-6}},
> which caps a plotted \eqn{\log_{10}} Bayes factor at 6. A curve running flat
> along 6.0 has reached that display cap; it is not evidence levelling off. The
> uncapped value at the chosen scale is `x$edges$chosen_scale_log10_bf`, which
> `print()` reports and which can be far larger. `x$edges$saturated` does not
> mark capped curves: it records that the edge's inclusion indicator never
> flipped in the chain, which is a different condition.

`man/plot.bgms_prior_sensitivity.Rd` regenerated.

---

## 2. Findings

### 05-1 · `minor` · every line of `.Rbuildignore`, comments included, is a live regex

New. `R CMD build` and roxygen2's `package_files()` compile each non-blank line
of `.Rbuildignore` as a PCRE pattern — there is no comment syntax, only patterns
that happen not to match. My first draft of the F-017 hunk split an existing
comment across two lines, which left a parenthesis unbalanced, and
`roxygen2::roxygenise()` died outright:

```
Error in FUN(X[[i]], ...) :
  invalid regular expression '# ---- Legacy test fixtures (GitHub CI only, not'
In addition: Warning message:
  PCRE pattern compilation error 'missing closing parenthesis'
```

The shipped file already contains six comment lines with parentheses that
happen to balance, so the hazard is latent rather than active. Fixed here by
keeping each comment's parentheses balanced, and the constraint is recorded in
the commit message. Nothing else to do unless the lead wants the comments
stripped.

### 05-2 · `note` · the additive kernel is now only reachable below shape 0.5

Not a defect, but it is a consequence of PR #193 worth having on the record,
because it is what made the F-048 gauge repair non-trivial. After #193 the
engine has three routes and the additive-counts saddle serves exactly one band:
`alpha < 0.5`. Above 10 the route is the isolated-edge ratio, not additive;
between 0.5 and 10 it is the surface. So any future test that wants to exercise
the additive kernel end to end has to sit below 0.5 — where the Z-ratio
*constants* are themselves outside their certified range `[0.5, 20]` and warn.
The coarse kernel and the uncertified-constants band now coincide exactly. That
is fine for a detector test (the detector is the thing under test), but it means
there is no certified cell in which the additive kernel can be exercised.

### 05-3 · `note` (positive) · the harm channel is intact

Stated because report 04 could only infer it. Measured on this build, the harm
channel fires when it should (shape 0.4: `harm_pred` 0.043 > 0.01) and stays
down when it should (shape 2 under the surface: 0.00794 < 0.01), on the *same*
fixture. Both directions are now covered by tests.

---

## 3. Evidence — the verification gate

All four gates run from the clean install of the branch head.

### Gate 1 — changed test files pass

`NOT_CRAN=true BGMS_RUN_SLOW_TESTS=true`, from the full slow run below:

| file | tests | pass | fail | warn | skip | secs |
|---|---:|---:|---:|---:|---:|---:|
| `test-zratio-surface-build.R` | 10 | 37 | 0 | 0 | 0 | 21.7 |
| `test-zratio-gauge.R` | 11 | 56 | 0 | 0 | 0 | 14.1 |
| `test-correction-tables.R` | 11 | 97 | 0 | 0 | 0 | 0.5 |
| `test-bgm-hier-spec.R` | 13 | 55 | 0 | 0 | 0 | 27.1 |
| `test-sbc-ggm.R` | 6 | 12 | 0 | **0** | 0 | — |
| `test-parameter-recovery-ggm.R` | 2 | 32 | 0 | **0** | 0 | — |
| `test-mixed-nuts.R` | 11 | 264 | 0 | **0** | 0 | — |
| `test-scaling-diagnostics.R` | 9 | 45 | 0 | **0** | 0 | — |

(Report 04 measured 1000 / 100 / 38 / 14 warnings on those last four.)

### Gate 2 — the full slow tier

```
=== TOTALS (slow) ===
files   : 77
tests   : 1152
pass    : 8762
fail    : 1
error   : 0
warning : 0
skip    : 7
elapsed : 702.5 s
```

The single failure is the documented F-049, unchanged and untouched:

```
══ Failed ══════════════════════════════════════════════════════════════════
── 1. Failure ('test-prior-sensitivity.R:432:3'): the difference-scale reweighti
Expected `max(abs(rw$pip[2, pairwise] - pip_of(f2)[pairwise]))` < `4 * noise`.
Actual comparison: 0.0255 >= 0.0202
Difference: 0.0053 >= 0
```

Against report 04's slow run on rc1:

| | rc1 (report 04) | this branch |
|---|---:|---:|
| tests | 1129 | 1152 |
| pass | 8650 | **8762** |
| **fail** | **4** | **1** (F-049) |
| error | 0 | 0 |
| **warning** | **1152** | **0** |
| skip | 7 | 7 |
| wall clock | 710.4 s | 702.5 s |

The three F-048 failures are gone, the 1152 deprecation warnings are gone, and
the tier is otherwise unchanged. The +23 tests are 19 from the merged
checking-layer batch (`adf87013`) and 4 new ones here.

### Gate 3 — the full CRAN-mode suite

```
=== TOTALS (cran) ===
files   : 77
tests   : 1152
pass    : 7342
fail    : 0
error   : 0
warning : 0
skip    : 236
elapsed : 67.8 s
```

### Gate 4 — `R CMD check --as-cran`

Report 01's recipe verbatim, on a `git archive` export of the branch head
(`RSTUDIO_PANDOC` set, `_R_CHECK_CRAN_INCOMING_REMOTE_=false`):

```
* checking for hidden files and directories ... OK
* checking whether package ‘bgms’ can be installed ... [67s/68s] OK
* checking installed package size ... OK
* checking examples ... OK
* checking examples with --run-donttest ... [501s/256s] OK
* checking tests ...
  Running ‘testthat.R’ [93s/73s] OK
* checking re-building of vignette outputs ... [97s/52s] OK
* checking CRAN incoming feasibility ... NOTE
    The Date field is over a month old.
* checking HTML version of manual ... NOTE
    Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy.

Status: 2 NOTEs
```

**The same two NOTEs as report 01's baseline, and no others** — the stale `Date`
field (F-003, unrelated to this branch) and the environment's old HTML Tidy.
No new NOTE, WARNING or ERROR. `checking for hidden files and directories` is
OK, which is the `^\.git$` line doing nothing here (this tarball came from
`git archive`, which carries no `.git`); its effect is shown separately in
gate 4b.

### Gate 4b — the tarball listing

Full tarballs, both built with vignettes, so this is like-for-like.

**vs `develop` (a `git archive` of each, same build command):**

```
=== DELTA (only in the develop tarball) ===
  37 tests/compliance/…
   2 tests/fixtures/…
=== DELTA (only in the branch tarball) ===
(empty)
```

Exactly the 39 entries the item asked for, and **nothing else changed in either
direction**.

**vs report 01's rc1 baseline tarball** (`~/bgms-review/bgms_0.2.0.0.tar.gz`):
the same 39 removals, plus one addition — `tests/testthat/_snaps/plot-methods.md`
— which comes from the merged checking-layer batch on `develop`, not from this
branch (it is present in the `develop` tarball above, where the branch delta is
empty). Entry count 396 → 359. Size **4,524,133 → 1,582,541 bytes**.

**The `^\.git$` line**, checked empirically rather than by inspection, since a
worktree's `.git` is a file and so survives `R CMD build`'s own hidden-file
handling. Same export, same build, with a `.git` *file* planted at the root:

```
develop  + .git file -> tarball contains  bgms/.git
branch   + .git file -> tarball contains  (nothing)
```

---

## 4. Left failing on purpose

`test-prior-sensitivity.R:432` (F-049) fails in the slow tier and I did not
touch it, per the brief. It is deterministic and reproduces to the same digits;
report 04 §04-2 has the diagnosis (the gate's `4 * noise` threshold is built
from a 2-sample spread estimate, so the threshold itself has large sampling
variance). It is MM's statistical call whether the reweighting is off by ~0.005
pip or the tolerance is too tight.

**Consequence for the deadline:** the Monday 2026-08-03 03:00 UTC nightly runs
`main`, not this branch. If this batch is merged to `main` first, that nightly
goes red on exactly one test — F-049 — instead of four. If F-049 is expected to
stay open past Monday, the lead may want to say so in advance so the red is not
read as a new regression.

---

## 5. Open questions

1. **F-049 is still open and still red.** Nothing here changes that; it is
   flagged only because it is what the Monday nightly will show.

2. **Should there be a certified cell that exercises the additive kernel?**
   (finding 05-2) After #193 the additive band and the uncertified-constants
   band coincide exactly below shape 0.5, so the harm channel's end-to-end
   firing test necessarily runs in a cell whose constants are unscored. That is
   acceptable for a detector test and is documented in the test, but if the
   lead wants the additive kernel exercised inside the certified envelope it
   would need a deliberate hook (e.g. a test-only option that withholds the
   surface), which I did not add on my own authority.

3. **`.Rbuildignore` comments** (finding 05-1) are live regexes. Left as
   comments with balanced parentheses. If the lead prefers, they could be
   stripped or rewritten to be trivially safe.

4. **The fence message's fourth branch** was not in the brief. I added it
   because the alternative was a verbose q ≤ 3 fit being told its surface build
   failed, which is untrue. Flagging it as scope I took rather than scope I was
   given.

# easybgm and bgms 0.2.0.0 — what to change, and why

**For:** Nikola Sekulovski (easybgm)
**From:** the bgms pre-CRAN review, 2026-08-03
**Subject package:** easybgm 0.4.0 (CRAN, published 2026-04-02)
**Against:** bgms 0.2.0.0 (release candidate), with bgms 0.1.6.3 as the control

This document is self-contained. You should not need to read the bgms source
or the bgms review report to act on it. Every line number is easybgm 0.4.0's
CRAN source. Every claim about behaviour was executed, not inferred; the
measurements are at the end.

---

## 0. The short version

**Nothing is on fire.** easybgm 0.4.0 keeps working against bgms 0.2.0.0:

* your test suite passes **0 failed, 0 errors** (185 expectations);
* `R CMD check --as-cran` returns the **identical** status against old and new
  bgms — 1 WARNING, 2 NOTEs, byte-identical text, and all three are artifacts
  of re-checking a published tarball locally (version not incremented, stale
  timestamp, no `pandoc`), none caused by bgms;
* every easybgm workflow that runs against 0.1.6.3 also runs against 0.2.0.0.

That is not luck. bgms 0.2.0.0 moves its fit objects to **S7**, which would
break easybgm, and ships a **compatibility shim** that detects
`packageVersion("easybgm") < "0.5.0"` and hands you the old S3 list instead.
You will see one warning per fit saying so. The shim is temporary; bgms
intends to remove it a release or two after an S7-native easybgm is on CRAN.

So there is no deadline. There are, however, **four separate things** worth
your time, in this order:

| | what | why now |
|---|---|---|
| **A** | Make easybgm S7-native (§1) | 4 lines + a version bump. Turns off the shim and the per-fit warning. |
| **B** | Fix `MCSE_BF` (§2) | It now silently mixes two estimators and reports an interval that is valid for neither. This is the only *correctness* item. |
| **C** | Decide what to do about the numbers that moved (§3) | Your users' edge weights halve and their Bayes factors change. Mostly not your bug — but it is your documentation. |
| **D** | Pin `warmup` in two `.Rd` files and the test suite (§4) | Your `R CMD check` goes 83 s → 196 s. One argument each. |

Plus two pre-existing bugs unrelated to this release (§6), flagged because
you will be in those files anyway.

---

## 1. (A) Making easybgm S7-native — the required change

### What bgms did

From 0.2.0.0, `bgm()` and `bgmCompare()` return **S7 objects**:

```r
class(fit)   # c("bgms", "S7_object")
typeof(fit)  # "object"     <- NOT "list"
```

All the accessors you use still work on them — `fit$arguments`,
`fit$posterior_summary_pairwise`, `extract_arguments()`,
`extract_indicators()`, and the rest all dispatch correctly on an S7 object.

**Exactly one idiom breaks: assigning to `class(fit)`.** There are four such
lines in easybgm:

| file | line | statement |
|---|---|---|
| `R/functions.bgms.R` | 54 | `class(fit) <- "bgms"` |
| `R/functions.bgms.R` | 61 | `class(fit) <- "bgms"` |
| `R/functions.bgmscompare.R` | 169 | `class(fit) <- c("bgmCompare")` |
| `R/functions.bgmscompare.R` | 242 | `class(fit) <- c("bgmCompare")` |

Assigning strips `"S7_object"` from the class attribute. The assignment itself
**raises nothing** — but the object's base type is `object`, not `list`, so the
very next `.subset2()` inside bgms's `$` method fails. Executed, verbatim:

```
>>> class(trapped) <- "bgms"
    OUTCOME: NO ERROR (returned bgms)

>>> extract_arguments(trapped)                          ERROR: this object class is not subsettable
>>> trapped$arguments                                   ERROR: this object class is not subsettable
>>> extract_pairwise_interactions(trapped)              ERROR: this object class is not subsettable
>>> extract_posterior_inclusion_probabilities(trapped)  ERROR: this object class is not subsettable
>>> trapped$posterior_summary_pairwise                  ERROR: this object class is not subsettable
>>> bgm_extract.package_bgms(fit = <S7 fit>, ...)       ERROR: this object class is not subsettable
>>> plot_network(<S7 fit>)                              ERROR: this object class is not subsettable
```

The good news is that it fails **loudly** — a user cannot get a silently wrong
number out of this. The bad news is that the message names neither easybgm nor
the cause.

There is no `.subset2` anywhere in easybgm 0.4.0, so the class assignment is
the *only* S7 incompatibility. (bgms's own comment claims otherwise; that is a
bug in bgms's comment, reported separately.)

### E1 — `R/functions.bgms.R`, lines 50-62

Currently:

```r
  # --- Ensure proper bgms object and variable names ---
  if (!inherits(fit, "bgms")) {
    varnames <- fit$var_names
    fit <- fit$packagefit
    class(fit) <- "bgms"                                   # <- line 54: DELETE
  } else {
    varnames <- fit$arguments$data_columnnames             # <- line 56
    if (is.null(varnames)) {
      varnames <- paste0("V", 1:fit$arguments$no_variables)}
  }
  if(packageVersion("bgms") > "0.1.4.2"){
    class(fit) <- "bgms"                                   # <- line 61: DELETE (whole block)
  }
```

Replace with:

```r
  # --- Ensure proper bgms object and variable names ---
  if (!inherits(fit, "bgms")) {
    varnames <- fit$var_names
    fit <- fit$packagefit
  } else {
    fit_args <- bgms::extract_arguments(fit)
    varnames <- fit_args$data_columnnames
    if (is.null(varnames)) {
      varnames <- paste0("V", seq_len(fit_args$no_variables))
    }
  }
```

`fit$packagefit` already carries the correct class, so line 54 was never doing
anything you needed, and line 61 was a no-op even under 0.1.6.3. Switching
line 56 from `fit$arguments` to `bgms::extract_arguments(fit)` is not strictly
required — `$` works on S7 — but it is the documented accessor and it works on
both object shapes, which keeps the function version-agnostic.

`no_variables` is still present in a `bgm()` fit's argument list in 0.2.0.0
(bgms keeps it explicitly for you), so that line needs no change.

### E9 — `R/functions.bgmscompare.R`, lines 169 and 242

Delete both `class(fit) <- c("bgmCompare")` lines outright. Nothing replaces
them: `fit <- fit$packagefit` at line 89 already yields an object of class
`c("bgmCompare", "S7_object")`, and `extract_arguments(fit)` on the very next
line dispatches on it correctly.

### E2 — `DESCRIPTION`: bump `Version` to `0.5.0`

**This must ship in the same release as E1 and E9.** bgms's shim keys on
`packageVersion("easybgm") < "0.5.0"` and nothing else. Bumping the version
without the code changes converts a working warning into a hard failure;
making the code changes without the bump leaves the shim on, which is harmless
but pointless.

### E3 — `DESCRIPTION`: bump the bgms dependency (recommended)

Line 26 currently reads `bgms (>= 0.1.4)`. Once E1/E9 land, set
`bgms (>= 0.2.0.0)`.

Worth knowing: easybgm's highest bgms version gate anywhere in the package is
`packageVersion("bgms") > "0.1.6"`. A bgms 0.2.0.0 fit therefore falls into
exactly the same branch a 0.1.6.3 fit did — easybgm currently has no place to
notice that the release happened. Once you require `>= 0.2.0.0` you can delete
the now-dead branches at `functions.bgms.R:20-26`, `:47-49`, `:162` and
`functions.bgmscompare.R:7-10`, `:19-34`, `:103-159`.

---

## 2. (B) `MCSE_BF` now mixes two estimators — the correctness item

This is the one place where 0.2.0.0 makes easybgm produce a number that is not
right, rather than merely different.

### What changed underneath you

bgms 0.2.0.0 made **Rao-Blackwellized (RB)** inclusion inference canonical:

* `extract_posterior_inclusion_probabilities(fit)` now defaults to
  `estimator = "rb"` (`"raw"` is still available and reproduces the old
  numbers);
* the `n_eff` column of `fit$posterior_summary_indicator` is now the **RB**
  effective sample size. Under 0.1.6.3 that column held the transition-based
  ESS of the binary indicator chain.

The column name did not change. The meaning did.

### Why that breaks the interval

`R/functions.bgms.R:250-256`:

```r
    bgms_res$convergence_parameter <-  fit$posterior_summary_pairwise$Rhat
    bgms_res$MCSE_BF <-BF_MCSE(gamma_mat = extract_indicators(fit),
                               BF_vec = bgms_res$inc_BF[lower.tri(bgms_res$inc_BF)],
                               ess = fit$posterior_summary_indicator$n_eff,
                               return = "ci",
                               smooth_bf = FALSE)
```

Under bgms 0.1.6.3 all three inputs described the *same* estimator:
`gamma_mat` gave the raw `p_hat`, `BF_vec` came from the raw inclusion
probability, and `n_eff` was the transition ESS of that same chain.

Under 0.2.0.0:

* `BF_vec` is now **RB**-derived (it comes from `bgms_res$inc_BF`, which comes
  from the RB `inc_probs` at line 153);
* `ess` is now the **RB** effective sample size;
* but `BF_MCSE()` still recomputes `p_hat` from `gamma_mat` on the **raw**
  scale and uses it for the variance —
  `AuxiliaryFunctions.R:286` (`p_raw <- apply(gamma_mat, 2, mean)`) and
  `:322-326` (`p_for_variance <- p_raw`; `var_p <- p_for_variance * (1 - p_for_variance) / ess_vec`).

So the delta-method variance is built from a raw probability, an RB sample
size, and is attached to an RB Bayes factor. It is a valid MCSE for neither
estimator.

### E6 — the fix

Simplest version: give `BF_MCSE()` an optional `p_hat` that overrides
`p_raw`/`p_for_variance` when supplied, and pass the RB probability.

In `R/AuxiliaryFunctions.R`, `BF_MCSE()` (line 273) gains an argument:

```r
BF_MCSE <- function(gamma_mat,
                    BF_vec,
                    ess = NULL,
                    p_hat = NULL,          # <- new
                    smooth_bf = FALSE,
                    return = c("mcse_log", "mcse_bf", "ci"),
                    conf_level = 0.95) {
```

and at line ~322, where `p_for_variance` is set:

```r
  if (!is.null(p_hat)) {
    p_for_variance <- p_hat
  } else if (smooth_bf) {
    ...
```

Then in `R/functions.bgms.R:251-256`:

```r
    rb_p <- extract_posterior_inclusion_probabilities(fit)   # RB by default
    bgms_res$MCSE_BF <- BF_MCSE(
      gamma_mat = extract_indicators(fit),
      BF_vec    = bgms_res$inc_BF[lower.tri(bgms_res$inc_BF)],
      p_hat     = rb_p[lower.tri(rb_p)],
      ess       = fit$posterior_summary_indicator$n_eff,
      return    = "ci",
      smooth_bf = FALSE)
```

**Alternative, if you would rather not own the arithmetic:** bgms now computes
this on the RB scale throughout. `fit$posterior_summary_indicator` has columns

```
mean  mcse  sd  n_eff  Rhat  n0->1  n1->0
```

so `mcse` is the Monte Carlo standard error of the RB inclusion probability,
and `bgms::extract_inclusion_bf(fit, log = TRUE)` gives the log Bayes factor
on a scale that stays finite where the BF scale saturates. Building the
interval from those two is less code and needs no `gamma_mat` at all.

**Either way, handle `NA`.** bgms reports `NA` in `mcse`, `n_eff` and `Rhat`
where the RB draws are constant to double precision (i.e. the inclusion
probability sits at a numerical bound). Real example from a live fit:

```
                 mean         mcse           sd n_eff     Rhat n0->1 n1->0
intrusion-dreams    1           NA 3.662205e-13    NA 1.004189     0     0
intrusion-flash     1 1.572953e-08 2.224491e-07   200 1.003107     0     0
```

`BF_MCSE()` already degrades those to `NA` CI rows rather than erroring
(`AuxiliaryFunctions.R:357-358`), which is fine — but it is now a routine
occurrence rather than a pathology, so it belongs in the documentation for
`summary(..., BF_uncertainty = TRUE)`.

---

## 3. (C) The numbers that moved

None of these is an easybgm bug. They are listed because your users will
notice them and will email you, and because two of them make easybgm's own
documentation wrong.

### C-1. Every edge weight halves

bgms 0.2.0.0 reports pairwise parameters on the **association scale**: the
coefficient `omega` entering each conditional as `2 * omega * x`. Version
0.1.6.3 stored `sigma = 2 * omega`. So `extract_pairwise_interactions()` and
the raw posterior draws are **about half** their old values.

Reaches easybgm at `functions.bgms.R:146-148` and `:194-195`, and through
those:

* `res$parameters`
* `res$samples_posterior`
* `res$centrality` (strength = row sums of `abs(samples_posterior)`,
  `AuxiliaryFunctions.R:100-112`)
* every edge weight drawn by `plot_network()` and `plot_parameterHDI()`

Same on the comparison side (`functions.bgmscompare.R:207`, `:229-231`,
`:233`, `:282`, `:284`, `:300`, `:309`): `parameters`, `parameters_g1`,
`parameters_g2`, `group_estimates`, `overall_estimate`, `samples_posterior`.

**No code change needed.** But a note in `summary.easybgm`'s documentation or
`NEWS.md` will save you a lot of email.

### C-2. The default prior changed, and easybgm's docs now describe the old one

`bgm()`'s default interaction prior moved from **Cauchy(2.5)** to
**`normal_prior(scale = 1)`**, on the new coordinate. easybgm never passes
`interaction_scale`, so it inherits the new default: a plain
`easybgm(type = "binary", package = "bgms")` is now a **different model** than
it was.

That makes `R/easybgm.R:86-87` factually wrong:

> `interaction_scale` the scale of the Cauchy distribution that is used as a
> prior for the pairwise interaction parameters. The default is 2.5.

**E8 — recommended replacement** for the `\strong{bgms}` block at
`easybgm.R:82-110`:

* `interaction_prior` — prior on the pairwise interactions. Default
  `normal_prior(scale = 1)`. `cauchy_prior(scale = 2.5)` reproduces the
  pre-0.2.0.0 default *family*, but note the scale now acts on the association
  coordinate, so it is not the same model as the old default.
* `threshold_prior` — prior on the category thresholds. Default
  `beta_prime_prior(0.5, 0.5)`. (Your current text says `threshold_alpha` and
  `threshold_beta` default to 1; that has been wrong since before 0.2.0.0 —
  bgms has used 0.5/0.5 throughout.)
* `edge_prior` — now a prior *object*: `bernoulli_prior()`,
  `beta_bernoulli_prior()`, `sbm_prior()`. The old character strings
  (`"Bernoulli"`, `"Beta-Bernoulli"`, `"Stochastic-Block"`) still work through
  a deprecation branch, which is why your test at
  `tests/testthat/test-easybgm.R:66` still passes.

All the scalar prior arguments you document — `interaction_scale`,
`inclusion_probability`, `beta_bernoulli_alpha`/`beta`,
`beta_bernoulli_alpha_between`/`beta_between`, `dirichlet_alpha`,
`threshold_alpha`/`beta` — **still work**. Each raises one `lifecycle`
deprecation warning per session and is translated to the corresponding prior
object. Nothing you or your users pass will error.

On the comparison side, `bgmCompare()` moved **both** its baseline interaction
prior and the prior on the group differences — the focal parameters, the ones
the difference Bayes factors are about — from Cauchy to Normal. Verified on
both builds rather than taken from the release notes, because this is the
change most likely to move a published number:

* **bgms 0.1.6.3** has no `difference_family` argument. `difference_scale = 1`
  is documented as "Scale of the **Cauchy** prior for difference parameters"
  (`bgmCompare.R:94`), and the C++ gradient confirms it: the difference term is
  `-2 * value / (value^2 + difference_scale^2)`
  (`src/bgmCompare/bgmCompare_logp_and_grad.cpp:454`, `:474`, `:729`, `:747`),
  which is the derivative of a Cauchy log-density.
* **bgms 0.2.0.0** adds `difference_family = c("Normal", "Cauchy")`, so
  `match.arg()` makes **Normal** the default. Fitting `bgmCompare()` at all
  defaults and reading the spec back:

  ```
  bgmCompare() at ALL DEFAULTS
    difference_prior_type  : normal
    difference_scale       : 1
    interaction_prior_type : normal
    pairwise_scale         : 1
  ```

  i.e. **Normal(0, 1)** on the differences, where 0.1.6.3 had **Cauchy(0, 1)**.
  Passing `difference_family = "Cauchy"` restores the old family (verified:
  `difference_prior_type : cauchy`).

Practical consequence for easybgm users: a Normal slab has far lighter tails
than a Cauchy, so it bounds a weakly-identified difference much more tightly.
Difference verdicts move most on threshold differences for categories one group
never observed — exactly the cases C-5 is about. Nothing for easybgm to change;
worth knowing when someone reports that their `easybgm_compare()` Bayes factors
moved.

### C-3. Inclusion probabilities and Bayes factors change — mostly for the better

`extract_posterior_inclusion_probabilities()` now returns the RB estimate by
default (see §2). Downstream in easybgm: `res$inc_probs`, `res$inc_BF`
(`functions.bgms.R:155-178`, `:203-226`), and `res$structure = 1*(inc_probs > 0.5)`
at `:184` — which means an edge can cross the median-probability threshold
differently than before.

The most visible practical effect: a raw indicator average saturates at 0 or 1
on a short chain, which sent easybgm's hand-computed `inc_BF` to `0` or `Inf`.
The RB average does not saturate, so **your Bayes factors are now finite where
they used to be infinite.** That is an improvement, but it is a change.

**E5 — your decision.** Either adopt the new default (recommended; then also
do §2), or pin `estimator = "raw"` at lines 153 and 201 to reproduce 0.4.0's
numbers exactly. What you should not do is leave it as it is, because that is
the half-and-half state §2 describes.

**E7 — optional, and it would delete a lot of code.** Lines 154-181 and
202-229 reimplement the prior-odds division by hand, including the
Stochastic-Block within/between case. bgms 0.2.0.0 exports
`extract_inclusion_bf(fit, log = FALSE)`, which divides out the exact per-edge
prior inclusion odds on the accumulator's log scale, stays finite down to
log-acceptances of about −745, handles all three edge priors uniformly, and
returns a p × p matrix with the same dimnames you already assign:

```
             intrusion       dreams      flash    upset
intrusion           NA 4.039866e+13 7.0626e+07    0.082
dreams    4.039866e+13           NA 1.4700e-01 8549.499
flash     7.062600e+07 1.470000e-01         NA   29.925
upset     8.200000e-02 8.549499e+03 2.9925e+01       NA
```

Both blocks would collapse to
`bgms_res$inc_BF <- bgms::extract_inclusion_bf(fit)`, which also removes the
`stop("Unknown edge prior type.")` dead ends at lines 180 and 228.

**One gotcha if you do this:** bgms puts `NA` on the diagonal, where
`vector2matrix()` puts `0`. `summary.easybgm` only ever reads
`[lower.tri(...)]` so it is unaffected, but check any other consumer.

### C-4. R-hat moved

`res$convergence_parameter` (`functions.bgms.R:250`;
`functions.bgmscompare.R:235`, `:308`) is now the classic split-R-hat.
bgms 0.1.6.3 applied the Brooks-Gelman degrees-of-freedom adjustment, which
collapsed to a data-independent ~1.29 on nearly-saturated binary indicators —
i.e. on exactly the *most decisive* edges, which were then flagged as
unconverged. Those now report near 1. You may also see `NA` (every chain
identical) and `+Inf` (chains stuck constant at different values).

No code change needed; worth knowing before someone reports it as a bug.

### C-5. Comparison results move a lot when groups saw different categories

bgms 0.2.0.0 fixes a real defect in `bgmCompare()`. Version 0.1.6.3 kept only
the categories a variable was observed in **every** group and silently merged
the rest into their neighbours — so with groups observing categories 0-2 and
1-3 of a four-category item, both groups were folded onto two categories. The
affected variable's pairwise parameters were then badly overestimated, which
also inflated its neighbours'.

0.2.0.0 keeps every category any group observes. What reaches easybgm:

* `res$parameters`, `parameters_g1`, `parameters_g2` change, sometimes a lot;
* `num_categories` changes, so the `(main)`/`(diff)` column count of
  `extract_indicators()` changes, so the width of the `structures` strings you
  paste at `functions.bgmscompare.R:216` and `:295` changes — and with it
  `structure_probabilities`, `graph_weights`, and `sample_graph`;
* a `message()` about renumbering, and a **classed warning** of class
  `bgms_group_support_warning`, now reach the user.
  `easybgm_compare()`'s `tryCatch` at `easybgm_compare.R:184-193` catches
  errors only, so both pass through.

The classed warning is deliberate on bgms's side — you can catch that one
condition specifically without muffling everything else, e.g.

```r
withCallingHandlers(
  bgm_fit(...),
  bgms_group_support_warning = function(w) { ...; invokeRestart("muffleWarning") }
)
```

if you want to reformat it for easybgm users rather than pass it through raw.

**This does not fire in your own examples.** I checked: in
`ADHD[1:10, 1:3]` and `ADHD[11:20, 1:3]` every column observes both categories
in both groups, so `easybgm_compare.Rd` is unaffected. Real datasets with
unequal category support are where it will show up.

---

## 4. (D) Runtime — your check got 2.4× slower, and the fix is one argument

### The measurement

`R CMD check --as-cran` on easybgm 0.4.0, same machine, same easybgm, one job
at a time, only bgms differing:

| | bgms 0.1.6.3 | bgms 0.2.0.0 | ratio |
|---|---|---|---|
| **total wall** | **83 s** | **196 s** | **2.36×** |
| `checking examples with --run-donttest` | [78 s CPU / 25 s elapsed] | [252 s / 128 s] | 3.2× CPU |
| `checking tests` | below timing threshold | [31 s / 17 s] | — |
| your test suite, standalone | 4.6 s | 8.6 s | 1.87× |

Per example, from the checks' own `easybgm-Ex.timings` (elapsed seconds):

| example | backend | 0.1.6.3 | 0.2.0.0 | ratio |
|---|---|---|---|---|
| `complexity_probs` | bgms | 9.065 | 58.129 | **6.41** |
| `structure_probs` | bgms | 8.924 | 58.013 | **6.50** |
| `structure` | bgms | 0.485 | 1.894 | 3.91 |
| `centrality` | bgms | 0.557 | 1.896 | 3.40 |
| `HDI` | bgms | 0.613 | 1.984 | 3.24 |
| `easybgm_compare` | bgms | 1.090 | 2.472 | 2.27 |
| `edgeevidence` | **BGGM** | 0.597 | 0.486 | 0.81 |
| `network` | **BGGM** | 0.237 | 0.236 | 1.00 |
| `easybgm` | **BGGM** | 0.026 | 0.026 | 1.00 |

The three BGGM-backed examples are unchanged to within 1%, which is what makes
the rest attributable to bgms rather than to machine noise.

### Why

`bgm()`'s and `bgmCompare()`'s `warmup` default moved from **1e3 to 2e3**, and
NUTS `target_accept` from 0.60/0.65 to **0.80** (smaller steps ⇒ more leapfrog
steps per iteration). easybgm passes `iter` everywhere and **never passes
`warmup`**, so it inherits the doubling in full.

Two details worth knowing:

* `--as-cran` runs `\donttest` as a **second** examples pass on top of the
  normal one, so your seven bgms-reaching example fits are paid for twice per
  check.
* The cost scales with model size. On a 3-variable, 20-observation
  `easybgm_compare` the new bgms is actually *faster* (0.19 s → 0.07 s) because
  its backend rewrite more than pays for the extra warmup at that size. The
  ratio climbs with the model: ~1.7× at 5 variables, 3.68× at 17. So this
  matters for `complexity_probs.Rd` and `structure_probs.Rd`, which fit the
  full 17-variable Wenchuan, and barely at all elsewhere.

### E13 — the fix

**In `man/complexity_probs.Rd` and `man/structure_probs.Rd`** — these two are
~60% of the whole check:

```r
fit <- easybgm(data, type = "ordinal", save = TRUE, edge_selection = TRUE,
               iter = 100,      # for demonstration only
               warmup = 100)    # <- add: bgms now defaults to 2000
```

**In `tests/testthat/test-easybgm.R`**, add `warmup = 10` beside each
`iter = itr` for the bgms combos (lines 45-82 and 236-260).

Optionally the same in `man/HDI.Rd`, `man/centrality.Rd`, `man/structure.Rd`,
`man/easybgm_compare.Rd`, though those are small.

This is entirely under your control — bgms is not going to move the default
back — and it is worth doing before bgms 0.2.0.0 reaches CRAN, so the
reverse-dependency check sees the cheaper version.

---

## 5. New: bgms fits that easybgm has never seen

bgms 0.2.0.0 adds two model classes that did not exist in 0.1.6.3:

* **GGM** — `bgm(x, variable_type = "continuous")`, a Gaussian graphical model
  with edge selection;
* **mixed MRF** — `bgm(x, variable_type = c("ordinal", "continuous", ...))`.

`easybgm()` itself cannot produce these: it reroutes `type = "continuous"` and
`type = "mixed"` to BDgraph/BGGM with a warning (`easybgm.R:226-232`,
`easybgm_compare.R:167-173`). I confirmed that still happens.

**But your plotting methods are registered on class `bgms`** —
`plot_network.bgms`, `plot_edgeevidence.bgms`, `plot_structure.bgms`,
`plot_structure_probabilities.bgms`, `plot_complexity_probabilities.bgms`,
`plot_parameterHDI.bgms`, `plot_centrality.bgms`, `plot_centrality.list` — and
each funnels a raw bgms fit into `bgm_extract.package_bgms`
(`plottingfunctions.bgms.R:21-25`, `plottingfunctions.easybgm.R:667-672`,
`:778-783`). So a user who fits a GGM with bgms directly and then calls
`plot_network(fit)` lands in your code.

**It runs.** No error, no warning, real `qgraph`/`ggplot` output. That is the
problem — but not because it contradicts your documented "partial associations"
(`easybgm.R:23`). That phrase covers a pairwise coupling perfectly well, and
easybgm has always put the unbounded coupling in `parameters` for ordinal bgms
fits, where no standard normalization exists. The problem is that easybgm uses
a **different convention for Gaussian models**, and a bgms GGM is a Gaussian
model arriving down the discrete path:

| backend | what goes into `res$parameters` | site | this edge |
|---|---|---|---|
| BGGM | `out_select$pcor_mat` — partial correlations | `functions.bggm.R:44` | — |
| BDgraph | `qgraph::wi2net(fit$K_hat)`, i.e. `-cov2cor(K)` — partial correlations | `functions.bdgraph.R:91`, `:131` | **0.5563** |
| bgms, ordinal | pairwise interaction parameter (unbounded) | `functions.bgms.R:147`, `:195` | n/a |
| **bgms, GGM** (new) | **same as ordinal** — association coupling, `-0.5 × K_ij` | `functions.bgms.R:147` | **0.8434** |

Both of your Gaussian backends report the partial correlation. Applying
`qgraph::wi2net` to this fit's precision matrix gives **0.5563**, identical to
`bgms::extract_partial_correlations()`. The bgms GGM path gives **0.8434** for
the same edge of the same fit:

```
easybgm res$parameters (association coupling):  0.8434
qgraph::wi2net(K)  ==  extract_partial_correlations(fit):  0.5563
bgms::extract_precision(fit) off-diagonal:    -1.6868
   check:  -0.5 * -1.6868                   =  0.8434
           1.6868 / sqrt(3.0893 * 2.9761)   =  0.5563
```

Same model class, same slot, two scales, and nothing marks which is on the
plot. A user comparing a BGGM network against a bgms GGM network in one paper
would be comparing 0.5563 against 0.8434.

(`extract_category_thresholds()` returns `NULL` for a GGM, which is harmless
since nothing reads `bgms_res$thresholds`.)

**E14 — pick one:**

**(a) Guard.** At the top of `bgm_extract.package_bgms`, after you have the
arguments:

```r
  mt <- bgms::extract_arguments(fit)$model_type
  if (!is.null(mt) && mt %in% c("ggm", "mixed_mrf")) {
    stop("easybgm does not yet support bgms Gaussian or mixed graphical models. ",
         "Use the bgms package's own plot() and summary() methods for these fits.",
         call. = FALSE)
  }
```

**(b) Support it properly.** Branch on `model_type` and use
`bgms::extract_partial_correlations(fit)` for the GGM case — it returns the
p × p partial correlation matrix directly, with dimnames, and agrees with
`qgraph::wi2net` to the digit — so that a bgms GGM lands on the same scale as
your BGGM and BDgraph backends.

(a) is a few minutes and is honest. (b) is a feature. Doing neither leaves a
plausible-looking network on a scale the user did not ask for.

---

## 6. Two pre-existing bugs, unrelated to this release

Neither is caused by bgms 0.2.0.0 — both are equally wrong against 0.1.6.3 —
but you will be in these files.

### E10 — `R/functions.bgmscompare.R:189-194`: wrong argument names

```r
      } else if (args$difference_prior[1] == "Beta-Bernoulli") {
        edge.prior <- args$beta_bernoulli_alpha /
          (args$beta_bernoulli_alpha + args$beta_bernoulli_beta)
```

This one is easy to get wrong, and the reason is worth stating, because the
name easybgm uses is a real bgms name — just not the one the value is *stored*
under.

`beta_bernoulli_alpha` **is** a formal argument of `bgmCompare()` in bgms
0.1.6.0 through 0.1.6.3 (`bgmCompare.R:189`, documented at `:94`), and remains
a deprecated formal in 0.2.0.0. So reading `args$beta_bernoulli_alpha` looks
correct. But `extract_arguments()` on a compare fit has never returned it under
that name. Checked across every version easybgm's branches can reach:

| bgms | name in a `bgmCompare()` fit's `$arguments` |
|---|---|
| 0.1.4.2 | `pairwise_beta_bernoulli_alpha` / `main_beta_bernoulli_alpha` (prefixed) — `bgmCompare.R:311-314` |
| 0.1.6.0 | `difference_selection_alpha` / `difference_selection_beta` |
| 0.1.6.3 | `difference_selection_alpha` / `difference_selection_beta` |
| 0.2.0.0 | `difference_selection_alpha` / `difference_selection_beta` |

So the formal is `beta_bernoulli_alpha`, the stored field is
`difference_selection_alpha`, and easybgm reads the formal's name off the
stored list. Your own **multi-group** branch at lines 254-256 already uses the
stored names correctly — this is a rename from bgms 0.1.6.0 that got followed
through in one sub-branch and not the other. (The pre-0.1.6.0 branch at lines
116-117 has the same slip against 0.1.4.2's prefixed names, so it is
long-standing rather than new.)

As written, `edge.prior` is `numeric(0)` and `inc_BF` at line 214 comes back a
zero-length matrix.

Fix — copy the multi-group branch's names:

```r
      } else if (args$difference_prior[1] == "Beta-Bernoulli") {
        edge.prior <- args$difference_selection_alpha /
          (args$difference_selection_alpha + args$difference_selection_beta)
```

This only fires on a two-group comparison with a Beta-Bernoulli difference
prior, which is why the test suite does not catch it — the suite only exercises
the default Bernoulli prior.

### E12 — `R/plottingfunctions.easybgm.R:663`: a guard that is missing in one place

```r
      if(!fit_args$save){
        stop("Samples of the posterior distribution required ...")
      }
```

`save` **was** a stored field of a bgms fit's `$arguments` up to and including
bgms 0.1.4.2 (`bgm.R:466`). bgms 0.1.6.0 stopped storing it — samples became
unconditional — and 0.2.0.0 keeps `save` only as a deprecated formal that warns
if passed (`bgm.R:526`). So from 0.1.6.0 onward `fit_args$save` is `NULL` and
`if(!NULL)` raises `argument is of length zero`.

You already handle this: every sibling method sets the field before reading it —
`plottingfunctions.bgms.R:10-12`:

```r
  if(packageVersion("bgms") > "0.1.4.2"){
    fit_args$save <- TRUE
  }
```

`plot_centrality.list` is the one place the guard was not added. Add the same
two lines before line 663, and at `:776-783` where `save = fit_args$save` is
passed into `bgm_extract`. Nothing about bgms 0.2.0.0 changes this — it has
been broken since 0.1.6.0 — but it will keep being broken.

---

## 7. Suggested release plan

You are not blocked and nothing needs to be simultaneous. bgms will ship
0.2.0.0 first, with the shim carrying easybgm 0.4.0.

**Before bgms 0.2.0.0 reaches CRAN, if convenient:** §4 (E13) only — the two
`.Rd` files and the test suite. One argument each; it makes the
reverse-dependency check cheaper for everyone.

**easybgm 0.5.0, whenever suits you:**

- [ ] **E1** — delete `class(fit) <- "bgms"` at `functions.bgms.R:54` and `:61`
- [ ] **E9** — delete `class(fit) <- c("bgmCompare")` at `functions.bgmscompare.R:169` and `:242`
- [ ] **E2** — `DESCRIPTION: Version: 0.5.0` (same release as E1/E9, not before)
- [ ] **E3** — `DESCRIPTION: bgms (>= 0.2.0.0)`; drop the dead `< 0.1.6` branches
- [ ] **E6** — fix `MCSE_BF` (§2) — the correctness item
- [ ] **E5** — decide RB vs `estimator = "raw"` (§3, C-3)
- [ ] **E4** — `extract_category_thresholds()` → `bgms::extract_main_effects()` at `functions.bgms.R:149`, `:196` (removes a deprecation warning; same return shape for ordinal fits, `NULL` for GGM)
- [ ] **E8** — update the documented bgms prior defaults at `easybgm.R:82-110`
- [ ] **E14** — guard or support GGM/mixed fits (§5)
- [ ] **E10, E12** — the two pre-existing bugs (§6)
- [ ] **E7** — optional: replace the hand-rolled `inc_BF` with `bgms::extract_inclusion_bf()`

**After easybgm 0.5.0 has been on CRAN a while**, bgms removes the shim. There
is no hard date; bgms will not do it while a meaningful number of users are
still on 0.4.0.

---

## 8. How to reproduce any of this

Two private libraries, one per bgms build, sharing one easybgm:

```sh
mkdir -p ~/lib-new ~/lib-old
R CMD INSTALL --library=~/lib-new    bgms_0.2.0.0.tar.gz
R CMD INSTALL --library=~/lib-old    bgms_0.1.6.3.tar.gz
Rscript -e 'install.packages("easybgm", lib="~/lib-new", repos="https://cloud.r-project.org")'
```

Put `~/lib-old` *first* on `.libPaths()` for the 0.1.6.3 runs and omit it for
the 0.2.0.0 runs; easybgm resolves from `~/lib-new` either way.

```sh
# the paired check
R_LIBS=~/lib-old:~/lib-new R CMD check --as-cran --no-manual easybgm_0.4.0.tar.gz
R_LIBS=~/lib-new           R CMD check --as-cran --no-manual easybgm_0.4.0.tar.gz
```

To see the S7 failure for yourself, in a fresh session — note the order,
`library(easybgm)` **after** the fit:

```r
library(bgms)
fit <- bgm(na.omit(Wenchuan)[1:60, 1:4], iter = 100, warmup = 200, chains = 2)
class(fit)                  # "bgms" "S7_object"
library(easybgm)
plot_network(fit)           # Error: this object class is not subsettable
```

(If you load easybgm *first*, bgms's shim gives you an S3 list and everything
works — that asymmetry is the shim, not a bug in your code.)

---

## 9. Questions I could not answer from here

1. Whether you want easybgm's inclusion Bayes factor to stay independent of
   bgms's (E7). Adopting `extract_inclusion_bf()` changes your numbers a second
   time on top of C-3; doing both in one release is kinder to users than doing
   them a release apart, but it is your call.
2. Whether GGM/mixed support (E14b) is on your roadmap at all, or whether the
   guard (E14a) is the right permanent answer.
3. C-5's real-world reach. The category-union fix does not fire in easybgm's
   own examples, so nothing here exercises the largest comparison-side change
   end-to-end through easybgm. A two-group dataset with unequal category
   support would be the honest test and I did not have one.

Happy to send patches for any of the mechanical items (E1, E9, E10, E12, E13)
if that is useful.

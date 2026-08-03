# bgms 0.2.0.0

This release rebuilds most of the package on top of 0.1.6.3. `bgm()` now fits
Gaussian graphical models and mixed discrete-continuous networks alongside
ordinal ones; priors are supplied as objects instead of loose scalars;
inclusion inference is Rao-Blackwellized throughout; `bgmCompare()` has been
reworked; and a set of checking tools (`verdicts()`, `calibration_check()`,
`prior_sensitivity_check()`, network and edge plots) is new. Several defaults
changed, some arguments were removed, and pairwise parameters are reported on
a different scale, so a script written against 0.1.6.3 can run and give
different numbers. Those changes are listed first.

Throughout this file and in the package's output, "log" means the natural
logarithm.

## Breaking changes

**Pairwise parameters are reported on the association scale.**

* `bgm()` stores and reports pairwise interaction parameters on the
  association scale: the coefficient `omega` that enters every conditional as
  `2 * omega * x`. In 0.1.6.3 the same coupling was stored as `sigma = 2 *
  omega`. Reported pairwise effects — `coef()`, `summary()`,
  `extract_pairwise_interactions()`, and the raw posterior samples — are
  therefore about half their 0.1.6.3 values for the same data and model. Main
  effects, indicators, and inclusion probabilities are unaffected in scale.
  Code that reads raw pairwise samples, or that compares numbers against a
  0.1.6.3 run, needs the factor of two applied; there is no automatic
  conversion.

* `bgmCompare()` moves to the association scale with `bgm()`, in the same
  release and by the same factor: the two entry points were mutually consistent
  in 0.1.6.3 and remain so. Pairwise effects reported by `coef()`, `summary()`,
  `extract_pairwise_interactions()`, `extract_group_params()`, and the raw
  posterior samples are therefore about half their 0.1.6.3 values, exactly as
  for `bgm()`. **A stored `bgmCompare()` fit from an earlier version should be
  refit rather than rescaled**, because the prior acted on the old, doubled
  coordinate: the same `difference_scale` or `interaction_prior` number is now a
  tighter prior than it was, and difference inclusion Bayes factors move
  accordingly. The calibration of the `difference_scale` default under the
  association-scale parameterization is still under study, so difference
  verdicts near a decision threshold should be read as scale-contingent.

* `simulate_mrf()` reads its `pairwise` argument on the association scale too:
  the discrete Gibbs conditional now accumulates `2 * (score - baseline) *
  pairwise` where 0.1.6.3 accumulated `(score - baseline) * pairwise`. **A
  `pairwise` matrix that a 0.1.6.3 script passed in now means twice the
  coupling it used to**, so simulated data are more strongly dependent than
  before at the same input. Halve a hand-specified `pairwise` matrix to
  reproduce 0.1.6.3 output; a matrix taken from a 0.2.0.0 `bgm()` fit needs no
  adjustment, because the fit is on the same scale as the simulator. The
  continuous (`variable_type = "continuous"`) branch takes a precision matrix
  and is unaffected.

**Changed defaults.**

* The default interaction prior is `normal_prior(scale = 1)` in `bgm()`, where
  0.1.6.3 used a Cauchy at `pairwise_scale = 2.5`; `sample_ggm_prior()` uses
  the same default, so the GGM prior chain matches `bgm()`. `bgmCompare()`
  takes the same default for its baseline pairwise interactions, so a compare
  fit and separate `bgm()` fits at stated defaults now price a baseline
  interaction alike. Pass `interaction_prior = cauchy_prior(scale = 2.5)` to `bgm()` for
  the 0.1.6.3 prior, remembering that the scale now acts on the association
  coordinate. Under the joint precision-graph specification the normalizer
  correction table is keyed on the interaction prior, so the first fit at the
  new default cell (with a Beta-Bernoulli or Stochastic-Block edge prior)
  builds and caches a fresh table once.

* `bgmCompare()` now defaults to `interaction_prior = normal_prior(scale = 1)`
  for the baseline pairwise interactions and to `difference_family = "Normal"`
  for the group differences, both of which were Cauchy in 0.1.6.3 and through
  0.2.0.0's development; the baseline default now matches `bgm()`, so the two
  entry points no longer ship different priors for an identically named
  argument. Results under the defaults change: a fit at the new defaults is a
  different model from a fit at the old ones, and difference verdicts on
  categories one group never observed move the most, because a Normal slab
  bounds an unidentified threshold difference far more tightly than a Cauchy
  one. `cauchy_prior(scale = 1)` and `difference_family = "Cauchy"` remain
  fully available and reproduce the previous behaviour.

* `iter` and `warmup` default to `2e3` in `bgm()` and `bgmCompare()`, up from
  `1e3`. A default call therefore runs about twice as long and returns twice as
  many draws as it did in 0.1.6.3. Set `iter` and `warmup` explicitly to
  recover the old budget.

* The NUTS `target_accept` default is `0.80` in both `bgm()` (was `0.60`) and
  `bgmCompare()` (was `0.65`). Higher target acceptance means smaller steps:
  expect more leapfrog steps and a longer run per iteration, in exchange for
  fewer divergences. The adaptive-Metropolis default is unchanged at `0.44`.
  Pass `target_accept` to restore the old value.

* On continuous and mixed models the prior on the precision diagonal defaults
  to `exponential_prior(eta = 1)`, where `eta` is the rate in the frame in
  which the pairwise slab has unit scale; the raw rate is derived at fit time
  as `eta / s` for interaction-prior scale `s`. At the default unit slab scale
  this is the same distribution as `gamma_prior(shape = 1, rate = 1)`, but with
  a non-unit slab scale the diagonal rate now co-scales. One consequence: a
  scale-free interaction prior (`beta_prime_prior()`) on a continuous model
  requires the diagonal prior to be given with `rate` rather than `eta`.

**Removed and reshaped arguments.**

* `update_method = "hamiltonian-mc"` has been removed; use
  `update_method = "nuts"`. NUTS adapts trajectory length dynamically and is
  more reliable, especially with edge selection on GGM models. The
  `hmc_num_leapfrogs` argument has been removed with it. A call passing either
  now fails. On continuous data the new conjugate `update_method = "gibbs"`
  (see New features) is a second option that needs no step-size tuning at all.

* The `standardize` argument of `bgm()` and `bgmCompare()` no longer does
  anything: pairwise interactions are on the association scale and share one
  prior scale, so the per-pair adjustment by the product of the two variables'
  maximum scores has been dropped. The argument stays as a deprecated formal so
  the failure is informative rather than an unused-argument error.
  `standardize = FALSE`, the old default, warns and proceeds, because it is
  what the sampler already does; `standardize = TRUE` errors and points to
  setting the scale directly through `interaction_prior` (and
  `difference_scale` for `bgmCompare()`). There is no per-pair equivalent, and
  none is planned: standardization belongs in the model rather than in an
  argument, and a forthcoming g-prior formulation of the interaction prior
  addresses it there. The per-pair maximum-score adjustment is retired in
  favour of that direction.

* Priors are now objects. The scalar arguments `pairwise_scale`, `main_alpha`,
  `main_beta`, `inclusion_probability`, `beta_bernoulli_alpha`,
  `beta_bernoulli_beta`, `beta_bernoulli_alpha_between`,
  `beta_bernoulli_beta_between`, `dirichlet_alpha`, and `lambda` are deprecated
  formals: passing one warns and is translated into the corresponding prior
  object, so old calls keep working for now. `edge_prior` takes a prior object
  (`bernoulli_prior()`, `beta_bernoulli_prior()`, `sbm_prior()`); the old
  character values are still accepted through a deprecation branch. The
  replacements are `interaction_prior`, `threshold_prior`, `means_prior`,
  `precision_scale_prior`, and `edge_prior`. Note that `pairwise_scale = s`
  translates to `cauchy_prior(scale = s)`, preserving the 0.1.6.3 family, so a
  deprecated call and a call at the new default are not the same model.

* `extract_category_thresholds()` is deprecated in favour of
  `extract_main_effects()`, which covers category thresholds, continuous means,
  and precision diagonal entries in one call.

**Changed output and semantics.**

* `extract_ess()` reports the Rao-Blackwellized effective sample size for the
  edge (or difference) indicators — the `n_eff` column of the fit summary's
  inclusion table — where 0.1.6.3 returned the indicator chain's
  transition-based ESS. A 0.1.6.3-era call therefore gets different numbers
  with no warning. The extractor family is now coherent: the Rao-Blackwellized
  inclusion probability (`extract_posterior_inclusion_probabilities()`), its
  R-hat (`extract_rhat()`), and its ESS all describe the same estimate. Fits
  made with bgms < 0.2.0.0 carry no Rao-Blackwellized draws, so a default call
  on such a fit falls back to the transition ESS without warning, while an
  explicit `estimator = "rb"` errors.

* The Rao-Blackwellized estimate is the canonical inclusion probability
  throughout a fit: `extract_posterior_inclusion_probabilities()` defaults to
  `estimator = "rb"` (with `"raw"` still available) and
  `fit$posterior_mean_indicator` reports it. Numbers quoted from a 0.1.6.3 fit
  will not reproduce exactly; the two estimate the same quantity, the
  Rao-Blackwellized one with lower variance.

* The fit summary's inclusion table has different columns. It now reports the
  Rao-Blackwellized inclusion probability with a continuous Monte Carlo
  standard error (`mcse`), effective sample size (`n_eff`), and split-R-hat
  (`Rhat`), beside the per-direction transition counts (`n0->1`, `n1->0`). The
  transition-based `n_eff_mixt` column is gone (see Deprecated). The three
  Rao-Blackwellized columns are `NA` where the Rao-Blackwellized draws are
  constant to double precision, which places the inclusion probability at its
  numerical bound, and are computed everywhere else.

* On a `bgmCompare()` fit the Rao-Blackwellized inclusion probabilities are
  `NA` for indicators the sampler never updated — the main-effect differences,
  unless `main_difference_selection = TRUE`. **A previously all-numeric matrix
  can now contain `NA`**, which breaks downstream code that assumes otherwise.
  The `NA` is honest, and the contrast with `bgm()` is the point: an
  always-included edge indicator in `bgm()` is still *proposed* at every
  iteration, so the Rao-Blackwellized accumulators see its conditional
  inclusion odds and report a finite Bayes factor. An unselected main-difference
  indicator is never proposed at all, so there is no Rao-Blackwellized quantity
  to report — the entry is `NA` exactly when the indicator was never visited,
  not because a number was lost. Use `estimator = "raw"` for the old
  all-numeric behaviour, or filter on `is.na()`.

* R-hat is the classic split-R-hat (Gelman et al. 2013 / Stan;
  `Rhat = sqrt(var_plus / W)` with `var_plus = (n-1)/n * W + B/n`). 0.1.6.3
  applied the Brooks-Gelman degrees-of-freedom adjustment used by
  `coda::gelman.diag`, which collapsed to a data-independent `sqrt(5/3) ~ 1.29`
  on nearly-saturated binary edge (or difference) indicators — the typical
  shape of a decisive edge, where one split sub-chain carries a brief excursion
  and the rest are constant — so the most decisive edges were flagged as
  unconverged. Reported R-hat changes slightly for all parameters (the
  adjustment factor was near 1 for healthy chains) and substantially for
  near-saturated indicators, which now report values near 1, or `NA` when every
  chain is identical. Chains stuck constant at different values report `+Inf`
  instead of `NA`. See the diagnostics vignette for how to read indicator R-hat
  alongside the transition counts and the inclusion probability.

* Fitted objects from `bgm()` and `bgmCompare()` are S7 class objects (new
  dependency: `S7`). All existing `$`, `[[`, and `names()` access patterns
  continue to work. When an incompatible `easybgm` version is loaded, bgms
  returns plain S3 lists for backwards compatibility; this shim will be removed
  in a future release. The compatibility mode is chosen when the fit object is
  constructed, so load `easybgm` before fitting when the fit is meant for use
  with an older `easybgm`.

* `summary()` on a `bgms` fit returns an additional `quadratic` element (the
  precision diagonal, reported on the residual-variance scale for continuous
  and mixed fits) and carries a `main_label` attribute.

## New features

**New model classes.**

* Gaussian graphical models: `bgm(x, variable_type = "continuous")` fits a GGM
  with Bayesian edge selection. Sampling uses NUTS on a free-element Cholesky
  (theta-space) parameterization of the precision matrix, which keeps the
  precision matrix positive-definite by construction; adaptive-metropolis is
  also available.

* A conjugate GGM Gibbs sampler: `update_method = "gibbs"` fits a GGM with a
  row-by-row update of the precision matrix, with or without edge selection. It
  needs no step-size or proposal tuning and supports a Normal or Cauchy prior
  on the edges. Continuous data only.

* With `update_method = "gibbs"` and edge selection, the first 15% of warmup
  runs the full model so the precision matrix settles, and edge selection is
  active for the remaining 85%. Both windows scale with the warmup budget.

* Mixed MRF models: `bgm()` accepts a per-variable `variable_type` vector that
  mixes `"ordinal"`, `"blume-capel"`, and `"continuous"` types, to estimate
  networks with both discrete and continuous variables. `simulate.bgms()` and
  `predict.bgms()` support mixed models.

* Missing-data imputation now covers the new model classes:
  `na_action = "impute"` integrates over missing values during MCMC sampling
  for continuous and mixed models as well as ordinal ones.

**Prior specification.**

* Priors are constructed by dedicated functions and passed as objects:
  `normal_prior()`, `cauchy_prior()`, and `beta_prime_prior()` for parameter
  priors; `gamma_prior()` and `exponential_prior()` for the precision diagonal;
  `bernoulli_prior()`, `beta_bernoulli_prior()`, and `sbm_prior()` for the edge
  prior. Each has a `print()` method. The corresponding `bgm()` formals are
  `interaction_prior` (pairwise slab), `threshold_prior` (category thresholds,
  default `beta_prime_prior(0.5, 0.5)` — the 0.1.6.3 `main_alpha`/`main_beta`
  pair), `means_prior` (continuous means, default `normal_prior(scale = 1)`),
  `precision_scale_prior` (precision diagonal), and `edge_prior`.

* `gamma_prior()` and `exponential_prior()` accept the rate as `eta`, the rate
  on the precision diagonal in the frame where the pairwise slab prior has unit
  scale. The raw rate is derived at fit time as `eta / s` for interaction-prior
  scale `s`, so the prior geometry is held fixed when the slab scale changes;
  at fixed `eta`, graph and partial-correlation inference is invariant to the
  slab scale. Applies to GGM and mixed MRF models.

* `bgm(delta =)` sets the exponent of the determinant tilt on the continuous
  precision block. It defaults to `NULL`, resolved at fit time to
  `0.5 * log(p)` over the continuous variables.

* `bgmCompare(difference_family =)` chooses the family of the prior on the
  group differences, `"Normal"` (the default) or `"Cauchy"` (0.1.6.3's fixed
  behaviour). It governs the pairwise-interaction differences and the
  main-effect threshold differences alike, and under
  `difference_selection = TRUE` it is the slab of the spike-and-slab.

**Inclusion inference.**

* Rao-Blackwellized inclusion probabilities and Bayes factors, through
  `extract_posterior_inclusion_probabilities(estimator = "rb")` and
  `extract_inclusion_bf()`. For each indicator update the sampler records the
  one-step draw `J = gamma + (1 - 2 gamma) alpha` (pre-move state `gamma`,
  birth/death acceptance probability `alpha`), stored per chain and per edge in
  `fit$raw_samples$rb_inclusion` alongside the raw indicator draws. Averaging
  `J` is a lower-variance estimator of the inclusion probability than the raw
  indicator average (`estimator = "raw"`). The sampler also accumulates the
  inclusion odds on the acceptance-probability scale, so `extract_inclusion_bf()`
  stays finite down to log-acceptances of about -745: edges whose raw average
  saturates at 0 or 1 because the chain never flipped in a finite run still
  receive a finite, data-driven Bayes factor. The prior inclusion odds are
  divided out (from `extract_prior_inclusion_probabilities()` for `bgm()`, from
  the exchangeable difference prior for `bgmCompare()`), so the reported value
  is a Bayes factor rather than posterior odds, and equals the posterior odds
  only at a prior inclusion probability of 1/2. Because the Rao-Blackwellized
  draw is continuous, the standard MCSE/ESS/split-R-hat machinery applies to
  inclusion inference. Available for `bgm()` (ordinal, continuous, and mixed)
  and `bgmCompare()` (difference indicators). The estimator inherits the
  chain's mixing and does not rescue a chain that failed to explore the model
  space.

* `extract_inclusion_bf(log =)` selects the scale. The default `log = FALSE`
  returns the inclusion Bayes factor itself, matching the function's name;
  `log = TRUE` returns its natural logarithm. The accumulators are exact on the
  log scale everywhere, while the Bayes-factor scale saturates at double
  precision: an entry whose log exceeds about 709.78 nats (a Bayes factor
  beyond about 1.8e308) is reported as `+Inf` under `log = FALSE` even though
  its log-scale value is finite, and an entry whose log is `-Inf` becomes `0`.
  Workflows that must separate such extreme evidence should use `log = TRUE`.
  The argument applies to both the `bgm()` and the `bgmCompare()` method,
  including the stochastic-block difference path, whose return remains posterior
  odds.

* `extract_prior_inclusion_probabilities()` returns prior edge-inclusion
  probabilities in the same matrix layout as
  `extract_posterior_inclusion_probabilities()`, for prior/posterior
  inclusion-odds computations. Under the joint spike-and-slab prior on a
  continuous block the graph marginal is reweighted by the per-graph normalizer
  (positive-definite-cone mass shaped by the determinant tilt), so
  continuous-continuous edges do not keep the edge-prior marginal at any
  `delta` — with a uniform Beta-Bernoulli hyperprior at 3 variables the prior
  edge probability is about 0.37 rather than 0.5, and a fixed
  `bernoulli_prior(0.5)` at `delta = 0` lands near 0.27. The values are read
  from the cached correction table (`bernoulli_prior()`,
  `beta_bernoulli_prior()`) or estimated by a prior-only chain with the fit's
  own correction settings (`sbm_prior()`, cached on the fit). Mixed models
  report per-edge-class values (discrete-discrete, continuous-continuous,
  cross); ordinal models use the analytic edge-prior marginals, including the
  exchangeable partition mixture for `sbm_prior()`.

**Checking a fitted model.**

* `verdicts()` reads each edge (or difference) indicator's inclusion Bayes
  factor as a three-way verdict — evidence of presence, undecided, evidence of
  absence, at a threshold and its reciprocal — and flags the verdicts that a
  rerun of the sampler could change. The flag is the union of two standard
  errors of the logit inclusion probability: the Jeffreys-smoothed two-state
  model of the binary indicator chain, which stays defined when the indicator
  never flips, and the delta-method transform of the Rao-Blackwellized Monte
  Carlo standard error. An edge is fragile when either places a verdict
  boundary within two standard errors of the estimated evidence. The rule comes
  from a known-truth calibration study of 37,010 graded edge-fits across
  ordinal, binary, and Gaussian graphical models, in which every one of the 66
  verdict errors sat within 0.58 of a threshold on the log Bayes factor scale
  and no edge further out was ever misclassified: over those fits the two-state
  standard error alone recalled 0.742 of the errors and the Rao-Blackwellized
  one 0.939, while the union recalled all of them at a 3.0% false-alarm rate
  and transferred across model types. Classification goes through the same
  internal routine as `prior_sensitivity_check()`, so the two classify
  identically. Printing reports the verdict tally and, when any verdict is
  fragile, says so. Every arm of that study fitted a single network, so the
  operating point belongs to the edge indicators of `bgm()`; on `bgmCompare()`
  difference indicators the flag still marks verdicts near a boundary, but its
  error-catching rate and false-alarm rate there are unmeasured, and the print
  method says so rather than borrowing the single-network numbers.

* `calibration_check()` fits the isotonic (pool-adjacent-violators) reliability
  curve of the model's conditional predictions, one panel per variable: pooling
  across variables lets a variable predicted too high and one predicted too low
  cancel, and the pooled curve then tracks the diagonal while neither variable
  does. The consistency band resamples each case's category from its own
  predicted distribution rather than resampling the cumulative threshold events
  independently, because those events are nested within a case and treating them
  as independent understates the band. Continuous variables are covered through
  the probability integral transform: `u = F(y | rest)` is uniform exactly when
  the conditional predictive distribution is right, so the panel is the
  empirical distribution function of the `u` values against the uniform
  diagonal. `F` is the predictive mixture over `ndraws` posterior draws rather
  than a plug-in Gaussian at the posterior mean, so parameter uncertainty sits
  inside the distribution the observation is transformed by; the two panel kinds
  therefore differ in what they condition on, and the `kind` column of both
  returned tables records which construction produced each row. The continuous
  band is not resampled from the model, because under the transform the null is
  `Uniform(0, 1)` whatever the conditional density was: it is the simultaneous
  envelope of the empirical distribution functions of `n` independent uniforms.
  Both panel kinds live on the unit square with the diagonal as the calibrated
  reference, so a mixed fit produces one figure and one summary table.
  Evaluated on the fitted data by default; `newdata` makes it out-of-sample.
  In-sample the observation also entered the parameters it is judged against,
  which makes the transform mildly under-dispersed and the check conservative.
  A `bgmCompare()` method is available and reports per-group panels.

* `prior_sensitivity_check()` recovers the continuous inclusion-Bayes-factor
  curve of each edge across the interaction slab scale, and classifies how
  every edge's verdict behaves along it. The curve is anchored at a handful of
  fixed-scale fits (log-spaced `anchors`, default `0.4 / 0.63 / 1 / 1.6 / 2.5`
  times the chosen scale) and filled in between anchors by importance
  reweighting: a fixed-scale fit is reweighted to a nearby scale with per-draw
  slab-density ratios over the included edges, the likelihood cancels, and no
  normalizing constant is involved. The `1x` anchor is the original fit itself
  and is never refit, so the chosen-scale verdicts the check reports are exactly
  the analysis already run. Each display point pools every anchor that clears
  `ess_floor` (default 400) there, weighting each anchor's inclusion-probability
  estimate by its inverse variance on the probability scale and transforming the
  pooled probability to the Bayes-factor scale; points where no anchor clears
  the floor are reported `NA` rather than extrapolated, and non-overlapping
  anchor radii warn to add anchors. Exactness is kept off the pooled curve and
  on the anchor fits themselves: the per-anchor verdict columns and the
  chosen-scale quantities come from those fits' own Rao-Blackwellized
  statistics, not from a pooled row. An edge is called scale-sensitive only if
  its verdict changes along the curve *and* its Bayes-factor swing exceeds a
  run-to-run noise band (the 95th percentile of the spread between one anchor
  refit and a repeat of it), so boundary edges that merely flip between reruns
  of the same prior are reported as `indistinguishable-from-wobble` rather than
  as moves; an edge that is threshold-relevant in one of the two refits and
  saturated in the other has a censored rather than an infinite spread and is
  left out of that percentile, with the printed report counting the exclusions.
  Each anchor passes a convergence gate (median continuous split-R-hat, bulk
  Rao-Blackwellized inclusion R-hat, and E-BFMI for NUTS) before it enters the
  curve. The check works on any `bgm()` fit with edge selection; on Wenchuan it
  costs about one original fit with warm-started NUTS refits (a
  `refit_sampler` argument, default `"same-as-fit"`, carries the adapted NUTS
  step size and mass matrix).

* `prior_sensitivity_check(vary =)` names which prior the curve moves, and the
  printed report says which. On a continuous or mixed fit the interaction slab
  scale and the prior on the precision diagonal are tied through the
  standardized frame (raw rate `= eta / s`), so a sweep of the slab has to hold
  one of the two fixed, and the two choices answer different questions:
  `vary = "slab"` holds the raw diagonal rate at the fitted value and moves the
  interaction prior alone, while `vary = "slab-and-diagonal"` holds `eta` fixed
  and lets the raw rate follow, rescaling the prior's overall size with its
  shape held fixed. The default `vary = "auto"` follows the frame the fit itself
  used — the joint sweep when the diagonal prior was given as `eta`, the slab
  alone when it was given as a raw `rate`. On a five-variable Gaussian graphical
  model the two modes give visibly different curves, so the choice is not a
  formality. Discrete fits have no precision diagonal; `vary` is accepted, has
  no effect there, and the report omits the line.

* `prior_sensitivity_check()` works on `bgmCompare()` fits, tracing each
  difference indicator's inclusion Bayes factor across the difference slab
  scale. `bgmCompare()` gives the pairwise and the main-effect difference
  families one `difference_scale`, so a single curve covers both. The anchored
  machinery, the convergence gate, the noise band, and the reporting are the
  same as for `bgm()`; what differs is that one difference indicator gates
  several parameters (a pairwise difference in every group contrast, a
  main-effect difference across that variable's whole threshold block), so the
  importance weight sums the difference-slab density ratio over the gated
  parameters while the inclusion probability is still reported per indicator.
  Reweighting a fit at one scale to another was checked against real refits at
  that scale: on Wenchuan at a doubling of the scale the reweighted inclusion
  probabilities land within 0.016 of a refit's, against a 0.007 run-to-run
  spread between two refits at that same scale, at an importance effective
  sample size of 2536 out of 12000 draws. Indicators the sampler never updated
  carry no verdict at any scale, and the printed report counts them out rather
  than reporting them as undecided. A stochastic-block difference prior has no
  single marginal inclusion probability and is refused, since the curve would
  report posterior odds rather than Bayes factors.

* New vignette "Checking your fitted model": `verdicts()` and the
  verdict-encoded network, then `calibration_check()`, then a worked example of
  building a predictive display directly on `simulate()` — the joint-level check
  is an analyst's reading of a display, not a pass/fail statistic, so the
  vignette shows how to roll one rather than calling a blessed function. The
  worked statistic is the sum-score distribution against its pointwise
  predictive band, chosen because it lies outside the model's sufficient set:
  pairwise dependence statistics sit close to that set, so how they read depends
  on how the estimator is anchored rather than on whether the model fits.

**Plots and centrality.**

* `plot()` on a `bgm()` fit draws the **edge evidence plot**: three panels on
  one shared layout, holding the pairs the data support, the pairs the data
  rule out, and the pairs the data cannot decide. A single network drawing has
  to make every pair either an edge or a blank, and a blank cannot say which of
  those two a missing edge is; because `bgm()` returns an inclusion Bayes
  factor for every pair, that choice does not have to be made. Each panel is
  titled with what it holds and how many pairs are in it, and with the rule
  that put them there, stated as a Bayes factor rather than its logarithm. Only
  the first panel is weighted: line width is the posterior mean pairwise
  association and colour carries its sign, because that is where the effect
  sizes are; the colour convention is documented in `?plot.bgms`. The other two
  are drawn at uniform width, dashed and dotted, because for those pairs the
  classification is the result. The layout is computed once from every pair, so
  a node sits in the same place in all three panels, and an all-absence fit is
  a result rather than an error: the absence panel fills and the presence panel
  is empty. The threshold is the one `verdicts()` uses, so the picture and the
  reporting table state the same thing. A fit run without edge selection has no
  inclusion Bayes factor and nothing to split its pairs by, so `plot()` draws
  one weighted panel titled `Edge weights` — the title says which channel the
  figure is drawn in, so a wide line is never ambiguous between a large
  association and strong evidence for one. `type = "centrality"` draws the
  centrality display instead. Layout comes from qgraph, which stays a suggested
  package; without it the method errors and points to `verdicts()`.

* `plot()` on a `bgmCompare()` fit draws the same three panels for difference
  evidence, with the supported panel weighted by the posterior mean difference.
  The display draws for any number of groups: a pair has a single inclusion
  indicator shared across all `K - 1` contrasts, so the three-way split is as
  well defined for three groups as for two; what a pair does not have beyond
  two groups is one magnitude, so the supported panel is weighted only for two
  groups and drawn at uniform width otherwise, where the classification is the
  whole of what it reports. A fit with no differences left is a result rather
  than an error — "the groups do not differ anywhere" is a common and correct
  finding. Main-effect differences are not edges, and when
  `main_difference_selection = TRUE` gives them their own indicators their
  evidence rides on the nodes: each node wears a ring filled to that difference
  indicator's posterior inclusion probability — a full ring is 1, half a ring
  0.5. Under the default `main_difference_selection = FALSE` those indicators
  are never updated and no ring is drawn. `type = "groups"` draws each group's
  own posterior mean network on the same layout and pages them, so a fit with
  more groups than `max_panels` is drawn a page at a time rather than squeezed
  into one row; those panels use the same colour-vision-safe sign pair as the
  rest of the package rather than qgraph's green/red default. With
  `difference_selection = FALSE`, weights are the display: two groups get one
  network of every pair, width the posterior mean difference and colour its
  sign, titled `Difference weights`; beyond two groups the weights of a pair
  are `K - 1` numbers whose honest weighted picture is the groups themselves,
  so `plot()` draws the `type = "groups"` panels without being asked for it.
  Read the magnitudes per group with `plot(fit, type = "groups")` or
  `extract_group_params()`.
  `type = "centrality"` draws the centrality display, with `group` passed
  through.

* `plot_edge_posterior()` draws one edge the way JASP draws a parameter: the
  posterior of the edge weight against the prior it was updated from, on a
  numbered density axis. Under edge selection the posterior drawn is the weight
  *given the edge is included* — a real density integrating to one, which is
  what puts numbers on the y axis — and the mass at zero moves off the plot
  region onto a probability wheel filled to the posterior inclusion probability,
  matching the ring encoding the `bgmCompare()` panels use for their main-effect
  differences. The prior curve is computed in closed form from the fit's own
  `interaction_prior` rather than assumed, so a `cauchy_prior(scale = 2.5)` fit
  and a fit at the default `normal_prior(scale = 1)` do not look alike. The
  median, the 95% credible interval, and the natural-log inclusion Bayes factor
  are printed as annotations in JASP's positions.

  What the panel shows follows the estimator the fit licenses. With edge
  selection, the evidence is the Rao-Blackwellized indicator Bayes factor, which
  is not a ratio of densities at zero, so no Savage-Dickey ordinates are drawn —
  dots there would assert an estimator the package does not use. Without edge
  selection there is no indicator, the posterior of the weight is continuous,
  and Savage-Dickey *is* the licensed estimator: the panel then draws the JASP
  figure exactly, with both grey ordinate dots at zero, the ratio printed as a
  log Bayes factor, and the wheel filled to `BF/(1 + BF)`; what the two shares
  are is documented in `?plot_edge_posterior` rather than tagged on the figure.
  The posterior ordinate is estimated by
  `stats::density(bw = "SJ")` at zero; JASP's own implementations use a
  logspline fit, which is not adopted here because it would add a dependency for
  one number. On a continuous block the drawn prior is the slab the sampler
  evaluates, which is the marginal prior of a precision entry only up to the
  positive-definite restriction; `?plot_edge_posterior` says so. An edge that
  the data rule out gets a figure rather than an error: the prior, the wheel
  (nearly all pale) and the evidence are drawn.
  `plot_edge_posterior()` refuses a `bgmCompare()` fit
  rather than drawing its baseline pairwise effect as though it were one
  network's edge; use `plot(fit)` and `verdicts(fit)` there.

* `extract_centrality()` evaluates a node centrality on every posterior draw of
  the network and returns a draws-by-nodes matrix, so the posterior distribution
  of a node's centrality comes out of the same model-averaged posterior as the
  edge weights, with structural uncertainty included: an edge excluded at a
  given iteration contributes exactly zero there. `summary()` reports the
  posterior mean, credible interval, and the posterior probability of being the
  most central node; `plot()` draws the ordered interval display.
  `measure = "strength"` is the one accepted value for now.

* `extract_centrality()` works on `bgmCompare()` fits. `group = 1` gives that
  group's centrality, one row per posterior draw; `group = c(1, 2)` gives the
  difference between two groups' centralities, again per draw, so the credible
  interval is the interval of the difference and answers whether the groups
  differ in a node's centrality directly, where two separately drawn intervals
  do not. Each group's network is rebuilt on every draw as
  `baseline + (P %*% differences)` with the fit's own contrast projection, not
  from posterior means, which is what carries the uncertainty through; averaged
  over draws it reproduces `extract_group_params()` to machine precision. A draw
  in which a difference indicator is excluded gives both groups the same edge
  weight and contributes exactly zero, so a difference carries a point mass at
  zero and `summary()` reports `p_positive` (the probability the first group's
  centrality is higher) rather than `p_most_central`, with the remaining mass on
  no difference at all. `plot()` marks zero.

**New extractors.**

* `extract_main_effects()`: main effect samples — category thresholds,
  continuous means, and the precision diagonal — replacing
  `extract_category_thresholds()`.
* `extract_precision()`: posterior precision matrix samples from GGM and mixed
  models.
* `extract_partial_correlations()`: posterior partial correlation samples from
  GGM and mixed models.
* `extract_log_odds()`: log-odds for discrete pairwise interactions.

**Prior-only sampling.**

* `sample_graph_prior()` draws edge-inclusion indicators, together with any
  edge-prior hyperparameters, from the graph level of the spike-and-slab prior.
  Under the hierarchical specification the draw is ancestral and exact; under
  the joint specification it runs the zero-data prior chain, so the draws carry
  the per-graph normalizer tilt. Hyperparameters can be fixed instead of sampled
  (`theta` for Bernoulli and Beta-Bernoulli priors, `allocations` plus
  `block_probs` for the Stochastic-Block prior).

* `sample_sbm_prior()` gives ancestral draws of block allocations and
  pair-inclusion probabilities from the MFM-SBM edge-prior hyperprior.

* `sample_ggm_prior()` runs a prior-only GGM chain. It accepts
  `update_method = "gibbs"` for the `spec = "joint"` prior chain, using the
  conjugate row and edge updates instead of adaptive Metropolis; `edge_prior`
  objects (`bernoulli_prior()`, `beta_bernoulli_prior()`, `sbm_prior()`),
  returning sampled allocations under the block model; and
  `spec = "hierarchical"`.

**Edge priors on continuous and mixed models.**

* Corrected inclusion-probability updates for `beta_bernoulli_prior()` on
  continuous (GGM) models: under the determinant-tilted precision prior, the
  conjugate Beta update omits a normalizing-constant factor and biases the
  sampled inclusion probability toward sparsity (at 5 variables with a uniform
  hyperprior its prior mean lands near 0.37 instead of 0.5). The update now
  draws from the corrected conditional using a table built from the prior
  distribution at the first fit of a model configuration and cached on disk
  (`tools::R_user_dir("bgms", "cache")`; a one-time cost of the order of
  minutes, announced when `verbose = TRUE`). The sampled inclusion probability
  is returned per chain in `fit$inclusion_parameter_samples`.

* Corrected stochastic block model updates for `sbm_prior()` on continuous (GGM)
  models, using the same cached table: the block-probability draws, the
  block-allocation weights, and the new-cluster weight all carry the
  normalizing-constant correction. Without it the sampled partition collapses
  toward one block (at 20 variables the prior mean number of blocks lands near
  1.0 instead of 1.87); with it a prior-only chain reproduces the model's
  partition prior (total variation 0.01-0.02 at 5-20 variables).

* The normalizing-constant correction extends to mixed models with
  `beta_bernoulli_prior()` or `sbm_prior()`: the determinant tilt acts on the
  continuous precision block, so the correction table is built for the
  continuous variables and the block-structure corrections read
  continuous-continuous edges only. With fewer than two continuous variables no
  edge is tilted and the plain conjugate updates apply unchanged; with exactly
  two, the single tilted pair supports the inclusion-probability correction but
  not the block-model slope curve, so `sbm_prior()` warns and keeps the plain
  conjugate updates there. Prior-only mixed chains reproduce the Beta hyperprior
  on the inclusion probability and the partition prior on the number of blocks.

* The correction table announces itself only when it actually builds (a cache
  hit is silent) and states that the build is one-time and cached. In
  interactive sessions it draws a progress bar for both serial and parallel
  builds. The bar follows `display_progress` (the same control as the sampler's
  bar), not the advisory `bgms.verbose` flag.

* Under `precision_graph_prior = "joint"`, a fit with edge
  selection on a continuous precision block reports once that the realized
  edge-inclusion prior is `pi(Gamma) * Z(Gamma)`, the edge prior reweighted by
  the per-graph normalizer, rather than the nominal edge prior. This holds for
  every edge prior: at three variables a uniform `beta_bernoulli_prior(1, 1)`
  realizes about 0.37 and a fixed `bernoulli_prior(0.5)` at `delta = 0` realizes
  about 0.27. With a learned inclusion probability (`beta_bernoulli_prior()`,
  `sbm_prior()`) the corrected hyperparameter update keeps the update coherent
  with the joint model but does not restore the nominal prior; with a fixed
  inclusion probability nothing absorbs the tilt. The message follows `verbose`
  (and the advisory `bgms.verbose` flag) and points to
  `extract_prior_inclusion_probabilities()` for the realized prior and to
  `precision_graph_prior = "hierarchical"` for the nominal one.

**The hierarchical precision-graph prior (continuous and mixed models).**

* `precision_graph_prior = "hierarchical"` is the default on continuous and
  mixed models. It composes the edge prior and the precision prior as
  `p(Gamma) p(K | Gamma)` with `p(K | Gamma)` normalized per graph, so the
  graph marginal is exactly the edge prior; under
  `precision_graph_prior = "joint"` it is that prior reweighted by the
  per-graph normalizer. Both specifications are fully supported; the default
  is the one whose graph marginal is the prior you wrote down. Each edge move evaluates the normalizer ratio with a fast
  per-component surface approximation, built once per analysis at the fixed
  `(eta, delta)` from block-Gibbs anchors. It requires `edge_selection = TRUE`
  and a `normal_prior()` or `cauchy_prior()` interaction prior; the precision
  scale prior may be `exponential_prior()` or a `gamma_prior()` of any shape. On
  mixed models the specification normalizes the continuous block
  `p(K_yy | Gamma_yy)`: the approximation enters the continuous-continuous edge
  moves only, and at least two continuous variables are required. The same
  specification is available in `sample_ggm_prior(spec = "hierarchical")`. In
  simulation checks with data (50 variables, n = 25, Beta-Bernoulli(2, 4)) the
  posterior inclusion probabilities are unbiased to within +0.002. Chains
  sampled from the prior alone, with a sampled inclusion probability and many
  variables, can fall outside the approximation's validated range; the trust
  gauge flags affected chains.

* The argument is accepted where the choice between the two specifications is
  vacuous rather than erroring. The specifications differ only in how
  `p(K | Gamma)` is normalized across graphs, so they differ only where the
  sampler moves between graphs. With `edge_selection = FALSE` the graph is fixed
  and the two coincide exactly; the argument is accepted silently, and the
  correction surface and the trust gauge are skipped rather than built and left
  unused. With no continuous precision block — an ordinal model, or mixed data
  with fewer than two continuous variables — there is no `K` for the argument to
  refer to; the fit is accepted, and a message reports this when the argument
  was named and `verbose = TRUE`. Inheriting the default is silent: since
  `"hierarchical"` is the default, the value reaches every fit, and the
  message reports a user's choice rather than a default nobody made. Both
  cases run the joint path, which applies no correction
  on exactly these configurations either, so the posterior is identical to a
  `"joint"` fit from the same seed. A `beta_prime_prior()` interaction slab is
  still rejected where the choice is meaningful — a continuous precision block
  under edge selection — and is tolerated where the argument has no referent.

* The per-edge normalizer ratio is corrected by a theta-independent
  absolute-moment surface. The surface predicts the per-component first and
  second spectral moments for both mediating-block families (common-neighbour
  clusters and bipartite bridges), and is deployed by decomposing each block
  into disjoint components. Its constants are built at unit slab scale, where
  the normalizer ratio is invariant to the slab scale, so results do not depend
  on the user's scale choice, and a `cauchy_prior()` slab gets its own constants
  rather than reusing the Normal ones. Validated against a per-component gold
  reference across a graph-density sweep: about 0.003 nats to gold, against
  about 0.036 for the additive kernel it replaces, up to about 1000x tighter in
  the dense regime. The surface deploys at every standardized rate `eta`, with
  gold gates recorded at `eta = 1` and `eta = 2`.

* Where the surface deploys, and what serves the cells it does not:

  - `gamma_prior()` diagonal shapes from 0.5 to 10 are served by the surface. It
    is scored against a block-Gibbs reference at shapes 0.5, 1, 2, 3 and 5 —
    every interior cell inside the 0.003-nat envelope, against 0.016 to 0.40
    nats for the additive kernel — and the interior of that range is
    interpolated, not measured. Anchors at a non-unit shape are run at a raised
    sweep budget, resolved by matching the measured anchor spread to the shape-1
    reference.
  - Above shape 10 the mediating correction is switched off and every edge is
    served the isolated-edge normalizer ratio, which is exact for an edge with
    no mediating structure, so the entire error of the route is the mediation it
    drops: at most 0.00028 nats over the scored band at shape 10 and 0.00022 at
    shape 12 (the largest value measured above 10, including a 100-variable
    block), against the 0.003-nat envelope the surface's own accuracy is stated
    in. The decay above shape 10 is not pointwise monotone, so the claim rests
    on the measured band maximum. The bound was measured at standardized rates
    (`eta`) of 1 and 2; past 2 the same route deploys and a note says plainly
    that the bound does not cover the cell.
  - Below shape 0.5, or when a surface build fails inside the scored range, the
    additive fallback kernel serves. That kernel returns exactly zero on
    common-neighbour mediating blocks above a size the cell fixes, discarding
    the whole normalizer ratio rather than approximating it coarsely: the
    common-neighbour edge constant is a difference (the clique-of-two moment net
    of its two endpoint nodes) and is negative in most cells, while the edge
    count grows as the square of the block size against a linear node count, so
    the moment sum crosses zero at a size the density and the cell set. This is
    a property of the kernel, not of the diagonal shape — it occurs at the
    exponential shape too, at a common-neighbour block of 28 variables at a
    standardized rate of 2. Bipartite blocks are unaffected. Measured at the
    default `delta` across 6 to 120 variables, the boundary runs from about 12
    to 32 common-neighbour variables; the collapse is out of reach below roughly
    16 variables and becomes reachable from about 20 to 30 upward. A fit on
    which it happens says so, driven by the counters below, and reports the
    largest common-neighbour block involved. `?summarize_zratio_gauge`
    documents the mechanism and the measurement.

* The anchored size range reaches 80 variables, so the reachable mediating block
  of an ordinary large fit sits inside it rather than in the extension zone; an
  anchor tier is placed at the size cap itself whenever the surviving top tier
  is more than two sizes below it. Past the anchored range the surface continues
  along its own boundary slope in log-size rather than freezing at the range's
  edge. Scored against a block-Gibbs reference on single-component blocks of 90
  to 150 variables against an 80-variable anchored range, the extension holds a
  median of 0.0006 nats and at most 0.0011 for common-neighbour blocks, and a
  median of 0.0043 and at most 0.0060 for bipartite ones, where freezing reaches
  0.060 and 0.095 and the additive fallback reaches 0.21 and 0.31. Inside the
  anchored range the surface is at 0.003 or better, so the extension is bounded
  and much better than the alternatives rather than validated to the in-range
  envelope, and the bipartite bound in particular sits about twice outside it.
  The absolute moments grow with component size, so a fitted slope pointing
  downward at the range edge is a fit pathology; it is floored at zero and
  counted in the per-chain `n_slope_floor` diagnostic.

* The certified range of the tabulated constants is [0.5, 20]; the fixed
  quadrature grids were scored at shapes 2, 10, 12, 15 and 20 at standardized
  rates 1 and 2 against references sharing no quadrature with them, worst
  deviation 8.6e-08 against a 1e-07 to 1e-06 working tolerance. A shape above
  20 warns.

* The surface build is cached (session and disk, keyed on the fit cell and size
  cap, following `options(bgms.correction_table_cache)`), so repeated fits of
  the same configuration reuse it. A first build uses the fit's own `cores` (the
  chains have not launched yet; the result is identical for any core count) and
  can be overridden with `options(bgms.zratio_surface_cores)`; on Windows a
  socket cluster stands in for forking on large builds. It announces itself when
  it runs, reporting the size cap and the number of cores it uses, and reports
  the elapsed time when it finishes — this accounts for the pause before the
  chains start on a first hierarchical fit. A cache hit is silent, and the
  announcement follows the advisory `bgms.verbose` flag. The build costs about
  24 s serial and 7 s on four cores at the full 80-variable range at the
  exponential shape (43 to 47 s serial at a non-unit shape); a fit on fewer
  variables sizes the build to its own variable count.

* A trust gauge runs by default on the deployed hierarchical path — which is
  to say on a default continuous or mixed fit with edge selection — with
  `options(bgms.zratio_gauge_sweeps = 0L)` as the off switch;
  `sample_ggm_prior()` keeps its own `zratio_diagnostics` argument. While the
  sampler runs, the gauge redoes a subset of each chain's edge add/remove
  decisions with the exact calculation and records `flip_rate`, the fraction of
  decisions that would have come out differently; a chain whose flip rate
  exceeds 1% is flagged. A second channel, `harm_pred`, projects the measured
  approximation error onto the inclusion-probability scale using the chain's own
  edge sensitivities and the edge-prior feedback amplification, and flags when
  the projected distortion exceeds 0.02 — this catches a consistent error that
  shifts the recovered network without changing individual edge decisions, which
  `flip_rate` cannot see at chains whose decisions are far from their
  accept/reject boundaries. It is computed for Bernoulli and Beta-Bernoulli edge
  priors. The gauge is silent when clean and prints only flagged chains; the
  per-chain summary is attached as `fit$zratio_diag` (and as
  `$zratio_diagnostics` on `sample_ggm_prior()` output). Its cost is fixed per
  chain rather than proportional to `iter`: two assessment sweeps audit a capped
  number of edge moves, which is nothing on a sparse posterior and about 5-10
  seconds per chain on a dense posterior at 100-200 variables. That is
  negligible against a production-length fit and noticeable on a short
  exploratory one, which is what the off switch is for. See
  `summarize_zratio_gauge()` and the diagnostics vignette, which documents the
  audit as a sample: it resolves coherent error, which is what `harm_pred`
  projects and what shifts a recovered network, and does not resolve rare
  edge-specific failures.

* The advice printed for a flagged chain is a ladder. It reports the measured
  error share, the mediating-block sizes the audit covered (`block_lo`,
  `block_hi`) and how many of the chain's non-trivial moves were audited; then
  advises raising `options(bgms.zratio_gauge_sweeps)` and refitting to resolve
  whether the signal is real; and only then mentions
  `precision_graph_prior = "joint"`, stating that it targets a different model,
  whose graph marginal is the edge prior reweighted by the per-graph normalizer
  rather than the edge prior itself. The gauge never switches specification on
  its own.

* Per-chain diagnostics in `fit$zratio_diag` carry `flip_rate`, `amplification`,
  `harm_pred`, `harm_flag`, the `se_mcse`/`se_se` uncertainty of the error
  estimate, the audited block sizes (`block_lo`, `block_hi`), the retained-sweep
  extrapolation counters (`n_pred_retained`, `n_extrap_retained`,
  `max_extrap_size_retained`), `n_slope_floor`, `n_isolated` (edge evaluations
  taking the isolated-edge route above shape 10) and `n_collapsed` /
  `max_collapse_size` (the additive kernel's zero-collapse, counted at the
  moment-sign condition itself and excluding trust-gauge sweeps). Post-fit notes
  are driven by those counters, so they report what the chains did rather than
  what the specification intended, and they reach a fit run with
  `verbose = FALSE`. The extrapolation notice separates warmup from retained
  sweeps — the sampler initializes from a complete graph, which puts every
  mediating block at roughly the number of variables for the first sweeps — and
  says plainly when no retained sweep extrapolated at all; trust-gauge sweeps
  run after sampling and enter neither tally.

## Other changes

* NUTS uses Stan's multinomial candidate weighting (log-sum-exp of `H0 - h` per
  leaf, biased progressive sampling at the top level) in place of the
  Hoffman-Gelman slice variable. The two schemes target the same posterior; the
  multinomial variant produces lower-variance candidate selection and has been
  Stan's default since 2017. User-facing output is unchanged apart from the new
  `accept_prob` diagnostic.

* NUTS Stage-2 warmup windowing matches Stan's
  `windowed_adaptation::compute_next_window`: when the window after the next
  would overshoot the Stage-3a boundary, the current next window is stretched to
  absorb the remaining Stage-2 budget instead of emitting a small trailing
  window. This eliminates a disruptive mass-matrix update and step-size reinit
  at the end of warmup and improves dual-averaging convergence.

* NUTS diagnostics include the per-iteration mean Metropolis acceptance
  probability (`fit$nuts_diag$accept_prob`, paralleling Stan's
  `accept_stat__`) and a per-chain `mean_accept_prob` summary. The diagnostics
  summary also prints the `warmup_incomplete` flag (energy not stationary) it
  already computed.

* `bgm(progress_callback =)` and `bgmCompare(progress_callback =)` take a
  function called as sampling proceeds, so an embedding application (JASP) can
  drive its own progress display.

* Dropped `coda` from Imports; ESS and R-hat are computed in C++ with on-demand
  (lazy) evaluation, replacing the eager R-based computation of 0.1.6.3. `$` and
  `[[` accessors on fitted objects trigger that computation on first access.

* Refactored the C++ backend: unified model hierarchy
  (`BaseModel` -> `GGMModel` / `OMRFModel` / `MixedMRFModel`), shared NUTS
  infrastructure, and fused log-posterior and gradient computation.

* On Windows with `RcppParallel` >= 6.0.0 (which bundles oneTBB 2022), a
  multi-core run and a single-core run started from the same `seed` are no
  longer bit-for-bit identical: the second and later chains follow different
  sampling trajectories under the new scheduler. Each chain still samples the
  same posterior — posterior summaries agree to within Monte Carlo error — so
  results are statistically equivalent, only not bitwise reproducible across
  core counts. Other platforms and earlier `RcppParallel` are unaffected.

* `extract_arguments()` now reports `main_effect_indices` for `bgmCompare()`
  fits: the per-variable column blocks of the baseline main-effect parameters.
  The layout was computed during fitting but kept internal, so there was no
  supported way to map main-effect parameters back to variables.

## Bug fixes

* Fixed the handling of ordinal categories in `bgmCompare()`, a defect present
  in 0.1.6.3. The comparison kept only the categories a variable was observed
  in *every* group, and silently merged the rest into their neighbours: with
  two groups observing categories 0-2 and 1-3 of a four-category item, both
  groups' data were folded onto two categories, however many observations each
  category carried. Because groups being compared routinely differ in the
  categories they use, this destroyed well-observed categories rather than
  unused ones — a five-level variable could be reduced to a binary one — and
  the affected variable's pairwise parameters were then badly overestimated,
  which also inflated the parameters of the variables it connects to.
  `bgmCompare()` now keeps every category that any group observes and merges
  none of them; only a category value that no group observes at all is
  dropped, with the remaining categories renumbered contiguously and a
  `message()` saying so. When a kept category is empty in some group, that
  group's threshold for it rests on the prior rather than on data, so
  `bgmCompare()` warns and names each variable, category, and group affected,
  and records the per-group category counts in the fit
  (`extract_arguments(fit)$category_support`). Blume-Capel variables were
  never affected and are still exempt, because their parameters are functions
  of the numeric category score. **Any earlier group comparison in which the
  groups did not observe exactly the same categories should be refit**; the
  pairwise estimates it produced may be substantially too large. Comparisons
  whose groups shared the same observed categories are unchanged.

* Fixed the number-of-blocks summary for stochastic-block edge priors, a
  mismatch a 0.1.6.3 user would have seen. The conditional p(K | t) behind
  `posterior_num_blocks` placed a zero-truncated Poisson prior on the number of
  components, while the sampler's own partition coefficients use the shifted
  Poisson (K - 1 ~ Poisson(lambda)). The mismatch reweighted the reported
  distribution by a factor lambda/K across K, so the reported number of blocks
  was systematically off from the one the sampler drew under. The summary now
  uses the sampler's convention and is unit-tested against the generative
  partition prior.

* Fixed a compilation failure on Alpine/musl that also affected 0.1.6.3:
  `mrf_simulation.cpp` used `tbb::global_control` while relying on a transitive
  include for `<tbb/global_control.h>`, which is not available on all platforms.
  The header is now included directly.

* Fixed a `target_accept` that was partly ignored, in 0.1.6.3 as well. The value
  reached dual averaging, but the step-size heuristic that reruns after a mass
  matrix update used a hard-coded `0.625` instead, so a fit set to a different
  `target_accept` was restarted at the wrong target every time the metric was
  re-estimated. The heuristic now uses the requested target.

* Fixed NUTS acceptance-probability accumulation, present in 0.1.6.3 too: the
  top-level trajectory loop overwrote the Metropolis contribution with the last
  subtree's value instead of summing across the full trajectory, so the signal
  driving dual-averaging step-size adaptation was biased and the adapted step
  size with it.

* `predict()` now recodes `newdata` through the category map stored with the
  fit rather than by subtracting each column's minimum. Under the old rule a
  sparse category coding — one with gaps, such as values 1, 2, 4, 5 — was
  miscoded silently: the fit collapses those values to categories 0, 1, 2, 3,
  but the minimum shift mapped them to 0, 1, 3, 4, so every value above a gap
  was attributed to the wrong category and the largest value fell outside the
  fitted category range entirely, with no error and no warning. `simulate()`
  correspondingly returns ordinal data on the original category scale, so the
  `simulate()` to `predict()` round trip stays on the scale the model was
  fitted on.

## Deprecated

* `bgm(standardize =)` and `bgmCompare(standardize =)` — see Breaking changes;
  `FALSE` warns, `TRUE` errors.
* The scalar prior arguments of `bgm()` and `bgmCompare()` — `pairwise_scale`,
  `main_alpha`, `main_beta`, `inclusion_probability`, `beta_bernoulli_*`,
  `dirichlet_alpha`, `lambda` — in favour of the prior objects
  `interaction_prior`, `threshold_prior`, `means_prior`, and `edge_prior`.
  Passing one warns and is translated.
* `bgmCompare(difference_probability =)` in favour of
  `difference_prior = bernoulli_prior()`; the argument keeps no default.
* `extract_category_thresholds()` in favour of `extract_main_effects()`.
* `extract_ess(estimator = "mixt")`, the transition-based effective sample size
  of the binary indicator chain: it warns and recomputes the quantity from the
  raw indicator draws, and the `n_eff_mixt` column is gone from the fit
  summary's inclusion table. Three grounds. It carries no decision value: in a
  known-truth calibration study its best error-catching threshold rejected 41%
  of correct edge verdicts to catch 85% of the wrong ones, while a rule built on
  the inclusion standard error caught all of them at a 2.3% cost. Its own model
  is violated: the conversion from flip counts to an effective sample size
  assumes a two-state first-order Markov chain, and the study rejects
  first-orderness for 41.3% of edge-fits on transitions out of the included
  state (3.2% out of the excluded state), asymmetrically. And it is redundant:
  the precision of the inclusion probability and its Bayes factor is carried by
  the Rao-Blackwellized `n_eff` and `mcse`, while the switching frequency
  mirrors the inclusion probability, so a low `n_eff_mixt` mostly restates that
  the inclusion probability is extreme. The raw directional flip counts stay in
  the table.
* `plot_edge_posterior(binwidth =)` is deprecated and ignored. The panel
  expresses the weight as a density rather than as probability per weight bin,
  so there is no bin to set a width for.

# bgms 0.1.6.3

## New features

* `extract_rhat()`: extract R-hat convergence diagnostics from fitted objects
* `extract_ess()`: extract effective sample size estimates from fitted objects
* `verbose` argument: control informational messages; set `options(bgms.verbose = FALSE)` to suppress globally
* `simulate_mrf()`: standalone MRF simulation with user-specified parameters
* `simulate.bgms()`: generate observations from fitted models (supports parallel processing)
* `predict.bgms()`: compute conditional probabilities P(X_j | X_{-j})
* `main_difference_selection` argument in `bgmCompare()`: control threshold difference selection
* `standardize` argument: scale Cauchy prior by response score range
* `baseline_category` now stored in fitted object for Blume-Capel simulation/prediction

## Bug fixes

* fixed matrix indexing for `posterior_mean_indicator`: now correctly maps C++ row-major order to R matrices (#77)
* fixed mass matrix adaptation: now correctly uses variance instead of precision
* fixed step size heuristic: re-runs after mass matrix updates, resamples momentum each iteration
* fixed E-BFMI diagnostic: now uses actual accepted trajectory momentum
* fixed Blume-Capel interaction: uses centered scores `(c - ref)` in pseudolikelihood denominator

## Other changes

* NUTS: implemented generalized U-turn criterion following Betancourt (2017) and STAN
* NUTS: fused log-posterior and gradient computation eliminates redundant probability evaluations
* bgmCompare: BLAS-vectorized gradient computation for improved performance
* expanded test suite: input validation, extractor functions, S3 methods, simulation, and numerical sanity tests
* improved warmup schedule: fixed buffers (75/25/50) with proportional fallback for short warmup
* edge selection warmup now within user budget: 85%/10%/5% split for stages 1-3a/3b/3c
* streamlined user messages: concise warnings, consolidated NUTS diagnostics
* E-BFMI threshold adjusted to 0.2 (standard)

## Deprecated

* `mrfSampler()` → use `simulate_mrf()`

# bgms 0.1.6.2

## New features

* added option to separately specify beta priors for the within- and between-cluster probabilities for the SBM prior.

## Other changes

* reparameterized the Blume-capel model to use (score-baseline) instead of score.
* implemented a new way to compute the denominators and probabilities. This made their computation both faster and more stable.
* refactored c++ code for better maintainability.
* removed the prepared_data field from bgm objects.

## Bug fixes

* fixed numerical problems with Blume-Capel variables using HMC and NUTS.
* fixed a reporting bug where category thresholds for ordinal variables with a single category were incorrectly expanded to two parameters, resulting in spurious NA values.

# bgms 0.1.6.1

## Other changes

* added extractor function for joint SBM output
* cleaned up documentation, and c++ files
* changed length of warmup phase I in warmup scheduler HMC / NUTS (15% → 7.5%)

## Bug fixes

* fixed a problem with warmup scheduling for adaptive-metropolis in bgmCompare()
* fixed stability problems with parallel sampling for bgm()
* fixed spurious output errors printing to console after user interrupt. 

# bgms 0.1.6.0

## New features

* added NUTS and HMC options for sampling `bgm()` and `bgmCompare()` models
* added support for running multiple chains in parallel
* added user interrupt handling for parallel sampling
* added Markov chain diagnostics (effective sample size and R-hat) for sampled parameters
* added `summary()`, `print()`, and `coef()` methods for fitted objects
* MCMC sampling in `bgm()` and `bgmCompare()` is now reproducible when a `seed` argument is specified

## Other changes

* improved progress bar for parallel sampling
* `summary()` now integrates the functionality of the old `summary_SBM()`
* removed options for modeling main differences; main differences are now always estimated or selected, equivalent to the previous `main_difference_model = "collapse"` setting

## Bug fixes

* fixed an out-of-bounds error in `bgmCompare()` when handling missing data
* fixed a bug in the SBM prior computation

## Deprecated

- In `bgm()`, the following arguments are deprecated:
  - `interaction_scale` → use `pairwise_scale`
  - `burnin` → use `warmup`
  - `save` → no longer needed (all outputs are returned by default)
  - `threshold_alpha`, `threshold_beta` → use `main_alpha`, `main_beta`

- In `bgmCompare()`, arguments related to difference models are deprecated:
  - `main_difference_model` (removed without replacement)
  - `reference_category` → use `baseline_category`
  - `pairwise_difference_*`, `main_difference_*` → use unified `difference_*` arguments
  - `pairwise_beta_bernoulli_*`, `main_beta_bernoulli_*` → use unified `beta_bernoulli_*` arguments
  - `interaction_scale` → use `pairwise_scale`
  - `threshold_alpha`, `threshold_beta` → use `main_alpha`, `main_beta`
  - `burnin` → use `warmup`
  - `save` → no longer needed

- Deprecated extractor functions:
  - `extract_edge_indicators()` → use `extract_indicators()`
  - `extract_pairwise_thresholds()` → use `extract_category_thresholds()`

- Deprecated object fields:
  - `$gamma` (pre-0.1.4) and `$indicator` (0.1.4–0.1.5) → replaced by `$raw_samples$indicator`
  - `$main_effects` (pre-0.1.4) and `$posterior_mean_main` (0.1.4–0.1.5) → replaced by `$raw_samples$main` (raw samples) and `$posterior_summary_main` (summaries)


# bgms 0.1.5.0 (GitHub only)

## New features

* The bgmCompare function now allows for network comparison for two or more groups.
* The new summary_sbm function can be used to summarize the output from the bgm function with the "Stochastic-Block" prior. 
* Two new data sets are included in the package: ADHD and Boredom.

## Other changes

* The bgm function with the "Stochastic-Block" prior can now also return the sampled allocations and block probabilities, and sample and return the number of blocks.
* The underlying R and c++ functions received a massive update to improve their efficiency and maintainance.
* Repository moved to the Bayesian Graphical Modelling Lab organization.
* Included custom c++ implementations for exp and log on Windows. 

## Bug fixes

* Fixed a bug in the bgmCompare function with selecting group differences of blume-capel parameters. Parameter differences that were not selected and should be fixed to zero were still updated.
* Fixed a bug in the bgmCompare function with handling the samples of blume-capel parameters. Output was not properly stored.
* Fixed a bug in the bgmCompare function with handling threshold estimation when missing categories and main_model = "Free". The sufficient statistics and number of categories were not computed correctly.
* Partially fixed a bug in which the bgms package is slower on Windows than on Linux or MacOS. This is because the computation of exp and log using the gcc compiler for Windows is really slow. With a custom c++ implementation, the speed is now closer to the speed achieved on Linux and MacOS.


# bgms 0.1.4.2

## Bug fixes
* fixed a bug with adjusting the variance of the proposal distributions
* fixed a bug with recoding data under the "collapse" condition

## Other changes
* when `selection = TRUE`, the burnin phase now runs `2 * burnin` iterations instead of `1 * burnin`. This ensures the chain starts with well-calibrated parameter values
* changed the maximum standard deviation of the adaptive proposal from 20 back to 2

# bgms 0.1.4.1

This is a minor release that adds some documentation and output bug fixes.

# bgms 0.1.4

## New features
* Comparing the category threshold and pairwise interaction parameters in two independent samples with bgmCompare().
* The Stochastic Block model is a new prior option for the network structure in bgm().

## Other changes
* Exported extractor functions to extract results from bgm objects in a safe way.
* Changed the maximum standard deviation of the adaptive proposal from 2 to 20.
* Some small bug fixes.

# bgms 0.1.3

## New features
* Added support for Bayesian estimation without edge selection to bgm().
* Added support for simulating data from a (mixed) binary, ordinal, and Blume-Capel MRF to mrfSampler()
* Added support for analyzing (mixed) binary, ordinal, and Blume-Capel variables to bgm()

## User level changes
* Removed support of optimization based functions, mple(), mppe(), and bgm.em()
* Removed support for the Unit-Information prior from bgm()
* Removed support to do non-adaptive Metropolis from bgm()
* Reduced file size when saving raw MCMC samples

# bgms 0.1.2

This is a minor release that adds some bug fixes.

# bgms 0.1.1

This is a minor release adding some new features and fixing some minor bugs.

## New features

* Missing data imputation for the bgm function. See the `na.action` option.
* Prior distributions for the network structure in the bgm function. See the `edge_prior` option.
* Adaptive Metropolis as an alternative to the current random walk Metropolis algorithm in the bgm function. See the `adaptive` option.

## User level changes

* Changed the default specification of the interaction prior from UnitInfo to Cauchy. See the `interaction_prior` option.
* Changed the default threshold hyperparameter specification from 1.0 to 0.5. See the `threshold_alpha` and `threshold_beta` options.
* Analysis output now uses the column names of the data.

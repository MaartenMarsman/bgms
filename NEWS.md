# bgms 0.2.0.0

## Changes since the 0.2.0.0 development build

* R-hat is now the classic split-R-hat (Gelman et al. 2013 / Stan; `Rhat = sqrt(var_plus / W)` with `var_plus = (n-1)/n * W + B/n`). The previous estimator applied the Brooks-Gelman degrees-of-freedom adjustment used by `coda::gelman.diag`, which collapsed to a data-independent `sqrt(5/3) ~ 1.29` on nearly-saturated binary edge (or difference) indicators — the typical shape of a decisive edge, where one split sub-chain carries a brief excursion and the rest are constant — so the most decisive edges were flagged as unconverged. Reported R-hat now changes slightly for all parameters (the adjustment factor was ~1 for healthy chains) and substantially for near-saturated indicators (the 1.29 artifact is removed; such indicators report values near 1, or `NA` when every chain is identical). Chains stuck constant at different values now report `+Inf` instead of `NA`. See the diagnostics vignette for how to read indicator R-hat alongside the transition counts and the inclusion probability.
* On Windows with `RcppParallel` >= 6.0.0 (which bundles oneTBB 2022), a multi-core run and a single-core run started from the same `seed` are no longer bit-for-bit identical: the second and later chains follow different sampling trajectories under the new scheduler. Each chain still samples the same posterior — posterior summaries agree to within Monte Carlo error — so results are statistically equivalent, only not bitwise reproducible across core counts. Other platforms and earlier `RcppParallel` are unaffected.
* The default `interaction_prior` is now `normal_prior(scale = 1)` in both `bgm()` (previously `cauchy_prior(scale = 1)`) and `sample_ggm_prior()` (previously `cauchy_prior(scale = 2.5)`), so the GGM prior chain matches `bgm()`. `bgmCompare()` keeps its Cauchy default. The joint-spec normalizer correction table is keyed on the interaction prior, so the first joint fit at the new default cell (with a Beta-Bernoulli or Stochastic-Block edge prior) builds and caches a fresh table once.
* The hierarchical prior trust gauge runs by default on the deployed path, with `options(bgms.zratio_gauge_sweeps = 0L)` as the off switch; `sample_ggm_prior()` keeps its own `zratio_diagnostics` argument. It is silent when clean and prints only flagged chains, and `fit$zratio_diag` carries the per-chain summary whenever it ran. Its cost is fixed per chain rather than proportional to `iter`: two assessment sweeps audit a capped number of edge moves, which is nothing on a sparse posterior (no mediating block is non-trivial there, the normalizer ratio is exact, and there is nothing to audit) and about 5-10 seconds per chain on a dense posterior at 100-200 variables. That is negligible against a production-length fit and noticeable on a short exploratory one, which is what the off switch is for. The advice printed for a flagged chain is now a ladder: it reports the measured error share, the mediating-block sizes the audit covered (`block_lo`, `block_hi`, new in `per_chain`) and how many of the chain's non-trivial moves were audited; then advises raising `options(bgms.zratio_gauge_sweeps)` and refitting to resolve whether the signal is real; and only then mentions `precision_graph_prior = "joint"`, stating that it targets a different model, whose graph marginal is the edge prior reweighted by the per-graph normalizer rather than the edge prior itself. The gauge never switches specification on its own. The diagnostics vignette documents the audit as a sample: it resolves coherent error, which is what `harm_pred` projects and what shifts a recovered network, and does not resolve rare edge-specific failures. The Option-B surface build is also cached (session and disk, keyed on the fit cell and size cap, following `options(bgms.correction_table_cache)`), so repeated fits of the same configuration reuse it instead of rebuilding. A first build uses the fit's own `cores` (the chains have not launched yet; the result is identical for any core count) and can be overridden with `options(bgms.zratio_surface_cores)`; on Windows a socket cluster stands in for forking on large builds. The build now announces itself when it runs, reporting the size cap and the number of cores it uses, and reports the elapsed time when it finishes. This accounts for the pause before the chains start on a first hierarchical fit. A cache hit is silent, and the announcement follows the advisory `bgms.verbose` flag.
* The hierarchical graph-prior per-edge normalizer ratio is now corrected by a theta-independent absolute-moment surface, built once per analysis at the fixed `(eta, delta)` from block-Gibbs anchors, replacing the online warm-up calibrator. The surface predicts the per-component first and second spectral moments for both mediating-block families (common-neighbour clusters and bipartite bridges), and is deployed by decomposing each block into disjoint components. It covers the Normal and Cauchy interaction slabs at the exponential (`alpha = 1`) precision diagonal; a non-unit Gamma diagonal shape falls back to the additive kernel. The surface deploys at every `eta` — the anchors are rebuilt at the analysis's own `(eta, delta)`, unlike the reference implementation's single `eta = 2` build with its `eta < 2` additive tier — with gold gates recorded at `eta = 1` and `eta = 2`. Validated against a per-component gold reference across a graph-density sweep (surface ~0.003 nats to gold vs ~0.036 for additive, up to ~1000x tighter in the dense regime). The `calibration_window` argument (which never left the development build) is removed, and the per-chain `zratio` diagnostics no longer carry the calibration anchors or the `n_clamp`/`n_oracle`/`n_anchors`/`frozen` counters. When a fit deploys the surface on a mediating block larger than its validated size range (dense regions of large graphs, where the moment is extrapolated beyond the trained hull), a single graceful note reports the fraction and the largest block size and points to the trust gauge; sparse graphs never trigger it.
* The hierarchical prior trust gauge gained a second alarm channel: `harm_pred` projects the measured approximation error onto the inclusion-probability scale, using the chain's own edge sensitivities and the edge-prior feedback amplification, and flags when the projected distortion exceeds 0.01. This catches a consistent error that shifts the recovered network without changing individual edge decisions, which the `flip_rate` channel cannot see at chains whose decisions are far from their accept/reject boundaries. Reported per chain in `fit$zratio_diag` (`amplification`, `harm_pred`, `harm_flag`, plus the `se_mcse`/`se_se` uncertainty of the error estimate); computed for Bernoulli and Beta-Bernoulli edge priors. See the diagnostics vignette.

* Past its anchored size range the hierarchical prior's absolute-moment surface now continues along its own boundary slope in log-size instead of freezing at the range's edge. Freezing left the prediction fixed while the truth kept moving, so its error grew without bound in block size. Scored against a block-Gibbs reference on single-component blocks of 90 to 150 variables against an 80-variable anchored range, over two anchor-build seeds: the extension holds a median of 0.0006 nats and at most 0.0011 for common-neighbour blocks, and a median of 0.0043 and at most 0.0060 for bipartite ones, where freezing reaches 0.060 and 0.095 and the additive fallback reaches 0.21 and 0.31. Inside the anchored range the surface is at 0.003 or better, so the extension is bounded and much better than the alternatives rather than validated to the in-range envelope, and the bipartite bound in particular sits about twice outside it. Continuing the fitted quadratic as its own extrapolant was measured and rejected: it degrades with distance (to 0.0026 and 0.016 over the same band) where the boundary-slope tail stays flat. The absolute moments grow with component size, so a fitted slope that points downward at the range edge is a fit pathology; it is floored at zero, which reproduces the old freeze, and counted in the per-chain `n_slope_floor` diagnostic.

* The hierarchical prior's anchored size range now reaches 80 variables, up from 44, with anchor tiers added at 52, 64 and 80 (common-neighbour) and 30, 38, 52, 64 and 80 (bipartite), at 600 block-Gibbs sweeps for sizes from 46 up. The reachable mediating block of an ordinary large fit therefore sits inside the anchored range rather than in the extension zone. The one-time cached build costs 24 s serial and 7 s on four cores at the full range; a fit on fewer variables sizes the build to its own variable count and is unaffected (2 s on four cores at 44). The wider fit was checked against the reference at the small end before adoption: 20 to 42 variables stay at 0.0010 nats or better for common-neighbour blocks and 0.0012 for bipartite ones, inside the 0.003 in-range envelope.

* The hierarchical prior's extrapolation notice now separates warmup from retained sweeps. The sampler initializes from a complete graph, which puts every mediating block at roughly the number of variables for the first sweeps, so the previous whole-run share could report an initial transient as though it described the posterior. The notice now quotes the retained share with the warmup share in brackets, sizes the reported block from the retained sweeps, and says plainly when no retained sweep extrapolated at all. It also quotes the measured accuracy of the extension rather than describing it as a possible slight reduction in accuracy. Trust-gauge sweeps run after sampling and enter neither tally. The per-chain counters gained `n_pred_retained`, `n_extrap_retained` and `max_extrap_size_retained`.

* The hierarchical prior's absolute-moment surface now anchors at the size cap. The anchor grid was filtered to sizes at or below the cap, so a cap that fell between two grid tiers trained the hull at the lower tier: at 40 variables the cap of 40 dropped the size-42 tier and left a hull of 36, and blocks of 37-38 were extrapolated. An anchor tier is now placed at the cap itself whenever the surviving top tier is more than two sizes below it (two, because a mediating block excludes the edge's own endpoints and so never exceeds `cap - 2`). At 40 variables the trained hull goes from 36 to 40 and the extrapolation notice, which reported 0.1% of evaluations at a largest block of 38, no longer fires. The one-time build cost rises by about a second (6.9 s serial, 2.0 s on four cores at the default cap). Surfaces cached by an earlier grid are invalidated by the cache key.

* Under the joint precision-graph specification, a fit with edge selection on a continuous precision block now reports once that the realized edge-inclusion prior is `pi(Gamma) * Z(Gamma)`, the edge prior reweighted by the per-graph normalizer, rather than the nominal edge prior. This holds for every edge prior: at three variables a uniform `beta_bernoulli_prior(1, 1)` realizes about 0.37 and a fixed `bernoulli_prior(0.5)` at `delta = 0` realizes about 0.27. With a learned inclusion probability (`beta_bernoulli_prior()`, `sbm_prior()`) the corrected hyperparameter update keeps the update coherent with the joint model but does not restore the nominal prior; with a fixed inclusion probability nothing absorbs the tilt. The message follows `verbose` (and the advisory `bgms.verbose` flag) and points to `extract_prior_inclusion_probabilities()` for the realized prior and to `precision_graph_prior = "hierarchical"` for the nominal one.

* `extract_inclusion_bf()` gained a `log` argument. With the default `log = FALSE` it returns the inclusion Bayes factor itself, matching its name; `log = TRUE` returns the natural logarithm, which is what the development build returned. The accumulators are exact on the log scale everywhere, while the Bayes factor scale saturates at double precision: an entry whose log exceeds about 709.78 nats (a Bayes factor beyond about 1.8e308) is reported as `+Inf` under `log = FALSE` even though its log-scale value is finite, and an entry whose log is `-Inf` (no inclusion evidence remains) becomes `0`. Workflows that must separate such extreme evidence should use `log = TRUE`. The argument applies to both the `bgm()` and the `bgmCompare()` method, including the stochastic-block difference path, whose return remains posterior odds.

  Scripts written against a development build of 0.2.0.0 before this change silently change meaning on 0.2.0.0: a call that returned the natural log of the inclusion Bayes factor now returns the Bayes factor itself, on the same code and the same call, with no error and no warning. The two scales are hard to tell apart by eye near the decision boundaries, and the difference is invisible in any analysis that only compares values to a threshold. A one-line assertion catches it: with `log = FALSE` an analysis carrying evidence of absence must contain values below 1, whereas the log-scale return is negative there. Set `log = TRUE` to recover the development build's semantics.

* The fit summary's inclusion table reports the Rao-Blackwellized `mcse`, `n_eff`, and `Rhat` on every edge whose Rao-Blackwellized draws vary, where the development build masked them to `NA` for every edge whose indicator never flipped. The masking criterion is now the one that decides whether the quantities exist: the three columns are `NA` when the Rao-Blackwellized draws are constant to double precision, which places the inclusion probability at its numerical bound and the verdict beyond any threshold. Zero flips do not imply constant Rao-Blackwellized draws — a decisive edge's one-step inclusion draw varies over a range of about 1e-3 while the indicator sits still — so the development build blanked a computable Monte Carlo standard error on exactly the near-boundary, short-run edges where it is needed to judge whether a verdict is stable. The directional flip counts (`n0->1`, `n1->0`) are unchanged and remain the record of the indicator's exploration.

* The `n_eff_mixt` column, the transition-based effective sample size of the binary indicator chain, is removed from the fit summary's inclusion table, and `extract_ess(estimator = "mixt")` is deprecated: it warns and recomputes the quantity from the raw indicator draws. Three grounds. It carries no decision value: in a known-truth calibration study its best error-catching threshold rejected 41% of correct edge verdicts to catch 85% of the wrong ones, while a rule built on the inclusion standard error caught all of them at a 2.3% cost. Its own model is violated: the conversion from flip counts to an effective sample size assumes a two-state first-order Markov chain, and the study rejects first-orderness for 41.3% of edge-fits on transitions out of the included state (3.2% out of the excluded state), asymmetrically. And it is redundant: the precision of the inclusion probability and its Bayes factor is carried by the Rao-Blackwellized `n_eff` and `mcse`, while the switching frequency mirrors the inclusion probability, so a low `n_eff_mixt` mostly restates that the inclusion probability is extreme. The raw directional flip counts stay in the table.

* `plot()` on a `bgmCompare()` fit draws the group differences the data settle. The default `type = "difference"` is one network of the difference indicators, encoded exactly as `plot()` on a `bgm()` fit encodes edges: solid with a sign-carrying colour and a width scaled by the size of the difference where the data settle one, thin dotted grey where they leave it undecided, omitted where they rule it out. So the default picture and `verdicts(fit)` state the same thing at the same threshold, which matters here because the difference indicators are what `bgmCompare()` parameterizes. A network with no edges left is a result rather than an error — "the groups do not differ anywhere" is a common and correct finding — so the nodes are drawn on their own and the subtitle says so. Main-effect differences are not edges, and when `main_difference_selection = TRUE` gave them their own indicators their verdicts ride on the nodes: a square node with an accented border where the data settle a main-effect difference, a circle where they do not, so the encoding does not rest on colour alone. Under the default `main_difference_selection = FALSE` those indicators are never updated and have no verdict, every node is drawn alike, and the subtitle says the channel is empty; `verdicts()` remains where main-effect differences are read precisely. `type = "groups"` draws each group's own network beside the difference panel on one shared layout, so a node sits in the same place throughout and a reader compares by position; those panels use the same colour-vision-safe sign pair as the rest of the package rather than qgraph's green/red default. `type = "centrality"` draws the centrality display, with `group` passed through. Layout comes from qgraph, which stays a suggested package.

* `extract_centrality()` works on `bgmCompare()` fits. `group = 1` gives that group's centrality, one row per posterior draw; `group = c(1, 2)` gives the difference between two groups' centralities, again per draw, so the credible interval is the interval of the difference and answers whether the groups differ in a node's centrality directly, where two separately drawn intervals do not. Each group's network is rebuilt on every draw as `baseline + (P %*% differences)` with the fit's own contrast projection, not from posterior means, which is what carries the uncertainty through; averaged over draws it reproduces `extract_group_params()` to machine precision. A draw in which a difference indicator is excluded gives both groups the same edge weight and contributes exactly zero, so a difference carries a point mass at zero and `summary()` reports `p_positive` (the probability the first group's centrality is higher) rather than `p_most_central`, with the remaining mass on no difference at all. `plot()` marks zero.

* `prior_sensitivity_check()` works on `bgmCompare()` fits, tracing each difference indicator's inclusion Bayes factor across the difference slab scale. `bgmCompare()` gives the pairwise and the main-effect difference families one `difference_scale`, so a single curve covers both. The anchored-curve machinery, the convergence gate, the run-to-run noise band and the reporting are the same as for `bgm()`; what differs is that one difference indicator gates several parameters (a pairwise difference in every group contrast, a main-effect difference across that variable's whole threshold block), so the importance weight sums the difference-slab density ratio over the gated parameters while the inclusion probability is still reported per indicator. Reweighting a fit at one scale to another was checked against real refits at that scale: on Wenchuan at a doubling of the scale the reweighted inclusion probabilities land within 0.016 of a refit's, against a 0.007 run-to-run spread between two refits at that same scale, at an importance effective sample size of 2536 out of 12000 draws. Indicators the sampler never updated (the main-effect differences, unless `main_difference_selection = TRUE`) carry no verdict at any scale, and the printed report counts them out rather than reporting them as undecided. A stochastic-block difference prior has no single marginal inclusion probability and is refused, since the curve would report posterior odds rather than Bayes factors.

* `prior_sensitivity_check()` gained a `vary` argument, and its printed report now names which prior the curve moved. On a continuous or mixed fit the interaction slab scale and the prior on the precision diagonal are tied through the standardized frame (raw rate `= eta / s`), so a sweep of the slab has to hold one of the two fixed, and the two choices answer different questions: `vary = "slab"` holds the raw diagonal rate at the fitted value and moves the interaction prior alone, while `vary = "slab-and-diagonal"` holds `eta` fixed and lets the raw rate follow, rescaling the prior's overall size with its shape held fixed. The default `vary = "auto"` follows the frame the fit itself used — the joint sweep when the diagonal prior was given as `eta`, the slab alone when it was given as a raw `rate`. Discrete fits have no precision diagonal; `vary` is accepted and has no effect there, and the report omits the line. This changes the curve reported for a continuous or mixed fit whose diagonal prior was specified in the standardized frame, which is the default: such a fit previously swept the slab alone, because the spec stores the diagonal rate already resolved to the raw frame, and now sweeps both. On a five-variable Gaussian graphical model the two modes give visibly different curves, so the choice is not a formality. Pass `vary = "slab"` to recover the previous behaviour.

* `prior_sensitivity_check()`'s run-to-run noise band no longer becomes infinite when an edge saturates in one of the two identical refits that define it. The band is the 95th percentile of the per-edge spread between those refits over threshold-relevant edges; an edge that is threshold-relevant in one and saturated (`log10 BF` of `-Inf` or `+Inf`) in the other gave a spread of `Inf`, which `na.rm = TRUE` did not remove, so the band itself became `Inf` and every verdict change along the curve was reported as `indistinguishable-from-wobble` — the mover rule switched off silently, with no warning and a normal-looking report. Such an edge has a censored rather than an infinite spread and is now left out of the percentile; the printed report counts the edges left out, and when no threshold-relevant edge has a measurable spread the band is `NA` and the mover rule falls back to `max(tolerance, 2 * MCSE)`, which the report also states. Any edge can trigger this, not only a Gaussian one; runs whose reported noise band was finite are unaffected.

## Breaking changes

* `update_method = "hamiltonian-mc"` has been removed. Use `update_method = "nuts"` instead. NUTS dynamically adapts trajectory length and is more reliable, especially with edge selection on GGM models.
* The `hmc_num_leapfrogs` argument has been removed along with pure HMC.
* The `standardize` argument of `bgm()` and `bgmCompare()` has been removed: pairwise interactions are on the association scale and share one prior scale, so the per-pair adjustment by the product of the variables' maximum scores no longer applies. The argument stays a deprecated formal so the failure is informative rather than an unused-argument error: `standardize = FALSE`, the old default, warns and proceeds because it is what the sampler already does, while `standardize = TRUE` errors and points to setting the scale directly through `interaction_prior` (and `difference_scale` for `bgmCompare()`).

* Pairwise interaction parameters for ordinal MRFs are now stored on association scale (half the sigma scale used in 0.1.6.3). Code that interprets raw pairwise posterior samples or sets `pairwise_scale` explicitly will need adjustment.
* Default `pairwise_scale` changed from 2.5 to 1 to match the association-scale reparameterization.

* `bgmCompare()` stores pairwise interaction parameters on the association scale, matching `bgm()`. Its sampler put `omega * x` in the pseudolikelihood's linear predictor where every other path in the package puts `2 * omega * x` — `bgm()`'s sampler, `simulate_mrf()`, both `predict()` methods, and the mixed-model cross terms — so a `bgmCompare()` fit's pairwise effects came out about twice a `bgm()` fit's for the same coupling, and the two could not be read on one scale. Pairwise effects reported by `coef()`, `summary()`, `extract_pairwise_interactions()`, `extract_group_params()`, and the raw posterior samples are now about half their previous values. Main effects, indicators, and inclusion probabilities are unaffected in scale.

  `predict()` and `simulate()` on a `bgmCompare()` fit were wrong as a consequence and are now correct. Both hand the fit's pairwise effects to code that applies the factor two itself, so both had been applying each interaction twice over: on Wenchuan, `predict()`'s category probabilities for a five-variable fit came out `0.037 0.164 0.134 0.216 0.448` against observed marginals of `0.079 0.383 0.223 0.208 0.107`, where the same data through `bgm()` gave `0.075 0.384 0.224 0.210 0.108`. They now reproduce the marginals as the `bgm()` methods do.

  `interaction_prior` and `difference_scale` keep their numeric defaults, and a scale of `s` now means for `bgmCompare()` what it already means for `bgm()`: the prior acts on the association coordinate in both. Because that scale previously acted on parameters twice as large, the same number is a tighter prior than before, and difference inclusion Bayes factors move accordingly. Analyses that set `difference_scale` or `pairwise_scale` explicitly, or that read raw pairwise samples, need re-running rather than rescaling. The calibration of the `difference_scale` default under the association-scale parameterization is under study, so difference verdicts near a decision threshold should be read as scale-contingent.
* `extract_category_thresholds()` is deprecated in favor of `extract_main_effects()`, which covers category thresholds, continuous means, and precision diagonal entries.
* `extract_ess()` reports the Rao-Blackwellized effective sample size for the edge (or difference) indicators, the `n_eff` column of the fit summary's inclusion table, where 0.1.6.3 returned the indicator chain's transition-based ESS. The extractor family is now coherent: the Rao-Blackwellized inclusion probability (`extract_posterior_inclusion_probabilities()`), its R-hat (`extract_rhat()`), and its ESS all describe the same estimate. The transition ESS is deprecated (see below); `extract_ess(fit, estimator = "mixt")` warns and recomputes it from the raw indicator draws.

## New features

* `verdicts()` reads each edge (or difference) indicator's inclusion Bayes factor as a three-way verdict — evidence of presence, undecided, evidence of absence, at a threshold and its reciprocal — and flags the verdicts that a rerun of the sampler could change. The flag is the union of two standard errors of the logit inclusion probability: the Jeffreys-smoothed two-state model of the binary indicator chain, which stays defined when the indicator never flips, and the delta-method transform of the Rao-Blackwellized Monte Carlo standard error. An edge is fragile when either places a verdict boundary within two standard errors of the estimated evidence. The rule comes from a known-truth calibration study of 37,010 graded edge-fits across ordinal, binary, and Gaussian graphical models, in which every one of the 66 verdict errors sat within 0.25 of a threshold on the log10 Bayes factor scale and no edge further out was ever misclassified: over those fits the two-state standard error alone recalled 0.742 of the errors and the Rao-Blackwellized one 0.939, while the union recalled all of them at a 3.0% false-alarm rate and transferred across model types. Classification goes through the same internal routine as `prior_sensitivity_check()`, so the two classify identically, and distances are measured on the log10 Bayes factor scale, which has the prior inclusion odds already divided out and stays finite where the Bayes factor saturates. Printing reports the verdict tally and, when any verdict is fragile, says so. Every arm of that study fitted a single network, so the operating point belongs to the edge indicators of `bgm()`; on `bgmCompare()` difference indicators the flag still marks verdicts near a boundary, but its error-catching rate and false-alarm rate there are unmeasured, and the print method says so rather than borrowing the single-network numbers.

* `plot()` on a `bgm()` fit draws the model-averaged network with edges encoded by what the data settle about them: solid with weight-scaled width and a sign-carrying colour for evidence of presence, thin dotted grey for undecided, and omitted for evidence of absence. The threshold is the one `verdicts()` uses, so the default picture and the reporting table state the same thing. `type = "centrality"` draws the centrality display instead. Layout comes from qgraph, which stays a suggested package; without it the method errors and points to `verdicts()`. The three-panel evidence display, the structure plots, and the other network displays remain easybgm's and are deliberately not duplicated.

* `plot_edge_posterior()` draws one edge's spike-and-slab posterior: a stem at zero whose height is the posterior probability that the edge is absent, with the printed probability, and the slab of the weight given presence, drawn as probability per weight bin so that the two parts divide the posterior mass as the model does.

* `extract_centrality()` evaluates a node centrality on every posterior draw of the network and returns a draws-by-nodes matrix, so the posterior distribution of a node's centrality comes out of the same model-averaged posterior as the edge weights, with structural uncertainty included: an edge excluded at a given iteration contributes exactly zero there. `summary()` reports the posterior mean, credible interval, and the posterior probability of being the most central node; `plot()` draws the ordered interval display. `measure = "strength"` is the one accepted value for now.


* `calibration_check()` fits the isotonic (pool-adjacent-violators) reliability curve of the model's conditional predictions, one panel per variable: pooling across variables lets a variable predicted too high and one predicted too low cancel, and the pooled curve then tracks the diagonal while neither variable does. The consistency band resamples each case's category from its own predicted distribution rather than resampling the cumulative threshold events independently, because those events are nested within a case and treating them as independent understates the band. Evaluated on the fitted data by default; `newdata` makes it out-of-sample.

  Continuous variables are covered too, through the probability integral transform: `u = F(y | rest)` is uniform exactly when the conditional predictive distribution is right, so the panel is the empirical distribution function of the `u` values against the uniform diagonal. `F` is the predictive mixture over `ndraws` posterior draws rather than a plug-in Gaussian at the posterior mean, so the parameter uncertainty sits inside the distribution the observation is transformed by; the two panel kinds therefore differ in what they condition on, and the `kind` column of both returned tables records which construction produced each row. The continuous band is not resampled from the model, because under the transform the null is `Uniform(0, 1)` whatever the conditional density was: it is the simultaneous envelope of the empirical distribution functions of `n` independent uniforms, computed once for every continuous variable in the fit. Both panels live on the unit square with the diagonal as the calibrated reference, so a mixed fit produces one figure and one summary table. In-sample the observation also entered the parameters it is judged against, which makes the transform mildly under-dispersed and the check conservative; `newdata` removes that. The `curves` table's `pav` column is renamed `curve`, since it now holds either construction's curve.


* New vignette "Checking your fitted model": `verdicts()` and the verdict-encoded network, then `calibration_check()`, then a worked example of building a predictive display directly on `simulate()` — the joint-level check is an analyst's reading of a display, not a pass/fail statistic, so the vignette shows how to roll one rather than calling a blessed function. The worked statistic is the sum-score distribution against its pointwise predictive band, chosen because it lies outside the model's sufficient set: pairwise dependence statistics sit close to that set, so how they read depends on how the estimator is anchored rather than on whether the model fits.

* `prior_sensitivity_check()` recovers the continuous inclusion-Bayes-factor curve of each edge across the interaction slab scale, and classifies how every edge's verdict behaves along it. The curve is anchored at a handful of fixed-scale fits (log-spaced `anchors`, default `0.4 / 0.63 / 1 / 1.6 / 2.5` times the chosen scale) and filled in between anchors by importance reweighting: a fixed-scale fit is reweighted to a nearby scale with per-draw slab-density ratios over the included edges, the likelihood cancels, and no normalizing constant is involved. The `1x` anchor is the original fit itself and is never refit, so the chosen-scale verdicts the check reports are exactly the analysis already run. Each display point pools every anchor that clears `ess_floor` (default 400) there, weighting each anchor's inclusion-probability estimate by its inverse variance and transforming the pooled probability to the Bayes-factor scale; points where no anchor clears the floor are reported `NA` rather than extrapolated, and non-overlapping anchor radii warn to add anchors. An edge is called scale-sensitive only if its verdict changes along the curve *and* its Bayes-factor swing exceeds a run-to-run noise band (the 95th percentile of the spread between one anchor refit and a repeat of it), so boundary edges that merely flip between reruns of the same prior are reported as `indistinguishable-from-wobble` rather than as moves. Each anchor passes a convergence gate (median continuous split-R-hat, bulk Rao-Blackwellized inclusion R-hat, and E-BFMI for NUTS) before it enters the curve. The check works on any `bgm()` fit with edge selection; on Wenchuan it costs about one original fit with warm-started NUTS refits (a `refit_sampler` argument, default `"same-as-fit"`, carries the adapted NUTS step size and mass matrix). The single-fit conditional-density reweighting curve of the earlier development build has been removed: reweighting one fit across a global scale collapses (its importance effective sample size is about 1 across the advertised band; see `dev/audit/2026-07-28-phase-b-review.md`).
* `sample_graph_prior()`: draws edge-inclusion indicators, together with any edge-prior hyperparameters, from the graph level of the spike-and-slab prior. Under the hierarchical specification the draw is ancestral and exact; under the joint specification it runs the zero-data prior chain, so the draws carry the per-graph normalizer tilt. Hyperparameters can be fixed instead of sampled (`theta` for Bernoulli and Beta-Bernoulli priors, `allocations` plus `block_probs` for the Stochastic-Block prior).
* `sample_sbm_prior()`: ancestral draws of block allocations and pair-inclusion probabilities from the MFM-SBM edge-prior hyperprior.
* The edge-prior correction table now announces itself only when it actually builds (a cache hit is silent) and states that the build is one-time and cached. In interactive sessions it draws a progress bar for both serial and parallel builds (a parallel build advances the bar as each batch of forked chains completes). The bar follows `display_progress` (the same control as the sampler's bar), not the advisory `bgms.verbose` flag.
* Gaussian graphical models (GGM): `bgm(x, variable_type = "continuous")` fits a GGM with Bayesian edge selection. Sampling uses NUTS on a free-element Cholesky (theta-space) parameterization of the precision matrix, which keeps the precision matrix positive-definite by construction; adaptive-metropolis is also available.
* GGM Gibbs sampler: `update_method = "gibbs"` fits a GGM with a conjugate row-by-row update of the precision matrix, with or without edge selection. It needs no step-size or proposal tuning and supports a Normal or Cauchy prior on the edges. Continuous data only.
* Gibbs warmup staging: with `update_method = "gibbs"` and edge selection, the first 15% of the warmup runs the full model so the precision matrix settles, and edge selection is active for the remaining 85%. Previously edge selection only started at the first retained iteration, so the graph's equilibration happened inside the retained samples. Both windows scale with the warmup budget; the warmup default is unchanged.
* Standardized-frame diagonal prior: `gamma_prior()` and `exponential_prior()` accept the rate as `eta`, the rate on the precision diagonal in the frame where the pairwise (slab) prior has unit scale. The raw rate is derived at fit time as `eta / s` for interaction-prior scale `s`, so the prior geometry is held fixed when the slab scale changes; at fixed `eta`, graph and partial-correlation inference is invariant to the slab scale. The default `precision_scale_prior` is now `exponential_prior(eta = 1)` (the same distribution as `gamma_prior(shape = 1, eta = 1)`), which is identical to the previous `gamma_prior(shape = 1, rate = 1)` at the default unit slab scale; with a non-unit slab scale the diagonal rate now co-scales. Applies to GGM and mixed MRF models. One consequence of the new default: a scale-free interaction prior (`beta_prime_prior()`) on a continuous model now requires specifying the diagonal prior with `rate` instead of `eta`.
* Mixed MRF models: `bgm()` accepts a per-variable `variable_type` vector that mixes `"ordinal"`, `"blume-capel"`, and `"continuous"` types to estimate networks with both discrete and continuous variables. `simulate.bgms()` and `predict.bgms()` also support mixed models.
* Missing data imputation: `na_action = "impute"` integrates over missing values during MCMC sampling for ordinal, continuous, and mixed models.
* `extract_precision()`: extract posterior precision matrix samples from GGM and mixed models.
* `extract_partial_correlations()`: extract posterior partial correlation samples from GGM and mixed models.
* `extract_log_odds()`: extract log-odds for discrete pairwise interactions.
* `extract_main_effects()`: extract main effect samples (category thresholds, continuous means, and precision diagonal).
* `extract_posterior_inclusion_probabilities(estimator = "rb")` and `extract_inclusion_bf()`: Rao-Blackwellized inclusion inference. For each indicator update the sampler records the one-step draw `J = gamma + (1 - 2 gamma) alpha` (pre-move state `gamma`, birth/death acceptance probability `alpha`), stored per chain and per edge in `fit$raw_samples$rb_inclusion` alongside the raw indicator draws. Averaging `J` (the `"rb"` estimator, now the default) is a lower-variance estimator of the inclusion probability than the raw indicator average (`estimator = "raw"`). It also accumulates the inclusion odds on the acceptance-probability scale, so `extract_inclusion_bf()` reports a log inclusion Bayes factor that stays finite down to log-acceptances of about -745; edges whose raw average saturates at 0 or 1 because the chain never flipped in a finite run still receive a finite, data-driven Bayes factor. The prior inclusion odds are divided out (from `extract_prior_inclusion_probabilities()` for `bgm()`, from the exchangeable difference prior for `bgmCompare()`), so the reported value is a Bayes factor rather than posterior odds and equals the posterior odds only at a prior inclusion probability of 1/2. (The averaged `J` itself still rounds to a machine 0 or 1 near the boundary, since `1 - alpha` underflows in double precision; the Bayes factor extractor avoids forming it.) The RB draw is continuous, so the standard MCSE/ESS/split-R-hat machinery applies to inclusion inference: the printed fit summary's inclusion table now reports the RB estimate with a continuous MCSE (`mcse`), ESS (`n_eff`) and split-R-hat (`Rhat`) beside the per-direction transition counts (`n0->1`, `n1->0`), in place of the binary indicator average alone. The three RB columns are `NA` where the `J` draws are constant to double precision, and computed everywhere else. The Rao-Blackwellized estimate is now the canonical inclusion probability throughout a fit: `extract_posterior_inclusion_probabilities()` defaults to it (`estimator = "rb"`, with `"raw"` still available), and `fit$posterior_mean_indicator` reports it. Available for `bgm()` (ordinal, continuous/GGM, and mixed) and `bgmCompare()` (difference indicators); the estimator inherits the chain's mixing and does not rescue a chain that failed to explore the model space.
* NUTS diagnostics now include the per-iteration mean Metropolis acceptance probability (`fit$nuts_diag$accept_prob`, paralleling Stan's `accept_stat__`) and a per-chain `mean_accept_prob` summary.
* `sample_ggm_prior()` accepts `update_method = "gibbs"` for the `spec = "joint"` prior chain, using the conjugate row and edge updates instead of adaptive Metropolis, and `edge_prior = "beta-bernoulli"` for sampling the inclusion probability under its Beta hyperprior.
* Corrected inclusion-probability updates for `beta_bernoulli_prior()` on continuous (GGM) models: under the determinant-tilted precision prior, the conjugate Beta update omits a normalizing-constant factor and biases the sampled inclusion probability toward sparsity (at 5 variables with a uniform hyperprior its prior mean lands near 0.37 instead of 0.5). The update now draws from the corrected conditional using a table built from the prior distribution at the first fit of a model configuration and cached on disk (`tools::R_user_dir("bgms", "cache")`; a one-time cost of the order of minutes, announced when `verbose = TRUE`). The sampled inclusion probability is returned per chain in `fit$inclusion_parameter_samples`.
* Corrected stochastic block model updates for `sbm_prior()` on continuous (GGM) models, using the same cached table: the block-probability draws, the block-allocation weights, and the new-cluster weight all carry the normalizing-constant correction. Without it the sampled partition collapses toward one block (at 20 variables the prior mean number of blocks lands near 1.0 instead of 1.87); with it a prior-only chain reproduces the model's partition prior (total variation 0.01-0.02 at 5-20 variables). `sample_ggm_prior()` accepts `edge_prior` objects (`bernoulli_prior()`, `beta_bernoulli_prior()`, `sbm_prior()`) and returns sampled allocations under the block model.
* The normalizing-constant correction extends to mixed models with `beta_bernoulli_prior()` or `sbm_prior()`: the determinant tilt acts on the continuous precision block, so the correction table is built for the continuous variables and the block-structure corrections read continuous-continuous edges only. With fewer than two continuous variables no edge is tilted and the plain conjugate updates apply unchanged; with exactly two, the single tilted pair supports the inclusion-probability correction but not the block-model slope curve, so `sbm_prior()` warns and keeps the plain conjugate updates there. Prior-only mixed chains reproduce the Beta hyperprior on the inclusion probability and the partition prior on the number of blocks; without the correction the inclusion probability biases toward sparsity and the partition toward fewer blocks.
* Hierarchical graph-prior specification for continuous (GGM) models: `bgm(precision_graph_prior = "hierarchical")` composes the edge prior and the precision prior as `p(Gamma) p(K | Gamma)` with `p(K | Gamma)` normalized per graph, so the graph marginal is exactly the edge prior (under the joint specification it is reweighted by the per-graph normalizer). Each edge move evaluates the normalizer ratio with a fast per-component surface approximation, built once per analysis at the fixed `(eta, delta)` from block-Gibbs anchors. Requires `edge_selection = TRUE` and a `normal_prior()` or `cauchy_prior()` interaction prior; the precision scale prior may be `exponential_prior()` or a `gamma_prior()` of any shape (the Z-ratio constants carry the diagonal Gamma shape through every channel and the block-Gibbs anchor oracle). The same specification is available in `sample_ggm_prior(spec = "hierarchical")`. In simulation checks with data (50 variables, n = 25, Beta-Bernoulli(2, 4)) the posterior inclusion probabilities are unbiased to within +0.002. Chains sampled from the prior alone, with a sampled inclusion probability and many variables, can fall outside the approximation's validated range; the trust gauge flags affected chains. On mixed models the specification normalizes the continuous block `p(K_yy | Gamma_yy)`: the approximation enters the continuous-continuous edge moves only, and at least two continuous variables are required.
* Trust gauge for the hierarchical specification: while the sampler runs, it redoes a subset of each chain's edge add/remove decisions with the exact calculation and records `flip_rate`, the fraction of decisions that would have come out differently. A chain whose flip rate exceeds 1% is flagged, and the warning prints like other sampler warnings. The per-chain summary is attached as `fit$zratio_diag` (and as `$zratio_diagnostics` on `sample_ggm_prior()` output); see `summarize_zratio_gauge()` and the diagnostics vignette.
* The approximation's constants are built at unit slab scale (diagonal rate `eta = pairwise_scale * scale_rate`), where the normalizer ratio is invariant to the slab scale, so results no longer depend on the user's scale choice. Previously the constants were built at the raw slab scale: a `normal_prior(scale = 2.5)` fit gave a biased approximation, with the `sample_ggm_prior()` hierarchical graph marginal at 6 variables sitting at 0.292 for a 0.30 edge prior; unit-scale fits were unaffected.
* The hierarchical specification supports `gamma_prior()` diagonals of any shape: the Z-ratio constants generalize the diagonal prior from the exponential to `Gamma(shape, rate)` (generalized Gauss-Laguerre pair and bridge integrals, a shape-shifted node channel, and an independence-Metropolis correction in the clique-2 and block-Gibbs anchor sweeps), and the constants cell extends to `(delta, eta, shape, slab)`. The graph-law identity holds at shapes 0.5 and 2 on both update methods; at shape 1 every code path is unchanged.
* A `cauchy_prior()` interaction prior under the hierarchical specification now gets its own constants instead of reusing the Normal-prior ones. The reuse broke the hierarchical graph law: at 6 variables the graph marginal read 0.247 for a 0.30 edge prior. It now sits within Monte Carlo error of the target on both the adaptive-metropolis and gibbs update methods.
* The NUTS diagnostics summary now prints the `warmup_incomplete` flag (energy not stationary) it already computed.
* `extract_prior_inclusion_probabilities()`: prior edge-inclusion probabilities in the same matrix layout as `extract_posterior_inclusion_probabilities()`, for prior/posterior inclusion-odds computations. Under the joint spike-and-slab prior on a continuous block the graph marginal is reweighted by the per-graph normalizer (positive-definite-cone mass shaped by the determinant tilt), so continuous-continuous edges do not keep the edge-prior marginal at any `delta` — with a uniform Beta-Bernoulli hyperprior at 3 variables the prior edge probability is about 0.37 rather than 0.5, and a fixed `bernoulli_prior(0.5)` at `delta = 0` lands near 0.27. The values are read from the cached correction table (`bernoulli_prior()`, `beta_bernoulli_prior()`) or estimated by a prior-only chain with the fit's own correction settings (`sbm_prior()`, cached on the fit). Mixed models report per-edge-class values (discrete-discrete, continuous-continuous, cross); ordinal models use the analytic edge-prior marginals, including the exchangeable partition mixture for `sbm_prior()`.

## Other changes

* Fitted objects from `bgm()` and `bgmCompare()` are now S7 class objects (new dependency: `S7`). All existing `$`, `[[`, and `names()` access patterns continue to work. When an incompatible `easybgm` version is loaded, bgms returns plain S3 lists for backwards compatibility; this shim will be removed in a future release.
* Refactored the C++ backend: unified model hierarchy (`BaseModel` → `GGMModel` / `OMRFModel` / `MixedMRFModel`), shared NUTS/HMC infrastructure, and fused log-posterior and gradient computation.
* NUTS now uses Stan's multinomial candidate weighting (log-sum-exp of `H0 - h` per leaf, biased progressive sampling at the top level) in place of the Hoffman-Gelman slice variable. The two schemes target the same posterior; the multinomial variant produces lower-variance candidate selection and has been Stan's default since 2017. User-facing output is unchanged apart from the new `accept_prob` diagnostic.
* NUTS Stage-2 warmup windowing now matches Stan's `windowed_adaptation::compute_next_window`: when the window after the next would overshoot the Stage-3a boundary, the current next window is stretched to absorb the remaining Stage-2 budget instead of emitting a small trailing window. This eliminates a disruptive mass-matrix update + step-size reinit at the end of warmup and improves dual-averaging convergence.
* Dropped `coda` from Imports; ESS and R-hat are now computed in C++ with on-demand (lazy) evaluation, replacing the eager R-based computation from 0.1.6.3.
* `$` and `[[` accessors on fitted objects trigger lazy computation of MCMC diagnostics on first access.

## Bug fixes

* Fixed a category-scale mismatch between `simulate()` and `predict()` for discrete variables whose observed values were not already coded `0, 1, 2, ...`. `bgm()` recodes each ordinal variable's values to internal 0-based categories; `simulate()` returned data on the internal scale while the ordinal (OMRF) and group-comparison `predict()` paths expect new data on the original scale, so `simulate()` followed by `predict()` mismatched, giving a "category values not observed in the training data" warning and a miscoded prediction context. Separately, the mixed-model `predict()` path did not recode discrete new data at all, so passing original-scale data to it was silently miscoded with no warning (and could index category parameters out of range). Now `simulate()` returns discrete data on the original category scale for all model types, and the mixed `predict()` path recodes discrete new data the same way the ordinal path does; category relabeling leaves predictions unchanged. Blume-Capel scores carry no recode map and are unchanged.
* Fixed `extract_posterior_inclusion_probabilities()` on mixed MRF fits: the indicator samples are stored in block order (discrete-discrete, continuous-continuous, cross) over internally reordered variables, but the extractor filled the matrix in global pair order, misplacing the entries. It now maps them through the same block filler as `fit$posterior_mean_indicator`.
* Fixed the number-of-blocks summary for `sbm_prior()` fits: the conditional p(K | t) behind `posterior_num_blocks` placed a zero-truncated Poisson prior on the number of components, while the sampler's partition coefficients use the shifted Poisson (K - 1 ~ Poisson(lambda)). The mismatch reweighted the reported distribution by a factor lambda/K across K. The summary now uses the sampler's convention and is unit-tested against the generative partition prior.
* Fixed an indicator bookkeeping asymmetry in mixed MRF models: cross (discrete-continuous) edge moves updated only the upper triangle of the internal edge-indicator matrix, while the stochastic block edge prior reads full columns. Under `edge_prior = sbm_prior()`, the block-allocation update for a discrete variable therefore saw its cross edges as always included. Discrete-discrete and continuous-continuous moves, other edge priors, and the reported indicator samples were not affected.
* Fixed the `delta = NULL` default for mixed MRF models: the determinant-tilt exponent applies to the continuous precision block, but the default counted blume-capel variables in its dimension. A mixed model with blume-capel variables now gets `0.5 * log(#continuous)` instead of `0.5 * log(#continuous + #blume-capel)`; models whose discrete variables are all ordinal are unchanged.
* A rank-1 Cholesky downdate that would leave the precision factor non-positive-definite is now detected and triggers a rebuild of the factor from the precision matrix. Previously the failure signal was never checked and a partially updated factor could silently corrupt the carried covariance in long GGM or mixed MRF edge-selection runs.
* Fixed compilation failure on Alpine/musl: `mrf_simulation.cpp` relied on a transitive include for `<tbb/global_control.h>` that is not available on all platforms.
* Fixed stale gradient cache after missing data imputation caused NUTS to use outdated cached values for leapfrog integration.
* Fixed stale observation transpose after missing data imputation caused the pairwise gradient to use stale data.
* Fixed NUTS acceptance probability: target_accept now correctly passed to lower-level NUTS functions.
* Fixed NUTS acceptance probability accumulation: the top-level trajectory loop overwrote the Metropolis contribution with the last subtree's value instead of summing across the full trajectory, biasing the signal used by dual-averaging step-size adaptation.

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
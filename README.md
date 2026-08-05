
<!-- README.md is generated from Readme.Rmd. Please edit that file -->

<p align="center">

<img src="man/figures/bgms-banner.svg" width="100%" alt="bgms">
</p>

<!-- badges: start -->

[![CRAN
Version](https://www.r-pkg.org/badges/version/bgms)](https://cran.r-project.org/package=bgms)
[![Downloads](https://cranlogs.r-pkg.org/badges/bgms)](https://cran.r-project.org/package=bgms)
[![Total](https://cranlogs.r-pkg.org/badges/grand-total/bgms)](https://cran.r-project.org/package=bgms)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![Codecov](https://codecov.io/gh/Bayesian-Graphical-Modelling-Lab/bgms/branch/main/graph/badge.svg)](https://app.codecov.io/gh/Bayesian-Graphical-Modelling-Lab/bgms)
<!-- badges: end -->

**Bayesian analysis of graphical models**

The **bgms** package provides Bayesian estimation and edge selection for
Markov random field models of mixed binary, ordinal, and continuous
variables. The variable types in the data determine the model: an
**ordinal MRF** for ordinal data, a **Gaussian graphical model** for
continuous data, or a **mixed MRF** combining both. Posterior inference
uses Markov chain Monte Carlo, combining a Metropolis approach for
between-model moves (i.e., edge selection) with the No-U-Turn sampler
for within-model parameter updates. The package supports both
single-threaded and parallel chains, and uses a C++ backend for
computational efficiency.

## Main functions

- `bgm()` — estimate a graphical model in a one-sample design.
- `bgmCompare()` — compare graphical models between groups.

Both functions support **edge selection** via spike-and-slab priors,
yielding posterior inclusion probabilities for each edge. `bgm()` can
additionally model **community structure**, and `bgmCompare()` can test
for **group differences** in individual parameters.

## Large Gaussian graphical models

Continuous data are where graphs get large, and large graphs are where
runtime becomes the constraint. `bgm()` has a faster route for that
case: a conjugate Gibbs sampler for the Gaussian graphical model,
paired with a precision-graph prior that is normalized once rather than
once per graph. Two arguments select it.

``` r
fit = bgm(y,
  variable_type = "continuous",
  precision_graph_prior = "joint",
  update_method = "gibbs"
)
```

On an Apple M5 Pro, four chains on four threads, at the package's
default chain length: a 50-variable graph in 3.1 seconds against 31.4
for the default route, and a 200-variable graph, 19,900 candidate
edges, in under 7 minutes.

The joint specification is a different model, not only a faster sampler:
its prior over graphs is the edge prior reweighted, so inclusion Bayes
factors are not interchangeable between the two routes. The vignette
`vignette("fast-ggm")`, also on the [package
website](https://bayesian-graphical-modelling-lab.github.io/bgms/articles/),
sets out that trade-off, the timings, and where the route applies.

Measured against the neighboring CRAN packages on the same data and
machine, the route's headline is time to a converged edge ranking:
seconds for `bgms` at every tested sample size, a median 18 to 25
times sooner than **ssgraph** at a matched prior — replicated, and
partly because a `bgms` sweep costs 6.6 times less at 200 variables.
**modelSelection** converges about as fast, by a design that trades
away some recall and exact reproducibility. The full comparison, with
every setting and seeded script, is on the
[comparison page](https://bayesian-graphical-modelling-lab.github.io/bgms/guide/fast-ggm.html).

## Installation

Install from CRAN:

``` r
install.packages("bgms")
```

Or install the development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("Bayesian-Graphical-Modelling-Lab/bgms@develop")
```

## Citation

If you use bgms in your research, please cite the software package:

``` r
citation("bgms")
toBibtex(citation("bgms"))
```

Additional citation formats are available on the
[package website](https://bayesian-graphical-modelling-lab.github.io/bgms/).

## Contributing

Contributions are welcome. See
[CONTRIBUTING.md](https://github.com/Bayesian-Graphical-Modelling-Lab/bgms/blob/main/CONTRIBUTING.md)
for how to get started.

## Code of Conduct

This project follows the [Contributor Covenant Code of
Conduct](https://github.com/Bayesian-Graphical-Modelling-Lab/bgms/blob/main/CODE_OF_CONDUCT.md).


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

## Installation

Install from CRAN:

``` r
install.packages("bgms")
```

Or install the development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("Bayesian-Graphical-Modelling-Lab/bgms")
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

Contributions are welcome. See [CONTRIBUTING.md](CONTRIBUTING.md) for
how to get started.

## Code of Conduct

This project follows the [Contributor Covenant Code of
Conduct](CODE_OF_CONDUCT.md).

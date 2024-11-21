# bgms 0.1.4.2

Fixes a bug with the adaptation of the proposal variances.

# bgms 0.1.4.1

This is a minor release that adds some documentation and bug fixes.

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
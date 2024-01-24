# bgms 0.1.3

* Added the "edge_selection" option to the *bgm* function to allow Bayesian estimation without edge selection.
* Reduced file size when saving raw MCMC samples by changing row names.
* It is now possible to simulate data from the Blume-Capel ordinal MRF.

# bgms 0.1.2

This is a minor release adding fixing some minor bugs.


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
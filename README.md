# bgms
Bayesian methods for analyzing Markov Random Fields for binary and ordinal variables.

The project comprises several methods for the Bayesian analysis of a Markov Random Field (MRF) model for binary and ordinal variables. First, we combine an Expectation-Maximization variable selection approach with the continuous spike-and-slab prior set-up developed by Marsman, Huth, Waldorp, and Ntzoufras (2022; https://doi.org/10.1007/s11336-022-09848-8) for an MRF for binary variables. This approach aims to find a single (locally) optimal model for the data (e.g., Ročková & George; Journal of the American Statistical Association, 109(506):828-846, 2014). Second, we combine Gibbs sampling with a discrete spike-and-slab prior, using Gottardo and Raftery's (Journal of Computational and Graphical Statistics, 17(4):949-975; 2008) mixture of mutually singular distributions formulation of the discrete spike-and-slab prior and corresponding Metropolis algorithm. This approach allows for Bayesian model averaging that, in turn, will enable us to quantify structure uncertainty and the statistical evidence for local aspects of the model (e.g., conditional independence). A manuscript containing the details of this approach is forthcoming. Third, we offer procedures for optimizing the pseudolikelihood and pseudoposterior.

You can install the development version from GitHub with:

install.packages('devtools')

devtools::install_github('MaartenMarsman/bgms/tree/manuscript-version')


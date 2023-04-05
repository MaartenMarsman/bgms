#' bgms: Bayesian Analysis of Graphical Models
#'
#' @description
#' The \code{R} package \strong{bgms} provides tools for Bayesian analysis of
#' graphical models describing networks of variables. The package uses Bayesian
#' variable selection methods to model the underlying network structure.
#'
#' The package is organized around two Bayesian variable selection approaches:
#' (1) EM variable selection, and (2) Gibbs variable selection. The key
#' distinction is that the former uses a continuous spike and slab prior
#' distribution on the pairwise interactions
#' \insertCite{MarsmanEtAl_2022_objective}{bgms} that allows us to use EM
#' variable selection \insertCite{RockovaGeorge_2014}{bgms}. The Gibbs variable
#' selection approach \insertCite{GeorgeMcCulloch_1993}{bgms}, on the other
#' hand, stipulates a discrete spike and slab prior on the pairwise
#' interactions, which allows us to set the interactions to exact zeroes. To
#' account for the discontinuity at zero, we embed a Metropolis approach for
#' mixtures of mutually singular distributions
#' \insertCite{GottardoRaftery_2008}{bgms} in a Gibbs sampler.
#'
#' While the goal is to provide the above tools for Markov Random Field (MRF)
#' models for a wide range of variable types in the \strong{bgms} package, it
#' currently provides tools for analyzing networks of binary and/or ordinal
#' variables \insertCite{MarsmanHaslbeck_2023_OrdinalMRF}{bgms}. MRFs are a
#' special class of graphical models whose graph structure reflects the
#' conditional associations between their variables. Bayes factor tests on the
#' inclusion or exclusion of an edge in the MRF thus compare the conditional
#' dependence and conditional independence hypotheses for the corresponding
#' network variables \insertCite{HuthEtAl_2023_intro}{bgms}.
#'
#' The \strong{bgms} package offers several tools for analyzing the structure of
#' the MRF:
#'
#' \enumerate{
#'  \item Simulate response data from the MRF using the Gibbs sampler.
#'  \itemize{
#'    \item Simulate \code{\link{mrfSampler}}.
#'  }
#'
#'  \item Estimate the structure of the MRF using EM variable selection.
#'  \itemize{
#'    \item Network estimation with \code{\link{bgm.em}}.
#'  }
#
#'  \item Estimate the posterior distribution of the MRF's network structure
#'  using Gibbs variable selection.
#'  \itemize{
#'  \item Structure learning with \code{\link{bgm}}.
#'  }
#' }
#'
#' \bold{Additional Features}:
#'
#' \enumerate{
#'
#'  \item Optimizing the joint pseudolikelihood \code{\link{mple}}.
#'
#'  \item Optimizing the joint pseudoposterior \code{\link{mppe}}.
#' }
#'
#' @docType package
#' @keywords internal
#' @useDynLib bgms, .registration=TRUE
#' @references
#' \insertAllCited{}
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL


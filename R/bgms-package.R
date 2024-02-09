#' bgms: Bayesian Analysis of Networks of Binary and/or Ordinal Variables
#'
#' @description
#' The \code{R} package \strong{bgms} provides tools for Bayesian analysis of
#' graphical models describing networks of variables. The package uses Markov
#' chain Monte Carlo methods combined with a pseudolikelihood approach to
#' estimate the posterior distribution of model parameters.
#'
#' Gibbs variable selection \insertCite{GeorgeMcCulloch_1993}{bgms} is used to
#' model the underlying network structure of the graphical model. By imposing a
#' discrete spike and slab prior on the pairwise interactions, it is possible to
#' shrink the interactions to exactly zero. The Gibbs sampler embeds a
#' Metropolis approach for mixtures of mutually singular distributions
#' \insertCite{GottardoRaftery_2008}{bgms} to account for the discontinuity at
#' zero. The goal is to provide these tools for Markov Random Field (MRF) models
#' for a wide range of variable types in the \strong{bgms} package, and
#' it currently provides them for analyzing networks of binary and/or ordinal
#' variables \insertCite{MarsmanHaslbeck_2023_OrdinalMRF}{bgms}.
#'
#' While the goal is to provide the above tools for Markov Random Field (MRF)
#' models for a wide range of variable types in the \strong{bgms} package, it
#' currently provides tools for analyzing networks of binary and/or ordinal
#' variables \insertCite{MarsmanHaslbeck_2023_OrdinalMRF}{bgms}.
#'
#' MRFs are a special class of graphical models whose graph structure reflects
#' the conditional associations between their variables, making them useful for
#' testing for conditional independence or dependence. For example, the
#' inclusion Bayes factor tests for conditional independence or dependence of a
#' pair of variables in the network by comparing the predictive adequacy of
#' models that include the edge between these variables and models that exclude
#' the edge. \insertCite{HuthEtAl_2023_intro,SekulovskiEtAl_2023}{bgms}.
#'
#' The \strong{bgms} package offers several tools for analyzing the structure of
#' the MRF:
#'
#' \enumerate{
#'  \item Simulate response data from the MRF using the Gibbs sampler.
#'  \itemize{
#'    \item Simulate \code{\link{mrfSampler}}.
#'  }
#
#'  \item Estimate the posterior distribution of the MRF's parameters and
#'  possibly its network structure using Gibbs variable selection.
#'  \itemize{
#'  \item Bayesian estimation or Bayesian edge selection with \code{\link{bgm}}.
#'  }
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


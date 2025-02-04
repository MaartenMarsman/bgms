#' bgms: Bayesian Analysis of Networks of Binary and/or Ordinal Variables
#'
#' @description
#' The \code{R} package \strong{bgms} provides tools for Bayesian analysis of
#' the ordinal Markov random field, a graphical model describing a network of
#' binary and/or ordinal variables \insertCite{MarsmanVandenBerghHaslbeck_2024}{bgms}.
#' A pseudolikelihood is used to approximate the likelihood of the graphical
#' model, and Markov chain Monte Carlo methods are used to simulate from the
#' corresponding pseudoposterior distribution of the graphical model parameters.
#'
#' The \strong{bgm} function can be used for a one-sample design and the
#' \strong{bgmCompare} function can be used for a two-independent-samples design
#' \insertCite{MarsmanWaldorpSekulovskiHaslbeck_2024}{bgms}. Both functions can
#' model the selection of effects. In one-sample designs, the \strong{bgm}
#' function models the presence or absence of edges between pairs of variables
#' in the network. The estimated posterior inclusion probability indicates how
#' plausible it is that a network with an edge between the two corresponding
#' variables produced the observed data, and can be converted into a Bayes
#' factor test for conditional independence.
#'
#' In two-independent-samples designs, the \strong{bgmCompare} function models
#' the selection of group differences in edge weights and possibly category
#' thresholds. The estimated posterior inclusion probability indicates how
#' plausible it is that graphical models with a difference in the corresponding
#' edge weight or category threshold generated the data at hand, and can be
#' converted to a Bayes factor test for parameter equivalence.
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
#'  possibly its network structure in one-sample designs.
#'  \itemize{
#'  \item Bayesian estimation or Bayesian edge selection with \code{\link{bgm}}.
#'  }
#'
#'  \item Estimate the posterior distribution of the MRF's parameters in a
#'  two-independent-sample design, and possibly perform selection on group
#'  differences in MRF parameters.
#'  \itemize{
#'  \item Bayesian estimation or Bayesian difference selection with \code{\link{bgmCompare}}.
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


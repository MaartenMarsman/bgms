<!-- badges: start -->
[![CRAN Version](http://www.r-pkg.org/badges/version/bgms)](https://cran.r-project.org/package=bgms)
[![Downloads](https://cranlogs.r-pkg.org/badges/bgms)](https://cran.r-project.org/package=bgms)
[![Total](https://cranlogs.r-pkg.org/badges/grand-total/bgms)](https://cran.r-project.org/package=bgms)
<!-- badges: end -->

# bgms: Bayesian Analysis of Graphical Models

The `R` package <strong>bgms</strong> provides tools for Bayesian
analysis of graphical models describing networks of variables. The
package uses Markov chain Monte Carlo methods combined with a
pseudolikelihood approach to estimate the posterior distribution of
model parameters. Gibbs variable selection (George and McCulloch 1993)
is used to model the underlying network structure of the graphical
model. By imposing a discrete spike and slab prior on the pairwise
interactions, it is possible to shrink the interactions to exactly zero.
The Gibbs sampler embeds a Metropolis approach for mixtures of mutually
singular distributions (Gottardo and Raftery 2008) to account for the
discontinuity at zero. The goal is to provide these tools for Markov
Random Field (MRF) models for a wide range of variable types in the
<strong>bgms</strong> package, and it currently provides them for
analyzing networks of binary and/or ordinal variables (Marsman and
Haslbeck 2023).

## Why use Markov Random Fields?

Multivariate analysis using graphical models has received much attention
in the recent psychological and psychometric literature (Robinaugh et
al. 2020; Marsman and Rhemtulla 2022; Contreras et al. 2019). Most of
these graphical models are Markov Random Field (MRF) models, whose graph
structure reflects the conditional associations between variables
(Kindermann and Snell 1980). In these models, a missing edge between two
variables in the network implies that these variables are independent,
given the remaining variables (Lauritzen 2004). In other words, the
remaining variables of the network fully account for the potential
association between the unconnected variables.

## Why use a Bayesian approach to analyze the MRF?

Testing the structure of the MRF requires us to determine the
plausibility of the opposing hypotheses of conditional dependence and
conditional independence. That is, how plausible is it that the observed
data come from a network with a structure that includes the edge between
two variables compared to a network structure that excludes that edge?
Frequentist approaches are limited in this regard because they can only
reject the conditional independence hypothesis, not support it
(Wagenmakers et al. 2018; Wagenmakers 2007). This leads to the problem
that if an edge is excluded, we do not know whether this is because the
edge is absent in the population or because we lack the power to reject
the null hypothesis of independence. To avoid this problem, we will use
a Bayesian approach using Bayes factors (Kass and Raftery 1995)). The
inclusion Bayes factor (Huth et al. 2023; Sekulovski et al. 2024) allows
us to quantify how much the data support both conditional dependence
-<em>evidence of edge presence</em>- or conditional independence
-<em>evidence of edge absence</em>. It also allows us to conclude that
there is limited support for either hypothesis (Dienes 2014)-an
<em>absence of evidence</em>.

## Installation

The current developmental version can be installed with

``` r
if (!requireNamespace("remotes")) { 
  install.packages("remotes")   
}   
remotes::install_github("MaartenMarsman/bgms")
```

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-ContrerasEtAl_2019" class="csl-entry">

Contreras, A., I. Nieto, C. Valiente, R. Espinosa, and C. Vazquez. 2019.
“The Study of Psychopathology from the Network Analysis Perspective: A
Systematic Review.” *Psychotherapy and Psychosomatics* 88 (2): 71–83.
<https://doi.org/10.1159/000497425>.

</div>

<div id="ref-Dienes_2014" class="csl-entry">

Dienes, Z. 2014. “Using Bayes to Get the Most Out of Non-Significant
Results.” *Frontiers in Psychology* 5 (781): 1–17.
<https://doi.org/10.3389/fpsyg.2014.00781>.

</div>

<div id="ref-GeorgeMcCulloch_1993" class="csl-entry">

George, E. I., and R. E. McCulloch. 1993. “Variable Selection via Gibbs
Sampling.” *Journal of the American Statistical Association* 88 (423):
881–89. <https://doi.org/10.1080/01621459.1993.10476353>.

</div>

<div id="ref-GottardoRaftery_2008" class="csl-entry">

Gottardo, R., and A. E. Raftery. 2008. “Markov Chain Monte Carlo with
Mixtures of Mutually Singular Distributions.” *Journal of Computational
and Graphical Statistics* 17 (4): 949–75.
<https://doi.org/10.1198/106186008X386102>.

</div>

<div id="ref-HuthEtAl_2023_intro" class="csl-entry">

Huth, K., J. de Ron, A. E. Goudriaan, K. Luigjes, R. Mohammadi, R. J.
van Holst, E.-J. Wagenmakers, and M. Marsman. 2023. “Bayesian Analysis
of Cross-Sectional Networks: A Tutorial in R and JASP.” *Advances in
Methods and Practices in Psychological Science* 6 (4): 1–18.
<https://doi.org/10.1177/25152459231193334>.

</div>

<div id="ref-KassRaftery_1995" class="csl-entry">

Kass, R. E., and A. E. Raftery. 1995. “Bayes Factors.” *Journal of the
American Statistical Association* 90 (430): 773–95.
<https://doi.org/10.2307/2291091>.

</div>

<div id="ref-KindermannSnell1980" class="csl-entry">

Kindermann, R., and J. L. Snell. 1980. *Markov Random Fields and Their
Applications*. Vol. 1. Contemporary Mathematics. Providence: American
Mathematical Society.

</div>

<div id="ref-Lauritzen2004" class="csl-entry">

Lauritzen, S. L.. 2004. *Graphical Models*. Oxford: Oxford University
Press.

</div>

<div id="ref-MarsmanHaslbeck_2023_OrdinalMRF" class="csl-entry">

Marsman, M., and J. M. B. Haslbeck. 2023. “Bayesian Analysis of the
Ordinal Markov Random Field.” *PsyArXiv*.
<https://doi.org/10.31234/osf.io/ukwrf>.

</div>

<div id="ref-MarsmanRhemtulla_2022_SIintro" class="csl-entry">

Marsman, M., and M. Rhemtulla. 2022. “Guest Editors’ Introduction to the
Special Issue ‘Network Psychometrics in Action’: Methodological
Innovations Inspired by Empirical Problems.” *Psychometrika* 87 (1):
1–11. <https://doi.org/10.1007/s11336-022-09861-x>.

</div>

<div id="ref-RobinaughEtAl_2020" class="csl-entry">

Robinaugh, D. J., R. H. A. Hoekstra, E. R. Toner, and D. Borsboom. 2020.
“The Network Approach to Psychopathology: A Review of the Literature
2008–2018 and an Agenda for Future Research.” *Psychological Medicine*
50: 353–66. <https://doi.org/10.1017/S0033291719003404>.

</div>

<div id="ref-SekulovskiEtAl_2023" class="csl-entry">

Sekulovski N, Keetelaar S, Huth K, Wagenmakers E.-J, van Bork R, van den Bergh D, Marsman M (2024). “Testing conditional independence in psychometric networks: An analysis of three bayesian methods.” *Multivariate Behavioral Research*, 1–21. <doi:10.1080/00273171.2024.2345915>.

</div>

<div id="ref-Wagenmakers_2007" class="csl-entry">

Wagenmakers, E.-J. 2007. “A Practical Solution to the Pervasive Problems
of p Values.” *Psychonomic Bulletin & Review* 14: 779–804.
<https://doi.org/10.3758/BF03194105>.

</div>

<div id="ref-WagenmakersEtAl_2018_BIP1" class="csl-entry">

Wagenmakers, E.-J., M. Marsman, T. Jamil, A. Ly, J. Verhagen, J. Love,
R. Selker, et al. 2018. “Bayesian Inference for Psychology. Part I:
Theoretical Advantages and Practical Ramifications.” *Psychonomic
Bulleting & Review* 25 (1): 58–76.
<https://doi.org/10.3758/s13423-017-1343-3>.

</div>

</div>

<!-- badges: start -->
[![CRAN Version](http://www.r-pkg.org/badges/version/bgms)](https://cran.r-project.org/package=bgms)
[![Average](https://cranlogs.r-pkg.org/badges/bgms)](https://cran.r-project.org/package=bgms)
[![Total](https://cranlogs.r-pkg.org/badges/grand-total/bgms)](https://cran.r-project.org/package=bgms)
<!-- badges: end -->

# bgms: Bayesian Analysis of Graphical Models

The `R` package <strong>bgms</strong> provides tools for a Bayesian
analysis of graphical models describing networks of variables. The
package uses Bayesian variable selection methods to model the underlying
network structure. The methods are organized around two general
approaches for Bayesian variable selection: (1) EM variable selection
and (2) Gibbs variable selection. The key distinction is that the former
uses a continuous spike and slab prior distribution on the pairwise
interactions (Marsman et al. 2022) that allows us to use EM variable
selection (Ročková and George 2014). The Gibbs variable selection
approach (George and McCulloch 1993), on the other hand, stipulates a
discrete spike and slab prior on the pairwise interactions, which allows
us to set the interactions to exact zeroes. To account for the
discontinuity at zero, we embed a Metropolis approach for mixtures of
mutually singular distributions (Gottardo and Raftery 2008) in a Gibbs
sampler. The goal is to provide these tools for Markov Random Field
(MRF) models for a wide range of variable types in the
<strong>bgms</strong> package, and it currently provides them for
analyzing networks of binary and/or ordinal variables (Marsman and
Haslbeck 2023).

## Why use Markov Random Fields?

Multivariate analysis using graphical models has received much attention
in the recent psychological and psychometric literature (Robinaugh et
al. 2020; Marsman and Rhemtulla 2022; Steinley 2021; Contreras et al.
2019). Most of these graphical models are Markov Random Field (MRF)
models, whose graph structure reflects the conditional associations
between variables (Kindermann and Snell 1980). In these models, a
missing edge between two variables in the network implies that these
variables are independent, given the remaining variables (Lauritzen
2004). In other words, the remaining variables of the network fully
account for the potential association between the unconnected variables.

## Why use a Bayesian approach to analyze the MRF?

Testing the structure of the MRF requires us to determine the
plausibility of the opposing hypotheses of conditional dependence and
conditional independence. For example, how plausible are network
structures that include the edge between variables 3 and 9 compared to
network structures that exclude this edge? Frequentist approaches are
limited in this respect, because they can only reject the conditional
independence hypothesis, but not support it (Wagenmakers et al. 2018;
Wagenmakers 2007). This creates the problem that, if an edge is
excluded, we do not know whether this is because the edge is absent in
the population, or because we lack the power to reject the null
hypothesis of independence. To avoid this problem, we will use a
Bayesian approach using Bayes factors (Kass and Raftery 1995). The
inclusion Bayes factor (Huth et al. 2023) allows us to quantify how much
the data support both conditional dependence —<em>evidence of edge
presence</em>— or conditional independence —<em>evidence of edge
absence</em>. It also allows us to conclude that there is only limited
support for either hypothesis (Dienes 2014) —an <em>absence of
evidence</em>.

## Installation

You can install the latest version from CRAN using:

``` r
install.packages("bgms")
```

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

Dienes, Z. 2014. “Using Bayes to Get the Most Our of Non-Significant
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
of Cross-Sectional Networks: A Tutorial in R and JASP.” *PsyArXiv*.
<https://doi.org/10.31234/osf.io/ub5tc>.

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

<div id="ref-MarsmanEtAl_2022_objective" class="csl-entry">

Marsman, M., K. B. S. Huth, L. J. Waldorp, and I. Ntzoufras. 2022.
“Objective Bayesian Edge Screening and Structure Selection for Ising
Networks.” *Psychometrika* 87 (1): 47–82.
<https://doi.org/10.1007/s11336-022-09848-8>.

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

<div id="ref-RockovaGeorge_2014" class="csl-entry">

Ročková, V., and E. I. George. 2014. “EMVS: The EM Approach to Bayesian
Variable Selection.” *Journal of the American Statistical Association*
109 (506): 828–46. <https://doi.org/10.1080/01621459.2013.869223>.

</div>

<div id="ref-Steinley_2021_SIintro" class="csl-entry">

Steinley, D. 2021. “Recent Advances in (Graphical) Network Models.”
*Multivariate Behavioral Research* 56 (2): 171–74.
<https://doi.org/10.1080/00273171.2021.1911777>.

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

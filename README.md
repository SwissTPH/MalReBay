
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MalReBay

<!-- badges: start -->

<!-- badges: end -->

MalReBay is an R package designed to analyze genotyping data from
Therapeutic Efficacy Studies (TES) and classify recurrent malaria
parasite infections. It determines the posterior probability that a
given infection is a recrudescence (a treatment failure) versus a new
infection. Unlike simple allele-matching methods, MalReBay implements a
Bayesian algorithm, solved using a Markov Chain Monte Carlo (MCMC)
engine, to provide a robust probabilistic classification for each
sample.

It supports analysis of both traditional length-polymorphic markers
(microsatellites, msp1, msp2, and glurp) and modern amplicon sequencing
data. For length-polymorphic data (microsatellites, msp), it employs a
distance-based likelihood model that accounts for PCR stutter and
scoring errors. For amplicon sequencing data, it uses an exact-match
categorical model, which is suitable for discrete haplotypes. The
framework robustly handles polyclonal infections and imputes missing
genotypes as part of the MCMC sampling, providing a powerful tool for
modern TES data analysis.

## Installation

You can install the development version of MalReBay from
[GitHub](https://github.com/SwissTPH/MalReBay) with:

``` r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("SwissTPH/MalReBay")
```

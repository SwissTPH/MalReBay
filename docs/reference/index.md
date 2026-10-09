# Package index

## Full pipeline

Run the complete analysis, from importing the genotyping data to saving
the posterior probabilities and plots, in a single call.

- [`MalReBay()`](https://swisstph.github.io/MalReBay/reference/MalReBay.md)
  : Run the MalReBay malaria recrudescence classification pipeline

## Pipeline steps

Run the individual stages of the analysis separately, for example to
inspect intermediate results or to rerun the classification with
different settings.

- [`import_data()`](https://swisstph.github.io/MalReBay/reference/import_data.md)
  : Import Genotyping Data from Excel
- [`classify_infections()`](https://swisstph.github.io/MalReBay/reference/classify_infections.md)
  : Classify malaria infections using Bayesian MCMC
- [`summarise_results()`](https://swisstph.github.io/MalReBay/reference/summarise_results.md)
  : Summarise MCMC Classification Results
- [`save_results()`](https://swisstph.github.io/MalReBay/reference/save_results.md)
  : Save Classification Results and Generate Plots

## Data quality checks

Check whether the imported data are suitable for the analysis.

- [`data_quality_check()`](https://swisstph.github.io/MalReBay/reference/data_quality_check.md)
  : Check dataset viability per site
- [`compute_locus_comparability()`](https://swisstph.github.io/MalReBay/reference/compute_locus_comparability.md)
  : Compute locus comparability per patient
- [`detect_msp_variants()`](https://swisstph.github.io/MalReBay/reference/detect_msp_variants.md)
  : Identify which locinames belong to the MSP1/MSP2 families

## Convergence diagnostics

Assess whether the MCMC sampling has converged before interpreting the
results.

- [`check_mcmc_diagnostics()`](https://swisstph.github.io/MalReBay/reference/check_mcmc_diagnostics.md)
  : Check MCMC sampler diagnostics and report plain-language guidance
- [`plot_likelihood_diagnostics()`](https://swisstph.github.io/MalReBay/reference/plot_likelihood_diagnostics.md)
  : Plot MCMC Likelihood Diagnostics

## Plots

Visualise the input data and the classification results.

- [`plot_probability_histogram()`](https://swisstph.github.io/MalReBay/reference/plot_probability_histogram.md)
  : Plot Posterior Probability Histogram
- [`plot_comparison_heatmap()`](https://swisstph.github.io/MalReBay/reference/plot_comparison_heatmap.md)
  : Plot Match-Counting vs MalReBay Comparison Heatmap
- [`plot_moi()`](https://swisstph.github.io/MalReBay/reference/plot_moi.md)
  : Plot Multiplicity of Infection (MOI)
- [`plot_markers_diversity()`](https://swisstph.github.io/MalReBay/reference/plot_markers_diversity.md)
  : Generate and Save Diversity Pie Charts
- [`plot_allele_distribution()`](https://swisstph.github.io/MalReBay/reference/plot_allele_distribution.md)
  : Plot Allele Distribution

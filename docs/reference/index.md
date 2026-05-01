# Package index

## Functions

Run the full Bayesian recrudescence classification pipeline, from
importing genotype data to saving posterior probabilities and diagnostic
plots.

- [`MalReBay()`](https://swisstph.github.io/MalReBay/reference/MalReBay.md)
  : Run the MalReBay malaria recrudescence classification pipeline

## Data import

Import and validate paired genotyping data from an Excel file,
automatically detecting whether the data is length-polymorphic or
amplicon sequencing.

- [`import_data()`](https://swisstph.github.io/MalReBay/reference/import_data.md)
  : Import Genotyping Data from Excel

## Plotting helpers

Visualize results to help you tune hyperparameters and choose model
methods.

- [`plot_likelihood_diagnostics()`](https://swisstph.github.io/MalReBay/reference/plot_likelihood_diagnostics.md)
  : Plot MCMC Likelihood Diagnostics
- [`plot_probability_histogram()`](https://swisstph.github.io/MalReBay/reference/plot_probability_histogram.md)
  : Plot Posterior Probability Histogram
- [`plot_moi()`](https://swisstph.github.io/MalReBay/reference/plot_moi.md)
  : Plot Multiplicity of Infection (MOI)
- [`plot_markers_diversity()`](https://swisstph.github.io/MalReBay/reference/plot_markers_diversity.md)
  : Generate and Save Diversity Pie Charts

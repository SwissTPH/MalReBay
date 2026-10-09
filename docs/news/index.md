# Changelog

## MalReBay 0.0.1

This is the first release version of MalReBay! 🎉

- Added a `NEWS.md` file to track changes to the package.

- Functions:

  - [`MalReBay()`](https://swisstph.github.io/MalReBay/reference/MalReBay.md)
  - [`import_data()`](https://swisstph.github.io/MalReBay/reference/import_data.md)
  - [`plot_likelihood_diagnostics()`](https://swisstph.github.io/MalReBay/reference/plot_likelihood_diagnostics.md)
  - [`plot_probability_histogram()`](https://swisstph.github.io/MalReBay/reference/plot_probability_histogram.md)
  - [`plot_moi()`](https://swisstph.github.io/MalReBay/reference/plot_moi.md)
  - [`plot_markers_diversity()`](https://swisstph.github.io/MalReBay/reference/plot_markers_diversity.md)

- Example data:

  - **`Angola_2021_TES_7NMS.xlsx`** — genotyping data from a Therapeutic
    Efficacy Study conducted in Angola in 2021. Contains 70 patients
    from three sites (Benguela, Lunda Sul, Zaire) genotyped at 7
    microsatellite markers.

  - **`Amplicon_Sequencing.xlsx`** — example amplicon sequencing dataset
    for demonstrating the ampseq workflow.

  - **`makers_details.xlsx`** — marker metadata file containing marker
    IDs, binning methods, and repeat lengths required by
    [`import_data()`](https://swisstph.github.io/MalReBay/reference/import_data.md).

  - **`default_mcmc_config.xlsx`** — default MCMC configuration file
    with recommended settings for `n_chains`, `iter`, `burn_in_frac`,
    `adapt_delta`, and `random_seed`.

  All files can be accessed via:

  ``` r
  system.file("extdata", "filename.xlsx", package = "MalReBay")
  ```

- Vignettes:

  - [Getting
    started](https://SwissTPH.github.io/MalReBay/articles/MalReBay.html)

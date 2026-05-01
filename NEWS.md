# MalReBay 0.0.1

This is the first release version of MalReBay! 🎉

- Added a `NEWS.md` file to track changes to the package.
- Functions:
    - `MalReBay()`
    - `import_data()`
    - `plot_likelihood_diagnostics()`
    - `plot_probability_histogram()`
    - `plot_moi()`
    - `plot_markers_diversity()`
- Example data:
    - **`Angola_2021_TES_7NMS.xlsx`** — genotyping data from a Therapeutic Efficacy 
  Study conducted in Angola in 2021. Contains 70 patients from three sites 
  (Benguela, Lunda Sul, Zaire) genotyped at 7 microsatellite markers.

  - **`Amplicon_Sequencing.xlsx`** — example amplicon sequencing dataset for 
  demonstrating the ampseq workflow.

  - **`makers_details.xlsx`** — marker metadata file containing marker IDs, 
  binning methods, and repeat lengths required by `import_data()`.

  - **`default_mcmc_config.xlsx`** — default MCMC configuration file with 
  recommended settings for `n_chains`, `iter`, `burn_in_frac`, `adapt_delta`, 
  and `random_seed`.

  All files can be accessed via:
  ```r
  system.file("extdata", "filename.xlsx", package = "MalReBay")
  ```
- Vignettes:
    - [Getting started](https://SwissTPH.github.io/MalReBay/articles/MalReBay.html)
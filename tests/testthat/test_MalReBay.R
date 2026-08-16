library(testthat)
library(MalReBay)

zaire_imported <- readRDS(system.file("extdata", "imported_data.rds", package = "MalReBay"))
marker_file    <- system.file("extdata", "makers_details.xlsx", package = "MalReBay")

zaire_file <- tempfile(fileext = ".xlsx")
writexl::write_xlsx(
  list(Sheet1 = zaire_imported$late_failures, Sheet2 = zaire_imported$additional),
  zaire_file
)

mcmc_config <- list(
  n_chains     = 2,
  iter         = 1000,
  burn_in_frac = 0.5,
  random_seed  = 42,
  adapt_delta  = 0.8
)

cmdstan_ok <- tryCatch({
  path <- cmdstanr::cmdstan_path()
  nzchar(path) && file.exists(path)
}, error = function(e) FALSE)

if (cmdstan_ok) {
  tmp_out <- tempfile()
  result  <- suppressMessages(suppressWarnings(
    MalReBay(
      filepath        = zaire_file,
      marker_filepath = marker_file,
      mcmc_config     = mcmc_config,
      output_folder   = tmp_out,
      n_workers       = 1,
      verbose         = FALSE
    )
  ))
}

test_that("MalReBay returns correct list structure for Zaire", {
  skip_if_not(cmdstan_ok, "CmdStan not installed")
  expect_named(result, c("posterior_probabilities", "comparison", "convergence", "mcmc_loglikelihoods"))
})

test_that("MalReBay probabilities are between 0 and 1", {
  skip_if_not(cmdstan_ok, "CmdStan not installed")
  probs <- result$posterior_probabilities$Probability
  expect_true(all(probs >= 0 & probs <= 1, na.rm = TRUE))
})

test_that("MalReBay saves all expected output files", {
  skip_if_not(cmdstan_ok, "CmdStan not installed")
  expect_true(file.exists(file.path(tmp_out, "posterior_probabilities.csv")))
  expect_true(file.exists(file.path(tmp_out, "bayesian_match_counting_comparison.csv")))
  expect_true(file.exists(file.path(tmp_out, "diversity_length_polymorphic_comparison.png")))
  expect_gt(length(list.files(tmp_out, pattern = "^moi_per_marker_.*\\.png$")), 0)
  expect_true(file.exists(file.path(tmp_out, "recrudescence_probability_histogram.png")))
  expect_true(dir.exists(file.path(tmp_out, "convergence_diagnosis")))
})

test_that("MalReBay non-existent filepath throws error", {
  expect_error(MalReBay(filepath = "nonexistent.xlsx"), regexp = "Cannot read file")
})

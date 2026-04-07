library(testthat)
library(MalReBay)
# Shared fixtures loaded automatically from helper.R:
#   angola_data_path, angola_marker_path, fast_mcmc

# Helper — runs the full MalReBay() pipeline on real Angola data.
# angola_data_path / angola_marker_path resolve via system.file() against
# the installed package. If the package is not installed these will be ""
# and the skip_on_cran() + path guard will skip the tests cleanly.
run_pipeline <- function(tmp) {
  suppressMessages(suppressWarnings(
    MalReBay(
      filepath        = angola_data_path,
      marker_filepath = angola_marker_path,
      mcmc_config     = fast_mcmc,
      output_folder   = tmp,
      n_workers       = 1,
      verbose         = FALSE
    )
  ))
}

# Guard — skips all pipeline tests when Angola files are not reachable.
# This happens during devtools::test() before devtools::install().
skip_if_angola_missing <- function() {
  skip_if(
    !nzchar(angola_data_path) || !file.exists(angola_data_path),
    "Angola data files not found — run devtools::install() first"
  )
}

# ============================================================
# Full pipeline — real Angola data
# ============================================================

test_that("MalReBay: runs on built-in data without error", {
  skip_on_cran()
  skip_if_angola_missing()
  tmp <- tempfile()
  on.exit(unlink(tmp, recursive = TRUE))
  
  expect_no_error(run_pipeline(tmp))
})

test_that("MalReBay: returns correct list structure", {
  skip_on_cran()
  skip_if_angola_missing()
  tmp <- tempfile()
  on.exit(unlink(tmp, recursive = TRUE))
  
  result <- run_pipeline(tmp)
  expect_named(result, c("posterior_probabilities", "comparison",
                         "convergence", "mcmc_loglikelihoods"))
})

test_that("MalReBay: probabilities are between 0 and 1", {
  skip_on_cran()
  skip_if_angola_missing()
  tmp <- tempfile()
  on.exit(unlink(tmp, recursive = TRUE))
  
  probs <- run_pipeline(tmp)$posterior_probabilities$Probability
  expect_true(all(probs >= 0 & probs <= 1, na.rm = TRUE))
})

# ============================================================
# Output files
# ============================================================

test_that("MalReBay: saves all expected output files", {
  skip_on_cran()
  skip_if_angola_missing()
  tmp <- tempfile()
  on.exit(unlink(tmp, recursive = TRUE))
  
  run_pipeline(tmp)
  
  expect_true(file.exists(file.path(tmp, "posterior_probabilities.csv")))
  expect_true(file.exists(file.path(tmp, "bayesian_match_counting_comparison.csv")))
  expect_true(file.exists(file.path(tmp, "diversity_length_polymorphic_comparison.png")))
  expect_true(file.exists(file.path(tmp, "moi_per_marker_by_site.png")))
  expect_true(file.exists(file.path(tmp, "recrudescence_probability_histogram.png")))
  expect_true(dir.exists( file.path(tmp, "convergence_diagnosis")))
})

# ============================================================
# Bad input
# ============================================================

test_that("MalReBay: non-existent filepath throws error", {
  expect_error(
    MalReBay(filepath = "nonexistent.xlsx"),
    regexp = "Cannot read file"
  )
})
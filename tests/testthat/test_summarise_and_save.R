library(testthat)
library(MalReBay)
# Shared fixtures loaded automatically from helper.R:
#   create_mock_imported_data(), fast_mcmc,
#   make_test_mcmc_results(), make_test_summary_results()
#
# Note: make_test_summary_results() always passes a real output_folder
# to summarise_results() because plot_likelihood_diagnostics() requires one.
# When called with output_folder = NULL it uses a disposable tempfile().
# When called with an explicit path it writes files there for inspection.

# ============================================================
# summarise_results() — bad input
# ============================================================

test_that("summarise_results: NULL mcmc_results throws error", {
  expect_error(
    summarise_results(NULL, create_mock_imported_data()),
    regexp = "mcmc_results.*is empty"
  )
})

# ============================================================
# summarise_results() — output structure
# ============================================================

test_that("summarise_results: returns correct list keys", {
  result <- make_test_summary_results()
  expect_named(result, c("posterior_probabilities", "comparison",
                         "convergence", "mcmc_loglikelihoods"))
})

test_that("summarise_results: posterior_probabilities has correct columns", {
  result <- make_test_summary_results()
  expect_true(all(c("Sample.ID", "Site", "Probability", "N_Comparable_Loci") %in%
                    colnames(result$posterior_probabilities)))
})

test_that("summarise_results: probabilities are between 0 and 1", {
  probs <- make_test_summary_results()$posterior_probabilities$Probability
  expect_true(all(probs >= 0 & probs <= 1, na.rm = TRUE))
})

test_that("summarise_results: one row per patient in posterior_probabilities", {
  result <- make_test_summary_results()
  expect_equal(nrow(result$posterior_probabilities), 3)
  expect_equal(
    sort(result$posterior_probabilities$Sample.ID),
    sort(c("BD21-002", "BD21-040", "BD21-041"))
  )
})

test_that("summarise_results: comparison table contains all late_failures rows", {
  result <- make_test_summary_results()
  # 3 pairs x 2 timepoints = 6 rows
  expect_equal(nrow(result$comparison), 6)
  expect_true("Probability" %in% colnames(result$comparison))
})

test_that("summarise_results: convergence is NULL or data.frame", {
  # make_test_summary_results() always provides an output_folder so
  # convergence diagnostics will always be attempted — result is either
  # a data.frame (converged) or NULL (too few samples to diagnose)
  result <- make_test_summary_results()
  expect_true(is.null(result$convergence) || is.data.frame(result$convergence))
})

test_that("summarise_results: convergence plots saved when output_folder provided", {
  skip_on_cran()
  tmp <- tempfile()
  on.exit(unlink(tmp, recursive = TRUE))
  
  make_test_summary_results(output_folder = tmp)
  
  diag_dir   <- file.path(tmp, "convergence_diagnosis", "Benguela")
  plot_files <- list.files(diag_dir, pattern = "\\.png$")
  expect_true(dir.exists(diag_dir))
  expect_gt(length(plot_files), 0)
})

# ============================================================
# save_results() — bad input
# ============================================================

test_that("save_results: invalid summary_results throws error", {
  expect_error(
    save_results(list(wrong = "structure")),
    regexp = "summary_results.*must be a valid list"
  )
})

# ============================================================
# save_results() — file output
# ============================================================

test_that("save_results: creates output folder if it does not exist", {
  tmp    <- file.path(tempdir(), "new_folder_test")
  result <- make_test_summary_results()
  on.exit(unlink(tmp, recursive = TRUE))
  
  save_results(result, imported_data = create_mock_imported_data(), output_folder = tmp, verbose = FALSE)
  expect_true(dir.exists(tmp))
})

test_that("save_results: writes posterior_probabilities.csv and comparison.csv", {
  tmp    <- tempfile()
  result <- make_test_summary_results()
  on.exit(unlink(tmp, recursive = TRUE))
  
  save_results(result, output_folder = tmp, verbose = FALSE)
  
  expect_true(file.exists(file.path(tmp, "posterior_probabilities.csv")))
  expect_true(file.exists(file.path(tmp, "bayesian_match_counting_comparison.csv")))
})

test_that("save_results: convergence CSV written only when convergence is not NULL", {
  tmp    <- tempfile()
  result <- make_test_summary_results()
  on.exit(unlink(tmp, recursive = TRUE))
  
  save_results(result, output_folder = tmp, verbose = FALSE)
  
  cv_path <- file.path(tmp, "mcmc_convergence_summary.csv")
  if (!is.null(result$convergence)) expect_true(file.exists(cv_path))
  else                               expect_false(file.exists(cv_path))
})

test_that("save_results: saved CSVs contain correct columns", {
  tmp    <- tempfile()
  result <- make_test_summary_results()
  on.exit(unlink(tmp, recursive = TRUE))
  
  save_results(result, output_folder = tmp, verbose = FALSE)
  
  pp <- utils::read.csv(file.path(tmp, "posterior_probabilities.csv"))
  ct <- utils::read.csv(file.path(tmp, "bayesian_match_counting_comparison.csv"))
  
  expect_true(all(c("Sample.ID", "Probability") %in% colnames(pp)))
  expect_true("Probability" %in% colnames(ct))
})

test_that("save_results: returns named character vector of existing file paths", {
  tmp    <- tempfile()
  result <- make_test_summary_results()
  on.exit(unlink(tmp, recursive = TRUE))
  
  paths <- save_results(result, output_folder = tmp, verbose = FALSE)
  
  expect_type(paths, "character")
  expect_true("posterior_probabilities" %in% names(paths))
  expect_true("comparison"              %in% names(paths))
  expect_true(all(file.exists(paths)))
})
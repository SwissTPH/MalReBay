
example_file <- system.file("extdata", "Angola_2021_TES_7NMS.xlsx", package = "MalReBay")
marker_file  <- system.file("extdata", "marker_details.xlsx", package = "MalReBay")

# Helper to ensure data is processed correctly for tests
# We wrap this in a helper or keep it clean to avoid failures before tests start
get_test_data <- function() {
  testthat::skip_if(example_file == "" || !file.exists(example_file), "Example data file not found.")
  testthat::skip_if(marker_file == "" || !file.exists(marker_file), "Marker file not found.")
  
  import_data(
    filepath = example_file, 
    marker_filepath = marker_file, 
    verbose = FALSE
  )
}

test_that("classify_infections runs and returns correctly structured results", {
  processed_data <- get_test_data()
  
  quick_mcmc_config <- list(
    n_chains = 2,
    max_iterations = 100, # Reduced for faster testing
    chunk_size = 100,
    rhat_threshold = 1.1,
    ess_threshold = 50
  )
  
  # Call the main function
  results_list <- classify_infections(
    imported_data = processed_data,
    mcmc_config = quick_mcmc_config,
    output_folder = tempdir(),
    verbose = FALSE 
  )
  
  # Assertions
  expect_type(results_list, "list")
  expect_named(results_list, c("summary", "marker_details", "comparison", "mcmc_loglikelihoods"))
  
  final_summary <- results_list$summary
  expect_s3_class(final_summary, "data.frame")
  expect_gt(nrow(final_summary), 0)
  expect_true(all(final_summary$Probability >= 0 & final_summary$Probability <= 1))
})

test_that("classify_infections removes pre-existing convergence_diagnosis folder contents", {
  processed_data <- get_test_data()
  
  # ARRANGE: Create a temporary directory and an old folder with a file in it
  temp_output_dir <- file.path(tempdir(), "test_diag_cleanup")
  if (!dir.exists(temp_output_dir)) dir.create(temp_output_dir)
  
  old_convergence_dir <- file.path(temp_output_dir, "convergence_diagnosis")
  if (!dir.exists(old_convergence_dir)) dir.create(old_convergence_dir)
  
  old_file <- file.path(old_convergence_dir, "old_plot.png")
  file.create(old_file)
  
  # Verify setup
  expect_true(file.exists(old_file))
  
  # Use a minimal config to make this test run very fast
  quick_mcmc_config <- list(n_chains = 1, max_iterations = 5, chunk_size = 5)
  
  # ACT
  classify_infections(
    imported_data = processed_data,
    mcmc_config = quick_mcmc_config,
    output_folder = temp_output_dir,
    verbose = FALSE
  )
  
  expect_false(file.exists(old_file))
  
  # Cleanup
  unlink(temp_output_dir, recursive = TRUE)
})
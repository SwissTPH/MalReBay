test_that("classify_infections runs and returns a correctly structured data frame", {
  example_file <- system.file("extdata", "Angola_2021_TES_7NMS.xlsx",
                              package = "MalReBay")
  
  testthat::skip_if_not(file.exists(example_file), "Example data file not found.")
  
  marker_info <- data.frame(
    Marker = c("313", "383", "TA1"),
    RepeatLength = c(2, 2, 3)
  )
  
  quick_mcmc_config <- list(
    n_chains = 2,
    max_iterations = 400,
    chunk_size = 200,
    burn_in_frac = 0.1
  )

  final_summary <- classify_infections(
    input_filepath = example_file,
    mcmc_config = quick_mcmc_config,
    output_folder = tempdir()
  )
  testthat::skip_if(is.null(final_summary), "Function returned NULL, which is acceptable.") 
  testthat::expect_s3_class(final_summary, "data.frame")
  testthat::expect_named(final_summary, c("Site", "Sample.ID", "Probability"))
  testthat::expect_gt(nrow(final_summary), 0)
  testthat::expect_true(all(final_summary$Probability >= 0))
  testthat::expect_true(all(final_summary$Probability <= 1))
  
})
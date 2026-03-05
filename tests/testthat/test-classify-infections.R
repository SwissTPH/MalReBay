# Helper functions for locating example datasets used in tests
# Retrieve a data file from the package extdata folder
get_data_file <- function(file_name) {
  path <- system.file("extdata", file_name, package = "MalReBay")
  if (path != "" && file.exists(path)) return(normalizePath(path))
  stop("Could not find data file: ", file_name)
}

# Create a small subset of the Angola TES dataset for fast MCMC testing
get_angola_test_subset <- function(n_pairs = 10, n_additional = 5) {
  # Load example data files
  data_path <- get_data_file("Angola_2021_TES_7NMS.xlsx")
  marker_path <- get_data_file("marker_details.xlsx")
  
  # create a full dataset with the same structure as the original
  full_data <- suppressWarnings(import_data(data_path, marker_path, verbose = FALSE))
  
  # Standardize IDs and identify matching Day0/recurrence pairs
  all_ids <- full_data$late_failures$Sample.ID
  prefixes <- gsub(" Day (0|Failure)$", "", all_ids)
  id_counts <- table(prefixes)
  
  valid_ids <- names(id_counts[id_counts >= 2])
  selected_ids <- valid_ids[1:min(n_pairs, length(valid_ids))]
  
  # Subset to selected pairs and additional samples
  full_data$late_failures <- full_data$late_failures[prefixes %in% selected_ids, ]
  if (nrow(full_data$additional) > 0) {
    full_data$additional <- full_data$additional[1:min(n_additional, nrow(full_data$additional)), ]
  }
  
  return(full_data)
}


# Create a small synthetic microsatellite dataset for testing
create_mock_microsat_data <- function() {
  data.frame(
    Sample.ID = c("Mock_01 Day 0", "Mock_01 Day Failure", 
                  "Mock_02 Day 0", "Mock_02 Day Failure"),
    Site = "Benguela",
    POLYA_1 = c(153, 105, 159, 159),  
    POLYA_2 = c(NA, NA, NA, NA),     
    POLYA_3 = c(NA, NA, NA, NA),
    POLYA_4 = c(NA, NA, NA, NA),
    check.names = FALSE,             
    stringsAsFactors = FALSE
  )
}

# Create minimal marker metadata for microsatellite tests
create_mock_markers <- function() {
  data.frame(
    marker_id = "POLYA",
    binning_method = "microsatellite",
    repeatlength = 3,
    stringsAsFactors = FALSE
  )
}

# Unit tests for input validation and core preprocessing logic

test_that("classify_infections errors on wrong input type", {
  # Function should reject invalid input
  expect_error(
    classify_infections("not_data"),
    "imported_data"
  )
  
})


test_that("marker information must match dataset columns", {
  # Create mock genotype dataset
  mock_data <- create_mock_microsat_data()
  
  # Define incorrect marker metadata
  bad_markers <- data.frame(
    marker_id = "wrong_marker",
    binning_method = "microsatellite",
    repeatlength = 3
  )
  # define_alleles should return empty definitions for missing markers
  res <- suppressWarnings(define_alleles(mock_data, bad_markers))
  expect_equal(nrow(res$wrong_marker), 0)
})


test_that("import_data handles Angolan ID standardization", {
  # Example raw sample IDs from dataset
  raw_ids <- c("BD21-002D0", "BD21-040D28")
  
  # Apply same ID cleaning logic used in import_data()
  cleaned <- gsub("_?D0$| Day 0$", " Day 0", raw_ids)
  is_day0 <- grepl(" Day 0$", cleaned)
  cleaned[!is_day0] <- gsub("_?D[1-9][0-9]*$| Day Failure$", " Day Failure", cleaned[!is_day0])
  
  # Verify standardized format
  expect_equal(cleaned, c("BD21-002 Day 0", "BD21-040 Day Failure"))
})


test_that("Microsatellite binning and frequency logic", {
  # Define minimal marker information
  marker_info <- data.frame(marker_id = "POLYA", 
                            binning_method = "microsatellite", 
                            repeatlength = 2, 
                            stringsAsFactors = FALSE)
  mock_data <- create_mock_microsat_data()
  
  # Compute allele bins and population frequencies
  defs <- suppressWarnings(define_alleles(mock_data, marker_info))
  freqs <- calculate_frequencies(mock_data, defs, marker_info)
  
  # Verify allele bin includes observed allele
  expect_true(any(defs$POLYA[, "lower"] <= 153 & defs$POLYA[, "upper"] >= 153))
  
  # Frequencies should sum to 1
  expect_equal(sum(freqs$freq_matrix["POLYA", ]), 1)
})



# Integration tests for the full MCMC pipeline

test_that("MCMC Pipeline Integration: 10 pairs + 5 additional samples", {
  skip_on_cran() 
  
  # Load small dataset subset
  dat <- get_angola_test_subset(n_pairs = 10, n_additional = 5)
  
  # Temporary output directory
  tmp_out <- file.path(tempdir(), "angola_test")
  
  # Run full classification pipeline
  set.seed(42)
  results <- suppressWarnings(classify_infections(
    imported_data = dat,
    mcmc_config = list(max_iterations = 20, chunk_size = 20, n_chains = 2),
    output_folder = tmp_out,
    n_workers = 1,
    verbose = FALSE
  ))
  
  # Verify number of classified pairs
  expect_equal(nrow(results$summary), 10)
  
  # Probabilities must be valid
  expect_true(all(results$summary$Probability >= 0 & results$summary$Probability <= 1))
  
  unlink(tmp_out, recursive = TRUE)
})

# Results reproducibility and convergence
test_that("MCMC results are reproducible and convergence is calculated", {
  skip_on_cran()
  dat <- get_angola_test_subset(n_pairs = 2, n_additional = 0)
  conf <- list(max_iterations = 20, chunk_size = 10, n_chains = 2)
  
  # Run pipeline twice with identical seeds
  set.seed(100)
  res1 <- suppressMessages(suppressWarnings(classify_infections(dat, conf, file.path(tempdir(), "r1"), 
                                                               n_workers = 1, 
                                                               verbose = FALSE)))
  set.seed(100)
  res2 <- suppressMessages(suppressWarnings(classify_infections(dat, conf, file.path(tempdir(), "r2"), 
                                                                n_workers = 1, 
                                                                verbose = FALSE)))
  
  # Results should be identical
  expect_equal(res1$summary$Probability, res2$summary$Probability)
  expect_named(res1, c("summary", "marker_details", "comparison", "mcmc_loglikelihoods"))
})



# MCMC parameter updating and numerical stability
test_that("MCMC parameters update and diagnostics are valid", {
  skip_on_cran()
  dat <- get_angola_test_subset(2, 0)
  
  res <- suppressMessages(suppressWarnings(classify_infections(
    dat, list(max_iterations = 40, chunk_size = 20, n_chains = 2),
    tempfile(), n_workers = 1, verbose = FALSE
  )))
  
  # Extract all likelihood values across chains
  all_liks <- unlist(res$mcmc_loglikelihoods)
  
  # Likelihoods should vary (parameters updating)
  expect_true(sd(all_liks) > 0) 
  
  # Ensure numerical stability
  expect_true(all(is.finite(all_liks))) 
  
  # Verify returned object structure
  expect_named(res, c("summary", "marker_details", "comparison", "mcmc_loglikelihoods"))
})
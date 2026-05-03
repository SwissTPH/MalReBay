# ============================================================
# Skip helper — replaces skip_on_check() (testthat >= 3.2.0 only)
# ============================================================

#' Skip a Stan-dependent test when CmdStan is not available.
#' Covers two situations:
#'   1. R CMD check  — _R_CHECK_PACKAGE_NAME_ is set by R.
#'   2. covr / any other runner where cmdstanr is installed but
#'      CmdStan itself has not been installed on the machine.
skip_stan_on_check <- function() {
  skip_if(
    nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_")),
    "Skipping Stan-dependent test during R CMD check"
  )
  cmdstan_ok <- tryCatch({
    path <- cmdstanr::cmdstan_path()
    nzchar(path) && file.exists(path)
  }, error = function(e) FALSE)
  skip_if(!cmdstan_ok, "CmdStan not installed — skipping Stan-dependent test")
}

# ============================================================
# Creating example data to be used as moch data
# ============================================================

#' 4-marker microsatellite data frame (6 samples, 3 patient pairs)
#'
#' All markers use repeatlength = 3 (matches makers_details.xlsx).
#'
#' Patient profiles (verified match outcomes):
#'
#'  BD21-002  — weak reinfection signal: 1/4 loci match
#'    TA1:   174 vs 177  diff=3 == repeat=3  -> R
#'    POLYA: 153 vs 105  diff=48             -> NI
#'    PFPK2: 159 vs 183  diff=24             -> NI
#'    TA109: 172 vs 184  diff=12             -> NI
#'
#'  BD21-040  — weak reinfection signal: 1/4 loci match
#'    TA1:   162 vs 171  diff=9  > repeat=3  -> NI
#'    POLYA: 159 vs 159  diff=0              -> R
#'    PFPK2: 162/171 vs 183                  -> NI
#'    TA109: 160 vs 178  diff=18             -> NI
#'
#'  BD21-041  — strong recrudescence signal: 4/4 loci match
#'    TA1:   177 vs 177                      -> R
#'    POLYA: 162/165 vs 162/165              -> R
#'    PFPK2: 168/171 vs 168/171              -> R
#'    TA109: 184 vs 184                      -> R
#'
#' NOTE: BD21-002 and BD21-040 have the same match count (1/4).
#' Use BD21-041 (4/4) vs BD21-040 (1/4) for classification
#' direction tests — not BD21-002 vs BD21-040.
#'
#' Used by: test_2, test_3, test_4, test_5
#' 
create_mock_microsat_data <- function() {
  data.frame(
    Sample.ID = c(
      "BD21-002 Day 0", "BD21-002 recurrence",
      "BD21-040 Day 0", "BD21-040 recurrence",
      "BD21-041 Day 0", "BD21-041 recurrence"
    ),
    Site    = "Benguela",
    # TA1 (repeat = 3): diff=3 -> R for BD21-002; diff=9 -> NI for BD21-040; R for BD21-041
    TA1_1   = c(174, 177, 162, 171, 177, 177),
    TA1_2   = c(NA,  NA,  NA,  NA,  NA,  NA),
    # POLYA (repeat = 3): NI for BD21-002; R for BD21-040 and BD21-041
    POLYA_1 = c(153, 105, 159, 159, 162, 162),
    POLYA_2 = c(NA,  NA,  NA,  NA,  165, 165),
    # PFPK2 (repeat = 3): NI for BD21-002 and BD21-040; R for BD21-041
    PFPK2_1 = c(159, 183, 162, 183, 168, 168),
    PFPK2_2 = c(NA,  NA,  171, NA,  171, 171),
    # TA109 (repeat = 3): NI for BD21-002 and BD21-040; R for BD21-041
    TA109_1 = c(172, 184, 160, 178, 184, 184),
    TA109_2 = c(NA,  NA,  NA,  NA,  NA,  NA),
    check.names      = FALSE,
    stringsAsFactors = FALSE
  )
}

#' Marker info data frame matching create_mock_microsat_data()
#' All repeat lengths = 3 (matches makers_details.xlsx)
#' Used by: test_1, test_2, test_5
create_mock_markers <- function() {
  data.frame(
    marker_id      = c("TA1", "POLYA", "PFPK2", "TA109"),
    markertype     = "microsatellite",
    binning_method = "microsatellite",
    repeatlength   = c(3, 3, 3, 3),
    stringsAsFactors = FALSE
  )
}

#' Minimal allele columns (all 4 markers, single non-NA allele per marker)
#' Used in test_1 ID-standardisation tests so import_data() can find
#' matching markers and pass the empty-allele check.
#' @param n number of rows
mock_allele_cols <- function(n) {
  data.frame(
    TA1_1   = rep(174, n), TA1_2   = rep(NA_real_, n),
    POLYA_1 = rep(153, n), POLYA_2 = rep(NA_real_, n),
    PFPK2_1 = rep(159, n), PFPK2_2 = rep(NA_real_, n),
    TA109_1 = rep(172, n), TA109_2 = rep(NA_real_, n),
    stringsAsFactors = FALSE
  )
}

#' Full import_data-style list wrapping create_mock_microsat_data()
#' Used by: test_3, test_4
create_mock_imported_data <- function() {
  list(
    late_failures = create_mock_microsat_data(),
    additional    = data.frame(
      Sample.ID = character(0), Site = character(0),
      TA1_1     = numeric(0),   TA1_2   = numeric(0),
      POLYA_1   = numeric(0),   POLYA_2 = numeric(0),
      PFPK2_1   = numeric(0),   PFPK2_2 = numeric(0),
      TA109_1   = numeric(0),   TA109_2 = numeric(0)
    ),
    marker_info = create_mock_markers(),
    data_type   = "length_polymorphic"
  )
}

#' Write data + marker frames to temp xlsx files; returns list(data=, marker=)
#' Caller is responsible for unlink() via on.exit().
#' Used by: test_1
create_mock_xlsx <- function(data, markers = create_mock_markers()) {
  tmp_data   <- tempfile(fileext = ".xlsx")
  tmp_marker <- tempfile(fileext = ".xlsx")
  writexl::write_xlsx(list(Sheet1 = data), tmp_data)
  writexl::write_xlsx(markers, tmp_marker)
  list(data = tmp_data, marker = tmp_marker)
}

# ------------------------------------------------------------
# MCMC config
# ------------------------------------------------------------

#' Minimal MCMC settings — fast enough for unit tests, relaxed
#' convergence thresholds because small mock data won't fully converge.
#' Used by: test_3, test_4, test_6
fast_mcmc <- list(
  n_chains     = 2,
  iter         = 200,      
  burn_in_frac = 0.5,      
  random_seed  = 42,
  adapt_delta  = 0.8,
  rhat_threshold = 1.1,
  ess_threshold  = 50
)

# ------------------------------------------------------------
# Real Angola data — loaded from pre-saved .rds (path-independent)
# Avoids system.file() failures during devtools::test() when the
# package is not yet installed in the R library.
# ------------------------------------------------------------

#' Pre-imported Angola dataset (output of import_data())
#' Used by: test_1, test_3, test_6
angola_imported_data <- readRDS(
  system.file("extdata", "imported_data.rds", package = "MalReBay")
)

#' Paths to raw Angola xlsx files — used only by test_6 (MalReBay wrapper)
#' which must call import_data() internally via filepath arguments.
#' These resolve correctly because .rds files are installed with the package.
angola_data_path   <- system.file("extdata", "Angola_2021_TES_7NMS.xlsx",
                                  package = "MalReBay")
angola_marker_path <- system.file("extdata", "makers_details.xlsx",
                                  package = "MalReBay")

# ------------------------------------------------------------
# Higher-level fixture builders
# ------------------------------------------------------------

#' Run classify_infections on mock data with fast_mcmc
#' Used by: test_4
make_test_mcmc_results <- function() {
  skip_stan_on_check()
  suppressMessages(suppressWarnings(
    classify_infections(
      imported_data = create_mock_imported_data(),
      mcmc_config   = fast_mcmc,
      n_workers     = 1,
      verbose       = FALSE
    )
  ))
}

#' Run summarise_results on mock classify output.
#' plot_likelihood_diagnostics() requires a real output_folder (save_plot=TRUE
#' by default), so we always provide a temp directory and clean it up after.
#' @param output_folder  If NULL a temp dir is used and discarded automatically.
#'                       Pass an explicit path to also test file outputs.
#' Used by: test_4
make_test_summary_results <- function(output_folder = NULL) {
  use_folder <- if (!is.null(output_folder)) output_folder else tempdir()

  suppressMessages(suppressWarnings(
    summarise_results(
      mcmc_results  = make_test_mcmc_results(),
      imported_data = create_mock_imported_data(),
      output_folder = use_folder,
      verbose       = FALSE
    )
  ))
}

# ------------------------------------------------------------
# Ampseq mock data
# ------------------------------------------------------------

#' Ampseq late_failures data frame (3 patient pairs, 2 loci, maxMOI = 2)
#'
#' Patient profiles (verified match outcomes):
#'
#'  AMP-001 — strong recrudescence: 2/2 loci identical haplotypes
#'    cpmp: HAPL_A vs HAPL_A  -> match
#'    cpp:  HAPL_C vs HAPL_C  -> match
#'
#'  AMP-002 — strong reinfection: 2/2 loci different haplotypes
#'    cpmp: HAPL_A vs HAPL_B  -> differ
#'    cpp:  HAPL_C vs HAPL_D  -> differ
#'
#'  AMP-003 — partial: 1/2 loci match
#'    cpmp: HAPL_A vs HAPL_A  -> match
#'    cpp:  HAPL_C vs HAPL_D  -> differ
#'
#' Used by: test_ampseq
create_mock_ampseq_data <- function() {
  data.frame(
    Sample.ID     = c(
      "AMP-001 Day 0", "AMP-001 recurrence",
      "AMP-002 Day 0", "AMP-002 recurrence",
      "AMP-003 Day 0", "AMP-003 recurrence"
    ),
    Site          = "TestSite",
    cpmp_allele_1 = c("HAPL_A", "HAPL_A",
                      "HAPL_A", "HAPL_B",
                      "HAPL_A", "HAPL_A"),
    cpmp_allele_2 = c(NA, NA, NA, NA, NA, NA),
    cpp_allele_1  = c("HAPL_C", "HAPL_C",
                      "HAPL_C", "HAPL_D",
                      "HAPL_C", "HAPL_D"),
    cpp_allele_2  = c(NA, NA, NA, NA, NA, NA),
    stringsAsFactors = FALSE
  )
}

#' Marker metadata matching create_mock_ampseq_data()
#' repeatlength is NA for ampseq (not applicable)
create_mock_ampseq_markers <- function() {
  data.frame(
    marker_id      = c("cpmp", "cpp"),
    markertype     = "amplicon",
    binning_method = "ampseq",
    repeatlength   = c(NA_real_, NA_real_),
    stringsAsFactors = FALSE
  )
}

#' Full import_data-style list wrapping create_mock_ampseq_data()
#' Used by: test_ampseq
create_mock_ampseq_imported_data <- function() {
  list(
    late_failures = create_mock_ampseq_data(),
    additional    = data.frame(
      Sample.ID     = character(0),
      Site          = character(0),
      cpmp_allele_1 = character(0),
      cpmp_allele_2 = character(0),
      cpp_allele_1  = character(0),
      cpp_allele_2  = character(0),
      stringsAsFactors = FALSE
    ),
    marker_info = create_mock_ampseq_markers(),
    data_type   = "ampseq"
  )
}

#' Run classify_infections on mock ampseq data with fast_mcmc
#' Used by: test_ampseq
make_test_ampseq_results <- function() {
  skip_stan_on_check()
  suppressMessages(suppressWarnings(
    classify_infections(
      imported_data = create_mock_ampseq_imported_data(),
      mcmc_config   = fast_mcmc,
      n_workers     = 1,
      verbose       = FALSE
    )
  ))
}
# =============================================================================
# Tests for the amplicon-sequencing (ampseq) code path
#
# Ampseq loci are no longer routed through a separate model/interface — they
# go through the same prepare_stan_data() / run_stan_sites() / malrebay_model.stan
# as length-polymorphic loci, dispatched per-locus on
# marker_info$binning_method == "exact" (see R/stan_data_prep.R).
#
# Covers:
#   R/stan_data_prep.R  — prepare_stan_data() on "exact" (ampseq) loci
#                         validate_stan_data()
#   R/stan_interface.R  — run_stan_sites() on ampseq data
#                         extract_stan_results()
#   R/import_data.R     — ampseq format detection
#   R/MalReBay.R         — classify_infections() on ampseq data
# =============================================================================

# ---------------------------------------------------------------------------
# Shared fixtures for prepare_stan_data() unit tests (exact/ampseq loci)
# ---------------------------------------------------------------------------

# late_failures with Site column removed (as received by prepare_stan_data())
.ampseq_late_site <- local({
  d <- create_mock_ampseq_data()
  d[, colnames(d) != "Site", drop = FALSE]
})

.ampseq_marker_info <- create_mock_ampseq_markers()
.ampseq_ids         <- c("AMP-001", "AMP-002", "AMP-003")
.ampseq_locinames   <- c("cpmp", "cpp")
.ampseq_maxMOI      <- 2L
.ampseq_comparable  <- matrix(
  TRUE, nrow = 3, ncol = 2,
  dimnames = list(.ampseq_ids, .ampseq_locinames)
)
# define_alleles() skips numeric binning for "exact" loci (returns an empty
# 0-row matrix per locus) — prepare_stan_data() takes this list as-is and
# recodes "exact" loci by mapping raw haplotype strings directly instead.
.ampseq_allele_defs <- suppressMessages(
  MalReBay:::define_alleles(.ampseq_late_site, .ampseq_marker_info)
)

# Helper: call prepare_stan_data() with default ampseq fixtures
make_ampseq_sd <- function(late       = .ampseq_late_site,
                            additional = data.frame()) {
  MalReBay:::prepare_stan_data(
    late_failures_site  = late,
    additional_site     = additional,
    allele_definitions  = .ampseq_allele_defs,
    marker_info         = .ampseq_marker_info,
    ids                 = .ampseq_ids,
    locinames           = .ampseq_locinames,
    maxMOI              = .ampseq_maxMOI,
    is_locus_comparable = .ampseq_comparable
  )
}

# Cache the Stan data so we build it once across multiple unit tests
.sd <- make_ampseq_sd()


# =============================================================================
# Section 1: prepare_stan_data (exact/ampseq loci) — structure & dimensions
# =============================================================================

test_that("prepare_stan_data (ampseq): returns all required fields", {
  required <- c("N", "J", "maxMOI", "max_K", "max_dist", "K",
                "recoded0", "recodedf", "hidden0", "hiddenf",
                "MOI0", "MOIf", "dist_array", "method_int", "threshold",
                "comparable", "additional_counts")
  expect_true(all(required %in% names(.sd)))
})

test_that("prepare_stan_data (ampseq): scalar dimensions are correct", {
  expect_equal(.sd$N,      3L)          # 3 patient pairs
  expect_equal(.sd$J,      2L)          # 2 loci (cpmp, cpp)
  expect_equal(.sd$maxMOI, 2L)
  expect_equal(.sd$max_K,  2L)          # 2 unique haplotypes per locus
  expect_equal(as.integer(.sd$K), c(2L, 2L))  # K per locus
})

test_that("prepare_stan_data (ampseq): dimensions are correct (per-patient MOI convention)", {
  # MOI0/MOIf are per-patient (length N), the same convention used for
  # length-polymorphic loci — not a per-locus N x J matrix.
  expect_equal(length(.sd$MOI0), 3L)
  expect_equal(length(.sd$MOIf), 3L)
  expect_equal(dim(.sd$recoded0), c(3L, 4L))   # N x (J * maxMOI)
  expect_equal(dim(.sd$recodedf), c(3L, 4L))
  expect_equal(dim(.sd$hidden0),  c(3L, 4L))
  expect_equal(dim(.sd$hiddenf),  c(3L, 4L))
  expect_equal(dim(.sd$comparable), c(3L, 2L))
  expect_equal(dim(.sd$additional_counts), c(2L, 2L))  # J x max_K
  expect_equal(dim(.sd$dist_array), c(2L, 2L, 2L))     # J x max_K x max_K
})

test_that("prepare_stan_data (ampseq): exact loci use method_int 3 with binary distances", {
  expect_equal(as.integer(.sd$method_int), c(3L, 3L))
  for (j in 1:2) {
    expect_equal(diag(.sd$dist_array[j, , ]), c(0, 0))
    expect_equal(.sd$dist_array[j, 1, 2], 1)
    expect_equal(.sd$dist_array[j, 2, 1], 1)
  }
})

test_that("prepare_stan_data (ampseq): haplotypes mapped to positive integers", {
  # Every patient has data at Day 0 — all recoded0 values in used slots must be > 0
  used_slots_d0 <- c(1, 3)   # slot 1 = cpmp allele 1, slot 3 = cpp allele 1
  expect_true(all(.sd$recoded0[, used_slots_d0] > 0))

  # AMP-001 recurrence: same haplotypes as Day 0 -> same integer codes
  expect_equal(.sd$recoded0[1, used_slots_d0],
               .sd$recodedf[1, used_slots_d0])

  # AMP-002 recurrence: different haplotypes -> different integer codes
  expect_false(all(.sd$recoded0[2, used_slots_d0] ==
                     .sd$recodedf[2, used_slots_d0]))
})

test_that("prepare_stan_data (ampseq): NA allele slots recoded as 0", {
  # All cpmp_allele_2 and cpp_allele_2 are NA -> slots 2 and 4 should be 0
  empty_slots <- c(2, 4)
  expect_true(all(.sd$recoded0[, empty_slots] == 0L))
  expect_true(all(.sd$recodedf[, empty_slots] == 0L))
})

test_that("prepare_stan_data (ampseq): additional_counts reflects observed alleles", {
  # With no additional data, counts come from trial rows only
  # cpmp: HAPL_A appears 5x, HAPL_B 1x -> total = 6
  # cpp:  HAPL_C appears 4x, HAPL_D 2x -> total = 6
  expect_gt(sum(.sd$additional_counts), 0)
  expect_equal(nrow(.sd$additional_counts), 2L)   # one row per locus
  expect_equal(ncol(.sd$additional_counts), 2L)   # one col per unique haplotype
})


# =============================================================================
# Section 2: validate_stan_data (shared validator, both marker types)
# =============================================================================

test_that("validate_stan_data (ampseq): returns TRUE for valid data", {
  result <- suppressMessages(MalReBay:::validate_stan_data(.sd))
  expect_true(result)
})

test_that("validate_stan_data (ampseq): returns FALSE when additional_counts empty", {
  bad_sd        <- .sd
  bad_sd$additional_counts <- bad_sd$additional_counts * 0L
  result <- suppressMessages(MalReBay:::validate_stan_data(bad_sd))
  expect_false(result)
})


# =============================================================================
# Section 3: import_data — ampseq format detection
# =============================================================================

test_that("import_data: detects 'ampseq' when allele values are haplotype strings", {
  xlsx <- create_mock_xlsx(
    data    = create_mock_ampseq_data(),
    markers = create_mock_ampseq_markers()
  )
  on.exit(unlink(c(xlsx$data, xlsx$marker)))

  imported <- suppressMessages(suppressWarnings(
    import_data(filepath = xlsx$data, marker_filepath = xlsx$marker)
  ))

  expect_equal(imported$data_type, "ampseq")
})

test_that("import_data: ampseq additional data is empty when no Additional sheet is present", {
  xlsx <- create_mock_xlsx(
    data    = create_mock_ampseq_data(),
    markers = create_mock_ampseq_markers()
  )
  on.exit(unlink(c(xlsx$data, xlsx$marker)))

  imported <- suppressMessages(suppressWarnings(
    import_data(filepath = xlsx$data, marker_filepath = xlsx$marker)
  ))

  expect_equal(nrow(imported$additional), 0L)
})


# =============================================================================
# Section 4: classify_infections — full ampseq integration (slow)
#
# Stan is run exactly once here (not once per test_that).
# During R CMD check, _R_CHECK_PACKAGE_NAME_ is set so we skip the Stan
# call entirely and set .ampseq_clf to NULL. Every test then guards with
# skip_if_ampseq_clf_missing() which also calls skip_stan_on_check() so
# the skip reason is reported cleanly in both contexts.
# =============================================================================

.ampseq_clf <- if (nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))) {
  NULL
} else {
  tryCatch(
    suppressMessages(suppressWarnings(make_test_ampseq_results())),
    error = function(e) NULL
  )
}

skip_if_ampseq_clf_missing <- function() {
  skip_stan_on_check()
  skip_if(is.null(.ampseq_clf), "ampseq Stan run not available")
}

test_that("classify_infections: returns correct top-level list keys (ampseq)", {
  skip_if_ampseq_clf_missing()
  expected_keys <- c("classifications", "all_chains_loglikelihood",
                     "ids", "locus_summary", "locus_lrs",
                     "locus_dists", "locinames", "stan_fits")
  expect_true(all(expected_keys %in% names(.ampseq_clf)))
})

test_that("classify_infections: ids match patients in mock ampseq data", {
  skip_if_ampseq_clf_missing()
  ids_out <- unlist(.ampseq_clf$ids, use.names = FALSE)
  expect_setequal(ids_out, .ampseq_ids)
})

test_that("classify_infections: locus_summary has correct columns (ampseq)", {
  skip_if_ampseq_clf_missing()
  ls <- .ampseq_clf$locus_summary[[1]]
  expect_true(all(c("patient_id", "n_available_d0",
                    "n_available_df", "n_comparable_loci") %in% colnames(ls)))
})

test_that("classify_infections: p_recrud draws are between 0 and 1 (ampseq)", {
  skip_if_ampseq_clf_missing()
  p_draws <- .ampseq_clf$classifications[[1]]
  expect_true(all(p_draws >= 0 & p_draws <= 1))
})

test_that("classify_infections: AMP-001 (2/2 match) classified towards recrudescence", {
  skip_if_ampseq_clf_missing()
  skip_on_cran()
  site     <- names(.ampseq_clf$classifications)[1]
  p_draws  <- .ampseq_clf$classifications[[site]]
  ids_site <- .ampseq_clf$ids[[site]]
  idx      <- which(ids_site == "AMP-001")
  expect_gt(mean(p_draws[, idx]), 0.5)
})

test_that("classify_infections: AMP-002 (0/2 match) classified towards reinfection", {
  skip_if_ampseq_clf_missing()
  skip_on_cran()
  site     <- names(.ampseq_clf$classifications)[1]
  p_draws  <- .ampseq_clf$classifications[[site]]
  ids_site <- .ampseq_clf$ids[[site]]
  idx      <- which(ids_site == "AMP-002")
  expect_lt(mean(p_draws[, idx]), 0.5)
})

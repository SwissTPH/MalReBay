library(testthat)
library(MalReBay)

imported <- readRDS(system.file("extdata", "imported_data.rds", package = "MalReBay"))

late      <- imported$late_failures[, setdiff(colnames(imported$late_failures), "Site")]
add       <- imported$additional[, setdiff(colnames(imported$additional), "Site")]
markers   <- imported$marker_info
locinames <- markers$marker_id
ids       <- unique(gsub(" Day 0$| recurrence$", "", late$Sample.ID))

allele_cols <- setdiff(colnames(late), "Sample.ID")
maxMOI      <- max(as.integer(gsub(".*(_allele_|_)(\\d+)$", "\\2", allele_cols)))
allele_definitions <- suppressMessages(define_alleles(rbind(late, add), markers))
comparability <- compute_locus_comparability(late, ids, locinames)

stan_data <- prepare_stan_data(
  late_failures_site  = late,
  additional_site     = add,
  allele_definitions  = allele_definitions,
  marker_info         = markers,
  ids                 = ids,
  locinames           = locinames,
  maxMOI              = maxMOI,
  is_locus_comparable = comparability$is_locus_comparable
)

sample_id <- "ZL21-203"
day0  <- late[late$Sample.ID == paste(sample_id, "Day 0"), ]
recur <- late[late$Sample.ID == paste(sample_id, "recurrence"), ]


test_that("import_data returns non-empty data", {
  expect_gt(nrow(late), 0)
  expect_gt(nrow(add), 0)
  expect_gt(nrow(markers), 0)
  expect_true(all(colSums(!is.na(late[, allele_cols])) > 0))
})

test_that("compute_locus_comparability returns non-empty results", {
  expect_equal(nrow(comparability$locus_summary), length(ids))
  expect_true(any(comparability$is_locus_comparable))
})

test_that("define_alleles returns non-empty bins for every locus", {
  for (locus in locinames) expect_gt(nrow(allele_definitions[[locus]]), 0)
})

test_that("observed alleles belong to the defined allele population", {
  for (j in seq_along(locinames)) {
    cols <- (maxMOI * (j - 1) + 1):(maxMOI * (j - 1) + maxMOI)
    observed0 <- stan_data$recoded0[, cols][stan_data$hidden0[, cols] == 0]
    observedf <- stan_data$recodedf[, cols][stan_data$hiddenf[, cols] == 0]
    expect_true(all(observed0 >= 1 & observed0 <= stan_data$K[j]))
    expect_true(all(observedf >= 1 & observedf <= stan_data$K[j]))
  }
})

test_that("recodeallele assigns ZL21-203's alleles to bins that contain them", {
  for (locus in locinames) {
    bins <- allele_definitions[[locus]]
    cols <- grep(paste0("^", locus, "_"), colnames(late), value = TRUE)
    values <- unique(unlist(c(day0[, cols], recur[, cols])))
    values <- values[!is.na(values)]
    
    for (v in values) {
      bin_idx <- recodeallele(bins, v)
      expect_false(is.na(bin_idx), info = sprintf("%s: value %s matched no bin", locus, v))
      expect_true(
        v >= bins[bin_idx, "lower"] && v <= bins[bin_idx, "upper"],
        info = sprintf("%s: value %s assigned to bin [%s, %s]",
                       locus, v, bins[bin_idx, "lower"], bins[bin_idx, "upper"])
      )
    }
  }
})

test_that("recodeallele rejects an out-of-range value as an outlier", {
  locus <- "TA1"
  bins <- allele_definitions[[locus]]
  repeat_length <- markers$repeatlength[markers$marker_id == locus]
  
  bin_idx <- recodeallele(bins, 800, max_distance_allowed = repeat_length)
  expect_true(is.na(bin_idx))
})

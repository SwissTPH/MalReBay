library(testthat)
library(MalReBay)

zaire_data <- readRDS(system.file("extdata", "imported_data.rds", package = "MalReBay"))
pooled <- rbind(zaire_data$late_failures, zaire_data$additional)
allele_definitions <- suppressMessages(define_alleles(pooled, zaire_data$marker_info))

sample_id <- "ZL21-203"
locinames <- zaire_data$marker_info$marker_id
day0  <- zaire_data$late_failures[zaire_data$late_failures$Sample.ID == paste(sample_id, "Day 0"), ]
recur <- zaire_data$late_failures[zaire_data$late_failures$Sample.ID == paste(sample_id, "recurrence"), ]

test_that("recodeallele assigns ZL21-203's alleles to bins that contain them", {
  for (locus in locinames) {
    bins <- allele_definitions[[locus]]
    cols <- grep(paste0("^", locus, "_"), colnames(zaire_data$late_failures), value = TRUE)
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
  repeat_length <- zaire_data$marker_info$repeatlength[zaire_data$marker_info$marker_id == locus]
  
  bin_idx <- recodeallele(bins, 800, max_distance_allowed = repeat_length)
  expect_true(is.na(bin_idx))
})

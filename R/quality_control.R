#' Compute locus comparability per patient
#'
#' @description For each patient, determines which loci have non-missing
#'   data at both Day 0 and recurrence timepoints ("comparable"), and tallies
#'   per-patient availability. A locus is comparable only if at least one
#'   non-NA allele is present at that locus for both timepoints.
#'
#' @param late_site The late failures data for a single site (Site column
#'   already removed), with Sample.ID and per-locus allele columns.
#' @param ids Character vector of patient IDs (without " Day 0"/" recurrence"
#'   suffixes) to check.
#' @param locinames Character vector of locus names to check.
#'
#' @return A list with two elements:
#'   \item{locus_summary}{A data frame with one row per patient: patient_id,
#'     n_available_d0, n_available_df, n_comparable_loci.}
#'   \item{is_locus_comparable}{A logical matrix (patients x loci) indicating
#'     which patient-locus pairs are comparable.}
#'
#' @export
compute_locus_comparability <- function(late_site, ids, locinames) {
  locus_summary <- data.frame(patient_id = ids,
                              n_available_d0 = 0L,
                              n_available_df = 0L,
                              n_comparable_loci = 0L)
  
  is_locus_comparable <- matrix(FALSE,
                                nrow = length(ids),
                                ncol = length(locinames),
                                dimnames = list(ids, locinames))
  
  locus_cols <- setNames(
    lapply(locinames, function(ln) 
      grep(paste0("^", ln, "_"), colnames(late_site), value = TRUE)),
    locinames
  )
  
  for (i in seq_along(ids)) {
    pid <- ids[i]
    d0_row <- late_site[grepl(paste0("\\b", pid, " Day 0\\b"), late_site$Sample.ID), ]
    df_row <- late_site[grepl(paste0("\\b", pid, " recurrence\\b"), late_site$Sample.ID), ]
    if (nrow(d0_row) == 0 || nrow(df_row) == 0) next
    for (ln in locinames) {
      lc <- locus_cols[[ln]]
      if (any(!is.na(d0_row[, lc]))) locus_summary$n_available_d0[i] <- locus_summary$n_available_d0[i] + 1L
      if (any(!is.na(df_row[, lc]))) locus_summary$n_available_df[i] <- locus_summary$n_available_df[i] + 1L
      if (any(!is.na(d0_row[, lc])) && any(!is.na(df_row[, lc]))) {
        locus_summary$n_comparable_loci[i] <- locus_summary$n_comparable_loci[i] + 1L
        is_locus_comparable[pid, ln] <- TRUE
      }
    }
  }
  
  list(
    locus_summary = locus_summary,
    is_locus_comparable = is_locus_comparable
  )
}

#' Check dataset viability per site
#'
#' @description For each site, checks whether there are enough paired
#'   samples (Day 0 + recurrence) with at least one comparable locus to be
#'   viable for analysis.
#'
#' @param imported_data Output list from import_data().
#' @param min_paired_samples Minimum paired samples (each with at least
#'   min_comparable_loci comparable loci) required per site. Default 10.
#' @param min_comparable_loci Minimum comparable loci a paired sample must
#'   have to count toward min_paired_samples. Default 1.
#' @param verbose Logical. Warn for sites that fail the check.
#' @return A data frame: Site, n_paired, n_viable_pairs, viable.
#' @examples
#' \dontrun{
#'   data_file <- system.file("extdata", 
#'                            "Angola_2021_TES_7NMS.xlsx", 
#'                            package = "MalReBay")
#'   marker_file <- system.file("extdata", 
#'                              "makers_details.xlsx", 
#'                              package = "MalReBay")
#'   imported <- import_data(filepath = data_file, 
#'                           marker_filepath = marker_file)
#'   quality  <- data_quality_check(imported)
#' }
#' @export
data_quality_check <- function(imported_data,
                               min_paired_samples = 10,
                               min_comparable_loci = 1,
                               verbose = TRUE) {
  late_failures <- imported_data$late_failures
  locinames     <- imported_data$marker_info$marker_id
  sites         <- unique(late_failures$Site)
  
  site_results <- lapply(sites, function(site) {
    late_site <- late_failures[late_failures$Site == site, ]
    late_site <- late_site[, colnames(late_site) != "Site", drop = FALSE]
    
    day0_ids   <- unique(gsub(" Day 0", "", late_site$Sample.ID[grepl("Day 0", late_site$Sample.ID)]))
    recur_ids  <- unique(gsub(" recurrence", "", late_site$Sample.ID[grepl("recurrence", late_site$Sample.ID)]))
    paired_ids <- intersect(day0_ids, recur_ids)
    n_paired   <- length(paired_ids)
    
    if (n_paired == 0) {
      n_viable_pairs <- 0L
    } else {
      comparability  <- compute_locus_comparability(late_site, paired_ids, locinames)
      n_viable_pairs <- sum(comparability$locus_summary$n_comparable_loci >= min_comparable_loci)
    }
    
    viable <- n_paired >= min_paired_samples && n_viable_pairs >= min_paired_samples
    
    if (!viable && verbose) {
      warning("WARNING: Site '", site, "' is not viable -- ", n_paired,
              " paired sample(s), ", n_viable_pairs,
              " with >= ", min_comparable_loci, " comparable locus/loci ",
              "(need >= ", min_paired_samples, " of each).", call. = FALSE)
    }
    
    data.frame(Site = site, n_paired = n_paired, n_viable_pairs = n_viable_pairs,
               viable = viable, stringsAsFactors = FALSE)
  })
  
  do.call(rbind, site_results)
}

#' Identify which locinames belong to the MSP1/MSP2 families
#'
#' @description Cross-checks locinames against the standard MSP1 (K1, MAD20,
#'   RO33) and MSP2 (3D7, FC27, IC) allelic family names to identify which
#'   markers in a dataset belong to which family. Digit "0" and letter "O"
#'   are treated as equivalent (e.g. "RO33" vs "R033"), matching ignores
#'   case/whitespace.
#'
#' @param locinames Character vector of marker/locus names present in the dataset.
#' @return A list with msp1 and msp2 elements: the subset of locinames
#'   (original spelling, not normalized) matching each family's known variants.
#' @export
detect_msp_variants <- function(locinames) {
  msp1_ref <- c("K1", "MAD20", "RO33")
  msp2_ref <- c("3D7", "FC27", "IC")
  
  normalize <- function(x) gsub("0", "O", toupper(trimws(x)))
  normalized_locinames <- normalize(locinames)
  
  list(
    msp1 = locinames[normalized_locinames %in% normalize(msp1_ref)],
    msp2 = locinames[normalized_locinames %in% normalize(msp2_ref)]
  )
}

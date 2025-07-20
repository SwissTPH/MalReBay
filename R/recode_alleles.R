
recodeallele <- function(alleles_definitions_subset, proposed) {
  if (is.na(proposed) || nrow(alleles_definitions_subset) == 0) {
    return(NA_integer_)
  }
  # Find which bin 'proposed' falls into
  ret <- which(proposed > alleles_definitions_subset[, 1] & proposed <= alleles_definitions_subset[, 2])
  if (length(ret) == 0) {
    ret <- NA_integer_
  }
  return(ret[1]) 
}

recode_alleles <- function(genotypedata, alleles_definitions, marker_info) {
  ids <- unique(gsub(" Day 0", "", genotypedata$Sample.ID[grepl(" Day 0", genotypedata$Sample.ID)]))
  locinames <- names(alleles_definitions)
  nids <- length(ids)
  nloci <- length(locinames)
  
  MOI0 <- rep(0, nids)
  MOIf <- rep(0, nids)
  for (i in 1:nids) {
    for (j in 1:nloci) {
      locicolumns <- grepl(paste0("^", locinames[j], "_"), colnames(genotypedata))
      nalleles0 <- sum(!is.na(genotypedata[grepl(paste(ids[i], "Day 0"), genotypedata$Sample.ID), locicolumns]))
      nallelesf <- sum(!is.na(genotypedata[grepl(paste(ids[i], "Day Failure"), genotypedata$Sample.ID), locicolumns]))
      MOI0[i] <- max(MOI0[i], nalleles0)
      MOIf[i] <- max(MOIf[i], nallelesf)
    }
  }
  
  observeddatamatrix <- list()
  for (locus_name in locinames) {
    
    # Get the binning method for this specific locus
    locus_marker_info <- marker_info[marker_info$marker_id == locus_name, ]
    binning_method <- locus_marker_info$binning_method[1]
    
    # Get the raw allele data for this locus
    locus_cols_pattern <- paste0("^", locus_name, "_")
    locus_cols_logical <- grepl(locus_cols_pattern, colnames(genotypedata))
    raw_alleles_matrix <- as.matrix(genotypedata[, locus_cols_logical, drop = FALSE])
    
    recoded_matrix <- NULL
    if (binning_method %in% c("microsatellite", "cluster")) {
      current_definitions <- alleles_definitions[[locus_name]]
      recoded_matrix <- apply(raw_alleles_matrix, MARGIN = c(1, 2), FUN = function(cell_value) {
        recodeallele(alleles_definitions_subset = current_definitions, proposed = cell_value)
      })
      
    } else if (binning_method == "exact") {
      recoded_matrix <- raw_alleles_matrix
      
    } else {
      stop("Unknown binning_method '", binning_method, "' for locus: ", locus_name)
    }
    
    if (is.null(dim(recoded_matrix))) {
      recoded_matrix <- matrix(recoded_matrix, nrow = nrow(raw_alleles_matrix))
    }
    
    tempobservedata <- character(nids)
    for (i_patient_idx in 1:nids) {
      patient_id_str <- ids[i_patient_idx]
      day0_logical_idx <- grepl(paste(patient_id_str, "Day 0"), genotypedata$Sample.ID)
      dayf_logical_idx <- grepl(paste(patient_id_str, "Day Failure"), genotypedata$Sample.ID)
      
      day0_alleles <- recoded_matrix[day0_logical_idx, , drop = FALSE]
      day0_unique_sorted <- sort(unique(day0_alleles[!is.na(day0_alleles)]))
      
      dayf_alleles <- recoded_matrix[dayf_logical_idx, , drop = FALSE]
      dayf_unique_sorted <- sort(unique(dayf_alleles[!is.na(dayf_alleles)]))
      
      day0_str <- paste(day0_unique_sorted, collapse = "-")
      dayf_str <- paste(dayf_unique_sorted, collapse = "-")
      tempobservedata[i_patient_idx] <- paste(day0_str, dayf_str, sep = "/")
    }
    
    observeddatamatrix[[locus_name]] <- tempobservedata
  }
  
  MOItemp <- cbind(MOI0, MOIf)
  rownames(MOItemp) <- ids
  list(observeddatamatrix = observeddatamatrix, MOI = MOItemp)
}
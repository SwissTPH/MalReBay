

# Recode allele function
recodeallele = function(alleles_definitions_subset,proposed) {
  
  ret = which(proposed > alleles_definitions_subset[,1] & proposed <= alleles_definitions_subset[,2])
  if (length(ret) == 0) {
    ret = NA
  }
  ret
}

recode_alleles = function(genotypedata, alleles_definitions) {
  
  ########### generate MOI for each sample
  
  ids = unique(unlist(strsplit(genotypedata$Sample.ID[grepl("Day 0",genotypedata$Sample.ID)]," Day 0")))
  locinames = unique(sapply(colnames(genotypedata)[-1],function(x) strsplit(x,"_")[[1]][1]))
  nids = length(ids)
  nloci = length(locinames)
  
  
  MOI0 = rep(0,nids)
  MOIf = rep(0,nids)
  
  # for each individual, cycle through each locus and count number of separate alleles
  
  for (i in 1:nids) {
    for (j in 1:nloci) {
      locicolumns = grepl(paste(locinames[j],"_",sep=""),colnames(genotypedata))
      nalleles0 = sum(!is.na(genotypedata[grepl(paste(ids[i],"Day 0"),genotypedata$Sample.ID),locicolumns]))
      nallelesf = sum(!is.na(genotypedata[grepl(paste(ids[i],"Day Failure"),genotypedata$Sample.ID),locicolumns]))
      
      MOI0[i] = max(MOI0[i],nalleles0)
      MOIf[i] = max(MOIf[i],nallelesf)
    }
  }
  
  
  observeddatamatrix = list()
  for (j in 1:nloci) { 
    locus = locinames[j]
    # Use more robust pattern matching for column names
    locus_allele_cols_pattern = paste0("^", locus, "_") # Match from the beginning of the name
    locicolumns_logical = grepl(locus_allele_cols_pattern, colnames(genotypedata))
    
    current_locus_raw_df = genotypedata[, locicolumns_logical, drop = FALSE]
    
    raw_allele_values_matrix_locus_j = matrix(numeric(0), 
                                              nrow = nrow(genotypedata), 
                                              ncol = ncol(current_locus_raw_df)) # Pre-allocate
    
    if (ncol(current_locus_raw_df) > 0) {
      for(col_idx in 1:ncol(current_locus_raw_df)){
        raw_allele_values_matrix_locus_j[,col_idx] = as.numeric(as.character(current_locus_raw_df[[col_idx]]))
      }
    }
    
    newalleles_recoded_matrix = matrix(NA_integer_, # Store recoded integers or NA
                                       nrow = nrow(raw_allele_values_matrix_locus_j), 
                                       ncol = ncol(raw_allele_values_matrix_locus_j))
    
    if (ncol(raw_allele_values_matrix_locus_j) > 0) { # Only proceed if there are allele columns
      for (k_col_idx in 1:ncol(raw_allele_values_matrix_locus_j)) {   # Iterate over allele columns (e.g., _1, _2, _3)
        for (i_row_idx in 1:nrow(raw_allele_values_matrix_locus_j)) { # Iterate over samples in genotypedata
          
          proposed_value = raw_allele_values_matrix_locus_j[i_row_idx, k_col_idx]
          
          if (!is.na(proposed_value)) {
            
            newalleles_recoded_matrix[i_row_idx, k_col_idx] = recodeallele(
              alleles_definitions_subset = alleles_definitions[[j]], # Pass the definitions for current locus j
              proposed = proposed_value
            )
          } else {
            newalleles_recoded_matrix[i_row_idx, k_col_idx] = NA # Keep NA as NA
          }
        }
      }
    }
    
    tempobservedata = character(nids) # Pre-allocate for efficiency
    
    for (i_patient_idx in 1:nids) { # Inner loop over unique patient IDs (1 to nids)
      patient_id_str = ids[i_patient_idx] # Get the actual patient ID string
      
      day0_logical_idx = genotypedata$Sample.ID == paste(patient_id_str, "Day 0")
      dayf_logical_idx = genotypedata$Sample.ID == paste(patient_id_str, "Day Failure")
      
      day0_recoded_alleles_for_patient_locus = newalleles_recoded_matrix[day0_logical_idx, , drop = FALSE]
      day0_unique_sorted = sort(unique(day0_recoded_alleles_for_patient_locus[!is.na(day0_recoded_alleles_for_patient_locus)]))
      
      dayf_recoded_alleles_for_patient_locus = newalleles_recoded_matrix[dayf_logical_idx, , drop = FALSE]
      dayf_unique_sorted = sort(unique(dayf_recoded_alleles_for_patient_locus[!is.na(dayf_recoded_alleles_for_patient_locus)]))
      
      day0_str = paste(day0_unique_sorted, collapse = "-")
      dayf_str = paste(dayf_unique_sorted, collapse = "-")
      
      if (length(day0_unique_sorted) == 0) day0_str = ""
      if (length(dayf_unique_sorted) == 0) dayf_str = ""
      
      tempobservedata[i_patient_idx] = paste(day0_str, dayf_str, sep = "/")
    }
    observeddatamatrix[[j]] = tempobservedata
  }
  
  
  
  
  MOItemp = cbind(MOI0,MOIf)
  rownames(MOItemp) = ids
  list(observeddatamatrix = observeddatamatrix, MOI = MOItemp)
}
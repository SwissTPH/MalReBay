run_all_sites <- function(
    genotypedata_latefailures,
    additional_genotypedata,
    marker_info_subset,
    mcmc_config,
    data_type,
    output_folder) {
  
  n_chains <- mcmc_config$n_chains
  rhat_threshold <- mcmc_config$rhat_threshold
  ess_threshold <- mcmc_config$ess_threshold
  chunk_size <- mcmc_config$chunk_size
  max_iterations <- mcmc_config$max_iterations
  burn_in_frac <- mcmc_config$burn_in_frac
  record_hidden_alleles <- mcmc_config$record_hidden_alleles
  site_names <- unique(genotypedata_latefailures$Site)
  local_sites_classification <- list()
  local_sites_parameters <- list()
  local_sites_alleles0 <- list()
  local_sites_allelesf <- list()
  local_sites_ids <- list()
  local_sites_locus_summary <- list()
  local_sites_locus_lrs <- list()
  local_sites_locus_dists <- list()
  local_sites_locinames <- list()
  
  for (site in site_names) {
    jobname <- site
      
    # Site-Specific Data Preparation (THE CORRECTED LOGIC)
    genotypedata_RR <- genotypedata_latefailures[genotypedata_latefailures$Site == site, ]
    additional_neutral <- additional_genotypedata[additional_genotypedata$Site == site, ]

    ids <- unique(gsub(" Day 0", "", genotypedata_RR$Sample.ID[grepl("Day 0", genotypedata_RR$Sample.ID)]))
    nids <- length(ids)
    if (nids == 0) {
        message("INFO: No valid Day 0/Failure pairs found for site '", site, "'. Skipping.")
        next
    }

    mcmc_engine_function <- NULL
    engine_args <- NULL
    
    if (data_type == "ampseq") {
      message("INFO: Preparing AMPSEQ MCMC engine for site: ", site)
      
      # For ampseq, we don't need 'alleles_definitions' or complex MOI calcs here
      marker_cols <- grep("_allele_\\d+$", colnames(genotypedata_RR), value = TRUE)
      maxMOI <- if (length(marker_cols) > 0) max(as.numeric(gsub(".*_allele_", "", marker_cols)), na.rm = TRUE) else 0
      locinames <- unique(gsub("_allele_\\d+$", "", marker_cols))
      nloci <- length(locinames)
      alleles_definitions_RR <- list()
      mcmc_engine_function <- run_one_chain_ampseq
      engine_args <- list(
        nids = nids, ids = ids, nloci = nloci, maxMOI = maxMOI, locinames = locinames,
        genotypedata_RR = genotypedata_RR[, -which(names(genotypedata_RR) == "Site")],
        additional_neutral = additional_neutral[, -which(names(additional_neutral) == "Site")],
        marker_info = marker_info_subset,
        record_hidden_alleles = record_hidden_alleles
      )
      
    } else { # This is the length_polymorphic path
      message("INFO: Preparing LENGTH-POLYMORPHIC MCMC engine for site: ", site)
      
      marker_cols <- grep("_\\d+$", colnames(genotypedata_RR), value = TRUE)
      maxMOI <- if (length(marker_cols) > 0) max(as.numeric(gsub(".*_", "", marker_cols)), na.rm = TRUE) else 0
      maxalleles <- 30 # Or some other default
      k <- rep(maxalleles, nrow(marker_info_subset))
      names(k) <- marker_info_subset$marker_id
      
      # Define alleles using our unified function
      alleles_definitions_RR <- define_alleles(
        genotypedata = rbind(genotypedata_RR, additional_neutral), 
        marker_info_subset = marker_info_subset, 
        maxk = k
      )
      
      locinames <- names(alleles_definitions_RR)
      local_sites_locinames[[site]] <- locinames
      nloci <- length(locinames)
      
      mcmc_engine_function <- run_one_chain # Your original, updated function
      engine_args <- list(
        nids = nids, ids = ids, nloci = nloci, maxMOI = maxMOI, locinames = locinames,
        genotypedata_RR = genotypedata_RR[, -which(names(genotypedata_RR) == "Site")],
        additional_neutral = additional_neutral[, -which(names(additional_neutral) == "Site")],
        alleles_definitions_RR = alleles_definitions_RR,
        marker_info = marker_info_subset,
        record_hidden_alleles = record_hidden_alleles
      )
    }

    if (nloci == 0) {
      message("INFO: No valid loci with data found for site '", site, "'. Skipping.")
      next
    }

    locus_summary <- data.frame(
      patient_id = ids,
      n_available_d0 = 0,
      n_available_df = 0,
      n_comparable_loci = 0
    )
    
    is_locus_comparable <- matrix(FALSE, nrow = nids, ncol = nloci, 
                                  dimnames = list(ids, locinames))

    for (i in 1:nrow(locus_summary)) {
      patient_id <- locus_summary$patient_id[i]
      pattern_d0 <- paste0("\\b", patient_id, " Day 0\\b")
      pattern_df <- paste0("\\b", patient_id, " Day Failure\\b")
      # d0_row <- genotypedata_RR[grepl(paste0(patient_id, " Day 0"), genotypedata_RR$Sample.ID), ]
      # df_row <- genotypedata_RR[grepl(paste0(patient_id, " Day Failure"), genotypedata_RR$Sample.ID), ]
      d0_row <- genotypedata_RR[grepl(pattern_d0, genotypedata_RR$Sample.ID), ]
      df_row <- genotypedata_RR[grepl(pattern_df, genotypedata_RR$Sample.ID), ]
      
      if (nrow(d0_row) == 0 || nrow(df_row) == 0) { next }
      
      for (locus_name in locinames) {
        locus_cols <- grep(paste0("^", locus_name, "_"), colnames(genotypedata_RR), value = TRUE)
        has_d0_allele <- any(!is.na(d0_row[, locus_cols]))
        has_df_allele <- any(!is.na(df_row[, locus_cols]))
        
        if (has_d0_allele) {
          locus_summary$n_available_d0[i] <- locus_summary$n_available_d0[i] + 1
        }
        if (has_df_allele) {
          locus_summary$n_available_df[i] <- locus_summary$n_available_df[i] + 1
        }
        if (has_d0_allele && has_df_allele) {
          locus_summary$n_comparable_loci[i] <- locus_summary$n_comparable_loci[i] + 1
          is_locus_comparable[patient_id, locus_name] <- TRUE
        }
      }
    }
    
    local_sites_locus_summary[[site]] <- locus_summary
    engine_args$is_locus_comparable <- is_locus_comparable

    # Automated MCMC Loop
    full_loglik_history <- list()
    full_chain_results <- list()  
    total_iterations <- 0
    converged <- FALSE
    
    while (!converged && total_iterations < max_iterations) {
      total_iterations <- total_iterations + chunk_size
      
      chunk_results <- future.apply::future_lapply(1:n_chains, function(id) {
        full_args_for_chunk <- c(
          list(chain_id = id, nruns = chunk_size, burnin = 0, record_interval = 10),
          engine_args 
        )
        
        result <- do.call(mcmc_engine_function, full_args_for_chunk)
        if (data_type == "ampseq") {
          names(result)[names(result) == "state_classification"] <- "classification"
          names(result)[names(result) == "state_parameters"] <- "parameters"
          names(result)[names(result) == "state_loglikelihood"] <- "loglikelihood"
          names(result)[names(result) == "state_recoded0"] <- "alleles0"
          names(result)[names(result) == "state_recodedf"] <- "allelesf"
        }
        return(result)
      }, future.seed = TRUE)
      
      if (length(full_chain_results) == 0) {
        full_chain_results <- chunk_results
        full_loglik_history <- lapply(chunk_results, `[[`, "loglikelihood")
      } else {
        for (i in 1:n_chains) {
          full_chain_results[[i]]$parameters <- cbind(full_chain_results[[i]]$parameters, chunk_results[[i]]$parameters)
          full_chain_results[[i]]$classification <- cbind(full_chain_results[[i]]$classification, chunk_results[[i]]$classification)
          full_loglik_history[[i]] <- c(full_loglik_history[[i]], chunk_results[[i]]$loglikelihood)
          full_chain_results[[i]]$locus_lrs <- abind::abind(full_chain_results[[i]]$locus_lrs, chunk_results[[i]]$locus_lrs, along = 3)
          full_chain_results[[i]]$locus_dists <- abind::abind(full_chain_results[[i]]$locus_dists, chunk_results[[i]]$locus_dists, along = 3)
          
          if (record_hidden_alleles) {
            full_chain_results[[i]]$alleles0 <- abind::abind(full_chain_results[[i]]$alleles0, chunk_results[[i]]$alleles0, along = 3)
            full_chain_results[[i]]$allelesf <- abind::abind(full_chain_results[[i]]$allelesf, chunk_results[[i]]$allelesf, along = 3)
          }
        }
      }
      
      if (length(full_loglik_history) == 0 || length(full_loglik_history[[1]]) == 0) {
        cat("Log-likelihood history is not populated yet. Running another chunk...\n")
        next
      }
      
      mcmc_list_loglik <- coda::mcmc.list(lapply(full_loglik_history, coda::mcmc))
      n_samples <- nrow(mcmc_list_loglik[[1]])
      burn_in_end <- floor(burn_in_frac * n_samples)
      
      if (is.null(n_samples) || (n_samples - burn_in_end < 50)) { next }
      
      post_burn_mcmc <- stats::window(mcmc_list_loglik, start = burn_in_end + 1)
      r_hat <- try(coda::gelman.diag(post_burn_mcmc)$psrf[1, 1], silent = TRUE)
      ess <- try(coda::effectiveSize(post_burn_mcmc)[1], silent = TRUE)
      
      if (inherits(r_hat, "try-error") || inherits(ess, "try-error")) { next }
      
      r_hat_ok <- !is.na(r_hat) && r_hat < R_hat_threshold
      ess_ok <- !is.na(ess) && ess > ESS_threshold
      
      if (r_hat_ok && ess_ok) {
        converged <- TRUE
      }
    }
    
    
    if (length(full_loglik_history) > 0 && length(full_loglik_history[[1]]) > 0) {
      cat(paste("\nGenerating convergence diagnostics for site:", site, "\n"))
      # Call the updated function here
      generate_likelihood_diagnostics(
        all_chains_loglikelihood = full_loglik_history,
        site_name = site,
        output_folder = output_folder # Pass the main output folder
      )
    } else {
      cat(paste("WARNING: No MCMC output was generated for site '", site, "'. Skipping diagnostics.\n", sep=""))
    }
    
   if (length(full_chain_results) > 0) {
      num_total_samples_per_chain <- ncol(full_chain_results[[1]]$parameters)
      burn_in_samples_per_chain <- floor(burn_in_frac * num_total_samples_per_chain)
      if (burn_in_samples_per_chain < num_total_samples_per_chain) {
        keep_indices <- (burn_in_samples_per_chain + 1):num_total_samples_per_chain
        final_classification <- do.call(cbind, lapply(full_chain_results, function(x) x$classification[, keep_indices, drop = FALSE]))
        final_lrs <- do.call(abind::abind, list(lapply(full_chain_results, function(x) x$locus_lrs[,, keep_indices, drop = FALSE]), along = 3))
        final_dists <- do.call(abind::abind, list(lapply(full_chain_results, function(x) x$locus_dists[,, keep_indices, drop = FALSE]), along = 3))
        
        local_sites_locus_lrs[[site]] <- final_lrs
        local_sites_locus_dists[[site]] <- final_dists
        local_sites_classification[[site]] <- final_classification
        local_sites_ids[[site]] <- ids    
      } else {
        message("INFO: The number of iterations for '", site, "', was less than or equal to the burn-in. No post-burn-in samples to summarize.")
      }
    }
  }
  
  return(
    list(
      classifications = local_sites_classification,
      ids = local_sites_ids,
      locus_summary = local_sites_locus_summary,
      locus_lrs = local_sites_locus_lrs,
      locus_dists = local_sites_locus_dists,
      locinames = local_sites_locinames
    )
  )
}
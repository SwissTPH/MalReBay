run_all_sites <- function(
    genotypedata_latefailures,
    additional_genotypedata,
    locirepeats,
    n_chains,
    R_hat_threshold,
    ESS_threshold,
    chunk_size,
    max_iterations,
    burn_in_frac) {
  
  site_names <- unique(genotypedata_latefailures$Site)
  
  local_sites_classification <- list()
  local_sites_parameters <- list()
  local_sites_alleles0 <- list()
  local_sites_allelesf <- list()
  local_sites_ids <- list()
  
  for (site in site_names) {
    jobname <- site
    
    # Site-Specific Data Preparation (THE CORRECTED LOGIC)
    genotypedata_RR <- genotypedata_latefailures[genotypedata_latefailures$Site == site, -c(2)]
    additional_neutral <- additional_genotypedata[additional_genotypedata$Site == site, -c(2)]
    
    ids <- unique(gsub(" Day 0", "", genotypedata_RR$Sample.ID[grepl("Day 0", genotypedata_RR$Sample.ID)]))
    nids <- length(ids)
    
    # Calculate maxMOI safely
    marker_cols <- grep("_\\d+$", colnames(genotypedata_RR), value = TRUE)
    maxMOI <- if (length(marker_cols) > 0) max(as.numeric(gsub(".*_", "", marker_cols)), na.rm = TRUE) else 0
    maxalleles <- 30
    k <- rep(maxalleles, length(locirepeats))
    names(k) <- names(locirepeats)
    alleles_definitions_RR <- define_alleles(rbind(genotypedata_RR, additional_neutral), locirepeats, k)
    locinames <- names(alleles_definitions_RR)
    nloci <- length(locinames)
    if (nloci == 0) {
      message("INFO: No valid loci with data found for site '", site, "'. Skipping.")
      next
    }
    
    # Automated MCMC Loop
    full_loglik_history <- list()
    full_chain_results <- list()  
    total_iterations <- 0
    converged <- FALSE
    
    while (!converged && total_iterations < max_iterations) {
      total_iterations <- total_iterations + chunk_size
      
      chunk_results <- future.apply::future_lapply(1:n_chains, function(id) {
        run_one_chain(
          chain_id = id, nruns = chunk_size, burnin = 0, record_interval = 10,
          nids = nids, ids = ids, nloci = nloci, maxMOI = maxMOI, locinames = locinames,
          genotypedata_RR = genotypedata_RR, 
          additional_neutral = additional_neutral, 
          alleles_definitions_RR = alleles_definitions_RR
        )
      }, future.seed = TRUE)
      
      if (length(full_chain_results) == 0) {
        full_chain_results <- chunk_results
        full_loglik_history <- lapply(chunk_results, `[[`, "loglikelihood")
      } else {
        for (i in 1:n_chains) {
          full_chain_results[[i]]$parameters <- cbind(full_chain_results[[i]]$parameters, chunk_results[[i]]$parameters)
          full_chain_results[[i]]$classification <- cbind(full_chain_results[[i]]$classification, chunk_results[[i]]$classification)
          full_chain_results[[i]]$alleles0 <- abind::abind(full_chain_results[[i]]$alleles0, chunk_results[[i]]$alleles0, along = 3)
          full_chain_results[[i]]$allelesf <- abind::abind(full_chain_results[[i]]$allelesf, chunk_results[[i]]$allelesf, along = 3)
          full_loglik_history[[i]] <- c(full_loglik_history[[i]], chunk_results[[i]]$loglikelihood)
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
    
    num_total_samples_per_chain <- ifelse(length(full_chain_results) > 0, ncol(full_chain_results[[1]]$parameters), 0)
    if (num_total_samples_per_chain == 0) { next }
    burn_in_samples_per_chain <- floor(burn_in_frac * num_total_samples_per_chain)
    if (burn_in_samples_per_chain >= num_total_samples_per_chain) { next }
    
    keep_indices <- (burn_in_samples_per_chain + 1):num_total_samples_per_chain
    
    final_classification <- do.call(cbind, lapply(full_chain_results, function(x) x$classification[, keep_indices, drop=FALSE]))
    
    local_sites_classification[[site]] <- final_classification
    local_sites_ids[[site]] <- ids
  }
  
  return(
    list(
      classifications = local_sites_classification,
      ids = local_sites_ids
    )
  )
}
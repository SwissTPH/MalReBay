# =============================================================================
# R/stan_data_prep.R
# Step 1 of Stan transition — pure R, no Stan dependency
#
# Purpose:
#   Packages allele_utils.R output into a named list that
#   cmdstanr::cmdstan_model()$sample() accepts as its `data` argument.
#
#   Handles both marker families in a single pass, branching per locus on
#   marker_info$binning_method:
#     - "microsatellite" / "msp_glurp": alleles recoded via recodeallele()
#       against the binned allele_definitions, distance = numeric bin gap.
#     - "exact" (amplicon sequencing haplotypes): alleles recoded by mapping
#       each unique raw haplotype string to an integer directly (no binning),
#       distance = 0 for an exact match, 1 otherwise.
#
# Depends on:  R/allele_utils.R  (unchanged)
# Called by:   R/stan_interface.R  (Step 3)
#
# Debugging workflow:
#   sd <- prepare_stan_data(...)
#   validate_stan_data(sd)   # must print ALL CHECKS PASSED before running Stan
# =============================================================================


# -----------------------------------------------------------------------------
# prepare_stan_data()
#
# Arguments:
#   late_failures_site   data.frame — site-filtered rows, Site column removed
#   additional_site      data.frame — site-filtered additional data
#   allele_definitions   list       — output of define_alleles()
#   marker_info          data.frame — marker metadata, pre-filtered to site loci
#   ids                  character  — patient ID vector (length nids)
#   locinames            character  — locus name vector (length nloci), ordered
#   maxMOI               integer    — maximum multiplicity of infection
#   is_locus_comparable  matrix     — logical [nids × nloci]
#
# Returns:
#    Named list ready for cmdstanr::cmdstan_model()$sample(data = ...).
# -----------------------------------------------------------------------------
prepare_stan_data <- function(late_failures_site,
                              additional_site,
                              allele_definitions,
                              marker_info,
                              ids,
                              locinames,
                              maxMOI,
                              is_locus_comparable) {

  nids  <- length(ids)
  nloci <- length(locinames)

  # 1. Calculate the MOI per patient (max across all loci, both marker types).
  # Missing alleles are already true NA by the time data reaches here
  # (import_data()'s clean_data() step runs regardless of data_type), so a
  # plain is.na() check is sufficient for exact (ampseq) loci too.
  MOI0 <- rep(1L, nids); MOIf <- rep(1L, nids)
  for (i in seq_len(nids)) {
    for (j in seq_len(nloci)) {
      loci_cols <- grepl(paste0(locinames[j], "_"), colnames(late_failures_site))
      n0 <- sum(!is.na(late_failures_site[grepl(paste0("\\b", ids[i], " Day 0\\b"), late_failures_site$Sample.ID), loci_cols]))
      nf <- sum(!is.na(late_failures_site[grepl(paste0("\\b", ids[i], " recurrence\\b"), late_failures_site$Sample.ID), loci_cols]))
      MOI0[i] <- max(MOI0[i], n0)
      MOIf[i] <- max(MOIf[i], nf)
    }
  }

  # 2. Site-specific allele pruning & recoding
  recoded0 <- matrix(0L, nids, maxMOI * nloci)
  recodedf <- matrix(0L, nids, maxMOI * nloci)

  K_site <- integer(nloci)
  dist_mats <- vector("list", nloci)

  # 500-allele scratch buffer: ampseq loci can carry far more distinct
  # haplotypes than a microsatellite has length bins.
  tmp_counts <- matrix(0L, nloci, 500)

  day0_rows <- grepl("Day 0",      late_failures_site$Sample.ID)
  dayf_rows <- grepl("recurrence", late_failures_site$Sample.ID)
  pid_idx   <- setNames(seq_along(ids), ids)

  for (j in seq_len(nloci)) {
    locus <- locinames[j]
    cols  <- grepl(paste0("^", locus, "_"), colnames(late_failures_site))

    binning_method_j <- marker_info$binning_method[match(locus, marker_info$marker_id)]
    if (is.na(binning_method_j)) binning_method_j <- "exact"

    obs_trial <- as.matrix(late_failures_site[, cols, drop = FALSE])
    obs_add   <- if (nrow(additional_site) > 0) as.matrix(additional_site[, cols, drop = FALSE]) else matrix(nrow = 0, ncol = 0)

    if (binning_method_j == "exact") {
      # ---------------------------------------------------------------------
      # Amplicon sequencing (haplotype) locus: alleles are already discrete
      # strings, so map each unique observed haplotype directly to an
      # integer — no numeric binning via recodeallele().
      # ---------------------------------------------------------------------
      all_raw <- c(as.character(obs_trial), as.character(obs_add))
      all_raw <- all_raw[!is.na(all_raw) & nchar(trimws(all_raw)) > 0 & all_raw != "NA"]
      unique_alleles <- sort(unique(all_raw))

      K_site[j] <- length(unique_alleles)
      if (K_site[j] == 0) { K_site[j] <- 1L; unique_alleles <- "UNKNOWN" }

      id_map <- setNames(seq_along(unique_alleles), unique_alleles)

      for (row_idx in seq_len(nrow(late_failures_site))) {
        p_id  <- gsub(" Day 0| recurrence", "", late_failures_site$Sample.ID[row_idx])
        p_idx <- pid_idx[p_id]
        if (is.na(p_idx)) next

        alleles_row <- as.character(obs_trial[row_idx, ])
        alleles_row <- alleles_row[!is.na(alleles_row) & nchar(trimws(alleles_row)) > 0 & alleles_row != "NA"]
        n_fill <- min(length(alleles_row), maxMOI)
        if (n_fill == 0) next

        slots <- (maxMOI * (j - 1) + 1):(maxMOI * (j - 1) + n_fill)
        coded <- as.integer(id_map[alleles_row[seq_len(n_fill)]])
        coded[is.na(coded)] <- 0L

        if (day0_rows[row_idx]) recoded0[p_idx, slots] <- coded
        if (dayf_rows[row_idx]) recodedf[p_idx, slots] <- coded
      }

      # additional_counts: all trial alleles (Day 0 + recurrence), matching
      # the length-polymorphic convention (see note below the else branch).
      trial_raw <- as.character(obs_trial)
      trial_raw <- trial_raw[!is.na(trial_raw) & nchar(trimws(trial_raw)) > 0 & trial_raw != "NA"]
      if (length(trial_raw) > 0) {
        mapped <- id_map[trial_raw]
        mapped <- mapped[!is.na(mapped)]
        if (length(mapped) > 0) {
          cnts <- table(factor(mapped, levels = seq_along(unique_alleles)))
          tmp_counts[j, seq_along(unique_alleles)] <-
            tmp_counts[j, seq_along(unique_alleles)] + as.integer(cnts)
        }
      }
      if (nrow(obs_add) > 0) {
        add_raw <- as.character(obs_add)
        add_raw <- add_raw[!is.na(add_raw) & nchar(trimws(add_raw)) > 0 & add_raw != "NA"]
        if (length(add_raw) > 0) {
          mapped <- id_map[add_raw]
          mapped <- mapped[!is.na(mapped)]
          if (length(mapped) > 0) {
            cnts <- table(factor(mapped, levels = seq_along(unique_alleles)))
            tmp_counts[j, seq_along(unique_alleles)] <-
              tmp_counts[j, seq_along(unique_alleles)] + as.integer(cnts)
          }
        }
      }

      # Exact-match distance: 0 for the same haplotype, 1 for any other one.
      # Feeds method_int == 3 ("exact match only") in malrebay_model.stan.
      dist_mats[[j]] <- matrix(1, K_site[j], K_site[j]) - diag(K_site[j])

    } else {
      # ---------------------------------------------------------------------
      # Length-polymorphic locus (microsatellite / msp_glurp): recode via
      # the binned allele_definitions.
      # ---------------------------------------------------------------------
      all_obs_raw <- c(as.vector(obs_trial), as.vector(obs_add))
      all_obs_ids <- unique(sapply(all_obs_raw, function(x) recodeallele(allele_definitions[[j]], x)))
      all_obs_ids <- sort(all_obs_ids[!is.na(all_obs_ids) & all_obs_ids > 0])

      K_site[j] <- length(all_obs_ids)
      if (K_site[j] == 0) { K_site[j] <- 1; all_obs_ids <- 1 }

      id_map     <- setNames(seq_along(all_obs_ids), all_obs_ids)
      pruned_def <- allele_definitions[[j]][all_obs_ids, , drop = FALSE]

      for (row_idx in seq_len(nrow(late_failures_site))) {
        p_id  <- gsub(" Day 0| recurrence", "", late_failures_site$Sample.ID[row_idx])
        p_idx <- pid_idx[p_id]
        if (is.na(p_idx)) next

        for (c_idx in seq_len(ncol(obs_trial))) {
          old_id <- recodeallele(allele_definitions[[j]], obs_trial[row_idx, c_idx])
          if (!is.na(old_id) && old_id > 0) {
            new_id <- id_map[as.character(old_id)]
            slot <- (maxMOI * (j - 1) + c_idx)
            if (day0_rows[row_idx]) recoded0[p_idx, slot] <- as.integer(new_id)
            if (dayf_rows[row_idx]) recodedf[p_idx, slot] <- as.integer(new_id)
          }
        }
      }

      # Count ALL trial alleles (Day 0 + recurrence) into additional_counts.
      # Stan marginalises the recrudescence indicator R_i out, so we cannot
      # condition on classification the way the Gibbs sampler did. Using all
      # trial rows keeps the Dirichlet posterior well-informed and prevents
      # freq[j] from being under-estimated, which would inflate SLR and produce
      # false-positive recrudescence calls.
      if (nrow(obs_trial) > 0) {
        trial_ids_raw <- sapply(as.vector(obs_trial), function(x) recodeallele(allele_definitions[[j]], x))
        trial_ids_raw <- trial_ids_raw[!is.na(trial_ids_raw) & trial_ids_raw > 0]
        if (length(trial_ids_raw) > 0) {
          new_trial_ids <- id_map[as.character(trial_ids_raw)]
          cnts <- table(factor(new_trial_ids, levels = seq_along(all_obs_ids)))
          tmp_counts[j, seq_along(all_obs_ids)] <-
            tmp_counts[j, seq_along(all_obs_ids)] + as.integer(cnts)
        }
      }

      # Background additional_site data: population survey Day-0 samples,
      # not trial recurrence samples — safe to include without classification.
      if (nrow(obs_add) > 0) {
        add_ids_raw <- sapply(obs_add, function(x) recodeallele(allele_definitions[[j]], x))
        add_ids_raw <- add_ids_raw[!is.na(add_ids_raw) & add_ids_raw > 0]
        if (length(add_ids_raw) > 0) {
          new_add_ids <- id_map[as.character(add_ids_raw)]
          cnts <- table(factor(new_add_ids, levels = seq_along(all_obs_ids)))
          tmp_counts[j, seq_along(all_obs_ids)] <-
            tmp_counts[j, seq_along(all_obs_ids)] + as.integer(cnts)
        }
      }

      dist_mats[[j]] <- as.matrix(stats::dist(rowMeans(pruned_def)))
    }
  }

  max_K <- max(K_site)
  additional_counts <- tmp_counts[, 1:max_K, drop = FALSE]

  # 3. Assemble the pruned per-locus distance matrices into one padded array.
  dist_array <- array(0, dim = c(nloci, max_K, max_K))
  for (j in seq_len(nloci)) {
    if (K_site[j] > 0) {
      dist_array[j, 1:K_site[j], 1:K_site[j]] <- dist_mats[[j]]
    }
  }

  # 4. Method & threshold logic
  method_int <- sapply(locinames, function(ln) {
    m <- marker_info$binning_method[match(ln, marker_info$marker_id)]
    if(is.na(m)) return(3L)
    switch(as.character(m), "microsatellite" = 1L, "msp_glurp" = 2L, "exact" = 3L, 3L)
  })

  threshold_vec <- sapply(seq_len(nloci), function(j) {
    if (method_int[j] == 2L) {
      val <- marker_info$repeatlength[match(locinames[j], marker_info$marker_id)]
      return(if(is.na(val)) 0 else val)
    } else return(0)
  })

  # 5. Assemble and return
  # max_dist is derived from the actual largest allele distance observed,
  # not a hardcoded constant. This keeps log_dvect as small as needed.
  computed_max_dist <- max(1L, ceiling(max(dist_array)))

  return(list(
    N        = nids,
    J        = nloci,
    maxMOI   = as.integer(maxMOI),
    max_K    = as.integer(max_K),
    max_dist = as.integer(computed_max_dist),
    K        = as.array(K_site),
    recoded0 = recoded0,
    recodedf = recodedf,
    hidden0  = matrix(as.integer(recoded0 == 0), nids, nloci * maxMOI),
    hiddenf  = matrix(as.integer(recodedf == 0), nids, nloci * maxMOI),
    MOI0     = as.array(MOI0),
    MOIf     = as.array(MOIf),
    dist_array = dist_array,
    method_int = as.array(method_int),
    threshold  = as.array(threshold_vec),
    comparable = matrix(as.integer(is_locus_comparable), nids, nloci),
    additional_counts = additional_counts
  ))
}

validate_stan_data <- function(sd) {
  ok  <- TRUE
  chk <- function(label, test) {
    if (!isTRUE(test)) {
      message("  [FAIL] ", label)
      ok <<- FALSE
    } else {
      message("  [PASS] ", label)
    }
  }

  message("=== validate_stan_data ===")
  chk("N >= 1",       sd$N >= 1)
  chk("J >= 1",       sd$J >= 1)
  chk("max_K >= 1",   sd$max_K >= 1)
  chk("length(K) == J", length(sd$K) == sd$J)
  chk("dim(recoded0) == [N, J*maxMOI]", all(dim(sd$recoded0) == c(sd$N, sd$J * sd$maxMOI)))
  chk("additional_counts exists", !is.null(sd$additional_counts))
  chk("additional_counts has data", sum(sd$additional_counts) > 0)

  message(if (ok) "=== ALL CHECKS PASSED ===" else "=== VALIDATION FAILED ===")
  return(ok)
}

stan_data_only <- function(sd) {
  sd[!startsWith(names(sd), ".")]
}

#' Prepare Stan data for the amplicon-sequencing model
#'
#' Converts site-filtered haplotype data into a named list accepted by
#' \code{rstan::sampling(data = ...)} for the ampseq Stan model.
#' Haplotype strings are mapped to integers; allele counts from both trial
#' and additional-site data are used to inform the Dirichlet frequency prior.
#'
#' @param late_failures_site   Data frame of late-failure rows for one site
#'   (Site column removed).
#' @param additional_site      Data frame of additional background rows for
#'   one site, or \code{NULL} / zero-row data frame if absent.
#' @param marker_info          Data frame of marker metadata.
#' @param ids                  Character vector of patient IDs.
#' @param locinames            Character vector of locus names.
#' @param maxMOI               Integer maximum multiplicity of infection.
#' @param is_locus_comparable  Logical matrix \code{[nids x nloci]}.
#'
#' @return A named list ready for \code{rstan::sampling}. Private fields
#'   (names starting with \code{.}) carry diagnostics and are stripped by
#'   \code{stan_data_only} before sampling.
#' @noRd
prepare_stan_data_ampseq <- function(late_failures_site,
                                      additional_site,
                                      marker_info,
                                      ids,
                                      locinames,
                                      maxMOI,
                                      is_locus_comparable) {
  nids  <- length(ids)
  nloci <- length(locinames)

  # ---- 1. Per-locus MOI and recoding ----
  # MOI0 and MOIf are N x J (per locus), unlike length-polymorphic which is N x 1
  MOI0     <- matrix(0L, nids, nloci)
  MOIf     <- matrix(0L, nids, nloci)
  recoded0 <- matrix(0L, nids, maxMOI * nloci)
  recodedf <- matrix(0L, nids, maxMOI * nloci)

  K_site      <- integer(nloci)
  id_maps     <- vector("list", nloci)
  tmp_counts  <- matrix(0L, nloci, 500)  # 500-allele buffer; ampseq can have many haplotypes

  for (j in seq_len(nloci)) {
    locus <- locinames[j]
    cols  <- grepl(paste0("^", locus, "_"), colnames(late_failures_site))

    obs_trial <- as.matrix(late_failures_site[, cols, drop = FALSE])
    obs_add   <- if (!is.null(additional_site) && nrow(additional_site) > 0)
                   as.matrix(additional_site[, cols, drop = FALSE])
                 else matrix(character(0), nrow = 0, ncol = 0)

    # Collect all unique non-NA haplotype strings observed at this site
    all_raw <- c(as.vector(obs_trial), as.vector(obs_add))
    all_raw <- all_raw[!is.na(all_raw) & nchar(trimws(all_raw)) > 0 & all_raw != "NA"]
    unique_alleles <- sort(unique(all_raw))

    K_site[j] <- length(unique_alleles)
    if (K_site[j] == 0) {
      K_site[j] <- 1L
      unique_alleles <- "UNKNOWN"
    }

    # String -> integer map (1-indexed)
    id_map      <- setNames(seq_along(unique_alleles), unique_alleles)
    id_maps[[j]] <- id_map

    day0_rows <- grepl("Day 0",      late_failures_site$Sample.ID)
    dayf_rows <- grepl("recurrence", late_failures_site$Sample.ID)
    pid_idx   <- setNames(seq_along(ids), ids)

    # Fill recoded0 / recodedf and compute per-locus MOI
    for (row_idx in seq_len(nrow(late_failures_site))) {
      p_id  <- gsub(" Day 0| recurrence", "", late_failures_site$Sample.ID[row_idx])
      p_idx <- pid_idx[p_id]
      if (is.na(p_idx)) next

      alleles_row <- as.character(obs_trial[row_idx, ])
      alleles_row <- alleles_row[!is.na(alleles_row) & nchar(trimws(alleles_row)) > 0 & alleles_row != "NA"]
      n_alleles   <- length(alleles_row)
      if (n_alleles == 0) next

      n_fill <- min(n_alleles, maxMOI)
      slots  <- (maxMOI * (j - 1) + 1):(maxMOI * (j - 1) + n_fill)
      coded  <- as.integer(id_map[alleles_row[1:n_fill]])
      coded[is.na(coded)] <- 0L  # safety: unknown haplotype -> 0

      if (day0_rows[row_idx]) {
        recoded0[p_idx, slots] <- coded
        MOI0[p_idx, j]         <- as.integer(n_fill)
      }
      if (dayf_rows[row_idx]) {
        recodedf[p_idx, slots] <- coded
        MOIf[p_idx, j]         <- as.integer(n_fill)
      }
    }

    # BUG 2 equivalent fix: count BOTH trial and additional alleles for freq prior
    # Trial data
    trial_raw <- as.vector(obs_trial)
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
    # Additional site data
    if (nrow(obs_add) > 0) {
      add_raw <- as.vector(obs_add)
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
  }

  max_K             <- max(K_site)
  additional_counts <- tmp_counts[, 1:max_K, drop = FALSE]

  # ---- 2. Hidden flags ----
  # For ampseq, all observed alleles are fully typed: hidden = 0 everywhere an
  # allele was recorded, and position remains 0 (beyond MOI) otherwise.
  # We still build the matrices for model generality.
  hidden0 <- matrix(0L, nids, maxMOI * nloci)
  hiddenf <- matrix(0L, nids, maxMOI * nloci)

  # ---- 3. Assemble ----
  return(list(
    N                 = nids,
    J                 = nloci,
    maxMOI            = as.integer(maxMOI),
    max_K             = as.integer(max_K),
    K                 = as.array(K_site),
    MOI0              = MOI0,          # N x J integer matrix
    MOIf              = MOIf,          # N x J integer matrix
    recoded0          = recoded0,
    recodedf          = recodedf,
    hidden0           = hidden0,
    hiddenf           = hiddenf,
    comparable        = matrix(as.integer(is_locus_comparable), nids, nloci),
    additional_counts = additional_counts,
    # id_maps kept for diagnostics / back-mapping (not passed to Stan)
    .id_maps          = id_maps,
    .locinames        = locinames
  ))
}

#' Validate Stan data list for the ampseq model
#'
#' Runs a series of dimension and range checks on the list returned by
#' \code{prepare_stan_data_ampseq}.
#'
#' @param sd Named list produced by \code{prepare_stan_data_ampseq}.
#' @return Logical scalar — \code{TRUE} if all checks pass.
#' @noRd
validate_stan_data_ampseq <- function(sd) {
  ok  <- TRUE
  chk <- function(label, test) {
    if (!isTRUE(test)) { message("  [FAIL] ", label); ok <<- FALSE }
    else                { message("  [PASS] ", label) }
  }
  message("=== validate_stan_data_ampseq ===")
  chk("N >= 1",                     sd$N >= 1)
  chk("J >= 1",                     sd$J >= 1)
  chk("max_K >= 1",                 sd$max_K >= 1)
  chk("length(K) == J",             length(sd$K) == sd$J)
  chk("dim(MOI0) == [N, J]",        all(dim(sd$MOI0) == c(sd$N, sd$J)))
  chk("dim(recoded0) == [N,J*maxMOI]", all(dim(sd$recoded0) == c(sd$N, sd$J * sd$maxMOI)))
  chk("additional_counts has data", sum(sd$additional_counts) > 0)
  chk("no recoded0 > max_K",        max(sd$recoded0) <= sd$max_K)
  chk("no recodedf > max_K",        max(sd$recodedf) <= sd$max_K)
  message(if (ok) "=== ALL CHECKS PASSED ===" else "=== VALIDATION FAILED ===")
  return(ok)
}

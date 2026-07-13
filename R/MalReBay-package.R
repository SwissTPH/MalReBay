#' @keywords internal
"_PACKAGE"

#' @importFrom rlang .data
#' @importFrom stats dist rbinom sd setNames
#' @import methods
NULL

utils::globalVariables(c(
  "patient_id", "n_available_d0", "n_available_df", "n_comparable_loci",
  "Mean_LR", "Probability", "N_Comparable_Loci", "Sample.ID", "Site",
  "Interpretation", "Mean_Distance", "Patient.ID", "Number_Matches",
  "Base.ID", "cls", "id", ".data",
  "ESS_threshold", "R_hat_threshold", "full_loglik_history", "site", "p_recrud"
))
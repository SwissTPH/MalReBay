#' @keywords internal
"_PACKAGE"

# Add stanmodels and p_recrud here
utils::globalVariables(c("ESS_threshold", "R_hat_threshold", "full_loglik_history", 
                         "site", "stanmodels", "p_recrud"))

#' @import Rcpp
#' @import methods
#' @importFrom rstan sampling
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
#' @useDynLib MalReBay, .registration = TRUE
NULL


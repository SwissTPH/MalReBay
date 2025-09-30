#' Generate and Save MCMC Convergence Diagnostics
#'
#' @description
#' Calculates standard MCMC convergence diagnostics (Gelman-Rubin R-hat,
#' Effective Sample Size) for the log-likelihood and saves diagnostic plots
#' (traceplot, Gelman plot, histogram, ACF).
#'
#' @details
#' This function takes the log-likelihood history from multiple MCMC chains,
#' cleans the data, and uses the `coda` package to compute and display diagnostics.
#' It prints the diagnostic summaries to the console (if `verbose = TRUE`) and saves
#' a series of PNG plots to a site-specific subdirectory within the output folder.
#' The function includes error handling to prevent crashes if `coda` functions fail.
#'
#' @param all_chains_loglikelihood A list where each element is a numeric vector
#'   representing the log-likelihood history of one MCMC chain.
#' @param site_name A character string for labeling plots and creating the output
#'   subdirectory.
#' @param output_folder A character string specifying the main directory for results.
#' @param verbose A logical. If `TRUE` (the default), prints the Gelman-Rubin
#'   and ESS summaries to the console.
#'
#' @return A list containing the calculated diagnostic results (`gelman`, `ess`, etc.).
#'
#' @importFrom coda as.mcmc.list mcmc varnames gelman.diag effectiveSize gelman.plot
#' @importFrom grDevices png dev.off rainbow n2mfrow
#' @importFrom graphics matplot legend title hist par
#' @importFrom stats acf var
#'
#' @keywords internal
#' @noRd
#'
generate_likelihood_diagnostics <- function(all_chains_loglikelihood,
                                            site_name,
                                            output_folder,
                                            verbose = TRUE) {
  
  # Create a dedicated folder for this site's diagnostics
  site_dir <- file.path(output_folder, "convergence_diagnosis", site_name)
  dir.create(site_dir, recursive = TRUE, showWarnings = FALSE)
  clean_chains <- lapply(all_chains_loglikelihood, function(x) {
    x[!is.finite(x)] <- NA
    x <- x[!is.na(x)]
    return(x)
  })
  
  if (any(sapply(clean_chains, length) == 0)) {
    stop("One or more chains contained no valid log-likelihood values after cleaning.")
  }
  
  loglikelihood_mcmc <- tryCatch({
    coda::as.mcmc.list(lapply(clean_chains, function(x) coda::mcmc(matrix(x, ncol = 1))))
  }, error = function(e) {
    warning("Could not create mcmc.list for site '", site_name, "': ", e$message, call. = FALSE)
    return(NULL)
  })
  
  if (is.null(loglikelihood_mcmc)) return(invisible(NULL))
  coda::varnames(loglikelihood_mcmc) <- "loglikelihood"
  
  # Gelman-Rubin Diagnostic
  gelman_result <- tryCatch({
    coda::gelman.diag(loglikelihood_mcmc, autoburnin = FALSE, multivariate = FALSE)
  }, error = function(e) {
    warning("Gelman-Rubin calculation failed for site '", site_name, "': ", e$message, call. = FALSE)
    return(NULL)
  })
  
  # Effective Sample Size
  ess_result <- tryCatch({
    coda::effectiveSize(loglikelihood_mcmc)
  }, error = function(e) {
    warning("ESS calculation failed for site '", site_name, "': ", e$message, call. = FALSE)
    return(NULL)
  })
  
  # Conditional Console Output
  if (verbose) {
    cat("Gelman-Rubin R-hat for Log-Likelihood:\n")
    if (!is.null(gelman_result)) {
      print(gelman_result)
    } else {
      cat("  -> Calculation failed.\n")
    }
    
    cat("\nEffective Sample Size (ESS):\n")
    if (!is.null(ess_result)) {
      print(ess_result)
    } else {
      cat("  -> Calculation failed.\n")
    }
  }
  
  # Diagnostic Plots (saved to file)
  plot_base <- function(filename, expr) {
    grDevices::png(file.path(site_dir, filename), width = 800, height = 600)
    on.exit(grDevices::dev.off()) 
    try(expr, silent = TRUE)
  }
  
  plot_base("loglikelihood_traceplot.png", {
    colors <- grDevices::rainbow(length(loglikelihood_mcmc))
    graphics::matplot(do.call(cbind, lapply(loglikelihood_mcmc, as.numeric)),
                      type = "l", lty = 1, col = colors,
                      main = paste(site_name, "- Traceplot of Log-Likelihood"),
                      ylab = "Log-Likelihood", xlab = "Iterations")
    graphics::legend("topright", legend = paste("Chain", seq_along(loglikelihood_mcmc)),
                     col = colors, lty = 1, cex = 0.8, box.lty = 0)
  })
  
  plot_base("gelman_loglikelihood.png", {
    coda::gelman.plot(loglikelihood_mcmc, autoburnin = FALSE, ask = FALSE)
    graphics::title(main = paste(site_name, "- Gelman-Rubin for Log-Likelihood"), line = 2.5)
  })
  
  plot_base("loglikelihood_histogram.png", {
    graphics::hist(unlist(clean_chains), breaks = 40, main = paste(site_name, "- Histogram of Log-Likelihood"),
                   xlab = "Log-Likelihood", col = "steelblue", border = "white")
  })
  
  plot_base("loglikelihood_acf.png", {
    n_chains <- length(clean_chains)
    plot_dims <- grDevices::n2mfrow(n_chains)
    graphics::par(mfrow = plot_dims)
    for (i in seq_along(clean_chains)) {
      stats::acf(clean_chains[[i]], main = paste("ACF - Chain", i), lag.max = 50)
    }
  })
  
  # Return summary invisibly
  invisible(list(
    gelman     = gelman_result,
    ess        = ess_result,
    n_valid    = sapply(clean_chains, length)
  ))
}
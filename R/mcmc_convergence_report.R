#' Generate and Plot MCMC Convergence Diagnostics
#'
#' @description
#' Calculates standard MCMC convergence diagnostics (Gelman-Rubin R-hat,
#' Effective Sample Size) and generates diagnostic plots (traceplot, Gelman plot,
#' histogram, ACF) for the log-likelihood.
#'
#' @details
#' This function takes the log-likelihood history from multiple MCMC chains,
#' cleans the data, and uses the `coda` package to compute diagnostics.
#' It can either display plots in the active R graphics device or save them as
#' PNG files to a specified directory.
#'
#' @param all_chains_loglikelihood A list where each element is a numeric vector
#'   representing the log-likelihood history of one MCMC chain.
#' @param site_name A character string for labeling plots and creating the output
#'   subdirectory when saving.
#' @param save_plot A logical. If `TRUE`, plots are saved as PNG files. If `FALSE`
#'   (the default), plots are displayed in the active R graphics device.
#' @param output_folder A character string specifying the path to a directory where
#'   plots should be saved. Only used if `save_plot = TRUE`.
#' @param verbose A logical. If `TRUE` (the default), prints the Gelman-Rubin
#'   and ESS summaries to the console.
#'
#' @return An invisible list containing the calculated diagnostic results:
#' \itemize{
#'   \item `gelman`: The Gelman-Rubin diagnostic result.
#'   \item `ess`: The Effective Sample Size result.
#'   \item `n_valid`: The number of valid (finite) samples per chain.
#' }
#'
#' @importFrom coda as.mcmc.list mcmc varnames gelman.diag effectiveSize gelman.plot
#' @importFrom grDevices png dev.off rainbow n2mfrow
#' @importFrom graphics matplot legend title hist par
#' @importFrom stats acf var
#'
#' @export
#' @examples
#' # --- Example 1: Interactive Plotting ---
#' # To view plots in the RStudio Plots pane (will ask before showing each plot)
#' \dontrun{
#'   generate_likelihood_diagnostics(all_chains, site_name = "Test Site")
#' }
#'
#' # --- Example 2: Saving Plots to a File ---
#' # Create a temporary directory for the output
#' temp_dir <- tempdir()
#' # Simulate some dummy MCMC chains for demonstration
#' all_chains <- list(
#'   chain1 = rnorm(1000, mean = -150, sd = 5),
#'   chain2 = rnorm(1000, mean = -150, sd = 5)
#' )
#' # Run the function to save plots
#' generate_likelihood_diagnostics(all_chains,
#'                                 site_name = "Test Site",
#'                                 save_plot = TRUE,
#'                                 output_folder = temp_dir)
#'
#' # Check that the files were created
#' list.files(file.path(temp_dir, "convergence_diagnosis", "Test Site"))
#'
generate_likelihood_diagnostics <- function(all_chains_loglikelihood,
                                            site_name,
                                            save_plot = FALSE,
                                            output_folder = NULL,
                                            verbose = TRUE) {
  
  if (save_plot && is.null(output_folder)) {
    stop("`output_folder` must be provided when `save_plot = TRUE`.", call. = FALSE)
  }
  
  site_dir <- NULL
  if (save_plot) {
    site_dir <- file.path(output_folder, "convergence_diagnosis", site_name)
    
    if (dir.exists(site_dir)) {
      if (verbose) message("INFO: Removing existing convergence diagnosis folder for site: ", site_name)
      unlink(site_dir, recursive = TRUE, force = TRUE)
    }
    
    dir.create(site_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Data Cleaning
  clean_chains <- lapply(all_chains_loglikelihood, function(x) {
    x[!is.finite(x)] <- NA
    return(x[!is.na(x)])
  })
  
  if (any(sapply(clean_chains, length) == 0)) {
    warning("For site '", site_name, "', one or more chains had no valid log-likelihood values after cleaning. Skipping plots.", call. = FALSE)
    return(invisible(NULL))
  }
  
  loglikelihood_mcmc <- tryCatch({
    coda::as.mcmc.list(lapply(clean_chains, function(x) coda::mcmc(matrix(x, ncol = 1))))
  }, error = function(e) {
    warning("Could not create mcmc.list for site '", site_name, "': ", e$message, call. = FALSE)
    return(NULL)
  })
  
  if (is.null(loglikelihood_mcmc)) return(invisible(NULL))
  coda::varnames(loglikelihood_mcmc) <- "loglikelihood"
  
  gelman_result <- tryCatch(
    coda::gelman.diag(loglikelihood_mcmc, autoburnin = FALSE, multivariate = FALSE),
    error = function(e) {
      warning("Gelman-Rubin calculation failed for site '", site_name, "': ", e$message, call. = FALSE)
      return(NULL)
    })
  
  ess_result <- tryCatch(
    coda::effectiveSize(loglikelihood_mcmc),
    error = function(e) {
      warning("ESS calculation failed for site '", site_name, "': ", e$message, call. = FALSE)
      return(NULL)
    })
  
  if (verbose) {
    cat("--- Convergence Diagnostics for:", site_name, "---\n")
    if (!is.null(gelman_result)) { cat("Gelman-Rubin R-hat:\n"); print(gelman_result) }
    if (!is.null(ess_result)) { cat("\nEffective Sample Size (ESS):\n"); print(ess_result) }
  }
  
  plot_trace <- function() {
    colors <- grDevices::rainbow(length(loglikelihood_mcmc))
    graphics::matplot(do.call(cbind, lapply(loglikelihood_mcmc, as.numeric)),
                      type = "l", lty = 1, col = colors,
                      main = paste(site_name, "- Traceplot"),
                      ylab = "Log-Likelihood", xlab = "Iterations")
    graphics::legend("topright", legend = paste("Chain", seq_along(loglikelihood_mcmc)),
                     col = colors, lty = 1, cex = 0.8, box.lty = 0)
  }
  
  plot_gelman <- function() {
    coda::gelman.plot(loglikelihood_mcmc, autoburnin = FALSE, ask = FALSE)
  }
  
  plot_hist <- function() {
    graphics::hist(unlist(clean_chains), breaks = 40, main = paste(site_name, "- Histogram"),
                   xlab = "Log-Likelihood", col = "steelblue", border = "white")
  }
  
  save_plot_with_feedback <- function(plot_function, file_path) {
    tryCatch({
      grDevices::png(file_path, width = 800, height = 600)
      plot_function()
    }, error = function(e) {
      warning("Failed to generate plot '", basename(file_path), "' for site '", site_name, "'.\n  Reason: ", e$message, call. = FALSE)
    }, finally = {
      if (grDevices::dev.cur() != 1) grDevices::dev.off()
    })
  }
  
  # Generate and Save Plots
  if (save_plot) {
    save_plot_with_feedback(plot_trace, file.path(site_dir, "loglikelihood_traceplot.png"))
    save_plot_with_feedback(plot_gelman, file.path(site_dir, "gelman_loglikelihood.png"))
    save_plot_with_feedback(plot_hist, file.path(site_dir, "loglikelihood_histogram.png"))
    
    acf_path <- file.path(site_dir, "loglikelihood_acf.png")
    tryCatch({
      grDevices::png(acf_path, width = 800, height = 600)
      n_chains <- length(clean_chains)
      old_par <- graphics::par(mfrow = grDevices::n2mfrow(n_chains))
      on.exit(graphics::par(old_par), add = TRUE)
      
      for (i in seq_along(clean_chains)) {
        stats::acf(clean_chains[[i]], main = paste("ACF - Chain", i), lag.max = 50)
      }
    }, error = function(e) {
      warning("Failed to generate ACF plot for site '", site_name, "'.\n  Reason: ", e$message, call. = FALSE)
    }, finally = {
      if (grDevices::dev.cur() != 1) grDevices::dev.off()
    })
    
  } else {
    old_par <- graphics::par(no.readonly = TRUE)   # save old par settings
    on.exit(graphics::par(old_par), add = TRUE)   # restore old par on exit
    
    graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))  # 2x2 grid with margins
    
    colors <- grDevices::rainbow(length(loglikelihood_mcmc))
    graphics::matplot(do.call(cbind, lapply(loglikelihood_mcmc, as.numeric)),
                      type = "l", lty = 1, col = colors,
                      main = paste(site_name, "- Traceplot"),
                      ylab = "Log-Likelihood", xlab = "Iterations")
    graphics::legend("topright", legend = paste("Chain", seq_along(loglikelihood_mcmc)),
                     col = colors, lty = 1, cex = 0.8, box.lty = 0)
    
    coda::gelman.plot(loglikelihood_mcmc, autoburnin = FALSE, ask = FALSE)
    
    graphics::hist(unlist(clean_chains), breaks = 40,
                   main = paste(site_name, "- Histogram"),
                   xlab = "Log-Likelihood", col = "steelblue", border = "white")
    
    stats::acf(clean_chains[[1]], main = paste(site_name, "- ACF (Chain 1)"), lag.max = 50)
  }
  
  invisible(list(
    gelman     = gelman_result,
    ess        = ess_result,
    n_valid    = sapply(clean_chains, length)
  ))
}
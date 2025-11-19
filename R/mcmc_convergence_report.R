#' Generate and Plot MCMC Convergence Diagnostics
#'
#' @description
#' Calculates standard MCMC convergence diagnostics and generates diagnostic plots
#' in a multi-stage process: 1) Gelman-Rubin, 2) Traceplot & Histogram, 3) All ACF plots.
#'
#' @param all_chains_loglikelihood A list where each element is a numeric vector
#'   representing the log-likelihood history of one MCMC chain.
#' @param site_name A character string for labeling plots.
#' @param save_plot A logical. If `TRUE`, plots are saved as PNG files.
#' @param output_folder A character string specifying the path to save plots.
#' @param verbose A logical. If `TRUE`, prints diagnostic summaries.
#'
#' @return An invisible list containing the calculated diagnostic results.
#'
#' @importFrom coda as.mcmc.list mcmc varnames gelman.diag effectiveSize gelman.plot
#' @importFrom grDevices png dev.off rainbow n2mfrow
#' @importFrom graphics matplot legend title hist par plot text abline lines layout plot.new
#' @importFrom stats acf
#'
#' @export
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
    if (dir.exists(site_dir)) unlink(site_dir, recursive = TRUE)
    dir.create(site_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  clean_chains <- lapply(all_chains_loglikelihood, function(x) x[is.finite(x)])
  clean_chains <- Filter(function(x) length(x) >= 2, clean_chains)
  
  if (length(clean_chains) == 0) {
    warning("No valid chains for site '", site_name, "'. Skipping.", call. = FALSE)
    return(invisible(NULL))
  }
  
  loglikelihood_mcmc <- tryCatch({
    mlist <- lapply(clean_chains, coda::mcmc)
    mclist <- coda::as.mcmc.list(mlist)
    coda::varnames(mclist) <- " "
    mclist
  }, error = function(e) { return(NULL) })
  
  if (is.null(loglikelihood_mcmc)) return(invisible(NULL))
  gelman_result <- NULL
  if (length(loglikelihood_mcmc) > 1) {
    gelman_result <- tryCatch(coda::gelman.diag(loglikelihood_mcmc, autoburnin = FALSE), error = function(e) NULL)
  }
  ess_result <- tryCatch(coda::effectiveSize(loglikelihood_mcmc), error = function(e) NULL)
  
  if (verbose) {
    cat("Convergence Diagnostics for:", site_name, "\n")
    if (!is.null(gelman_result)) { cat("Gelman-Rubin R-hat:\n"); print(gelman_result) }
    if (!is.null(ess_result)) { cat("\nEffective Sample Size (ESS):\n"); print(ess_result) }
  }
  
  plot_diagnostics <- function() {
    
    #  Gelman-Rubin Plot
    if (length(loglikelihood_mcmc) > 1) {
      tryCatch({
        coda::gelman.plot(loglikelihood_mcmc, autoburnin = FALSE,
                          col = c("black", "indianred"))
        graphics::title(main = "Gelman-Rubin Diagnostic")
      }, error = function(e) {
        plot(1, type="n", axes=FALSE, ann=FALSE); text(1, 1, "Gelman plot unavailable", col = "red")
      })
    } else {
      plot.new(); text(0.5, 0.5, "Gelman plot not applicable\n(requires > 1 chain)")
    }
    
    # Pause for user to see the first plot in an interactive session
    if (interactive()) readline(prompt="Press [enter] to continue to next plot...")
    
    # Traceplot and Histogram
    layout(matrix(c(1, 2), nrow = 2, byrow = TRUE))
    par(mar = c(4.1, 4.1, 3.5, 1.1))
    
    # Traceplot
    colors <- grDevices::rainbow(length(loglikelihood_mcmc))
    matplot(do.call(cbind, lapply(loglikelihood_mcmc, as.numeric)),
            type = "l", lty = 1, col = colors,
            main = "Traceplot", ylab = "Log-Likelihood", xlab = "Iterations")
    legend("topright", legend = paste("Chain", seq_along(loglikelihood_mcmc)),
           col = colors, lty = 1, bty = "n", cex = 0.8)
    
    # Histogram
    hist(unlist(clean_chains), breaks = 40,
         main = "Posterior Density",
         xlab = "Log-Likelihood", col = "steelblue", border = "white")
    

    if (interactive()) readline(prompt="Press [enter] to continue to ACF plots...")
    
    # ACF Plots
    n_chains <- length(clean_chains)
    plot_dims <- grDevices::n2mfrow(n_chains) 
    par(mfrow = plot_dims, mar = c(4, 4, 3, 1))
    
    for (i in seq_along(clean_chains)) {
      # Plot ACF for the current chain
      stats::acf(clean_chains[[i]], lag.max = 50, main = "")
      # Add a clear title
      graphics::title(main = paste("Autocorrelation - Chain", i), cex.main = 1.2)
    }
    
    # Reset plotting parameters to default
    layout(1)
    par(mfrow = c(1, 1))
  }
  
  if (save_plot) {
    
    # Gelman Plot
    gelman_path <- file.path(site_dir, paste0("1_gelman_plot_", site_name, ".png"))
    grDevices::png(gelman_path, width = 1000, height = 700, res = 120)
    if (length(loglikelihood_mcmc) > 1) {
      tryCatch(coda::gelman.plot(loglikelihood_mcmc, autoburnin = FALSE, main = "Gelman-Rubin Diagnostic", col = c("black", "indianred")),
               error = function(e){ plot.new(); text(0.5, 0.5, "Gelman plot failed", col="red") })
    } else {
      plot.new(); text(0.5, 0.5, "Gelman plot not applicable\n(requires > 1 chain)")
    }
    grDevices::dev.off()
    
    # File 2: Traceplot and Histogram
    trace_hist_path <- file.path(site_dir, paste0("2_trace_hist_plot_", site_name, ".png"))
    grDevices::png(trace_hist_path, width = 1000, height = 1000, res = 120)
    layout(matrix(c(1, 2), nrow = 2, byrow = TRUE)); par(mar = c(4.1, 4.1, 3.5, 1.1))
    colors <- grDevices::rainbow(length(loglikelihood_mcmc)); matplot(do.call(cbind, lapply(loglikelihood_mcmc, as.numeric)), type = "l", lty = 1, col = colors, main = "Traceplot", ylab = "Log-Likelihood", xlab = "Iterations"); legend("topright", legend = paste("Chain", seq_along(loglikelihood_mcmc)), col = colors, lty = 1, bty = "n", cex = 0.8)
    hist(unlist(clean_chains), breaks = 40, main = "Posterior Density", xlab = "Log-Likelihood", col = "steelblue", border = "white")
    grDevices::dev.off()
    
    # File 3: ACF Plots
    acf_path <- file.path(site_dir, paste0("3_acf_plots_", site_name, ".png"))
    grDevices::png(acf_path, width = 1200, height = 1000, res = 120)
    n_chains <- length(clean_chains); plot_dims <- grDevices::n2mfrow(n_chains); par(mfrow = plot_dims, mar = c(4, 4, 3, 1))
    for (i in seq_along(clean_chains)) { stats::acf(clean_chains[[i]], lag.max = 50, main = ""); graphics::title(main = paste("Autocorrelation - Chain", i), cex.main = 1.2) }
    grDevices::dev.off()
    
    layout(1); par(mfrow = c(1,1))
    
  } else {
   
    plot_diagnostics()
  }
  
  invisible(list(gelman = gelman_result, ess = ess_result))
}
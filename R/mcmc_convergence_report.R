#' Generate an Allele Frequency Distribution Plot
#'
#' @description
#' Creates and saves a faceted bar plot visualizing the frequency distribution
#' of raw allele lengths. The plot compares samples from "Day 0" against
#' "Day of Failure" for each genetic marker.
#'
#' @details
#' This function serves as an internal utility for exploratory data analysis. It
#' processes a wide-format data frame, calculates allele frequencies, and uses
#' `ggplot2` to generate a comprehensive visualization. It can operate in either
#' a site-specific or an overall mode. The resulting plot is saved as a PDF.
#'
#' @param raw_data_df A data frame in a 'wide' format. It must contain:
#'   \itemize{
#'     \item A `Sample.ID` column where entries contain the strings "Day 0" or "Day Failure".
#'     \item Marker columns named in the format `LocusName_AlleleNumber` (e.g.,
#'       `TA40_1`, `TA40_2`).
#'   }
#' @param site_name A character string. The name of the site for the plot title
#'   and filename. If `NULL`, an "Overall" plot is generated.
#' @param output_folder A character string specifying the directory where the
#'   plot PDF will be saved.
#' @param plots_per_row An integer setting the number of marker plots to display
#'   per row in the faceted plot.
#' @param plot_width,plot_height Numeric values for the width and height of the
#'   output PDF file in inches.
#'
#' @return Invisibly returns the `ggplot` object, allowing for further customization.
#'
#' @importFrom dplyr mutate across all_of if_else group_by summarise n ungroup
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_bar facet_wrap scale_fill_manual labs theme_bw theme element_rect element_text ggsave
#'
#' @keywords internal
#' @noRd
#'

generate_likelihood_diagnostics <- function(all_chains_loglikelihood = full_loglik_history,
                                            site_name = site,
                                            output_folder = output_folder) {
  # Create a dedicated folder for this site's diagnostics
  site_dir <- file.path(output_folder, "convergence_diagnosis", site_name)
  dir.create(site_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Clean chains
  clean_chains <- lapply(all_chains_loglikelihood, function(x) {
    x[!is.finite(x)] <- NA
    x <- x[!is.na(x)]
    return(x)
  })
  
  if (any(sapply(clean_chains, length) == 0)) {
    stop("One or more chains contained no valid log-likelihood values after cleaning.")
  }
  
  # Convert to mcmc.list
  loglikelihood_mcmc <- tryCatch({
    coda::as.mcmc.list(lapply(clean_chains, function(x) coda::mcmc(matrix(x, ncol = 1))))
  }, error = function(e) {
    message("Error creating mcmc.list: ", e$message)
    return(NULL)
  })
  
  if (is.null(loglikelihood_mcmc)) return(invisible(NULL))
  coda::varnames(loglikelihood_mcmc) <- "loglikelihood"
  
  # Gelman-Rubin Diagnostic
  gelman_result <- tryCatch({
    coda::gelman.diag(loglikelihood_mcmc, autoburnin = FALSE, multivariate = FALSE)
  }, error = function(e) {
    message("Gelman-Rubin failed: ", e$message)
    chain_means <- sapply(clean_chains, mean, na.rm = TRUE)
    chain_vars  <- sapply(clean_chains, stats::var, na.rm = TRUE)
    n <- min(sapply(clean_chains, length))
    B <- n * stats::var(chain_means)
    W <- mean(chain_vars)
    Rhat <- sqrt(((n - 1)/n * W + B/n) / W)
    list(
      psrf  = matrix(c(Rhat, NA), nrow = 1, dimnames = list("loglikelihood", c("Point est.", "Upper C.I."))),
      mpsrf = Rhat
    )
  })
  
  cat("Gelman-Rubin R-hat for Log-Likelihood:\n")
  print(gelman_result)
  
  # Effective Sample Size
  ess_result <- tryCatch({
    coda::effectiveSize(loglikelihood_mcmc)
  }, error = function(e) {
    message("ESS calculation failed: ", e$message)
    sapply(clean_chains, length)
  })
  
  cat("\nEffective Sample Size (ESS):\n")
  print(ess_result)
  
  
  # Diagnostic Plots

  plot_base <- function(filename, expr) {
    grDevices::png(file.path(site_dir, filename), width = 800, height = 600)
    try(expr, silent = TRUE)
    grDevices::dev.off()
  }
  
  # Traceplot with legend
  plot_base("loglikelihood_traceplot.png", {
    colors <- grDevices::rainbow(length(loglikelihood_mcmc))
    graphics::matplot(do.call(cbind, lapply(loglikelihood_mcmc, as.numeric)),
            type = "l", lty = 1, col = colors,
            main = paste(site_name, "- Traceplot of Log-Likelihood"),
            ylab = "Log-Likelihood", xlab = "Iterations")
    graphics::legend("topright", legend = paste("Chain", seq_along(loglikelihood_mcmc)),
           col = colors, lty = 1, cex = 0.8, box.lty = 0)
  })
  
  # Gelman-Rubin plot
  plot_base("gelman_loglikelihood.png", {
    coda::gelman.plot(loglikelihood_mcmc, autoburnin = FALSE, ask = FALSE)
    graphics::title(main = paste(site_name, "- Gelman-Rubin for Log-Likelihood"), line = 2.5)
  })
  
  
  # Histogram
  plot_base("loglikelihood_histogram.png", {
    graphics::hist(unlist(clean_chains), breaks = 40, main = paste(site_name, "- Histogram of Log-Likelihood"),
         xlab = "Log-Likelihood", col = "steelblue", border = "white")
  })
  
  # Autocorrelation plots per chain
  plot_base("loglikelihood_acf.png", {
    n_chains <- length(clean_chains)
    plot_dims <- grDevices::n2mfrow(n_chains) 
    graphics::par(mfrow = plot_dims)
    for (i in seq_along(clean_chains)) {
      stats::acf(clean_chains[[i]], main = paste("ACF - Chain", i), lag.max = 50)
    }
  })

  
  # Return summary
  return(list(
    gelman     = gelman_result,
    ess        = ess_result,
    n_valid    = sapply(clean_chains, length),
    range      = sapply(clean_chains, function(x) range(x, na.rm = TRUE))
  ))
}
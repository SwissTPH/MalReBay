generate_likelihood_diagnostics <- function(all_chains_loglikelihood, site_name) {
  # Create a dedicated folder for this site's diagnostics
  site_dir <- file.path("diagnostics", site_name)
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
    graphics::title(main = paste(site_name, "- Gelman-Rubin for "), line = 2.5)
  })
  
  
  # Histogram
  plot_base("loglikelihood_histogram.png", {
    graphics::hist(unlist(clean_chains), breaks = 40, main = paste(site_name, "- Histogram of Log-Likelihood"),
         xlab = "Log-Likelihood", col = "steelblue", border = "white")
  })
  
  # Autocorrelation plots per chain
  plot_base("loglikelihood_acf.png", {
    graphics::par(mfrow = c(2, 2))
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




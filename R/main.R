
# Set your working directory
setwd("D:/OneDrive/Masters Class notes/AIMS Rwanda - Mathematical Sciences - Epidemiology/Thesis Phase/R_Package_Folder/MSP, GLUR, Microsatellite, ampseq") 


# Load all the libraries needed for the analysis
library(rJava)
library(gtools)
library(readxl)
library(coda)
library(future)
library(future.apply)
library(abind)
library(dplyr)
library(kableExtra)

# Set Java memory options manually for the script
options(java.parameters = "-Xmx4096m")

# Load all of your custom R functions from the R/ folder
# This is the "normal R" equivalent of devtools::load_all()
source("R/calculate_frequencies.R") 
source("R/define_alleles.R")
source("R/MalReBay.R")
source("R/define_alleles.R")
source("R/findposterior_frequencies.R") 
source("R/import_data.R")
source("R/mcmc__haplotype_ampseq.R")
source("R/mcmc_alleles_length.R")
source("R/mcmc_convergence_report.R")
source("R/recode_alleles.R")
source("R/run_all_sites.R")
source("R/unobserved_alleles_length.R")
source("R/unobserved_haplotype_ampseq.R")
source("R/match_counting.R")
source("R/diversity_plot.R")

# Input Excel file containing cleaned TES genotype data
input_file <- "Angola_2021_TES_7NMS.xlsx"
output_dir <- "bayesian_algorithm_output"
imported_data <- import_data(filepath = input_file)


# MCMC configuration
mcmc_params <- list(
  n_chains = 4,                # Number of parallel chains
  chunk_size = 6000,           # Iterations per convergence check
  max_iterations = 18000,     # Max iterations before stopping, even if not converge
  record_hidden_alleles = FALSE # Record hidden alleles in the output
)



# Run the classification
plan(multisession, workers = 8)

# Run the full pipeline
final_summary <- classify_infections(
  input_filepath = input_file, 
  mcmc_config        = mcmc_params,
  output_folder      = output_dir
)

plan(sequential)
message("MCMC analysis complete.")

final_summary <- final_summary$summary
marker_details <- final_summary$marker_details
final_summary$Probability <- as.numeric(as.character(final_summary$Probability))

# Histogram of posterior probabilities
pdf_file_path <- file.path(output_dir, "recrudescence_probability_histogram.pdf")

pdf(file = pdf_file_path, width = 8, height = 6)
graphics::hist(final_summary$Probability, 
               breaks = 20,
               main = "Posterior Probability of Recrudescence",
               xlab = "Probability of Recrudescence",
               ylab = "Frequency (Number of Patients)",
               col = "skyblue",
               border = "white")
dev.off()


# markers diversisty plots 
if (imported_data$data_type == "length_polymorphic") {
  
  # Generate the OVERALL plot
  overall_plot_data <- rbind(imported_data$late_failures, imported_data$additional)
  generate_allele_frequency_plot(
    raw_data_df   = overall_plot_data,
    site_name     = NULL, 
    output_folder = output_dir
  )
  
  # Generate the site-specific plots
  site_names <- unique(imported_data$late_failures$Site)
  for (site in site_names) {
    site_data <- imported_data$late_failures[imported_data$late_failures$Site == site, ]
    additional_site_data <- imported_data$additional[imported_data$additional$Site == site, ]
    
    generate_allele_frequency_plot(
      raw_data_df = rbind(site_data, additional_site_data),
      site_name = site,
      output_folder = output_dir
    )
  }
}

# Load the match-counting function
match_results_df <- perform_match_counting(genotypedata_latefailures = imported_data$late_failures, marker_info = imported_data$marker_info)
mcmc_summary_for_merge <- final_summary[, c("Sample.ID", "Probability", "N_Comparable_Loci")]
colnames(mcmc_summary_for_merge) <- c("Sample.ID", "Prob_Recrudescence", "N_Comparable_Loci")
alleles_with_patient_id <- dplyr::mutate(imported_data$late_failures, Patient.ID = gsub(" (Day 0|Day Failure)$", "", Sample.ID))
final_table <- dplyr::left_join(alleles_with_patient_id, match_results_df, by = c("Patient.ID" = "Sample.ID"))
final_table <- dplyr::left_join(final_table, mcmc_summary_for_merge, by = c("Patient.ID" = "Sample.ID"))
final_table_sorted <- dplyr::arrange(final_table, Patient.ID, Sample.ID)
id_cols <- c("Sample.ID", "Site")
original_allele_cols <- grep("(_allele_\\d+|_\\d+)$", colnames(final_table_sorted), value = TRUE)
analysis_cols <- c("Number_Matches", "Number_Loci_Compared", setdiff(colnames(match_results_df), c("Sample.ID", "Number_Matches", "Number_Loci_Compared")),
                   "Prob_Recrudescence", "N_Comparable_Loci")
final_cols_order <- c(id_cols, original_allele_cols, analysis_cols)
final_comparison_table <- dplyr::select(final_table_sorted, all_of(final_cols_order))
final_comparison_table$Patient.ID <- NULL 
output_path <- file.path(output_dir, "algorithm_classification_comparison_table.csv")
write.csv(final_comparison_table, output_path, row.names = FALSE, na = "")

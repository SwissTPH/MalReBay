knitr::opts_chunk$set(
    echo = TRUE,
    collapse = TRUE,
    comment = "#>"
)

library(MalReBay)
library(future)
library(dplyr)

# Input Excel file containing cleaned TES genotype data
input_file <- system.file("extdata", "Angola_2021_TES_7NMS.xlsx", 
                          package = "MalReBay")

imported_data <- MalReBay:::import_data(filepath = input_file)

# MCMC configuration
mcmc_params <- list(
  n_chains = 4,
  chunk_size = 2000,
  max_iterations = 5000,
  burn_in_frac = 0.25,
  record_hidden_alleles = FALSE 
)

# Output folder
output_dir <- file.path(tempdir(), "infection_classification_results")
if (!dir.exists(output_dir)) dir.create(output_dir)

plan(multisession, workers = mcmc_params$n_chains)

# Run the full pipeline
classification_summary <- classify_infections(
  input_filepath = input_file, 
  mcmc_config        = mcmc_params,
  output_folder      = output_dir
)

plan(sequential)


final_summary <- classification_summary$summary
marker_details <- classification_summary$marker_details

# Ensure the Probability column is numeric for plotting
final_summary$Probability <- as.numeric(as.character(final_summary$Probability))

knitr::kable(head(final_summary), caption = "Top rows of the MCMC classification summary.")

# Histogram of posterior probabilities
graphics::hist(final_summary$Probability, breaks = 10,
     main = "Posterior Probability of Recrudescence",
     xlab = "Probability")

# markers diversisty plots 
if (imported_data$data_type == "length_polymorphic") {
  
  # OVERALL plot
  overall_plot_data <- rbind(imported_data$late_failures, imported_data$additional)
  MalReBay:::generate_allele_frequency_plot(
    raw_data_df   = overall_plot_data,
    site_name     = NULL, 
    output_folder = output_dir
  )
  
  # Site-specific plots
  site_names <- unique(imported_data$late_failures$Site)
  for (site in site_names) {
    site_data <- imported_data$late_failures[imported_data$late_failures$Site == site, ]
    additional_site_data <- imported_data$additional[imported_data$additional$Site == site, ]
    
    MalReBay:::generate_allele_frequency_plot(
      raw_data_df = rbind(site_data, additional_site_data),
      site_name = site,
      output_folder = output_dir
    )
  }
}

match_results_df <- MalReBay:::perform_match_counting(
  genotypedata_latefailures = imported_data$late_failures, 
  marker_info = imported_data$marker_info)

mcmc_summary_for_merge <- final_summary[, c("Sample.ID", "Probability", "N_Comparable_Loci")]
colnames(mcmc_summary_for_merge) <- c("Sample.ID", "Prob_Recrudescence", "N_Comparable_Loci")
alleles_with_patient_id <- dplyr::mutate(imported_data$late_failures, Patient.ID = gsub(" (Day 0|Day Failure)$", "", Sample.ID))
final_table <- dplyr::left_join(alleles_with_patient_id, match_results_df, by = c("Patient.ID" = "Sample.ID"))
final_table <- dplyr::left_join(final_table, mcmc_summary_for_merge, by = c("Patient.ID" = "Sample.ID"))
final_table_sorted <- dplyr::arrange(final_table, Patient.ID, Sample.ID)
id_cols <- c("Sample.ID", "Site")
original_allele_cols <- grep("(_allele_\\d+|_\\d+)$", colnames(final_table_sorted), value = TRUE)

analysis_cols <- c("Number_Matches", "Number_Loci_Compared", 
                   setdiff(colnames(match_results_df), 
                           c("Sample.ID", "Number_Matches", "Number_Loci_Compared")),
                   "Prob_Recrudescence", "N_Comparable_Loci")

final_cols_order <- c(id_cols, original_allele_cols, analysis_cols)
final_comparison_table <- dplyr::select(final_table_sorted, all_of(final_cols_order))
final_comparison_table$Patient.ID <- NULL 
output_path <- file.path(output_dir, "algorithm_classification_comparison_table.csv")
write.csv(final_comparison_table, output_path, row.names = FALSE, na = "")

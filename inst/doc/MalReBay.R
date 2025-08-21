knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(MalReBay)
library(future)
library(dplyr)
library(ggplot2)

# This command finds the example file within the installed MalReBay package
input_file <- system.file("extdata", "Angola_2021_TES_7NMS.xlsx",
                          package = "MalReBay")

print(input_file)

# NOTE: These settings are for a quick demonstration only.
quick_mcmc_config <- list(
  n_chains = 2,          # Run two parallel chains
  chunk_size = 1000,       # Check convergence every 1000 iterations
  max_iterations = 1000,   # Stop after a maximum of 1000 total iterations
  rhat_threshold = 1.1,  # A relaxed convergence threshold for the example
  ess_threshold = 50     # A relaxed effective sample size for the example
)

output_dir <- file.path(tempdir(), "MalReBay_tutorial_results")
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

print(output_dir)

plan(multisession, workers = 2)

classification_summary <- classify_infections(
  input_filepath = input_file,
  mcmc_config = quick_mcmc_config,
  output_folder = output_dir
)

plan(sequential)

summary_df <- classification_summary$summary

knitr::kable(head(summary_df), caption = "Top rows of the classification summary.")

marker_details_df <- classification_summary$marker_details

knitr::kable(head(marker_details_df), caption = "Top rows of the marker-level summary.")

ggplot(summary_df, aes(x = Probability)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue", color = "white", boundary = 0) +
  labs(
    title = "Distribution of Posterior Probabilities of Recrudescence",
    x = "Posterior Probability",
    y = "Number of Patients"
  ) +
  theme_bw()

imported_data <- MalReBay:::import_data(filepath = input_file)

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

imported_data <- MalReBay:::import_data(filepath = input_file)

match_results_df <- MalReBay:::perform_match_counting(
  genotypedata_latefailures = imported_data$late_failures, 
  marker_info = imported_data$marker_info)

mcmc_summary_for_merge <- summary_df[, c("Sample.ID", "Probability", "N_Comparable_Loci")]
colnames(mcmc_summary_for_merge) <- c("Sample.ID", "Prob_Recrudescence", "N_Comparable_Loci_MCMC")

alleles_with_patient_id <- dplyr::mutate(imported_data$late_failures, 
                                         Patient.ID = gsub(" (Day 0|Day Failure)$", "", Sample.ID))

final_table <- dplyr::left_join(alleles_with_patient_id, match_results_df, by = c("Patient.ID" = "Sample.ID"))
final_table <- dplyr::left_join(final_table, mcmc_summary_for_merge, by = c("Sample.ID" = "Sample.ID"))

final_table_sorted <- dplyr::arrange(final_table, Patient.ID, Sample.ID)
id_cols <- c("Sample.ID", "Site")
original_allele_cols <- grep("(_allele_\\d+|_\\d+)$", colnames(final_table_sorted), value = TRUE)
analysis_cols <- c("Number_Matches", "Number_Loci_Compared", 
                   setdiff(colnames(match_results_df), c("Sample.ID", "Number_Matches", "Number_Loci_Compared")),
                   "Prob_Recrudescence", "N_Comparable_Loci_MCMC")

final_cols_order <- intersect(c(id_cols, original_allele_cols, analysis_cols), colnames(final_table_sorted))

final_comparison_table <- dplyr::select(final_table_sorted, dplyr::all_of(final_cols_order))

output_path <- file.path(output_dir, "algorithm_classification_comparison_table.csv")
write.csv(final_comparison_table, output_path, row.names = FALSE, na = "")

knitr::kable(
  head(final_comparison_table[, 1:12]), 
  caption = "First 12 columns of the final comparison table."
)

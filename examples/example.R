###########################################
# Example of how to use MalReBay with length
# polymorphic and amplicon seuqencing data
###########################################

library(MalReBay)
library(ggplot2)

##############
# Example 1: length polymorphic data
#############

#----------USER FILE SYSTEM ------------------

# Input data
input_file_micr <- "~/GitRepos/MalReBay/inst/extdata/Angola_2021_TES_7NMS.xlsx"

# Folder where all results will be saved
output_dir_micr <- "~/Bayesian_tests/"

#------------------------------------------

quick_mcmc_config = list(
  n_chains = 2, 
  chunk_size = 100, #5000
  max_iterations = 500, #20000
  rhat_threshold = 1.1,
  ess_threshold = 400
)

mcmc_config = list(
  n_chains = 2, 
  chunk_size = 5000,
  max_iterations = 20000,
  rhat_threshold = 1.1,
  ess_threshold = 400
)

plan(multisession, workers = 2)

classification_summary = classify_infections(
  input_filepath = input_file_micr,
  mcmc_config = quick_mcmc_config,
  output_folder = output_dir_micr,
  verbose = FALSE
)

plan(sequential)

# Plot the distribution of the probabilities of recrudescence
ggplot(classification_summary$summary, aes(x = Probability)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue", color = "white", boundary = 0) +
  labs(
    title = "Distribution of Posterior Probabilities of Recrudescence",
    x = "Posterior Probability",
    y = "Number of Patients"
  ) +
  theme_bw()

##############
# Example 2: amplicon sequencing data
#############


#----------USER FILE SYSTEM ------------------
input_file_amp_seq <- "~/GitRepos/MalReBay/inst/extdata/Amplicon_Sequencing.xlsx"

n_cores <- 2
---------------------------------------------

# Main command of the package, parallel computation of each chain
plan(multisession, workers = n_cores)

classification_summary_ampseq <- classify_infections(
  input_filepath = input_file_amp_seq,
  mcmc_config = quick_mcmc_config,
  output_folder = output_dir,
  verbose = FALSE
)

plan(sequential)

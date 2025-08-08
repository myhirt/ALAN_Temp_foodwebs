#------------ Load all packages ------------
suppressPackageStartupMessages(library(tidyverse))

# Set the folder path
folder_path <- "../../ALAN_temp_fw_original"

# Get list of all CSV files in the folder
csv_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)

# Read and combine all CSV files into one dataframe
combined_df <- csv_files %>%
  lapply(read.csv) %>%
  bind_rows()


write.csv(combined_df, "../output/results_combined.csv")
# -------------------------------------------------------------------------
# loads narrow peak ATAC-seq data, omits sex chromosome data
# original scripts by Adina: filtering.R
# edited by sheryl 06.07.2026
# -------------------------------------------------------------------------

# Load libraries
library(dplyr)

# global variables - check these ------------------------------------------

# Input and output directories
input_dir <- "V:/MARIA/bfx2668_ATAC_E14.5/results/bwa/merged_library/macs2/narrow_peak"

output_dir <- "V:/MARIA/Sheryl/atac_seq_results/E14.5/filtered_noSex"

# Create output dir if it doesn't exist
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# List of files
files <- c(
  "KO_REP1.mLb.clN_peaks.narrowPeak",
  "KO_REP2.mLb.clN_peaks.narrowPeak",
  "KO_REP3.mLb.clN_peaks.narrowPeak",
  "KO_REP4.mLb.clN_peaks.narrowPeak",
  "WT_REP1.mLb.clN_peaks.narrowPeak",
  "WT_REP2.mLb.clN_peaks.narrowPeak",
  "WT_REP3.mLb.clN_peaks.narrowPeak",
  "WT_REP4.mLb.clN_peaks.narrowPeak",
  "WT_REP5.mLb.clN_peaks.narrowPeak",
  "WT_REP6.mLb.clN_peaks.narrowPeak"
)

# -------------------------------------------------------------------------

# Loop through and filter
for (file in files) {
  input_path <- file.path(input_dir, file)
  
  df <- read.table(input_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  
  # Keep only autosomes (chr1–chr19)
  df_noSex <- df %>% filter(grepl("^chr([1-9]|1[0-9])$", V1))
  
  # Create new filename with _noSex
  new_name <- sub("\\.narrowPeak$", "_noSex.narrowPeak", file)
  output_path <- file.path(output_dir, new_name)
  
  # Save
  write.table(df_noSex, output_path, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  print(paste("Finished for", file))
}

cat("Filtered noSex files saved to:", output_dir, "\n")


library(chromVAR)
library(betterChromVAR)
library(BiocParallel)

# Base directory containing dataset folders
data_dir <- "V:/MARIA/bfx2668_ATAC_E14.5/results/bedtools"

# Input folder containing filtered no-sex narrowPeak files
path_filtered <- file.path("V:/MARIA/Sheryl/atac_seq_results/E14.5/filtered_noSex")

wt_files <- list.files(path_filtered, pattern = "WT_.*narrowPeak$", full.names = TRUE)
ko_files <- list.files(path_filtered, pattern = "KO_.*narrowPeak$", full.names = TRUE)

wt_peaks <- sapply(wt_files, readNarrowpeaks)
ko_peaks <- sapply(ko_files, readNarrowpeaks)

wt_beds <- list.files(data_dir, pattern="WT_", full.names=TRUE)
ko_beds <- list.files(data_dir, pattern="KO_", full.names=TRUE)

all_peaks <- c(wt_peaks, ko_peaks)
all_beds <- c(wt_beds, ko_beds)

counts_all <- getCounts(
  alignment_files = all_beds,
  peaks = all_peaks,
  paired = TRUE,
  by_rg = FALSE,
  format = "bed"
)

colnames(counts_all) <- names(bam_files)
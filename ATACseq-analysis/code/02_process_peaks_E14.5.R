############################################################
# MERGED ATAC-seq PEAK ANALYSIS SCRIPT
# WT vs KO common/unique peaks + variable peak filtering
# merges original scripts by Adina: Unique_peaks.R and Variable_coordinates.R
# made by sheryl 06.07.2026 with chatGPT
# "unique peaks" = found in one group and not the other, "1bp" = discard a unique peak if it overlaps with a union peak by >1 bp
# "union peaks" = entire cumulative regions called across all samples
# "variable peaks" = not found in all replicates
# "common peaks" = use granges findoverlaps to find peaks that are shared between all WT samples / KO samples / WT+KO samples
############################################################

library(GenomicRanges)
library(rtracklayer)
library(dplyr)

############################################################
# 1. USER SETTINGS
############################################################

# Base directory containing dataset folders
data_dir <- "V:/MARIA/bfx2668_ATAC_E14.5/results/bwa/merged_library/macs2/narrow_peak"

# Input folder containing filtered no-sex narrowPeak files
path_filtered <- file.path("V:/MARIA/Sheryl/atac_seq_results/E14.5/filtered_noSex")

# Output folders
output_dir <- "V:/MARIA/Sheryl/atac_seq_results/E14.5/processed_peaks"
dir_unique <- file.path(output_dir, "unique_peaks")
dir_variable <- file.path(output_dir, "variable_peaks")
dir_common <- file.path(output_dir, "common_peaks")

# Create output folders if they do not exist
dir.create(dir_unique, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_variable, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_common, showWarnings = FALSE, recursive = TRUE)


############################################################
# 2. FUNCTION TO LOAD AND CLEAN PEAK FILES
############################################################

load_peak_file <- function(file) {
  df <- read.table(file, header = FALSE)
  
  colnames(df)[1:3] <- c("chr", "start", "end")
  
  # Keep only autosomal chromosomes, e.g. chr1–chr19
  df <- df[grepl("^chr[0-9]+$", df$chr), ]
  
  gr <- makeGRangesFromDataFrame(df[, 1:3])
  return(gr)
}


############################################################
# 3. FUNCTION TO CREATE UNION PEAKS
############################################################

make_union_peaks <- function(files) {
  peaks <- lapply(files, load_peak_file)
  
  cat("Number of peaks per replicate:\n")
  print(sapply(peaks, length))
  
  combined <- do.call(c, peaks)
  combined <- as(combined, "GRanges")
  
  union_peaks <- GenomicRanges::reduce(combined)
  
  return(list(
    replicate_peaks = peaks,
    union_peaks = union_peaks
  ))
}


############################################################
# 4. FUNCTION TO FIND VARIABLE PEAKS
############################################################

find_variable_peaks <- function(replicate_peaks, union_peaks) {
  
  peak_matrix <- do.call(cbind, lapply(replicate_peaks, function(gr) {
    as.integer(countOverlaps(union_peaks, gr) > 0)
  }))
  
  peak_overlap_counts <- rowSums(peak_matrix)
  
  # Variable peaks are peaks not found in all replicates
  variable_peaks <- union_peaks[peak_overlap_counts < length(replicate_peaks)]
  
  return(variable_peaks)
}


############################################################
# 5. FUNCTION FOR RECIPROCAL OVERLAP
############################################################

reciprocal_overlap <- function(gr1, gr2, threshold = 0.5) {
  
  hits <- findOverlaps(gr1, gr2)
  
  q <- gr1[queryHits(hits)]
  s <- gr2[subjectHits(hits)]
  
  overlap <- pintersect(q, s)
  ol_width <- width(overlap)
  
  frac_q <- ol_width / width(q)
  frac_s <- ol_width / width(s)
  
  keep <- which(frac_q >= threshold & frac_s >= threshold)
  
  list(
    gr1_hits = queryHits(hits)[keep],
    gr2_hits = subjectHits(hits)[keep]
  )
}


############################################################
# 6. FUNCTION TO EXPORT COMMON AND UNIQUE PEAKS
############################################################

export_common_unique <- function(union_wt, union_ko, threshold, label) {
  
  if (threshold == "1bp") {
    
    hits <- findOverlaps(union_wt, union_ko, minoverlap = 1)
    
    common_wt <- union_wt[queryHits(hits)]
    common_ko <- union_ko[subjectHits(hits)]
    
  } else {
    
    ol <- reciprocal_overlap(union_wt, union_ko, threshold = threshold)
    
    common_wt <- union_wt[ol$gr1_hits]
    common_ko <- union_ko[ol$gr2_hits]
  }
  
  # Create consensus common peaks
  common_peaks <- GenomicRanges::reduce(as(c(common_wt, common_ko), "GRanges"))
  names(common_peaks) <- paste0("peak_", seq_along(common_peaks))
  
  # Save common peaks
  common_file <- file.path(dir_common, paste0("common_peaks_", label, ".bed"))
  export(common_peaks, common_file, format = "BED")
  
  # Unique WT = WT union peaks that do not overlap common peaks
  wt_hits <- findOverlaps(union_wt, common_peaks, minoverlap = 1)
  unique_wt <- union_wt[-unique(queryHits(wt_hits))]
  
  # Unique KO = KO union peaks that do not overlap common peaks
  ko_hits <- findOverlaps(union_ko, common_peaks, minoverlap = 1)
  unique_ko <- union_ko[-unique(queryHits(ko_hits))]
  
  # Save unique peaks
  unique_wt_file <- file.path(dir_unique, paste0("unique_WT_", label, ".bed"))
  unique_ko_file <- file.path(dir_unique, paste0("unique_KO_", label, ".bed"))
  
  export(unique_wt, unique_wt_file, format = "BED")
  export(unique_ko, unique_ko_file, format = "BED")
  
  cat("\nOverlap:", label, "\n")
  cat("Common peaks:", length(common_peaks), "\n")
  cat("Unique WT:", length(unique_wt), "\n")
  cat("Unique KO:", length(unique_ko), "\n")
  
  return(list(
    common = common_peaks,
    unique_wt = unique_wt,
    unique_ko = unique_ko
  ))
}


############################################################
# 7. LOAD WT AND KO FILES
############################################################

wt_files <- list.files(path_filtered, pattern = "WT_.*narrowPeak$", full.names = TRUE)
ko_files <- list.files(path_filtered, pattern = "KO_.*narrowPeak$", full.names = TRUE)

cat("WT files found:\n")
print(wt_files)

cat("\nKO files found:\n")
print(ko_files)

if (length(wt_files) == 0) {
  stop("No WT narrowPeak files found.")
}

if (length(ko_files) == 0) {
  stop("No KO narrowPeak files found.")
}


############################################################
# 8. CREATE WT AND KO UNION PEAKS
############################################################

cat("\nLoading WT peaks...\n")
wt_data <- make_union_peaks(wt_files)

cat("\nLoading KO peaks...\n")
ko_data <- make_union_peaks(ko_files)

wt_peaks <- wt_data$replicate_peaks
ko_peaks <- ko_data$replicate_peaks

union_wt <- wt_data$union_peaks
union_ko <- ko_data$union_peaks

cat("\nUnion WT peaks:", length(union_wt), "\n")
cat("Union KO peaks:", length(union_ko), "\n")


############################################################
# 9. FIND AND SAVE VARIABLE PEAKS
############################################################

variable_wt <- find_variable_peaks(wt_peaks, union_wt)
variable_ko <- find_variable_peaks(ko_peaks, union_ko)

variable_wt_file <- file.path(dir_variable, "WT_variable_peaks.bed")
variable_ko_file <- file.path(dir_variable, "KO_variable_peaks.bed")

export(variable_wt, variable_wt_file, format = "BED")
export(variable_ko, variable_ko_file, format = "BED")

cat("\nVariable WT peaks:", length(variable_wt), "\n")
cat("Saved to:", variable_wt_file, "\n")

cat("\nVariable KO peaks:", length(variable_ko), "\n")
cat("Saved to:", variable_ko_file, "\n")


############################################################
# 10. EXPORT COMMON AND UNIQUE PEAKS
############################################################

# 1 bp overlap
results_1bp <- export_common_unique(
  union_wt = union_wt,
  union_ko = union_ko,
  threshold = "1bp",
  label = "1bp"
)

# 50% reciprocal overlap
results_50 <- export_common_unique(
  union_wt = union_wt,
  union_ko = union_ko,
  threshold = 0.5,
  label = "50percent"
)

# 80% reciprocal overlap
results_80 <- export_common_unique(
  union_wt = union_wt,
  union_ko = union_ko,
  threshold = 0.8,
  label = "80percent"
)

# 100% reciprocal overlap
results_100 <- export_common_unique(
  union_wt = union_wt,
  union_ko = union_ko,
  threshold = 1.0,
  label = "100percent"
)


############################################################
# 11. EXPORT ALL UNION PEAKS
############################################################

all_union <- GenomicRanges::reduce(c(union_wt, union_ko))

all_union_file <- file.path(output_dir, "all_union_peaks.bed")

export(all_union, all_union_file, format = "BED")

cat("\nAll union peaks:", length(all_union), "\n")
cat("Saved to:", all_union_file, "\n")


############################################################
# 12. FILTER UNIQUE PEAKS BY REMOVING VARIABLE PEAKS
############################################################

filter_unique_by_variable <- function(unique_file, variable_peaks, output_file) {
  
  unique_peaks <- import(unique_file, format = "BED")
  
  filtered_unique <- unique_peaks[!unique_peaks %over% variable_peaks]
  
  export(filtered_unique, output_file, format = "BED")
  
  return(list(
    original = unique_peaks,
    filtered = filtered_unique
  ))
}


# Filter only the 1 bp unique peaks, matching your original script
unique_wt_1bp_file <- file.path(dir_unique, "unique_WT_1bp.bed")
unique_ko_1bp_file <- file.path(dir_unique, "unique_KO_1bp.bed")

filtered_unique_wt_1bp_file <- file.path(dir_unique, "filtered_unique_WT_1bp.bed")
filtered_unique_ko_1bp_file <- file.path(dir_unique, "filtered_unique_KO_1bp.bed")

filtered_wt <- filter_unique_by_variable(
  unique_file = unique_wt_1bp_file,
  variable_peaks = variable_wt,
  output_file = filtered_unique_wt_1bp_file
)

filtered_ko <- filter_unique_by_variable(
  unique_file = unique_ko_1bp_file,
  variable_peaks = variable_ko,
  output_file = filtered_unique_ko_1bp_file
)

cat("\nWT unique peaks before filtering:", length(filtered_wt$original), "\n")
cat("WT unique peaks after filtering:", length(filtered_wt$filtered), "\n")

cat("\nKO unique peaks before filtering:", length(filtered_ko$original), "\n")
cat("KO unique peaks after filtering:", length(filtered_ko$filtered), "\n")


############################################################
# 13. OPTIONAL SUMMARY TABLE FOR E14.5 AND E14.5
############################################################

count_peaks <- function(file) {
  gr <- import(file, format = "BED")
  length(gr)
}

conditions <- c("WT", "KO")

summary_results <- data.frame()

for (cond in conditions) {
    
    folder <- file.path(output_dir, "unique_peaks")
    
    orig_file <- file.path(folder, paste0("unique_", cond, "_1bp.bed"))
    filt_file <- file.path(folder, paste0("filtered_unique_", cond, "_1bp.bed"))
    
    if (file.exists(orig_file) && file.exists(filt_file)) {
      
      orig_count <- count_peaks(orig_file)
      filt_count <- count_peaks(filt_file)
      
      summary_results <- rbind(
        summary_results,
        data.frame(
          Condition = cond,
          Original_Count = orig_count,
          Filtered_Count = filt_count,
          Difference = orig_count - filt_count,
          Percent_Retained = round((filt_count / orig_count) * 100, 2)
        )
      )
    }
  }

cat("\nSummary of unique peak filtering:\n")
print(summary_results)
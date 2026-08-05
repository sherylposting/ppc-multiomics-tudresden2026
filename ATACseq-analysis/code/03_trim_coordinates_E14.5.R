# -------------------------------------------------------------------------
# original script by Adina: Trimming_files.R
# edited by sheryl 06.07.2026 to use mm39 instead of mm10
# -------------------------------------------------------------------------

library(GenomicRanges)
library(rtracklayer)
library(GenomeInfoDb)

# global variables - check these ------------------------------------------

# Load mm39 chromosome sizes from UCSC
mm39_sizes <- read.table(
  "V:/MARIA/Sheryl/atac_seq_results/mm39.chrom.sizes.txt",
  col.names = c("chr", "size")
)

base_dir <- "V:/MARIA/Sheryl/atac_seq_results/E14.5/processed_peaks/unique_peaks"

# Load the BED files
KO_peaks <- import(file.path(base_dir, "filtered_unique_KO_1bp.bed"))
WT_peaks <- import(file.path(base_dir, "filtered_unique_WT_1bp.bed"))

# -------------------------------------------------------------------------

# Keep only peaks within valid chr sizes
valid_sizes <- setNames(mm39_sizes$size, mm39_sizes$chr)
KO_peaks <- KO_peaks[ end(KO_peaks) <= valid_sizes[as.character(seqnames(KO_peaks))] ]
WT_peaks <- WT_peaks[ end(WT_peaks) <= valid_sizes[as.character(seqnames(WT_peaks))] ]

# Save filtered BED
# export(KO_peaks, file.path(base_dir, "unique_KO_1bp_validrange.bed"), format="BED")
# export(WT_peaks, file.path(base_dir, "unique_WT_1bp_validrange.bed"), format="BED")

# -------------------------------------------------------------------------

mm39_seqinfo <- Seqinfo(seqnames = mm39_sizes$chr, seqlengths = mm39_sizes$size)

# Match seqlevels (chromosome names) with mm39
common_chr <- intersect(seqlevels(KO_peaks), seqlevels(mm39_seqinfo))
mm39_seqinfo_clean <- mm39_seqinfo[common_chr]
seqinfo(KO_peaks) <- mm39_seqinfo_clean

common_chr <- intersect(seqlevels(WT_peaks), seqlevels(mm39_seqinfo))
mm39_seqinfo_clean <- mm39_seqinfo[common_chr]
seqinfo(WT_peaks) <- mm39_seqinfo_clean

# Trim peaks to fit chromosome limits
KO_peaks_trimmed <- trim(KO_peaks)
WT_peaks_trimmed <- trim(WT_peaks)

# Export to a new trimmed BED file
export(KO_peaks_trimmed,
       file.path(base_dir, "unique_KO_1bp_trimmed.bed"),
       format = "BED")

export(WT_peaks_trimmed,
       file.path(base_dir, "unique_WT_1bp_trimmed.bed"),
       format = "BED")
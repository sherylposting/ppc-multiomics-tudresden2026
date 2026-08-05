# -------------------------------------------------------------------------
# original script by Adina: diffbind.R
# edited by sheryl 06.07.2026 to use mm39 instead of mm10
# input raw BAM and narrowPeak files into diffbind, subset diffpeaks in wt, ko, and both
# -------------------------------------------------------------------------

# === Load all required libraries ===
library(DiffBind)
library(BiocParallel)
library(org.Mm.eg.db)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)
library(dplyr)
library(ggplot2)
library(ggvenn)
library(readxl)

# Enable parallel processing
#register(MulticoreParam(8))  
SnowParam(8) # for windows

# global variables - check these ------------------------------------------

# needed to download files locally due to weird read permissions problems
bam_dir_E14.5 <- "C:/Users/shery/Documents/bioinfo_windows/ppc_multiomics_tudresden2026/ATACseq-analysis/bam_E14.5"

peak_dir_E14.5 <- "V:/MARIA/bfx2668_ATAC_E14.5/results/bwa/merged_library/macs2/narrow_peak"

sample_ids_E14.5 <- c("WT_REP1", "WT_REP2", "WT_REP3", "WT_REP4", "WT_REP5", "WT_REP6",
                      "KO_REP1", "KO_REP2", "KO_REP3", "KO_REP4")

samples_E14.5 <- data.frame(
  SampleID = sample_ids_E14.5,
  Condition = ifelse(grepl("KO", sample_ids_E14.5), "KO", "WT"),
  Replicate = c(1:4, 1:6), # CHECK THIS!
  bamReads = file.path(bam_dir_E14.5, paste0(sample_ids_E14.5, ".mLb.clN.sorted.bam")),
  Peaks = file.path(peak_dir_E14.5, paste0(sample_ids_E14.5, ".mLb.clN_peaks.narrowPeak")),
  PeakCaller = "narrow"
)

output_dir <- "V:/MARIA/Sheryl/atac_seq_results/E14.5/diffbind"

# -------------------------------------------------------------------------

# Create output folders if they do not exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# === DiffBind differential accessibility analysis ===
dba_object_E14.5 <- dba(sampleSheet = samples_E14.5)
dba_object_E14.5 <- dba.count(dba_object_E14.5)
dba_object_E14.5 <- dba.contrast(dba_object_E14.5, categories = DBA_CONDITION, minMembers = 2, reorderMeta=list(Condition="WT")) # explicitly set WT as baseline/reference level
dba_object_E14.5 <- dba.analyze(dba_object_E14.5)
diff_peaks_E14.5 <- dba.report(dba_object_E14.5)

diff_peaks_E14.5 <- diff_peaks_E14.5[!seqnames(diff_peaks_E14.5) %in% c("chrX", "chrY")]
diff_peaks_df_E14.5 <- as.data.frame(diff_peaks_E14.5)

# Save differential peaks 
write.table(diff_peaks_df_E14.5, file = file.path(output_dir, "DiffBind_differential_peaks.mm39.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# PCA plot
pdf(file.path(output_dir, "DiffBind_PCA_plot.E14.5.pdf"), width = 10, height = 8)
dba.plotPCA(dba_object_E14.5, DBA_CONDITION, label = DBA_ID)
dev.off()

# MA plot
pdf(file.path(output_dir, "DiffBind_MA_plot.E14.5.pdf"))
dba.plotMA(dba_object_E14.5)
dev.off()

# Heatmap
pdf(file.path(output_dir, "DiffBind_Heatmap.E14.5.pdf"))
dba.plotHeatmap(dba_object_E14.5)
dev.off()

save(dba_object_E14.5, diff_peaks_E14.5, file = file.path(output_dir, "dba_object_E14.5.RData"))

# wt vs ko separate annotation -------------------------------------------------------

# === Annotate differential peaks using ChIPseeker ===
txdb <- TxDb.Mmusculus.UCSC.mm39.knownGene

# KO-specific peaks (up in KO, high confidence and effect size)
ko_up_peaks_E14.5 <- diff_peaks_E14.5[diff_peaks_E14.5$Fold > 0 & diff_peaks_E14.5$FDR < 0.05]
wt_up_peaks_E14.5 <- diff_peaks_E14.5[diff_peaks_E14.5$Fold < 0 & diff_peaks_E14.5$FDR < 0.05]

wt_annot_E14.5 <- annotatePeak(wt_up_peaks_E14.5, TxDb = txdb, tssRegion = c(-3000, 3000), annoDb = "org.Mm.eg.db")
ko_annot_E14.5 <- annotatePeak(ko_up_peaks_E14.5, TxDb = txdb, tssRegion = c(-3000, 3000), annoDb = "org.Mm.eg.db")

wt_df_E14.5 <- as.data.frame(wt_annot_E14.5)
ko_df_E14.5 <- as.data.frame(ko_annot_E14.5)

wt_df_E14.5 <- wt_df_E14.5[!wt_df_E14.5$seqnames %in% c("chrX", "chrY"), ]
ko_df_E14.5 <- ko_df_E14.5[!ko_df_E14.5$seqnames %in% c("chrX", "chrY"), ]

# === Simplify annotation categories ===
simplify_annotation <- function(annotation) {
  dplyr::case_when(
    grepl("Promoter|promoter|TSS", annotation) ~ "Promoter",
    grepl("5' UTR", annotation) ~ "5' UTR",
    grepl("3' UTR", annotation) ~ "3' UTR",
    grepl("exon|CDS", annotation) ~ "Exon",
    grepl("intron", annotation) ~ "Intron",
    grepl("Downstream", annotation) ~ "Downstream",
    grepl("Intergenic", annotation) ~ "Intergenic",
    TRUE ~ "Other"
  )
}

wt_df_E14.5$Simplified_Annotation <- simplify_annotation(wt_df_E14.5$annotation)
ko_df_E14.5$Simplified_Annotation <- simplify_annotation(ko_df_E14.5$annotation)

# === Plot genomic distribution ===
wt_dist_E14.5 <- wt_df_E14.5 %>%
  count(Simplified_Annotation) %>%
  mutate(Percent = 100 * n / sum(n), Group = "WT_up")

ko_dist_E14.5 <- ko_df_E14.5 %>%
  count(Simplified_Annotation) %>%
  mutate(Percent = 100 * n / sum(n), Group = "KO_up")

combined_dist_E14.5 <- bind_rows(wt_dist_E14.5, ko_dist_E14.5)

pdf(file.path(output_dir, "genomic_distribution_WT_KO_diffbind.E14.5.pdf"), width = 10, height = 6)
ggplot(combined_dist_E14.5, aes(x = Simplified_Annotation, y = Percent, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = paste0(round(Percent, 1), "%")),
            position = position_dodge(width = 0.9),
            vjust = -0.5, size = 3) +
  labs(title = "Genomic Distribution of Differential Peaks (E14.5 WT-Up vs KO-Up)",
       y = "Percentage", x = "Genomic Feature") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# === Identify promoter genes from ATAC ===
wt_promoter_genes_E14.5 <- wt_df_E14.5 %>%
  filter(grepl("Promoter", annotation)) %>%
  pull(geneId) %>%
  unique()

ko_promoter_genes_E14.5 <- ko_df_E14.5 %>%
  filter(grepl("Promoter", annotation)) %>%
  pull(geneId) %>%
  unique()

wt_promoter_genes_symbols_E14.5 <- AnnotationDbi::mapIds(org.Mm.eg.db, keys = wt_promoter_genes_E14.5,
                                                       column = "SYMBOL", keytype = "ENTREZID", multiVals = "first") %>%
  na.omit()

ko_promoter_genes_symbols_E14.5 <- AnnotationDbi::mapIds(org.Mm.eg.db, keys = ko_promoter_genes_E14.5,
                                                       column = "SYMBOL", keytype = "ENTREZID", multiVals = "first") %>%
  na.omit()

# === Venn diagram of promoter gene overlap ===
gene_sets_E14.5 <- list(
  WT_up_promoter_genes = unname(wt_promoter_genes_symbols_E14.5),
  KO_up_promoter_genes = unname(ko_promoter_genes_symbols_E14.5)
)

pdf(file.path(output_dir, "Venn_promoter_genes_WT_KO.E14.5.pdf"), width = 6, height = 6)
ggvenn(gene_sets_E14.5,
       fill_color = c("#E69F00", "#56B4E9"),
       stroke_size = 0.5,
       set_name_size = 4,
       text_size = 4)
dev.off()

wt_promoter_peaks_E14.5 <- wt_df_E14.5 %>%
  filter(grepl("Promoter", annotation))

ko_promoter_peaks_E14.5 <- ko_df_E14.5 %>%
  filter(grepl("Promoter", annotation))

# === Save gene lists ===
write.table(wt_promoter_genes_symbols_E14.5, file = file.path(output_dir, "WT_up_promoter_gene_symbols.mm39.tsv"), sep = "\t", quote = FALSE, col.names = NA)
write.table(ko_promoter_genes_symbols_E14.5, file = file.path(output_dir, "KO_up_promoter_gene_symbols.mm39.tsv"), sep = "\t", quote = FALSE, col.names = NA)

write.table(intersect(wt_promoter_genes_symbols_E14.5, ko_promoter_genes_symbols_E14.5), file = file.path(output_dir, "Intersect_promoter_genes.mm39.tsv"), sep = "\t", quote = FALSE, col.names = NA)
write.table(setdiff(wt_promoter_genes_symbols_E14.5, ko_promoter_genes_symbols_E14.5), file = file.path(output_dir, "Unique_WT_promoter_genes.mm39.tsv"), sep = "\t", quote = FALSE, col.names = NA)
write.table(setdiff(ko_promoter_genes_symbols_E14.5, wt_promoter_genes_symbols_E14.5), file = file.path(output_dir, "Unique_KO_promoter_genes.mm39.tsv"), sep = "\t", quote = FALSE, col.names = NA)


write.table(
  wt_promoter_peaks_E14.5[, c("seqnames", "start", "end")],
  file = file.path(output_dir, "WT_promoter_peaks_E14.5.mm39.bed"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)

write.table(
  ko_promoter_peaks_E14.5[, c("seqnames", "start", "end")],
  file = file.path(output_dir, "KO_promoter_peaks_E14.5.mm39.bed"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)

# wt and ko shared annotation ----------------------------------------------------

diff_peaks_E14.5 <- diff_peaks_E14.5[diff_peaks_E14.5$FDR < 0.05]
diff_annot_E14.5 <- annotatePeak(diff_peaks_E14.5, TxDb = txdb, tssRegion = c(-3000, 3000), annoDb = "org.Mm.eg.db")
diff_annot_df_E14.5 <- as.data.frame(diff_annot_E14.5)
diff_annot_df_E14.5$Simplified_Annotation <- simplify_annotation(diff_annot_df_E14.5$annotation)
write.table(diff_annot_df_E14.5, file.path(output_dir, "all_ATAC_diffgenes.E14.5.tsv"))

save(diff_annot_df_E14.5, wt_df_E14.5, ko_df_E14.5, file = "ATAC_dfs.RData")

# -------------------------------------------------------------------------

# # === RNA-ATAC integration ===
rna_path_E14.5 <- "V:/MARIA/bfx2557_RNA_seq_E14.5_Katrin_no_trimming/bfx2557_all raw files/de-analysis_three_samples_each/data/bfx2557.deseq-results.separate.de-expl.xlsx"

rna_df_E14.5 <- read_excel(rna_path_E14.5, sheet = "3_Cond_knock_out_v_wild_type_F5")
rna_df_E14.5 <- rna_df_E14.5[!rna_df_E14.5$Chromosome %in% c("X", "Y"), ]

sig_rna_df_E14.5 <- rna_df_E14.5[rna_df_E14.5$padj < 0.05, ]
sig_rna_df_E14.5 <- sig_rna_df_E14.5[, c("Ensembl_ID", "Gene_Symbol", "log2FoldChange", "padj")]

wt_merged_E14.5 <- merge(wt_df_E14.5, sig_rna_df_E14.5, by.x = "ENSEMBL", by.y = "Ensembl_ID")
ko_merged_E14.5 <- merge(ko_df_E14.5, sig_rna_df_E14.5, by.x = "ENSEMBL", by.y = "Ensembl_ID")
diff_merged_E14.5 <- merge(diff_annot_df_E14.5, sig_rna_df_E14.5, by.x = "ENSEMBL", by.y = "Ensembl_ID")

write.table(wt_merged_E14.5, "V:/MARIA/Sheryl/atac_seq_results/E14.5/narrow_peak/diffbind/WT_up_overlap_with_RNAseq.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ko_merged_E14.5, "V:/MARIA/Sheryl/atac_seq_results/E14.5/narrow_peak/diffbind/KO_up_overlap_with_RNAseq.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(diff_merged_E14.5, "V:/MARIA/Sheryl/atac_seq_results/E14.5/narrow_peak/diffbind/all_ATAC_overlap_with_RNAseq.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# prepare subset of promoters-only
diff_merged_E14.5_promoters <- diff_merged_E14.5[diff_merged_E14.5$Simplified_Annotation == "Promoter",]

# integrated shared hits scatter plots ------------------------------------

# define helper function
integrated_plotter <- function(diff_merged_df, padj_label, title){
  if(padj_label == "ATAC-seq"){PADJ = "FDR"; COLOR="red"} # color-code according to either ATAC or RNA seq padj values
  if(padj_label == "RNA-seq"){PADJ = "padj"; COLOR="dodgerblue"}
  
  integrated_plot <- ggplot(diff_merged_df, aes(x=Fold, y=log2FoldChange, col=.data[[PADJ]])) +
    geom_point() +
    geom_text_repel(aes(label = SYMBOL), 
                    # data = subset(diff_merged_df, diff_merged_df[[PADJ]] < 0.03), # preferentially label higher padj
                    max.overlaps = 12 # fit more labels
                    #size = 3.2 # text a bit smaller
                    ) + 
    labs(
      x = "ATAC-seq log2fc",
      y = "RNA-seq log2fc",
      title = title,
      color=paste0("padj value (", padj_label, ")")
    ) +
    geom_hline(yintercept=0, linetype="dashed") +
    geom_vline(xintercept=0, linetype="dashed") +
    scale_colour_gradient(low = COLOR, high = "black")
  
  print(integrated_plot)
}

# plot and save ---------------------------------------------------

# full merged data
pdf(file.path(output_dir, "integrated_scatterplot.E14.5.pdf"), width = 10, height = 8)
integrated_plotter(diff_merged_E14.5, "RNA-seq", "E14.5 KO vs WT shared hits between ATAC-seq and RNA-seq (FDR<0.05)")
integrated_plotter(diff_merged_E14.5, "ATAC-seq", "E14.5 KO vs WT shared hits between ATAC-seq and RNA-seq (FDR<0.05)")
dev.off()

# promoters-only
pdf(file.path(output_dir, "integrated_scatterplot_promoters.E14.5.pdf"), width = 10, height = 8)
integrated_plotter(diff_merged_E14.5_promoters, "RNA-seq", "E14.5 KO vs WT shared hits between ATAC-seq and RNA-seq, promoters only (FDR<0.05)")
integrated_plotter(diff_merged_E14.5_promoters, "ATAC-seq", "E14.5 KO vs WT shared hits between ATAC-seq and RNA-seq, promoters only (FDR<0.05)")
dev.off()

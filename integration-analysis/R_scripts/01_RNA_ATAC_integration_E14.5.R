# -------------------------------------------------------------------------
# take diffbind results as well as previous RNA DEseq2 results, subset and plot the genes which are differential in both results. also consider the subset of diffpeaks in promoters only
# -------------------------------------------------------------------------

library(ggplot2)
library(ggrepel)
library(ggvenn)
library(readxl)

# global variables - check these ------------------------------------------

OUTPUT_DIR <- "V:/MARIA/Sheryl/integration_results/RNA_ATAC_integration_plots"

dir.create(OUTPUT_DIR)

# load RNA dfs
rna_path_E14.5 <- "V:/MARIA/bfx2557_RNA_seq_E14.5_Katrin_no_trimming/bfx2557_all raw files/de-analysis_three_samples_each/data/bfx2557.deseq-results.separate.de-expl.xlsx"
# select the sheet for 0.05 FDR
rna_df_E14.5 <- read_excel(rna_path_E14.5, sheet = "3_Cond_knock_out_v_wild_type_F5")

# load ATAC dfs
load("V:/MARIA/Sheryl/atac_seq_results/E14.5/R_scripts/ATAC_dfs.RData") # -> ATAC dfs for wt, ko, both, from diffbind.R

# merge RNA and ATAC dfs --------------------------------------------------

rna_df_E14.5 <- rna_df_E14.5[!rna_df_E14.5$Chromosome %in% c("X", "Y"), ]

sig_rna_df_E14.5 <- rna_df_E14.5[rna_df_E14.5$padj < 0.05, ]
sig_rna_df_E14.5 <- sig_rna_df_E14.5[, c("Ensembl_ID", "Gene_Symbol", "log2FoldChange", "padj")]

# merge on ENSEMBL ID
wt_merged_E14.5 <- merge(wt_df_E14.5, sig_rna_df_E14.5, by.x = "ENSEMBL", by.y = "Ensembl_ID")
ko_merged_E14.5 <- merge(ko_df_E14.5, sig_rna_df_E14.5, by.x = "ENSEMBL", by.y = "Ensembl_ID")
diff_merged_E14.5 <- merge(diff_annot_df_E14.5, sig_rna_df_E14.5, by.x = "ENSEMBL", by.y = "Ensembl_ID")

write.table(wt_merged_E14.5, file.path(OUTPUT_DIR, "WT_up_overlap_with_RNAseq.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ko_merged_E14.5, file.path(OUTPUT_DIR, "KO_up_overlap_with_RNAseq.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(diff_merged_E14.5, file.path(OUTPUT_DIR, "all_ATAC_overlap_with_RNAseq.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

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
pdf(file.path(OUTPUT_DIR, "integrated_scatterplot.E14.5.pdf"), width = 10, height = 8)

rna_atac_all_1 <- integrated_plotter(diff_merged_E14.5, "RNA-seq", "E14.5 KO vs WT shared hits between ATAC-seq and RNA-seq (FDR<0.05)")
rna_atac_all_2 <- integrated_plotter(diff_merged_E14.5, "ATAC-seq", "E14.5 KO vs WT shared hits between ATAC-seq and RNA-seq (FDR<0.05)")

dev.off()

# promoters-only
pdf(file.path(OUTPUT_DIR, "integrated_scatterplot_promoters.E14.5.pdf"), width = 10, height = 8)

rna_atac_promoters_1 <- integrated_plotter(diff_merged_E14.5_promoters, "RNA-seq", "E14.5 KO vs WT shared hits between ATAC-seq and RNA-seq, promoters only (FDR<0.05)")
rna_atac_promoters_2 <- integrated_plotter(diff_merged_E14.5_promoters, "ATAC-seq", "E14.5 KO vs WT shared hits between ATAC-seq and RNA-seq, promoters only (FDR<0.05)")

dev.off()

# save RData for shiny app
save(rna_atac_all_1, rna_atac_all_2, rna_atac_promoters_1, rna_atac_promoters_2, file = file.path(OUTPUT_DIR, "rna_atac_integration_plots.RData"))
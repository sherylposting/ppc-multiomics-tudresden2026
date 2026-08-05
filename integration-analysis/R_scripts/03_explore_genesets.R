# -------------------------------------------------------------------------
# creates volcanoplots of RNA and ATAC data, marks genes based on whether they are present in custom genesets
# has only been run on E14.5
# -------------------------------------------------------------------------

library(readxl)
library(ggplot2)
library(ggrepel)
library(dplyr)

# global variables - check these ------------------------------------------
RNA_PATH_E14.5 <- "V:/MARIA/bfx2557_RNA_seq_E14.5_Katrin_no_trimming/bfx2557_all raw files/de-analysis_three_samples_each/data/bfx2557.deseq-results.separate.de-expl.xlsx"

ATAC_PATH_E14.5 <- "V:/MARIA/Sheryl/atac_seq_results/E14.5/diffbind/all_ATAC_diffgenes.E14.5.tsv"

OUTPUT_DIR_E14.5 <- "V:/MARIA/Sheryl/atac_seq_results/E14.5/explore_genesets"

# load RNA seq results
rna_df_E14.5 <- read_excel(RNA_PATH_E14.5, sheet = "3_Cond_knock_out_v_wild_type_F5")
rna_df_E14.5 <- rna_df_E14.5[!rna_df_E14.5$Chromosome %in% c("X", "Y"), ]
rna_df_E14.5 <- rna_df_E14.5[rna_df_E14.5$padj < 0.05, ]

# load ATAC seq results
atac_df_E14.5 <- read.table(ATAC_PATH_E14.5)

# -------------------------------------------------------------------------

# rename ATAC df columns to match RNA
colnames(atac_df_E14.5)[colnames(atac_df_E14.5) == "Fold"] <- "log2FoldChange"
colnames(atac_df_E14.5)[colnames(atac_df_E14.5) == "FDR"] <- "padj"
colnames(atac_df_E14.5)[colnames(atac_df_E14.5) == "SYMBOL"] <- "Gene_Symbol"

# custom genesets ---------------------------------------------------------

# from https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/ZHOU_PANCREATIC_ENDOCRINE_PROGENITOR.html
geneset_zhou_endocrine <- c("Arx","Isl1","Mafa","Mafb","Mlxipl","Myt1","Neurod1","Neurog3","Nkx2-2","Nkx6-1","Pax4","Pax6","Pou3f4","Vdr")

# from https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/ZHOU_PANCREATIC_BETA_CELL.html
geneset_zhou_beta <- c("Foxa1","Foxo1","Hnf1a","Hnf4a","Isl1","Mafa","Neurod1","Nkx2-2","Nkx6-1","Pax6","Pdx1")


geneset_zhou_exocrine <- c("Foxa2","Hhex","Hnf1b","Hnf4a","Mnx1","Nr5a2","Onecut1","Pdx1","Prox1","Ptf1a","Sox9")

# from https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/HALLMARK_PANCREAS_BETA_CELLS.html
geneset_hallmark <- c("Abcc8","Akt3","Chga","Dcx","Dpp4","Elp4","Foxa2","Foxo1","G6pc2","Gcg","Gck","Hnf1a","Iapp","Ins2","Insm1","Isl1","Lmo2","Mafb","Neurod1","Neurog3","Nkx2-2","Nkx6-1","Pak3","Pax4","Pax6","Pcsk1","Pcsk2","Pdx1","Pklr","Scgn","Sec11a","Slc2a2","Spcs1","Srp14","Srprb","Sst","Stxbp1","Syt13","Vdr")


# check if hits are in custom genesets ------------------------------------

atac_df_E14.5 <- atac_df_E14.5 %>%
  mutate(geneset = case_when(
    Gene_Symbol %in% geneset_zhou_beta ~ "zhou_beta_cell",
    Gene_Symbol %in% geneset_zhou_endocrine ~ "zhou_endocrine",
    Gene_Symbol %in% geneset_hallmark ~ "hallmark",
    Gene_Symbol %in% geneset_zhou_exocrine ~ "zhou_exocrine"
  ))

rna_df_E14.5 <- rna_df_E14.5 %>%
  mutate(geneset = case_when(
    Gene_Symbol %in% geneset_zhou_beta ~ "zhou_beta_cell",
    Gene_Symbol %in% geneset_zhou_endocrine ~ "zhou_endocrine",
    Gene_Symbol %in% geneset_hallmark ~ "hallmark",
    Gene_Symbol %in% geneset_zhou_exocrine ~ "zhou_exocrine"
  ))

# one-by-one t/f labeling
# atac_df_E14.5 %>%
#   mutate(
#     zhou = SYMBOL %in% geneset_zhou_exocrine,
#     hallmark = SYMBOL %in% geneset_hallmark,
#   )

# volcano plot geneset labels ---------------------------------------------

# helper function to plot
volcanoplotter_geneset <- function(DF, TITLE, XLIM = NULL, YLIM = NULL) {
  p <- ggplot(data=DF, aes(x=log2FoldChange, y=-log10(padj), col=geneset))+
    geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
    geom_point(
      data = subset(DF, is.na(geneset)),
      colour = "gray80",
      size = 2
    ) +
    geom_point(
      data = subset(DF, !is.na(geneset)),
      aes(colour = geneset),
      size = 2
    ) +
    geom_text_repel(
      data = subset(DF, !is.na(geneset)),
      aes(label = Gene_Symbol, colour = geneset),
      max.overlaps = 37,
      show.legend = FALSE
    ) +
    labs(
      y = "-log10(padj)",
      x = "log2 fold change",
      title = TITLE,
    )
  
  if (!is.null(XLIM))
    p <- p + xlim(XLIM)
  
  if (!is.null(YLIM))
    p <- p + ylim(YLIM)
  
  print(p)
}

# volcano plot gene symbol labels (regular) -------------------------------

# helper function to plot
volcanoplotter_regular <- function(DF, TITLE, XLIM = NULL, YLIM = NULL, ...) {
  p <- ggplot(data=DF, aes(x=log2FoldChange, y=-log10(padj), col=padj))+
    geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
    geom_point() +
    geom_text_repel(
      data = subset(DF, DF$padj < 0.03), # preferentially label higher padj
      aes(label = Gene_Symbol),
      show.legend = FALSE,
      ... # extra options here just because for some reason the labels look crazy on the plot
    ) +
    labs(
      y = "-log10(padj)",
      x = "log2 fold change",
      title = TITLE,
    ) +
    scale_colour_gradient(low = "dodgerblue", high = "black")
  
  if (!is.null(XLIM))
    p <- p + xlim(XLIM)
  
  if (!is.null(YLIM))
    p <- p + ylim(YLIM)
  
  print(p)
}

# print to pdf ------------------------------------------------------------
pdf(file.path(OUTPUT_DIR_E14.5, "explore_genesets_volcanoplots.E14.5.pdf"), width=10, height=8)

# make geneset plots
volcanoplotter_geneset(atac_df_E14.5, "ATAC-seq hits in E14.5, custom genesets", c(-4,4))
volcanoplotter_geneset(atac_df_E14.5, "ATAC-seq hits in E14.5, custom genesets (zoom)", c(-1,-0.25), c(0, 7.6))
volcanoplotter_geneset(rna_df_E14.5, "RNA-seq hits in E14.5, custom genesets", c(-4,4))

# make regular plots
volcanoplotter_regular(atac_df_E14.5, "ATAC-seq hits in E14.5", c(-4,4), direction="x")
volcanoplotter_regular(atac_df_E14.5, "ATAC-seq hits in E14.5 (zoom)", c(-1,-0.25), c(0, 7.6))
volcanoplotter_regular(rna_df_E14.5, "RNA-seq hits in E14.5", c(-4,4))

dev.off()

# atac promoters only ----------------------------------------------------------
pdf(file.path(OUTPUT_DIR_E14.5, "explore_genesets_volcanoplots_promoters.E14.5.pdf"), width=10, height=8)

atac_E14.5_promoters <- atac_df_E14.5[atac_df_E14.5$Simplified_Annotation == "Promoter",]

# make geneset plots
volcanoplotter_geneset(atac_E14.5_promoters, "ATAC-seq hits in E14.5, promoters only, custom genesets", c(-4,4))
volcanoplotter_geneset(atac_E14.5_promoters, "ATAC-seq hits in E14.5, promoters only, custom genesets (zoom)", c(-1,-0.25), c(0, 7.6))

# make regular plots
volcanoplotter_regular(atac_E14.5_promoters, "ATAC-seq hits in E14.5, promoters only", c(-4,4), direction="x")
volcanoplotter_regular(atac_E14.5_promoters, "ATAC-seq hits in E14.5, promoters only (zoom)", c(-1,-0.25), c(0, 7.6))

dev.off()
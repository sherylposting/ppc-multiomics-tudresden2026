library(decoupleR)
library(dorothea)
library(dplyr)
library(tidyr)
library(pheatmap)
library(OmnipathR)
library(readxl)

# global variables - check these ------------------------------------------
RNA_PATH_E14.5 <- "V:/MARIA/bfx2557_RNA_seq_E14.5_Katrin_no_trimming/bfx2557_all raw files/de-analysis_three_samples_each/data/bfx2557.deseq-results.separate.de-expl.xlsx"

OUTPUT_DIR <- "V:/MARIA/Sheryl/integration_results/RNA_decoupleR_results"

MINSIZE <- 5

dir.create(OUTPUT_DIR)

# load RNA seq results
rna_df_E14.5 <- read_excel(RNA_PATH_E14.5, sheet = "3_Cond_knock_out_v_wild_type_F5")
rna_df_E14.5 <- rna_df_E14.5[!rna_df_E14.5$Chromosome %in% c("X", "Y"), ]

# prepare dorothea network ------------------------------------------------

# omnipath / collectri doesn't work, maybe need new version of R
data("dorothea_mm", package = "dorothea")

# filter higher confidence targets
dorothea_net <- dorothea_mm %>%
  filter(confidence %in% c("A", "B", "C")) %>%
  select(tf, target, mor)


# prepare rna df to decoupleR matrix --------------------------------------

rna_df_E14.5 <- rna_df_E14.5 %>%
  filter(
    !is.na(stat),
    !is.na(Gene_Symbol))

# mat for all genes
mat_all <- rna_df_E14.5 %>%
  select(stat) %>%
  as.matrix()

rownames(mat_all) <- rna_df_E14.5$Gene_Symbol

# mat for DE genes only
rna_df_diff <- rna_df_E14.5[rna_df_E14.5$padj < 0.05, ] %>%
  select(stat, Gene_Symbol) %>%
  drop_na()

mat_DE <- rna_df_diff %>%
  select(stat) %>%
  as.matrix()

rownames(mat_DE) <- rna_df_diff$Gene_Symbol

# mat for DE genes only, all samples
rna_df_diff <- rna_df_E14.5[rna_df_E14.5$padj < 0.05, ] %>%
  select(Gene_Symbol, contains(c("KO", "WT"))) %>%
  drop_na()

mat_DE_samples <- rna_df_diff %>%
  select(
    contains(c("KO", "WT"))
  ) %>%
  as.matrix()

rownames(mat_DE_samples) <- rna_df_diff$Gene_Symbol

# separate df, remove duplicate entries
rna_unique <- rna_df_E14.5 %>%
  filter(
    !is.na(Gene_Symbol),
    Gene_Symbol != "",
    !is.na(stat)
  ) %>%
  group_by(Gene_Symbol) %>%
  slice_max(
    order_by = abs(stat),
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

# run ULM -------------------------------------------

ulm_all <- decoupleR::run_ulm(mat = mat_all, 
                              network = dorothea_net, 
                              .source = 'tf', 
                              .target = 'target', 
                              .mor = 'mor', 
                              minsize = MINSIZE)

ulm_DE <- decoupleR::run_ulm(mat = mat_DE, 
                              network = dorothea_net, 
                              .source = 'tf', 
                              .target = 'target', 
                              .mor = 'mor', 
                              minsize = MINSIZE)

# run ULM and plot for DE genes, all individual samples -------------------------------------

ulm_DE_samples <- decoupleR::run_ulm(mat = mat_DE_samples, 
                              network = dorothea_net, 
                              .source = 'tf', 
                              .target = 'target', 
                              .mor = 'mor', 
                              minsize = MINSIZE)

# Transform to matrix
mat_ulm <- ulm_DE_samples %>%
  dplyr::filter(statistic=='ulm') %>%
  decoupleR::pivot_wider_profile(id_cols = source, 
                                 names_from = condition, 
                                 values_from = score) %>%
  as.matrix()


# ulm hits as sorted tibble -----------------------------------------------

# add FDR
ulm_all$padj <- p.adjust(ulm_all$p_value, method="BH")
ulm_DE$padj <- p.adjust(ulm_DE$p_value, method="BH")

ulm_all <- arrange(ulm_all, padj)
ulm_DE <- arrange(ulm_DE, padj)

ulm_all
ulm_DE


# a volcano plot?? --------------------------------------------------------

# helper function for plotting
volcanoplotter <- function(DF, TITLE) {
  p <- ggplot(data=DF, aes(x=score, y=-log10(padj), col=padj)) +
    geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
    geom_point() +
    geom_text_repel(
      aes(label = source),
      show.legend = FALSE
    ) +
    labs(
      y = "-log10(padj)",
      x = "score",
      title = TITLE,
    ) +
    scale_colour_gradient(low = "dodgerblue", high = "black") +
    ylim(c(0.5, NA))
  
  print(p)
}


# print to pdf ------------------------------------------------------------

pdf(file.path(OUTPUT_DIR, "predicted_TFs_rnaseq.pdf"), height=8, width=10)

decoupleR_plot1 <- volcanoplotter(ulm_DE, "decoupleR predicted TFs from RNAseq, E14.5, DE genes only")
decoupleR_plot2 <- volcanoplotter(ulm_all, "decoupleR predicted TFs from RNAseq, E14.5, all transcripts")

# pheatmap
n_top <- 50

top_tfs <- names(
  sort(
    apply(mat_ulm, 1, var, na.rm = TRUE),
    decreasing = TRUE
  )
)[seq_len(min(n_top, nrow(mat_ulm)))]

mat_subset <- mat_ulm[top_tfs, , drop = FALSE]

# Color scale
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use <- grDevices::colorRampPalette(colors = colors)(100)

# Heatmap
decoupleR_pheatmap <- pheatmap::pheatmap(mat = mat_subset,
                                         color = colors.use,
                                         border_color = "white",
                                         cluster_rows = TRUE,
                                         cluster_cols = TRUE,
                                         main="Predicted TFs driving DE genes, E14.5, minsize=5")


dev.off()

# save to RData for shiny app
save(
  dorothea_net,
  rna_unique,
  ulm_all,
  ulm_DE,
  decoupleR_plot1,
  decoupleR_plot2,
  decoupleR_pheatmap,
  file = file.path(OUTPUT_DIR, "decoupleR_results.RData")
)

# look up single TF in our dataset ----------------------------------------

mytf <- "Smad1"

tf_targets <- dorothea_net %>%
  filter(tf == mytf)

tf_targets

tf_res <- tf_targets %>%
  left_join(
    rna_unique %>%
      select(Gene_Symbol, stat),
    by = c("target" = "Gene_Symbol")
  ) %>%
  mutate(contribution = mor * stat) %>%
  arrange(desc(abs(contribution)))

tf_res

ulm_all[which(ulm_all$source == mytf),]

# look up a gene's TFs in our dataset----------------------------------------------------

mygene <- "Cox11"

gene_tfs <- dorothea_net %>%
  filter(target == mygene)

gene_tfs %>%
  left_join(
    rna_unique %>%
      select(Gene_Symbol, stat),
    by = c("target" = "Gene_Symbol")
  ) %>%
  mutate(contribution = mor * stat) %>%
  arrange(desc(abs(contribution)))
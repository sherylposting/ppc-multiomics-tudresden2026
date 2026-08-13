library(decoupleR)
library(dorothea)
library(dplyr)
library(tidyr)
library(pheatmap)
library(OmnipathR)
library(readxl)

# global variables - check these ------------------------------------------
RNA_PATH <- "V:/MARIA/Sheryl/integration_results/virtual_KO/virtual_KO_df.RData"

OUTPUT_DIR <- "V:/MARIA/Sheryl/integration_results/virtual_KO"

MINSIZE <- 5

dir.create(OUTPUT_DIR)

# load RNA seq results
load(RNA_PATH)

# prepare dorothea network ------------------------------------------------

# omnipath / collectri doesn't work, maybe need new version of R
data("dorothea_mm", package = "dorothea")

# filter higher confidence targets
dorothea_net <- dorothea_mm %>%
  filter(confidence %in% c("A", "B", "C")) %>%
  select(tf, target, mor)

# prepare rna df to decoupleR matrix --------------------------------------

res <- res %>%
  filter(
    !is.na(stat),
    !is.na(Gene_Symbol))

# mat for all genes
mat_all <- res %>%
  select(stat) %>%
  as.matrix()

rownames(mat_all) <- res$Gene_Symbol

# mat for DE genes only
rna_df_diff <- res[res$padj < 0.05, ] %>%
  select(stat, Gene_Symbol) %>%
  drop_na()

mat_DE <- rna_df_diff %>%
  select(stat) %>%
  as.matrix()

rownames(mat_DE) <- rna_df_diff$Gene_Symbol

# mat for DE genes only, all samples
rna_df_diff <- res[res$padj < 0.05, ] %>%
  select(Gene_Symbol, contains(c("KO", "WT"))) %>%
  drop_na()

# separate df, remove duplicate entries
rna_unique <- res %>%
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
volcanoplotter <- function(DF, TITLE, XLIM = NULL, YLIM = NULL) {
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
    scale_colour_gradient(low = "dodgerblue", high = "black")
  
  if (!is.null(XLIM))
    p <- p + xlim(XLIM)
  
  if (!is.null(YLIM))
    p <- p + ylim(YLIM)
  
  print(p)
}

# print to pdf ------------------------------------------------------------

pdf(file.path(OUTPUT_DIR, "predicted_TFs_virtualKO.pdf"), height=8, width=10)

decoupleR_plot1 <- volcanoplotter(ulm_all, "decoupleR predicted TFs from virtual KO, E14.5, all transcripts")
decoupleR_plot2 <- volcanoplotter(ulm_all, YLIM = c(0,5), XLIM = c(-5,5.1), "decoupleR predicted TFs from virtual KO, E14.5, all transcripts, zoomed")

dev.off()

# save to RData for shiny app
save(
  dorothea_net,
  rna_unique,
  ulm_all,
  ulm_DE,
  decoupleR_plot1,
  decoupleR_plot2,
  file = file.path(OUTPUT_DIR, "decoupleR_results.RData")
)

# write results to txt
write.table(ulm_all[which(ulm_all$padj < 0.05),], file = file.path(OUTPUT_DIR, "decoupleR_virtualKO.csv"))
write(ulm_all[which(ulm_all$padj < 0.05 & ulm_all$score > 0),]$source, file = file.path(OUTPUT_DIR, "decoupleR_virtualKO_upTFs.txt"))
write(ulm_all[which(ulm_all$padj < 0.05 & ulm_all$score < 0),]$source, file = file.path(OUTPUT_DIR, "decoupleR_virtualKO_downTFs.txt"))

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
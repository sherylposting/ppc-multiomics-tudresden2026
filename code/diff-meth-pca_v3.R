# these should be specified in the sbatch script
LIBPATH <- Sys.getenv("LIBPATH", unset = "/data/horse/ws/shli842i-p_dna15_1/rpacks")
WORKDIR <- Sys.getenv("WORKDIR", unset = "/home/shli842i/p_dna15")
VERSION <- Sys.getenv("VERSION", unset = "v3.x")

# link R packages installed in temp directory on cluster
.libPaths(c(LIBPATH, .libPaths()))

library(methylKit)
library(ggplot2)
library(ggrepel)

# set wd to home folder on cluster (writeable)
setwd(WORKDIR)

# -------------------------------------------------------------------------

# these were from v2 (min 2x cutoff, omit WT1)
load(paste0("code/RData/methBase_preliminary_", VERSION, ".RData")) # -> methBase
load(paste0("code/RData/myDiff_diff-meth_", VERSION, ".RData")) # -> myDiff, myDiff_df
TREATMENT <- c(0,0,0,0,0,0,0,1,1,1,1,1)

PCA_SAVENAME <- paste0("results/sig-diffmeth-pca_diff-meth_", VERSION, ".pdf")

# -------------------------------------------------------------------------

sig_diffmeth <- myDiff[myDiff$qvalue < 0.05]

sig_methBase <- methBase[which(methBase$start %in% sig_diffmeth$start)]

pdf(PCA_SAVENAME, width = 7, height = 6)

pca <- PCASamples(sig_methBase, obj.return = TRUE)

pca_plotdf <- as.data.frame(pca$x)
pca_plotdf$group <- c('WT', 'KO')[TREATMENT+1] # label on pca with group names
var_explained <- 100 * (pca$sdev^2 / sum(pca$sdev^2))

# pca plot
pca_plot <- ggplot(pca_plotdf, aes(PC1, PC2, color=group)) +
  geom_point() +
  geom_text_repel(aes(label = rownames(pca_plotdf))) +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2]),
    title = '3x cutoff, omit none, after differential methylation analysis'
  )

print(pca_plot)

dev.off()
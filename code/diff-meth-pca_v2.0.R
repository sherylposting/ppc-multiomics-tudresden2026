# these should be specified in the sbatch script
LIBPATH <- Sys.getenv("LIBPATH", unset = "/data/horse/ws/shli842i-p_dna15_1/rpacks")
WORKDIR <- Sys.getenv("WORKDIR", unset = "/home/shli842i/p_dna15")

# link R packages installed in temp directory on cluster
.libPaths(c(LIBPATH, .libPaths()))

library(methylKit)
library(ggplot2)
library(ggrepel)

# set wd to home folder on cluster (writeable)
setwd(WORKDIR)

# -------------------------------------------------------------------------

# these were from v2 (min 2x cutoff, omit WT1)
load("code/RData/methBase_preliminary_v2.4.RData") # -> methBase
load("code/RData/myDiff_diff-meth_v2.4.RData") # -> myDiff, myDiff_df

TREATMENT <- c(0,0,0,0,0,0,0,1,1,1,1,1)

# -------------------------------------------------------------------------

sig_diffmeth <- myDiff[myDiff$qvalue < 0.05]

sig_methBase <- methBase[which(methBase$start %in% sig_diffmeth$start)]

pdf("results/sig-diffmeth-pca2_diff-meth_v2.4.pdf", width = 7, height = 6)

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
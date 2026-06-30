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

# these were from v1 (min 3x cutoff)
load("code/RData/methBase_preliminary_v1.2.RData") # -> methBase
load("code/RData/myDiff_diff-meth_v1.1.RData") # -> myDiff, myDiff_df

TREATMENT <- c(0,0,0,0,0,0,0,1,1,1,1,1)
LONG_SAMPLE_NAMES <- c("L188015_WT1","L188016_WT2","L188017_WT3","L188018_WT4","L188019_WT5","L188020_WT6","L188021_WT7","L188022_KO1","L188023_KO2","L188024_KO3","L188025_KO4","L188026_KO5")

# -------------------------------------------------------------------------

sig_diffmeth <- myDiff[myDiff$qvalue < 0.05]

sig_methBase <- methBase[which(methBase$start %in% sig_diffmeth$start)]

# omit WT1
sig_methBase <- reorganize(sig_methBase,
           sample.ids = LONG_SAMPLE_NAMES[-1],
           treatment = TREATMENT[-1])

pdf("results/sig-diffmeth-pca2_diff-meth_v1.1.pdf", width = 7, height = 6)

pca <- PCASamples(sig_methBase, obj.return = TRUE)

pca_plotdf <- as.data.frame(pca$x)
pca_plotdf$group <- c('WT', 'KO')[TREATMENT[-1]+1] # label on pca with group names
var_explained <- 100 * (pca$sdev^2 / sum(pca$sdev^2))

# pca plot
pca_plot <- ggplot(pca_plotdf, aes(PC1, PC2, color=group)) +
  geom_point() +
  geom_text_repel(aes(label = rownames(pca_plotdf))) +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2]),
    title = '3x cutoff, omit WT1, after differential methylation analysis'
  )

print(pca_plot)

dev.off()
# link R packages installed in temp directory on cluster
.libPaths(c("/data/horse/ws/shli842i-p_dna15_1/rpacks", .libPaths()))

library(methylKit)
library(ggplot2)
library(ggrepel)

# set wd to home folder on cluster (writeable)
setwd('/home/shli842i/p_dna15')

# load workspace object from preliminary
load("code/methBase_preliminary_v1.2.RData")

cat('everything is loaded')

# differential methylation ------------------------------------------------

myDiff=calculateDiffMeth(methBase, mc.cores=8)

save(myDiff, file = "code/methylKit_diffmeth.RData")

# assign 'direction' column for coloring on plot
DE_df <- data.frame(myDiff, direction = c('negative','positive')[(myDiff$meth.diff>0)+1])
# cutoff q < 0.05 
DE_df[abs(DE_df$qvalue) > 0.05,]$direction <- 'none'

DE_volcano <- ggplot(data=DE_df, aes(x=meth.diff, y=-log10(qvalue), col=direction))+
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size = 2) +
  scale_color_manual(values = c("cyan3", "grey", "red3"))

ggsave(
  filename = "results/volcano-plot_diff-meth_v1.0.pdf",
  plot = DE_volcano,
  width = 8,
  height = 6,
  units = "in"
)

save(list = c("myDiff", "DE_df"), file = "code/methylKit_diffmeth.RData")

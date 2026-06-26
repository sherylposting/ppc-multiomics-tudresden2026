# these should be specified in the sbatch script
LIBPATH <- Sys.getenv("LIBPATH", unset = NULL)
WORKDIR <- Sys.getenv("WORKDIR", unset = NULL)

# link R packages installed in temp directory on cluster
.libPaths(c(LIBPATH, .libPaths()))

library(methylKit)
library(ggplot2)

# set wd to home folder on cluster (writeable)
setwd(WORKDIR)

# differential methylation ------------------------------------------------

load("code/RData/myDiff_diff-meth_v2.3.RData") # -> myDiff, myDiff_df

DE_volcano <- ggplot(data=myDiff_df, aes(x=meth.diff, y=-log10(qvalue), col=direction))+
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size = 2) +
  scale_color_manual(values = c("cyan3", "grey", "red3"))

ggsave(
  filename = "results/volcano-plot_diff-meth_v2.0.png",
  plot = DE_volcano,
  width = 8,
  height = 6,
  units = "in"
)
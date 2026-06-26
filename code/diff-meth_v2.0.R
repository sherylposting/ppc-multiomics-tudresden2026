# link R packages installed in temp directory on cluster
.libPaths(c("/data/horse/ws/shli842i-p_dna15_1/rpacks", .libPaths()))

library(methylKit)
library(ggplot2)
library(ggrepel)

# set wd to home folder on cluster (writeable)
setwd('/home/shli842i/p_dna15')

# load workspace object from preliminary
load("code/RData/methBase_preliminary_v2.2.RData")
cat('everything is loaded')

# differential methylation ------------------------------------------------

myDiff=calculateDiffMeth(methBase, mc.cores=8)
cat('diff meth analysis is done')

save(myDiff, file = "code/RData/myDiff_diff-meth_v2.0.RData")
cat('diff meth analysis is saved')

# assign 'direction' column for coloring on plot
myDiff_df <- data.frame(myDiff, direction = c('negative','positive')[(myDiff$meth.diff>0)+1])
# assign none direction to cutoff q < 0.05 and sort
myDiff_df[abs(myDiff_df$qvalue) > 0.05,]$direction <- 'none'
myDiff_df <- myDiff_df[order(myDiff_df$qvalue),]

save(list = c("myDiff", "myDiff_df"), file = "code/RData/myDiff_diff-meth_v2.0.RData")

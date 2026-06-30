### ---------------------- ###
# this code takes a methBase object and performs differential methylation analysis between KO - WT batches
# input: 
  # methBase
# output:
  # myDiff
  # myDiff_df - data frame with diff cpgs, sorted by q value, labeled with "direction" column
# note: calculateDiffMeth is very resource-intensive and methBase is very large. using 256GB RAM and 8 cores will take 2 hours
### ---------------------- ###

# these should be specified in the sbatch script
LIBPATH <- Sys.getenv("LIBPATH", unset = "/data/horse/ws/shli842i-p_dna15_1/rpacks")
WORKDIR <- Sys.getenv("WORKDIR", unset = "/home/shli842i/p_dna15")

# link R packages installed in temp directory on cluster
.libPaths(c(LIBPATH, .libPaths()))

library(methylKit)

# set wd to home folder on cluster (writeable)
setwd(WORKDIR)

# stuff for you to edit and check -----------------------------------------

# load workspace object from preliminary
load("code/RData/methBase_preliminary_v2.4.RData") # -> methBase
cat('everything is loaded')

DIFFMETH_SAVENAME <- "code/RData/myDiff_diff-meth_v2.4.RData"

# differential methylation ------------------------------------------------

myDiff=calculateDiffMeth(methBase, mc.cores = 8)
cat('diff meth analysis is done')

save(myDiff, file = DIFFMETH_SAVENAME)
cat('diff meth analysis is saved')

# assign 'direction' column for coloring on plot
myDiff_df <- data.frame(myDiff, direction = c('negative','positive')[(myDiff$meth.diff>0)+1])
# assign none direction to cutoff q < 0.05 and sort
myDiff_df[abs(myDiff_df$qvalue) > 0.05,]$direction <- 'none'
myDiff_df <- myDiff_df[order(myDiff_df$qvalue),]

save(myDiff, myDiff_df, file = DIFFMETH_SAVENAME)

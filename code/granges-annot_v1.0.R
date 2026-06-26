# link R packages installed in temp directory on cluster
.libPaths(c("/data/horse/ws/shli842i-p_dna15_1/rpacks", .libPaths()))

library(genomation)
library(methylKit)

# set wd to home folder on cluster (writeable)
setwd('/home/shli842i/p_dna15')

# -------------------------------------------------------------------------

refseq <- readTranscriptFeatures("data/GRCm39_RefSeq.bed.txt")

# annotate our dataset ----------------------------------------------------

load('code/myDiff_diff-meth_v1.1.RData')

load('code/methRawList_load-methraw_v1.0.RData')

sig_diffmeth <- read.table('results/sig_diffmeth.csv')

myDiff_feats <- annotateWithGeneParts(as(myDiff,"GRanges"), refseq)

promoters <- regionCounts(as(myDiff, "methylRawList"), refseq$promoters)

bed <- data.frame(
  chr = myDiff_df$chr,
  start = myDiff_df$start - 1,   # BED is 0-based
  end = myDiff_df$end,
  name = paste0("diffMeth_", round(myDiff_df$meth.diff, 2)),
  score = pmin(1000, -log10(myDiff_df$qvalue) * 100),
  strand = "."
)

write.table(
  bed,
  file = "results/myDiff_IGV.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
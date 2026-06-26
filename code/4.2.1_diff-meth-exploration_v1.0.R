# link R packages installed in temp directory on cluster - do not let it access packs in /home/shli842
.libPaths("/data/horse/ws/shli842i-p_dna15_1/rpacks_4.2.1")

library(methylKit)

# set wd to home folder on cluster (writeable)
setwd('/home/shli842i/p_dna15')

# load workspace object from diff-meth
load("code/myDiff_diff-meth_v1.1.RData")

top500 <- myDiff_df[1:500,]

rm(myDiff)
rm(myDiff_df)

load("code/methBase_preliminary_v1.2.RData")

methBase_df <- getData(methBase)
rm(methBase)

top500_temp <- methBase_df[as.numeric(rownames(top500)), 5:ncol(methBase_df)]
merged_coverage <- rowSums(top500_temp[, grep("^coverage", names(top500_temp), value = TRUE)], na.rm=TRUE)
top500$tot_coverage <- merged_coverage
sig_diffmeth <- top500[top500$qvalue < 0.05,]

write.table(sig_diffmeth, file='results/sig_diffmeth.csv')

cat(nrow(sig_diffmeth), 'differentially methylated CpGs') # 135
cat(nrow(sig_diffmeth[sig_diffmeth$direction=='positive',]), 'positively methylated CpGs') # 68
cat(nrow(sig_diffmeth[sig_diffmeth$direction=='negative',]), 'negatively methylated CpGs') # 67

bed <- data.frame(
  chr = sig_diffmeth$chr,
  start = sig_diffmeth$start - 1,   # BED is 0-based
  end = sig_diffmeth$end,
  name = paste0("diffMeth_", round(sig_diffmeth$meth.diff, 2)),
  score = pmin(1000, -log10(sig_diffmeth$qvalue) * 100),
  strand = "."
)

write.table(
  bed,
  file = "results/sig-diffmeth.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
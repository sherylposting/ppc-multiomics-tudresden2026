### ---------------------- ###
# this code takes a methylDiff object and provides regional annotation information (ex. sites in promoters, other gene features)
# input: 
  # myDiff
  # methRawList
  # GRCm39_RefSeq.bed.txt reference genome - see lab notebook for how to download
# output:
  # sig_diffmeth.csv - table of significantly differential sites
  # hyper-sig-diffmeth.bed, hypo-sig-diffmeth.bed - BED files to view significant sites on IGV viewer
### ---------------------- ###

# these should be specified in the sbatch script
LIBPATH <- Sys.getenv("LIBPATH", unset = NULL)
WORKDIR <- Sys.getenv("WORKDIR", unset = NULL)

# link R packages installed in temp directory on cluster
.libPaths(c(LIBPATH, .libPaths()))

library(methylKit)
library(genomation)

# set wd to home folder on cluster (writeable)
setwd(WORKDIR)

# stuff for you to edit and check -----------------------------------------

refseq <- readTranscriptFeatures("data/GRCm39_RefSeq.bed.txt")

load('code/myDiff_diff-meth_v1.1.RData')
load('code/methRawList_load-methraw_v1.0.RData')

SIG_DIFFMETH_SAVENAME <- 'results/sig-diffmeth_granges-annot_v2.0.csv'
BED_SAVENAME <- "results/sig-diffmeth_granges-annot_v2.0.bed"

# annotate our dataset ----------------------------------------------------

myDiff_feats <- annotateWithGeneParts(as(myDiff,"GRanges"), refseq)

promoters <- regionCounts(as(myDiff, "methylRawList"), refseq$promoters)

top500_temp <- methBase_df[as.numeric(rownames(top500)), 5:ncol(methBase_df)]
merged_coverage <- rowSums(top500_temp[, grep("^coverage", names(top500_temp), value = TRUE)], na.rm=TRUE)
top500$tot_coverage <- merged_coverage
sig_diffmeth <- top500[top500$qvalue < 0.05,]

write.table(sig_diffmeth, file=SIG_DIFFMETH_SAVENAME)

cat(nrow(sig_diffmeth), 'differentially methylated CpGs') # 135
cat(nrow(sig_diffmeth[sig_diffmeth$direction=='positive',]), 'positively methylated CpGs') # 68
cat(nrow(sig_diffmeth[sig_diffmeth$direction=='negative',]), 'negatively methylated CpGs') # 67


# export significant cpg as .bed files for IGV viewer ----------------------------

# create two bed files / tracks for hyper and hypo methylation
beds <- list(
  hyper = sig_diffmeth[sig_diffmeth$direction=='positive', ],
  hypo = sig_diffmeth[sig_diffmeth$direction=='negative', ]
)

for(name in names(beds)){
  diffmeth <- beds[[name]]
  
  bed <- data.frame(
    chr = diffmeth$chr,
    start = diffmeth$start - 1,   # BED is 0-based
    end = diffmeth$end,
    name = paste0("diffMeth_", round(diffmeth$meth.diff, 2)),
    score = pmin(1000, -log10(diffmeth$qvalue) * 100),
    strand = "."
  )
  
  write.table(
    bed,
    file = paste0(name, "-", BED_SAVENAME),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}
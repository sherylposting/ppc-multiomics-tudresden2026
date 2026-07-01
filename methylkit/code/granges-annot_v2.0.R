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
LIBPATH <- Sys.getenv("LIBPATH", unset = "/data/horse/ws/shli842i-p_dna15_1/rpacks")
WORKDIR <- Sys.getenv("WORKDIR", unset = "/home/shli842i/p_dna15/methylkit")

# link R packages installed in temp directory on cluster
.libPaths(c(LIBPATH, .libPaths()))

library(methylKit)
library(genomation)

# set wd to home folder on cluster (writeable)
setwd(WORKDIR)

# stuff for you to edit and check -----------------------------------------

refseq <- readTranscriptFeatures("data/GRCm39_RefSeq.bed.txt")

load('code/RData/myDiff_diff-meth_v2.4.RData')
load('code/RData/methBase_preliminary_v2.4.RData')

SIG_DIFFMETH_SAVENAME <- 'results/sig-diffmeth_granges-annot_v2.4.csv'
BED_SAVENAME <- "results/sig-diffmeth_granges-annot_v2.4.bed"

# annotate our dataset ----------------------------------------------------

myDiff_feats <- annotateWithGeneParts(as(myDiff,"GRanges"), refseq)
methBase_df <- getData(methBase)

top5000 <- myDiff_df[1:50000,]
top5000_temp <- methBase_df[as.numeric(rownames(top5000)), 5:ncol(methBase_df)]
merged_coverage <- rowSums(top5000_temp[, grep("^coverage", names(top5000_temp), value = TRUE)], na.rm=TRUE)
top5000$tot_coverage <- merged_coverage
sig_diffmeth <- top5000[top5000$qvalue < 0.05,]

write.table(sig_diffmeth, file=SIG_DIFFMETH_SAVENAME)

cat(nrow(sig_diffmeth), 'differentially methylated CpGs') # 543
cat(nrow(sig_diffmeth[sig_diffmeth$direction=='positive',]), 'positively methylated CpGs') # 204
cat(nrow(sig_diffmeth[sig_diffmeth$direction=='negative',]), 'negatively methylated CpGs') # 334


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

  BED_SAVENAME_temp <- sub(
    "results/",
    paste0("results/", name, "-"),
    BED_SAVENAME
  )
  
  write.table(
    bed,
    file = BED_SAVENAME_temp,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}
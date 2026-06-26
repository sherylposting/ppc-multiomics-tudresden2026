# link R packages installed in temp directory on cluster
.libPaths(c("/data/horse/ws/shli842i-p_dna15_1/rpacks", .libPaths()))

library(methylKit)
library(ggplot2)
library(ggrepel)

# following this methylKit tutorial:
# https://www.bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html
# input bismark .cov format:
# <chromosome>  <start position>  <end position>  <methylation percentage>  <count methylated>  <count unmethylated>

# set wd to home folder on cluster (writeable)
setwd('/home/shli842i/p_dna15')

# -------------------------------------------------------------------------

load("code/methRawList_load-methraw_v1.0.RData")
sample.names <- list("WT1","WT2","WT3","WT4","WT5","WT6","WT7","KO1","KO2","KO3","KO4","KO5")
  
# filtered: discard cpg that have <=3x coverage or <=30x
filtered.methRawList <- filterByCoverage(methRawList,lo.count=3,lo.perc=NULL,
                              hi.count=30, hi.perc=NULL)

# normalize: normalize coverage for overrepresented samples using a scaling factor (difference to median)
norm.methRawList <- normalizeCoverage(filtered.methRawList, method = "median", chunk.size = 1e+06, save.db=FALSE)

# # correlation matrix
# corr <- as.matrix(read.table("data/EM_seq_files/corr.txt")) # already done in previous run... for some reason the methylkit function doesn't return the object so i saved it as a .txt file

# coverage distribution data frame ----------------------------------------

methNormList.dfList <- lapply(norm.methRawList, getData)

# helper function to collect coverage vals + omit huge outliers for better plotting
omit_huge <- function(x){
  return(x$coverage[which(x$coverage<30)])
}

methNormList.dfList <- lapply(methNormList.dfList, omit_huge)
save(methNormList.dfList, file = "code/methRawList-dfs_preliminary_v2.0.RData")

# save summary df pre and post filtering
prefilt_summary <- lapply(methRawList, summary)
postfilt_summary <- lapply(norm.methRawList, summary)

nraw <- as.numeric(lapply(methRawList, nrow))
nfilt <- as.numeric(lapply(norm.methRawList, nrow))
filt_df <- data.frame('raw_sites' = nraw, 'filt_sites' = nfilt, 'perc_loss' = round(nfilt/nraw, 3))

save(prefilt_summary, postfilt_summary, filt_df, file = "code/filt-summarydf_preliminary_v2.0.RData")

# coverage after normalization barplot ---------------------------------

pdf('results/methylation_coverage.pdf', width=8, height=6)

# remove NAs because boxplot get angy
cov.list <- lapply(methNormList.dfList, function(x) {
  x[!is.na(x)]
})

boxplot(
  cov.list,
  names = sample.names,
  outline = FALSE,
  las = 2,
  ylab = "Coverage (x)",
  main = "Methylation coverage normalized (x)"
)


# coverage after normalization histogram ------------------------------

par(mfrow=c(3,4), mar=c(4,4,2,1), oma=c(0,0,3,0))

for(i in 1:length(methNormList.dfList)){
  hist(methNormList.dfList[[i]], 
       xlim=c(0,30), main=sample.names[i], xlab = "")
}

mtext("Methylation coverage normalized (x)", side=3,line=1,outer=TRUE,cex=1.2,font=2)

# sanity check plot on filtered un-norm

methFiltList.dfList <- lapply(filtered.methRawList, getData)
methFiltList.dfList <- lapply(methFiltList.dfList, omit_huge)

par(mfrow=c(3,4), mar=c(4,4,2,1), oma=c(0,0,3,0))

for(i in 1:length(methFiltList.dfList)){
  hist(methFiltList.dfList[[i]], 
       xlim=c(0,30), main=sample.names[i], xlab = "")
}

mtext("Methylation coverage filtered (x)", side=3,line=1,outer=TRUE,cex=1.2,font=2)


# also create plots for the old pre-normalized ------------------------------

methRawList.dfList <- lapply(methRawList, getData)
methRawList.dfList <- lapply(methRawList.dfList, omit_huge)

save(methNormList.dfList, methRawList.dfList, file = "code/methRawList-dfs_preliminary_v2.0.RData")

cov.list <- lapply(methRawList.dfList, function(x) {
  x[!is.na(x)]
})

par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), oma=c(0,0,0,0))

boxplot(
  cov.list,
  names = sample.names,
  outline = FALSE,
  las = 2,
  ylab = "Coverage (x)",
  main = "Methylation coverage un-normalized (x)"
)

par(mfrow=c(3,4), mar=c(4,4,2,1), oma=c(0,0,3,0))

for(i in 1:length(methRawList.dfList)){
  hist(methRawList.dfList[[i]], 
       xlim=c(0,30), main=sample.names[i], xlab = "")
}

mtext("Methylation coverage un-normalized (x)", side=3,line=1,outer=TRUE,cex=1.2,font=2)

dev.off()

# # unite and pca -----------------------------------------------------------
# 
# # unite: merge all samples into one object
# methBase <- unite(norm.methRawList, destrand=FALSE, min.per.group = 1L)
# save(methBase, file = "code/methBase_preliminary_v2.0.RData")
#
# # pca
# pca <- PCASamples(methBase, obj.return=TRUE)
# 
# # pca plot ----------------------------------------------------------------
# 
# # plotting df prep
# pca_plotdf <- as.data.frame(pca$x)
# pca_plotdf$group <- c('WT','WT','WT','WT','WT','WT','WT','KO','KO','KO','KO','KO')
# 
# # plot
# pca_plot <- ggplot(pca_plotdf, aes(PC1, PC2, color=group)) +
#   geom_point() +
#   geom_text_repel(aes(label = rownames(pca_plotdf)))
# 
# ggsave(
#   filename = "results/pca-plot_preliminary_v2.0.pdf",
#   plot = pca_plot,
#   width = 8,
#   height = 6,
#   units = "in"
# )

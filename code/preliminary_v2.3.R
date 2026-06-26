# link R packages installed in temp directory on cluster
.libPaths(c("/data/horse/ws/shli842i-p_dna15_1/rpacks", .libPaths()))

library(methylKit)
library(ggplot2)
library(ggrepel)

# set wd to home folder on cluster (writeable)
setwd('/home/shli842i/p_dna15')

# -------------------------------------------------------------------------

# load workspace objects
load("code/RData/methRawList_load-methraw_v1.0.RData") # -> methRawList

SAMPLE_NAMES <- list("WT1","WT2","WT3","WT4","WT5","WT6","WT7","KO1","KO2","KO3","KO4","KO5")
LONG_SAMPLE_NAMES <- list("L188015_WT1","L188016_WT2","L188017_WT3","L188018_WT4","L188019_WT5","L188020_WT6","L188021_WT7","L188022_KO1","L188023_KO2","L188024_KO3","L188025_KO4","L188026_KO5")
TREATMENT <- c(0,0,0,0,0,0,0,1,1,1,1,1)

PLOTS_FILENAME <- "results/coverage_preliminary_v2.3.pdf"
PCA_FILENAME <- "results/pca-plot_preliminary_v2.3.pdf"

FILT_LO_COUNT <- 2 # discard <=2x coverage
FILT_HI_COUNT <- 30 # discard >30x coverage


# function definitions ----------------------------------------------------

# helper function to collect coverage vals + omit huge outliers for better plotting
omit_huge <- function(x){
  return(x$coverage[which(x$coverage<30)])
}

# function to make boxplots
boxplotter <- function(dfList, samplenames, title, ...){
  par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), oma=c(0,0,0,0))
  
  # remove NAs because boxplot get angy
  cov.list <- lapply(dfList, function(x) {
    x[!is.na(x)]
  })
  
  boxplot(
    cov.list,
    names = samplenames,
    outline = FALSE,
    las = 2,
    ylab = "Coverage (x)",
    main = title,
    ...
  )
}

# function to make histograms
histogrammer <- function(dfList, samplenames, title, ...){
  par(mfrow=c(3,4), mar=c(4,4,2,1), oma=c(0,0,3,0))
  
  for(i in 1:length(dfList)){
    hist(dfList[[i]], 
         xlim=c(0,30), 
         main=samplenames, 
         xaxt = "n",
         breaks = seq(0.5, 30.5, by = 1),
         xlab = "")
    
    axis(1, at = 0:30)
  }
  
  mtext(title, side=3,line=1,outer=TRUE,cex=1.2,font=2)
}

# -------------------------------------------------------------------------

# filtered: omit outlier samples (WT1)
filt.SAMPLE_NAMES <- LONG_SAMPLE_NAMES[-1]
filt.TREATMENT <- TREATMENT[-1]
omitted.methRawList <- reorganize(methRawList, 
                                  sample.ids = filt.SAMPLE_NAMES, 
                                  TREATMENT = filt.TREATMENT)

# reassign short sample names for plot
filt.SAMPLE_NAMES <- SAMPLE_NAMES[-1]

# filtered: discard cpg that have <=2x coverage or >30x
filtered.methRawList <- filterByCoverage(omitted.methRawList,
                                         lo.count=FILT_LO_COUNT,
                                         lo.perc=NULL,
                                         hi.count=FILT_HI_COUNT, 
                                         hi.perc=NULL)


# normalize: normalize coverage for overrepresented samples using a scaling factor (difference to median)
norm.methRawList <- normalizeCoverage(filtered.methRawList, method = "median", chunk.size = 1e+06, save.db=FALSE)

# # correlation matrix
# corr <- as.matrix(read.table("data/EM_seq_files/corr.txt")) # already done in previous run... for some reason the methylkit function doesn't return the object so i saved it as a .txt file

save(filtered.methRawList, norm.methRawList, file = "code/RData/processed-methRawLists_preliminary_v2.3.RData")

# pre and post filtering summary data frame ----------------------------------------

# save summary df pre and post filtering + normalization
preproc_summary <- lapply(methRawList, summary)
postproc_summary <- lapply(norm.methRawList, summary)

nraw <- as.numeric(lapply(omitted.methRawList, nrow))
nproc <- as.numeric(lapply(norm.methRawList, nrow))
proc_df <- data.frame('raw_sites' = nraw, 'final_sites' = nproc, 'perc_loss' = round(nproc/nraw, 3))

save(preproc_summary, postproc_summary, proc_df, file = "code/RData/processed-summarydfs_preliminary_v2.3.RData")


# prep dfLists for plotting -----------------------------------------------

methNormList.dfList <- lapply(norm.methRawList, getData)
methNormList.dfList <- lapply(methNormList.dfList, omit_huge)
methRawList.dfList <- lapply(methRawList, getData)
methRawList.dfList <- lapply(methRawList.dfList, omit_huge)
methFiltList.dfList <- lapply(filtered.methRawList, getData)
methFiltList.dfList <- lapply(methFiltList.dfList, omit_huge)

save(methNormList.dfList, methRawList.dfList, methFiltList.dfList, file = "code/RData/methRawList-dfs_preliminary_v2.3.RData")

# plot coverage after normalization ---------------------------------

# open pdf device to save all the plots
pdf(PLOTS_FILENAME, width=8, height=6)

boxplotter(methNormList.dfList, SAMPLE_NAMES, title = "Methylation coverage normalized (x)")
histogrammer(methNormList.dfList, LONG_SAMPLE_NAMES, title = "Methylation coverage normalized (x)")

dev.off()

# unite and pca normalized -----------------------------------------------------------

# unite: merge all samples into one object
methBase <- unite(norm.methRawList, destrand=FALSE, min.per.group = 1L)
save(methBase, file = "code/RData/methBase_preliminary_v2.3.RData")

# pca
pca <- PCASamples(methBase, obj.return=TRUE)

# pca plot of normalized ----------------------------------------------------------------

# plotting df prep
pca_plotdf <- as.data.frame(pca$x)
pca_plotdf$group <- c('WT', 'KO')[filt.TREATMENT+1] # label on pca with group names. assumes filtered

# plot
pca_plot <- ggplot(pca_plotdf, aes(PC1, PC2, color=group)) +
  geom_point() +
  geom_text_repel(aes(label = rownames(pca_plotdf)))

ggsave(
  filename = PCA_FILENAME,
  plot = pca_plot,
  width = 8,
  height = 6,
  units = "in"
)

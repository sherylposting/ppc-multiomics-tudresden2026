### ---------------------- ###
# this code takes a methRawList object and outputs summary plots, to help decide on a filtering threshold that maximizes the variance
# input: 
  # methRawList
# output:
  # coverage boxplots, histograms, pca saved as pdf into a directory "results/filter-exploration"
### ---------------------- ###

# these should be specified in the sbatch script
LIBPATH <- Sys.getenv("LIBPATH", unset = NULL)
WORKDIR <- Sys.getenv("WORKDIR", unset = NULL)

# link R packages installed in temp directory on cluster
.libPaths(c(LIBPATH, .libPaths()))

library(methylKit)
library(ggplot2)
library(ggrepel)

# set wd to home folder on cluster (writeable)
setwd(WORKDIR)

# stuff for you to edit and check -----------------------------------------

# load workspace objects
load("code/RData/methRawList_load-methraw_v1.0.RData") # -> methRawList

SAMPLE_NAMES <- list("WT1","WT2","WT3","WT4","WT5","WT6","WT7","KO1","KO2","KO3","KO4","KO5")
LONG_SAMPLE_NAMES <- list("L188015_WT1","L188016_WT2","L188017_WT3","L188018_WT4","L188019_WT5","L188020_WT6","L188021_WT7","L188022_KO1","L188023_KO2","L188024_KO3","L188025_KO4","L188026_KO5")
TREATMENT <- c(0,0,0,0,0,0,0,1,1,1,1,1)

PLOTS_FILENAME <- "results/filter-exploration/filter-exploration-"
dir.create("results/filter-exploration", recursive = TRUE, showWarnings = FALSE)

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

# function to run the entire pipeline for each filtering threshold
pipeline <- function(rawList, filt.lo.count, omitted){
  
  # plot coverage of full data
  boxplotter(rawList.dfList, SAMPLE_NAMES, title = "Methylation coverage full data (x)")
  histogrammer(rawList.dfList, SAMPLE_NAMES, title = "Methylation coverage full data (x)")
  
  # filtered: omit outlier samples
  keep <- setdiff(seq_along(LONG_SAMPLE_NAMES), omitted)
  filt.SAMPLE_NAMES <- LONG_SAMPLE_NAMES[keep]
  filt.TREATMENT <- TREATMENT[keep]
  omitted.rawList <- reorganize(rawList, 
                                    sample.ids = filt.SAMPLE_NAMES, 
                                    treatment = filt.TREATMENT)
  
  # reassign short sample names for plot
  filt.SAMPLE_NAMES <- SAMPLE_NAMES[keep]
  
  # filtered: discard cpg that have <=2x coverage or >30x
  filtered.rawList <- filterByCoverage(omitted.rawList,
                                           lo.count=filt.lo.count,
                                           lo.perc=NULL,
                                           hi.count=FILT_HI_COUNT, 
                                           hi.perc=NULL)
  
  
  # normalize: normalize coverage for overrepresented samples using a scaling factor (difference to median)
  norm.rawList <- normalizeCoverage(filtered.rawList, method = "median", chunk.size = 1e+06, save.db=FALSE)
  
  normList.dfList <- lapply(norm.rawList, getData)
  normList.dfList <- lapply(normList.dfList, omit_huge)
  rawList.dfList <- lapply(rawList, getData)
  rawList.dfList <- lapply(rawList.dfList, omit_huge)
  filtList.dfList <- lapply(filtered.rawList, getData)
  filtList.dfList <- lapply(filtList.dfList, omit_huge)
  
  # open pdf device to save all the plots
  pdf(paste0(PLOTS_FILENAME, filt.lo.count, "x", omitted, "omit_v2.3.pdf"), width=8, height=6)
  
  omitname <- SAMPLE_NAMES[omitted]
  if(omitted == 0){omitname="none"}
  
  # plot coverage after filtering
  boxplotter(filtList.dfList, filt.SAMPLE_NAMES, title = paste("omit", omitname, filt.lo.count, "x cutoff:", "Methylation coverage filtered (x)"))
  histogrammer(filtList.dfList, filt.SAMPLE_NAMES, title = paste("omit", omitname, filt.lo.count, "x cutoff:", "Methylation coverage filtered (x)"))
  
  # plot coverage after normalization
  boxplotter(normList.dfList, filt.SAMPLE_NAMES, title = paste("omit", omitname, filt.lo.count, "x cutoff:", "Methylation coverage normalized (x)"))
  histogrammer(normList.dfList, filt.SAMPLE_NAMES, title = paste("omit", omitname, filt.lo.count, "x cutoff:", "Methylation coverage normalized (x)"))

# -------------------------------------------------------------------------
  
  # unite: merge all samples into one object
  methBase <- unite(norm.rawList, destrand=FALSE, min.per.group = 1L)
  
  # pca
  par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), oma=c(0,0,0,0))
  pca <- PCASamples(methBase, obj.return=TRUE)
  pca_plotdf <- as.data.frame(pca$x)
  pca_plotdf$group <- c('WT', 'KO')[filt.TREATMENT+1] # label on pca with group names. assumes filtered
  var_explained <- 100 * (pca$sdev^2 / sum(pca$sdev^2))
  
  # pca plot
  pca_plot <- ggplot(pca_plotdf, aes(PC1, PC2, color=group)) +
    geom_point() +
    geom_text_repel(aes(label = rownames(pca_plotdf))) +
    labs(
      x = sprintf("PC1 (%.1f%%)", var_explained[1]),
      y = sprintf("PC2 (%.1f%%)", var_explained[2])
    )
  
  print(pca_plot)
  
  dev.off()
}


# explore with various filtering cutoffs ------------------------------
pipeline(methRawList, filt.lo.count = 2, omitted = 0)
pipeline(methRawList, filt.lo.count = 3, omitted = 0)
pipeline(methRawList, filt.lo.count = 4, omitted = 0)
pipeline(methRawList, filt.lo.count = 2, omitted = 1)
pipeline(methRawList, filt.lo.count = 3, omitted = 1)
pipeline(methRawList, filt.lo.count = 4, omitted = 1)
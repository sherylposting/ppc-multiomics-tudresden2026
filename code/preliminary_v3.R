### ---------------------- ###
# this code takes a methRawList object, performs filtering, normalization, and quality control plots (coverage plot, PCA)
# input: 
  # methRawList
# output:
  # methBase - filtered, normalized, united CpGs
  # coverage_preliminary.pdf - coverage boxplots and histograms
  # pca_preliminary.pdf - pca plot
  # processed-summarydfs_preliminary.RData - summary statistics of the counts pre- and post- normalization
### ---------------------- ###

# these should be specified in the sbatch script
LIBPATH <- Sys.getenv("LIBPATH", unset = "/data/horse/ws/shli842i-p_dna15_1/rpacks")
WORKDIR <- Sys.getenv("WORKDIR", unset = "/home/shli842i/p_dna15")
VERSION <- Sys.getenv("VERSION", unset = "v3.x")
FILT_LO_COUNT <- Sys.getenv("FILT_LO_COUNT", unset = 2)
OMITTED <- Sys.getenv("OMITTED", unset = 0)

# link R packages installed in temp directory on cluster
.libPaths(c(LIBPATH, .libPaths()))

library(methylKit)
library(ggplot2)
library(ggrepel)

# set wd to home folder on cluster (writeable)
setwd(WORKDIR)

# global variables - check these -----------------------------------------

# load workspace objects
load("code/RData/methRawList_load-methraw_v1.0.RData") # -> methRawList

SAMPLE_NAMES <- list("WT1","WT2","WT3","WT4","WT5","WT6","WT7","KO1","KO2","KO3","KO4","KO5")
LONG_SAMPLE_NAMES <- list("L188015_WT1","L188016_WT2","L188017_WT3","L188018_WT4","L188019_WT5","L188020_WT6","L188021_WT7","L188022_KO1","L188023_KO2","L188024_KO3","L188025_KO4","L188026_KO5")
TREATMENT <- c(0,0,0,0,0,0,0,1,1,1,1,1)

METHBASE_SAVENAME <- paste0("code/RData/methBase_preliminary_", VERSION, ".RData")
SUMMARY_SAVENAME <- paste0("code/RData/processed-summarydfs_preliminary_", VERSION, ".RData")
PLOTS_SAVENAME <- paste0("results/preliminary-plots_", VERSION, ".pdf")

FILT_HI_COUNT <- 30 # discard >30x coverage
UNITE_MIN <- 1L # loosest setting, does not discard any sites

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
         main=samplenames[i], 
         xaxt = "n",
         breaks = seq(0.5, 30.5, by = 1),
         ylim = c(0, 8e6),
         xlab = "")
    
    axis(1, at = 0:30)
  }
  
  mtext(title, side=3,line=1,outer=TRUE,cex=1.2,font=2)
}

# function to run the entire pipeline for each filtering threshold
pipeline <- function(rawList, filt.lo.count, omitted){
  
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
  pdf(PLOTS_SAVENAME, width=8, height=6)
  
  # plot coverage of full data
  boxplotter(rawList.dfList, SAMPLE_NAMES, title = "Methylation coverage, full data (x)")
  histogrammer(rawList.dfList, SAMPLE_NAMES, title = "Methylation coverage, full data (x)")
  
  omitname <- SAMPLE_NAMES[omitted]
  if(all(omitted== 0)){omitname="none"}
  
  # plot coverage after filtering
  boxplotter(
    filtList.dfList,
    filt.SAMPLE_NAMES,
    title = paste(
      "Omit",
      paste(omitname, collapse = ", "),
      paste0(filt.lo.count, "x cutoff:"),
      "Methylation coverage, filtered (x)"
    )
  )
  
  histogrammer(
    filtList.dfList,
    filt.SAMPLE_NAMES,
    title = paste(
      "Omit",
      paste(omitname, collapse = ", "),
      paste0(filt.lo.count, "x cutoff:"),
      "Methylation coverage, filtered (x)"
    )
  )
  
  # plot coverage after normalization
  boxplotter(
    normList.dfList,
    filt.SAMPLE_NAMES,
    title = paste(
      "Omit",
      paste(omitname, collapse = ", "),
      paste0(filt.lo.count, "x cutoff:"),
      "Methylation coverage, normalized (x)"
    )
  )
  
  histogrammer(
    normList.dfList,
    filt.SAMPLE_NAMES,
    title = paste(
      "Omit",
      paste(omitname, collapse = ", "),
      paste0(filt.lo.count, "x cutoff:"),
      "Methylation coverage, normalized (x)"
    )
  )
  
  # -------------------------------------------------------------------------
  
  # unite: merge all samples into one object
  methBase <- unite(norm.rawList, destrand=FALSE, min.per.group = 1L)
  save(methBase, file = METHBASE_SAVENAME)
  
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
  
  # pre and post filtering summary data frame ----------------------------------------
  
  # save summary df pre and post filtering + normalization
  preproc_summary <- lapply(rawList, summary)
  postproc_summary <- lapply(norm.rawList, summary)
  
  nraw <- as.numeric(lapply(omitted.rawList, nrow))
  nproc <- as.numeric(lapply(norm.rawList, nrow))
  proc_df <- data.frame('raw_sites' = nraw, 'final_sites' = nproc, 'perc_loss' = round(1-(nproc/nraw), 3))
  
  save(preproc_summary, postproc_summary, proc_df, file = SUMMARY_SAVENAME)
  
}


# explore with various filtering cutoffs ------------------------------
pipeline(methRawList, filt.lo.count = FILT_LO_COUNT, omitted = OMITTED)
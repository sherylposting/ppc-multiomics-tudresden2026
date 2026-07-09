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
LIBPATH <- Sys.getenv("LIBPATH", unset = "/data/horse/ws/shli842i-p_dna15_1/rpacks_4.2.1")
WORKDIR <- Sys.getenv("WORKDIR", unset = "/home/shli842i/p_dna15/DMRichR")
DATADIR <- Sys.getenv("DATADIR", unset = "/projects/p_dna15/data/EM_seq_files/cytosine_reports")
VERSION <- Sys.getenv("VERSION", unset = "v1.x")

# link R packages installed in temp directory on cluster
.libPaths(LIBPATH)

#remotes::install_local(
#  "/home/shli842i/DMRichR",
#  dependencies = FALSE,
#  upgrade = "never",
#  force = TRUE
#)

BiocManager::install(
  c("Rhdf5lib", "rhdf5", "rhdf5filters", "HDF5Array", "DelayedArray"),
  force = TRUE,
  update = FALSE
)

library(DMRichR) # if this doesn't load due to rlang / bsseq then R is using the wrong libpath, restart R

# set wd to home folder on cluster (writeable)
setwd(WORKDIR)

# global variables - check these -----------------------------------------

DATADIR <- "/projects/p_dna15/data/EM_seq_files/cytosine_reports"
SAMPLEINFO_PATH <- "sample_info.xlsx"
RESPATH <- "/home/shli842i/p_dna15/data/"

# load data ---------------------------------------------------------------

# Individual smoothed values ----------------------------------------------
 
load("RData/bismark.RData")
load("RData/DMRs.RData")

cores=1
genome = "mm39"

  cat("\n[DMRichR] Smoothing individual methylation values \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  start_time <- Sys.time()
  
  bs.filtered.bsseq <- bsseq::BSmooth(bs.filtered,
                                      BPPARAM = BiocParallel::MulticoreParam(workers = cores,
                                                                             progressbar = TRUE))
  
  # Drop chrY in Rat only due to poor quality (some CpGs in females map to Y)
  if(genome == "rn6"){
    bs.filtered.bsseq <- GenomeInfoDb::dropSeqlevels(bs.filtered.bsseq,
                                                     "chrY",
                                                     pruning.mode = "coarse")
    GenomeInfoDb::seqlevels(bs.filtered.bsseq)
  }
  
  bs.filtered.bsseq
  
  print(glue::glue("Extracting individual smoothed methylation values of DMRs..."))
  bs.filtered.bsseq %>%
    DMRichR::smooth2txt(regions = sigRegions,
                        txt = "DMRs/DMR_individual_smoothed_methylation.txt")
  
  print(glue::glue("Extracting individual smoothed methylation values of background regions..."))
  bs.filtered.bsseq %>%
    DMRichR::smooth2txt(regions = regions,
                        txt = "DMRs/background_region_individual_smoothed_methylation.txt")
  
  print(glue::glue("Individual smoothing timing..."))
  end_time <- Sys.time()
  end_time - start_time
  
  print(glue::glue("Saving Rdata..."))
  save(bs.filtered.bsseq,
       file = "RData/bsseq.RData")
  #load("RData/bsseq.RData")


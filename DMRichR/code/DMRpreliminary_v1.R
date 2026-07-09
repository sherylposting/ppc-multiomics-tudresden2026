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

library(DMRichR) # if this doesn't load due to rlang / bsseq then R is using the wrong libpath, restart R

# set wd to home folder on cluster (writeable)
setwd(WORKDIR)

# global variables - check these -----------------------------------------

DATADIR <- "/projects/p_dna15/data/EM_seq_files/cytosine_reports"
SAMPLEINFO_PATH <- "sample_info.xlsx"
RESPATH <- "/home/shli842i/p_dna15/data/"

# load data ---------------------------------------------------------------

DM.R(
  genome = "mm39",
  coverage = 1,
  perGroup = 0.75,
  minCpGs = 5,
  maxPerms = 10,
  maxBlockPerms = 10,
  cutoff = 0.05,
  testCovariate = "Group",
  adjustCovariate = "Sex",
  matchCovariate = NULL,
  cores = 1,
  GOfuncR = TRUE,
  sexCheck = TRUE,
  EnsDb = FALSE,
  resPath = RESPATH,
  dataPath = DATADIR,
  sampleinfoPath = SAMPLEINFO_PATH
)

save.image(file = ".RData")

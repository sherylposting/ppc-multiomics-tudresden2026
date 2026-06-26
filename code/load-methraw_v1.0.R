### ---------------------- ###
# this code loads in bismark .cov files for differential methylation analysis of CpG islands using methylKit and genomation. it creates a SINGLE methRawList object as opposed to a methylDB
# input: 
  # folder of .cov files
# output: 
  # single methRawList object at "code/methRawList_load-methraw_v1.0.RData"
# warning: the methRawList object is extremely large and will need 100+ GB RAM to work with

# following this methylKit tutorial:
# https://www.bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html
# input bismark .cov format:
# <chromosome>  <start position>  <end position>  <methylation percentage>  <count methylated>  <count unmethylated>

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

# -------------------------------------------------------------------------

data_dir <- "/projects/p_dna15/data/EM_seq_files"

long.sample.names <- list("L188015_WT1","L188016_WT2","L188017_WT3","L188018_WT4","L188019_WT5","L188020_WT6","L188021_WT7","L188022_KO1","L188023_KO2","L188024_KO3","L188025_KO4","L188026_KO5")

treatment <- c(0,0,0,0,0,0,0,1,1,1,1,1)

# -------------------------------------------------------------------------

# load bismark .cov file names from folder
file.list <- as.list(list.files(data_dir, pattern = "\\.cov.gz$", full.names = TRUE)) # name list of all .cov files in folder

# load files as single object
methRawList=methRead(file.list,
                     sample.id=long.sample.names,
                     assembly="GRCm39", # just annotation
                     treatment=treatment,
                     context="CpG", # bismark .cov / bedgraph by default returns cpg context only
                     mincov = 2, # default is 10x coverage
                     pipeline = "bismarkCoverage"
)

save.image(file = "code/methRawList_load-methraw_v1.0.RData")
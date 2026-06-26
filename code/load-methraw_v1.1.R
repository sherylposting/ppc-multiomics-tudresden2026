### this code loads in bismark .cov files for differential methylation analysis of CpG islands using methylKit and genomation. it creates a METHYLDB directory + object as opposed to a single methRawList object.
# input: folder of .cov files
# output: methDB directory, methRawList methDB object at "code/methDB.RData"
# note: this will write a new directory called methylDB and initial database creation takes a long time. it seems to take longer than a single methRawList object to work with.
# note: this is depreciated and most of my analysis has been done with a single object made by methraw_v1.0, not methraw_v1.1

# following this methylKit tutorial:
# https://www.bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html
# input bismark .cov format:
# <chromosome>  <start position>  <end position>  <methylation percentage>  <count methylated>  <count unmethylated>

# link R packages installed in temp directory on cluster
.libPaths(c("/data/horse/ws/shli842i-p_dna15_1/rpacks", .libPaths()))

library(methylKit)
library(ggplot2)
library(ggrepel)

# set wd to home folder on cluster (writeable)
setwd('/home/shli842i/p_dna15')

# load bismark .cov file names from folder -------------------------------------

file.list <- as.list(list.files("data/EM_seq_files", pattern = "\\.cov.gz$", full.names = TRUE)) # name list of all .cov files in folder

sample.names <- list("L188015_WT1","L188016_WT2","L188017_WT3","L188018_WT4","L188019_WT5","L188020_WT6","L188021_WT7","L188022_KO1","L188023_KO2","L188024_KO3","L188025_KO4","L188026_KO5")

# load files as tabix
methRawList=methRead(file.list,
                     sample.id=sample.names,
                     assembly="mm9", # just annotation
                     treatment=c(0,0,0,0,0,0,0,1,1,1,1,1),
                     context="CpG", # bismark .cov / bedgraph by default returns cpg context only
                     mincov = 2, # default is 10x coverage
                     pipeline = "bismarkCoverage",
                     dbtype = "tabix",
                     dbdir = "methylDB" # creates tabix database directory called methylDB
)

save(methRawList, file = "code/methDB.RData")
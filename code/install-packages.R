# these should be specified in the sbatch script
LIBPATH <- Sys.getenv("LIBPATH", unset = "/data/horse/ws/shli842i-p_dna15_1/rpacks")
WORKDIR <- Sys.getenv("WORKDIR", unset = "/home/shli842i/p_dna15")

# link R packages installed in temp directory on cluster
.libPaths(c(LIBPATH, .libPaths()))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

pkgs <- c(
  "methylKit",
  "genomation",
  "rtracklayer",
  "TxDb.Mmusculus.UCSC.mm39.knownGene",
  "ggplot2",
  "ggrepel"
)

BiocManager::install(pkgs, ask = FALSE, update = TRUE)
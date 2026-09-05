# global variables - check these -----------------------------------------

# these should be specified in the sbatch script
LIBPATH <- Sys.getenv("LIBPATH", unset = "/data/horse/ws/shli842i-p_dna15_1/rpacks_4.2.1_DMRichR2")
WORKDIR <- Sys.getenv("WORKDIR", unset = "/home/shli842i/p_dna15/DMRichR/v1.3")
VERSION <- Sys.getenv("VERSION", unset = "v1.x")

# set wd to home folder on cluster (writeable)
setwd(WORKDIR)

SAMPLEINFO_PATH <- "sample_info.xlsx"

# link R packages installed in temp directory on cluster
.libPaths(LIBPATH)

library(DMRichR) # if this doesn't load due to ggplot2 then R is using the wrong libpath, restart R
library(BSgenome.Mmusculus.UCSC.mm39)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)

# set arguments ---------------------------------------------------------------

genome = "mm39"
cores = 8

# -------------------------------------------------------------------------

load("RData/DMRs.RData")

sigRegions %>% 
  DMRichR::prepareHOMER(regions = regions)

DMRichR::HOMER(genome = genome,
               cores = cores)

# -------------------------------------------------------------------------

homer_bin <- "/dev/shm/conda_envs_shli842i/homer/bin"

Sys.setenv(
  PATH = paste(homer_bin, Sys.getenv("PATH"), sep = .Platform$path.sep)
)

Sys.which("findMotifsGenome.pl")

cat("\n[DMRichR] HOMER known transcription factor motif analysis \t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
tryCatch({
  if(Sys.which("findMotifsGenome.pl") == ""){
    print(glue::glue("HOMER was not detected in PATH, skipping motif analysis. Did you load the module?"))
  }else{
    print(glue::glue("HOMER was detected in PATH, now performing motif enrichment for {genome} using {cores} cores"))
    system(paste(shQuote(system.file("exec/HOMER.sh", package = "DMRichR")),genome,cores))
  }
},
error = function(error_condition) {
  print(glue::glue("There was an error with running HOMER for {genome}. 
                  Have you confirmed that {genome} is avaiable in HOMER and installed using:
                  perl {path}/configureHomer.pl -list
                  perl {path}/configureHomer.pl -install {genome}
                  Note: see https://www.biostars.org/p/443759/ for HOMER path if using a conda install.",
                   path = dirname(Sys.which("findMotifsGenome.pl"))))
})

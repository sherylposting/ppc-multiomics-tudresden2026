# these should be specified in the sbatch script
LIBPATH <- Sys.getenv("LIBPATH", unset = "/data/horse/ws/shli842i-p_dna15_1/rpacks_4.2.1_DMRichR2")
WORKDIR <- Sys.getenv("WORKDIR", unset = "/home/shli842i/p_dna15/DMRichR")
DATADIR <- Sys.getenv("DATADIR", unset = "/projects/p_dna15/data/EM_seq_files/cytosine_reports")
VERSION <- Sys.getenv("VERSION", unset = "v1.x")

# link R packages installed in temp directory on cluster
.libPaths(LIBPATH)

library(DMRichR) # if this doesn't load due to ggplot2 then R is using the wrong libpath on HPC, restart R and try library again
library(dplyr)
library(BSgenome.Mmusculus.UCSC.mm39)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)

#.rs.restartR()

# set wd to home folder on cluster (writeable)
setwd(WORKDIR)

# global variables - check these -----------------------------------------

DATADIR <- "/projects/p_dna15/data/EM_seq_files/cytosine_reports"
SAMPLEINFO_PATH <- "sample_info.xlsx"
RESPATH <- "/home/shli842i/p_dna15/data/"

# args initially passed into DMRichR
genome = "mm39"
coverage = 1
perGroup = 0.75
minCpGs = 5
maxPerms = 10
maxBlockPerms = 10
cutoff = 0.05
testCovariate = "Group"
adjustCovariate = "Sex"
matchCovariate = NULL
cores = 1
GOfuncR = TRUE
sexCheck = TRUE
EnsDb = FALSE
resPath = RESPATH
dataPath = DATADIR
sampleinfoPath = SAMPLEINFO_PATH

# -------------------------------------------------------------------------

load("RData/bsseq.RData")
load("RData/DMRs.RData")

assign("goi", BSgenome.Mmusculus.UCSC.mm39, envir = parent.frame())
assign("annoDb", "org.Mm.eg.db", envir = parent.frame())

# make bioconductor 3.16 compatible TxDb
txdb_mm39 <- GenomicFeatures::makeTxDbFromGFF(
  file = "/projects/p_dna15/data/mm39.ncbiRefSeq.gtf.gz",
  format = "gtf",
  organism = "Mus musculus"
)
assign("TxDb", TxDb.Mmusculus.UCSC.mm39.knownGene, envir = parent.frame())

AnnotationDbi::saveDb(
  txdb_mm39,
  file = "TxDb.mm39.ncbiRefSeq.sqlite"
)

# txdb_mm39 <- AnnotationDbi::loadDb(
#   "TxDb.mm39.ncbiRefSeq.sqlite"
# )

old_gene_id <- as.character(gtf$gene_id)

symbol_to_entrez <- mapIds(
  org.Mm.eg.db,
  keys = unique(na.omit(old_gene_id)),
  keytype = "SYMBOL",
  column = "ENTREZID",
  multiVals = "first"
)

gtf$original_gene_symbol <- old_gene_id
gtf$gene_id <- unname(symbol_to_entrez[old_gene_id])

# Smoothed global, chromosomal, and CGi methylation statistics ------------

dir.create("Global")

bs.filtered.bsseq %>%
  DMRichR::globalStats(genome = genome,
                       testCovariate = testCovariate,
                       adjustCovariate = adjustCovariate,
                       matchCovariate = matchCovariate) %>%
  openxlsx::write.xlsx("Global/smoothed_globalStats.xlsx")

# Global plots ------------------------------------------------------------

windows <- bs.filtered.bsseq %>%
  DMRichR::windows(goi = goi)

CpGs <- bs.filtered.bsseq %>%
  DMRichR::CpGs()

plots <- c("windows", "CpGs")

CGi <- bs.filtered.bsseq %>%
  DMRichR::CGi(genome = genome, resPath=resPath)

plots <- c("windows", "CpGs", "CGi")

purrr::walk(plots,
            function(plotMatrix,
                     group =  bs.filtered.bsseq %>%
                       pData() %>%
                       dplyr::as_tibble() %>%
                       dplyr::pull(!!testCovariate) %>%
                       forcats::fct_rev()){

              title <- dplyr::case_when(plotMatrix == "windows" ~ "20Kb Windows",
                                        plotMatrix == "CpGs" ~ "Single CpG",
                                        plotMatrix == "CGi" ~ "CpG Island")

              plotMatrix %>%
                get() %>%
                DMRichR::PCA(testCovariate = testCovariate,
                             bs.filtered.bsseq = bs.filtered.bsseq) %>%
                ggplot2::ggsave(glue::glue("Global/{title} PCA.pdf"),
                                plot = .,
                                device = NULL,
                                width = 11,
                                height = 8.5)

              plotMatrix %>%
                get() %>%
                DMRichR::densityPlot(group = group) %>%
                ggplot2::ggsave(glue::glue("Global/{title} Density Plot.pdf"),
                                plot = .,
                                device = NULL,
                                width = 11,
                                height = 4)

              Glimma::glMDSPlot(plotMatrix %>%
                                  get(),
                                groups = cbind(bsseq::sampleNames(bs.filtered.bsseq),
                                               pData(bs.filtered.bsseq)) %>%
                                  dplyr::as_tibble() %>%
                                  dplyr::select(-col) %>%
                                  dplyr::rename(Name = bsseq..sampleNames.bs.filtered.bsseq.),
                                path = getwd(),
                                folder = "interactiveMDS",
                                html = glue::glue("{title} MDS plot"),
                                launch = FALSE)
            })

# Heatmap -----------------------------------------------------------------

sigRegions %>%
  DMRichR::smoothPheatmap(bs.filtered.bsseq = bs.filtered.bsseq,
                          testCovariate = testCovariate)

# CpG and genic enrichment testing ----------------------------------------

cat("\n[DMRichR] Performing DMRichments \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

DMRich <- function(x){
  
  if(genome %in% c("mm39", "hg38", "hg19", "mm10", "mm9", "rheMac10", "rheMac8", "rn6", "danRer11", "galGal6", "bosTau9", "panTro6", "dm6", "susScr11", "canFam3")){
    print(glue::glue("Running CpG annotation enrichments for {names(dmrList)[x]}"))
    dmrList[x] %>% 
      DMRichR::DMRichCpG(regions = regions,
                         genome = genome,
                         resPath = resPath) %T>%
      openxlsx::write.xlsx(file = glue::glue("DMRichments/{names(dmrList)[x]}_CpG_enrichments.xlsx")) %>% 
      DMRichR::DMRichPlot(type = "CpG") %>% 
      ggplot2::ggsave(glue::glue("DMRichments/{names(dmrList)[x]}_CpG_enrichments.pdf"),
                      plot = ., 
                      width = 4,
                      height = 3)
  }
  
  print(glue::glue("Running gene region annotation enrichments for {names(dmrList)[x]}"))
  dmrList[x] %>% 
    DMRichR::DMRichGenic(regions = regions,
                         TxDb = TxDb,
                         annoDb = annoDb,
                         genome = genome) %T>%
    openxlsx::write.xlsx(file = glue::glue("DMRichments/{names(dmrList)[x]}_genic_enrichments.xlsx")) %>% 
    DMRichR::DMRichPlot(type = "genic") %>% 
    ggplot2::ggsave(glue::glue("DMRichments/{names(dmrList)[x]}_genic_enrichments.pdf"),
                    plot = ., 
                    width = 4,
                    height = 4)
}

dmrList <- sigRegions %>% 
  DMRichR::dmrList()

dir.create("DMRichments")

purrr::walk(seq_along(dmrList),
            DMRich)

purrr::walk(dplyr::case_when(genome %in% c("mm39", "hg38", "hg19", "mm10", "mm9", "rn6") ~ c("CpG", "genic"),
                             TRUE ~ "genic") %>%
              unique(),
            function(type){
              
              print(glue::glue("Creating DMRichMultiPlots for {type} annotations"))
              
              DMRichR::DMparseR(direction =  c("All DMRs",
                                               "Hypermethylated DMRs",
                                               "Hypomethylated DMRs"),
                                type = type) %>%
                DMRichR::DMRichPlot(type = type,
                                    multi = TRUE) %>% 
                ggplot2::ggsave(glue::glue("DMRichments/{type}_multi_plot.pdf"),
                                plot = .,
                                device = NULL,
                                height = dplyr::case_when(type == "genic" ~ 5,
                                                          type == "CpG" ~ 3.5),
                                width = 7)
            })

# Machine learning --------------------------------------------------------
tryCatch({
  methylLearnOutput <- DMRichR::methylLearn(bs.filtered.bsseq = bs.filtered.bsseq,
                                            sigRegions = sigRegions,
                                            testCovariate = testCovariate,
                                            TxDb = TxDb,
                                            annoDb = annoDb,
                                            topPercent = 1,
                                            output = "all",
                                            saveHtmlReport = TRUE)
  
  if(!dir.exists("./Machine_learning")) {
    dir.create("./Machine_learning")
  } 
  
  if(length(methylLearnOutput) == 1) {
    openxlsx::write.xlsx(list(Annotations_Common_DMRs = methylLearnOutput), 
                         file = "./Machine_learning/Machine_learning_output_one.xlsx") 
  } else {
    openxlsx::write.xlsx(list(Annotations_Common_DMRs = methylLearnOutput$`Annotated common DMRs`,
                              RF_Ranking_All_DMRs = methylLearnOutput$`RF ranking`,
                              SVM_Ranking_All_DMRs = methylLearnOutput$`SVM ranking`),
                         file = "./Machine_learning/Machine_learning_output_all.xlsx") 
  }
  
  print(glue::glue("Saving RData..."))
  save(methylLearnOutput, file = "RData/machineLearning.RData")
  #load("RData/machineLearning.RData")
},
error = function(error_condition) {
  print(glue::glue("Warning: methylLearn did not finish. \\
                      You may have not had enough top DMRs across algrothims."))
})
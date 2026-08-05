library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)
library(org.Mm.eg.db)

OUTPUT_DIR <- "V:/MARIA/Sheryl/integration_results"

DMRICHR_PATH <- "V:/MARIA/Sheryl/DMRichR_results/260724_v1.3-loosecutoff/DMRs/DMRs_annotated.xlsx"
METHYLKIT_PATH <- "V:/MARIA/Sheryl/methylkit_results/260701_v2.4_3xcutoff-results/sig-diffmeth_granges-annot_v2.4.csv"

methylkit_df <- read.table(METHYLKIT_PATH)
methylkit_grange <- makeGRangesFromDataFrame(methylkit_df, keep.extra.columns=TRUE)

dmrichr_df <- read_xlsx(DMRICHR_PATH)
dmrichr_grange <- makeGRangesFromDataFrame(dmrichr_df, keep.extra.columns=TRUE)

hits <- findOverlaps(dmrichr_grange, methylkit_grange)

dmrichr_grange_overlap <- dmrichr_grange[unique(queryHits(hits))]
methylkit_grange_overlap <- methylkit_grange[unique(subjectHits(hits))]

write.csv(dmrichr_grange_overlap, file=file.path(OUTPUT_DIR, "dmrichr_methylkit_overlap.csv"))


# -------------------------------------------------------------------------

annotated <- annotatePeak(
  methylkit_grange,
  TxDb = TxDb.Mmusculus.UCSC.mm39.knownGene,
  annoDb = "org.Mm.eg.db"
)

gr_annotated <- as.GRanges(annotated)

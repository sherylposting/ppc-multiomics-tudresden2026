library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)
library(org.Mm.eg.db)
library(readxl)
library(GenomicRanges)
library(rtracklayer)

OUTPUT_DIR <- "V:/MARIA/Sheryl/integration_results"

DMRICHR_PATH <- "V:/MARIA/Sheryl/DMRichR_results/260724_v1.3-loosecutoff/DMRs/DMRs_annotated.xlsx"
METHYLKIT_PATH <- "V:/MARIA/Sheryl/methylkit_results/260701_v2.4_3xcutoff-results/sig-diffmeth-annotated_v2.4.csv"

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


# sanity check - figure out overlap between atac seq and either hyper or hypomethylated DMRichR (mislabeled?) --------

DMRichR_hypo <- import("V:/MARIA/Sheryl/DMRichR_results/260724_v1.3-loosecutoff/HOMER/DMRs_hypo.bed", format = "BED")
DMRichR_hyper <- import("V:/MARIA/Sheryl/DMRichR_results/260724_v1.3-loosecutoff/HOMER/DMRs_hyper.bed", format = "BED")
atac_WT <- import("V:/MARIA/Sheryl/atac_seq_results/E14.5/processed_peaks/unique_peaks/unique_WT_1bp_trimmed.bed", format = "BED")

hypo_hits <- findOverlaps(
  atac_WT,
  DMRichR_hypo,
  maxgap = 1,
  ignore.strand = TRUE
)

hyper_hits <- findOverlaps(
  atac_WT,
  DMRichR_hyper,
  maxgap = 1,
  ignore.strand = TRUE
)

n_atac_hypo <- length(unique(queryHits(hypo_hits)))
n_atac_hyper <- length(unique(queryHits(hyper_hits)))

n_hypo_dmrs_with_atac <- length(unique(subjectHits(hypo_hits)))
n_hyper_dmrs_with_atac <- length(unique(subjectHits(hyper_hits)))

pct_hypo_dmrs_with_atac <-
  100 * n_hypo_dmrs_with_atac / length(DMRichR_hypo)

pct_hyper_dmrs_with_atac <-
  100 * n_hyper_dmrs_with_atac / length(DMRichR_hyper)

pct_hypo_dmrs_with_atac
pct_hyper_dmrs_with_atac

# i think the hyper and hypo might be swapped. rather than hypermethylated in WT, might be hypomethylated. but need to ultimately check by re-running DMRichR

# check diffbind ----------------------------------------------------------

diffbind <- read.table("V:/MARIA/Sheryl/atac_seq_results/E14.5/diffbind/all_ATAC_diffgenes.E14.5.tsv")
diffbind_grange <- makeGRangesFromDataFrame(diffbind, keep.extra.columns=TRUE)

hits <- findOverlaps(dmrichr_grange, diffbind_grange)

dmrichr_grange_overlap <- dmrichr_grange[unique(queryHits(hits))]
diffbind_grange_overlap <- diffbind_grange[unique(subjectHits(hits))]

data.frame("DMRichR diff fixed" = -dmrichr_grange_overlap$difference, "diffbind fold" = diffbind_grange_overlap$Fold)

# DMRichR is WT - KO, diffbind is KO - WT
# flip DMRichR to get KO - WT to match diffbind
# DMRichR regions that are hypomethylated in the KO (accessible) are inaccessible in the KO (silenced)
# however, most diffbind peaks were inaccessible, so this could be by chance. only 1.7% of diffbind peaks were positive

length(which(diffbind$Fold > 0)) / length(diffbind$Fold) * 100
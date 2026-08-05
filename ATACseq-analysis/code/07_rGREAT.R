# this script currently outputs results for just promoters

library(rGREAT) #install TxDb mm39 too
library(clusterProfiler)

library(GenomicFeatures)
library(GenomicRanges)
library(AnnotationDbi)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)
library(org.Mm.eg.db)
library(msigdbr)
library(readxl)
library(dplyr)
library(tidyr)
library(BiocParallel)

SnowParam(8)

# global variables - check these ------------------------------------------

RNA_PATH <- "V:/MARIA/bfx2557_RNA_seq_E14.5_Katrin_no_trimming/bfx2557_all raw files/de-analysis_three_samples_each/data/bfx2557.deseq-results.separate.de-expl.xlsx"

ATAC_PATH <- "V:/MARIA/Sheryl/atac_seq_results/E14.5/diffbind/all_ATAC_diffgenes.E14.5.tsv"

EM_PATH <- "V:/MARIA/Sheryl/DMRichR_results/260724_v1.3-loosecutoff/DMRs/DMRs_annotated.xlsx"

OUTPUT_DIR <- "V:/MARIA/Sheryl/integration_results/rGREAT_results"

ANTHONY_GSEA_PATH <- "V:/MARIA/Sheryl/integration_results/R_scripts/custom_gene_sets_GSEA_Anthony_FIXED.xlsx"

# -------------------------------------------------------------------------

# load ATAC seq results
atac_df <- read.table(ATAC_PATH)

# load EM seq (DMRichR) results
em_df <- read_xlsx(EM_PATH)

# load RNA seq results
rna_df <- read_xlsx(RNA_PATH, sheet = "3_Cond_knock_out_v_wild_type_F5")

# limit to promoters
atac_promoters <- atac_df %>%
  filter(Simplified_Annotation == "Promoter")

em_promoters <- em_df %>%
  filter(annotation == "Promoter")

# prepare custom gene module granges --------------------------------------

anthony_gsea <- read_excel(ANTHONY_GSEA_PATH)

# NOTE: in the original excel, i deleted egr7, etv7, cela2b as they are deleted in mice. renamed hla-a and hla-e -> h2-k, h2-d, h2-l but there are no homologs. some others i directly renamed to match newest gene symbol

# convert from tibble to list
anthony_gsea <- anthony_gsea %>%
  as.list()

# load ref genome
txdb <- TxDb.Mmusculus.UCSC.mm39.knownGene

# prepare gene set granges object -----------------------------------------

# Gene ranges; names are Entrez gene IDs
gene_gr <- genes(txdb)

# Convert each gene to a one-base TSS range
tss_gr <- promoters(
  gene_gr,
  upstream = 0,
  downstream = 1
)

# Map Entrez IDs to mouse gene symbols
symbol_map <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = names(tss_gr),
  keytype = "ENTREZID",
  columns = c("ENTREZID", "SYMBOL")
)

# Remove duplicate Entrez mappings, if present
symbol_map <- symbol_map[!duplicated(symbol_map$ENTREZID), ]

# Match symbols back to the GRanges object
mcols(tss_gr)$gene_id <- symbol_map$SYMBOL[
  match(names(tss_gr), symbol_map$ENTREZID)
]

# Remove the Entrez IDs stored as GRanges names
names(tss_gr) <- NULL

# Optionally remove genes without a mapped symbol
tss_gr <- tss_gr[!is.na(tss_gr$gene_id)]

anthony_gsea_entrez <- lapply(
  anthony_gsea,
  function(symbols) {
    # Remove existing missing values
    symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
    
    entrez <- AnnotationDbi::mapIds(
      org.Mm.eg.db,
      keys = symbols,
      keytype = "SYMBOL",
      column = "ENTREZID",
      multiVals = "first"
    )
    
    # Remove genes that could not be mapped
    unname(entrez[!is.na(entrez)])
  }
)

# double check list is legit
anthony_gsea_entrez


# prep mgsigdbr and kegg sets ------------------------------------------------------------------

msigdbr_all <- msigdbr(species="Mus musculus")

genes = getGenomeDataFromNCBI(getKEGGGenome("mmu"))
kegg_gene_sets = getKEGGPathways("mmu")

# c2_gene_sets = split(c2_gene_sets$entrez_gene, c2_gene_sets$gs_name)
# c2_gene_sets = lapply(c2_gene_sets, as.character)  # just to make sure gene IDs are all in character.

# run all GREAT analyses on ATAC --------------------------------------------------

# GO:BP
go.bp_great_atac <- great(makeGRangesFromDataFrame(atac_promoters), "GO:BP", "txdb:mm39", cores=8)

# GO:MF
go.mf_great_atac <- great(makeGRangesFromDataFrame(atac_promoters), "GO:MF", "txdb:mm39", cores=8)

# GO:CC
go.cc_great_atac <- great(makeGRangesFromDataFrame(atac_promoters), "GO:CC", "txdb:mm39", cores=8)

# # KEGG
kegg_great_atac <- great(makeGRangesFromDataFrame(atac_promoters), kegg_gene_sets, "mm39", cores=8)

# anthony custom genesets
anthony_great_atac <- great(makeGRangesFromDataFrame(atac_promoters), anthony_gsea_entrez, "mm39", cores=8)

save(go.bp_great_atac, go.mf_great_atac, go.cc_great_atac, kegg_great_atac, anthony_great_atac, file = file.path(OUTPUT_DIR, "rGREAT_results_ATAC_promoters.RData"))

# run all GREAT analyses on EM --------------------------------------------------

# GO:BP
go.bp_great_em <- great(makeGRangesFromDataFrame(em_promoters), "GO:BP", "txdb:mm39", cores=8)

# GO:MF
go.mf_great_em <- great(makeGRangesFromDataFrame(em_promoters), "GO:MF", "txdb:mm39", cores=8)

# GO:CC
go.cc_great_em <- great(makeGRangesFromDataFrame(em_promoters), "GO:CC", "txdb:mm39", cores=8)

# # KEGG
kegg_great_em <- great(makeGRangesFromDataFrame(em_promoters), kegg_gene_sets, "mm39", cores=8)

# anthony custom genesets
anthony_great_em <- great(makeGRangesFromDataFrame(em_promoters), anthony_gsea_entrez, "mm39", cores=8)

save(go.bp_great_em, go.mf_great_em, go.cc_great_em, kegg_great_em, anthony_great_em, file = file.path(OUTPUT_DIR, "rGREAT_results_EM_promoters.RData"))

# -------------------------------------------------------------------------

atac_tables <- list(
  GO_BP = go.bp_great_atac@table,
  GO_MF = go.mf_great_atac@table,
  GO_CC = go.cc_great_atac@table,
  KEGG  = kegg_great_atac@table,
  ANTHONY = anthony_great_atac@table
)

em_tables <- list(
  GO_BP = go.bp_great_em@table,
  GO_MF = go.mf_great_em@table,
  GO_CC = go.cc_great_em@table,
  KEGG  = kegg_great_em@table,
  ANTHONY = anthony_great_em@table
)

intersection_results <- list()

for (analysis_name in names(atac_tables)) {
  intersection_results[[analysis_name]] <- atac_tables[[analysis_name]] %>%
    inner_join(
      em_tables[[analysis_name]],
      by = "id",
      suffix = c("_ATAC", "_EM")
    ) %>%
    arrange(p_adjust_ATAC)
  
  write.csv(
    atac_tables[[analysis_name]],
    file = file.path(OUTPUT_DIR, paste0(analysis_name, "_atac.csv")),
    row.names = FALSE
  )
  
  write.csv(
    em_tables[[analysis_name]],
    file = file.path(OUTPUT_DIR, paste0(analysis_name, "_em.csv")),
    row.names = FALSE
  )
  
  write.csv(
    intersection_results[[analysis_name]],
    file = file.path(OUTPUT_DIR, paste0(analysis_name, "_intersection.csv")),
    row.names = FALSE
  )
}

save(intersection_results, file=file.path(OUTPUT_DIR, "rGREAT_intersection_results_promoters.RData"))

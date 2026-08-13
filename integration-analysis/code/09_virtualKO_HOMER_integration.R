library(dplyr)
library(stringr)
library(RobustRankAggreg)
library(dorothea)
library(eulerr)
library(tidyr)

# global variables - check these ------------------------------------------

DATA_DIR <- "V:/MARIA/Sheryl/individual_hits_tables/"

load(file.path(DATA_DIR, "individual_dfs.RData")) # -> ATAC_diffgenes, DMRichR_diffgenes, ATAC_HOMER, DMRichR_HOMER, decoupleR_all. this was previously custom-saved

# new "corrected" analysis? use "hypo" as referring to hypomethylation / accessible in the WT (rather than the KO)
DMRichR_HOMER <- read.delim("V:/MARIA/Sheryl/DMRichR_results/260724_v1.3-loosecutoff/HOMER/hypo/knownResults.txt")
OUTPUT_DIR <- "V:/MARIA/Sheryl/integration_results/virtual_KO/HOMER_integrated"

# new virtual KO decoupleR results
load("V:/MARIA/Sheryl/integration_results/virtual_KO/E16/decoupleR_results.RData")
decoupleR_all <- ulm_all

dir.create(OUTPUT_DIR)

# match HOMER results to dorothea simple TF names -------------------------

# filter only downregulated hits from decoupleR
decoupleR_down <- decoupleR_all[which(decoupleR_all$score < 0),]

# helper function definition
match_homer_to_decoupler <- function(
    homer_df,
    decoupler_df,
    motif_col = "Motif.Name",
    source_col = "tf"
) {
  motif_values <- homer_df[[motif_col]]
  
  source_values <- decoupler_df[[source_col]] %>%
    as.character() %>%
    unique()
  
  # Remove AR before matching because it frequently creates false matches
  source_values <- source_values[
    !is.na(source_values) &
      tolower(source_values) != "ar"
  ]
  
  matches <- crossing(
    motif = unique(motif_values),
    tf = source_values
  ) %>%
    filter(
      !is.na(motif),
      str_detect(
        tolower(motif),
        fixed(tolower(tf))
      )
    )
  
  matched_df <- homer_df %>%
    left_join(
      matches,
      by = setNames("motif", motif_col)
    ) %>%
    mutate(
      tf = coalesce(
        tf,
        str_trim(sub("\\(.*", "", .data[[motif_col]]))
      ),
      match_type = if_else(
        is.na(match(.data[[motif_col]], matches$motif)),
        "HOMER_name",
        "decoupleR_match"
      )
    )
  
  print(paste(sum(matched_df$match_type == "decoupleR_match"), "names matched to dorothea database,", sum(matched_df$match_type == "HOMER_name"), "names didn't match"))
  
  return(matched_df)
}

# load dorothea. omnipath / collectri doesn't work, maybe need new version of R
data("dorothea_mm", package = "dorothea")

# filter higher confidence targets
dorothea_net <- dorothea_mm %>%
  filter(confidence %in% c("A", "B", "C")) %>%
  dplyr::select(tf, target, mor)

# run matching annotations
ATAC_HOMER_TFnamed <- match_homer_to_decoupler(ATAC_HOMER, dorothea_net)
DMRichR_HOMER_TFnamed <- match_homer_to_decoupler(DMRichR_HOMER, dorothea_net)

# make a rank aggregation score to combine ranks across 3 datasets --------

# standardize columns, remove duplicates and keep the highest ranked one
ATAC_HOMER_ranked <- ATAC_HOMER_TFnamed %>%
  arrange(P.value) %>%
  mutate(rank = row_number()) %>%
  dplyr::select(tf, pval = P.value, rank) %>%
  distinct(tf, .keep_all = TRUE)

DMRichR_HOMER_ranked <- DMRichR_HOMER_TFnamed %>%
  arrange(P.value) %>%
  mutate(rank = row_number()) %>%
  dplyr::select(tf, pval = P.value, rank) %>%
  distinct(tf, .keep_all = TRUE)

decoupleR_ranked <- decoupleR_down %>%
  arrange(p_value) %>%
  mutate(rank = row_number()) %>%
  dplyr::select(tf = source, pval = p_value, rank) %>%
  distinct(tf, .keep_all = TRUE)

# aggregate ranks, get rra score
rra_all3 <- aggregateRanks(
  glist = list(
    ATAC_HOMER_ranked$tf, 
    DMRichR_HOMER_ranked$tf, 
    decoupleR_ranked$tf),
  method = "RRA"
)

# run another one with just EM and ATAC
# NOTE that rra_merge2 is using dorothea annotations, which will collapse some HOMER annotations into one TF name
rra_merge2 <- aggregateRanks(
  glist = list(
    ATAC_HOMER_ranked$tf, 
    DMRichR_HOMER_ranked$tf), 
  method = "RRA"
)

# append individual pval columns
rra_all3 <- rra_all3 %>%
  rename(
    tf = Name,
    rra_score = Score
  ) %>%
  left_join(
    ATAC_HOMER_ranked %>%
      dplyr::select(tf, pval_ATAC = pval),
    by = "tf"
  ) %>%
  left_join(
    DMRichR_HOMER_ranked %>%
      dplyr::select(tf, pval_DMRichR = pval),
    by = "tf"
  ) %>%
  left_join(
    decoupleR_ranked %>%
      dplyr::select(tf, pval_decoupleR = pval),
    by = "tf"
  ) %>%
  mutate(
    n_present = rowSums(
      !is.na(across(starts_with("pval_")))
    )
  ) %>%
  arrange(rra_score)

rra_merge2 <- rra_merge2 %>%
  rename(
    tf = Name,
    rra_score = Score
  ) %>%
  left_join(
    ATAC_HOMER_ranked %>%
      dplyr::select(tf, pval_ATAC = pval),
    by = "tf"
  ) %>%
  left_join(
    DMRichR_HOMER_ranked %>%
      dplyr::select(tf, pval_DMRichR = pval),
    by = "tf"
  ) %>%
  arrange(rra_score)

# merge ATAC and EM HOMER along raw HOMER labels -------------------------------------------------

# filter only q < 0.05 results
ATAC_HOMER <- ATAC_HOMER[which(ATAC_HOMER$q.value..Benjamini. < 0.05),]
DMRichR_HOMER <- DMRichR_HOMER[which(DMRichR_HOMER$q.value..Benjamini. < 0.05),]
decoupleR_down <- decoupleR_down[which(decoupleR_down$padj < 0.05),]

# merge overlapping TF hits between HOMER analysis for atac and EM seq (DMRichR)
HOMER_merge2 <- merge(ATAC_HOMER_TFnamed, DMRichR_HOMER, by="Consensus", suffixes=c(".ATAC", ".EM"))
HOMER_merge2 <- subset(HOMER_merge2, select=-Motif.Name.EM)
HOMER_merge2 <- HOMER_merge2[HOMER_merge2$P.value.ATAC < 0.05 & HOMER_merge2$P.value.EM < 0.05,]

# sort along the best pvalue in either column. this means hits towards the top may or may not drastically disagree between the two datasets
HOMER_merge2 <- HOMER_merge2 %>%
  mutate(best_p = pmax(P.value.ATAC, P.value.EM, na.rm = TRUE)) %>%
  arrange(best_p) %>%
  distinct(Consensus, .keep_all = TRUE) %>%
  dplyr::select(-best_p)

# filter only q < 0.05 HOMER rra_merge2 hits, for venn counting
rra_merge2_sig <- rra_merge2 %>%
  inner_join(HOMER_merge2, by = "tf")

# make TF predictions venn diagram -------------------------------------------------------

# count all exact intersections with RNA dataset
HOMER_ATAC_decoupleR_consensus <- ATAC_HOMER_TFnamed %>%
  filter(tf %in% decoupleR_down$source)

HOMER_DMRichR_decoupleR_consensus <- DMRichR_HOMER_TFnamed %>%
  filter(tf %in% decoupleR_down$source)

# all 3
HOMER_all3_consensus <- HOMER_merge2 %>%
  distinct(tf, .keep_all = TRUE) %>%
  filter(tf %in% decoupleR_down$source)

# directional disagreements with decoupleR
HOMER_all3_disagree <- HOMER_merge2 %>%
  distinct(tf, .keep_all = TRUE) %>%
  filter(tf %in% decoupleR_all$source) %>%
  setdiff(HOMER_all3_consensus)

# counting using deduplicated TFnamed (mix of decoupler and simplified HOMER names)
manual_counts <- c(
  "ATAC"         = nrow(ATAC_HOMER_ranked),  # ATAC only
  "EM"           = nrow(DMRichR_HOMER_ranked),   # EM only
  "RNA"          = nrow(decoupleR_all),  # RNA only
  "ATAC&EM"      = nrow(rra_merge2_sig),   # ATAC + EM, but not RNA
  "ATAC&RNA"     = nrow(HOMER_ATAC_decoupleR_consensus),   # ATAC + RNA, but not EM
  "EM&RNA"       = nrow(HOMER_DMRichR_decoupleR_consensus),   # EM + RNA, but not ATAC
  "ATAC&EM&RNA"  = nrow(HOMER_all3_consensus)    # shared by all three
)

homer_venn <- euler(manual_counts)

pdf(file.path(OUTPUT_DIR, "venn_TF_predictions.pdf"), width = 6, height = 6)

plot(
  homer_venn,
  quantities = TRUE,
  labels = TRUE,
  fills = c("#E69F00", "#56B4E9", "brown1"),
  edges = list(lwd = 1),
  main = list(
    label = "Overlapping TF predictions (HOMER and decoupleR, p < 0.05)",
    cex = 1)
)

dev.off()

# save all merged comparisons ---------------------------------------------

merged_results <- list(rra_all3,
                       rra_merge2)
merged_names <- list("rra_all3.tsv",
                     "rra_merge2.tsv"
)

for(i in seq_along(merged_results)){
  write.table(merged_results[[i]], file=file.path(OUTPUT_DIR, merged_names[[i]]), sep="\t", quote=FALSE, row.names=FALSE)
}

# export neat copy-pastable TF names from all3
write(HOMER_all3_consensus$tf, file = file.path(OUTPUT_DIR, "all3_TFnames.txt"))

rra_nominal <- rra_all3 %>%
  filter(rra_score < 0.05)

rra_nominal$tf %>%
  write(file = file.path(OUTPUT_DIR, "all3_sigRRA_TFnames.txt"))

# save RData for app
save(HOMER_merge2,
     HOMER_ATAC_decoupleR_consensus,
     HOMER_DMRichR_decoupleR_consensus,
     HOMER_all3_consensus,
     rra_all3,
     rra_merge2,
     homer_venn,
     file = file.path(OUTPUT_DIR, "HOMER_merged.RData"))

# print summary message ---------------------------------------------------

write(
  paste(
    "There were ", nrow(HOMER_all3_consensus),
    " explicitly shared TFs across all datasets: ",
    paste(HOMER_all3_consensus$tf, collapse = ", "),
    ".\n\nThere were ", nrow(rra_nominal),
    " TFs with nominal (unadjusted) RRA scores < 0.05: ",
    paste(rra_nominal$tf, collapse = ", "),
    ".\n\nThere were ", length(HOMER_all3_disagree),
    " TFs which were significant, but ATAC-EM and decoupleR predicted opposite directionality: ",
    paste(HOMER_all3_disagree$tf, collapse = ", ")
  ),
  file = file.path(OUTPUT_DIR, "summary_msg.txt")
)

message(
  "There were ", nrow(HOMER_all3_consensus),
  " explicitly shared, silenced TFs across all datasets: ",
  paste(HOMER_all3_consensus$tf, collapse = ", "),
  ".\n\nThere were ", nrow(rra_nominal),
  " TFs with nominal (unadjusted) RRA scores < 0.05: ",
  paste(rra_nominal$tf, collapse = ", "),
  ".\n\nThere were ", length(HOMER_all3_disagree),
  " TFs which were significant, but ATAC-EM and decoupleR predicted opposite directionality: ",
  paste(HOMER_all3_disagree$tf, collapse = ", ")
)



# ==========================================================================
# global.R
# made by sheryl 06.08.2026, vibe-coded with claude
#
# Loaded once when the app starts, before ui.R and server.R. Holds package
# library() calls, file paths, and all data loading, so ui.R/server.R stay
# focused on UI layout and reactive logic.
#
# With this file layout (global.R + ui.R + server.R, no app.R) and an R/
# folder next to them, Shiny (>= 1.5.0) automatically sources every .R file
# in R/ first, then global.R, then launches ui.R + server.R as the app.
# ==========================================================================

library(shiny)
library(DT)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(plotly)
library(eulerr)
library(rGREAT)   # needed by R/great_report_module.R

# --- paths - check these ---------------------------------------------------

DATA_PATH  <- "V:/MARIA/Sheryl/integration_results/rGREAT_results"
PLOTS_PATH <- "V:/MARIA/Sheryl/integration_results/RNA_ATAC_integration_plots"
PLOTS_FILE <- "rna_atac_integration_plots.RData"

DECOUPLER_PATH <- "V:/MARIA/Sheryl/integration_results/virtual_KO"
DECOUPLER_FILE <- "decoupleR_results.RData"

# swap these labels due to backwards DMRichR contrast ("hyper" is hypermethylation in the WT)
HOMER_EM_HYPER_PATH <- "V:/MARIA/Sheryl/DMRichR_results/260724_v1.3-loosecutoff/HOMER/hypo"
HOMER_EM_HYPO_PATH <- "V:/MARIA/Sheryl/DMRichR_results/260724_v1.3-loosecutoff/HOMER/hyper"
HOMER_ATAC_PATH <- "V:/MARIA/Sheryl/atac_seq_results/E14.5/HOMER/homer_filtered_unique_WT_1bp.mm39.bed"
HOMER_MERGED_PATH <- "V:/MARIA/Sheryl/integration_results/HOMER_merged_new/HOMER_merged.RData" # -> HOMER_merge2, HOMER_ATAC_decoupleR_consensus, HOMER_DMRichR_decoupleR_consensus, rra_all3, rra_merge2, homer_venn
HOMER_FILE <- "knownResults.html"

ATAC_TABLE_PATH <- "V:/MARIA/Sheryl/atac_seq_results/E14.5/diffbind"
ATAC_TABLE_FILE <- "all_ATAC_diffgenes.E14.5.tsv"
ATAC_GENE_COL   <- "SYMBOL"

# --- HOMER results -----------------------------------------------------------

# serve the report folders as static files so the .html reports can be
# embedded via <iframe>
addResourcePath("homer_em_hyper_path", HOMER_EM_HYPER_PATH)
addResourcePath("homer_em_hypo_path", HOMER_EM_HYPO_PATH)
addResourcePath("homer_atac_path", HOMER_ATAC_PATH)

load(HOMER_MERGED_PATH)

# move most relevant columns of HOMER dfs to the front
HOMER_merge2 <- HOMER_merge2 %>% 
  relocate(tf)

HOMER_ATAC_decoupleR_consensus <- HOMER_ATAC_decoupleR_consensus %>%
  relocate(tf)

HOMER_DMRichR_decoupleR_consensus <- HOMER_DMRichR_decoupleR_consensus %>%
  relocate(tf)

HOMER_all3_consensus <- HOMER_all3_consensus %>%
  relocate(tf)

# pre-prepared pairwise intersection tables for the "Intersection tables" tab
intersection_tables <- list(
  "Rank-aggregated all-3 overlap" = rra_all3,
  "Rank-aggregated ATAC & EM overlap" = rra_merge2,
  "ATAC & EM exact overlap"  = HOMER_merge2,
  "ATAC & RNA exact overlap" = HOMER_ATAC_decoupleR_consensus,
  "EM & RNA exact overlap"   = HOMER_DMRichR_decoupleR_consensus,
  "All-3 exact overlap"   = HOMER_all3_consensus
)

# --- GREAT results -------------------------------------------------------

load(file.path(DATA_PATH, "rGREAT_results_EM_promoters.RData"))
load(file.path(DATA_PATH, "rGREAT_results_ATAC_promoters.RData"))
load(file.path(DATA_PATH, "rGREAT_intersection_results_promoters.RData"))

# static, pre-extracted tables for the filterable "GREAT Results" tab
results <- list(
  GO_BP_intersect = intersection_results$GO_BP,
  GO_MF_intersect = intersection_results$GO_MF,
  GO_CC_intersect = intersection_results$GO_CC,
  GO_KEGG_intersect = intersection_results$KEGG
)

# the raw GreatObject objects, kept separately from `results` because the
# interactive report module (R/great_report_module.R) needs the object
# itself, not just its extracted @table. Only the primary per-ontology
# objects are "explorable" this way - the pre-merged intersection tables
# in `results` are plain data frames, not GreatObjects, so they're left out.
great_objects <- list(
  "ATAC - GO Biological Process" = go.bp_great_atac,
  "ATAC - GO Molecular Function" = go.mf_great_atac,
  "ATAC - GO Cellular Component" = go.cc_great_atac,
  "ATAC - KEGG"                  = kegg_great_atac,
  "EM - GO Biological Process"   = go.bp_great_em,
  "EM - GO Molecular Function"   = go.mf_great_em,
  "EM - GO Cellular Component"   = go.cc_great_em
)

# possible p-value column names across all GREAT result types, with readable labels
p_col_candidates <- c(
  "Both"          = "p_adjust",
  "Hyper"         = "p_adjust_hyper",
  "Hypo"          = "p_adjust_hypo",
  "Both (ATAC)"   = "p_adjust_ATAC",
  "Hyper (ATAC)"  = "p_adjust_hyper_ATAC",
  "Hypo (ATAC)"   = "p_adjust_hypo_ATAC",
  "Both (EM)"     = "p_adjust_EM",
  "Hyper (EM)"    = "p_adjust_hyper_EM",
  "Hypo (EM)"     = "p_adjust_hypo_EM"
)

# --- plots from the other analysis ------------------------------------------
# Loads whatever objects are in the .RData file into a temporary environment,
# then keeps only ggplot objects ("gg") or base-R recorded plots
# ("recordedplot", made with recordPlot()) and puts them in a named list.
#
# EASIEST OPTION: if you control that RData file, just save a single named
# list called `plot_list` in it, e.g.:
#     plot_list <- list(Volcano = p1, PCA = p2, Heatmap = p3)
#     save(plot_list, file = "other_analysis_plots.RData")
# then the code below will just pick it up directly (see the check at the end).

plot_env <- new.env()
load(file.path(PLOTS_PATH, PLOTS_FILE), envir = plot_env)

loaded_objs <- mget(ls(plot_env), envir = plot_env)

if (!is.null(loaded_objs$plot_list) && is.list(loaded_objs$plot_list)) {
  # a ready-made named list was found - use it as-is
  plot_list <- loaded_objs$plot_list
} else {
  # otherwise auto-detect individual plot objects and name them after their variable names
  plot_list <- Filter(function(x) inherits(x, "gg") || inherits(x, "recordedplot"), loaded_objs)
}

if (length(plot_list) == 0) {
  warning("No plot objects (ggplot or recordedplot) were found in ", PLOTS_FILE)
}

# --- decoupleR TF prediction results ----------------------------------------
# Expects an .RData file containing: dorothea_net, rna_unique, ulm_all, ulm_DE
# (saved from your analysis script - see the save() call to add there)

load(file.path(DECOUPLER_PATH, DECOUPLER_FILE))

# move most relevant columns to the front
rna_unique <- rna_unique %>%
  relocate(log2FoldChange) %>%
  relocate(padj) %>%
  relocate(Gene_Symbol)

ulm_results <- list(
  "All transcripts" = ulm_all,
  "DE genes only"    = ulm_DE
)

# --- ATAC-seq results table --------------------------------------------------
# Expects an .RData file containing a data frame with ATAC-seq differential
# results, including a gene-name column (see ATAC_GENE_COL above).
#
# If the file contains exactly one data frame, that's used automatically.
# Otherwise, save it as an object literally named `atac_df`, e.g.:
#     save(atac_df, file = "atac_seq_results.RData")

atac_env <- new.env()
atac_env$atac_df <- read.table(file.path(ATAC_TABLE_PATH, ATAC_TABLE_FILE))

# move most relevant columns of ATAC df to the front
atac_env$atac_df <- atac_env$atac_df %>%
  relocate(Fold) %>%
  relocate(FDR) %>%
  relocate(SYMBOL)

atac_loaded_objs <- mget(ls(atac_env), envir = atac_env)

if (!is.null(atac_loaded_objs$atac_df)) {
  atac_df <- atac_loaded_objs$atac_df
} else {
  atac_df_candidates <- Filter(is.data.frame, atac_loaded_objs)
  if (length(atac_df_candidates) == 1) {
    atac_df <- atac_df_candidates[[1]]
  } else if (length(atac_df_candidates) > 1) {
    stop(
      "Multiple data frames found in ", ATAC_TABLE_FILE,
      " - save just one as `atac_df` so the app knows which to use."
    )
  } else {
    stop("No data frame found in ", ATAC_TABLE_FILE)
  }
}

if (!ATAC_GENE_COL %in% names(atac_df)) {
  warning(
    "ATAC_GENE_COL ('", ATAC_GENE_COL, "') not found in atac_df columns: ",
    paste(names(atac_df), collapse = ", ")
  )
}

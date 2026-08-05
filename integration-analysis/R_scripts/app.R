library(shiny)
library(DT)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(plotly)

# global variables - check these ------------------------------------------

DATA_PATH  <- "V:/MARIA/Sheryl/integration_results/rGREAT_results"
PLOTS_PATH <- "V:/MARIA/Sheryl/integration_results/RNA_ATAC_integration_plots"
PLOTS_FILE <- "rna_atac_integration_plots.RData"          

DECOUPLER_PATH <- "V:/MARIA/Sheryl/integration_results/RNA_decoupleR_results"
DECOUPLER_FILE <- "decoupleR_results.RData"

HOMER_EM_HYPER_PATH <- "V:/MARIA/Sheryl/DMRichR_results/260724_v1.3-loosecutoff/HOMER/hyper"
HOMER_EM_HYPO_PATH <- "V:/MARIA/Sheryl/DMRichR_results/260724_v1.3-loosecutoff/HOMER/hypo"
HOMER_ATAC_PATH <- "V:/MARIA/Sheryl/atac_seq_results/E14.5/HOMER/homer_filtered_unique_WT_1bp.mm39.bed"
HOMER_MERGED_PATH <- "V:/MARIA/Sheryl/integration_results/HOMER_merged/HOMER_merged.RData" # -> HOMER_merge2, HOMER_ATAC_decoupleR_consensus, HOMER_DMRichR_decoupleR_consensus, HOMER_all3_consensus_down, homer_venn
HOMER_FILE <- "knownResults.html"

ATAC_TABLE_PATH <- "V:/MARIA/Sheryl/atac_seq_results/E14.5/diffbind"
ATAC_TABLE_FILE <- "all_ATAC_diffgenes.E14.5.tsv"
ATAC_GENE_COL   <- "SYMBOL"

# load HOMER results ------------------------------------------------------

# serve the report folder as static files so the .html can be embedded via <iframe>
addResourcePath("homer_em_hyper_path", HOMER_EM_HYPER_PATH)
addResourcePath("homer_em_hypo_path", HOMER_EM_HYPO_PATH)
addResourcePath("homer_atac_path", HOMER_ATAC_PATH)

load(HOMER_MERGED_PATH)

# load GREAT results ------------------------------------------------------

load(file.path(DATA_PATH, "rGREAT_results_EM_promoters.RData"))
load(file.path(DATA_PATH, "rGREAT_results_ATAC_promoters.RData"))
load(file.path(DATA_PATH, "rGREAT_intersection_results_promoters.RData"))

results <- list(
  ATAC_GO_BP = go.bp_great_atac@table,
  ATAC_GO_MF = go.mf_great_atac@table,
  ATAC_GO_CC = go.cc_great_atac@table,
  ATAC_KEGG = kegg_great_atac@table,
  EM_GO_BP   = go.bp_great_em@table,
  EM_GO_MF   = go.mf_great_em@table,
  EM_GO_CC   = go.cc_great_em@table,
  GO_BP_intersect = intersection_results$GO_BP,
  GO_MF_intersect = intersection_results$GO_MF,
  GO_CC_intersect = intersection_results$GO_CC,
  GO_KEGG_intersect = intersection_results$KEGG
)

# --- Load plots from the other analysis -----------------------------------
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

# --- Load decoupleR TF prediction results ----------------------------------
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

# --- Load ATAC-seq results table --------------------------------------------
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

# ui ----------------------------------------------------------------------

ui <- fluidPage(
  titlePanel("E14.5 PPC Multiomics Explorer (Summer 2026)"),
  
  tabsetPanel(
    tabPanel(
      "GREAT Results",
      br(),
      selectInput(
        "result_name",
        "Result",
        choices = names(results)
      ),
      
      selectInput(
        "p_column",
        "Adjusted p-value",
        choices = character(0)
      ),
      
      numericInput(
        "cutoff",
        "Cutoff",
        value = 0.05,
        min = 0,
        max = 1,
        step = 0.01
      ),
      
      downloadButton("download", "Download CSV"),
      
      br(),
      br(),
      
      DTOutput("table")
    ),
    
    tabPanel(
      "RNA ATAC integrated plots",
      br(),
      selectInput(
        "plot_choice",
        "Select Plot",
        choices = names(plot_list)
      ),
      plotOutput("selected_plot", height = "600px")
    ),
    
    tabPanel(
      "decoupleR TF prediction (RNA-seq)",
      br(),
      
      h4("decoupleR TF prediction - volcano plot"),
      selectInput(
        "ulm_choice",
        "Gene set",
        choices = names(ulm_results)
      ),
      plotlyOutput("volcano_plot", height = "600px"),
      
      hr(),
      
      fluidRow(
        column(
          6,
          h4("Look up a TF"),
          selectizeInput(
            "tf_choice",
            "Transcription factor",
            choices = sort(unique(dorothea_net$tf)),
            options = list(placeholder = "Start typing a TF name...")
          ),
          strong("Predicted activity (from selected gene set above):"),
          verbatimTextOutput("tf_ulm_summary"),
          strong("Target genes driving this TF's score:"),
          DTOutput("tf_targets_table")
        ),
        
        column(
          6,
          h4("Look up a gene"),
          textInput(
            "gene_choice",
            "Gene",
            placeholder = "Type a gene name, e.g. Cox11"
          ),
          strong("TFs predicted to regulate this gene:"),
          DTOutput("gene_tfs_table")
        )
      ),
      
      hr(),
      
      h4("RNA-seq results"),
      textInput(
        "rnaseq_gene_search",
        "Search gene",
        placeholder = "Type a gene name, e.g. Cox11"
      ),
      DTOutput("rnaseq_table"),
      
      hr(),
      
      h4("ATAC-seq results"),
      textInput(
        "atacseq_gene_search",
        "Search gene",
        placeholder = "Type a gene name, e.g. Cox11"
      ),
      DTOutput("atacseq_table")
    ),
    
    tabPanel(
      "HOMER TF predictions (ATAC, EM) and all 3 merged",
      br(),
      selectInput(
        "homer_choice",
        "Dataset",
        choices = c("EM (hyper in WT)" = "homer_em_hyper_path", "EM (hypo in KO)" = "homer_em_hypo_path", "ATAC" = "homer_atac_path", "EM and ATAC overlap (silenced in KO)" = "homer_merge2", "All 3 merged" = "homer_all3")
      ),
      uiOutput("homer_output")
    )
  )
)

# server ------------------------------------------------------------------

server <- function(input, output, session) {
  
  # possible p-value column names across all result types, with readable labels
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
  
  # update the p_column dropdown to only show columns that exist in the
  # currently selected result data frame (fires on load and whenever
  # result_name changes)
  observe({
    df <- results[[input$result_name]]
    available <- p_col_candidates[p_col_candidates %in% names(df)]
    
    updateSelectInput(
      session,
      "p_column",
      choices = available
    )
  })
  
  filtered_results <- reactive({
    df <- results[[input$result_name]]
    
    req(input$p_column %in% names(df))
    
    df %>%
      filter(.data[[input$p_column]] < input$cutoff) %>%
      arrange(.data[[input$p_column]])
  })
  
  output$table <- renderDT({
    filtered_results()
  })
  
  output$download <- downloadHandler(
    filename = function() {
      paste0(
        input$result_name,
        "_",
        input$p_column,
        "_filtered.csv"
      )
    },
    content = function(file) {
      write.csv(
        filtered_results(),
        file,
        row.names = FALSE
      )
    }
  )
  
  output$selected_plot_plotly <- renderPlotly({
    req(input$plot_choice)
    p <- plot_list[[input$plot_choice]]
    req(inherits(p, "gg"))
    
    ggplotly(p) %>%
      layout(dragmode = "zoom")
  })
  
  output$selected_plot <- renderPlot({
    req(input$plot_choice)
    p <- plot_list[[input$plot_choice]]
    
    if (inherits(p, "recordedplot")) {
      replayPlot(p)
    } else {
      print(p)
    }
  })
  
  # --- decoupleR TF prediction tab ------------------------------------------------------
  
  output$volcano_plot <- renderPlotly({
    req(input$ulm_choice)
    df <- ulm_results[[input$ulm_choice]]
    
    p <- ggplot(data = df, aes(x = score, y = -log10(p_value), col = p_value, label = source)) +
      geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = "dashed") +
      geom_hline(yintercept = -log10(0.05), col = "gray", linetype = "dashed") +
      geom_point() +
      geom_text(nudge_y=0.02, size = 3, show.legend = FALSE) +
      labs(
        y = "-log10(pval)",
        x = "score",
        title = paste("decoupleR predicted TFs -", input$ulm_choice)
      ) +
      scale_colour_gradient(low = "dodgerblue", high = "black") +
      ylim(c(0.5, NA))
    
    ggplotly(p, tooltip = c("label", "x", "y")) %>%
      layout(dragmode = "zoom")
  })
  
  # target genes + contribution for the selected TF
  tf_targets_reactive <- reactive({
    req(input$tf_choice)
    
    dorothea_net %>%
      filter(tf == input$tf_choice) %>%
      left_join(
        rna_unique %>% select(Gene_Symbol, stat),
        by = c("target" = "Gene_Symbol")
      ) %>%
      mutate(contribution = mor * stat) %>%
      arrange(desc(abs(contribution)))
  })
  
  output$tf_targets_table <- renderDT({
    tf_targets_reactive()
  })
  
  # the selected TF's own predicted activity score, from whichever gene set is chosen above
  output$tf_ulm_summary <- renderPrint({
    req(input$tf_choice, input$ulm_choice)
    
    df <- ulm_results[[input$ulm_choice]]
    hit <- df[df$source == input$tf_choice, ]
    
    if (nrow(hit) == 0) {
      cat("No result for this TF in the selected gene set (may not meet minsize threshold).")
    } else {
      print(hit)
    }
  })
  
  # TFs predicted to regulate the selected gene
  # debounce so we're not re-filtering on every single keystroke
  gene_choice_debounced <- reactive({
    trimws(input$gene_choice)
  }) %>% debounce(400)
  
  gene_tfs_reactive <- reactive({
    gene_input <- gene_choice_debounced()
    req(nzchar(gene_input))
    
    # case-insensitive match against the actual gene symbol
    matched_gene <- rna_unique$Gene_Symbol[
      tolower(rna_unique$Gene_Symbol) == tolower(gene_input)
    ]
    
    validate(
      need(length(matched_gene) > 0, paste0("No gene found matching '", gene_input, "'"))
    )
    
    dorothea_net %>%
      filter(target == matched_gene[1]) %>%
      left_join(
        rna_unique %>% select(Gene_Symbol, stat),
        by = c("target" = "Gene_Symbol")
      ) %>%
      mutate(contribution = mor * stat) %>%
      arrange(desc(abs(contribution)))
  })
  
  output$gene_tfs_table <- renderDT({
    gene_tfs_reactive()
  })
  
  # --- RNA-seq / ATAC-seq tables at bottom of decoupleR TF prediction tab --------------
  
  # RNA-seq table: uses rna_unique, already loaded from the decoupleR results file
  rnaseq_search_debounced <- reactive({
    trimws(input$rnaseq_gene_search)
  }) %>% debounce(400)
  
  output$rnaseq_table <- renderDT({
    search_term <- rnaseq_search_debounced()
    
    if (nzchar(search_term)) {
      rna_unique %>%
        filter(grepl(search_term, Gene_Symbol, ignore.case = TRUE))
    } else {
      rna_unique
    }
  })
  
  # ATAC-seq table: uses atac_df, loaded from the ATAC-seq results file above
  atacseq_search_debounced <- reactive({
    trimws(input$atacseq_gene_search)
  }) %>% debounce(400)
  
  output$atacseq_table <- renderDT({
    search_term <- atacseq_search_debounced()
    
    if (nzchar(search_term) && ATAC_GENE_COL %in% names(atac_df)) {
      atac_df %>%
        filter(grepl(search_term, .data[[ATAC_GENE_COL]], ignore.case = TRUE))
    } else {
      atac_df
    }
  })
  
  # --- HOMER tab --------------------------------------------------------
  
  output$homer_output <- renderUI({
    req(input$homer_choice)
    
    if (input$homer_choice == "homer_merge2") {
      
      DTOutput("homer_table_merge2")
      
    } else if (input$homer_choice == "homer_all3") {
      
      tagList(
        plotOutput("homer_venn_plot", height = "600px"),
        br(),
        DTOutput("homer_table_all3")
      )
      
    } else {
      
      tags$iframe(
        src = file.path(input$homer_choice, HOMER_FILE),
        width = "100%",
        height = "900px",
        frameborder = 0,
        style = "border: none;"
      )
    }
  })
  
  output$homer_table_merge2 <- renderDT({
    req(input$homer_choice == "homer_merge2")
    
    datatable(
      HOMER_merge2,
      filter = "top",
      rownames = FALSE,
      options = list(
        pageLength = 25,
        scrollX = TRUE
      )
    )
  })
  
  output$homer_table_all3 <- renderDT({
    req(input$homer_choice == "homer_all3")
    
    datatable(
      HOMER_all3_consensus_down,
      filter = "top",
      rownames = FALSE,
      options = list(
        pageLength = 25,
        scrollX = TRUE
      )
    )
  })
  
  output$homer_venn_plot <- renderPlot({
    req(input$homer_choice == "homer_all3")
    
    plot(
      homer_venn,
      quantities = TRUE,
      labels = TRUE
    )
  })
}

shinyApp(ui, server)
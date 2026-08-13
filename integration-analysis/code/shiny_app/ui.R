ui <- fluidPage(
  titlePanel("E14.5 PPC Multiomics Explorer (Summer 2026)"),
  
  tabsetPanel(
    
    tabPanel(
      "GREAT Interactive Report",
      br(),
      selectInput(
        "great_report_choice",
        "Result",
        choices = names(great_objects)
      ),
      hr(),
      greatReportUI("great_report")
    ),
    
    tabPanel(
      "GREAT Intersections",
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
      "HOMER TF predictions (ATAC, EM)",
      br(),
      selectInput(
        "homer_choice",
        "Dataset",
        choices = c(
          "EM (hyper, silenced in KO)" = "homer_em_hyper_path",
          "EM (hypo, accessible in KO)"  = "homer_em_hypo_path",
          "ATAC"             = "homer_atac_path"
        )
      ),
      uiOutput("homer_output")
    ),
    
    tabPanel(
      
      "TF predictions, all 3 datasets",
      br(),
      
      h4("Venn diagram: exact overlap of TF predictions across EM, ATAC, and RNA-seq"),
      plotOutput("homer_venn_plot", height = "600px"),
      br(),
      
      h4("Overlapping TF predictions tables: Robust Rank-Aggregated (RRA) overlaps and exact overlaps"),
      h5("RRA scores weigh the TF prediction rankings according to p-values across all datasets. Exact overlaps include hits whose names explicitly match between datasets. Some hits may be absent from RNA results due to different annotation databases (ATAC and EM use HOMER motif names, whereas RNA uses Dorothea names)."),
      
      selectInput(
        "intersection_choice",
        "Comparison",
        choices = names(intersection_tables)
      ),
      
      DTOutput("intersection_table"),
    ),
  )
)

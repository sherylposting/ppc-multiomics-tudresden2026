# ==========================================================================
# great_report_module.R
#
# A Shiny module version of rGREAT:::shinyReport(<GreatObject>). The
# original function builds and launches its own standalone shinyApp(),
# which can't be nested inside another running app. This reimplements the
# same UI/behavior as a namespaced module (greatReportUI / greatReportServer)
# so it can live inside a tabPanel() of a larger app instead.
#
# Files in R/ are auto-sourced by Shiny before ui.R/server.R run, as long as
# there's no app.R in the same directory (Shiny >= 1.5.0).
#
# Usage in ui.R:
#   greatReportUI("great_report")
#
# Usage in server.R:
#   greatReportServer("great_report", great_object = reactive({
#     req(input$great_report_choice)
#     great_objects[[input$great_report_choice]]
#   }))
#
# `great_object` must be a reactive() that returns a single GreatObject
# (as produced by rGREAT::great()).
# ==========================================================================

greatReportUI <- function(id) {
  ns <- NS(id)

  tagList(
    tags$style(HTML(paste0(
      "#", id, " pre { width: 800px; padding: 20px; }
       #", id, " .global_control td { padding-right: 20px; }
       #", id, " .fake_link { color: #337ab7; text-decoration: none; cursor: pointer; }
       #", id, " .fake_link:hover { text-decoration: underline; }
       #", id, " .error { color: red; }
       #", id, " .message { padding: 20px 0; }"
    ))),

    div(id = id,

      h4("Job description"),
      verbatimTextOutput(ns("job_desc")),

      h4("Global region-gene associations"),
      tags$pre("plotRegionGeneAssociations(object)"),
      plotOutput(ns("global_plot"), width = "1000px", height = "400px"),

      hr(),
      h4("Controls"),
      div(class = "global_control",
        tags$table(tags$tr(
          tags$td(textInput(ns("padj_cutoff"), "Cutoff for adjusted p-values (Binomial test)", value = "0.05", width = 380)),
          tags$td(textInput(ns("observed_hits_cutoff"), "Cutoff for observed region hits", value = "5", width = 400))
        ))
      ),
      hr(),

      tabsetPanel(type = "tabs",
        tabPanel("Enrichment table",
          htmlOutput(ns("error")),
          htmlOutput(ns("enrichment_table"))
        ),
        tabPanel("Volcano plot",
          br(),
          tags$pre("plotVolcano(object)"),
          radioButtons(ns("volcano_x_values"), "Values on x-axis",
            c("Fold enrichment: log2(obs/exp)" = "fold_enrichment",
              "z-score: (obs-exp)/sd" = "z-score"),
            selected = "fold_enrichment", inline = TRUE),
          radioButtons(ns("volcano_y_values"), "Values on y-axis",
            c("Raw p-values" = "p_value",
              "Adjusted p-values" = "p_adjust"),
            selected = "p_value", inline = TRUE),
          plotOutput(ns("volcano_plot"), width = "600px", height = "600px")
        )
      )
    )
  )
}

greatReportServer <- function(id, great_object) {
  # great_object: reactive() returning the currently-selected GreatObject
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # builds the enrichment-table datatable, with a clickable term-name
    # column that pops up the region-gene association modal for that term
    format_table <- function(tb) {
      tb$id <- paste0(
        "<a class='fake_link' onclick=\"Shiny.onInputChange('", ns("select_term"),
        "', '');Shiny.onInputChange('", ns("select_term"), "', '", tb$id, "');false;\">",
        tb$id, "</a>"
      )

      has_desc <- "description" %in% colnames(tb)
      if (has_desc) {
        offset <- 1
        tb <- tb[, c("id", "description", "mean_tss_dist", "p_value", "p_adjust",
                     "fold_enrichment", "observed_region_hits", "genome_fraction",
                     "p_value_hyper", "p_adjust_hyper", "fold_enrichment_hyper",
                     "observed_gene_hits", "gene_set_size")]
        colnames(tb) <- c("Term Name", "Term Description", "Mean Abs Dist to TSS (bp)",
                           "Binom Raw P-value", "Binom Adjusted P-value", "Binom Fold Enrichment",
                           "Binom Observed Region Hits", "Genome Fraction",
                           "Hyper Raw P-value", "Hyper Adjusted P-value", "Hyper Fold Enrichment",
                           "Observed Gene Hits", "Total Genes in Gene Set")
      } else {
        offset <- 0
        tb <- tb[, c("id", "mean_tss_dist", "p_value", "p_adjust",
                     "fold_enrichment", "observed_region_hits", "genome_fraction",
                     "p_value_hyper", "p_adjust_hyper", "fold_enrichment_hyper",
                     "observed_gene_hits", "gene_set_size")]
        colnames(tb) <- c("Term Name", "Mean Abs Dist to TSS (bp)",
                           "Binom Raw P-value", "Binom Adjusted P-value", "Binom Fold Enrichment",
                           "Binom Observed Region Hits", "Genome Fraction",
                           "Hyper Raw P-value", "Hyper Adjusted P-value", "Hyper Fold Enrichment",
                           "Observed Gene Hits", "Total Genes in Gene Set")
      }

      dt <- datatable(
        tb, escape = FALSE, rownames = FALSE, selection = "none",
        width = "100%", height = "auto",
        options = list(
          searching = FALSE,
          rowCallback = JS(paste0(
            "function(row, data) {",
            "  $(this.api().cell(row, ", 2 + offset, ").node()).html(data[", 2 + offset, "].toExponential(3));",
            "  $(this.api().cell(row, ", 3 + offset, ").node()).html(data[", 3 + offset, "].toExponential(3));",
            "  $(this.api().cell(row, ", 7 + offset, ").node()).html(data[", 7 + offset, "].toExponential(3));",
            "  $(this.api().cell(row, ", 8 + offset, ").node()).html(data[", 8 + offset, "].toExponential(3));",
            "}"
          ))
        )
      )
      dt <- formatRound(dt, "Binom Fold Enrichment", 3)
      dt <- formatPercentage(dt, "Genome Fraction", 3)
      dt <- formatPercentage(dt, "Hyper Fold Enrichment", 3)
      dt
    }

    output$job_desc <- renderPrint({
      req(great_object())
      show(great_object())
      cat("\n")
      cat("Cutoff for adjusted p-values (Binomial test): ", input$padj_cutoff, "\n", sep = "")
      cat("Cutoff for observed region hits: ", input$observed_hits_cutoff, "\n", sep = "")
    })

    observe({
      req(great_object())
      object <- great_object()

      suppressWarnings(padj_cutoff <- as.numeric(input$padj_cutoff))
      suppressWarnings(observed_hits_cutoff <- as.numeric(input$observed_hits_cutoff))

      output$volcano_plot <- renderPlot({
        plotVolcano(object, min_region_hits = observed_hits_cutoff,
                    x_values = input$volcano_x_values, y_values = input$volcano_y_values)
      })

      if (is.na(padj_cutoff) || is.na(observed_hits_cutoff)) {
        output$error <- renderUI(HTML("<p class='error message'>Wrong format for cutoffs.</p>"))
        output$enrichment_table <- renderUI(HTML(""))
        return(NULL)
      }

      tb <- getEnrichmentTable(object, min_region_hits = observed_hits_cutoff)
      tb <- tb[tb$p_adjust <= padj_cutoff, , drop = FALSE]

      if (nrow(tb) == 0) {
        output$error <- renderUI(HTML(""))
        output$enrichment_table <- renderUI(HTML("<p class='message'>No significant term under current cutoffs.</p>"))
      } else {
        output$error <- renderUI(HTML(""))
        output$enrichment_table <- renderUI({
          div(
            HTML(paste0("<h4>Enrichment table (", nrow(tb), " significant terms)</h4>")),
            format_table(tb)
          )
        })
      }
    })

    # clicking a term name in the enrichment table opens a modal with the
    # region-gene association plot + table for that specific term
    observeEvent(input$select_term, {
      req(input$select_term, great_object())
      term <- input$select_term
      object <- great_object()

      tb <- getRegionGeneAssociations(object, term_id = term)
      tb <- as.data.frame(tb)
      colnames(tb) <- c("Chromosome", "Start", "End", "Width", "Strand", "Annotated Genes", "Distance to TSSs")
      tb <- tb[, -5]

      output$select_term_plot <- renderPlot({
        plotRegionGeneAssociations(object, term_id = term)
      })

      output$select_term_table <- DT::renderDT({
        datatable(tb, escape = FALSE, rownames = FALSE, selection = "none",
                  options = list(searching = FALSE))
      })

      showModal(modalDialog(
        title = paste0("Region-gene associations for term: ", term),
        tags$pre(paste0("plotRegionGeneAssociations(object, term_id = '", term, "')")),
        plotOutput(ns("select_term_plot"), width = "1000px", height = "400px"),
        hr(),
        tags$pre(paste0("getRegionGeneAssociations(object, term_id = '", term, "')")),
        DTOutput(ns("select_term_table")),
        easyClose = TRUE,
        size = "l"
      ))
    })

    output$global_plot <- renderPlot({
      req(great_object())
      plotRegionGeneAssociations(great_object())
    }, res = 100)
  })
}

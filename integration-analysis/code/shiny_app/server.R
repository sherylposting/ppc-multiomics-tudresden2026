server <- function(input, output, session) {
  
  # --- GREAT Results tab (static filterable tables) --------------------------
  
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
  
  # --- GREAT Interactive Report tab -------------------------------------------
  
  greatReportServer("great_report", great_object = reactive({
    req(input$great_report_choice)
    great_objects[[input$great_report_choice]]
  }))
  
  # --- RNA ATAC integrated plots tab ------------------------------------------
  
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
  
  # --- decoupleR TF prediction tab --------------------------------------------
  
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
  
  # --- HOMER tab (individual dataset HTML reports) ----------------------------
  
  output$homer_output <- renderUI({
    req(input$homer_choice)
    
    tags$iframe(
      src = file.path(input$homer_choice, HOMER_FILE),
      width = "100%",
      height = "800px",
      style = "border: none;"
    )
  })
  
  output$homer_venn_plot <- renderPlot({
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
  })
  
  intersection_data <- reactive({
    req(input$intersection_choice)
    intersection_tables[[input$intersection_choice]]
  })
  
  output$intersection_table <- renderDT({
    req(intersection_data())
    
    datatable(
      intersection_data(),
      filter = "top",
      rownames = FALSE,
      options = list(
        pageLength = 25,
        scrollX = TRUE
      )
    )
  })
  
  # --- TF predictions, all 3 datasets tab -------------------------------------
  
  output$rra_merge2_table <- renderDT({
    datatable(
      rra_merge2,
      filter = "top",
      rownames = FALSE,
      options = list(
        pageLength = 25,
        scrollX = TRUE
      )
    )
  })
  
  output$rra_all3_table <- renderDT({
    datatable(
      rra_all3,
      filter = "top",
      rownames = FALSE,
      options = list(
        pageLength = 25,
        scrollX = TRUE
      )
    )
  })
  
  # --- Intersection tables tab -------------------------------------------------
  
  intersection_data <- reactive({
    req(input$intersection_choice)
    intersection_tables[[input$intersection_choice]]
  })
  
  output$intersection_table <- renderDT({
    req(intersection_data())
    
    datatable(
      intersection_data(),
      filter = "top",
      rownames = FALSE,
      options = list(
        pageLength = 25,
        scrollX = TRUE
      )
    )
  })
}
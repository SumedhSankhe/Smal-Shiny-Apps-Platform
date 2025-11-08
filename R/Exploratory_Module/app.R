# Exploratory Module: Boxplots, Heatmap, PCA
#
# Purpose: Allow users to upload a sample metadata CSV (with at least column
#          sampleName and one additional covariate column) and an expression/sample matrix CSV where
#          rows are targets (features) and columns are sample names. Provides
#          interactive visualization: single-target boxplot (expression by covariate),
#          heatmap of top variable targets, and PCA colored by covariate.
#
# Contract:
# Inputs:
#   - metadata_file: CSV with columns: sampleName plus >=1 covariate column (user chooses which to use)
#   - matrix_file:   CSV, first column optional; if first column has header or implicit rownames
#                    each row = target, columns = sampleName. Will be read with row.names=1.
# Outputs (rendered):
#   - Boxplot (selected target vs chosen covariate)
#   - Heatmap (top N targets by variance across samples)
#   - PCA scatter (PC1 vs PC2, samples colored by chosen covariate)
# Error modes:
#   - Missing required columns
#   - Sample name mismatch between metadata and matrix
#   - Non-numeric matrix values (after coercion)
# Edge Cases handled:
#   - Duplicate sampleNames in metadata (deduplicated / flagged)
#   - Targets with zero variance (excluded from variance ranking)
#   - Fewer targets than requested top N (fallback to available)
#
# Dependencies: shiny, bslib, ggplot2
# (All imported explicitly via box::use). If ggplot2 not present, install via: renv::install("ggplot2", lock=TRUE)

box::use(
  shiny[NS, moduleServer, reactive, req, validate, need, renderPlot, renderUI,
        fileInput, selectInput, numericInput, plotOutput, uiOutput, h4, p, div, tags, tagList],
  bslib[layout_sidebar, sidebar, card, card_title, card_body, navset_card_tab, nav_panel],
  ggplot2[ggplot, aes, geom_boxplot, geom_point, theme_minimal, labs],
  rlang[sym],
  stats[prcomp],
  utils[read.csv],
  plotly[plotlyOutput, renderPlotly, plot_ly],
  shinycssloaders[withSpinner],
  logger[log_info, log_debug, log_warn, log_error],
  htmltools[HTML]
)

exploratory_module_ui <- function(id) {
  ns <- NS(id)
  # Consistent font styling (Nunito or fallback) applied within module scope
  font_css <- paste0(
    "#", ns("main_nav"), " { font-family: 'Nunito','Helvetica Neue','Arial','sans-serif'; } ",
    "#", ns("main_nav"), " .card-title, #", ns("main_nav"), " .nav-link, #", ns("main_nav"), " .plotly, #", ns("main_nav"), " .plot-container { font-family: inherit; }"
  )
  layout_sidebar(
    sidebar = sidebar(
      position = "right",
      card(
        card_title("Data Upload"),
        card_body(
          fileInput(ns("metadata"), "Metadata CSV", accept = c(".csv")),
          fileInput(ns("matrix"), "Expression Matrix CSV", accept = c(".csv")),
          div(style="font-size:0.75em;", p("Metadata: must contain 'sampleName' + ≥1 covariate column."))
        )
      ),
      card(
        card_title("Options"),
        card_body(
          uiOutput(ns("sidebar_options"))
        )
      )
    ),
    navset_card_tab(id = ns("main_nav"),
                    nav_panel("Boxplot", withSpinner(plotOutput(ns("boxplot"), height = "420px"), type = 4)),
                    nav_panel("Heatmap", withSpinner(plotlyOutput(ns("heatmap"), height = "420px"), type = 4)),
                    nav_panel("PCA", withSpinner(plotOutput(ns("pca_plot"), height = "420px"), type = 4)),
                    footer = div(style = "margin-top:0.5rem;", uiOutput(ns("status")))
    ),
    # Inject font CSS
    tags$style(font_css)
  )
}

exploratory_module_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive: read metadata
    metadata <- reactive({
      req(input$metadata)
      df <- try(read.csv(input$metadata$datapath, stringsAsFactors = FALSE), silent = TRUE)
      validate(need(!inherits(df, "try-error"), "Failed to read metadata CSV."))
      validate(need("sampleName" %in% names(df), "Metadata must have column: sampleName"))
      # Require at least one other column to act as covariate
      validate(need(length(setdiff(names(df), "sampleName")) >= 1,
                    "Metadata needs at least one additional covariate column besides sampleName."))
      # Deduplicate sampleName if repeated
      if (any(duplicated(df$sampleName))) {
        dupes <- unique(df$sampleName[duplicated(df$sampleName)])
        warning(sprintf("Duplicate sampleNames removed: %s", paste(dupes, collapse=", ")))
        df <- df[!duplicated(df$sampleName), ]
      }
      log_info(sprintf("Metadata loaded: %d samples; covariate columns: %s", nrow(df), paste(setdiff(names(df), 'sampleName'), collapse=',')))
      df
    })
    
    # Reactive: read matrix
    expr_matrix <- reactive({
      req(input$matrix)
      mat <- try(read.csv(input$matrix$datapath, row.names = 1, check.names = FALSE), silent = TRUE)
      validate(need(!inherits(mat, "try-error"), "Failed to read matrix CSV."))
      validate(need(nrow(mat) > 1, "Matrix must have at least 2 targets."))
      # Coerce numeric
      num_mat <- suppressWarnings(apply(mat, 2, as.numeric))
      # If coercion fails (all NA) for a column, keep original as character -> validation fails
      na_cols <- colSums(is.na(num_mat)) == nrow(num_mat)
      validate(need(!any(na_cols), paste("Non-numeric columns found:", paste(names(mat)[na_cols], collapse=", "))))
      rownames(num_mat) <- rownames(mat)
      log_info(sprintf("Expression matrix loaded: %d targets x %d samples", nrow(num_mat), ncol(num_mat)))
      as.data.frame(num_mat, check.names = FALSE)
    })
    
    # Validation: sample name alignment
    output$status <- renderUI({
      if (is.null(input$metadata) || is.null(input$matrix)) return(p("Awaiting uploads..."))
      mdf <- metadata(); mat <- expr_matrix()
      meta_samples <- mdf$sampleName
      matrix_samples <- colnames(mat)
      missing_in_matrix <- setdiff(meta_samples, matrix_samples)
      missing_in_meta <- setdiff(matrix_samples, meta_samples)
      msgs <- list()
      if (length(missing_in_matrix)) msgs <- c(msgs, sprintf("Samples missing in matrix: %s", paste(missing_in_matrix, collapse=", ")))
      if (length(missing_in_meta)) msgs <- c(msgs, sprintf("Samples missing in metadata: %s", paste(missing_in_meta, collapse=", ")))
      if (!length(msgs)) {
        popover_id <- ns("alignment_info")
        detail_html <- HTML(paste0(
          "<strong>All samples aligned</strong>:<br>",
          "Every 'sampleName' in the metadata file matches a column in the expression matrix, ",
          "with no missing or extra samples. This guarantees one-to-one pairing for boxplots, heatmap, and PCA.")
        )
        msgs <- list(
          tags$span(
            "All samples aligned ",
            tags$a(
              href = "#", id = popover_id, `data-bs-toggle` = "popover", `data-bs-trigger` = "focus", tabindex = "0",
              `data-bs-html` = "true", `data-bs-content` = detail_html, `data-bs-placement` = "right",
              style = "text-decoration:none; cursor:help;",
              tags$span(class = "badge bg-info", "?")
            )
          )
        )
      }
      scripts <- tags$script(HTML("document.querySelectorAll('[data-bs-toggle=popover]').forEach(function(el){if(window.bootstrap){window.bootstrap.Popover.getOrCreateInstance(el);} });"))
      div(lapply(msgs, p), scripts)
    })
    
    # Target selector UI
    output$target_selector <- renderUI({
      req(expr_matrix())
      selectInput(ns("target"), "Target for boxplot", choices = rownames(expr_matrix()), selected = rownames(expr_matrix())[1])
    })
    
    # Covariate selector UI (choose any non-sampleName column)
    output$covariate_selector <- renderUI({
      req(metadata())
      mdf <- metadata()
      covariate_cols <- setdiff(names(mdf), "sampleName")
      if (!length(covariate_cols)) return(p("No covariate columns available."))
      selectInput(ns("covariate_col"), "Covariate column", choices = covariate_cols, selected = covariate_cols[1])
    })
    
    # Boxplot
    output$boxplot <- renderPlot({
      req(metadata(), expr_matrix(), input$target, input$covariate_col)
      mdf <- metadata(); mat <- expr_matrix()
      validate(need(input$target %in% rownames(mat), "Selected target not found."))
      validate(need(input$covariate_col %in% names(mdf), "Covariate column not found."))
      # Build data frame for plotting
      df <- data.frame(
        expression = as.numeric(mat[input$target, mdf$sampleName]),
        covariate = mdf[[input$covariate_col]]
      )
      ggplot(df, aes(x = covariate, y = expression, fill = covariate)) +
        geom_boxplot(outlier.alpha = 0.6) +
        theme_minimal() +
        labs(title = paste("Boxplot:", input$target), x = input$covariate_col, y = "Expression")
    })
    
    # Interactive heatmap with hover using pure plotly
    output$heatmap <- renderPlotly({
      req(metadata(), expr_matrix())
      mdf <- metadata(); mat <- expr_matrix()
      shared_samples <- intersect(mdf$sampleName, colnames(mat))
      validate(need(length(shared_samples) > 1, "Need at least 2 matching samples for heatmap."))
      sub_mat <- mat[, shared_samples, drop = FALSE]
      vars <- apply(sub_mat, 1, var, na.rm = TRUE)
      vars[is.na(vars)] <- 0
      ord <- order(vars, decreasing = TRUE)
      top_n_val <- if (!is.null(input$top_n)) input$top_n else 25
      top_n <- min(top_n_val, length(ord))
      top_idx <- ord[seq_len(top_n)]
      heat_mat <- sub_mat[top_idx, , drop = FALSE]
      # Ensure unique target names to avoid hover text expansion issues
      if (any(duplicated(rownames(heat_mat)))) {
        dup_ct <- sum(duplicated(rownames(heat_mat)))
        rownames(heat_mat) <- make.unique(rownames(heat_mat), sep = "_dup")
        log_warn(sprintf("Made %d duplicated target names unique for heatmap", dup_ct))
      }
      scaled <- t(scale(t(as.matrix(heat_mat))))
      log_debug(sprintf("Heatmap render: top_n=%d; scaled range=%.3f..%.3f", top_n, min(scaled, na.rm=TRUE), max(scaled, na.rm=TRUE)))
      # Build hover text matrix with exact same dimensions as 'scaled'
      hover_text <- matrix(
        paste0(
          "Target: ", rep(rownames(scaled), times = ncol(scaled)),
          "<br>Sample: ", rep(colnames(scaled), each = nrow(scaled)),
          "<br>Z-score: ", round(c(scaled), 3)
        ),
        nrow = nrow(scaled),
        ncol = ncol(scaled),
        byrow = FALSE,
        dimnames = list(rownames(scaled), colnames(scaled))
      )
      plot_ly(
        x = colnames(scaled),
        y = rownames(scaled),
        z = scaled,
        type = "heatmap",
        text = hover_text,
        hoverinfo = "text",
        colorscale = list(
          list(0, '#313695'), list(0.1,'#4575B4'), list(0.2,'#74ADD1'), list(0.3,'#ABD9E9'),
          list(0.4,'#E0F3F8'), list(0.5,'#FFFFBF'), list(0.6,'#FEE090'), list(0.7,'#FDAE61'),
          list(0.8,'#F46D43'), list(0.9,'#D73027'), list(1,'#A50026')
        )
      ) |> plotly::layout(
        title = sprintf("Heatmap: Top %d variable targets", top_n),
        xaxis = list(title = "Sample", tickangle = -45),
        yaxis = list(title = "Target", autorange = "reversed"),
        margin = list(l=80,r=20,b=80,t=50)
      )
    })
    
    # Dynamic sidebar options based on selected tab
    output$sidebar_options <- renderUI({
      current <- input$main_nav
      if (is.null(current)) return(NULL)
      switch(current,
             "Boxplot" = tagList(
               uiOutput(ns("target_selector")),
               uiOutput(ns("covariate_selector"))
             ),
             "Heatmap" = tagList(
               numericInput(ns("top_n"), "Top variable targets (heatmap)", value = 25, min = 5, max = 200, step = 5)
             ),
             "PCA" = tagList(
               uiOutput(ns("covariate_selector")),
               uiOutput(ns("pca_axes"))
             ),
             NULL)
    })
    
    # Reactive PCA data frame (computed once for plotting + axis selectors)
    pca_data <- reactive({
      req(metadata(), expr_matrix(), input$covariate_col)
      mdf <- metadata(); mat <- expr_matrix()
      shared_samples <- intersect(mdf$sampleName, colnames(mat))
      validate(need(length(shared_samples) > 2, "Need >=3 matching samples for PCA."))
      sub_mat <- mat[, shared_samples, drop = FALSE]
      vars <- apply(sub_mat, 1, var, na.rm = TRUE)
      keep <- vars > 0
      validate(need(sum(keep) > 2, "Insufficient variable targets for PCA."))
      pca <- prcomp(t(sub_mat[keep, ]), scale. = TRUE)
      pcs <- as.data.frame(pca$x)
      pcs$sampleName <- rownames(pcs)
      validate(need(input$covariate_col %in% names(mdf), "Covariate column not found."))
      pcs <- merge(pcs, mdf[, c("sampleName", input$covariate_col)], by = "sampleName", all.x = TRUE)
      colnames(pcs)[ncol(pcs)] <- "covariate"
      log_info(sprintf("PCA computed: %d samples, %d components", nrow(pcs), sum(grepl('^PC[0-9]+$', names(pcs)))))
      pcs
    })
    
    # PCA axis selectors (appear only when PCA tab active)
    output$pca_axes <- renderUI({
      req(pca_data())
      pc_cols <- grep("^PC[0-9]+$", names(pca_data()), value = TRUE)
      if (length(pc_cols) < 2) return(p("Not enough principal components."))
      # Defaults
      default_x <- if (!is.null(input$pca_x) && input$pca_x %in% pc_cols) input$pca_x else pc_cols[1]
      default_y <- if (!is.null(input$pca_y) && input$pca_y %in% pc_cols) input$pca_y else pc_cols[2]
      tagList(
        selectInput(ns("pca_x"), "X Axis PC", choices = pc_cols, selected = default_x),
        selectInput(ns("pca_y"), "Y Axis PC", choices = pc_cols, selected = default_y)
      )
    })
    
    # PCA plot
    output$pca_plot <- renderPlot({
      pcs <- pca_data()
      req(pcs)
      pc_cols <- grep("^PC[0-9]+$", names(pcs), value = TRUE)
      validate(need(length(pc_cols) >= 2, "Insufficient principal components."))
      x_pc <- if (!is.null(input$pca_x) && input$pca_x %in% pc_cols) input$pca_x else pc_cols[1]
      y_pc <- if (!is.null(input$pca_y) && input$pca_y %in% pc_cols) input$pca_y else pc_cols[min(2, length(pc_cols))]
      validate(need(x_pc != y_pc, "Choose two different PCs for X and Y."))
      log_debug(sprintf("PCA plot render: x=%s y=%s", x_pc, y_pc))
      ggplot(pcs, aes(x = !!sym(x_pc), y = !!sym(y_pc), color = covariate)) +
        geom_point(size = 3, alpha = 0.85) +
        theme_minimal() +
        labs(title = sprintf("PCA: %s vs %s", x_pc, y_pc), color = input$covariate_col, x = x_pc, y = y_pc)
    })
  })
}

# Standalone mode for development
if (sys.nframe() == 0) {
  shinyApp(exploratory_module_ui("explore1"), function(input, output, session) {
    exploratory_module_server("explore1")
  })
}

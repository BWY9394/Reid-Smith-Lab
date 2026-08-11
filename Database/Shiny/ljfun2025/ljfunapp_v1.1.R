library(shiny)
library(DBI)
library(RSQLite)
library(dplyr)
library(ggplot2)
library(DT)
library(sortable)

ui <- fluidPage(
  titlePanel("ljfun RNA-seq browser"),
  
  sidebarLayout(
    sidebarPanel(
      #textInput("gene_search", "Search locusName", value = ""),
      
      textAreaInput(
        "paste_loci",
        "Paste locusName IDs",
        value = "",
        rows = 5,
        placeholder = "One locusName per line, or comma/space separated"
      ),
      
      actionButton("select_pasted", "Select pasted loci"),
      actionButton("clear_selection", "Clear all selections"),
      
      actionButton("reset_order", "Reset x-axis order"),
      
      hr(),
      h4("Drag to reorder x-axis"),
      uiOutput("x_order_ui"),
      
      hr(),
      h4("Columns to display"),
      uiOutput("column_selector"),
      
      width = 3
    ),
    
    mainPanel(
      h3("Expression plot"),
      checkboxInput("custom_height", "Use custom plot height (px)", value = FALSE),
      conditionalPanel(
        condition = "input.custom_height == true",
        numericInput("plot_height", "Plot height (px)", value = 1000, min = 200, max = 4000, step = 50)
      ),
      uiOutput("plot_ui"),
      uiOutput("plot_info"),
      hr(),
      h3("Selected gene boxplots"),
      uiOutput("gene_boxplot_selector"),
      uiOutput("gene_boxplot_ui"),
      selectInput("download_format", "Download format:",
                  choices = c("PNG", "PDF", "SVG"), selected = "PNG"),
      h5("Download dimensions (inches)"),
      fluidRow(
        column(6,
               numericInput("download_width", "Width", value = 12, min = 4, max = 20, step = 0.5)
        ),
        column(6,
               numericInput("download_height", "Height", value = 8, min = 4, max = 20, step = 0.5)
        )
      ),
      downloadButton("download_plot", "Download plot"),
      
      hr(),
      
      h3("Gene metadata"),
      DTOutput("gene_table")
    )
  )
)

server <- function(input, output, session) {
  
  con <- dbConnect(SQLite(), "rnaseq_lj_test_copy.sqlite")
  onStop(function() dbDisconnect(con))
  
  metadata_tbl <- dbReadTable(con, "lotusmetadata")
  counts_tbl   <- dbReadTable(con, "ljnormcounts")
  atlast_tbl  <- dbReadTable(con, "gifuatlas")
  
  default_cols <- c(
    "locusName",
    "Gene Symbol",
    "Function",
    "Best.hit.arabi.name",
    #"",
    "Best.hit.arabi.defline"
  )
  default_cols <- default_cols[default_cols %in% names(metadata_tbl)]
  
  all_combids <- sort(unique(counts_tbl$CombID))
  
  output$column_selector <- renderUI({
    checkboxGroupInput(
      "display_cols",
      label = NULL,
      choices = names(metadata_tbl),
      selected = default_cols
    )
  })
  
  pasted_loci <- reactive({
    pasted <- unlist(strsplit(input$paste_loci, "[,;[:space:]]+"))
    pasted <- trimws(pasted)
    pasted[pasted != ""]
  })
  
  filtered_metadata <- reactive({
    dat <- metadata_tbl
    
    #if (input$gene_search != "") {
    #  dat <- dat %>%
    #    filter(grepl(input$gene_search, locusName, ignore.case = TRUE))
    #  }
    
    if (length(pasted_loci()) > 0) {
      dat <- dat %>%
        filter(locusName %in% pasted_loci())
    }
    
    dat
  })
  
  output$gene_table <- renderDT({
    req(input$display_cols)
    
    display_tbl <- filtered_metadata() %>%
      select(all_of(input$display_cols))
    
    datatable(
      display_tbl,
      selection = list(mode = "multiple", target = "row"),
      filter = "top",
      options = list(
        pageLength = 10,
        scrollX = TRUE
      ),
      rownames = FALSE
    )
  })
  
  selected_loci <- reactiveVal(character(0))
  
  observeEvent(input$gene_table_rows_selected, {
    rows <- input$gene_table_rows_selected
    
    if (!is.null(rows) && length(rows) > 0) {
      selected_loci(filtered_metadata()$locusName[rows])
    }
  })
  
  observeEvent(input$clear_selection, {
    selected_loci(character(0))
    
    proxy <- dataTableProxy("gene_table")
    selectRows(proxy, NULL)
  })
  
  observeEvent(input$select_pasted, {
    valid <- pasted_loci()[pasted_loci() %in% metadata_tbl$locusName]
    selected_loci(unique(valid))
  })
  
  selected_genes <- reactive({
    selected_loci()
  })
  
  boxplot_genes <- reactive({
    sel <- selected_genes()
    if (length(sel) == 0) {
      return(character(0))
    }
    
    if (!is.null(input$boxplot_gene_selection) && length(input$boxplot_gene_selection) > 0) {
      intersect(sel, input$boxplot_gene_selection)
    } else {
      character(0)
    }
  })
  
  ordered_selected_genes <- reactive({
    sel <- selected_genes()
    pasted <- unique(pasted_loci())
    if (length(sel) == 0 || length(pasted) == 0) {
      return(sel)
    }
    ordered <- intersect(pasted, sel)
    extra <- setdiff(sel, ordered)
    c(ordered, extra)
  })
  
  current_order <- reactiveVal(all_combids)
  
  observeEvent(input$reset_order, {
    current_order(all_combids)
  })
  
  output$x_order_ui <- renderUI({
    rank_list(
      text = "Drag sample groups into the order you want",
      labels = current_order(),
      input_id = "combid_rank",
      options = sortable_options()
    )
  })
  
  observe({
    req(input$combid_rank)
    current_order(input$combid_rank)
  })
  
  # dynamic UI for plot height (toggleable)
  output$plot_ui <- renderUI({
    h <- if (isTRUE(input$custom_height)) {
      # guard: ensure numeric value
      ph <- as.integer(input$plot_height)
      if (is.na(ph) || ph < 100) ph <- 1000
      ph
    } else {
      500
    }
    plotOutput("counts_plot", height = paste0(h, "px"))
  })
  
  output$plot_info <- renderUI({
    req(selected_genes())
    info <- metadata_tbl %>%
      filter(locusName %in% selected_genes()) %>%
      arrange(match(locusName, ordered_selected_genes())) %>%
      select(locusName, Best.hit.arabi.defline)
    if (nrow(info) == 0) return(NULL)
    
    tagList(
      h4("Best.hit.arabi.defline"),
      lapply(seq_len(nrow(info)), function(i) {
        div(
          strong(info$locusName[i]), ": ",
          span(info$Best.hit.arabi.defline[i])
        )
      })
    )
  })
  
  output$gene_boxplot_selector <- renderUI({
    sel <- selected_genes()
    if (length(sel) == 0) {
      return(NULL)
    }
    
    checkboxGroupInput(
      inputId = "boxplot_gene_selection",
      label = "Choose genes to plot",
      choices = sel,
      selected = sel
    )
  })
  
  output$gene_boxplot_ui <- renderUI({
    plotOutput("gene_boxplot", height = "500px", width = "1500px")
  })
  
  output$gene_boxplot <- renderPlot({
    plot_genes <- boxplot_genes()
    
    if (length(plot_genes) == 0) {
      plot.new()
      text(0.5, 0.5, "Select one or more genes from the metadata table", cex = 1.2)
      return()
    }
    
    filtergenes <- atlast_tbl %>%
      filter(locusName %in% plot_genes)
    
    if (nrow(filtergenes) == 0) {
      plot.new()
      text(0.5, 0.5, "No data found for the plotted genes", cex = 1.2)
      return()
    }
    
    filtergenes <- filtergenes %>%
        mutate(
          Condition = factor(Condition, levels = unique(Condition)),
          locusName = factor(locusName, levels = plot_genes)
        )
    
    ggplot(filtergenes, aes(x = count, y = Condition, color = locusName)) +
      geom_boxplot() +
      theme_bw() +
      theme(axis.text.y = element_text(size = 10),
            legend.position = "top",
            legend.justification = "left") +
      labs(title = "Gene atlas", x = "Count", y = "Sample")
  })
  
  output$counts_plot <- renderPlot({
    
    if (length(selected_genes()) == 0) {
      plot.new()
      text(0.5, 0.5, "Select one or more genes from the metadata table", cex = 1.2)
      return()
    }
    
    plot_data <- counts_tbl %>%
      filter(locusName %in% selected_genes())
    
    if (nrow(plot_data) == 0) {
      plot.new()
      text(0.5, 0.5, "No data found for selected genes", cex = 1.2)
      return()
    }
    
    plot_data <- plot_data %>%
      mutate(
        Rep = factor(Rep),
        CombID = factor(CombID, levels = current_order()),
        locusName = factor(locusName, levels = ordered_selected_genes())
      )
    
    # build ggplot object so it can be reused for download
    plot_obj <- ggplot(plot_data, aes(x = CombID, y = normcounts)) +
      geom_boxplot(fatten=NULL,fill = NA, color = "gray70", outlier.shape = NA) +
      geom_point(aes(color = Rep), size = 2) +
      facet_wrap(~ locusName, nrow = 4,ncol=5, scales = "free") +
      labs(
        title = "Normalized counts of selected genes",
        y = "Normalized counts",
        x = "Genotypic response to N treatments at D1 and D3"
      ) +
      theme_bw() +
      theme(
        axis.title.x = element_text(color = "black", size = 15, face = "bold"),
        axis.title.y = element_text(color = "black", size = 15, face = "bold"),
        axis.text.x = element_text(size = 12, angle = -40, hjust = 0, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        strip.text.x = element_text(size = 13, face = "bold"),
        strip.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()
      ) +
      stat_summary(
        fun = mean,
        geom = "errorbar",
        aes(ymax = after_stat(y), ymin = after_stat(y)),
        color = "gray70",
        width = 0.7,
        size = 0.4,
        linetype = "solid"
      )
    
    
    print(plot_obj)
  })
  
  # allow downloading the last plot in the selected format
  output$download_plot <- downloadHandler(
    filename = function() {
      fmt <- tolower(input$download_format %||% "png")
      paste0("counts_plot_", Sys.Date(), ".", fmt)
    },
    content = function(file) {
      if (length(selected_genes()) == 0) {
        stop("No genes selected")
      }
      
      plot_data <- counts_tbl %>%
        filter(locusName %in% selected_genes()) %>%
        mutate(
          Rep = factor(Rep),
          CombID = factor(CombID, levels = current_order()),
          locusName = factor(locusName, levels = ordered_selected_genes())
        )
      
      p <- ggplot(plot_data, aes(x = CombID, y = normcounts)) +
        geom_boxplot(fatten=NULL,fill = NA, color = "gray70", outlier.shape = NA) +
        geom_point(aes(color = Rep), size = 2) +
        facet_wrap(~ locusName, nrow = 4,ncol=5, scales = "free") +
        labs(
          title = "Normalized counts of selected genes",
          y = "Normalized counts",
          x = "Genotypic response to N treatments at D1 and D3"
        ) +
        theme_bw() +
        theme(
          axis.title.x = element_text(color = "black", size = 15, face = "bold"),
          axis.title.y = element_text(color = "black", size = 15, face = "bold"),
          axis.text.x = element_text(size = 12, angle = -40, hjust = 0, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold"),
          strip.text.x = element_text(size = 13, face = "bold"),
          strip.background = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank()
        )+
        stat_summary(
          fun = mean,
          geom = "errorbar",
          aes(ymax = after_stat(y), ymin = after_stat(y)),
          color = "gray70",
          width = 0.7,
          size = 0.4,
          linetype = "solid"
        )
      
      fmt <- tolower(input$download_format %||% "png")
      dl_width <- input$download_width %||% 12
      dl_height <- input$download_height %||% 8
      
      # choose device/extension for ggsave
      if (fmt == "png") {
        ggplot2::ggsave(file, plot = p, device = "png", width = dl_width, height = dl_height, dpi = 300)
      } else if (fmt == "pdf") {
        ggplot2::ggsave(file, plot = p, device = "pdf", width = dl_width, height = dl_height)
      } else if (fmt == "svg") {
        ggplot2::ggsave(file, plot = p, device = "svg", width = dl_width, height = dl_height)
      } else {
        ggplot2::ggsave(file, plot = p, device = "png", width = dl_width, height = dl_height, dpi = 300)
      }
    }
  )
}

shinyApp(ui, server)
library(shiny)
library(DBI)
library(RSQLite)
library(dplyr)
library(ggplot2)
library(DT)
library(sortable)

ui <- fluidPage(
  titlePanel("VuFun RNA-seq browser"),
  
  sidebarLayout(
    sidebarPanel(
      textInput("gene_search", "Search locusName", value = ""),
      
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
      plotOutput("counts_plot", height = "500px"),
      
      hr(),
      
      h3("Gene metadata"),
      DTOutput("gene_table")
    )
  )
)

server <- function(input, output, session) {
  
  con <- dbConnect(SQLite(), "rnaseq_vufun_test.sqlite")
  onStop(function() dbDisconnect(con))
  
  metadata_tbl <- dbReadTable(con, "cowpeametadata")
  counts_tbl   <- dbReadTable(con, "vunormcounts")
  
  default_cols <- c(
    "locusName",
    "Cluster",
    "Gene Symbol...10",
    "Function...11",
    "RNAseqDEG",
    "plusNresponse"
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
  
  filtered_metadata <- reactive({
    if (input$gene_search == "") {
      metadata_tbl
    } else {
      metadata_tbl %>%
        filter(grepl(input$gene_search, locusName, ignore.case = TRUE))
    }
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
  
  selected_genes <- reactive({
    rows <- input$gene_table_rows_selected
    
    if (is.null(rows) || length(rows) == 0) {
      return(character(0))
    }
    
    filtered_metadata()$locusName[rows]
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
        CombID = factor(CombID, levels = current_order())
      )
    
    ggplot(plot_data, aes(x = CombID, y = normcounts)) +
      geom_boxplot(fill = NA, color = "gray70", outlier.shape = NA) +
      geom_point(aes(color = Rep), size = 2) +
      facet_wrap(~ locusName, nrow = 2, scales = "free") +
      labs(
        title = "Normalized counts of selected genes",
        y = "Normalized counts",
        x = "ZBF line"
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
      )
  })
}

shinyApp(ui, server)

library(shiny)
library(DBI)
library(RSQLite)
library(dplyr)
library(ggplot2)
library(DT)

ui <- fluidPage(
  titlePanel("RNA-seq Gene Browser"),
  sidebarLayout(
    sidebarPanel(
      textInput("gene_search", "Gene ID or symbol"),
      selectInput("cluster", "Cluster", choices = c("All")),
      actionButton("run_filter", "Search")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Genes", DTOutput("gene_table")),
        tabPanel("Gene details", tableOutput("gene_details")),
        tabPanel("Expression plot", plotOutput("counts_plot"))
      )
    )
  )
)

server <- function(input, output, session) {
  con <- dbConnect(SQLite(), "rnaseq_lab.sqlite")

  genes_data <- eventReactive(input$run_filter, {
    dbReadTable(con, "genes") %>%
      filter(
        grepl(input$gene_search, locusName, ignore.case = TRUE) |
        grepl(input$gene_search, gene_symbol, ignore.case = TRUE)
      )
  })

  output$gene_table <- renderDT({
    genes_data()
  }, selection = "multiple")

  selected_genes <- reactive({
    req(input$gene_table_rows_selected)
    genes_data()$locusName[input$gene_table_rows_selected]
  })

  output$gene_details <- renderTable({
    req(selected_genes())
    dbReadTable(con, "genes") %>%
      filter(locusName %in% selected_genes())
  })

  output$counts_plot <- renderPlot({
    req(selected_genes())

    counts <- dbReadTable(con, "normalized_counts")
    samples <- dbReadTable(con, "samples")

    plot_data <- counts %>%
      filter(locusName %in% selected_genes()) %>%
      left_join(samples, by = "sample_id")

    ggplot(plot_data, aes(x = CombID, y = normcounts)) +
      geom_boxplot(fill = NA, color = "gray90", outlier.shape = NA) +
      geom_point(aes(color = factor(Rep))) +
      facet_wrap(~locusName, scales = "free") +
      theme_bw()
  })
}

shinyApp(ui, server)
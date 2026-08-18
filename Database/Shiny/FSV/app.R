#Author: Benjamin Wee Yang
#Email: B.weeyang@latrobe.edu.au

# Agricultural Trial Spatial Variation Visualizer

library(shiny)
library(bslib)
library(ggplot2)
library(plotly)
library(DT)
library(dplyr)

options(shiny.maxRequestSize = 30 * 1024^2)

# App metadata and required column definitions.
APP_TITLE <- "FieldSpatialViewer"
REQUIRED_COLUMNS <- c("row", "column")
NON_TRAIT_COLUMNS <- c(
  "trial_id", "plot_id", "row", "column", "replicate", "block",
  "treatment", "genotype", "x", "y", "latitude", "longitude"
)

# Convert values to numeric while suppressing warnings from coercion.
# Non-numeric values become NA, so they are treated as missing rather than as a hard rejection.
safe_numeric <- function(x) {
  y <- suppressWarnings(as.numeric(x))
  y[is.na(y)] <- NA_real_
  y
}

# Read a trial CSV and validate that row/column coordinates exist and are numeric.
read_trial_csv <- function(path) {
  dat <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                  na.strings = c("", " ", "NA", "N/A", "n/a", "NULL", "null"))
  names(dat) <- make.names(tolower(trimws(names(dat))), unique = TRUE)

  missing <- setdiff(REQUIRED_COLUMNS, names(dat))
  if (length(missing) > 0) {
    stop("CSV is missing required column(s): ", paste(missing, collapse = ", "))
  }

  dat$row <- safe_numeric(dat$row)
  dat$column <- safe_numeric(dat$column)
  if (any(!is.finite(dat$row)) || any(!is.finite(dat$column))) {
    stop("Columns 'row' and 'column' must contain numeric values only.")
  }
  if (anyDuplicated(dat[c("row", "column")])) {
    stop("Each row-column coordinate must identify one unique plot.")
  }

  if (!"plot_id" %in% names(dat)) {
    dat$plot_id <- sprintf("R%s-C%s", dat$row, dat$column)
  }
  dat
}

# Identify numeric columns that are likely trait measurements for visualization.
# Ignore the previous 50% completeness requirement and only require a reasonable
# number of valid numeric observations.
numeric_traits <- function(dat) {
  candidates <- setdiff(names(dat), NON_TRAIT_COLUMNS)
  candidates[vapply(dat[candidates], function(x) {
    y <- safe_numeric(x)
    valid <- y[is.finite(y)]
    length(valid) >= 4 && length(unique(valid)) > 2
  }, logical(1))]
}

# Fit the spatial trend and residuals for the selected trait.
# Optional treatment and replicate adjustment is performed first.
fit_spatial_components <- function(dat, response, treatment_col = NULL, replicate_col = NULL) {
  y <- safe_numeric(dat[[response]])
  keep <- is.finite(y) & is.finite(dat$row) & is.finite(dat$column)
  adjusted <- rep(NA_real_, nrow(dat))
  trend <- rep(NA_real_, nrow(dat))
  residual <- rep(NA_real_, nrow(dat))

  model_data <- dat[keep, , drop = FALSE]
  model_data$.response <- y[keep]

  fixed_terms <- character(0)
  if (!is.null(treatment_col) && treatment_col %in% names(model_data)) {
    model_data[[treatment_col]] <- factor(model_data[[treatment_col]])
    if (nlevels(model_data[[treatment_col]]) > 1) fixed_terms <- c(fixed_terms, treatment_col)
  }
  if (!is.null(replicate_col) && replicate_col %in% names(model_data)) {
    model_data[[replicate_col]] <- factor(model_data[[replicate_col]])
    if (nlevels(model_data[[replicate_col]]) > 1) fixed_terms <- c(fixed_terms, replicate_col)
  }

  if (length(fixed_terms) > 0) {
    fixed_formula <- as.formula(paste(".response ~", paste(fixed_terms, collapse = " + ")))
    fixed_fit <- lm(fixed_formula, data = model_data)
    adjusted_values <- residuals(fixed_fit) + mean(model_data$.response, na.rm = TRUE)
  } else {
    adjusted_values <- model_data$.response
  }
  adjusted[keep] <- adjusted_values

  spatial_data <- data.frame(
    adjusted = adjusted_values,
    row = model_data$row,
    column = model_data$column
  )

  # A robust, dependency-light 2D smoother. Span increases slightly for small trials.
  span <- if (nrow(spatial_data) < 120) 0.65 else 0.45
  smooth_fit <- try(
    loess(adjusted ~ column + row, data = spatial_data, span = span, degree = 2,
          surface = "direct", control = loess.control(surface = "direct")),
    silent = TRUE
  )

  if (inherits(smooth_fit, "try-error")) {
    fallback <- lm(adjusted ~ poly(column, 2, raw = TRUE) + poly(row, 2, raw = TRUE) + column:row,
                   data = spatial_data)
    trend_values <- predict(fallback, newdata = spatial_data)
  } else {
    trend_values <- predict(smooth_fit, newdata = spatial_data)
  }

  trend[keep] <- trend_values
  residual[keep] <- adjusted_values - trend_values

  list(raw = y, adjusted = adjusted, trend = trend, residual = residual)
}

# Format numeric values for display in tooltips and labels.
format_value <- function(x) {
  ifelse(is.finite(x), scales::number(x, accuracy = 0.01), "NA")
}

# Format axis labels for row/column coordinates: keep integers as whole numbers,
# but preserve decimal precision when the input coordinate is non-integer.
format_axis_label <- function(x) {
  ifelse(abs(x - round(x)) < 1e-6, as.integer(round(x)), scales::number(x, accuracy = 0.01))
}

# Build the app user interface: sidebar controls, theme, and main tab panels.
ui <- page_sidebar(
  title = div(
    class = "brand-wrap",
    div(class = "brand-mark", "FSV"),
    div(
      div(class = "brand-title", APP_TITLE),
      div(class = "brand-subtitle", "Agricultural Field trial Spatial Variation viewer")
    )
  ),
  theme = bs_theme(
    version = 5,
    bootswatch = "minty",
    bg = "#f5f7f6",
    fg = "#17211b",
    primary = "#176b4d",
    secondary = "#5e7468",
    success = "#23855d",
    base_font = font_google("Inter"),
    heading_font = font_google("Manrope")
  ),
  fillable = TRUE,
  sidebar = sidebar(
    width = 330,
    open = "desktop",
    div(class = "sidebar-intro",
        h5("Explore your field"),
        p("Upload plot-level data, isolate environmental variation, and inspect spatial structure.")),
    fileInput(
      "file", "Upload trial CSV", accept = c(".csv"),
      buttonLabel = "Choose CSV", placeholder = "Using built-in demo data"
    ),
    uiOutput("trait_ui"),
    selectInput(
      "view_mode", "Map layer",
      choices = c(
        "Observed values" = "raw",
        "Treatment-adjusted values" = "adjusted",
        "Estimated environmental trend" = "trend",
        "Residual after spatial trend" = "residual"
      ),
      selected = "trend"
    ),
    uiOutput("treatment_ui"),
    uiOutput("replicate_ui"),
    selectInput(
      "palette", "Colour palette",
      choices = c("Viridis" = "viridis", "Magma" = "magma", "Plasma" = "plasma", "Cividis" = "cividis"),
      selected = "viridis"
    ),
    checkboxInput("reverse_rows", "Display field row 1 at top", TRUE),
    checkboxInput("show_labels", "Show plot labels", FALSE),
    hr(),
    downloadButton("download_processed", "Download processed data", class = "btn-primary w-100"),
    div(class = "help-note",
        tags$strong("Required columns:"), " row, column.", tags$br(),
        "Add numeric trait columns plus optional treatment, replicate, block, genotype and plot_id fields.")
  ),
  tags$head(
    tags$style(HTML("
      :root { --fl-green:#176b4d; --fl-ink:#17211b; --fl-muted:#66756d; }
      body { letter-spacing:-0.01em; }
      .bslib-page-sidebar > .navbar { box-shadow:0 1px 0 rgba(23,33,27,.08); }
      .brand-wrap { display:flex; align-items:center; gap:.8rem; }
      .brand-mark { width:38px; height:38px; border-radius:12px; display:grid; place-items:center;
        color:white; background:linear-gradient(145deg,#176b4d,#2a9270); font-weight:800; font-size:.8rem;
        box-shadow:0 8px 18px rgba(23,107,77,.24); }
      .brand-title { font-weight:800; line-height:1.05; }
      .brand-subtitle { color:#6b7b72; font-size:.75rem; margin-top:.15rem; }
      .sidebar-intro { padding:.2rem 0 .65rem; }
      .sidebar-intro h5 { font-weight:750; margin-bottom:.25rem; }
      .sidebar-intro p, .help-note { color:var(--fl-muted); font-size:.84rem; line-height:1.45; }
      .help-note { margin-top:.8rem; padding:.8rem; border-radius:12px; background:rgba(23,107,77,.06); }
      .card { border:1px solid rgba(23,33,27,.08); border-radius:18px; box-shadow:0 10px 30px rgba(33,54,43,.06); overflow:hidden; }
      .card-header { background:transparent; border-bottom:1px solid rgba(23,33,27,.07); font-weight:750; }
      .value-box { border-radius:18px; box-shadow:0 10px 30px rgba(33,54,43,.06); }
      .value-box-title { font-size:.77rem; letter-spacing:.04em; text-transform:uppercase; opacity:.82; }
      .value-box-value { font-weight:800; }
      .nav-tabs { border:0; gap:.4rem; }
      .nav-tabs .nav-link { border:0; border-radius:10px; color:#526159; font-weight:650; }
      .nav-tabs .nav-link.active { color:var(--fl-green); background:rgba(23,107,77,.09); }
      .control-label, .form-label { font-size:.83rem; font-weight:700; color:#34473d; }
      .form-control, .form-select, .selectize-input { border-radius:11px !important; border-color:#d6dfda !important; }
      .plotly.html-widget { border-radius:14px; }
      .plot-detail { padding:.85rem 1rem; border-radius:14px; background:#f6f9f7; min-height:105px; }
      .plot-detail .plot-id { font-size:1.05rem; font-weight:800; color:var(--fl-green); }
      .plot-detail-grid { display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); gap:.3rem .8rem; margin-top:.55rem; font-size:.85rem; }
      .status-pill { display:inline-flex; align-items:center; gap:.35rem; padding:.3rem .55rem; border-radius:999px;
        background:#e8f5ee; color:#176b4d; font-size:.76rem; font-weight:750; }
      .status-dot { width:7px; height:7px; border-radius:50%; background:#22a06b; }
      @media (max-width: 768px) { .brand-subtitle { display:none; } }
    "))),
  layout_columns(
    col_widths = c(3, 3, 3, 3),
    value_box(title = "Plots", value = textOutput("kpi_plots"), showcase = bsicons::bs_icon("grid-3x3-gap"), theme = "primary"),
    value_box(title = "Trait mean", value = textOutput("kpi_mean"), showcase = bsicons::bs_icon("activity"), theme = "success"),
    value_box(title = "Field range", value = textOutput("kpi_range"), showcase = bsicons::bs_icon("arrows-expand"), theme = "info"),
    value_box(title = "Spatial signal", value = textOutput("kpi_spatial"), showcase = bsicons::bs_icon("broadcast-pin"), theme = "warning")
  ),
  navset_card_tab(
    id = "main_tabs",
    nav_panel(
      "Field map",
      layout_columns(
        col_widths = c(8, 4),
        card(
          full_screen = TRUE,
          card_header(div(class = "d-flex justify-content-between align-items-center",
                          span("Interactive plot map"),
                          uiOutput("data_status"))),
          plotlyOutput("field_map", height = "620px")
        ),
        card(
          card_header("Selected plot"),
          uiOutput("plot_details"),
          hr(),
          card_header("Layer distribution"),
          plotlyOutput("distribution", height = "340px")
        )
      )
    ),
    nav_panel(
      "Profiles",
      layout_columns(
        col_widths = c(6, 6),
        card(full_screen = TRUE, card_header("Mean by field row"), plotlyOutput("row_profile", height = "500px")),
        card(full_screen = TRUE, card_header("Mean by field column"), plotlyOutput("column_profile", height = "500px"))
      )
    ),
    nav_panel(
      "Diagnostics",
      layout_columns(
        col_widths = c(7, 5),
        card(full_screen = TRUE, card_header("Observed vs spatial trend"), plotlyOutput("observed_fitted", height = "500px")),
        card(card_header("Model summary"), uiOutput("model_summary"), plotOutput("residual_qq", height = "330px"))
      )
    ),
    nav_panel(
      "Data",
      card(full_screen = TRUE, card_header("Trial data and calculated layers"), DTOutput("data_table"))
    )
  )
)

server <- function(input, output, session) {
  # Path to the bundled demo data used when no file is uploaded.
  demo_path <- "mock_agricultural_trial.csv"

  # Reactive loaded trial data; uses the uploaded file if available.
  trial_data <- reactive({
    path <- if (is.null(input$file)) demo_path else input$file$datapath
    validate(need(file.exists(path), paste0("Demo file not found. Place mock_agricultural_trial.csv beside app.R or upload a CSV.")))
    tryCatch(
      read_trial_csv(path),
      error = function(e) stop(conditionMessage(e), call. = FALSE)
    )
  })

  # Populate the trait selector based on detected numeric trait columns.
  observeEvent(trial_data(), {
    traits <- numeric_traits(trial_data())
    selected <- if ("yield_t_ha" %in% traits) "yield_t_ha" else traits[[1]]
    updateSelectInput(session, "trait", choices = traits, selected = selected)
  }, ignoreInit = FALSE)

  # Render the trait selection UI dynamically.
  output$trait_ui <- renderUI({
    traits <- numeric_traits(trial_data())
    validate(need(length(traits) > 0, "No numeric trait columns were found."))
    selectInput("trait", "Trait", choices = traits, selected = if ("yield_t_ha" %in% traits) "yield_t_ha" else traits[[1]])
  })

  output$treatment_ui <- renderUI({
    dat <- trial_data()
    choices <- intersect(c("treatment", "genotype"), names(dat))
    selectInput("treatment_col", "Adjust for treatment factor", choices = c("None" = "", choices),
                selected = if ("treatment" %in% choices) "treatment" else "")
  })

  output$replicate_ui <- renderUI({
    dat <- trial_data()
    choices <- intersect(c("replicate", "block"), names(dat))
    selectInput("replicate_col", "Adjust for design factor", choices = c("None" = "", choices),
                selected = if ("replicate" %in% choices) "replicate" else "")
  })

  analysed <- reactive({
    req(input$trait)
    dat <- trial_data()
    validate(need(input$trait %in% names(dat), "Choose a valid trait."))
    treatment_col <- if (is.null(input$treatment_col) || input$treatment_col == "") NULL else input$treatment_col
    replicate_col <- if (is.null(input$replicate_col) || input$replicate_col == "") NULL else input$replicate_col
    components <- fit_spatial_components(dat, input$trait, treatment_col, replicate_col)
    dat$.raw <- components$raw
    dat$.adjusted <- components$adjusted
    dat$.trend <- components$trend
    dat$.residual <- components$residual
    dat
  })

  # Add the selected display layer to the dataset for plotting.
  mapped_data <- reactive({
    dat <- analysed()
    layer_col <- paste0(".", input$view_mode)
    dat$.display <- dat[[layer_col]]
    dat
  })

  selected_plot <- reactiveVal(NULL)

  # Track which plot was clicked in the interactive field map.
  observeEvent(event_data("plotly_click", source = "fieldmap"), {
    click <- event_data("plotly_click", source = "fieldmap")
    if (!is.null(click$key)) selected_plot(as.character(click$key))
  })

  # Render the main field grid map with the selected display layer.
  output$field_map <- renderPlotly({
    dat <- mapped_data()
    validate(need(any(is.finite(dat$.display)), "No finite values are available for this layer."))

    label_text <- paste0(
      "<b>", dat$plot_id, "</b>",
      "<br>Row ", dat$row, " · Column ", dat$column,
      "<br>", input$trait, ": ", format_value(dat$.raw),
      "<br>Displayed: ", format_value(dat$.display),
      if ("treatment" %in% names(dat)) paste0("<br>Treatment: ", dat$treatment) else "",
      if ("genotype" %in% names(dat)) paste0("<br>Genotype: ", dat$genotype) else ""
    )

    p <- ggplot(dat, aes(x = column, y = row, fill = .display, key = plot_id, text = label_text)) +
      geom_tile(colour = "white", linewidth = 0.35) +
      coord_equal(expand = FALSE) +
      scale_x_continuous(breaks = function(x) unique(floor(pretty(c(0, max(x)), n = 10))), labels = format_axis_label) +
      scale_y_continuous(breaks = function(x) unique(floor(pretty(c(0, max(x)), n = 10))), labels = format_axis_label) +
      scale_fill_viridis_c(option = input$palette, na.value = "#e6ebe8") +
      labs(x = "Field column", y = "Field row", fill = input$view_mode) +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid = element_blank(),
        axis.title = element_text(face = "bold", colour = "#46564d"),
        axis.text = element_text(colour = "#5f6f66"),
        legend.position = "bottom",
        legend.key.width = grid::unit(2.2, "cm"),
        plot.margin = margin(12, 14, 8, 10)
      )
    if (isTRUE(input$reverse_rows)) p <- p + scale_y_reverse()
    if (isTRUE(input$show_labels) && nrow(dat) <= 600) {
      p <- p + geom_text(aes(label = plot_id), size = 2.2, colour = "#17211b", check_overlap = TRUE)
    }

    ggplotly(p, tooltip = "text", source = "fieldmap") %>%
      layout(
        dragmode = "pan",
        hoverlabel = list(bgcolor = "white", bordercolor = "#d6dfda", font = list(color = "#17211b")),
        margin = list(l = 55, r = 20, b = 65, t = 10)
      ) %>%
      config(displaylogo = FALSE, scrollZoom = TRUE, modeBarButtonsToRemove = c("lasso2d", "select2d"))
  })

  # Show metadata for the currently selected plot on the map.
  output$plot_details <- renderUI({
    dat <- mapped_data()
    id <- selected_plot()
    if (is.null(id) || !id %in% dat$plot_id) {
      return(div(class = "plot-detail", p(class = "text-muted mb-0", "Click a plot in the field map to inspect its values and trial metadata.")))
    }
    x <- dat[match(id, dat$plot_id), , drop = FALSE]
    fields <- c("treatment", "genotype", "replicate", "block")
    fields <- fields[fields %in% names(x)]
    meta <- lapply(fields, function(f) tags$div(tags$span(class = "text-muted", paste0(tools::toTitleCase(f), ": ")), tags$strong(as.character(x[[f]]))))
    div(
      class = "plot-detail",
      div(class = "plot-id", x$plot_id),
      div(class = "text-muted small", paste("Row", x$row, "· Column", x$column)),
      div(class = "plot-detail-grid",
          tags$div(tags$span(class = "text-muted", "Observed: "), tags$strong(format_value(x$.raw))),
          tags$div(tags$span(class = "text-muted", "Displayed: "), tags$strong(format_value(x$.display))),
          meta)
    )
  })

  # Render a distribution plot for the currently selected layer.
  # Render a distribution plot for the currently selected layer.
  output$distribution <- renderPlotly({
    dat <- mapped_data() %>% filter(is.finite(.display))
    p <- ggplot(dat, aes(x = .display)) +
      geom_histogram(aes(y = after_stat(density)), bins = 24, fill = "#23855d", colour = "white", linewidth = 0.3) +
      geom_density(linewidth = 0.8, colour = "#17211b", adjust = 1.1) +
      labs(x = NULL, y = "Density") +
      theme_minimal(base_size = 11) +
      theme(panel.grid.minor = element_blank(), plot.margin = margin(8, 8, 8, 8))
    ggplotly(p, tooltip = c("x", "y")) %>% config(displaylogo = FALSE)
  })

  # Helper for row/column profile plots that show mean layer values along one axis.
  profile_plot <- function(dat, group_col, axis_label) {
    summary <- dat %>%
      filter(is.finite(.display)) %>%
      group_by(.data[[group_col]]) %>%
      summarise(mean = mean(.display), se = sd(.display) / sqrt(dplyr::n()), .groups = "drop")
    names(summary)[1] <- "position"
    p <- ggplot(summary, aes(position, mean, text = paste0(axis_label, ": ", position, "<br>Mean: ", format_value(mean)))) +
      geom_ribbon(aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se), alpha = 0.16, fill = "#23855d") +
      geom_line(linewidth = 0.9, colour = "#176b4d") +
      geom_point(size = 2, colour = "#176b4d") +
      labs(x = axis_label, y = "Mean displayed value") +
      theme_minimal(base_size = 12) +
      theme(panel.grid.minor = element_blank())
    ggplotly(p, tooltip = "text") %>% config(displaylogo = FALSE)
  }

  output$row_profile <- renderPlotly(profile_plot(mapped_data(), "row", "Field row"))
  output$column_profile <- renderPlotly(profile_plot(mapped_data(), "column", "Field column"))

  # Render the observed vs fitted trend diagnostic scatter plot.
  output$observed_fitted <- renderPlotly({
    dat <- analysed() %>% filter(is.finite(.raw), is.finite(.trend))
    p <- ggplot(dat, aes(.trend, .raw, key = plot_id,
                         text = paste0("<b>", plot_id, "</b><br>Observed: ", format_value(.raw), "<br>Trend: ", format_value(.trend)))) +
      geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#93a39a") +
      geom_point(alpha = 0.75, size = 2.4, colour = "#176b4d") +
      labs(x = "Estimated environmental trend", y = paste("Observed", input$trait)) +
      theme_minimal(base_size = 12) +
      theme(panel.grid.minor = element_blank())
    ggplotly(p, tooltip = "text") %>% config(displaylogo = FALSE)
  })

  # Render a QQ plot of residuals from the spatial model.
  output$residual_qq <- renderPlot({
    dat <- analysed()
    r <- dat$.residual[is.finite(dat$.residual)]
    validate(need(length(r) >= 4, "Not enough residuals."))
    qqnorm(r, pch = 19, col = "#176b4d", main = "Residual normal Q-Q", xlab = "Theoretical quantiles", ylab = "Residual quantiles")
    qqline(r, col = "#66756d", lwd = 2)
  }, res = 110)

  output$model_summary <- renderUI({
    dat <- analysed()
    ok <- is.finite(dat$.adjusted) & is.finite(dat$.trend)
    spatial_r2 <- if (sum(ok) > 2) cor(dat$.adjusted[ok], dat$.trend[ok])^2 else NA_real_
    rmse <- sqrt(mean(dat$.residual^2, na.rm = TRUE))
    div(
      class = "px-2 pb-2",
      p("The environmental surface is estimated with a two-dimensional LOESS smoother after optional treatment and design-factor adjustment."),
      tags$div(class = "d-flex justify-content-between py-1", span(class = "text-muted", "Spatial R²"), tags$strong(scales::percent(spatial_r2, accuracy = 0.1))),
      tags$div(class = "d-flex justify-content-between py-1", span(class = "text-muted", "Residual RMSE"), tags$strong(format_value(rmse))),
      tags$div(class = "d-flex justify-content-between py-1", span(class = "text-muted", "Complete plots"), tags$strong(sum(ok)))
    )
  })

  # KPI summaries shown at the top of the app.
  output$kpi_plots <- renderText(format(nrow(trial_data()), big.mark = ","))
  output$kpi_mean <- renderText(format_value(mean(mapped_data()$.raw, na.rm = TRUE)))
  output$kpi_range <- renderText({
    x <- mapped_data()$.display
    paste0(format_value(min(x, na.rm = TRUE)), " – ", format_value(max(x, na.rm = TRUE)))
  })
  output$kpi_spatial <- renderText({
    dat <- analysed()
    ok <- is.finite(dat$.adjusted) & is.finite(dat$.trend)
    if (sum(ok) > 2) scales::percent(cor(dat$.adjusted[ok], dat$.trend[ok])^2, accuracy = 1) else "NA"
  })

  output$data_status <- renderUI({
    label <- if (is.null(input$file)) "Demo data" else input$file$name
    div(class = "status-pill", span(class = "status-dot"), label)
  })

  output$data_table <- renderDT({
    dat <- analysed()
    names(dat)[names(dat) == ".raw"] <- paste0(input$trait, "_observed")
    names(dat)[names(dat) == ".adjusted"] <- paste0(input$trait, "_adjusted")
    names(dat)[names(dat) == ".trend"] <- paste0(input$trait, "_spatial_trend")
    names(dat)[names(dat) == ".residual"] <- paste0(input$trait, "_spatial_residual")
    datatable(
      dat,
      rownames = FALSE,
      filter = "top",
      extensions = c("Buttons", "Scroller"),
      options = list(
        dom = "Bfrtip", buttons = c("copy", "csv"), pageLength = 15,
        scrollX = TRUE, deferRender = TRUE, scroller = TRUE, scrollY = 520
      )
    )
  })

  output$download_processed <- downloadHandler(
    filename = function() paste0("FSV_", input$trait, "_processed.csv"),
    content = function(file) {
      dat <- analysed()
      names(dat)[names(dat) == ".raw"] <- paste0(input$trait, "_observed")
      names(dat)[names(dat) == ".adjusted"] <- paste0(input$trait, "_adjusted")
      names(dat)[names(dat) == ".trend"] <- paste0(input$trait, "_spatial_trend")
      names(dat)[names(dat) == ".residual"] <- paste0(input$trait, "_spatial_residual")
      write.csv(dat, file, row.names = FALSE, na = "")
    }
  )
}

shinyApp(ui, server)

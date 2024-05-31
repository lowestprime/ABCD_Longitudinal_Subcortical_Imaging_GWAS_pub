library(shiny)
library(dplyr)
library(gt)

server <- function(input, output) {
  smri.R5.1.baseline.y2 <- readRDS("smri.R5.1.baseline.y2.lite.rds")
  smri.R5.1.all <- readRDS("smri.R5.1.all.lite.rds")
  data_selected <- reactive({
    data <- if (input$dataset == "all") {
      smri.R5.1.all
    } else {
      smri.R5.1.baseline.y2
    }
    if (!is.null(input$ethnicityFilter) && !all(input$ethnicityFilter == "All")) {
      data <- data %>% filter(ethnicity %in% input$ethnicityFilter)
    } else {
      data <- data %>% filter(ethnicity %in% c("AFR", "AMR", "EUR"))
    }
    data
  })
  output$plot <- renderPlotly({
    data_summary <- data_selected() %>%
      distinct(src_subject_id, .keep_all = TRUE) %>%
      group_by(ethnicity, sex) %>%
      summarise(n = n_distinct(src_subject_id), .groups = 'drop') %>%
      mutate(Total_by_ethnicity = sum(n), .by = ethnicity) %>%
      ungroup() %>%
      group_by(ethnicity) %>%
      mutate(Total = sum(n))
    dataset_label <- if (input$dataset == "all") {"All Timepoints"} else {"Baseline Y2"}
    fig <- plot_ly(data_summary, x = ~ethnicity, y = ~n, type = 'bar', color = ~sex,
                   colors = c('F' = '#FF69B4', 'M' = '#1E90FF', 'NA' = '#4B0082'), 
                   hoverinfo = 'text',
                   hovertemplate = 'Count: %{y}<br>Total by Ethnicity: %{customdata}<extra></extra>',
                   customdata = ~Total_by_ethnicity,
                   textposition = 'none') %>%
      layout(title = paste("Ethnicity and Sex Counts (", dataset_label, ")", sep = ""),
             xaxis = list(title = "Ethnicity", titlefont = list(size = 16)),
             yaxis = list(title = "Count", titlefont = list(size = 16)),
             barmode = 'stack',
             margin = list(l = 50, r = 50, t = 70, b = 50),  # Adjusted margins for better use of space
             legend = list(x = 1.05, y = 0.5, orientation = "v", font = list(size = 14)),  # Vertically oriented legend
             plot_bgcolor = 'white')
    fig
  })
  output$table <- render_gt({
    data_summary <- data_selected() %>%
      distinct(src_subject_id, .keep_all = TRUE) %>%
      group_by(ethnicity, sex) %>%
      summarise(n = n_distinct(src_subject_id), .groups = 'drop')
    gt(data_summary) %>%
      grand_summary_rows(
        columns = c(n),
        fns = list(Total = ~sum(.))
      ) %>%
      tab_style(
        style = cell_fill(color = "transparent"),
        locations = cells_body(columns = c(n))
      )
  })
}
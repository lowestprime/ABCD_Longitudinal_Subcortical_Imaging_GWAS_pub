library(shiny)
library(dplyr)
library(gt)

smri.R5.1.baseline.y2.ROC.long <- readRDS("smri.R5.1.baseline.y2.ROC.long.rds")
smri.R5.1.all.ROC.long <- readRDS("smri.R5.1.all.ROC.long.rds")

server <- function(input, output) {
  data_selected <- reactive({
    if (input$dataset == "all") {
      data <- smri.R5.1.all.ROC.long
      data <- data %>% filter(roi == input$region & Time_Comparison == input$timeComparison)
    } else {
      data <- smri.R5.1.baseline.y2.ROC.long
      data <- data %>% filter(roi == paste0(input$region, "_ROC0_2"))
    }
    data
  })
  
  output$summary_table <- render_gt({
    data_summary <- data_selected() %>%
      summarise(
        avg = if (all(is.na(Value))) NA else mean(Value, na.rm = TRUE),
        max = if (all(is.na(Value))) NA else max(Value, na.rm = TRUE),
        min = if (all(is.na(Value))) NA else min(Value, na.rm = TRUE),
        std_dev = if (all(is.na(Value))) NA else sd(Value, na.rm = TRUE),
        count = n()
      )
    
    gt(data_summary) %>%
      tab_header(
        title = "Summary Statistics"
      ) %>%
      cols_label(
        avg = "Average",
        max = "Maximum",
        min = "Minimum",
        std_dev = "Standard Deviation",
        count = "Count"
      ) %>%
      tab_style(
        style = cell_fill(color = "transparent"),
        locations = cells_body(columns = c(avg, max, min, std_dev, count))
      )
  })
  
  output$detailed_table <- render_gt({
    if (input$dataset == "all") {
      data_sorted <- smri.R5.1.all.ROC.long %>%
        group_by(roi, Time_Comparison) %>%
        summarise(
          avg_value = mean(Value, na.rm = TRUE),
          .groups = 'drop'
        ) %>%
        arrange(if (input$sortBy == "asc") avg_value else desc(avg_value))
      
      gt(data_sorted) %>%
        opt_interactive() %>%
        tab_header(
          title = "Detailed Table of Region ROCs"
        ) %>%
        cols_label(
          roi = "Region",
          Time_Comparison = "Time Comparison",
          avg_value = "Average Rate of Change"
        ) %>%
        tab_style(
          style = cell_fill(color = "transparent"),
          locations = cells_body(columns = c(avg_value))
        )
    } else {
      data_sorted <- smri.R5.1.baseline.y2.ROC.long %>%
        group_by(roi) %>%
        summarise(
          avg_value = mean(Value, na.rm = TRUE),
          .groups = 'drop'
        ) %>%
        arrange(if (input$sortBy == "asc") avg_value else desc(avg_value))
      
      gt(data_sorted) %>%
        opt_interactive() %>%
        tab_header(
          title = "Detailed Table of Region ROCs"
        ) %>%
        cols_label(
          roi = "Region",
          avg_value = "Average Rate of Change"
        ) %>%
        tab_style(
          style = cell_fill(color = "transparent"),
          locations = cells_body(columns = c(avg_value))
        )
    }
  })
}

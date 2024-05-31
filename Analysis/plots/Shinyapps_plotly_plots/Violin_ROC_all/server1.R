library(shiny)
library(plotly)

server1 <- function(input, output) {
  # Load data
  smri.R5.1.all.ROC.long <- readRDS("smri.R5.1.all.ROC.long.rds");
  # run
  output$roi_plot <- renderPlotly({
    filtered_data <- smri.R5.1.all.ROC.long[smri.R5.1.all.ROC.long$roi == input$roi_select, ]
    
    p <- plot_ly(data = filtered_data, type = "violin",
                 x = ~Time_Comparison, y = ~Value,
                 color = ~Time_Comparison,
                 hoverinfo = 'text',
                 text = ~paste('ROI:', roi, '<br>Value:', Value, '<br>Time Comparison:', Time_Comparison)) %>%
      layout(
        title = paste("Distribution of Rate of Change for", input$roi_select, "(all timepoint subjects)"),
        yaxis = list(title = "Rate of Change (mm³/yr)", titlefont = list(size = 16), tickfont = list(size = 12), standoff = 20),
        xaxis = list(title = "Time Comparison", tickangle = 0, tickfont = list(size = 12)),
        legend = list(title = "Time Comparison", orientation = "v", x = 1.05, y = 0.5, xanchor = 'left', yanchor = 'middle', font = list(size = 16)),
        margin = list(l = 40, r = 40, t = 50, b = 50)
      )
    return(p)
  })
}
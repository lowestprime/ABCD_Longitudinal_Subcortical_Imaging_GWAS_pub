library(shiny)
library(plotly)

server1 <- function(input, output) {
  # Load data
  smri.R5.1.baseline.y2.ROC.long <- readRDS("smri.R5.1.baseline.y2.ROC.long.rds");
  # run
  output$roi_plot <- renderPlotly({
    filtered_data <- smri.R5.1.baseline.y2.ROC.long[smri.R5.1.baseline.y2.ROC.long$roi == input$roi_select, ]
    
    p <- plot_ly(data = filtered_data, type = "violin",
                 x = ~roi,  # Using ROI as the x-axis
                 y = ~Value,  # Value of the measurement
                 hoverinfo = 'text',
                 text = ~paste('ROI:', roi, '<br>Value:', Value)) %>%
      layout(
        title = paste("Distribution of Rate of Change for", input$roi_select, "(baseline + y2 only subjects)"),
        yaxis = list(title = "Rate of Change (mm³/yr)", titlefont = list(size = 16), tickfont = list(size = 12), standoff = 20),
        xaxis = list(title = "ROI", tickangle = 0, tickfont = list(size = 12)),  # Adjusting for better readability
        margin = list(l = 40, r = 40, t = 50, b = 50)
      )
    return(p)
  })
}
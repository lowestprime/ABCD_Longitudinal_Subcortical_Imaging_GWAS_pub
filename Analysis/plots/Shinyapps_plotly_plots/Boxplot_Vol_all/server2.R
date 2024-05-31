library(shiny)
library(plotly)

server2 <- function(input, output) {
  # Load data
  smri.R5.1.all.long <- readRDS("smri.R5.1.all.long.rds");
  # run
  output$boxplot <- renderPlotly({
    filtered_data <- smri.R5.1.all.long %>%
      filter(volume_type == input$volume_type_select)
    p <- plot_ly(data = filtered_data, type = 'box',
                 x = ~factor(timepoint), y = ~volume,
                 color = ~factor(timepoint),
                 colors = c('#636EFA', '#EF553B', '#00CC96')) %>%
      layout(
        title = paste("Distribution of Volume Measurements for", input$volume_type_select, "(all timepoint subjects)"),
        yaxis = list(title = "Volume (mm³)", titlefont = list(size = 16), standoff = 30),
        xaxis = list(title = "Timepoint", titlefont = list(size = 16)),
        margin = list(l = 50, r = 50, t = 50, b = 60),
        legend = list(font = list(size = 16))
      )
    return(p)
  })
}
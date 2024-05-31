library(shiny)
library(plotly)

# Load data
smri.R5.1.baseline.y2.long <- readRDS("smri.R5.1.baseline.y2.long.rds");

# ui
ui2 <- fluidPage(
  tags$head(
    tags$style(HTML("
            body { text-align: center; margin: auto; }
            .shiny-input-container { margin: 20px auto; width: 80%; }
            .shiny-plot-output { margin: auto; width: 90%; }
            .js-plotly-plot .plotly .legend { font-size: 16px; } /* Increase legend font size */
        "))
  ),
  titlePanel("Interactive Boxplots of SCS ROI Volume Measurements", windowTitle = "Volume Measurements"),
  plotlyOutput("boxplot", height = "650px"),
  selectInput("volume_type_select", "Select Volume Type:",
              choices = unique(smri.R5.1.baseline.y2.long$volume_type),
              selected = unique(smri.R5.1.baseline.y2.long$volume_type)[1])
)
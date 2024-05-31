library(shiny)
library(plotly)

# Load data
smri.R5.1.all.ROC.long <- readRDS("smri.R5.1.all.ROC.long.rds");

# ui
ui1 <- fluidPage(
  tags$head(
    tags$style(HTML("
            body { text-align: center; margin: auto; }
            .shiny-input-container { margin: 20px auto; width: 50%; } /* Center and adjust width of dropdown */
            .shiny-plot-output { margin: auto; width: 90%; } /* Adjust plot width and centering */
            .js-plotly-plot .plotly .legend { font-size: 16px; } /* Increase legend font size */
        "))
  ),
  titlePanel("Interactive Violin Plots of SCS ROI Rate of Change", windowTitle = "Rate of Change Analysis"),
  selectInput("roi_select", "Select ROI:", choices = unique(smri.R5.1.all.ROC.long$roi)),
  plotlyOutput("roi_plot", height = "600px")  # Adjust the plot height
)
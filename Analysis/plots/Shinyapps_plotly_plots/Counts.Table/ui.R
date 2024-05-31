library(shiny)
library(plotly)
library(gt)

# ui
ui <- fluidPage(
  titlePanel("Data Summary Dashboard"),
  sidebarLayout(
    sidebarPanel(
      selectInput("dataset", "Choose a dataset:",
                  choices = c("All Timepoints" = "all", "Baseline Y2" = "baseline.y2")),
      checkboxGroupInput("ethnicityFilter", "Filter by Ethnicity:",
                         choices = c("AFR", "AMR", "EUR"))
    ),
    mainPanel(
      div(style = "display: flex; flex-direction: row;",  # Using flexbox for layout
          div(plotlyOutput("plot", height = "600px"), style = "flex: 75%;"),  # Adjusted for 75% width
          div(gt_output("table"), style = "flex: 25%;")  # Adjusted for 25% width
      )
    )
  )
)
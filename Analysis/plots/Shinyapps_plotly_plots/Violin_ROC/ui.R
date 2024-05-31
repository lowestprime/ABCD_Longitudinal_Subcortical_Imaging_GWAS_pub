smri.R5.1.baseline.y2.ROC.long <- readRDS("smri.R5.1.baseline.y2.ROC.long.rds")
smri.R5.1.all.ROC.long <- readRDS("smri.R5.1.all.ROC.long.rds")

ui <- fluidPage(
  titlePanel("Summary of Percent Rates of Change"),
  sidebarLayout(
    sidebarPanel(
      selectInput("dataset", "Select Dataset", 
                  choices = list("All" = "all", "Baseline" = "baseline")),
      selectInput("region", "Select Region", 
                  choices = unique(gsub("_ROC0_2", "", smri.R5.1.baseline.y2.ROC.long$roi))),
      conditionalPanel(
        condition = "input.dataset == 'all'",
        selectInput("timeComparison", "Select Time Comparison", 
                    choices = unique(smri.R5.1.all.ROC.long$Time_Comparison))
      ),
      selectInput("sortBy", "Sort By", 
                  choices = list("Ascending" = "asc", "Descending" = "desc"),
                  selected = "asc")
    ),
    mainPanel(
      titlePanel("Interactive Violin Plots of SCS ROI Rate of Change", windowTitle = "Rate of Change Analysis"),
      selectInput("roi_select", "Select ROI:", choices = unique(smri.R5.1.all.ROC.long$roi)),
      plotlyOutput("roi_plot", height = "600px")  # Adjust the plot height
    )
  )
)
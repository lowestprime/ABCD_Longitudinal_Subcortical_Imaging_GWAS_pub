library(shiny)

source("ui2.R")
source("server2.R")

# Load data
smri.R5.1.baseline.y2.long <- readRDS("smri.R5.1.baseline.y2.long.rds");

shinyApp(ui = ui2, server = server2)
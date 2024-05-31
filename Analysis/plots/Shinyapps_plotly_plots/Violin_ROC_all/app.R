library(shiny)

source("ui1.R")
source("server1.R")

# Load data
smri.R5.1.all.ROC.long <- readRDS("smri.R5.1.all.ROC.long.rds");

shinyApp(ui = ui1, server = server1)
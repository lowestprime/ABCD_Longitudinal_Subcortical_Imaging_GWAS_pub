library(shiny)

source("ui2.R")
source("server2.R")

# Load data
smri.R5.1.all.long <- readRDS("smri.R5.1.all.long.rds");

shinyApp(ui = ui2, server = server2)
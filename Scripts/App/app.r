# nolint start

library(shiny)
library(ggplot2)

PROJECT <- "/home/daniele/Desktop/IV_course/II_semester/TF_analysis/"
source(paste0(PROJECT, "Scripts/App/ui.R"))
source(paste0(PROJECT, "Scripts/App/server.R"))
shinyApp(ui = ui, server = server)

# nolint end
library(shiny)
library(data.table)

server <- function(input, output, session) {
  output$files <- renderTable(input$upload)
  file <- input$upload
  read.csv(file$datapath, header = input$header)
}
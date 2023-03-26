# nolint start

ui <- fluidPage( 
  titlePanel("Transcription factor analysis dashboard"),
  sidebarPanel(
    fileInput(
        inputId = "bigbed",
        label = "Provide BigBed file(-s):",
        multiple = TRUE,
        buttonLabel = "Browse files",
        placeholder = "No file selected"
    ),
    selectInput(
        inputId = "tf_options",
        label = "Specify transcription factor:",
        choices = c("Tbx5", "GATA3", "Tcf21"),
        selected = "GATA3"
    ),
    fileInput(
        inputId = "pwm",
        label = "Provide PWM matrix:",
        multiple = FALSE,
        buttonLabel = "Browse files",
        placeholder = "No file selected"
    )
  ),
  mainPanel(
    plotOutput(outputId = "table")
  )
    

)

# nolint end
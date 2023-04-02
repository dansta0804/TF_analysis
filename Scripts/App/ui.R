# nolint start
library(shinydashboard)


ui <- dashboardPage(
  
  dashboardHeader(
    title = "Transcription factor analysis dashboard",
    titleWidth = 350
  ),
  dashboardSidebar(
    width = 350,
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
    radioButtons(
        inputId = "genome",
        label = "Specify genome:",
        choices = c("Homo sapiens" = "hg", "Mus musculus" = "mm",
                    "Danio rerio" = "dr")
    ),
    fileInput(
        inputId = "pwm",
        label = "Provide PWM matrix:",
        multiple = FALSE,
        buttonLabel = "Browse files",
        placeholder = "No file selected"
  )),
  dashboardBody(
     tags$head(tags$style(HTML("
      .skin-blue .main-header .navbar {
        background-color: #cf882a;
      }
      .skin-blue .main-header .logo {
        background-color: #dda153;
        color: black;
        font-weight: bolder;
        border-bottom: 0 solid transparent;
      }
    "))),
    fluidRow(
      box(plotOutput(outputId = "plot1")),
      box(plotOutput(outputId = "plot2")),
      box(plotOutput(outputId = "plot3")),
      box(plotOutput(outputId = "plot4")),
      box(plotOutput(outputId = "plot5"))
    )
  )
)




# nolint end
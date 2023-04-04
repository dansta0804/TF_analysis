# nolint start
library(shinydashboard)

# VARIANTAS 1:

# sidebar <- dashboardSidebar(
#   width = 350,
#   sidebarMenu(
#     menuItem(
#       "Duomenų kokybė", tabName = "data_quality", icon = icon("circle-check", class = NULL, lib = "font-awesome"),
#       fileInput(
#         inputId = "bigbed", label = "Provide BigBed file(-s):", multiple = TRUE,
#         buttonLabel = "Browse files", placeholder = "No file selected"
#       )
#     ),
#     menuItem(
#       "Analizės", tabName = "analysis", icon = icon("magnifying-glass-chart"),
#       selectInput(
#         inputId = "tf_options", label = "Specify transcription factor:",
#         choices = c("Tbx5", "GATA3", "Tcf21"), selected = "GATA3"
#       ),
#       fileInput(
#         inputId = "pwm", label = "Provide PWM matrix:", multiple = FALSE,
#         buttonLabel = "Browse files", placeholder = "No file selected"
#       )
#     ),
#     menuItem(
#       "Taikinių spėjimas", tabName = "target_prediction",
#       icon = icon("question", class = NULL, lib = "font-awesome"),
#       radioButtons(
#         inputId = "genome", label = "Specify genome:",
#         choices = c("Homo sapiens" = "hg", "Mus musculus" = "mm",
#                     "Danio rerio" = "dr")
#       )
#     )
#   )
# )

# body <- dashboardBody(
#   tags$head(tags$style(HTML("
#       .skin-blue .main-header .navbar {
#         background-color: #cf882a;
#       }
#       .skin-blue .main-header .logo {
#         background-color: #dda153;
#         color: black;
#         font-weight: bolder;
#         border-bottom: 0 solid transparent;
#       }
#       .main-sidebar {
#         font-size: 17px;
#       }
#     "))),
#   tabItems(
#     tabItem(tabName = "data_quality",
#     plotOutput(outputId = "plot1")
    
#     )
#   )
# )

# ui <- dashboardPage(
#   dashboardHeader(
#     title = "Transcription factor analysis dashboard",
#     titleWidth = 350),
#   sidebar,body
# )


# VARIANTAS 2:  

ui <- dashboardPage(
  dashboardHeader(
    title = "Transcription factor analysis dashboard",
    titleWidth = 350),

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
      box(
        title = "Pikų skaičius mėginiuose", solidHeader = TRUE,
        status = "primary", plotOutput(outputId = "plot1")
      ),
      box(
        title = "Pikų skaičius chromosomose", solidHeader = TRUE,
        status = "primary", plotOutput(outputId = "plot2")
      ),
      box(
        title = "Mėginių panašumas", solidHeader = TRUE,
        status = "primary", plotOutput(outputId = "plot3")
      ),
      box(
        title = "Transkripcijos faktoriaus motyvas", solidHeader = TRUE,
        status = "primary", plotOutput(outputId = "plot4")
      ),
      box(
        title = "PWM matricos atitikimai", solidHeader = TRUE,
        status = "primary", plotOutput(outputId = "plot5")
      )
    )
  )
)


# nolint end
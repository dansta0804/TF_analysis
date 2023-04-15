# nolint start
library(shinydashboard)

PROJECT <- "/home/daniele/Desktop/IV_course/II_semester/TF_analysis/"
source(paste0(PROJECT, "Scripts/App/data_quality.R"))
source(paste0(PROJECT, "Scripts/App/server.R"))


# PRIE DUOMENŲ KOKYBĖS TURI BŪTI:
#   PIKŲ SKAIČIUS (PATIKIMŲ IR NEPATIKIMŲ);
#   IŠSIDĖSTYMAS APIE TSS;
#   GENOMINĖ DISTRIBUCIJA;
#   KIEK TF MOTYVŲ NUSTATYTA;
#   PAVAIZDUOTI, AR MOTYVAI IŠSIDĖSTĘ PIKO CENTRO REGIONE;
#   PIKŲ SKAIČIUS CHROMOSOMOSE;
#   MĖGINIŲ PANAŠUMAS;

# VARIANTAS 2:

# body <- dashboardBody(
#     fluidRow(
#         tabBox(
#             title = NULL, width = 12,
#             # The id lets us use input$tabset1 on the server to find the current tab
#             id = "tabset1", height = "250px",
#             tabPanel("Victim", "Victim tab"),
#             tabPanel("Trafficker", "Trafficker tab")
#         )
#     ),
#     fluidRow(infoBoxOutput("tabset1Selected"))
# )

ui <- navbarPage("ChIP sekoskaitos analizės", theme = shinytheme("cosmo"),
tags$style(HTML("
  body {
    font-family: 'Source Sans Pro',sans-serif;
    font-size: 15px;
    line-height: 1.42857143;
    color: #333333;
    background-color: #ffffff;
  }
  .navbar-default {
    background-color: #683030;
  }
  .navbar-default .navbar-nav>.active>a, .navbar-default
  .navbar-nav>.active>a:hover, .navbar-default .navbar-nav>.active>a:focus {
    color: #ffffff;
    background-color: #340a0a;
  }
  .btn-default {
    color: #ffffff;
    background-color: #a75f21;
    border-color: #723e0d;
  }
  a {
    color: #700c0c;
    font-weight: bold;
  }
  a:hover, a:focus {
    color: #a75f21;
    text-decoration: underline;
  }
  p {
    margin: 20px;
    text-align: justify;
    font-size: 22px;
  }
")),
  tabPanel("Duomenų kokybė",
    sidebarLayout(
        sidebarPanel(
            width = 4,
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
          )
        ),
        mainPanel(
          width = 8,
      #         fluidRow(
      #   box(
      #     title = "Pikų skaičius mėginiuose", solidHeader = TRUE,
      #     status = "primary", plotOutput(outputId = "plot1")
      #   ),
      #   box(
      #     title = "Pikų skaičius chromosomose", solidHeader = TRUE,
      #     status = "primary", plotOutput(outputId = "plot2")
      #   ),
      #   box(
      #     title = "Mėginių panašumas", solidHeader = TRUE,
      #     status = "primary", plotOutput(outputId = "plot3")
      #   ),
      #   box(
      #     title = "Transkripcijos faktoriaus motyvas", solidHeader = TRUE,
      #     status = "primary", plotOutput(outputId = "plot4")
      #   ),
      #   box(
      #     title = "PWM matricos atitikimai", solidHeader = TRUE,
      #     status = "primary", plotOutput(outputId = "plot5")
      #   ),
      #   box(
      #     width = 10, title = "Genominis pasiskirstymas", solidHeader = TRUE,
      #     status = "primary", DT::dataTableOutput(outputId = "table1"),
      #     downloadButton("downloadData", "Download")
      #   )
      # )
      tabsetPanel(
        
        tabPanel("Pikų skaičius mėginiuose", 
          p("Description goes here..."),
          box(
            width = 12, 
            withLoader(plotOutput("plot1"), type = "html", loader = "dnaspin")
          )
        ),
        tabPanel("Pikų skaičius chromosomose",
          p("Description goes here..."),
          box(
            width = 12,
            withLoader(plotOutput("plot2"), type = "html", loader = "dnaspin")
          )
        ),
        tabPanel("Mėginių panašumas",
          p("Description goes here..."),
          box(
            width = 12,
            withLoader(plotOutput("plot3"), type = "html", loader = "dnaspin")
          ))
      )
        )
    )
  ),
  tabPanel("Analizės",
          h2("Trafficker tab")
  ),
  tabPanel("Taikinių spėjimas")
)

# ui <- dashboardPage(
#   dashboardHeader(
#     title = "Transcription factor analysis dashboard",
#     titleWidth = 350),

#   dashboardSidebar(
#     width = 350,
#     fileInput(
#         inputId = "bigbed",
#         label = "Provide BigBed file(-s):",
#         multiple = TRUE,
#         buttonLabel = "Browse files",
#         placeholder = "No file selected"
#     ),
#     selectInput(
#         inputId = "tf_options",
#         label = "Specify transcription factor:",
#         choices = c("Tbx5", "GATA3", "Tcf21"),
#         selected = "GATA3"
#     ),
#     radioButtons(
#         inputId = "genome",
#         label = "Specify genome:",
#         choices = c("Homo sapiens" = "hg", "Mus musculus" = "mm",
#                     "Danio rerio" = "dr")
#     ),
#     fileInput(
#         inputId = "pwm",
#         label = "Provide PWM matrix:",
#         multiple = FALSE,
#         buttonLabel = "Browse files",
#         placeholder = "No file selected"
#   )),
#   dashboardBody(
#      tags$head(tags$style(HTML("
#       .skin-blue .main-header .navbar {
#         background-color: #cf882a;
#       }
#       .skin-blue .main-header .logo {
#         background-color: #dda153;
#         color: black;
#         font-weight: bolder;
#         border-bottom: 0 solid transparent;
#       }
#     "))),
#     fluidRow(
#       box(
#         title = "Pikų skaičius mėginiuose", solidHeader = TRUE,
#         status = "primary", plotOutput(outputId = "plot1")
#       ),
#       box(
#         title = "Pikų skaičius chromosomose", solidHeader = TRUE,
#         status = "primary", plotOutput(outputId = "plot2")
#       ),
#       box(
#         title = "Mėginių panašumas", solidHeader = TRUE,
#         status = "primary", plotOutput(outputId = "plot3")
#       ),
#       box(
#         title = "Transkripcijos faktoriaus motyvas", solidHeader = TRUE,
#         status = "primary", plotOutput(outputId = "plot4")
#       ),
#       box(
#         title = "PWM matricos atitikimai", solidHeader = TRUE,
#         status = "primary", plotOutput(outputId = "plot5")
#       ),
#       box(
#         width = 10, title = "Genominis pasiskirstymas", solidHeader = TRUE,
#         status = "primary", DT::dataTableOutput(outputId = "table1"),
#         downloadButton("downloadData", "Download")
#       )
#     )


    
#   )
# )

# ui <- shinyUI(  
#   navbarPage(title = "Transkripcijos faktorių analizė",
#   tabPanel(
#     "Duomenų kokybė",
#     box(
#       title = "Pikų skaičius mėginiuose", solidHeader = TRUE,
#       status = "primary", plotOutput(outputId = "plot1")
#     ),
#     box(
#       title = "Pikų skaičius chromosomose", solidHeader = TRUE,
#       status = "primary", plotOutput(outputId = "plot2")
#     ),
#     box(
#       title = "Mėginių panašumas", solidHeader = TRUE,
#       status = "primary", plotOutput(outputId = "plot3")
#     ),
#     box(
#       title = "Transkripcijos faktoriaus motyvas", solidHeader = TRUE,
#       status = "primary", plotOutput(outputId = "plot4")
#     ),
#     box(
#       title = "PWM matricos atitikimai", solidHeader = TRUE,
#       status = "primary", plotOutput(outputId = "plot5")
#     )
#   ),
#   tabPanel(
#     "Analizės", p("labasĄ")
#   ),
#   tabPanel(
#     "Taikinių spėjimas", p("labas!")
#   )
#   )
# )


# nolint end
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

# Pridėti lentelę, kur naudotojas gali pasirinkti mėginį, kurio duomenų kokybę
# nori įvertinti. Toliau išsidėstymą apie TSS, genominę distribuciją vaizduoti
# tik pasirinktam mėginiui! Antradienį!

ui <- navbarPage("ChIP sekoskaitos analizės", theme = shinytheme("cosmo"),
  tags$style(HTML("
    body {
      font-family: sans-serif;
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
      margin: 18px;
      text-align: justify;
      font-size: 22px;
    }
    table.dataTable.hover>tbody>tr.selected:hover>*,table.dataTable.display>tbody>tr.selected:hover>* {
      box-shadow: inset 0 0 0 9999px rgb(77 6 6 / 90%)
    }
    table.dataTable.stripe>tbody>tr.odd.selected>*, table.dataTable.display>tbody>tr.odd.selected>* {
      box-shadow: inset 0 0 0 9999px rgb(77 6 6 / 90%)
    }
    table.dataTable.stripe>tbody>tr.odd.selected>*,table.dataTable.display>tbody>tr.odd.selected>* {
      box-shadow: inset 0 0 0 9999px rgb(77 6 6 / 90%)
      color: white
      font-weight: bold
    }
    table.dataTable tbody tr.selected>* {
      box-shadow: inset 0 0 0 9999px rgb(77 6 6 / 90%);
      color: white
      text-decoration: none
      font-weight: bold
    }
    table.dataTable tbody tr.selected>* {
      color: white;
      font-weight: bold;
    }
    .DNA:after {
      letter-spacing: 5px;
      display: inline-block;
      padding-top: 100px;
      color: black;
      font-family: arial;
      text-align: center;
      font-size: 36pt;
      content: 'Skaičiuojama...';
    }
    .strain:before {
      content: '';
      display: block;
      position: absolute;
      width: 18px;
      height: 18px;
      border-radius: 50%;
      background-color: #4f0707;
      left: 0;
      margin-top: -9px;
      margin-left: -9px;
      -webkit-animation-name: strain-size-left;
              animation-name: strain-size-left;
      -webkit-animation-duration: 1.7s;
              animation-duration: 1.7s;
      -webkit-animation-timing-function: linear;
              animation-timing-function: linear;
      -webkit-animation-iteration-count: infinite;
              animation-iteration-count: infinite;
      -webkit-animation-play-state: running;
              animation-play-state: running;
    }
    .strain:after {
      content: '';
      display: block;
      position: absolute;
      width: 18px;
      height: 18px;
      border-radius: 50%;
      background-color: #d68325;
      right: 0;
      margin-top: -9px;
      margin-right: -9px;
      -webkit-animation-name: strain-size-right;
              animation-name: strain-size-right;
      -webkit-animation-duration: 1.7s;
              animation-duration: 1.7s;
      -webkit-animation-timing-function: linear;
              animation-timing-function: linear;
      -webkit-animation-iteration-count: infinite;
              animation-iteration-count: infinite;
      -webkit-animation-play-state: running;
              animation-play-state: running;
    }
  ")),
  # GALBŪT ČIA REIKIA NAUJO TAB'O, KURIAME REIKĖTŲ ĮKELTI VISUS DUOMENIS?
  tabPanel("Duomenų kokybė",
    sidebarLayout(
      sidebarPanel(
        width = 4,
        fileInput(
          inputId = "bigbed",
          label = "Įkelkite BED formato failą (-us):",
          multiple = TRUE,
          buttonLabel = "Ieškoti failo",
          placeholder = "Failas nepasirinktas"
        ),
        radioButtons(
          inputId = "genome",
          label = "Pasirinkite genomą:",
          choices = c("Homo sapiens" = "hg", "Mus musculus" = "mm",
                      "Danio rerio" = "dr")
        )
      ),
      mainPanel(
        width = 8,
        tabsetPanel(
          tabPanel("Pikų skaičius mėginiuose", 
            p("Pateiktoje stulpelinėje diagramoje pavaizduota, kiek yra pikų
              kiekviename pateiktame BED formato faile:"),
              DT::dataTableOutput('samples'),
            box(
              width = 12, 
              withLoader(plotOutput("plot1"), type = "html", loader = "dnaspin")
            )
          ),
          tabPanel("Pikų skaičius chromosomose",
            p("Pateiktose stulpelinėse diagramose pavaizduota, kaip pikų skaičius
              pasiskirstęs skirtingose chromosomose:"),
            box(
              width = 12,
              withLoader(plotOutput("plot2"), type = "html", loader = "dnaspin")
            )
          ),
          tabPanel("Mėginių panašumas",
            p("Pateiktame spalvų intensyvumo grafike pavaizduota, kokia pikų
              dalis (procentiškai) sutampa tarp skirtingų mėginių:"),
            box(
              width = 12,
              withLoader(plotOutput("plot3"), type = "html", loader = "dnaspin")
            )
          ),
          tabPanel("Genominė distribucija",
            p("DESCRIPTION GOES HERE..."),
            box(
              width = 12,
              withLoader(DT::dataTableOutput("table2"), type = "html",
                         loader = "dnaspin")
            ) 
          )
        )
      )
    )
  )
  # tabPanel("Analizės",
  #   sidebarLayout(
  #     sidebarPanel(
  #       width = 4,
  #       selectInput(
  #         inputId = "tf_options",
  #         label = "Pasirinkite transkripcijos faktorių:",
  #         choices = c("Tbx5", "GATA3", "Tcf21"),
  #         selected = "GATA3"
  #       ),
  #       fileInput(
  #         inputId = "pwm",
  #         label = "Įkelkite transkripcijos faktoriaus PWM matricą:",
  #         multiple = FALSE,
  #         buttonLabel = "Ieškoti failo",
  #         placeholder = "Failas nepasirinktas"
  #       ),
  #       radioButtons(
  #         inputId = "genome",
  #         label = "Pasirinkite genomą:",
  #         choices = c("Homo sapiens" = "hg", "Mus musculus" = "mm",
  #                     "Danio rerio" = "dr")
  #       )
  #     ),
  #     mainPanel(
  #       width = 8,
  #       tabsetPanel(
  #         tabPanel("Transkripcijos faktoriaus motyvas", 
  #           p("Description goes here..."),
  #           box(
  #             width = 12, 
  #             withLoader(plotOutput("plot4"), type = "html", loader = "dnaspin")
  #           )
  #         ),
  #         tabPanel("PWM matricos atitikimai",
  #           p("Description goes here..."),
  #           box(
  #             width = 12,
  #             withLoader(plotOutput("plot5"), type = "html", loader = "dnaspin")
  #           )
  #         ),
  #         tabPanel("Genominis pasiskirstymas",
  #           p("Description goes here..."),
  #           box(
  #             width = 12,
  #             withLoader(DT::dataTableOutput(outputId = "table1"),
  #                        type = "html", loader = "dnaspin"),
  #             downloadButton("downloadData", "Download"),
  #           #   DT::dataTableOutput('x3'),
  #           # verbatimTextOutput('x4')
  #           )
  #         )
  #       )
  #     )
  #   )
  # ),
  # tabPanel("Taikinių spėjimas",
  #   sidebarLayout(
  #     sidebarPanel(
  #       width = 4,
  #       radioButtons(
  #         inputId = "genome",
  #         label = "Pasirinkite genomą:",
  #         choices = c("Homo sapiens" = "hg", "Mus musculus" = "mm",
  #                     "Danio rerio" = "dr")
  #       )
  #     ),
  #     mainPanel(
  #       width = 8,
        
  #     )
  #   )
  # )
)

# nolint end
# nolint start
library(pacman)
p_load(shinydashboard)

PROJECT <- "./"

ui <- navbarPage("ChIP sekoskaitos analizės", theme = shinytheme("cosmo"),
  includeCSS(paste0(PROJECT, "Scripts/styles.css")),
  # tabPanel("Formatų konversija",
  #   sidebarLayout(
  #     sidebarPanel(
  #       width = 4,
  #       fileInput(
  #         inputId = "provided_file",
  #         label = "Įkelkite ChIP sekoskaitos duomenų failą (-us)*:",
  #         multiple = TRUE,
  #         buttonLabel = "Ieškoti failo",
  #         placeholder = "Failas nepasirinktas"
  #       ),
  #       radioButtons(
  #         inputId = "file_format",
  #         label = "Nurodykite įvestų duomenų formatą:",
  #         choices = c("BED" = "bed", "BigWig" = "bw", "BigBed" = "bb")
  #       ),
  #     ),
  #     mainPanel(
  #       width = 8
  #     )
  #   )
  # ),
  tabPanel("Duomenų įkėlimas",
    sidebarLayout(
      sidebarPanel(
        width = 4,
        fileInput(
          inputId = "bigbed",
          label = "Įkelkite BED formato failą (-us)*:",
          multiple = TRUE,
          buttonLabel = "Ieškoti failo",
          placeholder = "Failas nepasirinktas"
        ),
        selectInput(
          inputId = "organism",
          label = "Nurodykite, iš kokio organizmo išgauti mėginiai:",
          choices = c(
            "Mus musculus", "Homo sapiens", "Rattus norvegicus", "Danio rerio",
            "Bos taurus", "Drosophila melanogaster", "Gallus gallus",
            "Macaca mulatta", "Pan troglodytes", "Sus scrofa", "Nenurodyta"
          ),
          selected = "Nenurodyta"
        ),
        fileInput(
          inputId = "pwm",
          label = "Įkelkite transkripcijos faktoriaus PWM matricą:",
          multiple = FALSE,
          buttonLabel = "Ieškoti failo",
          placeholder = "Failas nepasirinktas"
        ),
        # selectInput(
        #   inputId = "tf_options",
        #   label = "Pasirinkite transkripcijos faktorių:",
        #   choices = c(
        #     "Tbx5", "GATA3", "GATA4", "Tbx5", "Tcf21", "CTCF", "FOXA2",
        #     "Nenurodyta"),
        #   selected = "Nenurodyta"
        # ),
        textInput(
            inputId = "tf_options", 
            label = "Įveskite transkripcijos faktoriaus pavadinimą:",
            value = "Nenurodyta", 
            width = "100%",
            placeholder = "Transkripcijos faktoriaus pavadinimas..."),
        p("* - privalomas įvesties laukas", class = "info_text"),
        br(),
        br(),
        actionButton(
          inputId = "sample_data",
          label = "Pavyzdiniai duomenys",
          icon = icon("th"), 
          onclick = "window.open('https://vult-my.sharepoint.com/:f:/g/personal/daniele_stasiunaite_mif_stud_vu_lt/Eu4Rfi8-aqJBq-B_XcmeU_QBcdlQ7_Epb4yfy0VIIa5K_g?e=CFFJSG', '_blank')"
        )
      ),
      mainPanel(
        width = 8,
        tableOutput("table0"),
        shinydashboard::box(
            width = 12,
            withLoader(plotOutput("plot7", width = "40%", height = "250px"),
                       type = "html", loader = "dnaspin")
        ),
        # p("Mėginių informacija:"),
        DT::dataTableOutput("samples")
      )
    )
  ),
  tabPanel("Kokybės vertinimas",
    sidebarLayout(
      sidebarPanel(
        width = 4,
        p("Pasirinkite vieną arba kelis pateiktus mėginius, kurių duomenų
          kokybę vertinsite:", style = "font-weight:bold; font-size:17px;
          margin-left:0px"),
        DT::dataTableOutput("samples2")
      ),
      mainPanel(
        width = 8,
        tabsetPanel(
          tabPanel("Pikų skaičius mėginiuose", 
            p("Pateiktoje stulpelinėje diagramoje pavaizduota, kiek yra pikų
              kiekviename pateiktame BED formato faile:"),
            shinydashboard::box(
              width = 12, 
              withLoader(plotOutput("plot1"), type = "html",
                         loader = "dnaspin"),
              downloadButton("download1", "Atsisiųsti vaizdą")
            )
          ),
          tabPanel("Pikų skaičius chromosomose",
            p("Pateiktose stulpelinėse diagramose pavaizduota, kaip pikų
              skaičius pasiskirstęs skirtingose chromosomose:"),
            shinydashboard::box(
              width = 12,
              withLoader(plotOutput("plot2"), type = "html",
                         loader = "dnaspin"),
              downloadButton("download2", "Atsisiųsti vaizdą")
            )
          ),
          tabPanel("Mėginių panašumas",
            p("Pateiktame spalvų intensyvumo grafike pavaizduota, kokia pikų
              dalis (procentiškai) sutampa tarp skirtingų mėginių:"),
            shinydashboard::box(
              width = 12,
              withLoader(plotOutput("plot3"), type = "html",
                         loader = "dnaspin"),
              downloadButton("download3", "Atsisiųsti vaizdą")
            )
          ),
          tabPanel("Genominė distribucija",
            p("Pateiktame grafike pavaizduota pasirinktų mėginių genominių
              elementų procentinė dalis:"),
            shinydashboard::box(
              width = 12,
              withLoader(plotOutput("plot4"), type = "html",
                         loader = "dnaspin"),
              downloadButton("download4", "Atsisiųsti vaizdą")
            )
          ),
          tabPanel("Atstumas iki TSS",
            p("Pateiktame grafike pavaizduota pasirinktų mėginių anotuotų
              pikų atstumai iki TSS (angl. Transcription Start Site):"),
            shinydashboard::box(
              width = 12,
              withLoader(plotOutput("plot5"), type = "html",
                         loader = "dnaspin"),
              downloadButton("download5", "Atsisiųsti vaizdą")
            )
          ),
          tabPanel("Pikų profilis",
            p("Pateiktuose grafikuose pavaizduoti DNR nuskaitymų dažniai
               atitinkamose genominėse pozicijose:"),
            shinydashboard::box(
              width = 12,
              withLoader(plotOutput("plot6"), type = "html",
                         loader = "dnaspin")
            )
          )
        )
      )
    )
  ),
  tabPanel("Biologinės analizės",
    sidebarLayout(
      sidebarPanel(
        width = 4,
        p("Pasirinkite vieną arba kelis pateiktus mėginius, kuriems atliksite
          analizes:", style = "font-weight:bold; font-size:17px;
          margin-left:0px"),
        DT::dataTableOutput("samples3")
      ),
      mainPanel(
        width = 8,
        tabsetPanel(
          tabPanel("PWM matricos atitikimai",
            p("Pateiktoje stulpelinėje diagramoje pavaizduota, kokią procentinę
               dalį sudaro įkeltą transkripcijos faktoriaus pozicinę svorių
               matricą atitinkantys sekų fragmentai, palyginus su bendru pikų
               skaičiumi:"),
            shinydashboard::box(
              width = 12,
              withLoader(plotOutput("plot8"), type = "html",
                         loader = "dnaspin"),
              downloadButton("download8", "Atsisiųsti vaizdą")
            )
          ),
          tabPanel("GO analizė",
            tabsetPanel(
              id = "go",
              tabPanel("GO lentelė",
                id = "go_table",
                p("Biologinių procesų (BP) genų ontologijos rezultatai:"),
                width = 12,
                withLoader(reactableOutput("table1"), type = "html",
                           loader = "dnaspin"),
                p("Molekulinių funkcijų (MF) genų ontologijos rezultatai:"),
                width = 12,
                withLoader(reactableOutput("table2"), type = "html",
                           loader = "dnaspin"),
                p("Ląstelinių komponentų (CC) genų ontologijos rezultatai:"),
                width = 12,
                withLoader(reactableOutput("table3"), type = "html",
                           loader = "dnaspin")
              ),
              tabPanel("GO aciklinis grafas",
                id = "go_graph",
                p("Biologinių procesų (BP) kryptinis aciklinis grafas:"),
                width = 12,
                withLoader(plotOutput("plot9"), type = "html",
                           loader = "dnaspin"),
                p("Molekulinių funkcijų (MF) kryptinis aciklinis grafas:"),
                width = 12,
                withLoader(plotOutput("plot10"), type = "html",
                           loader = "dnaspin"),
                p("Ląstelinių komponentų (CC) kryptinis aciklinis grafas:"),
                width = 12,
                withLoader(plotOutput("plot11"), type = "html",
                           loader = "dnaspin")
              ),
              tabPanel("GO medžio struktūra",
                id = "go_tree",
                p("Biologinių procesų (BP) hierarchinis klasterizavimas:"),
                width = 12,
                withLoader(plotOutput("plot12"), type = "html",
                           loader = "dnaspin"),
                p("Molekulinių funkcijų (MF) hierarchinis klasterizavimas:"),
                width = 12,
                withLoader(plotOutput("plot13"), type = "html",
                           loader = "dnaspin"),
                p("Ląstelinių komponentų (CC) hierarchinis klasterizavimas:"),
                width = 12,
                withLoader(plotOutput("plot14"), type = "html",
                           loader = "dnaspin")
              )
            )
          ),
          tabPanel("Motyvų paieška De novo",
            p("Pateiktoje lentelėje pateikti identifikuoti pasirinkto mėginio
               motyvai, atlikus De novo motyvų paiešką."),
            p("Pastaba: priklausomai nuo pasirinkto mėginio dydžio De novo
               motyvų paieška gali trukti ilgiau nei 10 minučių.",
               style = "font-weight:bold; color:red"),
            shinydashboard::box(
              width = 12,
              withLoader(DT::dataTableOutput(outputId = "table4"),
                         type = "html", loader = "dnaspin"),
              downloadButton("downloadData", "Atsisiųsti lentelę")
            )
          )
        )
      )
    )
  ),
  tabPanel("Taikinių spėjimas",
    sidebarLayout(
      sidebarPanel(
        width = 4
      ),
      mainPanel(
        width = 8
      )
    )
  )
)
  

# nolint end
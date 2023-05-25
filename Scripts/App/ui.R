# nolint start
library(pacman)
p_load(shiny, shinythemes, shinydashboard, shinycustomloader, reactable)

PROJECT <- "./"

ui <- navbarPage(
  "ChIP sekoskaitos analizės",
  theme = shinytheme("cosmo"),
  includeCSS(paste0(PROJECT, "Scripts/styles.css")),
  tabPanel(
    "Duomenų įkėlimas",
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
        textInput(
          inputId = "tf_options", 
          label = "Įveskite transkripcijos faktoriaus pavadinimą:",
          value = "Nenurodyta", 
          width = "100%",
          placeholder = "Transkripcijos faktoriaus pavadinimas..."
        ),
        p("* - privalomas įvesties laukas", class = "info_text"),
        br(),
        br(),
        actionButton(
          inputId = "sample_data",
          label = "Pavyzdiniai duomenys",
          icon = icon("th"), 
          onclick = "window.open('https://vult-my.sharepoint.com/:f:/g/personal/daniele_stasiunaite_mif_stud_vu_lt/EtEGQ8POkapLhPv6eHvl48cB-jmes81M0JPW8PVWTz2QgA?e=wjSSKJ', '_blank')"
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
        DT::dataTableOutput("samples")
      )
    )
  ),
  tabPanel(
    "Kokybės vertinimas",
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
          tabPanel(
            "Pikų skaičius mėginiuose", 
            p("Pateiktoje stulpelinėje diagramoje pavaizduota, kiek yra pikų
              kiekviename pateiktame BED formato faile:"),
            shinydashboard::box(
              width = 12, 
              withLoader(
                plotOutput("plot1"), type = "html", loader = "dnaspin"
              ),
              downloadButton("download1", "Atsisiųsti vaizdą")
            )
          ),
          tabPanel(
            "Pikų skaičius chromosomose",
            p("Pateiktose stulpelinėse diagramose pavaizduota, kaip pikų
              skaičius pasiskirstęs skirtingose chromosomose:"),
            shinydashboard::box(
              width = 12,
              withLoader(
                plotOutput("plot2"), type = "html", loader = "dnaspin"
              ),
              downloadButton("download2", "Atsisiųsti vaizdą")
            )
          ),
          tabPanel(
            "Mėginių panašumas",
            p("Pateiktame spalvų intensyvumo grafike pavaizduota, kokia pikų
              dalis (procentiškai) sutampa tarp skirtingų mėginių:"),
            shinydashboard::box(
              width = 12,
              withLoader(
                plotOutput("plot3"), type = "html", loader = "dnaspin"
              ),
              downloadButton("download3", "Atsisiųsti vaizdą")
            )
          ),
          tabPanel(
            "Genominė distribucija",
            p("Pateiktame grafike pavaizduota pasirinktų mėginių genominių
              elementų procentinė dalis:"),
            shinydashboard::box(
              width = 12,
              withLoader(
                plotOutput("plot4"), type = "html", loader = "dnaspin"
              ),
              downloadButton("download4", "Atsisiųsti vaizdą")
            )
          ),
          tabPanel(
            "Atstumas iki TSS",
            p("Pateiktame grafike pavaizduota pasirinktų mėginių anotuotų
              pikų atstumai iki TSS (angl. Transcription Start Site):"),
            shinydashboard::box(
              width = 12,
              withLoader(
                plotOutput("plot5"), type = "html", loader = "dnaspin"
              ),
              downloadButton("download5", "Atsisiųsti vaizdą")
            )
          ),
          tabPanel(
            "Pikų profilis",
            p("Pateiktuose grafikuose pavaizduoti DNR nuskaitymų dažniai
               atitinkamose genominėse pozicijose:"),
            shinydashboard::box(
              width = 12,
              withLoader(
                plotOutput("plot6"), type = "html", loader = "dnaspin"
              ),
              downloadButton("download6", "Atsisiųsti vaizdą")
            )
          )
        )
      )
    )
  ),
  tabPanel(
    "Biologinės analizės",
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
          tabPanel(
            "PWM matricos atitikimai",
            p("Pateiktoje stulpelinėje diagramoje pavaizduota, kokią procentinę
               dalį sudaro įkeltą transkripcijos faktoriaus pozicinę svorių
               matricą atitinkantys sekų fragmentai, palyginus su bendru pikų
               skaičiumi:"),
            shinydashboard::box(
              width = 12,
              withLoader(
                plotOutput("plot8"), type = "html", loader = "dnaspin"
              ),
              downloadButton("download8", "Atsisiųsti vaizdą")
            )
          ),
          tabPanel(
            "GO analizė",
            tabsetPanel(
              id = "go",
              tabPanel(
                "GO lentelė",
                id = "go_table",
                p("Biologinių procesų (BP) genų ontologijos rezultatai:"),
                width = 12,
                withLoader(
                  reactableOutput("table1"), type = "html", loader = "dnaspin"
                ),
                p("Molekulinių funkcijų (MF) genų ontologijos rezultatai:"),
                width = 12,
                withLoader(
                  reactableOutput("table2"), type = "html", loader = "dnaspin"
                ),
                p("Ląstelinių komponentų (CC) genų ontologijos rezultatai:"),
                width = 12,
                withLoader(
                  reactableOutput("table3"), type = "html", loader = "dnaspin"
                )
              ),
              tabPanel(
                "GO aciklinis grafas",
                id = "go_graph",
                p("Biologinių procesų (BP) kryptinis aciklinis grafas:"),
                width = 12,
                withLoader(
                  plotOutput("plot9"), type = "html", loader = "dnaspin"
                ),
                p("Molekulinių funkcijų (MF) kryptinis aciklinis grafas:"),
                width = 12,
                withLoader(
                  plotOutput("plot10"), type = "html", loader = "dnaspin"
                ),
                p("Ląstelinių komponentų (CC) kryptinis aciklinis grafas:"),
                width = 12,
                withLoader(
                  plotOutput("plot11"), type = "html", loader = "dnaspin"
                )
              ),
              tabPanel(
                "GO medžio struktūra",
                id = "go_tree",
                p("Biologinių procesų (BP) hierarchinis klasterizavimas:"),
                width = 12,
                withLoader(
                  plotOutput("plot12"), type = "html", loader = "dnaspin"
                ),
                p("Molekulinių funkcijų (MF) hierarchinis klasterizavimas:"),
                width = 12,
                withLoader(
                  plotOutput("plot13"), type = "html", loader = "dnaspin"
                ),
                p("Ląstelinių komponentų (CC) hierarchinis klasterizavimas:"),
                width = 12,
                withLoader(
                  plotOutput("plot14"), type = "html", loader = "dnaspin"
                )
              )
            )
          ),
          tabPanel(
            "Motyvų paieška de novo",
            p("Pateiktoje lentelėje pateikti identifikuoti pasirinkto mėginio
               motyvai, atlikus De novo motyvų paiešką."),
            p("Pastaba: priklausomai nuo pasirinkto mėginio dydžio De novo
               motyvų paieška gali trukti ilgiau nei 10 minučių.",
               style = "font-weight:bold; color:red"),
            shinydashboard::box(
              width = 12,
              withLoader(
                DT::dataTableOutput(outputId = "table4"),
                type = "html",
                loader = "dnaspin"
              )
            )
          )
        )
      )
    )
  ),
  tabPanel(
    "Taikinių spėjimas",
    sidebarLayout(
      sidebarPanel(
        width = 4,
        p("Pasirinkite vieną arba kelis pateiktus mėginius, kuriuos naudosite
           transkripcijos faktoriaus taikinių spėjimui:",
           style = "font-weight:bold; font-size:17px; margin-left:0px"),
        DT::dataTableOutput("samples4"),
        br(),
        selectInput(
          inputId = "organism_predict",
          label = "Nurodykite, kokiame organizme spėsite transkripcijos
                  faktorių taikinius:",
          choices = c(
            "Mus musculus", "Homo sapiens", "Rattus norvegicus", "Danio rerio",
            "Bos taurus", "Drosophila melanogaster", "Gallus gallus",
            "Macaca mulatta", "Pan troglodytes", "Sus scrofa", "Nenurodyta"
          ),
          selected = "Nenurodyta"
        ),
        br(),
        selectInput(
          inputId = "min_score",
          label = "Nurodykite minimalų transkripcijos faktoriaus PWM matricos
                  ir sekų fragmentų atitikimo procentą:",
          choices = c("70%", "75%", "80%", "85%", "90%", "95%"),
          selected = "70%"
        ),
        br(),
        br(),
        p("Pastaba: priklausomai nuo pasirinktų mėginių dydžių metodo
          realizacija gali užtrukti ilgiau nei 20 minučių.",
          style = "font-weight:bold; color:red"),
      ),
      mainPanel(
        width = 8,
        tabsetPanel(
          tabPanel(
            "PWM atitikimai", 
            p("Pateiktoje sklaidos diagramoje pavaizduota, kiek kiekviename
              gene nustatyta pozicinės svorių matricos
              atitikimų:"),
            # DT::dataTableOutput(outputId = "table5")
            shinydashboard::box(
              width = 12,
              withLoader(
                plotOutput(outputId = "plot15"),
                type = "html",
                loader = "dnaspin"
              ),
              downloadButton("download15", "Atsisiųsti vaizdą")
            )
          ),
          tabPanel(
            "PWM atitikimų skirtumai mėginiuose",
            p("Pateiktoje sudėtinėje stulpelinėje diagramoje pavaizduota,
              kokią procentinę dalį sudaro genai, kuriuose nustatytas didesnis
              PWM matricos atitikimų skaičius užklausos sekose:"),
            shinydashboard::box(
              width = 12,
              withLoader(
                plotOutput(outputId = "plot16"),
                type = "html",
                loader = "dnaspin"
              )
            )
          ),
          tabPanel(
            "Rezultatų atsisiuntimas", 
            p("Atsisiųsti išgautas anotuotų pikų genų aminorūgščių sekas:"),
            downloadButton("download18", "Atsisiųsti failą"),
            p("Atsisiųsti Blastp rezultatų lentelę:"),
            downloadButton("download19", "Atsisiųsti failą")
          )
        )
      )
    )
  )
)

# nolint end
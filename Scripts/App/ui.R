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
  includeCSS(paste0(PROJECT, "Scripts/styles.css")),
  # GALBŪT ČIA REIKIA NAUJO TAB'O, KURIAME REIKĖTŲ ĮKELTI VISUS DUOMENIS?
  tabPanel("Duomenų įkėlimas",
    sidebarLayout(
      sidebarPanel(
        width = 6,
        fileInput(
          inputId = "bigbed",
          label = "Įkelkite BED formato failą (-us)*:",
          multiple = TRUE,
          buttonLabel = "Ieškoti failo",
          placeholder = "Failas nepasirinktas"
        ),
        fileInput(
          inputId = "pwm",
          label = "Įkelkite transkripcijos faktoriaus PWM matricą:",
          multiple = FALSE,
          buttonLabel = "Ieškoti failo",
          placeholder = "Failas nepasirinktas"
        ),
        p("* - privalomas įvesties laukas", class = "info_text")
      ),
      mainPanel()
    )
  ),
  tabPanel("Duomenų kokybė",
    sidebarLayout(
      sidebarPanel(
        width = 4,
        p("Pasirinkite vieną arba kelis pateiktus mėginius, kurių duomenų
          kokybę vertinsite:", style = "font-weight:bold; font-size:17px;
          margin-left:0px"),
        DT::dataTableOutput("samples"),
        p("Nurodykite, ar norite naudoti apjungtų mėginių duomenis:",
          style = "font-weight:bold; font-size:17px; margin-left:0px; padding-top: 80px"),
        selectInput(
          inputId = "sample_overlap",
          label = "",
          choices = c(
            "Netaikyti mėginių apjungimo" = "no_join",
            "Apjungti pažymėtus mėginius" = "marked_samples",
            "Apjungti visus įkeltus mėginius" = "all_samples",
            "Rasti bendrus regionus pažymėtuose mėginiuose" = "common_marked",
            "Rasti bendrus regious visuose mėginiuose" = "common_all")
        )
      ),
      mainPanel(
        width = 8,
        tabsetPanel(
          tabPanel("Pikų skaičius mėginiuose", 
            p("Pateiktoje stulpelinėje diagramoje pavaizduota, kiek yra pikų
              kiekviename pateiktame BED formato faile:"),
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
            p("Pateiktame grafike pavaizduota pasirinktų mėginių genominių
              elementų procentinė dalis:"),
            box(
              width = 12,
              withLoader(plotOutput("plot4"), type = "html", loader = "dnaspin")
            )
          ),
          tabPanel("Atstumas iki TSS",
            p("Pateiktame grafike pavaizduota pasirinktų mėginių anotuotų
              pikų atstumai iki TSS (angl. Transcription Start Site):"),
            box(
              width = 12,
              withLoader(plotOutput("plot5"), type = "html", loader = "dnaspin")
            )
          ),
          tabPanel("Pikų profilis",
            p("Pateiktuose grafikuose pavaizduoti DNR nuskaitymų dažniai
               atitinkamose genominėse pozicijose:"),
            box(
              width = 12,
              withLoader(plotOutput("plot6"), type = "html", loader = "dnaspin")
            )
          )
        )
      )
    )
  ),
  tabPanel("Analizės",
    sidebarLayout(
      sidebarPanel(
        width = 4,
        selectInput(
          inputId = "tf_options",
          label = "Pasirinkite transkripcijos faktorių:",
          choices = c("Tbx5", "GATA3", "Tcf21"),
          selected = "GATA3"
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
          tabPanel("Transkripcijos faktoriaus motyvas", 
            p("Pateiktame paveiksle pavaizduotas įkeltą pozicinę svorių matricą
               atitinkančio transkripcijos faktoriaus motyvo sekos logotipas:"),
            box(
              width = 12, 
              withLoader(plotOutput("plot7", width = "70%"), type = "html",
                         loader = "dnaspin")
            )
          ),
          tabPanel("PWM matricos atitikimai",
            p("Pateiktoje stulpelinėje diagramoje pavaizduota, kokią procentinę
               dalį sudaro įkeltą transkripcijos faktoriaus pozicinę svorių
               matricą atitinkantys sekų fragmentai, palyginus su bendru pikų
               skaičiumi:"),
            box(
              width = 12,
              withLoader(plotOutput("plot8"), type = "html", loader = "dnaspin")
            )
          ),
          tabPanel("GO analizė",
            p("Description goes here..."),
            box(
              width = 12,
              withLoader(plotOutput("plot66"), type = "html", loader = "dnaspin")
            )
          ),
          tabPanel("Motyvų paieška De novo",
            p("Pateiktoje lentelėje pateikti identifikuoti pasirinkto mėginio
               motyvai, atlikus De novo motyvų paiešką."),
            p("Pastaba: priklausomai nuo pasirinkto mėginio dydžio De novo
               motyvų paieška gali trukti ilgiau nei 10 minučių.",
               style = "font-weight:bold; color:red"),
            box(
              width = 12,
              withLoader(DT::dataTableOutput(outputId = "table1"),
                         type = "html", loader = "dnaspin"),
              downloadButton("downloadData", "Download")
            )
          )
        )
      )
    )
  ),
  tabPanel("Taikinių spėjimas",
    sidebarLayout(
      sidebarPanel(
        width = 4,
        radioButtons(
          inputId = "genome",
          label = "Pasirinkite genomą:",
          choices = c("Homo sapiens" = "hg", "Mus musculus" = "mm",
                      "Danio rerio" = "dr")
        )
      ),
      mainPanel(
        width = 8
        
      )
    )
  )
)

# nolint end
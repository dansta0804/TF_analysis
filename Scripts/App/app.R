# nolint start
library(pacman)
p_load(shiny, data.table, rtracklayer, ggplot2, ggthemes, plyranges, ggpubr,
       BRGenomics, reshape2, plotly, heatmaply, dplyr, gplots, genomation,
       Biostrings, scales, GenomicRanges, DT, shinythemes, shinycustomloader,
       TxDb.Mmusculus.UCSC.mm10.knownGene)

PROJECT <- "/home/daniele/Desktop/IV_course/II_semester/TF_analysis/"
# setwd("/home/daniele/Desktop/IV_course/II_semester/TF_analysis/Scripts/App/") 
# source(paste0("/home/daniele/Desktop/IV_course/II_semester/TF_analysis/Scripts/App/ui.R"))
# source(paste0(PROJECT, "Scripts/App/server.R"))

options(scipen = 100)
options(shiny.maxRequestSize = 300 * 1024 ^ 2)

FUNCTIONS   <- paste0(PROJECT, "Scripts/functions.R")
INPUTS      <- paste0(PROJECT, "Analyses/Tbx5_analysis_I/Inputs/BED/")
INPUT       <- paste0(PROJECT, "Input/")
RESULTS     <- paste0(PROJECT, "Results/")
source(FUNCTIONS)

server <- function(input, output, session) {
  # Klaida yra tame, kad vienodi pavadinimai 'output$samples', todėl
  # išsiaiškinti, kaip padaryti, kad būtų galima naudoti tuos pačius
  # pavadinimus.
  observe({
    req(input$file_format)
    if (input$file_format == "bed") {
      output$samples <-
        DT::renderDataTable({input$bigbed[, c("name", "size")]}, server = TRUE)
    } else if (input$file_format == "bw") {
      output$samples <- DT::renderDataTable({
        converted_table <- data.frame(matrix(ncol = 3, nrow = 0))
        colnames(converted_table) <-
          c("Old file name", "New file name", "Graph name")
        
        for(sample in 1:length(input$bigbed[, c("name")])) {
          names <- input$bigbed[sample, 'name']
          system(paste("/home/daniele/Tools/bigWigToBedGraph",
                 paste0(INPUT, "Mus_musculus/", names),
                 paste0(INPUT, "Mus_musculus/", names, ".bedGraph")))
          row <- c(names, paste0(names, ".bedGraph"),
                   paste0("Sample", sample, "_", input$organism))

          converted_table[nrow(converted_table) + 1, ] <- row
        }
        converted_table
      })
    }
  })

  observe({
    req(input$samples_rows_selected)
    selRow <- input$bigbed[input$samples_rows_selected,]
    names <- selRow[[1]]
    samples <- selRow[[4]]

    
    # if (input$sample_overlap == "no_join") {
    #     selRow <- input$bigbed[input$samples_rows_selected, ]
    #     names <- selRow[[1]]
    #     samples <- selRow[[4]]
    # }

    ########################### GENOMIC DATA QUALITY #######################
    # GENOMIC DATA QUALITY - PEAK COUNT (PLOT 1):
    output$plot1 <- renderPlot({
      bigbed_files <- list()

      


      # for(sample in 1:length(samples)) {
      #   names <- input$bigbed[[sample, 'name']]
      #   system(paste("/home/daniele/Tools/bigWigToBedGraph",
      #                paste0(INPUTS, "Homo_sapiens/", names),
      #                paste0(INPUTS, "Homo_sapiens/", names, ".bedGraph")))
      #   bigbed_files[[sample]] <-
      #     read.table(file = paste0(INPUTS, "Homo_sapiens/", names, ".bedGraph"))
      #       #  names(bigbed_files)[sample] <-
      #       #  substring(names[sample], 1, 11)
      # }




      for(sample in 1:length(samples)) {
        bigbed_files[[sample]] <- read.table(file = samples[sample])
        names(bigbed_files)[sample] <- substring(names[sample], 1, 11)
      }

      if (length(input$samples_rows_selected) == 0) {
        return()
      } else {
        peaks <- data.frame(matrix(ncol = 2, nrow = 0))
        colnames(peaks) <- c("Experiment", "Peak_count")

        for (file in 1:length(bigbed_files)) {
          nrows <- length(rownames(bigbed_files[[file]]))
          row <- c(names(bigbed_files[file]), nrows)
          peaks[nrow(peaks) + 1, ] <- row
        }

        ggplot(peaks, aes(x = Experiment, y = as.numeric(Peak_count))) +
          geom_bar(stat = "identity", position = "dodge", width = 0.5,
                   color = "black", fill = "#930d1f") +
          labs(x = "", y = "Pikų skaičius", size = 5) +
          ylim(0, 200000) +
          scale_y_continuous(labels = label_number(suffix = " K", scale = 1e-3)) +
          theme(panel.background = element_rect(fill = "#eeeef1",
                                                colour = "#4c0001"),
                panel.grid.major.y = element_line(colour = "#cab5b5",
                                                  size = 0.3,
                                                  linetype = "dashed"),
                panel.grid.minor.y = element_line(colour = "#cab5b5",
                                                  size = 0.3,
                                                  linetype = "dashed"),
                panel.grid.major.x = element_line(colour = "#cab5b5",
                                                  size = 0.2,
                                                  linetype = "longdash"),
                panel.grid.minor.x = element_line(colour = "#cab5b5",
                                                  size = 0.2,
                                                  linetype = "longdash"),
                axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                                           size = 11, face = "bold",
                                           color = "black"),
                axis.text.y = element_text(size = 11, face = "bold",
                                           color = "black"),
                axis.title.x = element_text(size = 2),
                axis.title.y = element_text(size = 16),
                plot.title = element_text(hjust = 0.5, face = "bold"))
      }
    })

    # GENOMIC DATA QUALITY - PEAK COUNT DISTRIBUTION BY CHROMOSOME (PLOT 2):
    output$plot2 <- renderPlot({
      bigbed_files <- list()
      grl <- GRangesList()

      for(sample in 1:length(samples)) {
        bigbed_files[[sample]] <- read.table(file = samples[sample])

        colnames(bigbed_files[[sample]]) <-
          c("chrom", "start", "end", "name", rep(paste0("Other", 1:7)))
        grl[[sample]] <- makeGRangesFromDataFrame(bigbed_files[[sample]],
                                                  keep.extra.columns = TRUE)
        names(grl)[sample] <- substring(names[sample], 1, 11)
      }

      if (length(input$samples_rows_selected) == 0) {
        return()
      } else {
        peak_counts <-
          lapply(names(grl), count_peaks, objects = grl) %>%
          bind_rows()

          # Calling factor() function in order to maintain certain Chromosome
          # and Name order:
          unique_chr <- unique(peak_counts$Chromosome)
          peak_counts$Chromosome <-
            factor(peak_counts$Chromosome, levels = c(unique_chr))
          peak_counts$Name <-
            factor(peak_counts$Name, levels = unique(peak_counts$Name))

          # Creating barplots that visualize peak differences between
          # different chromosomes:
          ggplot(peak_counts, aes(x = Name, y = as.numeric(Peak_count))) +
            geom_bar(stat = "identity", color = "black", fill = "#930d1f") +
            ylab("Pikų skaičius") +
            facet_wrap(~ Chromosome, ncol = 7) +
            xlab("") +
            scale_y_continuous(labels = label_number(suffix = " K",
                                                     scale = 1e-3)) +
            theme_linedraw() +
            theme(axis.text.x = element_text(angle = 90, size = 12,
                                             vjust = 0.5),
                  legend.position = "none",
                  axis.text.y = element_text(size = 14, face = "bold"),
                  axis.title.x = element_text(size = 14, colour = "black"),
                  axis.title.y = element_text(size = 20, colour = "black"),
                  strip.background = element_rect(fill = "white"),
                  strip.text = element_text(colour = "black", face = "bold",
                                            size = 16))
        }
    })

    # GENOMIC DATA QUALITY - OVERLAPS BETWEEN SAMPLES (JACCARD) (PLOT 3):
    output$plot3 <- renderPlot({
      bigbed_files <- list()

      for(sample in 1:length(samples)) {
        bigbed_files[[sample]] <- read.table(file = samples[sample])
        names(bigbed_files)[sample] <- substring(names[sample], 1, 11)
        colnames(bigbed_files[[sample]]) <-
              c("chrom", "start", "end", "name", rep(paste0("Other", 1:7)))
        bigbed_files[[sample]] <-
          makeGRangesFromDataFrame(bigbed_files[[sample]],
                                   keep.extra.columns = TRUE)
      }

      if (length(input$samples_rows_selected) == 0) {
        return()
      } else {
        coef_matrix <- matrix(nrow = length(bigbed_files),
                              ncol = length(bigbed_files))

        # Calculating Jaccard coefficient for sample pair:
        for (i in 1:length(bigbed_files)) {
          for (y in 1:length(bigbed_files)) {
            coef_matrix[i, y] <- jaccard(bigbed_files, i, y)
          }
        }

        # Setting colnames and rownames for the matrix:
        colnames(coef_matrix) <- names(bigbed_files)
        rownames(coef_matrix) <- names(bigbed_files)

        coef_mat1 <- coef_matrix
        coef_mat2 <- coef_matrix

        # Passing Jaccard coefficients to matrix except for the diagonal -
        # it contains 'NA':
        coef_mat1[lower.tri(coef_mat1, diag = TRUE)] <- NA
        coef_mat2[upper.tri(coef_mat2, diag = TRUE)] <- NA

        # Binding two matrixes:
        coef_mat_bind <- rbind(data.matrix(coef_mat1), data.matrix(coef_mat2))

        # Translating matrix to dataframe using melt() function:
        melt_coef_mat <- melt(coef_mat_bind, na.rm = TRUE)

        # Creating a heatmap that shows similarity between samples:
        ggplot(melt_coef_mat, aes(x = Var2, y = Var1, fill = value)) +
          geom_tile(color = "black") +
          geom_text(aes(label = round(value, digits = 3)), size = 4.5,
                        color = "#030101", fontface = "bold") +
          labs(x = "", y = "") +
          scale_fill_gradient(low = "#ffee8e", high = "#ab1f1f") +
          guides(fill = guide_colourbar(title = "Koeficientas",
                                        face = "bold")) +
          theme(axis.text = element_text(size = 12, colour = "black",
                                         face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title.x = element_text(size = 14, colour = "black"),
          axis.title.y = element_text(size = 14, colour = "black"),
          panel.grid.major = element_line(color = "#eeeeee"),
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          legend.position = "bottom")
      }
    })

    # GENOMIC DATA QUALITY - GENOMIC DISTRIBUTION (PLOT 4):
    output$plot4 <- renderPlot({
      grl <- GRangesList()
      bigbed_files <- list()


      for(sample in 1:length(samples)) {
        bigbed_files[[sample]] <- read.table(file = samples[sample])

        colnames(bigbed_files[[sample]]) <-
          c("chrom", "start", "end", "name", rep(paste0("Other", 1:7)))
        grl[[sample]] <-
          makeGRangesFromDataFrame(bigbed_files[[sample]],
                                    keep.extra.columns = TRUE)
        names(grl)[sample] <- substring(names[sample], 1, 11)
      }

      if (length(input$samples_rows_selected) == 0) {
        return()
      } else {
        mm_known_genes <- TxDb.Mmusculus.UCSC.mm10.knownGene
        peak_annotations <- lapply(grl, annotatePeak, TxDb = mm_known_genes,
                                   tssRegion = c(-3000, 3000), verbose = FALSE)

        plotAnnoBar(peak_annotations, ylab = "Procentinė dalis (%)", title = "")
      }
    })

    # GENOMIC DATA QUALITY - DISTANCE TO TSS (PLOT 5):
    output$plot5 <- renderPlot({
      grl <- GRangesList()
      bigbed_files <- list()

      for(sample in 1:length(samples)) {
        bigbed_files[[sample]] <- read.table(file = samples[sample])

        colnames(bigbed_files[[sample]]) <-
          c("chrom", "start", "end", "name", rep(paste0("Other", 1:7)))
        grl[[sample]] <-
          makeGRangesFromDataFrame(bigbed_files[[sample]],
                                   keep.extra.columns = TRUE)
        names(grl)[sample] <- substring(names[sample], 1, 11)
      }

      if (length(input$samples_rows_selected) == 0) {
        return()
      } else {
        mm_known_genes <- TxDb.Mmusculus.UCSC.mm10.knownGene
        peak_annotations <- lapply(grl, annotatePeak, TxDb = mm_known_genes,
                                   tssRegion = c(-3000, 3000), verbose = FALSE)
        plotDistToTSS(peak_annotations, ylab = "Atstumas", title = "")
      }
    })

    # GENOMIC DATA QUALITY - PEAK PROFILE (PLOT 6):
    output$plot6 <- renderPlot({
      grl <- GRangesList()
      plots <- list()
      bigbed_files <- list()

      for(sample in 1:length(samples)) {
        bigbed_files[[sample]] <- read.table(file = samples[sample])

        colnames(bigbed_files[[sample]]) <-
          c("chrom", "start", "end", "name", rep(paste0("Other", 1:7)))
        grl[[sample]] <-
          makeGRangesFromDataFrame(bigbed_files[[sample]],
                                   keep.extra.columns = TRUE)
        names(grl)[sample] <- substring(names[sample], 1, 11)
      }

      if (length(input$samples_rows_selected) == 0) {
        return()
      } else {
        mm_known_genes <- TxDb.Mmusculus.UCSC.mm10.knownGene
        promoter <- getPromoters(TxDb = mm_known_genes, upstream = 3000,
                                 downstream = 3000)
        
        for (peak_file in 1:length(grl)) {
          remove_y <- theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank()
          )

          peak <- makeGRangesFromDataFrame(grl[[peak_file]],
                                           keep.extra.columns = TRUE)
          
          tagMatrix <- getTagMatrix(peak, windows = promoter)
          
          plots[[peak_file]] <- plotAvgProf(tagMatrix, xlim = c(-3000, 3000),
                                            xlab = "Genomic Region (5'->3')",
                                            ylab = "Read Count Frequency") +
          remove_y
        }
        ggarrange(plotlist = plots, nrow = 1, common.legend = TRUE,
                  widths = 5, legend = "bottom")
      }
    })
  })

  output$samples2 <- DT::renderDataTable({input$bigbed[, c("name", "size")]},
                                         server = TRUE)

  ########################## GENOMIC DATA ANALYSIS #######################
  observe({       
    # GENOMIC DATA ANALYSIS - TF MOTIF (PLOT 7):
    output$plot7 <- renderPlot({
      req(input$pwm)
      mpwm <- read.table(file = input$pwm$datapath)
      mpwm <- t(mpwm)

      # Setting matrix rownames:
      rownames(mpwm) <- c("A", "C", "G", "T")
      ggseqlogo(mpwm)
    })
                     
    req(input$samples2_rows_selected)
    selRow <- input$bigbed[input$samples2_rows_selected, ]
    names <- selRow[[1]]
    samples2 <- selRow[[4]]

    # GENOMIC DATA ANALYSIS - PWM MATRIX MATCHES (PLOT 8):
    output$plot8 <- renderPlot({
      req(input$bigbed)
      req(input$pwm)
      bigbed_files <- list()

      for(sample in 1:length(samples2)) {
        bigbed_files[[sample]] <-
          read.table(file = input$bigbed[[sample, 'datapath']])
        
        colnames(bigbed_files[[sample]]) <-
          c("chrom", "start", "end", "name", rep(paste0("Other", 1:7)))
        bigbed_files[[sample]] <-
          makeGRangesFromDataFrame(bigbed_files[[sample]],
                                   keep.extra.columns = TRUE)
        names(bigbed_files)[sample] <- input$bigbed[[sample, 'name']]
      }

      if (length(input$samples2_rows_selected) == 0) {
        return()
      } else {
        # Čia turės būti naudojama nebūtinai šita matrica. Reikia pateikti visų
        # galimų TF PWM matricų sąrašą?
        mpwm <- read.table(file = input$pwm$datapath)
        mpwm <- t(mpwm)
        rownames(mpwm) <- c("A", "C", "G", "T")

        peak_sequences <- list()

        for (file in 1:length(bigbed_files)) {
          peak_sequences[[file]] <-
            getSeq(BSgenome.Mmusculus.UCSC.mm10, bigbed_files[[file]])
          names(peak_sequences)[file] <- names(bigbed_files)[file]
        }

        tbx5_motifs <- data.frame(matrix(nrow = 0, ncol = 4))
        colnames(tbx5_motifs) <-
          c("Sample", "Motif_count", "Peak_count", "Percentage")

        for (i in 1:length(peak_sequences)) {
          filename <- input$bigbed[[i, 'datapath']]
          name <- substring(names(peak_sequences[i]), 1, 11)

          motif <- find_motif_hits(peak_sequences[[i]], mpwm)
          peaks <- calculate_peaks(filename)
          percentage <- round((motif / peaks) * 100, 2)

          data_row <- c(name, motif, peaks, paste0(percentage, "%"))
          tbx5_motifs[nrow(tbx5_motifs) + 1, ] <- data_row
        }

        # Reading a file that stores information about Tbx5 motif counts and
        # peak percentages:
        motif_data <- tbx5_motifs

        # Calling factor() function in order to maintain certain Sample
        # order:
        motif_data$Sample <-
                factor(motif_data$Sample, levels = motif_data$Sample)

        # Subseting data to extract all columns except for 'X' column
        # (column 1):
        subset_df <- motif_data[, 1:3]

        # 'Melting' the dataframe:
        melted_df <- melt(subset_df, id = c("Sample"))

        ggplot(data = melted_df, aes(x = Sample, y = as.numeric(value),
                                     fill = variable, label = value)) +
          geom_bar(stat = "identity", colour = "#35170450", size = 0.5,
                   width = 0.8) +
          scale_fill_manual(values = c("#e3a15e", "#c7633b"),
                            labels = c("Tbx5 motyvų skaičius",
                                       "Pikų skaičius")) +
          scale_y_continuous(labels = label_number(suffix = " K",
                                                   scale = 1e-3)) +
          guides(fill = guide_legend(title = "Spalvų paaiškinimas", size = 6)) +
          labs(title = "", x = "", y = "TF/Pikų skaičius") +
          theme(axis.text = element_text(size = 10, colour = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                                       size = 12, face = "bold"),
            axis.text.y = element_text(size = 12, face = "bold"),
            axis.title.x = element_text(size = 14, colour = "black"),
            axis.title.y = element_text(size = 14, colour = "black"),
            panel.grid.major = element_line(color = "#eeeeee"),
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            panel.background = element_rect(fill = "#eeeef1",
                                            colour = "#4c0001"),
            panel.grid.major.y = element_line(colour = "#cab5b5",
                                              size = 0.3, linetype = "dashed"),
            panel.grid.minor.y = element_line(colour = "#cab5b5",
                                              size = 0.3, linetype = "dashed"),
            panel.grid.major.x = element_line(colour = "#cab5b5",
                                              size = 0.2, linetype = "longdash"),
            panel.grid.minor.x = element_line(colour = "#cab5b5",
                                              size = 0.2, linetype = "longdash"),
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 11))
      }
    })

    # GENOMIC DATA ANALYSIS - DE NOVO MOTIFS (TABLE 1):
    output$table1 <- DT::renderDataTable({
      req(input$bigbed)
      bigbed_files <- list()
      names <- c()

      for(sample in 1:length(input$bigbed[, 1])) {
        names <- c(names, input$bigbed[[sample, 'name']])
        bigbed_files[[sample]] <-
          read.table(file = input$bigbed[[sample, 'datapath']])
        names(bigbed_files)[sample] <-
          substring(input$bigbed[[sample, 'name']], 1, 11)
      }

      system(paste("findMotifsGenome.pl", paste0(INPUTS, names), "mm10",
             paste0(PROJECT, RESULTS, names, "/"), "-size 200 -mask"),
             intern = FALSE)

      homer_motifs <-
        read.table(paste0(PROJECT, RESULTS, names, "/knownResults.txt"),
                   sep = "\t", header = FALSE, skip = 1)

      colnames(homer_motifs) <-
        c("Motif Name", "Consensus", "P value", "Log P value",
          "q value (Benjamini)", "Target Sequences with Motif (#)",
          "Target Sequences with Motif (%)",
          "Background Sequences with Motif (#)",
          "Background Sequences with Motif (%)")

      for (motif in 1:nrow(homer_motifs)) {
        homer_motifs$`Motif Name`[motif] <-
            strsplit(homer_motifs$`Motif Name`, "/")[[motif]][1]
      }
      homer_motifs
    })
  })
}

ui <- navbarPage("ChIP sekoskaitos analizės", theme = shinytheme("cosmo"),
  includeCSS(paste0(PROJECT, "Scripts/styles.css")),
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
        fileInput(
          inputId = "pwm",
          label = "Įkelkite transkripcijos faktoriaus PWM matricą:",
          multiple = FALSE,
          buttonLabel = "Ieškoti failo",
          placeholder = "Failas nepasirinktas"
        ),
        radioButtons(
          inputId = "file_format",
          label = "Nurodykite įvestų duomenų formatą:",
          choices = c("BED" = "bed", "BigWig" = "bw", "BigBed" = "bb")
        ),
        selectInput(
          inputId = "organism",
          label = "Nurodykite, iš kokio organizmo išgauti mėginiai:",
          choices = c(
            "Mus musculus" = "mm", "Homo sapiens" = "hs", "Danio rerio" = "dr"
          )
        ),
        p("* - privalomas įvesties laukas", class = "info_text")
      ),
      mainPanel(
        width = 8,
        p("Esančioje lentelėje pateikta įkeltus mėginius apibendrinanti
          informacija: originalus mėginio pavadinimas, failo tipas bei
          naujas mėginiui priskirtas pavadinimas, kuris bus naudojamas
          duomenų kokybės bei analizių grafikuose."),
        DT::dataTableOutput("samples")
      )
    )
  ),
  tabPanel("Duomenų kokybė",
    sidebarLayout(
      sidebarPanel(
        width = 4,
        p("Pasirinkite vieną arba kelis pateiktus mėginius, kurių duomenų
          kokybę vertinsite:", style = "font-weight:bold; font-size:17px;
          margin-left:0px"),
        # DT::dataTableOutput("samples"),
        p("Nurodykite, ar norite naudoti apjungtų mėginių duomenis:",
          style = "font-weight:bold; font-size:17px; margin-left:0px;
                   padding-top: 80px"),
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
            shinydashboard::box(
              width = 12, 
              withLoader(plotOutput("plot1"), type = "html", loader = "dnaspin")
            )
          ),
          tabPanel("Pikų skaičius chromosomose",
            p("Pateiktose stulpelinėse diagramose pavaizduota, kaip pikų skaičius
              pasiskirstęs skirtingose chromosomose:"),
            shinydashboard::box(
              width = 12,
              withLoader(plotOutput("plot2"), type = "html", loader = "dnaspin")
            )
          ),
          tabPanel("Mėginių panašumas",
            p("Pateiktame spalvų intensyvumo grafike pavaizduota, kokia pikų
              dalis (procentiškai) sutampa tarp skirtingų mėginių:"),
            shinydashboard::box(
              width = 12,
              withLoader(plotOutput("plot3"), type = "html", loader = "dnaspin")
            )
          ),
          tabPanel("Genominė distribucija",
            p("Pateiktame grafike pavaizduota pasirinktų mėginių genominių
              elementų procentinė dalis:"),
            shinydashboard::box(
              width = 12,
              withLoader(plotOutput("plot4"), type = "html", loader = "dnaspin")
            )
          ),
          tabPanel("Atstumas iki TSS",
            p("Pateiktame grafike pavaizduota pasirinktų mėginių anotuotų
              pikų atstumai iki TSS (angl. Transcription Start Site):"),
            shinydashboard::box(
              width = 12,
              withLoader(plotOutput("plot5"), type = "html", loader = "dnaspin")
            )
          ),
          tabPanel("Pikų profilis",
            p("Pateiktuose grafikuose pavaizduoti DNR nuskaitymų dažniai
               atitinkamose genominėse pozicijose:"),
            shinydashboard::box(
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
        p("Pasirinkite vieną arba kelis pateiktus mėginius, kuriems atliksite
          analizes:", style = "font-weight:bold; font-size:17px;
          margin-left:0px"),
        DT::dataTableOutput("samples2"),
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
            shinydashboard::box(
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
            shinydashboard::box(
              width = 12,
              withLoader(plotOutput("plot8"), type = "html", loader = "dnaspin")
            )
          ),
          tabPanel("GO analizė",
            p("Description goes here..."),
            shinydashboard::box(
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
            shinydashboard::box(
              width = 12,
              withLoader(DT::dataTableOutput(outputId = "table1"),
                         type = "html", loader = "dnaspin"),
              downloadButton("downloadData", "Download")
            )
          )
        )
      )
    )
  )
)
  






shinyApp(ui = ui, server = server)

# nolint end
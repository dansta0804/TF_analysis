# nolint start
library(pacman)
p_load(shiny, data.table, rtracklayer, ggplot2, ggthemes, plyranges, ggpubr,
       BRGenomics, reshape2, plotly, heatmaply, dplyr, gplots, genomation,
       Biostrings, scales, GenomicRanges, DT, shinythemes, shinycustomloader,
       ggseqlogo, ChIPseeker, tools, reactable, annotables, enrichplot,
       clusterProfiler, shinyalert, rjson)

PROJECT <- "/home/daniele/Desktop/IV_course/II_semester/TF_analysis/"
# setwd("/home/daniele/Desktop/IV_course/II_semester/TF_analysis/Scripts/App/") 
# source(paste0("/home/daniele/Desktop/IV_course/II_semester/TF_analysis/Scripts/App/ui.R"))
# source(paste0(PROJECT, "Scripts/App/server.R"))

options(scipen = 100)
options(shiny.maxRequestSize = 300 * 1024 ^ 2)

FUNCTIONS     <- paste0(PROJECT, "Scripts/functions.R")
INPUTS        <- paste0(PROJECT, "Analyses/Tbx5_analysis_I/Inputs/BED/")
INPUT         <- paste0(PROJECT, "Input/")
RESULTS       <- paste0(PROJECT, "Results/")
HOMER_RESULTS <- paste0(RESULTS, "HOMER/")
source(FUNCTIONS)

server <- function(input, output, session) {
  observe({
    result <- fromJSON(file = paste0(INPUT, "genome_data.json"))
    print(input$organism)
    print(result[[1]][5])
    for (genome in 1:length(result)) {
      if (input$organism == "Nenurodyta") {
        # shinyalert("Nurodykite organizmą!", type = "error",
        #            confirmButtonText = "Rinktis organizmą")
        return()
      } else if (input$organism == result[[genome]][5]) {
        if (system.file(package = result[[genome]][1]) != 0 &&
            system.file(package = result[[genome]][2]) != 0 &&
            system.file(package = result[[genome]][3]) != 0) {

              library(pacman)
          p_load(result[[genome]][1], result[[genome]][2], result[[genome]][3])
            } else {
              BiocManager::install(result[[genome]][1], update = FALSE)
              BiocManager::install(result[[genome]][2], update = FALSE)
              BiocManager::install(result[[genome]][3], update = FALSE)
            }
      } else { next }
    }
  })

  observe({
    output$table0 <- renderTable({
      input_data <- data.frame(matrix(ncol = 2, nrow = 0), row.names = NULL)
      colnames(input_data) <- c("Reikalinga informacija", "Reikšmė")

      checkable_values <-
        c(input$tf_options, input$organism, input$pwm)
      col1_text <-
        c("Transkripcijos faktorius:", "Organizmas, iš kurio išgauti duomenys:",
          "PWM motyvo logotipas:")

      for (value in 1:length(checkable_values)) {
        if (checkable_values[value] != "Nenurodyta") {
          input_data[nrow(input_data) + 1, ] <-
            c(col1_text[value], checkable_values[value])
        } else {
          input_data[nrow(input_data) + 1, ] <-
            c(col1_text[value], "Trūksta informacijos!")
        }
      }
      input_data
    })
  })

  ############################# FORMAT CONVERSION #############################
  observe({
    req(input$provided_file)
    input_files <- input$provided_file[, c("name")]
    print(input_files)

    # test_path <- system.file("tests", package = "rtracklayer")
  # format_detection(input_files)

  ## Returns ranges with all fields
  


  

    # if (input$file_format == "bed") {
    #   colnames(converted_table) <-
    #     c("Old file name", "New file name", "Graph name")
        
    #   for(sample in 1:length(input$bigbed[, c("name")])) {
    #     names <- input$bigbed[sample, 'name']
    #     row <- c(names, "-", paste0("Sample", sample, "_", input$organism))
    #     converted_table[nrow(converted_table) + 1, ] <- row
    #   }
    #   output$samples <- DT::renderDataTable({converted_table})
    # } else if (input$file_format == "bw") {
    #     colnames(converted_table) <-
    #       c("Old file name", "New file name", "Graph name")
        
    #     for(sample in 1:length(input$bigbed[, c("name")])) {
    #       names <- input$bigbed[sample, 'name']
    #       system(paste("/home/daniele/Tools/bigWigToBedGraph",
    #              paste0(BIGWIG, "Mus_musculus/", names),
    #              paste0(BEDGRAPH, "Mus_musculus/", names, ".bedGraph")))
    #       row <- c(names, paste0(names, ".bedGraph"),
    #                paste0("Sample", sample, "_", input$organism))
    #       converted_table[nrow(converted_table) + 1, ] <- row
    #     }
    #     output$samples <- DT::renderDataTable({converted_table})
    # }
  })

  converted_table <- data.frame(matrix(ncol = 3, nrow = 0))
  observe({
    req(input$bigbed)
    
    colnames(converted_table) <-
      c("Originalus pavadinimas", "Grafikų pavadinimas", "Mėginio dydis (Mb)")
    
    for(sample in 1:length(input$bigbed[, c("name")])) {
      names <- input$bigbed[sample, 'name']
      size <- as.numeric(input$bigbed[sample, 'size']) / 1000000
      row <- c(names, paste0("Sample", sample, "_", input$organism), size)
      converted_table[nrow(converted_table) + 1, ] <- row
    }
      
    output$samples <- DT::renderDataTable({converted_table}, server = TRUE)
    output$samples2 <-
      DT::renderDataTable({
        converted_table[, c("Grafikų pavadinimas", "Mėginio dydis (Mb)")]}, server = TRUE)
    output$samples3 <-
      DT::renderDataTable({
        converted_table[, c("Grafikų pavadinimas", "Mėginio dydis (Mb)")]}, server = TRUE)
  
    ########################### GENOMIC DATA QUALITY #######################
    # GENOMIC DATA QUALITY - PEAK COUNT (PLOT 1):
    output$plot1 <- renderPlot({
      req(input$samples2_rows_selected)
      bigbed_files <- list()
      row_index <- input$samples2_rows_selected

      for (rows in 1:length(row_index)) {
        bigbed_files[[rows]] <-
          read.table(file = paste0(INPUTS, converted_table[row_index[rows],
                                                           "Originalus pavadinimas"]))
        names(bigbed_files)[rows] <-
          converted_table[row_index[rows], "Grafikų pavadinimas"]
      }

      if (length(row_index) == 0) {
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
          coord_cartesian(ylim = c(0, as.numeric(max(peaks$Peak_count)) + 10000)) +
          geom_text(aes(label = as.numeric(Peak_count)),
              color = "#030101", size = 5, vjust = -1, fontface = 2) +
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
      req(input$samples2_rows_selected)
      bigbed_files <- list()
      grl <- GRangesList()
      row_index <- input$samples2_rows_selected

      for(rows in 1:length(row_index)) {
        bigbed_files[[rows]] <-
          read.table(file = paste0(INPUTS, converted_table[row_index[rows],
                                                            "Originalus pavadinimas"]))
        if (length(colnames) > 4) {
          colnames(bigbed_files[[rows]]) <-
            c("chrom", "start", "end", "name",
              rep(paste0("other", 1:length(colnames) - 4)))
        } else {
          colnames(bigbed_files[[rows]]) <-
            c("chrom", "start", "end", "other")
        }

        grl[[rows]] <- makeGRangesFromDataFrame(bigbed_files[[rows]],
                                                keep.extra.columns = TRUE)
        names(grl)[rows] <- converted_table[row_index[rows], "Grafikų pavadinimas"]
      }

      if (length(row_index) == 0) {
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
      req(input$samples2_rows_selected)
      bigbed_files <- list()
      grl <- list()
      row_index <- input$samples2_rows_selected

      for(rows in 1:length(row_index)) {
        bigbed_files[[rows]] <-
          read.table(file = paste0(INPUTS, converted_table[row_index[rows],
                                                            "Originalus pavadinimas"]))
        if (length(colnames) > 4) {
          colnames(bigbed_files[[rows]]) <-
            c("chrom", "start", "end", "name",
              rep(paste0("other", 1:length(colnames) - 4)))
        } else {
          colnames(bigbed_files[[rows]]) <-
            c("chrom", "start", "end", "other")
        }

        grl[[rows]] <- makeGRangesFromDataFrame(bigbed_files[[rows]],
                                                keep.extra.columns = TRUE)
        names(grl)[rows] <- converted_table[row_index[rows], "Grafikų pavadinimas"]
      }


      if (length(row_index) == 0) {
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
      req(input$samples2_rows_selected)
      bigbed_files <- list()
      grl <- GRangesList()
      row_index <- input$samples2_rows_selected

      for(rows in 1:length(row_index)) {
        bigbed_files[[rows]] <-
          read.table(file = paste0(INPUTS, converted_table[row_index[rows],
                                                            "Originalus pavadinimas"]))
        if (length(colnames) > 4) {
          colnames(bigbed_files[[rows]]) <-
            c("chrom", "start", "end", "name",
              rep(paste0("other", 1:length(colnames) - 4)))
        } else {
          colnames(bigbed_files[[rows]]) <-
            c("chrom", "start", "end", "other")
        }

        grl[[rows]] <- makeGRangesFromDataFrame(bigbed_files[[rows]],
                                                keep.extra.columns = TRUE)
        names(grl)[rows] <- converted_table[row_index[rows], "Grafikų pavadinimas"]
      }


      if (length(row_index) == 0) {
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
      req(input$samples2_rows_selected)
      bigbed_files <- list()
      grl <- GRangesList()
      row_index <- input$samples2_rows_selected

      for(rows in 1:length(row_index)) {
        bigbed_files[[rows]] <-
          read.table(file = paste0(INPUTS, converted_table[row_index[rows],
                                                            "Originalus pavadinimas"]))
        if (length(colnames) > 4) {
          colnames(bigbed_files[[rows]]) <-
            c("chrom", "start", "end", "name",
              rep(paste0("other", 1:length(colnames) - 4)))
        } else {
          colnames(bigbed_files[[rows]]) <-
            c("chrom", "start", "end", "other")
        }

        grl[[rows]] <- makeGRangesFromDataFrame(bigbed_files[[rows]],
                                                keep.extra.columns = TRUE)
        names(grl)[rows] <- converted_table[row_index[rows], "Grafikų pavadinimas"]
      }

      if (length(row_index) == 0) {
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
      req(input$samples2_rows_selected)
      bigbed_files <- list()
      plots <- list()
      grl <- GRangesList()
      row_index <- input$samples2_rows_selected

      for(rows in 1:length(row_index)) {
        bigbed_files[[rows]] <-
          read.table(file = paste0(INPUTS, converted_table[row_index[rows],
                                                            "Originalus pavadinimas"]))
        if (length(colnames) > 4) {
          colnames(bigbed_files[[rows]]) <-
            c("chrom", "start", "end", "name",
              rep(paste0("other", 1:length(colnames) - 4)))
        } else {
          colnames(bigbed_files[[rows]]) <-
            c("chrom", "start", "end", "other")
        }

        grl[[rows]] <- makeGRangesFromDataFrame(bigbed_files[[rows]],
                                                keep.extra.columns = TRUE)
        names(grl)[rows] <- converted_table[row_index[rows], "Grafikų pavadinimas"]
      }

      if (length(row_index) == 0) {
        return()
      } else {
        mm_known_genes <- TxDb.Mmusculus.UCSC.mm10.knownGene
        promoter <- getPromoters(TxDb = mm_known_genes, upstream = 3000,
                                 downstream = 3000)
        
        for (peak_file in 1:length(grl)) {
          remove_y <- theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(face = "bold")
          )

          peak <- makeGRangesFromDataFrame(grl[[peak_file]],
                                           keep.extra.columns = TRUE)
          
          tagMatrix <- getTagMatrix(peak, windows = promoter)
          
          plots[[peak_file]] <- plotAvgProf(tagMatrix, xlim = c(-3000, 3000),
                                            xlab = "Genomic Region (5'->3')",
                                            ylab = "Read Count Frequency") +
          remove_y +
          ggtitle(names(grl)[peak_file])
        }
        ggarrange(plotlist = plots, nrow = 1, common.legend = TRUE,
                  widths = 5, legend = "bottom")
      }
    })


    # if (length(row_index) > 1) {
    #     shinyalert("Pasirinkite tik vieną mėginį!", type = "error")
    #     return()
    #   } else if (length(row_index) == 0) {
        
    #   }

    

  ########################## GENOMIC DATA ANALYSIS #######################
    output$text <- renderText({
      req(input$tf_options)
      paste0(input$tf_options)
    })

    # GENOMIC DATA ANALYSIS - TF MOTIF (PLOT 7):
    output$plot7 <- renderPlot({
      req(input$pwm)
      mpwm <- read.table(file = input$pwm$datapath)
      mpwm <- t(mpwm)

      # Setting matrix rownames:
      rownames(mpwm) <- c("A", "C", "G", "T")
      ggseqlogo(mpwm)
    })

    # GENOMIC DATA ANALYSIS - PWM MATRIX MATCHES (PLOT 8):
    output$plot8 <- renderPlot({
      req(input$samples3_rows_selected)
      req(input$pwm)
      bigbed_files <- list()
      row_index <- input$samples3_rows_selected

      for(rows in 1:length(row_index)) {
        bigbed_files[[rows]] <-
          read.table(file = paste0(INPUTS, converted_table[row_index[rows],
                                                            "Originalus pavadinimas"]))
        
        if (length(colnames) > 4) {
          colnames(bigbed_files[[rows]]) <-
            c("chrom", "start", "end", "name",
              rep(paste0("other", 1:length(colnames) - 4)))
        } else {
          colnames(bigbed_files[[rows]]) <- c("chrom", "start", "end", "other")
        }

        bigbed_files[[rows]] <-
          makeGRangesFromDataFrame(bigbed_files[[rows]],
                                   keep.extra.columns = TRUE)
        names(bigbed_files)[rows] <- converted_table[row_index[rows],
                                                     "Grafikų pavadinimas"]
      }

      if (length(input$samples3_rows_selected) == 0) {
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
    
    go_data <- reactive({
      options(scipen = 0)
      req(input$samples3_rows_selected)
      row_index <- input$samples3_rows_selected
      bigbed_files <- list()

      if (length(row_index) > 1) {
        shinyalert("Pasirinkite tik vieną mėginį!", type = "error",
                   confirmButtonText = "Rinktis mėginį")
        return()
      }

      for(rows in 1:length(row_index)) {
        bigbed_files[[rows]] <-
          read.table(file = paste0(INPUTS, converted_table[row_index[rows],
                                                            "Originalus pavadinimas"]))
        
        if (length(colnames(bigbed_files[[rows]])) > 4) {
          colnames(bigbed_files[[rows]]) <-
          c("seqnames", "start", "end", "name", "uniprot_id", "abs_summit", "pileup",
              "p_value", "fold_enrichment", "q_value", "ID")
            # c("chrom", "start", "end", "name",
            #   rep(paste0("other", 1:length(colnames) - 4)))
        } else {
          colnames(bigbed_files[[rows]]) <- c("seqnames", "start", "end", "other")
        }

        bigbed_files[[rows]] <-
          makeGRangesFromDataFrame(bigbed_files[[rows]],
                                   keep.extra.columns = TRUE)
        names(bigbed_files)[rows] <- converted_table[row_index[rows],
                                                     "Grafikų pavadinimas"]
      }

      if (length(input$samples3_rows_selected) == 0) {
        return()
      } else {
        annotated_peaks <-
          annotate_peaks(bigbed_files, grcm38,
                         TxDb.Mmusculus.UCSC.mm10.knownGene, "org.Mm.eg.db")

        mmgene_ranges <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)

        # Adding 'gene_symbol' column to gene_ranges GRange:
        mmgene_ranges$gene_symbol <-
            mapIds(org.Mm.eg.db, keys = mmgene_ranges$gene_id, column = "SYMBOL",
                  keytype = "ENTREZID", multiVals = "first")

        mm_genes <- GRangesList()

        for (sample in seq_along(annotated_peaks)) {
            selected_genes <-
                tolower(unique(annotated_peaks[[sample]]$gene_symbol)) %in%
                tolower(mmgene_ranges$gene_symbol)
                                    
            genes <- mmgene_ranges[selected_genes]
            mm_genes[[sample]] <- genes 
        }

        merged_expanded <- list()
        plots <- list()

        for (l in 1:length(bigbed_files)) {
            merged_data <- merge(annotated_peaks[[l]], bigbed_files[[l]])
            expanded_merged <-
                bitr(merged_data$gene_symbol, fromType = "SYMBOL",
                    toType = c("ENTREZID"), OrgDb = org.Mm.eg.db)
            merged_expanded[[l]] <-
                merge(merged_data, expanded_merged, by.x = "gene_symbol",
                      by.y = "SYMBOL", all.x = TRUE) %>%
                as.data.frame() %>%
                makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
        }
        merged_expanded
      }
    })

    # GENOMIC DATA ANALYSIS - GO ANALYSIS (TABLE 1):
    output$table1 <- renderReactable({
        pl <- enrichGO(gene = go_data()[[1]]$ENTREZID,
                                OrgDb = org.Mm.eg.db, ont = "CC",
                                pAdjustMethod = "BH", pvalueCutoff  = 0.01,
                                qvalueCutoff  = 0.05, readable = TRUE)

        pl <- as.data.frame(pl)
        pl$Status <- "Peržiūrėti genų sąrašą"
        pl <-
          pl %>%
          dplyr::select(c("ID", "Description", "GeneRatio", "Count",
                          "Status", "geneID")) %>%
          rename("GO ID" = "ID", "Apibūdinimas" = "Description",
                 "Genų santykis" = "GeneRatio", "Genų skaičius" = "Count",
                 "Peržiūra" = "Status")
        
        reactable(
          pl, searchable = FALSE, showSortable = TRUE, rownames = FALSE,
          pagination = FALSE, highlight = TRUE, height = 750,
          defaultColDef = colDef(
            align = "center",
            minWidth = 70
            # headerStyle = list(background = "#f7f7f8")
          ),
          columns = list(
            `Peržiūra` = colDef(
              style = function(value) {
                list(color = "#500909", fontWeight = "bold")
              },
              details = function(value) {
                genes <- unlist(strsplit(pl$geneID, split = "/")[value])
                add_NA <- 8 - (length(genes) - 8 * (length(genes) %/% 8))
                genes <- c(genes, rep(" ", add_NA))
                tb <- data.frame(matrix(genes, ncol = 8), ncol = 8)
                colnames(tb) <- c(rep(paste0("Genai ", 1:8)))

                htmltools::div(
                  style = "padding: 1rem",
                  reactable(tb[, 1:8], outlined = TRUE, fullWidth = TRUE,
                    defaultColDef = colDef(
                      align = "center",
                      minWidth = 70
                      # headerStyle = list(background = "#f7f7f8")
                    )
                  )
                )
              }
            ),
            geneID = colDef(show = FALSE)
          )
        )   
    })

    # GENOMIC DATA ANALYSIS - GO ANALYSIS (ACYCLIC GRAPH) (PLOT 9):
    output$plot9 <- renderPlot({
      pl <- enrichGO(gene = go_data()[[1]]$ENTREZID,
                                OrgDb = org.Mm.eg.db, ont = "CC",
                                pAdjustMethod = "BH", pvalueCutoff  = 0.01,
                                qvalueCutoff  = 0.05, readable = TRUE)
      goplot(pl)
    })

    # GENOMIC DATA ANALYSIS - GO ANALYSIS (TREE POT) (PLOT 10):
    output$plot10 <- renderPlot({
      pl <- enrichGO(gene = go_data()[[1]]$ENTREZID,
                                OrgDb = org.Mm.eg.db, ont = "CC",
                                pAdjustMethod = "BH", pvalueCutoff  = 0.01,
                                qvalueCutoff  = 0.05, readable = TRUE)
      pl_modified <- setReadable(pl, 'org.Mm.eg.db', 'ENTREZID')
      pl_modified <- pairwise_termsim(pl_modified)
      p <- treeplot(pl_modified, cluster.params = list(method = "average"),
                    xlim = c(0, 30))
      p
    })

    # GENOMIC DATA ANALYSIS - DE NOVO MOTIFS (TABLE 1):
    output$table2 <- DT::renderDataTable({
      req(input$samples3_rows_selected)
      bigbed_files <- list()
      names <- c()
      row_index <- input$samples3_rows_selected

      if (length(row_index) > 1) {
        shinyalert("Pasirinkite tik vieną mėginį!", type = "error",
                   confirmButtonText = "Rinktis mėginį")
        return()
      }

      for (rows in 1:length(row_index)) {
        bigbed_files[[rows]] <-
          read.table(file = paste0(INPUTS, converted_table[row_index[rows],
                                                           "Originalus pavadinimas"]))
        names(bigbed_files)[rows] <-
          converted_table[row_index[rows], "Grafikų pavadinimas"]
      }

      if (length(input$samples3_rows_selected) == 0) {
        return()
      } else {
        homer_motifs <- c()
        for (file in 1:length(bigbed_files)) {
          folder_name <-
            file_path_sans_ext(basename(converted_table[row_index[file],
                                                        "Originalus pavadinimas"]))
          filename <- converted_table[row_index[file], "Originalus pavadinimas"]
          output_path <- paste0(HOMER_RESULTS, folder_name)
          
          
          if (!file.exists(output_path)) {
            dir.create(output_path)
          }
          
          if (!file.exists(paste0(output_path, "/knownResults.txt"))) {
            system(paste("findMotifsGenome.pl",
                       paste0(INPUTS, filename), "mm10",
                       paste0(output_path, "/"),
                 "-size 200 -mask"), intern = FALSE)
          } else {
            homer_motifs <- read.table(paste0(output_path, "/knownResults.txt"),
                                       sep = "\t", header = FALSE, skip = 1)
          }

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
        }
        homer_motifs
      }
    })

  #   url <- a("Google Homepage", href="https://vult-my.sharepoint.com/:f:/g/personal/daniele_stasiunaite_mif_stud_vu_lt/Eu4Rfi8-aqJBq-B_XcmeU_QBcdlQ7_Epb4yfy0VIIa5K_g?e=CFFJSG")
  #   output$link <- renderUI({
  #     tagList("URL link:", url)
  #   })

  #   observeEvent(input$do, {
  #   session$sendCustomMessage(type = 'testmessage',
  #     message = 'Thank you for clicking')
  # })
  })
}

ui <- navbarPage("ChIP sekoskaitos analizės", theme = shinytheme("cosmo"),
  useShinyalert(),
  includeCSS(paste0(PROJECT, "Scripts/styles.css")),
  tabPanel("Formatų konversija",
    sidebarLayout(
      sidebarPanel(
        width = 4,
        fileInput(
          inputId = "provided_file",
          label = "Įkelkite ChIP sekoskaitos duomenų failą (-us)*:",
          multiple = TRUE,
          buttonLabel = "Ieškoti failo",
          placeholder = "Failas nepasirinktas"
        ),
        radioButtons(
          inputId = "file_format",
          label = "Nurodykite įvestų duomenų formatą:",
          choices = c("BED" = "bed", "BigWig" = "bw", "BigBed" = "bb")
        ),
      ),
      mainPanel(
        width = 8
      )
    )
  ),
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
        selectInput(
          inputId = "tf_options",
          label = "Pasirinkite transkripcijos faktorių:",
          choices = c("Tbx5", "GATA3", "Tcf21", "Nenurodyta"),
          selected = "Nenurodyta"
        ),
        selectInput(
          inputId = "organism",
          label = "Nurodykite, iš kokio organizmo išgauti mėginiai:",
          choices = c(
            "Mus musculus" = "mm", "Homo sapiens" = "hs",
            "Rattus norvegicus" = "rn", "Danio rerio" = "dr",
            "Bos taurus" = "bt", "Drosophila melanogaster" = "dm",
            "Gallus gallus" = "gg", "Macaca mulatta" = "mmul",
            "Pan troglodytes" = "pt", "Sus scrofa" = "ss", "Nenurodyta" = "nn"
          ),
          selected = "nn"
        ),
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
        p("Pateiktose lentelėse apibendrinta mėginių bei kita įvesta
           informacija:"),
        tableOutput("table0"),
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
        DT::dataTableOutput("samples2"),
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
        DT::dataTableOutput("samples3")
        # ,
        # radioButtons(
        #   inputId = "genome",
        #   label = "Pasirinkite genomą:",
        #   choices = c("Homo sapiens" = "hg", "Mus musculus" = "mm",
        #               "Danio rerio" = "dr")
        # )
      ),
      mainPanel(
        width = 8,
        tabsetPanel(
          tabPanel("Transkripcijos faktoriaus motyvas", 
            p("Pateiktame paveiksle pavaizduotas įkeltą pozicinę svorių matricą
               atitinkančio transkripcijos faktoriaus motyvo sekos logotipas."),
            h4(strong("Transkripcijos faktorius: "), textOutput("text")),
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
            tabsetPanel(
              id = "go",
              tabPanel("GO lentelė",
                id = "go_table",
                p("Pateiktoje lentelėje nurodyti gauti genų ontologijos
                   rezultatai:"),
                width = 12,
                withLoader(reactableOutput("table1"), type = "html",
                           loader = "dnaspin")
              ),
              tabPanel("GO aciklinis grafas",
                id = "go_graph",
                p("Gauto GO analizės rezultato pavaizdavimas kryptiniu
                   acikliniu grafu:"),
                width = 12,
                withLoader(plotOutput("plot9"), type = "html",
                           loader = "dnaspin")
              ),
              tabPanel("GO medžio struktūra",
                id = "go_tree",
                p("Gauto GO analizės rezultato pavaizdavimas, pritaikius
                  hierarchinį klasterizavimą:"),
                width = 12,
                withLoader(plotOutput("plot10"), type = "html",
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
              withLoader(DT::dataTableOutput(outputId = "table2"),
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
        width = 4
      ),
      mainPanel(
        width = 8
      )
    )
  )
)
  






shinyApp(ui = ui, server = server)

# nolint end
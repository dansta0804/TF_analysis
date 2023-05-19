# nolint start
library(pacman)
p_load(shiny, data.table, rtracklayer, ggplot2, ggthemes, plyranges, ggpubr,
       BRGenomics, reshape2, plotly, heatmaply, dplyr, gplots, genomation,
       Biostrings, scales, GenomicRanges, DT, shinythemes, shinycustomloader,
       ggseqlogo, ChIPseeker, tools, reactable, annotables, enrichplot,
       clusterProfiler, shinyalert, rjson)

options(scipen = 100)
options(shiny.maxRequestSize = 300 * 1024 ^ 2)

PROJECT <- "./"
FUNCTIONS <- paste0(PROJECT, "Scripts/functions.R")
INPUT <- paste0(PROJECT, "Input/")
RESULTS <- paste0(PROJECT, "Results/")
HOMER_RESULTS <- paste0(RESULTS, "HOMER/")

source(FUNCTIONS)

server <- function(input, output, session) {
  observe({
    result <- fromJSON(file = paste0(INPUT, "genome_data.json"))
    req(input$organism)



    genome <- input$organism
    bsgenome <- unlist(result[[genome]][1])
    txdb_g <- unlist(result[[genome]][2])
    orgdb <- unlist(result[[genome]][3])

    if (length(bsgenome) != 0 && length(txdb_g) != 0 && length(orgdb) != 0) {
        if (system.file(package = bsgenome) != "" &&
            system.file(package = txdb_g) != "" &&
            system.file(package = orgdb) != "" ) {
                library(bsgenome, character.only = TRUE)
                library(txdb_g, character.only = TRUE)
                library(orgdb, character.only = TRUE)
        } else {
            BiocManager::install(bsgenome, update = FALSE)
            BiocManager::install(txdb, update = FALSE)
            BiocManager::install(orgdb, update = FALSE)
        }
    }
  })

  json <- fromJSON(file = paste0(INPUT, "genome_data.json"))

  output$table0 <- renderTable({
    req(input$organism)
    req(input$tf_options)
    
    input_data <- data.frame(matrix(ncol = 2, nrow = 0), row.names = NULL)
    colnames(input_data) <- c("Reikalinga informacija", "Reikšmė")

    checkable_values <-
      c(input$tf_options, input$organism)
    col1_text <-
      c("Transkripcijos faktorius:", "Organizmas, iš kurio išgauti duomenys:")

    for (value in 1:length(checkable_values)) {
      if (checkable_values[value] != "Nenurodyta") {
        input_data[nrow(input_data) + 1, ] <-
          c(col1_text[value], checkable_values[value])
      } else {
        input_data[nrow(input_data) + 1, ] <-
          c(col1_text[value], checkable_values[value])
      }
    }
    input_data
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

  converted_table <- data.frame(matrix(ncol = 4, nrow = 0))
  observe({
    req(input$bigbed)
    
    colnames(converted_table) <-
      c("Originalus pavadinimas", "Grafikų pavadinimas",
        "Mėginio dydis (pikais)", "Kelias")
    
    for(sample in 1:length(input$bigbed[, c("name")])) {
      names <- input$bigbed[sample, 'name']
      path <- input$bigbed[sample, 'datapath']
      # size <- as.numeric(input$bigbed[sample, 'size']) / 1000000
      file <- read.table(file = path)
      size <- spaces(length(rownames(file)))
      row <- c()
      organism <- input$organism

      if (organism == "Nenurodyta") {
        row <- c(names, paste0("Sample", sample, "_", "nn"), size, path)
      } else {
        abbrv <- unlist(json[[organism]][5])
        row <- c(names, paste0("Sample", sample, "_", abbrv), size, path)
      }
      converted_table[nrow(converted_table) + 1, ] <- row
    }
      
    output$samples <-
      DT::renderDataTable({
        converted_table[, c("Originalus pavadinimas", "Grafikų pavadinimas",
                            "Mėginio dydis (pikais)")]}, server = TRUE)
    output$samples2 <-
      DT::renderDataTable({
        converted_table[, c("Grafikų pavadinimas", "Mėginio dydis (pikais)")]},
        server = TRUE)
    output$samples3 <-
      DT::renderDataTable({
        converted_table[, c("Grafikų pavadinimas", "Mėginio dydis (pikais)")]},
        server = TRUE)
  
    ########################### GENOMIC DATA QUALITY #######################
    # GENOMIC DATA QUALITY - PEAK COUNT (PLOT 1):
    plot1 <- reactive({
      req(input$samples2_rows_selected)
      bigbed_files <- list()
      row_index <- input$samples2_rows_selected

      for (rows in 1:length(row_index)) {
        bigbed_files[[rows]] <-
          read.table(file = converted_table[row_index[rows], "Kelias"])
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
    output$plot1 <- renderPlot({plot1()})
    output$download1 <- downloadHandler(
      filename = function() {
        paste("Pikų_skaičius_mėginiuose", "png", sep = ".")
      },
      content = function(file){
        ggsave(file, plot1(), device = "png")
      }
    )

    # GENOMIC DATA QUALITY - PEAK COUNT DISTRIBUTION BY CHROMOSOME (PLOT 2):
    plot2 <- reactive({
      req(input$samples2_rows_selected)
      bigbed_files <- list()
      grl <- GRangesList()
      row_index <- input$samples2_rows_selected

      genome <- input$organism
      unlisted_json <- unlist(json[[genome]][6])
      start_end <- unlist(strsplit(unlisted_json[1], ":"))

      chromosomes <-
        c(rep(paste0("chr", seq(start_end[1], start_end[2]))),
          rep(paste0("chr", unlisted_json[2:(length(unlisted_json))])))

      for(rows in 1:length(row_index)) {
        bigbed_files[[rows]] <-
          read.table(file = converted_table[row_index[rows], "Kelias"])

        col_len <- length(colnames(bigbed_files[[rows]]))

        if (col_len > 4) {
          colnames(bigbed_files[[rows]]) <-
            c("chrom", "start", "end", "name",
              rep(paste0("other", 1:(col_len - 4))))
        } else {
          colnames(bigbed_files[[rows]]) <- c("chrom", "start", "end", "other")
        }

        filtered_df <-
          bigbed_files[[rows]] %>%
          filter(chrom %in% chromosomes)

        grl[[rows]] <- makeGRangesFromDataFrame(filtered_df,
                                                keep.extra.columns = TRUE)
        names(grl)[rows] <-
          converted_table[row_index[rows], "Grafikų pavadinimas"]
      }

      if (length(row_index) == 0) {
        return()
      } else {
        chr_len <- unlist(json[[genome]][9])

        peak_counts <-
          lapply(names(grl), count_peaks, objects = grl,
                 chr_abbr = chromosomes, chr_lengths = chr_len) %>%
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
    output$plot2 <- renderPlot({plot2()})
    output$download2 <- downloadHandler(
      filename = function() {
        paste("Pikų_pasiskirstymas_chromosomose", "png", sep = ".")
      },
      content = function(file){
        ggsave(file, plot2(), device = "png")
      }
    )

    # GENOMIC DATA QUALITY - OVERLAPS BETWEEN SAMPLES (JACCARD) (PLOT 3):
    plot3 <- reactive({
      req(input$samples2_rows_selected)
      bigbed_files <- list()
      grl <- GRangesList()
      row_index <- input$samples2_rows_selected

      for(rows in 1:length(row_index)) {
        bigbed_files[[rows]] <-
          read.table(file = converted_table[row_index[rows], "Kelias"])
        col_len <- length(colnames(bigbed_files[[rows]]))

        if (col_len > 4) {
          colnames(bigbed_files[[rows]]) <-
            c("chrom", "start", "end", "name",
              rep(paste0("other", 1:(col_len - 4))))
        } else {
          colnames(bigbed_files[[rows]]) <- c("chrom", "start", "end", "other")
        }

        grl[[rows]] <- makeGRangesFromDataFrame(bigbed_files[[rows]],
                                                keep.extra.columns = TRUE)
        names(grl)[rows] <-
          converted_table[row_index[rows], "Grafikų pavadinimas"]
      }

      if (length(row_index) == 0) {
        return()
      } else {
        coef_matrix <- matrix(nrow = length(grl), ncol = length(grl))

        # Calculating Jaccard coefficient for sample pair:
        for (i in 1:length(grl)) {
            for (y in 1:length(grl)) {
                coef_matrix[i, y] = jaccard(grl, i, y)
            }
        }

        # Setting colnames and rownames for the matrix:
        colnames(coef_matrix) <- names(grl)
        rownames(coef_matrix) <- names(grl)

        coef_mat1 <- coef_matrix
        coef_mat2 <- coef_matrix

        # Passing Jaccard coefficients to matrix except for the diagonal - it
        # contains 'NA':
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
          guides(fill = guide_colourbar(title = "Koeficientas", face = "bold")) +
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
    output$plot3 <- renderPlot({plot3()})
    output$download3 <- downloadHandler(
      filename = function() {
        paste("Mėginių_panašumas", "png", sep = ".")
      },
      content = function(file){
        ggsave(file, plot3(), device = "png")
      }
    )

    # GENOMIC DATA QUALITY - GENOMIC DISTRIBUTION (PLOT 4):
    plot4 <- reactive({
      req(input$samples2_rows_selected)
      bigbed_files <- list()
      grl <- GRangesList()
      row_index <- input$samples2_rows_selected

      for(rows in 1:length(row_index)) {
        bigbed_files[[rows]] <-
          read.table(file = converted_table[row_index[rows], "Kelias"])

        col_len <- length(colnames(bigbed_files[[rows]]))

        if (col_len > 4) {
          colnames(bigbed_files[[rows]]) <-
            c("chrom", "start", "end", "name",
              rep(paste0("other", 1:(col_len - 4))))
        } else {
          colnames(bigbed_files[[rows]]) <- c("chrom", "start", "end", "other")
        }

        grl[[rows]] <- makeGRangesFromDataFrame(bigbed_files[[rows]],
                                                keep.extra.columns = TRUE)
        names(grl)[rows] <-
          converted_table[row_index[rows], "Grafikų pavadinimas"]
      }

      if (length(row_index) == 0) {
        return()
      } else {
        genome <- input$organism
        mm_known_genes <- unlist(json[[genome]][2])
        peak_annotations <- lapply(grl, annotatePeak, TxDb = get(mm_known_genes),
                                   tssRegion = c(-3000, 3000), verbose = FALSE)

        plotAnnoBar(peak_annotations, ylab = "Procentinė dalis (%)", title = "")
      }
    })
    output$plot4 <- renderPlot({plot4()})
    output$download4 <- downloadHandler(
      filename = function() {
        paste("Genominė_distribucija", "png", sep = ".")
      },
      content = function(file){
        ggsave(file, plot4(), device = "png")
      }
    )

    # GENOMIC DATA QUALITY - DISTANCE TO TSS (PLOT 5):
    plot5 <- reactive({
      req(input$samples2_rows_selected)
      bigbed_files <- list()
      grl <- GRangesList()
      row_index <- input$samples2_rows_selected

      for(rows in 1:length(row_index)) {
        bigbed_files[[rows]] <-
          read.table(file = converted_table[row_index[rows], "Kelias"])

        col_len <- length(colnames(bigbed_files[[rows]]))

        if (col_len > 4) {
          colnames(bigbed_files[[rows]]) <-
            c("chrom", "start", "end", "name",
              rep(paste0("other", 1:(col_len - 4))))
        } else {
          colnames(bigbed_files[[rows]]) <- c("chrom", "start", "end", "other")
        }

        grl[[rows]] <- makeGRangesFromDataFrame(bigbed_files[[rows]],
                                                keep.extra.columns = TRUE)
        names(grl)[rows] <-
          converted_table[row_index[rows], "Grafikų pavadinimas"]
      }

      if (length(row_index) == 0) {
        return()
      } else {
        genome <- input$organism
        mm_known_genes <- unlist(json[[genome]][2])
        peak_annotations <- lapply(grl, annotatePeak, TxDb = get(mm_known_genes),
                                   tssRegion = c(-3000, 3000), verbose = FALSE)
        plotDistToTSS(peak_annotations, ylab = "Atstumas", title = "")
      }
    })
    output$plot5 <- renderPlot({plot5()})
    output$download5 <- downloadHandler(
      filename = function() {
        paste("Atstumas_iki_TSS", "png", sep = ".")
      },
      content = function(file){
        ggsave(file, plot5(), device = "png")
      }
    )

    # GENOMIC DATA QUALITY - PEAK PROFILE (PLOT 6):
    output$plot6 <- renderPlot({
      req(input$samples2_rows_selected)
      bigbed_files <- list()
      plots <- list()
      grl <- GRangesList()
      row_index <- input$samples2_rows_selected

      for(rows in 1:length(row_index)) {
        bigbed_files[[rows]] <-
          read.table(file = converted_table[row_index[rows], "Kelias"])

        col_len <- length(colnames(bigbed_files[[rows]]))

        if (col_len > 4) {
          colnames(bigbed_files[[rows]]) <-
            c("chrom", "start", "end", "name",
              rep(paste0("other", 1:(col_len - 4))))
        } else {
          colnames(bigbed_files[[rows]]) <- c("chrom", "start", "end", "other")
        }

        grl[[rows]] <- makeGRangesFromDataFrame(bigbed_files[[rows]],
                                                keep.extra.columns = TRUE)
        names(grl)[rows] <-
          converted_table[row_index[rows], "Grafikų pavadinimas"]
      }

      if (length(row_index) == 0) {
        return()
      } else {
        genome <- input$organism
        mm_known_genes <- unlist(json[[genome]][2])
        promoter <- getPromoters(TxDb = get(mm_known_genes), upstream = 3000,
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
                                            xlab = "Regionas (5'->3')",
                                            ylab = "Read Count Frequency") +
          remove_y +
          ggtitle(names(grl)[peak_file])
        }
        ggarrange(plotlist = plots, nrow = 1, common.legend = TRUE,
                  widths = 5, legend = "bottom")
      }
    })

  ########################## GENOMIC DATA ANALYSIS #######################
    output$text <- renderText({
      req(input$tf_options)
      paste0(input$tf_options)
    })

    # GENOMIC DATA ANALYSIS - TF MOTIF (PLOT 7):
     output$plot7 <- renderPlot({
      req(input$pwm)
      mpwm <- read.table(file = input$pwm$datapath,  skip = 1)
      mpwm <- t(mpwm)

      # Setting matrix rownames:
      rownames(mpwm) <- c("A", "C", "G", "T")
      ggseqlogo(mpwm)
    })
    
    plot8 <- reactive({
      req(input$samples3_rows_selected)
      bigbed_files <- list()
      row_index <- input$samples3_rows_selected

      if (length(input$pwm) == 0) {
        shinyalert("PWM matrica neįkelta!", type = "error",
                   confirmButtonText = "Įkelti matricą")
        return()
      }

      if (input$organism == "Nenurodyta") {
        shinyalert("Nenurodytas organizmas!", type = "error",
                   confirmButtonText = "Nurodyti organizmą")
        return()
      }

      for (rows in 1:length(row_index)) {
        bigbed_files[[rows]] <-
          read.table(file = converted_table[row_index[rows], "Kelias"])
        
        col_len <- length(colnames(bigbed_files[[rows]]))

        if (col_len > 4) {
          colnames(bigbed_files[[rows]]) <-
            c("chrom", "start", "end", "name",
              rep(paste0("other", 1:(col_len - 4))))
        } else {
          colnames(bigbed_files[[rows]]) <- c("chrom", "start", "end", "other")
        }

        genome <- input$organism
        unlisted_json <- unlist(json[[genome]][6])
        start_end <- unlist(strsplit(unlisted_json[1], ":"))

        chromosomes <-
          c(rep(paste0("chr", seq(start_end[1], start_end[2]))),
            rep(paste0("chr", unlisted_json[2:(length(unlisted_json))])))

        filtered_df <-
          bigbed_files[[rows]] %>%
          filter(chrom %in% chromosomes)

        bigbed_files[[rows]] <-
          makeGRangesFromDataFrame(filtered_df, keep.extra.columns = TRUE)
        names(bigbed_files)[rows] <- converted_table[row_index[rows],
                                                     "Grafikų pavadinimas"]
      }

      if (length(input$samples3_rows_selected) == 0) {
        return()
      } else {
        mpwm <- read.table(file = input$pwm$datapath, skip = 1)
        mpwm <- t(mpwm)
        rownames(mpwm) <- c("A", "C", "G", "T")

        peak_sequences <- list()
        genome <- input$organism
        
        for (file in 1:length(bigbed_files)) {
          peak_sequences[[file]] <-
            getSeq(get(unlist(json[[genome]][1])), bigbed_files[[file]])
          names(peak_sequences)[file] <- names(bigbed_files)[file]
        }

        tbx5_motifs <- data.frame(matrix(nrow = 0, ncol = 4))
        colnames(tbx5_motifs) <-
          c("Sample", "Motif_count", "Peak_count", "Percentage")

        for (i in 1:length(peak_sequences)) {
          filename <- input$bigbed[[i, 'datapath']]
          name <- substring(names(peak_sequences[i]), 1, 11)

          motif <- find_motif_hits(as.character(peak_sequences[[i]]), mpwm)
          peaks <- calculate_peaks(filename)
          percentage <- round((motif / peaks) * 100, 2)

          data_row <- c(name, motif, peaks, as.numeric(percentage))
          tbx5_motifs[nrow(tbx5_motifs) + 1, ] <- data_row
        }

        # Reading a file that stores information about Tbx5 motif counts and
        # peak percentages:
        motif_data <- tbx5_motifs

        # Calling factor() function in order to maintain certain Sample
        # order:
        motif_data$Sample <-
                factor(motif_data$Sample, levels = motif_data$Sample)

        # 'Melting' the dataframe:
        melted_df <- melt(motif_data, id = c("Sample", "Percentage"))

        ggplot(melted_df,
               aes(fill = variable, y = as.numeric(value), x = Sample)) + 
          geom_bar(width = 0.4, size = 0.2, colour = "#3f2704",
                   stat = "identity", position = position_dodge(0.4)) +
          scale_fill_manual(values = c("#e3a15e", "#c7633b"),
                            labels = c(paste0(input$tf_options,
                                              " motyvų skaičius"),
                                       "Pikų skaičius")) +
          scale_y_continuous(labels = label_number(suffix = " K",
                                                   scale = 1e-3)) +
          geom_text(aes(label = ifelse(variable == "Motif_count",
                                       paste0(round(as.numeric(Percentage),
                                                    digits = 2), "%"), ""),
                        fontface = 2), vjust = 3.2, hjust = -0.3, size = 5) +
          guides(fill = guide_legend(title = "Spalvų paaiškinimas", size = 6)) +
          coord_cartesian(ylim = c(0, as.numeric(max(melted_df$value)) + 200000)) +
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
            panel.grid.major.y = element_line(colour = "#cab5b5", size = 0.3,
                                              linetype = "dashed"),
            panel.grid.minor.y = element_line(colour = "#cab5b5", size = 0.3,
                                              linetype = "dashed"),
            panel.grid.major.x = element_line(colour = "#cab5b5", size = 0.2,
                                              linetype = "longdash"),
            panel.grid.minor.x = element_line(colour = "#cab5b5",  size = 0.2,
                                              linetype = "longdash"),
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 11)) +
          coord_flip()
      }
    })

    output$plot8 <- renderPlot({plot8()})
    output$download8 <- downloadHandler(
      filename = function() {
        paste("PWM_atitikimai", "png", sep = ".")
      },
      content = function(file){
        ggsave(file, plot8(), device = "png")
      }
    )
    
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
          read.table(file = converted_table[row_index[rows], "Kelias"])
        
        col_len <- length(colnames(bigbed_files[[rows]]))

        if (col_len > 4) {
          colnames(bigbed_files[[rows]]) <-
            c("chrom", "start", "end", "name",
              rep(paste0("other", 1:(col_len - 4))))
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
        genome <- input$organism
        mm_known_genes <- unlist(json[[genome]][2])
        orgdb <- unlist(json[[genome]][3])
        ensembl_annot <- unlist(json[[genome]][7])

        annotated_peaks <-
          annotate_peaks(bigbed_files, get(ensembl_annot),
                         get(mm_known_genes), orgdb)
        # print(annotated_peaks)                 
        mmgene_ranges <- genes(get(mm_known_genes))

        # Adding 'gene_symbol' column to gene_ranges GRange:
        mmgene_ranges$gene_symbol <-
            mapIds(get(orgdb), keys = mmgene_ranges$gene_id, column = "SYMBOL",
                  keytype = "ENTREZID", multiVals = "first")

        # print(mmgene_ranges)
        mm_genes <- GRangesList()

        for (sample in seq_along(annotated_peaks)) {
            selected_genes <-
                tolower(unique(annotated_peaks[[sample]]$gene_symbol)) %in%
                tolower(mmgene_ranges$gene_symbol)
            # print(selected_genes)                        
            genes <- mmgene_ranges[selected_genes]
            mm_genes[[sample]] <- genes 
        }

        merged_expanded <- list()
        plots <- list()

        for (l in 1:length(bigbed_files)) {
            merged_data <- merge(annotated_peaks[[l]], bigbed_files[[l]])
            expanded_merged <-
                bitr(merged_data$gene_symbol, fromType = "SYMBOL",
                    toType = c("ENTREZID"), OrgDb = get(orgdb))
            merged_expanded[[l]] <-
                merge(merged_data, expanded_merged, by.x = "gene_symbol",
                      by.y = "SYMBOL", all.x = TRUE) %>%
                as.data.frame() %>%
                makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
        }
        merged_expanded
      }
    })

    genome <- reactive({unlist(json[[input$organism]][3])})

    # GENOMIC DATA ANALYSIS - GO ANALYSIS (TABLE 1):
    output$table1 <-
      renderReactable({find_ontologies_table(go_data(), genome(), "BP")})
    output$table2 <-
      renderReactable({find_ontologies_table(go_data(), genome(), "MF")})
    output$table3 <-
      renderReactable({find_ontologies_table(go_data(), genome(), "CC")})

    # GENOMIC DATA ANALYSIS - GO ANALYSIS (ACYCLIC GRAPH) (PLOT 9):
    output$plot9 <-
      renderPlot({find_ontologies_graph(go_data(), genome(), "BP")})
    output$plot10 <-
      renderPlot({find_ontologies_graph(go_data(), genome(), "MF")})
    output$plot11 <-
      renderPlot({find_ontologies_graph(go_data(), genome(), "CC")})

    # GENOMIC DATA ANALYSIS - GO ANALYSIS (TREE POT) (PLOT 10):
    output$plot12 <-
      renderPlot({find_ontologies_tree(go_data(), genome(), "BP")})
    output$plot13 <-
      renderPlot({find_ontologies_tree(go_data(), genome(), "MF")})
    output$plot14 <-
      renderPlot({find_ontologies_tree(go_data(), genome(), "CC")})

    # GENOMIC DATA ANALYSIS - DE NOVO MOTIFS (TABLE 1):
    output$table4 <- DT::renderDataTable({
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
          read.table(file = converted_table[row_index[rows], "Kelias"])
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
          filename <- converted_table[row_index[file], "Kelias"]
          output_path <- paste0(HOMER_RESULTS, folder_name)
          
          if (!file.exists(output_path)) {
            dir.create(output_path)
          }
          
          if (!file.exists(paste0(output_path, "/knownResults.txt"))) {
            system(paste("findMotifsGenome.pl",
                       filename, unlist(json[[input$organism]][8]),
                       paste0(output_path, "/"),
                 "-size 200 -mask"), intern = FALSE)
          } else {
            homer_motifs <- read.table(paste0(output_path, "/knownResults.txt"),
                                       sep = "\t", header = FALSE, skip = 1)
          }

          print(homer_motifs)
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
  })
}






# nolint end

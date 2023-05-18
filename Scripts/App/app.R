# nolint start
library(pacman)
p_load(shiny, data.table, rtracklayer, ggplot2, ggthemes, plyranges, ggpubr,
       BRGenomics, reshape2, plotly, heatmaply, dplyr, gplots, genomation,
       Biostrings, scales, GenomicRanges, DT, shinythemes, shinycustomloader,
       ggseqlogo, ChIPseeker, tools, reactable, annotables, enrichplot,
       clusterProfiler, shinyalert, rjson, BSgenome.Mmusculus.UCSC.mm10,
       TxDb.Mmusculus.UCSC.mm10.knownGene, org.Mm.eg.db)

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
  # genome <- reactive({
  #   result <- fromJSON(file = paste0(INPUT, "genome_data.json"))
  #   # print(input$organism)
  #   # print(result[[1]][5])
  #   # for (genome in 1:length(result)) {
  #   #   if (input$organism == "Nenurodyta") {
  #   #     # shinyalert("Nurodykite organizmą!", type = "error",
  #   #     #            confirmButtonText = "Rinktis organizmą")
  #   #     return()
  #   #   } else if (input$organism == result[[genome]][5]) {
  #   #     if (system.file(package = result[[genome]][1], character.only = TRUE) != 0 &&
  #   #         system.file(package = result[[genome]][2], character.only = TRUE) != 0 &&
  #   #         system.file(package = result[[genome]][3], character.only = TRUE) != 0) {

  #   #           library(pacman)
  #   #       p_load(result[[genome]][1], result[[genome]][2], result[[genome]][3], character.only = TRUE)
  #   #         } else {
  #   #           BiocManager::install(result[[genome]][1], update = FALSE)
  #   #           BiocManager::install(result[[genome]][2], update = FALSE)
  #   #           BiocManager::install(result[[genome]][3], update = FALSE)
  #   #         }
  #   #     return(result)
  #   #   } else { next }
  #   # }
  # })

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
      
    output$samples <- DT::renderDataTable({converted_table}, server = TRUE)
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

         genome <- input$organism
        unlisted_json <- unlist(json[[genome]][6])
        start_end <- unlist(strsplit(unlisted_json[1], ":"))

        chromosomes <-
          c(rep(paste0("chr", seq(start_end[1], start_end[2]))),
            rep(paste0("chr", unlisted_json[2:(length(unlisted_json))])))

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
        mm_known_genes <- TxDb.Mmusculus.UCSC.mm10.knownGene
        peak_annotations <- lapply(grl, annotatePeak, TxDb = mm_known_genes,
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
        mm_known_genes <- TxDb.Mmusculus.UCSC.mm10.knownGene
        peak_annotations <- lapply(grl, annotatePeak, TxDb = mm_known_genes,
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
      print("labas")
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
          coord_cartesian(ylim = c(0, as.numeric(max(melted_df$value)) + 20000)) +
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
    output$table1 <- renderReactable({find_ontologies_table(go_data(), "BP")})
    output$table2 <- renderReactable({find_ontologies_table(go_data(), "MF")})
    output$table3 <- renderReactable({find_ontologies_table(go_data(), "CC")})

    # GENOMIC DATA ANALYSIS - GO ANALYSIS (ACYCLIC GRAPH) (PLOT 9):
    output$plot9 <- renderPlot({find_ontologies_graph(go_data(), "BP")})
    output$plot10 <- renderPlot({find_ontologies_graph(go_data(), "MF")})
    output$plot11 <- renderPlot({find_ontologies_graph(go_data(), "CC")})

    # GENOMIC DATA ANALYSIS - GO ANALYSIS (TREE POT) (PLOT 10):
    output$plot12 <- renderPlot({find_ontologies_tree(go_data(), "BP")})
    output$plot13 <- renderPlot({find_ontologies_tree(go_data(), "MF")})
    output$plot14 <- renderPlot({find_ontologies_tree(go_data(), "CC")})

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
        selectInput(
          inputId = "tf_options",
          label = "Pasirinkite transkripcijos faktorių:",
          choices = c(
            "Tbx5", "GATA3", "GATA4", "Tbx5", "Tcf21", "CTCF", "FOXA2",
            "Nenurodyta"),
          selected = "Nenurodyta"
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
        tableOutput("table0"),
        # p("Transkripcijos faktoriaus motyvo logotipas:"),
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
              withLoader(plotOutput("plot1"), type = "html", loader = "dnaspin"),
              downloadButton("download1", "Atsisiųsti vaizdą")
            )
          ),
          tabPanel("Pikų skaičius chromosomose",
            p("Pateiktose stulpelinėse diagramose pavaizduota, kaip pikų skaičius
              pasiskirstęs skirtingose chromosomose:"),
            shinydashboard::box(
              width = 12,
              withLoader(plotOutput("plot2"), type = "html", loader = "dnaspin"),
              downloadButton("download2", "Atsisiųsti vaizdą")
            )
          ),
          tabPanel("Mėginių panašumas",
            p("Pateiktame spalvų intensyvumo grafike pavaizduota, kokia pikų
              dalis (procentiškai) sutampa tarp skirtingų mėginių:"),
            shinydashboard::box(
              width = 12,
              withLoader(plotOutput("plot3"), type = "html", loader = "dnaspin"),
              downloadButton("download3", "Atsisiųsti vaizdą")
            )
          ),
          tabPanel("Genominė distribucija",
            p("Pateiktame grafike pavaizduota pasirinktų mėginių genominių
              elementų procentinė dalis:"),
            shinydashboard::box(
              width = 12,
              withLoader(plotOutput("plot4"), type = "html", loader = "dnaspin"),
              downloadButton("download4", "Atsisiųsti vaizdą")
            )
          ),
          tabPanel("Atstumas iki TSS",
            p("Pateiktame grafike pavaizduota pasirinktų mėginių anotuotų
              pikų atstumai iki TSS (angl. Transcription Start Site):"),
            shinydashboard::box(
              width = 12,
              withLoader(plotOutput("plot5"), type = "html", loader = "dnaspin"),
              downloadButton("download5", "Atsisiųsti vaizdą")
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
              withLoader(plotOutput("plot8"), type = "html", loader = "dnaspin"),
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
  






shinyApp(ui = ui, server = server)

# nolint end
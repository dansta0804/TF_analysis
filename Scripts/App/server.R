# nolint start
library(dplyr)
library(data.table)
library(rtracklayer)
library(randomcoloR)
library(ggplot2)
library(scales)
library(GenomicRanges)
options(shiny.maxRequestSize=30*1024^2)

source("/home/daniele/Desktop/IV_course/II_semester/TF_analysis/Scripts/functions.R")


server <- function(input, output, session) {
  ####################### 1st PLOT FROM TBX5_ANALYSIS_I ####################### <- WORKS
  # output$table <- renderPlot({

  #   req(input$bigbed)
  #   bigbed_files <- list()

  #     for(sample in 1:length(input$bigbed[, 1])) {
  #       bigbed_files[[sample]] <-
  #         read.table(file = input$bigbed[[sample, 'datapath']])
  #       names(bigbed_files)[sample] <-
  #         substring(input$bigbed[[sample, 'name']],1, 11)
  #     }

  #   peaks <- data.frame(matrix(ncol = 2, nrow = 0))
  #   colnames(peaks) <- c("Experiment", "Peak_count")

  #   for (file in 1:length(bigbed_files)) {
  #       nrows <- length(rownames(bigbed_files[[file]]))
  #       row <- c(names(bigbed_files[file]), nrows)
  #       peaks[nrow(peaks) + 1, ] <- row
  #   }

  #   ggplot(peaks, aes(x = Experiment, y = as.numeric(Peak_count))) +
  #     geom_bar(stat = "identity", position = "dodge", width = 0.5,
  #               color = "black", fill = "#930d1f") +
  #     labs(x = "", y = "Pikų skaičius", size = 5) +
  #     ylim(0, 200000) +
  #     scale_y_continuous(labels = label_number(suffix = " K", scale = 1e-3)) +
  #     theme(panel.background = element_rect(fill = "#eeeef1",
  #                                           colour = "#4c0001"),
  #           panel.grid.major.y = element_line(colour = "#cab5b5", size = 0.3,
  #                                             linetype = "dashed"),
  #           panel.grid.minor.y = element_line(colour = "#cab5b5", size = 0.3,
  #                                             linetype = "dashed"),
  #           panel.grid.major.x = element_line(colour = "#cab5b5", size = 0.2,
  #                                             linetype = "longdash"),
  #           panel.grid.minor.x = element_line(colour = "#cab5b5", size = 0.2,
  #                                             linetype = "longdash"),
  #           axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
  #                                       size = 11, face = "bold",
  #                                       color = "black"),
  #           axis.text.y = element_text(size = 11, face = "bold",
  #                                       color = "black"),
  #           axis.title.x = element_text(size = 2),
  #           axis.title.y = element_text(size = 16),
  #           plot.title = element_text(hjust = 0.5, face = "bold"))
  #   })

    ####################### 3rd PLOT FROM TBX5_ANALYSIS_I ####################### <- WORKS
    # output$plot2 <- renderPlot({
    #   req(input$bigbed)
    #   bigbed_files <- list()

    #   for(sample in 1:length(input$bigbed[, 1])) {
    #     bigbed_files[[sample]] <-
    #       read.table(file = input$bigbed[[sample, 'datapath']])
    #     names(bigbed_files)[sample] <-
    #       substring(input$bigbed[[sample, 'name']],1, 11)
    #     colnames(bigbed_files[[sample]]) <-
    #       c("chrom", "start", "end", "name", rep(paste0("Other", 1:7)))
    #     bigbed_files[[sample]] <-
    #       makeGRangesFromDataFrame(bigbed_files[[sample]],
    #                                keep.extra.columns = TRUE)
    #   }

    #   coef_matrix <- matrix(nrow = length(bigbed_files),
    #                         ncol = length(bigbed_files))

    #   # Calculating Jaccard coefficient for sample pair:
    #   for (i in 1:length(bigbed_files)) {
    #       for (y in 1:length(bigbed_files)) {
    #           coef_matrix[i, y] = jaccard(bigbed_files, i, y)
    #       }
    #   }

    #   # Setting colnames and rownames for the matrix:
    #   colnames(coef_matrix) <- names(bigbed_files)
    #   rownames(coef_matrix) <- names(bigbed_files)

    #   coef_mat1 <- coef_matrix
    #   coef_mat2 <- coef_matrix

    #   # Passing Jaccard coefficients to matrix except for the diagonal - it
    #   # contains 'NA':
    #   coef_mat1[lower.tri(coef_mat1, diag = TRUE)] <- NA
    #   coef_mat2[upper.tri(coef_mat2, diag = TRUE)] <- NA

    #   # Binding two matrixes:
    #   coef_mat_bind <- rbind(data.matrix(coef_mat1), data.matrix(coef_mat2))

    #   # Translating matrix to dataframe using melt() function:
    #   melt_coef_mat <- melt(coef_mat_bind, na.rm = TRUE)

    #   # Creating a heatmap that shows similarity between samples:
    #   ggplot(melt_coef_mat, aes(x = Var2, y = Var1, fill = value)) +
    #     geom_tile(color = "black") +
    #     geom_text(aes(label = round(value, digits = 3)), size = 4.5,
    #                   color = "#030101", fontface = "bold") +
    #     labs(x = "", y = "") +
    #     scale_fill_gradient(low = "#ffee8e", high = "#ab1f1f") +
    #     guides(fill = guide_colourbar(title = "Koeficientas", face = "bold")) +
    #     theme(axis.text = element_text(size = 12, colour = "black",
    #                                    face = "bold"),
    #           axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    #           axis.title.x = element_text(size = 14, colour = "black"),
    #           axis.title.y = element_text(size = 14, colour = "black"),
    #           panel.grid.major = element_line(color = "#eeeeee"),
    #           plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    #           legend.position = "bottom")
    # })

    ####################### 4th PLOT FROM TBX5_ANALYSIS_I #######################
    output$plot3 <- renderPlot({
      req(input$bigbed)
      bigbed_files <- list()

      for(sample in 1:length(input$bigbed[, 1])) {
        bigbed_files[[sample]] <-
          read.table(file = input$bigbed[[sample, 'datapath']])
        names(bigbed_files)[sample] <-
          substring(input$bigbed[[sample, 'name']],1, 11)
        colnames(bigbed_files[[sample]]) <-
          c("chrom", "start", "end", "name", rep(paste0("Other", 1:7)))
        bigbed_files[[sample]] <-
          makeGRangesFromDataFrame(bigbed_files[[sample]],
                                   keep.extra.columns = TRUE)
      }

      # Čia turės būti naudojama nebūtinai šita matrica. Reikia pateikti visų
      # galimų TF PWM matricų sąrašą?
      mpwm <- read.table(paste0(INTERMEDIATE_FILES, "TBX5_MOUSE.H11MO.0.D.pwm"))
      mpwm <- t(mpwm)

      # Setting matrix rownames:
      rownames(mpwm) <- c("A", "C", "G", "T")

      # Generating Tbx5 sequence logo and saving the image:
      tbx5_motif <- ggseqlogo(mpwm)

      # # Calculating Tbx5 hits within fetched fasta files of peaks:
      # peak_sequences <- list.files(path = paste0(INPUTS, "FASTA/"), "*fasta")

      # hit_vec <- c()

      # # Creating an empty dataframe:
      # tbx5_hits <- data.frame(matrix(ncol = 3, nrow = 0)) %>%
      #                 setNames(., c("Sample", "Peak_count", "Tbx5_hits"))

      # # Declaring a function that calculates Tbx5 motif hits in each sample:
      # find_motif_hits <- function(sequences) {
      #     for (i in 1:length(sequences)) {
      #         hits <- countPWM(as.matrix(mpwm), sequences[[i]], min.score = "75%")
      #         if (hits == 0) { next }
      #         else { hit_vec <- c(hit_vec, hits)}
      #     }
      #     return(sum(hit_vec))
      # }

      # # Declaring a function that calculates total peak count in each sample:
      # calculate_peaks <- function(filename) {
      #     name <- paste0(tools::file_path_sans_ext(filename), ".bb")
      #     file <- import(paste0(INPUTS, "BigBed/", name)) %>%
      #                     filter(seqnames %in% chr_abr)
      #     region_count <- length(file)
      #     return(region_count)
      # }

      # # Creating a new dataframe that stores data:
      # #   - generated sample name (used in plots);
      # #   - Tbx5 motif hit count;
      # #   - Total peak count;
      # #   - Tbx5 motif percentage.

      # # TEMPORARY COMMENTED CHUNK
      # # columns <- c("Sample", "Motif_count", "Peak_count", "Percentage")
      # # tbx5_motifs <- data.frame(matrix(nrow = 0, ncol = length(columns)))

      # # colnames(tbx5_motifs) <- columns

      # # for (i in 1:length(peak_sequences)) {
      # #     filename <- peak_sequences[i]
      # #     bb_file <- paste0(tools::file_path_sans_ext(filename), ".bb")
      # #     samples_plot_name <- samples$Graph_names[samples$Filename == bb_file]

      # #     seq <- readDNAStringSet(paste0(pth_fasta, filename))
      # #     motif <- find_motif_hits(seq)
      # #     peaks <- calculate_peaks(filename)
      # #     percentage <- round((motif / peaks) * 100, 2)

      # #     data_row <- c(samples_plot_name, motif, peaks, paste0(percentage, "%"))
      # #     tbx5_motifs[nrow(tbx5_motifs) + 1, ] <- data_row
      # # }

      # # write.csv(tbx5_motif, paste0(INTERMEDIATE_FILES, "tbx5_motif_data.csv"))
      # # THE END OF THE COMMENT

      # # Reading a file that stores information about Tbx5 motif counts and
      # # peak percentages:
      # motif_data <- read.csv(file = paste0(INTERMEDIATE_FILES, "tbx5_motif_data.csv"))

      # # Calling factor() function in order to maintain certain Sample order:
      # motif_data$Sample <- factor(motif_data$Sample, levels = motif_data$Sample)

      # # Subseting data to extract all columns except for 'X' column (column 1):
      # subset_df <- motif_data[, 2:5]

      # # 'Melting' the dataframe:
      # melted_df <- melt(subset_df)

      # # Creating a barplot that visualizes Tbx5 motif distribution within samples:
      # plot4 <- ggplot(data = melted_df, aes(x = Sample, y = value,
      #                                       fill = variable, label = value)) +
      #             geom_bar(stat = "identity", colour = "#35170450", size = 0.5,
      #                     width = 0.8) +
      #             scale_fill_manual(values = c("#e3a15e", "#c7633b"),
      #                               labels = c("Tbx5 motyvų skaičius",
      #                                         "Pikų skaičius")) +
      #             scale_y_continuous(labels = label_number(suffix = " K",
      #                                                     scale = 1e-3)) +
      #             guides(fill = guide_legend(title = "Spalvų paaiškinimas",
      #                                       size = 6)) +
      #             labs(title = "", x = "", y = "TF/Pikų skaičius") +
      #             theme(axis.text = element_text(size = 10, colour = "black"),
      #                   axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
      #                                             size = 12, face = "bold"),
      #                   axis.text.y = element_text(size = 12, face = "bold"),
      #                   axis.title.x = element_text(size = 14, colour = "black"),
      #                   axis.title.y = element_text(size = 14, colour = "black"),
      #                   panel.grid.major = element_line(color = "#eeeeee"),
      #                   plot.title = element_text(hjust = 0.5, size = 16,
      #                                             face = "bold"),
      #                   panel.background = element_rect(fill = "#eeeef1",
      #                                                   colour = "#4c0001"),
      #                   panel.grid.major.y = element_line(colour = "#cab5b5",
      #                                             size = 0.3, linetype = "dashed"),
      #                   panel.grid.minor.y = element_line(colour = "#cab5b5",
      #                                             size = 0.3, linetype = "dashed"),
      #                   panel.grid.major.x = element_line(colour = "#cab5b5",
      #                                             linetype = "longdash", size = 0.2),
      #                   panel.grid.minor.x = element_line(colour = "#cab5b5",
      #                                             linetype = "longdash", size = 0.2),
      #                   legend.position = c(0.79, 0.9),
      #                   legend.title = element_text(size = 12),
      #                   legend.text = element_text(size = 11))








    })
}
# nolint end

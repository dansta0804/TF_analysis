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
    })
}
# nolint end
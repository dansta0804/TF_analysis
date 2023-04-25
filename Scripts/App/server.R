# nolint start

library(pacman)
p_load(dplyr, data.table, rtracklayer, randomcoloR, ggplot2, scales,
       GenomicRanges, ggseqlogo, BSgenome.Mmusculus.UCSC.mm10)

options(scipen = 100)
options(shiny.maxRequestSize = 30 * 1024 ^ 2)

PROJECT <- "/home/daniele/Desktop/IV_course/II_semester/TF_analysis/"
FUNCTIONS <- paste0(PROJECT, "Scripts/functions.R")
source(FUNCTIONS)
source(paste0(PROJECT, "Scripts/App/data_quality.R"))

server <- function(input, output, session) {
#   output$plot1 <- renderPlot({
#     req(input$bigbed)
#     bigbed_files <- list()

#     for(sample in 1:length(input$bigbed[, 1])) {
#       bigbed_files[[sample]] <-
#         read.table(file = input$bigbed[[sample, 'datapath']])
#       names(bigbed_files)[sample] <-
#         substring(input$bigbed[[sample, 'name']],1, 11)
#     }

#     peaks <- data.frame(matrix(ncol = 2, nrow = 0))
#     colnames(peaks) <- c("Experiment", "Peak_count")

#     for (file in 1:length(bigbed_files)) {
#         nrows <- length(rownames(bigbed_files[[file]]))
#         row <- c(names(bigbed_files[file]), nrows)
#         peaks[nrow(peaks) + 1, ] <- row
#     }

#     ggplot(peaks, aes(x = Experiment, y = as.numeric(Peak_count))) +
#       geom_bar(stat = "identity", position = "dodge", width = 0.5,
#                 color = "black", fill = "#930d1f") +
#       labs(x = "", y = "Pikų skaičius", size = 5) +
#       ylim(0, 200000) +
#       scale_y_continuous(labels = label_number(suffix = " K", scale = 1e-3)) +
#       theme(panel.background = element_rect(fill = "#eeeef1",
#                                             colour = "#4c0001"),
#             panel.grid.major.y = element_line(colour = "#cab5b5", size = 0.3,
#                                               linetype = "dashed"),
#             panel.grid.minor.y = element_line(colour = "#cab5b5", size = 0.3,
#                                               linetype = "dashed"),
#             panel.grid.major.x = element_line(colour = "#cab5b5", size = 0.2,
#                                               linetype = "longdash"),
#             panel.grid.minor.x = element_line(colour = "#cab5b5", size = 0.2,
#                                               linetype = "longdash"),
#             axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
#                                         size = 11, face = "bold",
#                                         color = "black"),
#             axis.text.y = element_text(size = 11, face = "bold",
#                                         color = "black"),
#             axis.title.x = element_text(size = 2),
#             axis.title.y = element_text(size = 16),
#             plot.title = element_text(hjust = 0.5, face = "bold"))
# })

# ####################### 2nd PLOT FROM TBX5_ANALYSIS_I ######################
# output$plot2 <- renderPlot({
#     req(input$bigbed)
#     bigbed_files <- list()
#     grl <- GRangesList()

#     for(sample in 1:length(input$bigbed[, 1])) {
#     bigbed_files[[sample]] <-
#         read.table(file = input$bigbed[[sample, 'datapath']])
    
#     colnames(bigbed_files[[sample]]) <-
#         c("chrom", "start", "end", "name", rep(paste0("Other", 1:7)))
#     grl[[sample]] <-
#         makeGRangesFromDataFrame(bigbed_files[[sample]],
#                                 keep.extra.columns = TRUE)
#     names(grl)[sample] <-
#         substring(input$bigbed[[sample, 'name']], 1, 11)
#     }

#     peak_counts <-
#         lapply(names(grl), count_peaks, objects = grl) %>%
#         bind_rows()

#     # Calling factor() function in order to maintain certain Chromosome
#     # and Name order:
#     unique_chr <- unique(peak_counts$Chromosome)
#     peak_counts$Chromosome <-
#         factor(peak_counts$Chromosome, levels = c(unique_chr))
#     peak_counts$Name <-
#         factor(peak_counts$Name, levels = unique(peak_counts$Name))

#     # Creating barplots that visualize peak differences between different
#     # chromosomes:
#     ggplot(peak_counts, aes(x = Name, y = as.numeric(Peak_count))) +
#     geom_bar(stat = "identity", color = "black", fill = "#930d1f") +
#     ylab("Pikų skaičius") +
#     facet_wrap(~ Chromosome, ncol = 7) +
#     xlab("") +
#     scale_y_continuous(labels = label_number(suffix = " K", scale = 1e-3)) +
#     theme_linedraw() +
#     theme(axis.text.x = element_text(angle = 90, size = 12, vjust = 0.5),
#         legend.position = "none",
#         axis.text.y = element_text(size = 14, face = "bold"),
#         axis.title.x = element_text(size = 14, colour = "black"),
#         axis.title.y = element_text(size = 20, colour = "black"),
#         strip.background = element_rect(fill = "white"),
#         strip.text = element_text(colour = "black", face = "bold", size = 16))
# })

# ####################### 3rd PLOT FROM TBX5_ANALYSIS_I ######################
# output$plot3 <- renderPlot({
#     req(input$bigbed)
#     bigbed_files <- list()

#     for(sample in 1:length(input$bigbed[, 1])) {
#     bigbed_files[[sample]] <-
#         read.table(file = input$bigbed[[sample, 'datapath']])
#     names(bigbed_files)[sample] <-
#         substring(input$bigbed[[sample, 'name']],1, 11)
#     colnames(bigbed_files[[sample]]) <-
#         c("chrom", "start", "end", "name", rep(paste0("Other", 1:7)))
#     bigbed_files[[sample]] <-
#         makeGRangesFromDataFrame(bigbed_files[[sample]],
#                                 keep.extra.columns = TRUE)
#     }

#     coef_matrix <- matrix(nrow = length(bigbed_files),
#                         ncol = length(bigbed_files))

#     # Calculating Jaccard coefficient for sample pair:
#     for (i in 1:length(bigbed_files)) {
#         for (y in 1:length(bigbed_files)) {
#             coef_matrix[i, y] = jaccard(bigbed_files, i, y)
#         }
#     }

#     # Setting colnames and rownames for the matrix:
#     colnames(coef_matrix) <- names(bigbed_files)
#     rownames(coef_matrix) <- names(bigbed_files)

#     coef_mat1 <- coef_matrix
#     coef_mat2 <- coef_matrix

#     # Passing Jaccard coefficients to matrix except for the diagonal - it
#     # contains 'NA':
#     coef_mat1[lower.tri(coef_mat1, diag = TRUE)] <- NA
#     coef_mat2[upper.tri(coef_mat2, diag = TRUE)] <- NA

#     # Binding two matrixes:
#     coef_mat_bind <- rbind(data.matrix(coef_mat1), data.matrix(coef_mat2))

#     # Translating matrix to dataframe using melt() function:
#     melt_coef_mat <- melt(coef_mat_bind, na.rm = TRUE)

#     # Creating a heatmap that shows similarity between samples:
#     ggplot(melt_coef_mat, aes(x = Var2, y = Var1, fill = value)) +
#     geom_tile(color = "black") +
#     geom_text(aes(label = round(value, digits = 3)), size = 4.5,
#                   color = "#030101", fontface = "bold") +
#     labs(x = "", y = "") +
#     scale_fill_gradient(low = "#ffee8e", high = "#ab1f1f") +
#     guides(fill = guide_colourbar(title = "Koeficientas", face = "bold")) +
#     theme(axis.text = element_text(size = 12, colour = "black",
#                                     face = "bold"),
#         axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#         axis.title.x = element_text(size = 14, colour = "black"),
#         axis.title.y = element_text(size = 14, colour = "black"),
#         panel.grid.major = element_line(color = "#eeeeee"),
#         plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
#         legend.position = "bottom")
# })

# output$plot4 <- renderPlot({
#     req(input$pwm)
#     mpwm <- read.table(file = input$pwm$datapath)
#     mpwm <- t(mpwm)

#     # Setting matrix rownames:
#     rownames(mpwm) <- c("A", "C", "G", "T")
#     ggseqlogo(mpwm)
# })

# ####################### 5th PLOT FROM TBX5_ANALYSIS_I ######################
# output$plot5 <- renderPlot({
#     req(input$bigbed)
#     req(input$pwm)
#     bigbed_files <- list()

#     for(sample in 1:length(input$bigbed[, 1])) {
#         bigbed_files[[sample]] <-
#             read.table(file = input$bigbed[[sample, 'datapath']])
        
#         colnames(bigbed_files[[sample]]) <-
#             c("chrom", "start", "end", "name", rep(paste0("Other", 1:7)))
#         bigbed_files[[sample]] <-
#             makeGRangesFromDataFrame(bigbed_files[[sample]],
#                                     keep.extra.columns = TRUE)
#         names(bigbed_files)[sample] <-
#             input$bigbed[[sample, 'name']]
#     }

#     # Čia turės būti naudojama nebūtinai šita matrica. Reikia pateikti visų
#     # galimų TF PWM matricų sąrašą?
#     mpwm <- read.table(file = input$pwm$datapath)
#     mpwm <- t(mpwm)

#     # Setting matrix rownames:
#     rownames(mpwm) <- c("A", "C", "G", "T")

#     peak_sequences <- list()

#     for (file in 1:length(bigbed_files)) {
#         peak_sequences[[file]] <-
#             getSeq(BSgenome.Mmusculus.UCSC.mm10, bigbed_files[[file]])
#         names(peak_sequences)[file] <- names(bigbed_files)[file]
#     }

#     tbx5_motifs <- data.frame(matrix(nrow = 0, ncol = 4))
#     colnames(tbx5_motifs) <-
#         c("Sample", "Motif_count", "Peak_count", "Percentage")

#     for (i in 1:length(peak_sequences)) {
#         filename <- input$bigbed[[i, 'datapath']]
#         name <- substring(names(peak_sequences[i]), 1, 11)

#         motif <- find_motif_hits(peak_sequences[[i]], mpwm)
#         peaks <- calculate_peaks(filename)
#         percentage <- round((motif / peaks) * 100, 2)

#         data_row <- c(name, motif, peaks, paste0(percentage, "%"))
#         tbx5_motifs[nrow(tbx5_motifs) + 1, ] <- data_row
#     }

#     # Reading a file that stores information about Tbx5 motif counts and
#     # peak percentages:
#     motif_data <- tbx5_motifs

#     # Calling factor() function in order to maintain certain Sample order:
#     motif_data$Sample <- factor(motif_data$Sample, levels = motif_data$Sample)

#     # Subseting data to extract all columns except for 'X' column (column 1):
#     subset_df <- motif_data[, 1:3]

#     # 'Melting' the dataframe:
#     melted_df <- melt(subset_df, id = c("Sample"))

#     ggplot(data = melted_df, aes(x = Sample, y = as.numeric(value),
#                                 fill = variable, label = value)) +
#         geom_bar(stat = "identity", colour = "#35170450", size = 0.5,
#                  width = 0.8) +
#         scale_fill_manual(values = c("#e3a15e", "#c7633b"),
#                           labels = c("Tbx5 motyvų skaičius", "Pikų skaičius")) +
#         scale_y_continuous(labels = label_number(suffix = " K", scale = 1e-3)) +
#         guides(fill = guide_legend(title = "Spalvų paaiškinimas", size = 6)) +
#         labs(title = "", x = "", y = "TF/Pikų skaičius") +
#         theme(axis.text = element_text(size = 10, colour = "black"),
#             axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
#                                        size = 12, face = "bold"),
#             axis.text.y = element_text(size = 12, face = "bold"),
#             axis.title.x = element_text(size = 14, colour = "black"),
#             axis.title.y = element_text(size = 14, colour = "black"),
#             panel.grid.major = element_line(color = "#eeeeee"),
#             plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
#             panel.background = element_rect(fill = "#eeeef1",
#                                             colour = "#4c0001"),
#             panel.grid.major.y = element_line(colour = "#cab5b5", size = 0.3,
#                                               linetype = "dashed"),
#             panel.grid.minor.y = element_line(colour = "#cab5b5", size = 0.3,
#                                               linetype = "dashed"),
#             panel.grid.major.x = element_line(colour = "#cab5b5", size = 0.2,
#                                               linetype = "longdash"),
#             panel.grid.minor.x = element_line(colour = "#cab5b5", size = 0.2,
#                                               linetype = "longdash"),
#             legend.position = c(0.79, 0.9),
#             legend.title = element_text(size = 12),
#             legend.text = element_text(size = 11))
# })

# ####################### 6th PLOT FROM TBX5_ANALYSIS_I ######################
# output$table1 = DT::renderDataTable({
#     req(input$bigbed)
#     bigbed_files <- list()
#     names <- c()

#     for(sample in 1:length(input$bigbed[, 1])) {
#     names <- c(names, input$bigbed[[sample, 'name']])
#       bigbed_files[[sample]] <-
#         read.table(file = input$bigbed[[sample, 'datapath']])
#       names(bigbed_files)[sample] <-
#         substring(input$bigbed[[sample, 'name']],1, 11)
#     }

# system(paste("findMotifsGenome.pl", paste0(INPUTS, "BED/", names), "mm10",
#        paste0(PROJECT, "results2/"), "-size 200 -mask"), intern = FALSE)

#     homer_motifs <- read.table(paste0(PROJECT, "results2/knownResults.txt"),
#                                sep = "\t", header = FALSE, skip = 1)



#     colnames(homer_motifs) <-
#         c("Motif Name", "Consensus", "P value",
#             "Log P value", "q value (Benjamini)",
#             "Target Sequences with Motif (#)",
#             "Target Sequences with Motif (%)",
#             "Background Sequences with Motif (#)",
#             "Background Sequences with Motif (%)")

#     for (motif in 1:nrow(homer_motifs)) {
#         homer_motifs$`Motif Name`[motif] <-
#             strsplit(homer_motifs$`Motif Name`, "/")[[motif]][1]
#     }
#     homer_motifs
# })

# output$downloadData <- downloadHandler(
#     filename = function() {
#         paste("De_novo_motifs_", Sys.Date(), ".csv")
#     },
#     content = function(file) {
#         homer_motifs <- read.table(paste0(PROJECT, "results/knownResults.txt"),
#                                    sep = "\t", header = FALSE, skip = 1)

#         colnames(homer_motifs) <-
#             c("Motif Name", "Consensus", "P value",
#             "Log P value", "q value (Benjamini)",
#             "Target Sequences with Motif (#)",
#             "Target Sequences with Motif (%)",
#             "Background Sequences with Motif (#)",
#             "Background Sequences with Motif (%)")

#         for (motif in 1:nrow(homer_motifs)) {
#             homer_motifs$`Motif Name`[motif] <-
#                 strsplit(homer_motifs$`Motif Name`, "/")[[motif]][1]
#         }
#         write.csv(homer_motifs, file, row.names = FALSE)
#     }
# )

# Galbūt čia reikia pridėti mygtuką, kurį paspaudus pasirinkti mėginiai yra
# apjungiami į vieną ir tada atliekamos analizės tik su apjungu duomenų rinkiniu?

    output$samples <-
        DT::renderDataTable({input$bigbed[, c("name", "size")]}, server = TRUE)

    observe({
        req(input$samples_rows_selected)
        selRow <- input$bigbed[input$samples_rows_selected,]
        print(selRow[[1]])
        print(selRow[[4]])
        names <- selRow[[1]]
        samples <- selRow[[4]]

        bigbed_files <- list()
        output$plot1 <- renderPlot({
            for(sample in 1:length(samples)) {
                bigbed_files[[sample]] <-
                    read.table(file = samples[sample])
                names(bigbed_files)[sample] <-
                    substring(names[sample], 1, 11)
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
                scale_y_continuous(labels = label_number(suffix = " K",
                                                         scale = 1e-3)) +
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
    })
}
# nolint end

# nolint start

PROJECT <- "/home/daniele/Desktop/IV_course/II_semester/Baigiamasis_darbas/TF_analysis/"

server <- function(input, output, session) {
######################## 4th PLOT PWM SEQUENCE LOGO ############################
output$plot4 <- renderPlot({
    req(input$pwm)
    mpwm <- read.table(file = input$pwm$datapath)
    mpwm <- t(mpwm)

    # Setting matrix rownames:
    rownames(mpwm) <- c("A", "C", "G", "T")
    ggseqlogo(mpwm)
})

####################### 5th PLOT FROM TBX5_ANALYSIS_I ######################
output$plot5 <- renderPlot({
    req(input$bigbed)
    req(input$pwm)
    bigbed_files <- list()

    for(sample in 1:length(input$bigbed[, 1])) {
        bigbed_files[[sample]] <-
            read.table(file = input$bigbed[[sample, 'datapath']])
        
        colnames(bigbed_files[[sample]]) <-
            c("chrom", "start", "end", "name", rep(paste0("Other", 1:7)))
        bigbed_files[[sample]] <-
            makeGRangesFromDataFrame(bigbed_files[[sample]],
                                    keep.extra.columns = TRUE)
        names(bigbed_files)[sample] <-
            input$bigbed[[sample, 'name']]
    }

    # Čia turės būti naudojama nebūtinai šita matrica. Reikia pateikti visų
    # galimų TF PWM matricų sąrašą?
    mpwm <- read.table(file = input$pwm$datapath)
    mpwm <- t(mpwm)

    # Setting matrix rownames:
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

    # Calling factor() function in order to maintain certain Sample order:
    motif_data$Sample <- factor(motif_data$Sample, levels = motif_data$Sample)

    # Subseting data to extract all columns except for 'X' column (column 1):
    subset_df <- motif_data[, 1:3]

    # 'Melting' the dataframe:
    melted_df <- melt(subset_df, id = c("Sample"))

    ggplot(data = melted_df, aes(x = Sample, y = as.numeric(value),
                                fill = variable, label = value)) +
        geom_bar(stat = "identity", colour = "#35170450", size = 0.5,
                 width = 0.8) +
        scale_fill_manual(values = c("#e3a15e", "#c7633b"),
                          labels = c("Tbx5 motyvų skaičius", "Pikų skaičius")) +
        scale_y_continuous(labels = label_number(suffix = " K", scale = 1e-3)) +
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
            panel.grid.major.y = element_line(colour = "#cab5b5", size = 0.3,
                                              linetype = "dashed"),
            panel.grid.minor.y = element_line(colour = "#cab5b5", size = 0.3,
                                              linetype = "dashed"),
            panel.grid.major.x = element_line(colour = "#cab5b5", size = 0.2,
                                              linetype = "longdash"),
            panel.grid.minor.x = element_line(colour = "#cab5b5", size = 0.2,
                                              linetype = "longdash"),
            legend.position = c(0.79, 0.9),
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 11))
})

####################### 6th PLOT FROM TBX5_ANALYSIS_I ######################
output$table1 = DT::renderDataTable({
    # req(input$bigbed)
    # bigbed_files <- list()

    # for(sample in 1:length(input$bigbed[, 1])) {
    #   bigbed_files[[sample]] <-
    #     read.table(file = input$bigbed[[sample, 'datapath']])
    #   names(bigbed_files)[sample] <-
    #     substring(input$bigbed[[sample, 'name']],1, 11)
    # }

# system(paste0("findMotifsGenome.pl ~/Desktop/IV_course/II_semester/TF_analysis/Analyses/Tbx5_analysis_I/Inputs/BED/", bigbed_files[3], " mm10 ~/Desktop/IV_course/II_semester/TF_analysis/results/ -size 200 -mask"), intern = FALSE)

    homer_motifs <- read.table(paste0(PROJECT, "results/knownResults.txt"),
                               sep = "\t", header = FALSE, skip = 1)

    colnames(homer_motifs) <-
        c("Motif Name", "Consensus", "P value",
            "Log P value", "q value (Benjamini)",
            "Target Sequences with Motif (#)",
            "Target Sequences with Motif (%)",
            "Background Sequences with Motif (#)",
            "Background Sequences with Motif (%)")

    for (motif in 1:nrow(homer_motifs)) {
        homer_motifs$`Motif Name`[motif] <-
            strsplit(homer_motifs$`Motif Name`, "/")[[motif]][1]
    }
    homer_motifs
})

output$downloadData <- downloadHandler(
    filename = function() {
        paste("De_novo_motifs_", Sys.Date(), ".csv")
    },
    content = function(file) {
        homer_motifs <- read.table(paste0(PROJECT, "results/knownResults.txt"),
                                   sep = "\t", header = FALSE, skip = 1)

        colnames(homer_motifs) <-
            c("Motif Name", "Consensus", "P value",
            "Log P value", "q value (Benjamini)",
            "Target Sequences with Motif (#)",
            "Target Sequences with Motif (%)",
            "Background Sequences with Motif (#)",
            "Background Sequences with Motif (%)")

        for (motif in 1:nrow(homer_motifs)) {
            homer_motifs$`Motif Name`[motif] <-
                strsplit(homer_motifs$`Motif Name`, "/")[[motif]][1]
        }
        write.csv(homer_motifs, file, row.names = FALSE)
    }
)
}
# nolint end
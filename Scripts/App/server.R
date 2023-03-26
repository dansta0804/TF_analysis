# nolint start
library(dplyr)
library(data.table)
library(rtracklayer)
library(randomcoloR)
library(ggplot2)
library(scales)
options(shiny.maxRequestSize=30*1024^2)


server <- function(input, output, session) {
  output$table <- renderPlot({
    req(input$bigbed)
    bigbed_files <- list()

    for(sample in 1:length(input$bigbed[, 1])) {
      bigbed_files[[sample]] <-
        read.table(file = input$bigbed[[sample, 'datapath']])
      names(bigbed_files)[sample] <-
        substring(input$bigbed[[sample, 'name']],1, 11)
    }

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
            panel.grid.major.y = element_line(colour = "#cab5b5", size = 0.3,
                                              linetype = "dashed"),
            panel.grid.minor.y = element_line(colour = "#cab5b5", size = 0.3,
                                              linetype = "dashed"),
            panel.grid.major.x = element_line(colour = "#cab5b5", size = 0.2,
                                              linetype = "longdash"),
            panel.grid.minor.x = element_line(colour = "#cab5b5", size = 0.2,
                                              linetype = "longdash"),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                                        size = 11, face = "bold",
                                        color = "black"),
            axis.text.y = element_text(size = 11, face = "bold",
                                        color = "black"),
            axis.title.x = element_text(size = 2),
            axis.title.y = element_text(size = 16),
            plot.title = element_text(hjust = 0.5, face = "bold"))
    })
}
# nolint end
# nolint start

show_gene_counts <- function(data_frame) {
    plot <-
        ggplot(data_frame, aes(x = Sample, y = as.numeric(Gene_count))) +
        geom_bar(stat = "identity", fill = "#cc8b12", width = .6,
                 color = "#8d5c00") +
        labs(title = "", x = "",
             y = expression(paste("Sutampančių genų skaičius ",
                            italic("M. musculus")))) +
        geom_text(aes(label = Gene_count), color = "#030101",
                  size = 5,vjust = -1, fontface = 2) +
        theme(axis.text = element_text(size = 10, colour = "black"),
                axis.text.x =
                    element_text(angle = 45, hjust = 1, size = 12, vjust = 1,
                                 face = "bold"),
                axis.text.y = element_text(size = 12, face = "bold"),
                axis.title.x = element_text(size = 14, colour = "black"),
                axis.title.y = element_text(size = 14, colour = "black"),
                panel.grid.major = element_line(color = "#eeeeee"),
                plot.title =
                    element_text(hjust = 0.5, size = 16, face = "bold"),
                panel.background =
                    element_rect(fill = "#eeeef1", colour = "#3a1010"),
                panel.grid.major.y =
                    element_line(colour = "#cab5b5", size = 0.3,
                                 linetype = "dashed"),
                panel.grid.minor.y =
                    element_line(colour = "#cab5b5", size = 0.3,
                                 linetype = "dashed"),
                panel.grid.major.x =
                    element_line(colour = "#cab5b5", size = 0.2,
                                 linetype = "longdash"),
                panel.grid.minor.x =
                    element_line(colour = "#cab5b5", size = 0.2,
                                 linetype = "longdash"),
                legend.position = c(0.777, 0.788))
return(plot)
}


# nolint end
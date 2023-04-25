# nolint start

library(pacman)
p_load(shiny, data.table, rtracklayer, ggplot2, ggthemes, plyranges, ggpubr,
       BRGenomics, reshape2, plotly, heatmaply, dplyr, gplots, genomation,
       Biostrings, scales, GenomicRanges, DT, shinythemes, shinycustomloader)

PROJECT <- "/home/daniele/Desktop/IV_course/II_semester/TF_analysis/"
source(paste0("/home/daniele/Desktop/IV_course/II_semester/",
              "TF_analysis/Scripts/App/ui.R"))
source(paste0(PROJECT, "Scripts/App/server.R"))
shinyApp(ui = ui, server = server)

# nolint end
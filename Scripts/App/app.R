# nolint start
library(pacman)
p_load(shiny, data.table, rtracklayer, ggplot2, ggthemes, plyranges, ggpubr,
       BRGenomics, reshape2, plotly, heatmaply, dplyr, gplots, genomation,
       Biostrings, scales, GenomicRanges, DT, shinythemes, shinycustomloader,
       ggseqlogo, ChIPseeker, tools, reactable, annotables, enrichplot,
       clusterProfiler, shinyalert, rjson)

PROJECT <- "./"
source(paste0(PROJECT, "Scripts/App/ui.R"))
source(paste0(PROJECT, "Scripts/App/server.R"))

shinyApp(ui = ui, server = server)

# nolint end
# nolint start
library(pacman)

# Library loading:
p_load(shiny, data.table, rtracklayer, ggplot2, ggthemes, plyranges, ggpubr,
       BRGenomics, reshape2, plotly, heatmaply, dplyr, gplots, genomation,
       Biostrings, scales, GenomicRanges, DT, shinythemes, shinycustomloader,
       ggseqlogo, ChIPseeker, tools, reactable, annotables, enrichplot,
       clusterProfiler, shinyalert, rjson, shiny, shinythemes, shinydashboard)

library(BiocManager)
options(repos = BiocManager::repositories())

# Path declarations:
PROJECT <- "./"
source(paste0(PROJECT, "Scripts/App/ui.R"))
source(paste0(PROJECT, "Scripts/App/server.R"))

# A function that starts an app:
shinyApp(ui = ui, server = server)
# nolint end
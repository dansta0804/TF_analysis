
library(pacman)
p_load(GenomicRanges, ggplot2, rtracklayer, dplyr, plyranges, ChIPseeker,
TxDb.Mmusculus.UCSC.mm10.knownGene, TxDb.Hsapiens.UCSC.hg38.knownGene,
BSgenome.Mmusculus.UCSC.mm10, BSgenome.Hsapiens.UCSC.hg38, org.Mm.eg.db,
org.Hs.eg.db, annotables, rBLAST, rentrez)

# nolint start
PROJECT     <- "/home/daniele/Desktop/IV_course/II_semester/TF_analysis/"
BIGBEDS     <- paste0(PROJECT, "Input/")
MMUSCULUS   <- paste0(BIGBEDS, "Mus_musculus/")
HSAPIENS    <- paste0(BIGBEDS, "Homo_sapiens/")
INTER_FILES <- paste0(PROJECT, "Intermediate_data/")
SCRIPTS     <- paste0(PROJECT, "Scripts/")
RESULTS     <- paste0(INTER_FILES, "Generated_files")
FASTAS      <- paste0(INTER_FILES, "FASTA/")
FIGURES     <- paste0(PROJECT, "Figures/")

SAMPLES     <- read.csv(file = paste0(INTER_FILES, "sample_key.csv"))


# nolint start
# Listing sample BigBed files for each analyzed organism:
bbfiles_mm <- list.files(path = MMUSCULUS, "*bb")
bbfiles_hg <- list.files(path = HSAPIENS, "*bb")

grl_mm <- GRangesList()
for (object in 1:length(grl)) {
    name <- gsub(" ", "_", names(grl[object]))
    peak <- grl[[object]]
    seqlengths(peak) <- seqlengths(peak) - 10001
    peak_annotation <- annotatePeak(peak, tssRegion = c(-3000, 3000),
                         TxDb = mm_known_genes, annoDb = "org.Mm.eg.db")

    mm_annot <- as.data.frame(peak_annotation@anno)
    entrezids <- unique(mm_annot$geneId)
    entrez2gene <- grcm38 %>% filter(entrez %in% entrezids) %>%
                                dplyr::select(entrez, symbol)

    m <- match(mm_annot$geneId, entrez2gene$entrez)
    mm_annot <- cbind(mm_annot[, 1:14], gene_symbol = entrez2gene$symbol[m],
                      mm_annot[, 15:16])

    # Defining the same grl_smlr object that has two extra columns with
    # gene id and gene symbol:
    grl_mm[[object]] <- mm_annot
    names(grl_mm)[object] <- names(grl)[object]
}
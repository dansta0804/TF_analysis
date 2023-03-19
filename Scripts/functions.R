# nolint start
# A function that reads PWM of certain transcription factor:
get_PWM <- function(pwm_name) {
    TF_pwm <-
        read.table(file = paste0(INTER_FILES, pwm_name)) %>%
        t() %>%
        `rownames<-`(c("A", "C", "G", "T"))
    return(TF_pwm)
}

# A function that annotates peaks by assigning gene id and gene symbol:
annotate_peaks <- function(peak_set, genome, known_genes, anno_db) {
    grl_annotation <- GRangesList()
    for (object in 1:length(peak_set)) {
        name <- gsub(" ", "_", names(peak_set[object]))
        peak <- peak_set[[object]]
        seqlengths(peak) <- seqlengths(peak) - 10001
        peak_annotation <-
            annotatePeak(peak, tssRegion = c(-3000, 3000), TxDb = known_genes,
            annoDb = anno_db)

        annotation <- as.data.frame(peak_annotation@anno)
        entrezids <- unique(annotation$geneId)
        entrez2gene <-
            genome %>%
            filter(entrez %in% entrezids) %>%
            dplyr::select(entrez, symbol)

        m <- match(annotation$geneId, entrez2gene$entrez)
        organism_annot <-
            cbind(annotation[, 1:14], gene_symbol = entrez2gene$symbol[m],
                  annotation[, 15:16])

        # Defining an object that has two extra columns with gene id and symbol:
        grl_annotation[[object]] <- organism_annot
        names(grl_annotation)[object] <- names(peak_set)[object]
    }
    return(grl_annotation)
}






# nolint end
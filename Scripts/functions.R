# nolint start
PROJECT     <- "/home/daniele/Desktop/IV_course/II_semester/TF_analysis/"
INTER_FILES <- paste0(PROJECT, "Intermediate_data/")
RESULTS     <- paste0(INTER_FILES, "Generated_files")



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

count_pwm_hits <- function(sequence_dataset, pwm) {
    # This is a list that stores pwm hit counts for each sample:
    gene_hits_total <- list()

    # Predicting the number of transcription factor binding hits in Homo
    # sapiens by using PWM matrix:
    for (sample in 1:length(sequence_dataset)) {
        gene_hits <- data.frame(matrix(ncol = 3, nrow = 0))
        colnames(gene_hits) <- c("Gene", "Hits", "Length")
        for (gene in 1:length(sequence_dataset[[sample]])) {
            pwm_hits <- countPWM(as.matrix(pwm),
                                as.character(sequence_dataset[[sample]][gene]),
                                min.score = "75%")
            row <- c(names(sequence_dataset[[sample]][gene]), pwm_hits,
                    length(sequence_dataset[[sample]][[gene]]))
            gene_hits[nrow(gene_hits) + 1, ] <- row
        }
        gene_hits_total[[sample]] <- gene_hits

        # Writing the results into a file:
        write.table(gene_hits_total[[sample]],
                    file = paste0(RESULTS, "/Hg_pwm_hits.txt"),
                    sep = "\t", quote = FALSE, row.names = FALSE)
    }
    return(gene_hits_total)
}






# nolint end
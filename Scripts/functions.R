# nolint start
library(imager)

PROJECT     <- "./"
INTER_FILES <- paste0(PROJECT, "Intermediate_data/")
RESULTS     <- paste0(INTER_FILES, "Generated_files/")
FIGURES     <- paste0(PROJECT, "Figures/")

# A function that formats number by adding spaces (e.g., 5000 -> 5 000):
spaces <- function(number) {
  format(number, big.mark = " ")
}

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
    # seqlengths(peak) <- seqlengths(peak) - 10001
    peak_annotation <-
      annotatePeak(peak, tssRegion = c(-3000, 3000),
      TxDb = known_genes, annoDb = anno_db)

    annotation <- as.data.frame(peak_annotation@anno)
    entrezids <- unique(annotation$geneId)
    entrez2gene <-
      genome %>%
      dplyr::filter(entrez %in% entrezids) %>%
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

# A function that is used to create a sequence dataset be retrieving specified
# sequences from certain genome:
get_sequences <- function(genome, genes_of_interest) {
  sequence_dataset <- c()

  # Creating a sequence dataset that stores certain genome gene sequences:
  for (gene in 1:length(genes_of_interest)) {
    gene_sequences <- getSeq(genome, genes_of_interest[[gene]])
    names(gene_sequences) <-
      as.data.frame(genes_of_interest[[gene]])$gene_symbol

    # Creating a sequence vector for all the samples:
    sequence_dataset <- c(sequence_dataset, gene_sequences)
  }
  return(sequence_dataset)
}

# State is one of the two: ("predicted", "control")
count_pwm_hits <- function(sequence_dataset, pwm, state) {
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
    # print(gene_hits)

    # Writing the results into a file:
    write.table(gene_hits_total[[sample]],
                file = paste0(RESULTS, paste0("Hg_pwm_hits_", state, "_",
                sample, ".txt")), sep = "\t", quote = FALSE,
                row.names = FALSE)
  }
  return(gene_hits_total)
}

###################### FUNCTIONS FROM TBX5 ANALYSIS PART I #####################
# A function that calculates how many peaks are in each chromosome:
count_peaks <- function(name, objects, chr_abbr, chr_lengths) {
  # Defining chromosome abbreviations and lenghts (Mbp):
  chr_abr <- chr_abbr
  chr_mbp <- chr_lengths

  # Calculating chromosome lenghts in bp:
  chr_bp <- chr_mbp * 1000000

  # Creating a dataframe that stores data about peak counts in each
  # chromosome for different samples:
  peak_counts <-
    data.frame(matrix(ncol = 3, nrow = 0)) %>%
    setNames(., c("Name", "Chromosome", "Peak_count"))

  for (chr in 1:length(chr_abr)) {
    peaks <-
      objects[[name]] %>%
      filter(seqnames == chr_abr[chr]) %>%
      length() / chr_bp[chr] * 100000000
    peak_counts[nrow(peak_counts) + 1, ] <- c(name, chr_abr[chr], peaks)
  }
  return(peak_counts)
}

# Declaring a function that calculates modified Jaccard coefficient:
jaccard <- function(granges, a, b) {
  len <-
    reduce(c(granges[[a]], granges[[a]])) %>%
    length()

  return((length(GenomicRanges::intersect(granges[[a]],
                                          granges[[b]])) / len) * 100)
}

find_motif_hits <- function(sequences, mpwm) {
  hit_vec <- c()
  for (i in 1:length(sequences)) {
    hits <- countPWM(as.matrix(mpwm), sequences[[i]], min.score = "75%")
    if (hits == 0) { next }
    else {hit_vec <- c(hit_vec, hits)}
  }
  return(sum(hit_vec))
}

# A function that calculates total peak count in each sample:
calculate_peaks <- function(filename) {
  file <- read.table(filename)
  region_count <- length(rownames(file))
  return(region_count)
}

html <- function(x, inline = FALSE) {
  container <- if (inline) htmltools::span else htmltools::div
  container(dangerouslySetInnerHTML = list("__html" = x))
}

find_ontologies_table <- function(data, genome, subontology) {
  pl <- enrichGO(gene = data[[1]]$ENTREZID, OrgDb = get(genome),
                 ont = subontology, pAdjustMethod = "BH", pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05, readable = TRUE)

  if (length(rownames(as.data.frame(pl))) == 0) {
    error <-
      data.frame(`Klaida` = paste0("GO analizės rezultate nebuvo gauti genų ",
                                   "GO subontologijų apibūdinimai!"))
    reactable(
      error, searchable = FALSE, showSortable = FALSE, rownames = FALSE,
      pagination = FALSE, highlight = FALSE,
      defaultColDef = colDef(
        align = "center",
        minWidth = 70
      ),
      columns = list(
        `Klaida` = colDef(
          style = function(value) {
            list(color = "#8f2222", fontWeight = "bold", fontSize = "20pt")
          }
        )
      )
    )
  } else {
    pl <- as.data.frame(pl)
    pl$Status <- "Peržiūrėti genų sąrašą"
    pl <-
      pl %>%
      dplyr::select(c("ID", "Description", "GeneRatio", "Count",
                      "Status", "geneID")) %>%
      rename("GO ID" = "ID", "Apibūdinimas" = "Description",
              "Genų santykis" = "GeneRatio", "Genų skaičius" = "Count",
              "Peržiūra" = "Status")
    
    reactable(
      pl, searchable = FALSE, showSortable = TRUE, rownames = FALSE,
      pagination = TRUE, highlight = TRUE,
      defaultColDef = colDef(
        align = "center",
        minWidth = 70
      ),
      columns = list(
        `Peržiūra` = colDef(
          style = function(value) {
            list(color = "#500909", fontWeight = "bold")
          },
          details = function(value) {
            genes <- unlist(strsplit(pl$geneID, split = "/")[value])
            add_NA <- 8 - (length(genes) - 8 * (length(genes) %/% 8))
            genes <- c(genes, rep(" ", add_NA))
            tb <- data.frame(matrix(genes, ncol = 8), ncol = 8)
            colnames(tb) <- c(rep(paste0("Genai ", 1:8)))

            htmltools::div(
              style = "padding: 1rem",
              reactable(tb[, 1:8], outlined = TRUE, fullWidth = TRUE,
                defaultColDef = colDef(
                  align = "center",
                  minWidth = 70
                )
              )
            )
          }
        ),
        geneID = colDef(show = FALSE)
      )
    )
  }
}

find_ontologies_graph <- function(data, genome, subontology) {
  pl <- enrichGO(gene = data[[1]]$ENTREZID, OrgDb = get(genome),
                 ont = subontology, pAdjustMethod = "BH", pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05, readable = TRUE)
  if (length(rownames(as.data.frame(pl))) == 0) {
    plot(load.image(paste0(FIGURES, "Error_message.png")))
  } else {
    goplot(pl)
  }
}

find_ontologies_tree <- function(data, genome, subontology) {
  pl <- enrichGO(gene = data[[1]]$ENTREZID, OrgDb = get(genome),
                 ont = subontology, pAdjustMethod = "BH", pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05, readable = TRUE)
  if (length(rownames(as.data.frame(pl))) == 0) {
    plot(load.image(paste0(FIGURES, "Error_message.png")))
  } else {
    pl_modified <- setReadable(pl, 'org.Mm.eg.db', 'ENTREZID')
    pl_modified <- pairwise_termsim(pl_modified)
    treeplot(pl_modified, cluster.params = list(method = "average"),
            xlim = c(0, 30))
  }
}
# nolint end

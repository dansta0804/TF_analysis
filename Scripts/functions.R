# nolint start
# Library loading:
library(pacman)
p_load(imager, GenomicRanges, dplyr)

# Path declarations:
PROJECT     <- "./"
FIGURES     <- paste0(PROJECT, "Figures/")

########################### FUNCTION FOR SAMPLE TABLE ##########################
# A function that formats number by adding spaces (e.g., 5000 -> 5 000):
spaces <- function(number) {
  format(number, big.mark = " ")
}

########################## FUNCTIONS FOR DATA QUALITY ##########################
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
      dplyr::filter(seqnames == chr_abr[chr]) %>%
      length() / chr_bp[chr] * 100000000
    peak_counts[nrow(peak_counts) + 1, ] <- c(name, chr_abr[chr], peaks)
  }
  return(peak_counts)
}

# A function that calculates modified Jaccard coefficient:
jaccard <- function(granges, a, b) {
  len <-
    reduce(c(granges[[a]], granges[[a]])) %>%
    length()

  return(
    (length(GenomicRanges::intersect(granges[[a]], granges[[b]])) / len) * 100
  )
}

# A function that calculates positional weight matrix hits:
find_motif_hits <- function(sequences, mpwm) {
  hit_vec <- c()
  for (i in 1:length(sequences)) {
    hits <- countPWM(as.matrix(mpwm), sequences[[i]], min.score = "75%")
    if (hits == 0) { next }
    else { hit_vec <- c(hit_vec, hits) }
  }
  return(sum(hit_vec))
}

# A function that calculates total peak count in each sample:
calculate_peaks <- function(filename) {
  file <- read.table(filename)
  region_count <- length(rownames(file))
  return(region_count)
}

#################### FUNCTIONS FOR BIOLOGICAL DATA ANALYSIS ####################
# A function that reads PWM of certain transcription factor:
get_PWM <- function(pwm_name) {
  TF_pwm <-
    read.table(file = pwm_name, skip = 1) %>%
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
    peak_annotation <-
      annotatePeak(
        peak,
        tssRegion = c(-3000, 3000),
        TxDb = known_genes,
        annoDb = anno_db
      )

    annotation <- as.data.frame(peak_annotation@anno)
    entrezids <- unique(annotation$geneId)
    entrez2gene <-
      genome %>%
      dplyr::filter(entrez %in% entrezids) %>%
      dplyr::select(entrez, symbol)

    m <- match(annotation$geneId, entrez2gene$entrez)
    organism_annot <-
      cbind(
        annotation[, 1:14],
        gene_symbol = entrez2gene$symbol[m],
        annotation[, 15:16]
      )

    # Defining an object that has two extra columns with gene id and symbol:
    grl_annotation[[object]] <- organism_annot
    names(grl_annotation)[object] <- names(peak_set)[object]
  }
  return(grl_annotation)
}

# A function that creates a table after performing GO analysis for the
# specified subontology:
find_ontologies_table <- function(data, genome, subontology) {
  pl <-
    enrichGO(
      gene = data[[1]]$ENTREZID,
      OrgDb = get(genome),
      ont = subontology,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.01,
      qvalueCutoff  = 0.05,
      readable = TRUE
    )

  if (length(rownames(as.data.frame(pl))) == 0) {
    error <-
      data.frame(`Klaida` = paste0("GO analizės rezultate nebuvo gauti genų ",
                                   "GO subontologijų apibūdinimai!"))
    reactable(
      error,
      searchable = FALSE,
      showSortable = FALSE,
      rownames = FALSE,
      pagination = FALSE,
      highlight = FALSE,
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
      dplyr::select(
        c("ID", "Description", "GeneRatio", "Count", "Status", "geneID")
      ) %>%
      rename(
        "GO ID" = "ID",
        "Apibūdinimas" = "Description",
        "Genų santykis" = "GeneRatio",
        "Genų skaičius" = "Count",
        "Peržiūra" = "Status"
      )
    
    reactable(
      pl,
      searchable = FALSE,
      showSortable = TRUE,
      rownames = FALSE,
      pagination = TRUE,
      highlight = TRUE,
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
              reactable(
                tb[, 1:8],
                outlined = TRUE,
                fullWidth = TRUE,
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

# A function that plots a graph after performing GO analysis for the
# specified subontology:
find_ontologies_graph <- function(data, genome, subontology) {
  pl <-
    enrichGO(
      gene = data[[1]]$ENTREZID,
      OrgDb = get(genome),
      ont = subontology,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.01,
      qvalueCutoff  = 0.05,
      readable = TRUE
    )
  if (length(rownames(as.data.frame(pl))) == 0) {
    plot(load.image(paste0(FIGURES, "Error_message.png")))
  } else {
    goplot(pl)
  }
}

# A function that creates a treeplot with highlighted clusters after
# performing GO analysis for the specified subontology:
find_ontologies_tree <- function(data, genome, subontology) {
  pl <-
    enrichGO(
      gene = data[[1]]$ENTREZID,
      OrgDb = get(genome),
      ont = subontology,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.01,
      qvalueCutoff  = 0.05,
      readable = TRUE
    )
  if (length(rownames(as.data.frame(pl))) == 0) {
    plot(load.image(paste0(FIGURES, "Error_message.png")))
  } else {
    pl_modified <- setReadable(pl, 'org.Mm.eg.db', 'ENTREZID')
    pl_modified <- pairwise_termsim(pl_modified)
    treeplot(
      pl_modified,
      cluster.params = list(method = "average"),
      xlim = c(0, 30)
    )
  }
}

###################### FUNCTIONS FOR TF TARGET PREDICTION ######################
# A function that counts PWM hits in query and subject sequences:
find_PWM_hits <- function (q_seq_list, s_seq_list, pwm, min_score) {
  percentages_qs <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(percentages_qs) <- c("Sample", "Percentage")

  for (seq in 1:length(q_seq_list)) {
    percentages_query <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(percentages_query) <- c("Gene", "QueryPWMHits")

    percentages_subject <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(percentages_subject) <- c("Gene", "SubjectPWMHits")

    for (gene in 1:length(q_seq_list[[seq]])) {
      q_hits <-
        as.numeric(countPWM(as.matrix(pwm),
                            as.character(q_seq_list[[seq]][gene]),
                            min.score = min_score))

      row <- c(names(q_seq_list[[seq]][gene]), q_hits)
      percentages_query[nrow(percentages_query) + 1, ] <- row            
    }

    # write.csv(percentages_query,
    #           paste0(RESULTS, paste0("PWM_hits/percent_mm_", seq), ".csv"),
    #           row.names = FALSE)

    percentages_query$Gene <- toupper(percentages_query$Gene)

    for (gene in 1:length(s_seq_list[[seq]])) {
      s_hits <-
        as.numeric(countPWM(as.matrix(pwm),
                            as.character(s_seq_list[[seq]][gene]),
                            min.score = min_score))

      row <- c(names(s_seq_list[[seq]][gene]), hg_hits)
      percentages_subject[nrow(percentages_subject) + 1, ] <- row            
    }

    # write.csv(percentages_subject,
    #           paste0(RESULTS, paste0("PWM_hits/percent_hg_", seq), ".csv"),
    #           row.names = FALSE)

    percent_merge <- merge(percentages_query, percentages_subject)
    percent_merge["Percentage"] <- ""

    print(head(percent_merge))

    for (row in 1:(length(rownames(percent_merge)))) {
      query_hits <- percent_merge[row, "QueryPWMHits"]
      subject_hits <- percent_merge[row, "SubjectPWMHits"]
      percentage <- (as.numeric(query_hits) / as.numeric(subject_hits)) * 100

      if (percentage %in% c("Inf", "NaN", 0)) { next }
      else { percent_merge[row, "Percentage"] <- as.numeric(percentage) }    
    }

    percent_merges[[seq]] <- percent_merge
    names(percent_merges)[seq] <- paste0("Sample", seq)

    # write.csv(percent_merge, paste0(RESULTS, paste0("percentages_", seq), ".csv"),
    #         row.names = FALSE)

    percentages_qs[nrow(percentages_qs) + 1, ] <-
        c(paste("Sample", seq),
          mean(as.numeric(percent_merge$Percentage), na.rm = TRUE))
}}
# nolint end

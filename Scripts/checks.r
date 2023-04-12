# nolint start

library(pacman)
p_load(dplyr, data.table, rtracklayer, randomcoloR, ggplot2, scales,
       GenomicRanges, ggseqlogo, BSgenome.Mmusculus.UCSC.mm10)

source("/home/daniele/Desktop/IV_course/II_semester/TF_analysis/Scripts/functions.R")
INPUTS <- "/home/daniele/Desktop/IV_course/II_semester/TF_analysis/Analyses/Tbx5_analysis_I/Inputs/"
INTERMEDIATE_FILES <- "/home/daniele/Desktop/IV_course/II_semester/TF_analysis/Analyses/Tbx5_analysis_I/Intermediate_data_files/"

bed_files <- list.files(path = paste0(INPUTS, "BED/"), "*bed")
grl <- list()

for (i in 1:length(bed_files)) {
    grl[[i]] <- read.table(paste0(INPUTS, "BED/", bed_files[i]))
    names(grl)[i] <- bed_files[i]
}

# Čia turės būti naudojama nebūtinai šita matrica. Reikia pateikti visų
# galimų TF PWM matricų sąrašą?
mpwm <- read.table(paste0(INTERMEDIATE_FILES, "TBX5_MOUSE.H11MO.0.D.pwm"))
mpwm <- t(mpwm)

# Setting matrix rownames:
rownames(mpwm) <- c("A", "C", "G", "T")

# Generating Tbx5 sequence logo and saving the image:
tbx5_motif <- ggseqlogo(mpwm)


colnames(grl[[1]])
colnames(grl[[1]]) <- c("chrom", "start", "end", "name", rep(paste0("Other", 1:6)))

new_list <- GRangesList()

for (file in 1:length(grl)) {
    colnames(grl[[file]]) <- c("chrom", "start", "end", "name", rep(paste0("Other", 1:6)))
    new_list[[file]] <- makeGRangesFromDataFrame(grl[[file]], keep.extra.columns = TRUE)
    names(new_list)[file] <- paste0("FILE", file)
}
a <- makeGRangesFromDataFrame(grl[[2]], keep.extra.columns = TRUE)



?? rGADEM


library(rGADEM)
library(BSgenome.Hsapiens.UCSC.hg19)
library(rtracklayer)

pwd<-"" #INPUT FILES- BedFiles, FASTA, etc.
path <- system.file("extdata/Test_100.bed",package="rGADEM")
BedFile<-paste(pwd,path,sep="")
Sequences<-import(BedFile)


gadem<-GADEM(Sequences,verbose=1,genome=Hsapiens)
getPWM(gadem)
showClass("gadem")








options(scipen = 100)

# Creating a dataframe that stores data about peak counts in each
# chromosome for different samples:
peak_counts <- data.frame(matrix(ncol = 3, nrow = 0)) %>%
                    setNames(., c("Name", "Chromosome", "Peak_count"))

# Defining chromosome lenghts (Mbp) based on MGI data:
chr_mbp <- c(195, 182, 160, 157, 152, 150, 145, 130, 124, 131,
             122, 120, 121, 125, 104, 98, 95, 91, 61, 169, 91)

# Calculating chromosome lenghts in bp:
chr_bp <- chr_mbp * 1000000

# Declaring a function that calculates how many peaks are in each chromosome:
count_peaks <- function(name, objects) {
    for (chr in 1:21) {
            peaks <- objects[[name]] %>%
                    filter(seqnames == chr_abr[chr]) %>%
                    length() / chr_bp[chr] * 100000000
            peak_counts[nrow(peak_counts) + 1,] = c(name, chr_abr[chr], peaks)
        }
        return(peak_counts)
    }

peak_counts <- lapply(names(new_list), count_peaks, objects = new_list) %>% bind_rows()




cas <- unique(peak_counts$Chromosome)
peak_counts$Chromosome <- factor(peak_counts$Chromosome, levels = c(cas))
peak_counts$Name <- factor(peak_counts$Name, levels = unique(peak_counts$Name))

class(as.numeric(peak_counts$Peak_count))
# Note: it is essential to add 'as.numeric'. Otherwise, Peak_count is
# considered as 'string' type - column values are represented in
# lexicographic order which leads to wrong visualization.

# Creating barplots that visualize peak differences between different
# chromosomes:
ggplot(peak_counts, aes(x = Name, y = as.numeric(Peak_count))) +
            geom_bar(stat = "identity", color = "black", fill = "red") +
            ylab("Pikų skaičius") +
            facet_wrap(~ Chromosome, ncol = 7) +
            xlab("") +
            scale_y_continuous(labels = label_number(suffix = " K",
                               scale = 1e-3)) +
            theme_linedraw() +
            theme(axis.text.x = element_text(angle = 90, size = 12,
                  vjust = 0.5),
                  legend.position = "none",
                  axis.text.y = element_text(size = 14, face = "bold"),
                  axis.title.x = element_text(size = 14, colour = "black"),
                  axis.title.y = element_text(size = 20, colour = "black"),
                  strip.background = element_rect(fill = "white"),
                  strip.text = element_text(colour = "black", face = "bold",
                                            size = 16))










library(BSgenome.Mmusculus.UCSC.mm10)
chr_abr <- c(paste0("chr", 1:19), "chrX", "chrY")

seq = getSeq(BSgenome.Mmusculus.UCSC.mm10, a)

new_list[1]
peak_sequences <- list()

for (file in 1:length(new_list)) {
        peak_sequences[[file]] <- getSeq(BSgenome.Mmusculus.UCSC.mm10, new_list[[file]])
        names(peak_sequences)[file] <- names(new_list)[file]
      }

peak_sequences[1]

find_motif_hits <- function(sequences) {
    hit_vec <- c()
          for (i in 1:length(sequences)) {
              hits <- countPWM(as.matrix(mpwm), sequences[[i]], min.score = "75%")
              if (hits == 0) { next }
              else { hit_vec <- c(hit_vec, hits)}
          }
          return(sum(hit_vec))
      }

# # Declaring a function that calculates total peak count in each sample:
calculate_peaks <- function() {
    name <- "PEAKS054910_Tbx5_P70326_MACS2_3053.bed"
    file <- read.table(paste0(INPUTS, "BED/", name))
    region_count <- length(rownames(file))
    return(region_count)
}


# TEMPORARY COMMENTED CHUNK
columns <- c("Sample", "Motif_count", "Peak_count", "Percentage")
tbx5_motifs <- data.frame(matrix(nrow = 0, ncol = length(columns)))

colnames(tbx5_motifs) <- columns

samples_plot_name <- "NAME2"

motif <- find_motif_hits(seq)
peaks <- calculate_peaks()
percentage <- round((motif / peaks) * 100, 2)

data_row <- c(samples_plot_name, motif, peaks, paste0(percentage, "%"))
tbx5_motifs[nrow(tbx5_motifs) + 1, ] <- data_row





























peak_sequences <- BSgenome::getSeq(BSgenome.Mmusculus.UCSC.mm39, grl[[4]])
seqlengths(grl[[1]]) <- seqlengths(grl[[2]]) - 10001
peak_sequences <- Biostrings::writeXStringSet(peak_sequences , "my.fasta")

hit_vec <- c()

# Creating an empty dataframe:
tbx5_hits <- data.frame(matrix(ncol = 3, nrow = 0)) %>%
                setNames(., c("Sample", "Peak_count", "Tbx5_hits"))

# Declaring a function that calculates Tbx5 motif hits in each sample:
find_motif_hits <- function(sequences) {
    for (i in 1:length(sequences)) {
        hits <- countPWM(as.matrix(mpwm), sequences[[i]], min.score = "75%")
        if (hits == 0) { next }
        else { hit_vec <- c(hit_vec, hits)}
    }
    return(sum(hit_vec))
}

# # Declaring a function that calculates total peak count in each sample:
calculate_peaks <- function(filename) {
    name <- paste0(tools::file_path_sans_ext(filename), ".bb")
    file <- import(paste0(INPUTS, "BigBed/", name)) %>%
                    filter(seqnames %in% chr_abr)
    region_count <- length(file)
    return(region_count)
}

# Creating a new dataframe that stores data:
#   - generated sample name (used in plots);
#   - Tbx5 motif hit count;
#   - Total peak count;
#   - Tbx5 motif percentage.

# TEMPORARY COMMENTED CHUNK
# columns <- c("Sample", "Motif_count", "Peak_count", "Percentage")
# tbx5_motifs <- data.frame(matrix(nrow = 0, ncol = length(columns)))

# colnames(tbx5_motifs) <- columns

# for (i in 1:length(peak_sequences)) {
#     filename <- peak_sequences[i]
#     bb_file <- paste0(tools::file_path_sans_ext(filename), ".bb")
#     samples_plot_name <- samples$Graph_names[samples$Filename == bb_file]

#     seq <- readDNAStringSet(paste0(pth_fasta, filename))
#     motif <- find_motif_hits(seq)
#     peaks <- calculate_peaks(filename)
#     percentage <- round((motif / peaks) * 100, 2)

#     data_row <- c(samples_plot_name, motif, peaks, paste0(percentage, "%"))
#     tbx5_motifs[nrow(tbx5_motifs) + 1, ] <- data_row
# }

# write.csv(tbx5_motif, paste0(INTERMEDIATE_FILES, "tbx5_motif_data.csv"))
# THE END OF THE COMMENT

# Reading a file that stores information about Tbx5 motif counts and
# peak percentages:
motif_data <- tbx5_motifs
# Calling factor() function in order to maintain certain Sample order:
motif_data$Sample <- factor(motif_data$Sample, levels = motif_data$Sample)

# Subseting data to extract all columns except for 'X' column (column 1):
subset_df <- motif_data[, 1:3]

# 'Melting' the dataframe:
melted_df <- melt(subset_df, id = c("Sample"))

# Creating a barplot that visualizes Tbx5 motif distribution within samples:
plot4 <- ggplot(data = melted_df, aes(x = Sample, y = as.numeric(value),
                                      fill = variable, label = value)) +
            geom_bar(stat = "identity", colour = "#35170450", size = 0.5,
                    width = 0.8) +
            scale_fill_manual(values = c("#e3a15e", "#c7633b"),
                              labels = c("Tbx5 motyvų skaičius",
                                        "Pikų skaičius")) +
            # scale_y_continuous(labels = label_number(suffix = " K",
            #                                         scale = 1e-3)) +
            guides(fill = guide_legend(title = "Spalvų paaiškinimas",
                                      size = 6)) +
            labs(title = "", x = "", y = "TF/Pikų skaičius") +
            theme(axis.text = element_text(size = 10, colour = "black"),
                  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                                            size = 12, face = "bold"),
                  axis.text.y = element_text(size = 12, face = "bold"),
                  axis.title.x = element_text(size = 14, colour = "black"),
                  axis.title.y = element_text(size = 14, colour = "black"),
                  panel.grid.major = element_line(color = "#eeeeee"),
                  plot.title = element_text(hjust = 0.5, size = 16,
                                            face = "bold"),
                  panel.background = element_rect(fill = "#eeeef1",
                                                  colour = "#4c0001"),
                  panel.grid.major.y = element_line(colour = "#cab5b5",
                                            size = 0.3, linetype = "dashed"),
                  panel.grid.minor.y = element_line(colour = "#cab5b5",
                                            size = 0.3, linetype = "dashed"),
                  panel.grid.major.x = element_line(colour = "#cab5b5",
                                            linetype = "longdash", size = 0.2),
                  panel.grid.minor.x = element_line(colour = "#cab5b5",
                                            linetype = "longdash", size = 0.2),
                  legend.position = c(0.79, 0.9),
                  legend.title = element_text(size = 12),
                  legend.text = element_text(size = 11))














library(BSgenome.Hsapiens.UCSC.hg19)
library(rGADEM)
library(rtracklayer)
library(GenomicRanges)

library(parseR)

? 
BedFile <- bed_files[[3]]

bed_files <- list.files(path = paste0(INPUTS, "BED/"), "*bed")

BedFilee <- import(paste0(INPUTS, "BED/", BedFile))

colnames(BedFilee) <-
          c("chrom", "start", "end", "name", rep(paste0("Other", 1:7)))
grl <- makeGRangesFromDataFrame(BedFilee, keep.extra.columns = TRUE)
Sequences <- rtracklayer::import()
gadem <- GADEM(BedFilee,verbose=1,genome=BSgenome.Mmusculus.UCSC.mm39)

BedFilee[15,]

pwd<-"" #INPUT FILES- BedFiles, FASTA, etc.
path<- system.file("extdata","Test_100.bed",package="rGADEM")
BedFile<-paste(pwd,path,sep="")



BED<-read.table(paste0(INPUTS, "BED/", bed_files[4]),header=FALSE,sep="\t")



BED<-data.frame(chr=as.factor(BED[,1]),start=as.numeric(BED[,2]),end=as.numeric(BED[,3]))
#Create RD files

library(BSgenome.Mmusculus.UCSC.mm39)
library(BSgenome.Hsapiens.UCSC.hg38)
Sequences<-makeGRangesFromDataFrame(BED, keep.extra.columns = TRUE)
gadem<-GADEM(Sequences,verbose=1,genome=BSgenome.Mmusculus.UCSC.mm39)









library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(BSgenome.Mmusculus.UCSC.mm10)


library(ChIPseeker)
mm_known_genes <- TxDb.Mmusculus.UCSC.mm10.knownGene

grll <- grl[3]

names(grll)


name <- gsub(" ", "_", names(grll[1]))
peak <- grll[[1]]
class(peak)


colnames(peak) <-
          c("chrom", "start", "end", "name", rep(paste0("Other", 1:7)))
        peakk <-
          makeGRangesFromDataFrame(peak,
                                   keep.extra.columns = TRUE)

peakk <- makeGRangesFromDataFrame(peak)

# seqlengths(peak) <- seqlengths(peak) - 10001
peak_annotation <- annotatePeak(peakk, tssRegion = c(-3000, 3000),
                        TxDb = mm_known_genes, annoDb = "org.Mm.eg.db")

as.data.frame(peak_annotation)
write(as.data.frame(peak_annotation), paste0(PROJECT,"aaa.txt"))

colnames(as.data.frame(peak_annotation))

annotation(peak_annotation)

# mm_annot <- as.data.frame(peak_annotation@anno)
# entrezids <- unique(mm_annot$geneId)
# entrez2gene <- grcm38 %>% filter(entrez %in% entrezids) %>%
#                             dplyr::select(entrez, symbol)

# m <- match(mm_annot$geneId, entrez2gene$entrez)
# mm_annot <- cbind(mm_annot[, 1:14], gene_symbol = entrez2gene$symbol[m],
#                   mm_annot[, 15:16])

# # Defining the same grl_smlr object that has two extra columns with
# # gene id and gene symbol:
# grl_mm[[object]] <- mm_annot
# names(grl_mm)[object] <- names(grl)[object]



bed_files <- list.files(path = paste0(INPUTS, "BED/"), "*bed")
fil <- read.table(paste0(INPUTS, "BED/", bed_files[3]))



system(paste0("findMotifsGenome.pl ~/Desktop/IV_course/II_semester/TF_analysis/Analyses/Tbx5_analysis_I/Inputs/BED/", bed_files[3], " mm10 ~/Desktop/IV_course/II_semester/TF_analysis/results/ -size 200 -mask"), intern = FALSE)

# nolint end
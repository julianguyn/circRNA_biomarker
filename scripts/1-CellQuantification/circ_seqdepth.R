# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggpubr)
    library(stringr)
    library(viridis)
    library(patchwork)
})

options(stringsAsFactors = FALSE)
set.seed(123)

source("utils/palettes.R")

args <- "circ"

# -----------------------------------------------------------
# Parse args
# -----------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
analysis <- args[1]

valid <- c("cpm", "unnorm")
if (is.na(analysis) || !analysis %in% valid) {
  stop(
    sprintf("Invalid analysis argument '%s'. Must be one of: %s",
            analysis, paste(valid, collapse = ", ")),
    call. = FALSE
  )
}

# get input dirs
dir <- switch(
  analysis,
  cpm = "../data/processed_cellline/common_samples/",
  unnorm = "../data/seqdepth_cellline/"
)


############################################################
# Load in circRNA expression data
############################################################

# helper function to load in data
load_counts <- function(dir, filename) {
    counts <- fread(paste0(dir, filename), data.table = FALSE)
    rownames(counts) <- counts$sample
    counts$sample <- NULL
    return(counts)
}

ciri_gcsi <- load_counts(dir, "CIRI2/ciri_gcsi_counts.tsv")
ciri_gdsc <- load_counts(dir, "CIRI2/ciri_gdsc_counts.tsv")
ciri_ccle <- load_counts(dir, "CIRI2/ciri_ccle_counts.tsv")

circ_gcsi <- load_counts(dir, "CIRCexplorer2/circ_gcsi_counts.tsv")
circ_gdsc <- load_counts(dir, "CIRCexplorer2/circ_gdsc_counts.tsv")
circ_ccle <- load_counts(dir, "CIRCexplorer2/circ_ccle_counts.tsv")

cfnd_gcsi <- load_counts(dir, "circRNA_finder/cfnd_gcsi_counts.tsv")
cfnd_gdsc <- load_counts(dir, "circRNA_finder/cfnd_gdsc_counts.tsv")
cfnd_ccle <- load_counts(dir, "circRNA_finder/cfnd_ccle_counts.tsv")

fcrc_gcsi <- load_counts(dir, "find_circ/fcrc_gcsi_counts.tsv")
fcrc_gdsc <- load_counts(dir, "find_circ/fcrc_gdsc_counts.tsv")
fcrc_ccle <- load_counts(dir, "find_circ/fcrc_ccle_counts.tsv")

############################################################
# Load in circRNA read counts
############################################################

gcsi_reads <- fread(paste0("../data/seqdepth_cellline/readCounts/gcsi.tsv"), data.table = FALSE)
ccle_reads <- fread(paste0("../data/seqdepth_cellline/readCounts/ccle.tsv"), data.table = FALSE)
gdsc_reads <- fread(paste0("../data/seqdepth_cellline/readCounts/gdsc.tsv"), data.table = FALSE)

############################################################
# Count all transcripts
############################################################

# count reads
count_reads <- function(counts, reads, label) {
    count <- rowSums(counts) |> as.data.frame()
    colnames(count) <- "Count"
    count$Seq_Depth <- reads$avg_counts[match(rownames(count), reads$sample)]
    count$Label <- label
    count$Sample <- rownames(count)
    rownames(count) <- NULL
    return(count)
}

ciri_gcsi <- count_reads(ciri_gcsi, gcsi_reads, "ciri_gcsi")
ciri_gdsc <- count_reads(ciri_gdsc, gdsc_reads, "ciri_gdsc")
ciri_ccle <- count_reads(ciri_ccle, ccle_reads, "ciri_ccle")

circ_gcsi <- count_reads(circ_gcsi, gcsi_reads, "circ_gcsi")
circ_gdsc <- count_reads(circ_gdsc, gdsc_reads, "circ_gdsc")
circ_ccle <- count_reads(circ_ccle, ccle_reads, "circ_ccle")

cfnd_gcsi <- count_reads(cfnd_gcsi, gcsi_reads, "cfnd_gcsi")
cfnd_gdsc <- count_reads(cfnd_gdsc, gdsc_reads, "cfnd_gdsc")
cfnd_ccle <- count_reads(cfnd_ccle, ccle_reads, "cfnd_ccle")

fcrc_gcsi <- count_reads(fcrc_gcsi, gcsi_reads, "fcrc_gcsi")
fcrc_gdsc <- count_reads(fcrc_gdsc, gdsc_reads, "fcrc_gdsc")
fcrc_ccle <- count_reads(fcrc_ccle, ccle_reads, "fcrc_ccle")

save(
    ciri_gcsi, ciri_gdsc, ciri_ccle,
    circ_gcsi, circ_gdsc, circ_ccle,
    cfnd_gcsi, cfnd_gdsc, cfnd_ccle,
    fcrc_gcsi, fcrc_gdsc, fcrc_ccle,
    file = "../results/data/seqdepth_temp.RData"
)

############################################################
# Create toPlot
############################################################

toPlot <- rbind(
    ciri_gcsi, ciri_gdsc, ciri_ccle,
    circ_gcsi, circ_gdsc, circ_ccle,
    cfnd_gcsi, cfnd_gdsc, cfnd_ccle,
    fcrc_gcsi, fcrc_gdsc, fcrc_ccle
)
toPlot$Pipeline <- sub("_.*", "", toPlot$Label)
toPlot$PSet <- sub(".*_", "", toPlot$Label)

toPlot$Pipeline <- factor(toPlot$Pipeline, levels = c("ciri", "circ", "cfnd", "fcrc"), labels = names(pipeline_pal)[1:4])
toPlot$PSet <- factor(toPlot$PSet, levels = c("gcsi", "ccle", "gdsc"), labels = names(pset_pal))

############################################################
# Correlations (CPM and seqdepth)
############################################################

res <- data.frame(matrix(nrow=0, ncol=3))

for (pset in unique(toPlot$PSet)) {
    pset_df <- toPlot[toPlot$PSet == pset,]
    for (pipeline in unique(pset_df$Pipeline)) {
        subset_df <- pset_df[pset_df$Pipeline == pipeline,]
        df <- data.frame(
            pset = pset,
            pipeline = pipeline,
            corr = suppressWarnings(cor.test(subset_df$Count, subset_df$Seq_Depth, method = "spearman")$estimate)
        )
        res <- rbind(res, df)
    }
}

res$pipeline <- factor(res$pipeline, levels = rev(names(pipeline_pal)[1:4]))
res$pset <- factor(res$pset, levels = names(pset_pal)[1:3])
p <- ggplot(res, aes(y = pipeline, x = pset, fill = corr)) +
    geom_tile(color = "gray") +
    geom_text(aes(label = round(corr, 2)), size = 4) +
    scale_fill_gradient2("Spearman\nCorrelation", high = "blue", low = "red", mid = "white", limits = c(-1, 1)) +
    theme_minimal() +
    theme(
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 11)
    )

filename <- paste0("results/figures/suppfig5/", analysis, "_pset_correlations.png")
ggsave(filename, p, width = 5, height = 3)



############################################################
# Load in lung seqdepth
############################################################

# load in metadata
polyA_meta <- read.table("data/rnaseq_meta/lung_polyA.tsv", header = TRUE)
ribo0_meta <- read.csv("data/rnaseq_meta/lung_ribozero.csv")

# keep common 51 (and order)
polyA_meta <- polyA_meta[match(ribo0_meta$TB_id, polyA_meta$TB_id), ]

# load in seqdepth
polyA_lib <- read.table("data/raw_lung/polyA.tsv", header = TRUE)
polyA_lib <- polyA_lib[polyA_lib$sample %in% polyA_meta$sample,]
ribo0_lib <- read.table("data/raw_lung/ribo0.tsv", header = TRUE)

############################################################
# Plot
############################################################

toPlot <- unique(toPlot[,c(2,6)])
toPlot <- rbind(
    toPlot,
    data.frame(Seq_Depth = polyA_lib$avg_counts, PSet = "polyA"),
    data.frame(Seq_Depth = ribo0_lib$avg_counts, PSet = "ribo0")
)

ggplot(toPlot, aes(x = PSet, y = Seq_Depth)) +
    geom_violin() +
    geom_boxplot(width = 0.1) +
    geom_jitter(width = 0.2, alpha = 0.7) +
    theme_bw() +
    labs(y = "Average Read Counts", x = "Dataset")

ggplot(toPlot[toPlot$PSet %in% c("polyA", "ribo0"),], aes(x = PSet, y = Seq_Depth)) +
    geom_violin() +
    geom_boxplot(width = 0.1) +
    geom_jitter(width = 0.2, alpha = 0.7) +
    theme_bw() +
    labs(y = "Average Read Counts", x = "Dataset")


# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(matrixStats)
    library(reshape2)
    library(ggplot2)
    library(ggpubr)
})

options(stringsAsFactors = FALSE)
set.seed(123)

source("utils/palettes.R")
source("utils/compute_stability_index.R")

# -----------------------------------------------------------
# Parse args
# -----------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
analysis <- args[1]

valid <- c("circ", "GE")
if (is.na(analysis) || !analysis %in% valid) {
  stop(
    sprintf("Invalid analysis argument '%s'. Must be one of: %s",
            analysis, paste(valid, collapse = ", ")),
    call. = FALSE
  )
}

############################################################
# Load in circRNA expression data
############################################################

# helper function to load and filter data
load_counts <- function(path, filename) {
    counts <- fread(paste0(path, filename), data.table = FALSE)
    # filter circRNA transcripts with low detection rates
    counts <- counts[,-which(colSums(counts == 0) > 45)]
    return(counts)
}

# load circRNA expression data
path <- switch(
    analysis,
    circ = "../data/processed_cellline/common_samples/",
    GE = "../data/processed_cellline/GE_common_samples/"
)

ciri_gcsi <- load_counts(path, "CIRI2/ciri_gcsi_counts.tsv")
ciri_gdsc <- load_counts(path, "CIRI2/ciri_gdsc_counts.tsv")
ciri_ccle <- load_counts(path, "CIRI2/ciri_ccle_counts.tsv")

circ_gcsi <- load_counts(path, "CIRCexplorer2/circ_gcsi_counts.tsv")
circ_gdsc <- load_counts(path, "CIRCexplorer2/circ_gdsc_counts.tsv")
circ_ccle <- load_counts(path, "CIRCexplorer2/circ_ccle_counts.tsv")

cfnd_gcsi <- load_counts(path, "circRNA_finder/cfnd_gcsi_counts.tsv")
cfnd_gdsc <- load_counts(path, "circRNA_finder/cfnd_gdsc_counts.tsv")
cfnd_ccle <- load_counts(path, "circRNA_finder/cfnd_ccle_counts.tsv")

fcrc_gcsi <- load_counts(path, "find_circ/fcrc_gcsi_counts.tsv")
fcrc_gdsc <- load_counts(path, "find_circ/fcrc_gdsc_counts.tsv")
fcrc_ccle <- load_counts(path, "find_circ/fcrc_ccle_counts.tsv")


############################################################
# Get common transcripts
############################################################

# get common circRNA transcripts found in all psets for each pipeline
ciri_common <- intersect(intersect(colnames(ciri_gcsi), colnames(ciri_ccle)), colnames(ciri_gdsc))
circ_common <- intersect(intersect(colnames(circ_gcsi), colnames(circ_ccle)), colnames(circ_gdsc))
cfnd_common <- intersect(intersect(colnames(cfnd_gcsi), colnames(cfnd_ccle)), colnames(cfnd_gdsc))
fcrc_common <- intersect(intersect(colnames(fcrc_gcsi), colnames(fcrc_ccle)), colnames(fcrc_gdsc))

print(length(ciri_common)) #31, 50
print(length(circ_common)) #67, 122
print(length(cfnd_common)) #85, 269
print(length(fcrc_common)) #504, 176

############################################################
# Keep only common transcripts across all PSets
############################################################

#' Helper function to order circRNA dataframes and keep
#' only transcripts found in all psets for each pipeline
subset_df <- function(df, common_transcripts) {
    df <- df[,which(colnames(df) %in% common_transcripts)]
    df <- df[order(df$sample),order(colnames(df))]
    rownames(df) <- df$sample
    df$sample <- NULL
    return(df)
}

ciri_gcsi <- subset_df(ciri_gcsi, ciri_common)
ciri_ccle <- subset_df(ciri_ccle, ciri_common)
ciri_gdsc <- subset_df(ciri_gdsc, ciri_common)

circ_gcsi <- subset_df(circ_gcsi, circ_common)
circ_ccle <- subset_df(circ_ccle, circ_common)
circ_gdsc <- subset_df(circ_gdsc, circ_common)

cfnd_gcsi <- subset_df(cfnd_gcsi, cfnd_common)
cfnd_ccle <- subset_df(cfnd_ccle, cfnd_common)
cfnd_gdsc <- subset_df(cfnd_gdsc, cfnd_common)

fcrc_gcsi <- subset_df(fcrc_gcsi, fcrc_common)
fcrc_ccle <- subset_df(fcrc_ccle, fcrc_common)
fcrc_gdsc <- subset_df(fcrc_gdsc, fcrc_common)

# save subsetted dataframes for feature importance (circ only)
#save(ciri_gcsi, ciri_ccle, ciri_gdsc, circ_gcsi, circ_ccle, circ_gdsc,
#     cfnd_gcsi, cfnd_ccle, cfnd_gdsc, fcrc_gcsi, fcrc_ccle, fcrc_gdsc,
#     file = "../results/data/temp/circ_stability_subsetdf.RData")

############################################################
# Compute pairwise spearman correlations
############################################################

print("Starting CIRI2 stability indices")
ciri_stability <- rbind(
    compute_stability_index(ciri_gcsi, ciri_ccle, "gcsi", "ccle", "ciri"),
    compute_stability_index(ciri_gcsi, ciri_gdsc, "gcsi", "gdsc", "ciri"),
    compute_stability_index(ciri_gdsc, ciri_ccle, "gdsc", "ccle", "ciri")
)

print("Starting CIRCexplorer2 stability indices")
circ_stability <- rbind(
    compute_stability_index(circ_gcsi, circ_ccle, "gcsi", "ccle", "circ"),
    compute_stability_index(circ_gcsi, circ_gdsc, "gcsi", "gdsc", "circ"),
    compute_stability_index(circ_gdsc, circ_ccle, "gdsc", "ccle", "circ")
)

print("Starting circRNA_finder stability indices")
cfnd_stability <- rbind(
    compute_stability_index(cfnd_gcsi, cfnd_ccle, "gcsi", "ccle", "cfnd"),
    compute_stability_index(cfnd_gcsi, cfnd_gdsc, "gcsi", "gdsc", "cfnd"),
    compute_stability_index(cfnd_gdsc, cfnd_ccle, "gdsc", "ccle", "cfnd")
)

print("Starting find_circ stability indices")
fcrc_stability <- rbind(
    compute_stability_index(fcrc_gcsi, fcrc_ccle, "gcsi", "ccle", "fcrc"),
    compute_stability_index(fcrc_gcsi, fcrc_gdsc, "gcsi", "gdsc", "fcrc"),
    compute_stability_index(fcrc_gdsc, fcrc_ccle, "gdsc", "ccle", "fcrc")
)

save(ciri_stability, circ_stability, cfnd_stability, fcrc_stability, 
     file = "../results/data/temp/circ_stability.RData")

############################################################
# Compute pairwise spearman correlations with randomization
############################################################

print("Starting CIRI2 random stability indices")
ciri_stability_random <- rbind(
    compute_stability_index(ciri_gcsi, ciri_ccle, "gcsi", "ccle", "ciri", random = TRUE, iter = 100),
    compute_stability_index(ciri_gcsi, ciri_gdsc, "gcsi", "gdsc", "ciri", random = TRUE, iter = 100),
    compute_stability_index(ciri_gdsc, ciri_ccle, "gdsc", "ccle", "ciri", random = TRUE, iter = 100)
)

print("Starting CIRCexplorer2 random stability indices")
circ_stability_random <- rbind(
    compute_stability_index(circ_gcsi, circ_ccle, "gcsi", "ccle", "circ", random = TRUE, iter = 100),
    compute_stability_index(circ_gcsi, circ_gdsc, "gcsi", "gdsc", "circ", random = TRUE, iter = 100),
    compute_stability_index(circ_gdsc, circ_ccle, "gdsc", "ccle", "circ", random = TRUE, iter = 100)
)

print("Starting circRNA_finder random stability indices")
cfnd_stability_random <- rbind(
    compute_stability_index(cfnd_gcsi, cfnd_ccle, "gcsi", "ccle", "cfnd", random = TRUE, iter = 100),
    compute_stability_index(cfnd_gcsi, cfnd_gdsc, "gcsi", "gdsc", "cfnd", random = TRUE, iter = 100),
    compute_stability_index(cfnd_gdsc, cfnd_ccle, "gdsc", "ccle", "cfnd", random = TRUE, iter = 100)
)

print("Starting find_circ random stability indices")
fcrc_stability_random <- rbind(
    compute_stability_index(fcrc_gcsi, fcrc_ccle, "gcsi", "ccle", "fcrc", random = TRUE, iter = 100),
    compute_stability_index(fcrc_gcsi, fcrc_gdsc, "gcsi", "gdsc", "fcrc", random = TRUE, iter = 100),
    compute_stability_index(fcrc_gdsc, fcrc_ccle, "gdsc", "ccle", "fcrc", random = TRUE, iter = 100)
)

save(ciri_stability_random, circ_stability_random, cfnd_stability_random, fcrc_stability_random, 
     file = "../results/data/temp/circ_stability_random.RData")

############################################################
# SI for isoform and gene expression data
############################################################

# only do if analysis == "cells" (do once)

if (analysis == "cells") {

    # load in data
    load("../results/data/isoform_expression.RData")
    load("../results/data/gene_expression.RData")

    print("Starting gene expression stability indices")
    gene_stability <- rbind(
        compute_stability_index(expr_gcsi_p, expr_ccle_p, "gcsi", "ccle", "ciri", random = TRUE, iter = 100),
        compute_stability_index(expr_gcsi_p, expr_gdsc_p, "gcsi", "gdsc", "ciri", random = TRUE, iter = 100),
        compute_stability_index(expr_gdsc_p, expr_ccle_p, "gdsc", "ccle", "ciri", random = TRUE, iter = 100)
    )

    print("Starting isoform exprssion stability indices")
    transcript_stability <- rbind(
        compute_stability_index(expr_gcsi_i, expr_ccle_i, "gcsi", "ccle", "circ", random = TRUE, iter = 100),
        compute_stability_index(expr_gcsi_i, expr_gdsc_i, "gcsi", "gdsc", "circ", random = TRUE, iter = 100),
        compute_stability_index(expr_gdsc_i, expr_ccle_i, "gdsc", "ccle", "circ", random = TRUE, iter = 100)
    )

    print("Starting gene expression random  stability indices")
    gene_stability_random <- rbind(
        compute_stability_index(expr_gcsi_p, expr_ccle_p, "gcsi", "ccle", "cfnd", random = TRUE, iter = 100),
        compute_stability_index(expr_gcsi_p, expr_gdsc_p, "gcsi", "gdsc", "cfnd", random = TRUE, iter = 100),
        compute_stability_index(expr_gdsc_p, expr_ccle_p, "gdsc", "ccle", "cfnd", random = TRUE, iter = 100)
    )

    print("Starting isoform expression random stability indices")
    transcript_stability_random <- rbind(
        compute_stability_index(expr_gcsi_i, expr_ccle_i, "gcsi", "ccle", "fcrc", random = TRUE, iter = 100),
        compute_stability_index(expr_gcsi_i, expr_gdsc_i, "gcsi", "gdsc", "fcrc", random = TRUE, iter = 100),
        compute_stability_index(expr_gdsc_i, expr_ccle_i, "gdsc", "ccle", "fcrc", random = TRUE, iter = 100)
    )

    save(gene_stability, transcript_stability, 
     file = "../results/data/temp/gene_isoform_stability.RData")

    save(gene_stability_random, isof_stability_random, 
     file = "../results/data/temp/gene_isoform_stability_random.RData")

} else if (analysis == "lungs") {
    load("../results/data/temp/gene_isoform_stability.RData")
    load("../results/data/temp/gene_isoform_stability_random.RData")
}

colnames(transcript_stability) <- c("gcsi_ccle_spearman", "gcsi_gdsc_spearman", "gdsc_ccle_spearman")  
colnames(gene_stability) <- c("gcsi_ccle_spearman", "gcsi_gdsc_spearman", "gdsc_ccle_spearman")  


############################################################
# Compute median expression for isoform stability
############################################################

# format matrices
gcsi_matrix <- as.matrix(expr_gcsi_i)                                 
ccle_matrix <- as.matrix(expr_ccle_i)  
gdsc_matrix <- as.matrix(expr_gdsc_i)

# compute median
gcsi_median <- as.numeric(rowMedians(gcsi_matrix))                                 
ccle_median <- as.numeric(rowMedians(ccle_matrix))                                  
gdsc_median <- as.numeric(rowMedians(gdsc_matrix))  

# compile isoforms into one data frame                                 
transcript_stability$transcript_id <- transcripts
rownames(transcript_stability) <- transcript_stability$transcript_id
transcript_stability$gcsi_median <- gcsi_median
transcript_stability$ccle_median <- ccle_median
transcript_stability$gdsc_median <- gdsc_median

#remove transcripts that have median expression of 0 across all datasets
transcript_stability <- transcript_stability[-which(transcript_stability$gcsi_median == 0 & transcript_stability$ccle_median == 0 & transcript_stability$gdsc_median == 0),]  
   

############################################################
# Compute median expression for gene exp stability
############################################################

# format matrices
gcsi_matrix <- as.matrix(expr_gcsi_p)                                 
ccle_matrix <- as.matrix(expr_ccle_p)  
gdsc_matrix <- as.matrix(expr_gdsc_p)

# compute median      
gcsi_median <- as.numeric(rowMedians(gcsi_matrix))                                 
ccle_median <- as.numeric(rowMedians(ccle_matrix))                                  
gdsc_median <- as.numeric(rowMedians(gdsc_matrix)) 
        
# compile gene expression into one data frame                                 
gene_stability$gene_id <- genes
rownames(gene_stability) <- gene_stability$gene_id
gene_stability$gcsi_median <- gcsi_median
gene_stability$ccle_median <- ccle_median
gene_stability$gdsc_median <- gdsc_median

############################################################
# Compute average spearman correlations (Table 2)
############################################################

colMeans(gene_stability[,1:3])
colMeans(transcript_stability[,1:3])
colMeans(ciri_stability[,1:3])
colMeans(circ_stability[,1:3])
colMeans(cfnd_stability[,1:3])
colMeans(fcrc_stability[,1:3])

colMeans(gene_stability_random)
colMeans(isof_stability_random)
colMeans(ciri_stability_random)
colMeans(circ_stability_random)
colMeans(cfnd_stability_random)
colMeans(fcrc_stability_random)

############################################################
# Compute t-test (non-random > random) (Table 2)
############################################################

ttest_si <- function(nonrandom, random) {
    nonrandom <- nonrandom[,1:3]
    nonrandom <- c(nonrandom[,1], nonrandom[,2], nonrandom[,3])
    random <- c(random[,1], random[,2], random[,3])

    t.test(nonrandom, random, alternative = "greater")
}

ttest_si(gene_stability, gene_stability_random)
ttest_si(transcript_stability, isof_stability_random)
ttest_si(ciri_stability, ciri_stability_random)
ttest_si(circ_stability, circ_stability_random)
ttest_si(cfnd_stability, cfnd_stability_random)
ttest_si(fcrc_stability, fcrc_stability_random)

############################################################
# Format stability index matrices for plotting
############################################################

# function to format stability dataframes
format_df <- function(df, label, random = "NonRandom") {
    if (ncol(df) == 4) {df <- df[,c(1:3)]}
    colnames(df) <- c("gCSI/CCLE", "gCSI/GDSC2", "GDSC2/CCLE")
    toPlot <- melt(df)
    colnames(toPlot) <- c("PSet", "Stability")
    toPlot$label <- label
    toPlot$random <- random
    return(toPlot)
}

ciri_stability <- format_df(ciri_stability, "CIRI2", "NonRandom")
circ_stability <- format_df(circ_stability, "CIRCexplorer2", "NonRandom")
cfnd_stability <- format_df(cfnd_stability, "circRNA_finder", "NonRandom")
fcrc_stability <- format_df(fcrc_stability, "find_circ", "NonRandom")

ciri_stability_random <- format_df(ciri_stability_random, "CIRI2", "Random")
circ_stability_random <- format_df(circ_stability_random, "CIRCexplorer2", "Random")
cfnd_stability_random <- format_df(cfnd_stability_random, "circRNA_finder", "Random")
fcrc_stability_random <- format_df(fcrc_stability_random, "find_circ", "Random")


transcript_stability <- format_df(transcript_stability[,c("gcsi_ccle_spearman", "gcsi_gdsc_spearman", "gdsc_ccle_spearman")], "Isoforms")
gene_stability <- format_df(gene_stability[,c("gcsi_ccle_spearman", "gcsi_gdsc_spearman", "gdsc_ccle_spearman")], "Gene Expression")
gene_stability[is.na(gene_stability)] <- 0

gene_stability_random <- format_df(gene_stability_random, "Gene Expression", "Random")
transcript_stability_random <- format_df(isof_stability_random, "Isoforms", "Random")


############################################################
# Plot stability index distribution
############################################################

# merge results for plotting
toPlot <- rbind(
    gene_stability,
    transcript_stability,
    ciri_stability,
    circ_stability,
    cfnd_stability,
    fcrc_stability,
    gene_stability_random,
    transcript_stability_random,
    ciri_stability_random,
    circ_stability_random,
    cfnd_stability_random,
    fcrc_stability_random
)
toPlot <- reshape2::melt(toPlot)
toPlot$label <- factor(toPlot$label, levels = names(dataset_pal))

# plot distribution of stability indices
filename <- paste0("../results/figures/figure2/stability_", analysis, ".png")
png(filename, width=250, height=200, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = label, y = value, fill = label)) +
    geom_boxplot(data = toPlot, aes(alpha = random)) +
    scale_fill_manual(values = dataset_pal) +
    scale_alpha_manual(values = c(1, 0.2)) +
    facet_grid(factor(PSet)~.) +
    theme_classic() + labs(x = "", fill = "", y = "Stability Index", alpha = "Randomization") +
    theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
        legend.key.size = unit(0.7, 'cm')) +
    geom_hline(yintercept = 0, linetype = "dotted")
dev.off()

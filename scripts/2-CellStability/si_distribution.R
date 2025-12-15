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

valid <- c("circ", "circ_fi", "GE")
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
    circ_fi = "../data/processed_cellline/common_samples/",
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

if (analysis == "circ_fi") {
    # helper function to keep only transcripts in 2/3 of the psets
    get_common <- function(gcsi, ccle, gdsc) {
        all_transcripts <- c(colnames(gcsi), colnames(ccle), colnames(gdsc))
        freq <- table(all_transcripts)
        to_keep <- names(freq[freq >= 2])
        return(to_keep)
    }
    # get common circRNA transcripts found in all psets for each pipeline
    ciri_common <- get_common(ciri_gcsi, ciri_ccle, ciri_gdsc)
    circ_common <- get_common(circ_gcsi, circ_ccle, circ_gdsc)
    cfnd_common <- get_common(cfnd_gcsi, cfnd_ccle, cfnd_gdsc)
    fcrc_common <- get_common(fcrc_gcsi, fcrc_ccle, fcrc_gdsc)
} else {
    # keep transcripts found in all psets for each pipeline
    ciri_common <- intersect(intersect(colnames(ciri_gcsi), colnames(ciri_ccle)), colnames(ciri_gdsc))
    circ_common <- intersect(intersect(colnames(circ_gcsi), colnames(circ_ccle)), colnames(circ_gdsc))
    cfnd_common <- intersect(intersect(colnames(cfnd_gcsi), colnames(cfnd_ccle)), colnames(cfnd_gdsc))
    fcrc_common <- intersect(intersect(colnames(fcrc_gcsi), colnames(fcrc_ccle)), colnames(fcrc_gdsc))
}

############################################################
# Keep only common transcripts across all PSets
############################################################

#' Helper function to order circRNA dataframes and keep common transccripts
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

# save subsetted dataframes for get_features.R (feature importance)
if (analysis == "circ_fi") {
save(ciri_gcsi, ciri_ccle, ciri_gdsc, circ_gcsi, circ_ccle, circ_gdsc,
     cfnd_gcsi, cfnd_ccle, cfnd_gdsc, fcrc_gcsi, fcrc_ccle, fcrc_gdsc,
     file = "../results/data/temp/circ_stability_subsetdf2.RData")
}

############################################################
# Compute pairwise spearman correlations
############################################################

print("Starting CIRI2 stability indices")
ciri_stability <- rbind(
    compute_stability_index(ciri_gcsi, ciri_ccle, "gCSI", "CCLE", "CIRI2"),
    compute_stability_index(ciri_gcsi, ciri_gdsc, "gCSI", "GDSC2", "CIRI2"),
    compute_stability_index(ciri_gdsc, ciri_ccle, "GDSC2", "CCLE", "CIRI2")
)

print("Starting CIRCexplorer2 stability indices")
circ_stability <- rbind(
    compute_stability_index(circ_gcsi, circ_ccle, "gCSI", "CCLE", "CIRCexplorer2"),
    compute_stability_index(circ_gcsi, circ_gdsc, "gCSI", "GDSC2", "CIRCexplorer2"),
    compute_stability_index(circ_gdsc, circ_ccle, "GDSC2", "CCLE", "CIRCexplorer2")
)

print("Starting circRNA_finder stability indices")
cfnd_stability <- rbind(
    compute_stability_index(cfnd_gcsi, cfnd_ccle, "gCSI", "CCLE", "circRNA_finder"),
    compute_stability_index(cfnd_gcsi, cfnd_gdsc, "gCSI", "GDSC2", "circRNA_finder"),
    compute_stability_index(cfnd_gdsc, cfnd_ccle, "GDSC2", "CCLE", "circRNA_finder")
)

print("Starting find_circ stability indices")
fcrc_stability <- rbind(
    compute_stability_index(fcrc_gcsi, fcrc_ccle, "gCSI", "CCLE", "find_circ"),
    compute_stability_index(fcrc_gcsi, fcrc_gdsc, "gCSI", "GDSC2", "find_circ"),
    compute_stability_index(fcrc_gdsc, fcrc_ccle, "GDSC2", "CCLE", "find_circ")
)

save(ciri_stability, circ_stability, cfnd_stability, fcrc_stability,
     file = paste0("../results/data/stability/", analysis, "_stability.RData"))

############################################################
# Compute pairwise spearman correlations with randomization
############################################################

print("Starting CIRI2 random stability indices")
ciri_stability_random <- rbind(
    compute_stability_index(ciri_gcsi, ciri_ccle, "gCSI", "CCLE", "CIRI2", random = TRUE, iter = 100),
    compute_stability_index(ciri_gcsi, ciri_gdsc, "gCSI", "GDSC2", "CIRI2", random = TRUE, iter = 100),
    compute_stability_index(ciri_gdsc, ciri_ccle, "GDSC2", "CCLE", "CIRI2", random = TRUE, iter = 100)
)

print("Starting CIRCexplorer2 random stability indices")
circ_stability_random <- rbind(
    compute_stability_index(circ_gcsi, circ_ccle, "gCSI", "CCLE", "CIRCexplorer2", random = TRUE, iter = 100),
    compute_stability_index(circ_gcsi, circ_gdsc, "gCSI", "GDSC2", "CIRCexplorer2", random = TRUE, iter = 100),
    compute_stability_index(circ_gdsc, circ_ccle, "GDSC2", "CCLE", "CIRCexplorer2", random = TRUE, iter = 100)
)

print("Starting circRNA_finder random stability indices")
cfnd_stability_random <- rbind(
    compute_stability_index(cfnd_gcsi, cfnd_ccle, "gCSI", "CCLE", "circRNA_finder", random = TRUE, iter = 100),
    compute_stability_index(cfnd_gcsi, cfnd_gdsc, "gCSI", "GDSC2", "circRNA_finder", random = TRUE, iter = 100),
    compute_stability_index(cfnd_gdsc, cfnd_ccle, "GDSC2", "CCLE", "circRNA_finder", random = TRUE, iter = 100)
)

print("Starting find_circ random stability indices")
fcrc_stability_random <- rbind(
    compute_stability_index(fcrc_gcsi, fcrc_ccle, "gCSI", "CCLE", "find_circ", random = TRUE, iter = 100),
    compute_stability_index(fcrc_gcsi, fcrc_gdsc, "gCSI", "GDSC2", "find_circ", random = TRUE, iter = 100),
    compute_stability_index(fcrc_gdsc, fcrc_ccle, "GDSC2", "CCLE", "find_circ", random = TRUE, iter = 100)
)

save(ciri_stability_random, circ_stability_random, cfnd_stability_random, fcrc_stability_random, 
     file = paste0("../results/data/stability/", analysis, "_stability_random.RData"))

############################################################
# SI for isoform and gene expression data
############################################################

# only do if analysis == "circ" (do once)

if (analysis == "circ") {

    # load in data
    load("../results/data/isoform_expression.RData")
    load("../results/data/gene_expression.RData")

    print("Starting gene expression stability indices")
    gene_stability <- rbind(
        compute_stability_index(expr_gcsi_p, expr_ccle_p, "gCSI", "CCLE", "Gene_Expression"),
        compute_stability_index(expr_gcsi_p, expr_gdsc_p, "gCSI", "GDSC2", "Gene_Expression"),
        compute_stability_index(expr_gdsc_p, expr_ccle_p, "GDSC2", "CCLE", "Gene_Expression")
    )

    print("Starting isoform exprssion stability indices")
    transcript_stability <- rbind(
        compute_stability_index(expr_gcsi_i, expr_ccle_i, "gCSI", "CCLE", "Isoform_Expression"),
        compute_stability_index(expr_gcsi_i, expr_gdsc_i, "gCSI", "GDSC2", "Isoform_Expression"),
        compute_stability_index(expr_gdsc_i, expr_ccle_i, "GDSC2", "CCLE", "Isoform_Expression")
    )

    print("Starting gene expression random  stability indices")
    gene_stability_random <- rbind(
        compute_stability_index(expr_gcsi_p, expr_ccle_p, "gCSI", "CCLE", "Gene_Expression", random = TRUE, iter = 100),
        compute_stability_index(expr_gcsi_p, expr_gdsc_p, "gCSI", "GDSC2", "Gene_Expression", random = TRUE, iter = 100),
        compute_stability_index(expr_gdsc_p, expr_ccle_p, "GDSC2", "CCLE", "Gene_Expression", random = TRUE, iter = 100)
    )

    print("Starting isoform expression random stability indices")
    transcript_stability_random <- rbind(
        compute_stability_index(expr_gcsi_i, expr_ccle_i, "gCSI", "CCLE", "Isoform_Expression", random = TRUE, iter = 100),
        compute_stability_index(expr_gcsi_i, expr_gdsc_i, "gCSI", "GDSC2", "Isoform_Expression", random = TRUE, iter = 100),
        compute_stability_index(expr_gdsc_i, expr_ccle_i, "GDSC2", "CCLE", "Isoform_Expression", random = TRUE, iter = 100)
    )

    save(gene_stability, transcript_stability,
     file = "../results/data/stability/gene_isoform_stability.RData")

    save(gene_stability_random, transcript_stability_random,
     file = "../results/data/stability/gene_isoform_stability_random.RData")

} else if (analysis == "lungs") {
    load("../results/data/stability/gene_isoform_stability.RData")
    load("../results/data/stability/gene_isoform_stability_random.RData")
}

############################################################
# Compute average spearman correlations (Table 2)
############################################################

# helper function to compute mean of given dataset pair
get_mean <- function(stability_df) {
    pairs <- unique(stability_df$label)
    for (pair in pairs) {
        m <- mean(stability_df$stability[stability_df$label == pair])
        print(paste("Mean of", pair, ":", m))
    }
}

############################################################
# Compute t-test (Table 2)
############################################################

#' Helper function to compute t-test
#' 
#' @param greater_df data.frame. The expected larger mean (e.g. non-random, gene expression)
#' @param less_df data.frame. The expected smaller mean (e.g. random, circRNA)
#' 
ttest_si <- function(greater_df, less_df) {
    t.test(greater_df$stability, less_df$stability, alternative = "greater")
}

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
toPlot$dataset <- factor(toPlot$dataset, levels = names(dataset_pal))

# plot distribution of stability indices
filename <- paste0("../results/figures/figure2/stability_", analysis, ".png")
png(filename, width=250, height=200, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = dataset, y = stability, fill = dataset, alpha = random)) +
    scale_fill_manual(values = dataset_pal) +
    scale_alpha_manual(values = c(1, 0.2)) +
    facet_grid(factor(label)~.) +
    theme_classic() + labs(x = "", fill = "", y = "Stability Index", alpha = "Randomization") +
    theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
        legend.key.size = unit(0.7, 'cm')) +
    geom_hline(yintercept = 0, linetype = "dotted")
dev.off()

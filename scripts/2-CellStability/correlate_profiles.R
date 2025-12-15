# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(plyr)
    library(dplyr)
    library(reshape2)
    library(tidyverse)
    library(umap)
    library(ggplot2)
    library(ggpubr)
})

options(stringsAsFactors = FALSE)
set.seed(101)

source("utils/palettes.R")
source("utils/correlate_profiles.R")

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
# Filter matrices by transcripts
############################################################

# get common transcripts
transcripts <- data.frame(c(colnames(ciri_gcsi), colnames(ciri_ccle), colnames(ciri_gdsc),
                            colnames(circ_gcsi), colnames(circ_ccle), colnames(circ_gdsc),
                            colnames(cfnd_gcsi), colnames(cfnd_ccle), colnames(cfnd_gdsc),
                            colnames(fcrc_gcsi), colnames(fcrc_ccle), colnames(fcrc_gdsc)))
colnames(transcripts) <- "transcriptID"
transcript_counts <- transcripts %>% count(transcriptID)

# remove transcripts only in one method_dataset object (keep any in 2 or more)
transcript_counts <- transcript_counts[-which(transcript_counts$n == 1),]
common_transcripts <- c(transcript_counts$transcriptID) #n: 29764

############################################################
# Spearman corr of transcript exp across biological reps
############################################################

# compute circRNA profile correlations
ciri_corr <- correlate_profiles(ciri_gcsi, ciri_ccle, ciri_gdsc)
circ_corr <- correlate_profiles(circ_gcsi, circ_ccle, circ_gdsc)
cfnd_corr <- correlate_profiles(cfnd_gcsi, cfnd_ccle, cfnd_gdsc)
fcrc_corr <- correlate_profiles(fcrc_gcsi, fcrc_ccle, fcrc_gdsc)

# only do if analysis == "circ" (do once)
if (analysis == "circ") {

    # load in data
    load("../results/data/isoform_expression.RData")
    load("../results/data/gene_expression.RData")

    gexpr_corr <- correlate_profiles(expr_gcsi_p, expr_ccle_p, expr_gdsc_p, circ = FALSE)
    isoform_corr <- correlate_profiles(expr_gcsi_i, expr_ccle_i, expr_gdsc_i, circ = FALSE)

    save(gexpr_corr, isoform_corr, file = "../results/data/gene_isof_corr_expr.RData")

} else {
    load("../results/data/gene_isof_corr_expr.RData")
}

filename <- paste0("../results/data/", analysis, "_corr_expr.RData")
save(gexpr_corr, isoform_corr, ciri_corr, circ_corr, cfnd_corr, fcrc_corr, file = filename)


############################################################
# Compute average spearman correlations (Table 3)
############################################################

# helper function to get average spearman
get_avg <- function(corr_df) {
    print(colMeans(corr_df[,2:4], na.rm = TRUE))
}

############################################################
# Format dataframe for plotting
############################################################

# function to format correlation dataframes
format_df <- function(df, label) {
    colnames(df) <- c("Cells", "gCSI/CCLE", "gCSI/GDSC2", "GDSC2/CCLE")
    toPlot <- reshape2::melt(df)
    colnames(toPlot) <- c("Cells", "PSet", "Correlation")
    toPlot$label <- label
    return(toPlot)
}

gexpr_corr <- format_df(gexpr_corr, "Gene Expression")
isoform_corr <- format_df(isoform_corr, "Isoforms")
ciri_corr <- format_df(ciri_corr, "CIRI2")
circ_corr <- format_df(circ_corr, "CIRCexplorer2")
cfnd_corr <- format_df(cfnd_corr, "circRNA_finder")
fcrc_corr <- format_df(fcrc_corr, "find_circ")

# merge results for plotting
toPlot <- rbind(
    gexpr_corr, isoform_corr,
    ciri_corr, circ_corr,
    cfnd_corr, fcrc_corr
)
toPlot[is.na(toPlot)] <- 0
toPlot$label <- factor(toPlot$label, levels = names(dataset_pal))

############################################################
# Create plots
############################################################

filename <- paste0("../results/figures/figure2/correlate_reps", analysis, ".png")
png(filename, width=200, height=75, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = Correlation, y = label)) + 
    geom_violin(aes(fill = label, alpha = PSet), 
                    scale = "width", width = 1.1, 
                    position = position_dodge(width = 0.8)) + 
    theme_classic() + 
    labs(x = "Spearman Correlation", fill = "", alpha = "", y = "") +
    scale_fill_manual(values = dataset_pal) + 
    scale_alpha_manual(values = c(0.2, 0.5, 0.9)) +
    scale_y_discrete(limits = rev(names(dataset_pal))) + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3), legend.key.size = unit(0.5, 'cm')) 
dev.off()

filename <- paste0("../results/figures/figure5/corr_dist_density_", analysis, ".png")
png(filename, width=150, height=50, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = Correlation)) + 
    geom_density(aes(fill = label), alpha = 0.4, size = 0.5) +
    theme_classic() +  
    scale_fill_manual(values = dataset_pal) +
    labs(x = "Spearman Correlation", y = "Density", fill = "") +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3), legend.key.size = unit(0.4, 'cm'))
dev.off()


############################################################
# Create UMAP projections
############################################################

ciri_umap <- umap_cells(ciri_gcsi, ciri_ccle, ciri_gdsc)
circ_umap <- umap_cells(circ_gcsi, circ_ccle, circ_gdsc)
cfnd_umap <- umap_cells(cfnd_gcsi, cfnd_ccle, cfnd_gdsc)
fcrc_umap <- umap_cells(fcrc_gcsi, fcrc_ccle, fcrc_gdsc)

gexpr_umap <- umap_cells(expr_gcsi_p, expr_ccle_p, expr_gdsc_p, circ = FALSE)
isoform_umap <- umap_cells(expr_gcsi_i, expr_ccle_i, expr_gdsc_i, circ = FALSE)

############################################################
# Plot UMAP projections
############################################################

# function to plot umap
plot_umap <- function(umap_df, title) {
  p <- ggplot(data = umap_df, aes(x = UMAP1, y = UMAP2, group = cell_line)) +
    geom_line(show.legend = FALSE) +
    geom_point(aes(color = dataset, shape = dataset), size = 3) +
    scale_color_manual(
        guide = guide_legend(reverse = FALSE, title = "Dataset"),
        values = pset_pal) + 
    guides(shape = guide_legend(ncol = 1)) + #make legend one column instead of 3
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
    labs(title = title) +
    scale_shape_discrete(name = "Dataset")
  return(p)
}

p1 <- plot_umap(gexpr_umap, "Gene Expression")
p2 <- plot_umap(isoform_umap, "Isoform Expression")
p3 <- plot_umap(ciri_umap, "CIRI2 circRNA Expression")
p4 <- plot_umap(circ_umap, "CIRCexplorer2 circRNA Expression")
p5 <- plot_umap(cfnd_umap, "circRNA_finder circRNA Expression")
p6 <- plot_umap(fcrc_umap, "find_circ circRNA Expression")

filename <- paste0("../results/figures/suppfig2/umaps_", analysis, ".png")
png(filename, width=400, height=225, units='mm', res = 600, pointsize=80)
ggarrange(
    p1, p3, p5, p2, p4, p6,
    ncol = 3, nrow = 2,
    common.legend = TRUE,
    legend = "right")
dev.off()
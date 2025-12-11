# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(stats)
    library(ggplot2)
    library(ComplexHeatmap)
    library(ggh4x)
    library(stringr)
    library(RColorBrewer)
    library(ggrepel)
    library(ggsignif)
})

source("utils/biomarker_analysis.R")
source("utils/palettes.R")


############################################################
# Load in circRNA expression data
############################################################

# helper function to load in expression matrices
load_counts <- function(path, filename) {
    counts <- fread(paste0(path, filename), data.table = FALSE)
    counts[is.na(counts)] <- 0
    rownames(counts) <- counts$V1
    counts$V1 <- NULL
    # filter
    counts <- counts[,colnames(counts)[colSums(counts != 0) >= thres], drop = FALSE]
    # zero to NA for downstream analysis
    counts[counts == 0] <- NA
    return(counts)
}

# load circRNA expression data
path <- "../data/processed_cellline/merged_all_samples/"
thres <- 15

gcsi_df <- load_counts(path, "gcsi_counts.tsv")
ccle_df <- load_counts(path, "ccle_counts.tsv")
gdsc_df <- load_counts(path, "gdsc_counts.tsv")

save(
    gcsi_df, ccle_df, gdsc_df,
    file = "../results/data/biomarker_analysis_circ.RData"
)

############################################################
# Load in drug response data and subset
############################################################

# load in drug sensitivity from PSets
load("../data/temp/sensitivity_data.RData")

# drugs of interest
gcsi_ctrp <- intersect(rownames(gcsi_sen), rownames(ctrp_sen))
gcsi_gdsc <- intersect(rownames(gcsi_sen), rownames(gdsc_sen))
ctrp_gdsc <- intersect(rownames(ctrp_sen), rownames(gdsc_sen))

# keep drugs of interest
gcsi_sen <- gcsi_sen[rownames(gcsi_sen) %in% c(gcsi_ctrp, gcsi_gdsc), ]
ctrp_sen <- ctrp_sen[rownames(ctrp_sen) %in% c(gcsi_ctrp, ctrp_gdsc), ]
gdsc_sen <- gdsc_sen[rownames(gdsc_sen) %in% c(gcsi_gdsc, ctrp_gdsc), ]

############################################################
# Drug response associations
############################################################

gcsi_bin_dr <- binary_dr(gcsi_df, gcsi_sen)
ccle_bin_dr <- binary_dr(ccle_df, ctrp_sen)
gdsc_bin_dr <- binary_dr(gdsc_df, gdsc_sen)

# save drug response associations
dr_out <- "../results/data/bin_dr/circ_"

save(gcsi_bin_dr, file = paste0(dr_out, "gcsi_bin0.RData"))
save(ccle_bin_dr, file = paste0(dr_out, "ccle_bin0.RData"))
save(gdsc_bin_dr, file = paste0(dr_out, "gdsc_bin0.RData"))


############################################################
# Subset for significant associations
############################################################

# helper function to subset for significant associations
subset_dr <- function(dr, dataset) {
    dr <- dr[which(dr$FDR_drug < 0.1), ]
    dr$dataset <- dataset
    return(dr)
}

gcsi_fdr_b <- subset_dr(gcsi_bin_dr, "gCSI")
gdsc_fdr_b <- subset_dr(gdsc_bin_dr, "CCLE")
ccle_fdr_b <- subset_dr(ccle_bin_dr, "GDSC2")

############################################################
# Save Table 4
############################################################

table4 <- rbind(gcsi_fdr_b, ccle_fdr_b, gdsc_fdr_b)
colnames(table4) <- c(
    "Drug", "circRNA", "W", "pvalue", "diff.in.AAC",
    "Exp", "No.Samples", "rm", "rm2", "pair", "FDR", "dataset"
)
table4 <- table4[,-which(colnames(table4) %in% c("rm", "rm2"))]
write.csv(table4, "../results/data/table4.csv", quote = FALSE, row.names = FALSE)

############################################################
# Upset plot of overlapping biomarkers
############################################################

# create list object of transcripts for upset plot
comb_mat <- make_comb_mat(list(
    gCSI = gcsi_fdr_b$pair,
    CCLE = ccle_fdr_b$pair,
    GDSC = gdsc_fdr_b$pair
))

# upset plot
filename <- "../results/figures/figure3/upset_dr_bin0.png"
png(filename, width = 8, height = 3, res = 600, units = "in")
UpSet(comb_mat, set_order = names(pset_pal),
    top_annotation = upset_top_annotation(comb_mat, add_numbers = TRUE),
    comb_order = order(-comb_size(comb_mat)),
    right_annotation = upset_right_annotation(comb_mat, add_numbers = TRUE))
dev.off()

############################################################
# Stats of top biomarkers across all datasets
############################################################

# save top 20 associations from each PSet
gcsi_top <- gcsi_fdr_b[order(abs(gcsi_fdr_b$diff), decreasing = T),][1:20,]
ccle_top <- ccle_fdr_b[order(abs(ccle_fdr_b$diff), decreasing = T),][1:20,]
gdsc_top <- gdsc_fdr_b[order(abs(gdsc_fdr_b$diff), decreasing = T),][1:20,]

top_biomarkers <- unique(c(gcsi_top$pair, ccle_top$pair, gdsc_top$pair))

# annotate transcript and drugs
toPlot <- rbind(
    format_bin(gcsi_bin_dr, gcsi_df, gcsi_sen),
    format_bin(ccle_bin_dr, ccle_df, ctrp_sen),
    format_bin(gdsc_bin_dr, gdsc_df, gdsc_sen)
)
toPlot$dataset <- factor(toPlot$dataset, levels = c(names(pset_pal)))

# order by drug per dataset
order_toplot <- c(
    top_biomarkers[1:20][order(top_biomarkers[1:20])],
    top_biomarkers[21:40][order(top_biomarkers[21:40])],
    top_biomarkers[41:60][order(top_biomarkers[41:60])]
)
toPlot$Pair <- factor(toPlot$Pair, levels = rev(order_toplot))

# create labels
toPlot$label <- paste(toPlot$Drug, sub("^([^\\.]+\\.)([0-9]{3}).*$", "\\1\\2", toPlot$Feature), sep = ":    ")

# set up values to colour 
toPlot$Diff[is.na(toPlot$Status)] <- NA
toPlot$Diff[which(toPlot$Status == "ND")] <- 0

# plot overlapping biomarkers
png("../results/figures/figure3/top_biomarkers.png", width = 6, height = 10, res = 600, units = "in")
ggplot(toPlot, aes(x = dataset, y = Pair, fill = Diff)) + 
    geom_tile() +
    geom_hline(yintercept = c(20.5, 40.5), linetype = "dashed") +
    geom_text(aes(label = Status), size = 3) +
    scale_y_discrete(labels = setNames(toPlot$label, toPlot$Pair)) +
    theme_classic() + 
    scale_fill_gradient2(gradient_pal) +
    theme(
        panel.border = element_rect(color = "black", fill = NA, size = 0.5), 
        legend.title = element_text(size = 9)) +
    labs(fill = "Δ Mean AAC \n(Expression vs\nNon-expression)", y = "Drug-CircRNA Pair", x = "Dataset")
dev.off()

############################################################
# Plot Vorinostat_chr8.61680967.61684188 association
############################################################

# annotate by circ expression (high vs low)
toPlot <- rbind(
    label_exp(gcsi_df, gcsi_sen, "gCSI"),
    label_exp(ccle_df, ctrp_sen, "CCLE"),
    label_exp(gdsc_df, gdsc_sen, "GDSC")
)
toPlot$Dataset <- factor(toPlot$Dataset, levels = names(pset_pal))

# plot AAC by circ expression
png("../results/figures/figure3/vorinostat.png", width = 7, height = 5.5, res = 600, units = "in")
ggplot(toPlot, aes(x = Dataset, y = AAC, fill = Label)) +
    geom_boxplot() +
    theme_classic() +
    scale_fill_manual(values = vorinostat_pal) +
    ylim(c(-0.01, 0.6)) +
    geom_signif(
        y_position = c(0.47, 0.49, 0.54), 
        xmin = c(0.8, 1.8, 2.8), xmax = c(1.2, 2.2, 3.2),
        annotation = c("*", "**", "**"),
        tip_length = 0.01,
        textsize = 4) +
    theme(
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        legend.key.size = unit(0.5, 'cm'),
        legend.title = element_text(size = 9)) +
    labs(y = "Vorinostat Response (AAC)", fill = "Expression of\nchr8.61680967.61684188")
dev.off()

############################################################
# Pair-wise correlation plots across PSets
############################################################

# only doing CCLE and GDSC since there is only 20 fdr sig in gCSI
plot_corr_pset(ccle_bin_dr, gdsc_bin_dr, "../results/figures/figure3/CCLE_GDSC_bin0")

############################################################
# Plot distribution of 0 expression
############################################################

# function to plot the distribution of zero expression cross cell lines
plot_zeros <- function(bin, filename) {

    # format dataframe for plotting
    bin <- bin[,colnames(bin) %in% c("Feature", "num_high", "num_samples")] |> unique()
    bin <- bin[order(bin$num_high, decreasing = TRUE),]
    toPlot <- melt(bin)

    # get circRNAs order by proportion of >0 expression
    circs <- bin$Feature
    toPlot$Feature <- factor(toPlot$Feature, levels = unique(circs))

    png(paste0("../results/figures/suppfig3/dist_zeros_", filename, ".png"), width = 4, height = 5, res = 600, units = "in")
    print({ggplot(toPlot, aes(y = Feature, x = value, fill = variable)) + 
        geom_bar(position="dodge", stat="identity") +
        scale_fill_manual("Expression", values = zero_pal, label = c(">0", "0")) +
        theme(axis.text.y=element_blank(),
              axis.ticks.y=element_blank()) +
        labs(x = "Number of Cell Lines", y = "CircRNA Transcript")
    })
    dev.off()
}

plot_zeros(gcsi_bin_dr, "gcsi_bin0")
plot_zeros(ccle_bin_dr, "ccle_bin0")
plot_zeros(gdsc_bin_dr, "gdsc_bin0")
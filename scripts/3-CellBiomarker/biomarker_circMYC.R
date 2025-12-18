# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggrepel)
})

set.seed(101)

source("utils/biomarker_analysis.R")
source("utils/palettes.R")

############################################################
# Load in circRNA expression data
############################################################

# from circ_to_myc.R
load("../results/data/circ_myc.RData")

############################################################
# Robust normalization of circRNA expression
############################################################

# function to count and remove low expressed transcripts
# function from combine_pipeline_counts.R
robust_norm <- function(df) {

    df <- na.omit(df)
    for (circ in colnames(df)) {

        exp <- df[[circ]]
        non_zero <- exp != 0
        
        # median and IQR on non-zero values only
        med <- median(exp[non_zero])
        iqr <- IQR(exp[non_zero])
        
        # avoid division by zero in case IQR is 0
        if (iqr == 0 || is.na(iqr)) {
            df[non_zero, circ] <- 0  # or NA if preferred
            if (is.na(iqr)) {
                print(paste("No expression for", circ))
            } else {
                print(paste("Issue with IQR for", circ))
            }
        } else {
            df[non_zero, circ] <- (exp[non_zero] - med) / iqr
        }
    }
    return(df)
}

ciri_gcsi <- robust_norm(ciri_gcsi)
ciri_gdsc <- robust_norm(ciri_gdsc)
ciri_ccle <- robust_norm(ciri_ccle)

circ_gcsi <- robust_norm(circ_gcsi)
circ_gdsc <- robust_norm(circ_gdsc)
circ_ccle <- robust_norm(circ_ccle)

cfnd_gcsi <- robust_norm(cfnd_gcsi)
cfnd_gdsc <- robust_norm(cfnd_gdsc)
cfnd_ccle <- robust_norm(cfnd_ccle)

fcrc_gcsi <- robust_norm(fcrc_gcsi)
fcrc_gdsc <- robust_norm(fcrc_gdsc)
fcrc_ccle <- robust_norm(fcrc_ccle)

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
# Binarize transcript expression by zero expression
############################################################

ciri_gcsi_bin <- binary_dr(ciri_gcsi, gcsi_sen, "CIRI_gCSI") #14, 0 FDR<0.05
ciri_gdsc_bin <- binary_dr(ciri_gdsc, gdsc_sen, "CIRI_GDSC") #73, 12 FDR<0.05
ciri_ccle_bin <- binary_dr(ciri_ccle, ctrp_sen, "CIRI_CCLE") #72, 3 FDR<0.05

circ_gcsi_bin <- binary_dr(circ_gcsi, gcsi_sen, "CIRC_gCSI") #14, 1 FDR<0.05
circ_gdsc_bin <- binary_dr(circ_gdsc, gdsc_sen, "CIRC_GDSC") #73, 13 FDR<0.05
circ_ccle_bin <- binary_dr(circ_ccle, ctrp_sen, "CIRC_CCLE") #72, 1 FDR<0.05

cfnd_gcsi_bin <- binary_dr(cfnd_gcsi, gcsi_sen, "CFND_gCSI") #14, 0 FDR<0.05
cfnd_gdsc_bin <- binary_dr(cfnd_gdsc, gdsc_sen, "CFND_GDSC") #73, 11 FDR<0.05
cfnd_ccle_bin <- binary_dr(cfnd_ccle, ctrp_sen, "CFND_CCLE") #72, 33 FDR<0.05

fcrc_gcsi_bin <- binary_dr(fcrc_gcsi, gcsi_sen, "FCRC_gCSI") #14, 4 FDR<0.05
fcrc_gdsc_bin <- binary_dr(fcrc_gdsc, gdsc_sen, "FCRC_GDSC") #73, 9 FDR<0.05
fcrc_ccle_bin <- binary_dr(fcrc_ccle, ctrp_sen, "FCRC_CCLE") #72, 50 FDR<0.05

# save dataframes
save(ciri_gcsi_bin, ciri_gdsc_bin, ciri_ccle_bin,
     circ_gcsi_bin, circ_gdsc_bin, circ_ccle_bin,
     cfnd_gcsi_bin, cfnd_gdsc_bin, cfnd_ccle_bin,
     fcrc_gcsi_bin, fcrc_gdsc_bin, fcrc_ccle_bin,
     file = "../results/data/circ_myc_norm_dr.RData")

############################################################
# Volcano Plots
############################################################

plot_volcano <- function(bin_dr) {

    # get x limits
    x <- max(c(abs(min(bin_dr$diff)), max(bin_dr$diff)))

    # make volcano plot
    p <- ggplot() +
        geom_point(data = bin_dr, aes(x = diff, y = -log(FDR)), color = "gray") +
        geom_point(data = bin_dr[!bin_dr$to_label == "", ], aes(x = diff, y = -log(FDR), color = label)) +
        geom_text_repel(
            data = bin_dr[!bin_dr$to_label == "", ],
            aes(x = diff, y = -log(FDR), label = to_label),
            force_pull = 0.5, max.overlaps = Inf, force = 5,
            box.padding = 0.4) +
        geom_hline(yintercept = -log(0.1), linetype = "dotted") +
        scale_color_manual(
            values = c(myc_pal), 
            labels = c(
                "CFND_CCLE" = "CCLE (circRNA_finder)", 
                "FCRC_CCLE" = "CCLE (find_circ)", "FCRC_GDSC" = 
                "GDSC2 (find_circ)"),
            name = "Dataset (Pipeline)") +
        theme_classic() + 
        xlim(-x, x) +
        theme(
            panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            legend.key.size = unit(0.5, 'cm')) +
        labs(x = "Δ Mean AAC (circMYC Expression vs Non-Expression)", y = "-log(FDR)")
    return(p)

}

toPlot <- rbind(
    ciri_gcsi_bin, ciri_gdsc_bin, ciri_ccle_bin,
    circ_gcsi_bin, circ_gdsc_bin, circ_ccle_bin,
    cfnd_gcsi_bin, cfnd_gdsc_bin, cfnd_ccle_bin,
    fcrc_gcsi_bin, fcrc_gdsc_bin, fcrc_ccle_bin
)
toPlot <- toPlot[!is.na(toPlot$diff),]

# label top 20 associations
toPlot$temp <- paste(toPlot$pair, toPlot$label, sep = "_")
toPlot <- toPlot[order(abs(toPlot$diff), decreasing = TRUE),]
top20 <- toPlot$temp[toPlot$FDR < 0.1][1:20]
toPlot$to_label <- ifelse(toPlot$temp %in% top20, as.character(toPlot$Drug), "")

# label by significance
toPlot$Label <- ifelse(toPlot$FDR < 0.05, "FDR\nSignificant", "Not FDR\nSignificant")
toPlot$Label <- factor(toPlot$Label, levels = c("FDR\nSignificant", "Not FDR\nSignificant"))
toPlot$label <- factor(
    toPlot$label, 
    levels = c(
        "CIRI_gCSI", "CIRI_CCLE", "CIRI_GDSC",
        "CIRC_gCSI", "CIRC_CCLE", "CIRC_GDSC",
        "CFND_gCSI", "CFND_CCLE", "CFND_GDSC",
        "FCRC_gCSI", "FCRC_CCLE", "FCRC_GDSC"
    )
)

png("../results/figures/figure3/myc_volcano.png", width = 6, height = 4, res = 600, units = "in")
plot_volcano(toPlot)
dev.off()


############################################################
# Stats of MYC drug associations in other PSets/Pipelines
############################################################

# save significant MYC drug associations
drugs <- toPlot[!toPlot$to_label == "",]$Drug |> as.character() |> unique()

sig_drug <- rbind(
    formatMYC(gcsi_sen, "CIRI_gCSI"),
    formatMYC(ctrp_sen, "CIRI_CCLE"),
    formatMYC(gdsc_sen, "CIRI_GDSC"),
    formatMYC(gcsi_sen, "CIRC_gCSI"),
    formatMYC(ctrp_sen, "CIRC_CCLE"),
    formatMYC(gdsc_sen, "CIRC_GDSC"),
    formatMYC(gcsi_sen, "CFND_gCSI"),
    formatMYC(ctrp_sen, "CFND_CCLE"),
    formatMYC(gdsc_sen, "CFND_GDSC"),
    formatMYC(gcsi_sen, "FCRC_gCSI"),
    formatMYC(ctrp_sen, "FCRC_CCLE"),
    formatMYC(gdsc_sen, "FCRC_GDSC")
)
sig_drug$PSet <- gsub(".*_", "", sig_drug$Label)
sig_drug$Pipeline <- rep(names(pipeline_pal), each = nrow(sig_drug)/4)
sig_drug$Pipeline <- factor(sig_drug$Pipeline, levels = names(pipeline_pal))
sig_drug$PSet[sig_drug$PSet=="GDSC"] <- "GDSC2"
sig_drug$PSet <- factor(sig_drug$PSet, levels = names(pset_pal))
sig_drug$Drug <- factor(sig_drug$Drug, levels = unique(rev(sig_drug$Drug[order(sig_drug$Drug)])))

# plot overlapping biomarkers
png("../results/figures/figure3/myc_top_biomarkers.png", width = 8.5, height = 5, res = 600, units = "in")
ggplot(sig_drug, aes(x = PSet, y = Drug, fill = Diff)) +
  geom_tile() +
  geom_text(aes(label = Status), size = 3) +
  facet_grid(. ~ Pipeline, scales = "free_x", space = "free_x") +
  labs(fill = "Δ Mean AAC\n(circMYC Expression\nvs Non-Expression)", y = "Drug", x = "Dataset") +
  theme_classic() +
  scale_fill_gradient2(low = "#9D3737", mid = "white", high = "#3670A0", midpoint = 0, na.value = "grey") +
  theme(
    strip.background = element_rect(fill = "#f0f0f0"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    legend.title = element_text(size = 9)
)
dev.off()
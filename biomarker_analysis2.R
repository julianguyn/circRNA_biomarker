# script to binarize circ expression by 0 expression

# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(stats)
    library(ggplot2)
    library(ComplexHeatmap)
    library(ggh4x)
    library(ggVennDiagram)
    library(stringr)
})


############################################################
# Specify analysis
############################################################

# analysis = c("circ", "GE")  
# circ for circRNA counts, GE for circ mapped to GE

if (analysis == "circ") {                                               
    path <- "../data/processed_cellline/merged_all_samples/"
    thres = 15
    dr_out <- "../results/data/bin_dr/circ_"
    lm_out <- "../results/data/lm_dr/circ_"
}
if (analysis == "GE") {                                                 
    path <- "../data/processed_cellline/mergedGE_common_samples/"
    thres = 10       
    dr_out <- "../results/data/bin_dr/GE_"    
    lm_out <- "../results/data/lm_dr/GE_"                           
}


############################################################
# Load in circRNA expression data
############################################################

# load circRNA expression data
# load from biomarker_analysis.R
load("../results/data/biomarker_analysis_circ.RData")     


############################################################
# Load in drug response data and subset
############################################################

# load in drug sensitivity from PSets
load("../data/temp/sensitivity_data.RData")

# drugs of interest
gcsi_ctrp <- intersect(rownames(gcsi_sen), rownames(ctrp_sen))
gcsi_gdsc <- intersect(rownames(gcsi_sen), rownames(gdsc_sen))
ctrp_gdsc <- intersect(rownames(ctrp_sen), rownames(gdsc_sen))

#drugs <- intersect(intersect(rownames(gcsi_sen), rownames(ctrp_sen)), rownames(gdsc_sen))
# "Bortezomib", "Crizotinib", "Docetaxel", "Erlotinib", "Pictilisib", "Gemcitabine", "Lapatinib", "Entinostat", "Paclitaxel", "Sirolimus", "Vorinostat" 

# keep drugs of interest
gcsi_sen <- gcsi_sen[rownames(gcsi_sen) %in% c(gcsi_ctrp, gcsi_gdsc), ]
ctrp_sen <- ctrp_sen[rownames(ctrp_sen) %in% c(gcsi_ctrp, ctrp_gdsc), ]
gdsc_sen <- gdsc_sen[rownames(gdsc_sen) %in% c(gcsi_gdsc, ctrp_gdsc), ]


############################################################
# Binarize transcript expression by zero expression
############################################################

# function to compute drug response associations
binary_dr <- function(counts_df, drug_df) {

    features <- colnames(counts_df)

    # create data frame to hold results
    combinations <- expand.grid(Drug = rownames(drug_df), Feature = features)
    combinations$num_samples <- combinations$num_high <- combinations$pval <-  combinations$W <- NA
    combinations$drug <- combinations$feature <- NA
    
    # initiate row count
    row = 1
    print(paste("Number of features:", length(features)))

    for (i in seq_along(features)) {
        
        print(i)
        feature <- features[i]

        # binarize transcript expression by median
        subset <- counts_df[,i, drop = FALSE]

        high_exp <- rownames(subset[complete.cases(subset), , drop = FALSE])
        low_exp <- rownames(subset)[-which(rownames(subset) %in% high_exp)]

        # loop through each drug

        for (j in seq_along(rownames(drug_df))) {
            
            high = drug_df[j,colnames(drug_df) %in% high_exp] |> as.numeric()
            low = drug_df[j,colnames(drug_df) %in% low_exp] |> as.numeric()

            if (length(high[!is.na(high)]) > 1 & length(low[!is.na(low)]) > 1) {

                # wilcoxon rank sum test
                res <- wilcox.test(high, low, alternative = "two.sided", exact = FALSE)

                # save results
                combinations$W[row] <- res$statistic
                combinations$pval[row] <- res$p.value
                combinations$num_samples[row] <- nrow(subset)
                combinations$num_high[row] <- length(high_exp)
                combinations$feature[row] <- feature
                combinations$drug[row] <- rownames(drug_df)[j]

            } else {
                # save results
                combinations$W[row] <- NA
                combinations$pval[row] <- NA
                combinations$num_samples[row] <- nrow(subset)
                combinations$num_high[row] <- NA
                combinations$feature[row] <- feature
                combinations$drug[row] <- rownames(drug_df)[j]
            }

            row <- row + 1
        }
    }
    combinations$FDR <- p.adjust(combinations$pval, method = "BH")

    # format dataframe for plotting 
    combinations <- combinations[order(combinations$W, decreasing = T),]
    combinations$rank <- seq_len(nrow(combinations))
    combinations$pair <- paste(combinations$Drug, combinations$Feature, sep = "_")
    return(combinations)
}

gcsi_bin_dr <- binary_dr(gcsi_df, gcsi_sen)
ccle_bin_dr <- binary_dr(ccle_df, ccle_sen)
gdsc_bin_dr <- binary_dr(gdsc_df, gdsc_sen)


############################################################
# Save drug response associations
############################################################

save(gcsi_bin_dr, file = paste0(dr_out, "gcsi_bin0.RData"))
save(ccle_bin_dr, file = paste0(dr_out, "ccle_bin0.RData"))
save(gdsc_bin_dr, file = paste0(dr_out, "gdsc_bin0.RData"))


############################################################
# Quantify (TO REMOVE)
############################################################

load(paste0(dr_out, "gcsi.RData"))
load(paste0(dr_out, "ccle.RData"))
load(paste0(dr_out, "gdsc.RData"))

# number of biomarker associations
# gCSI: 5488
# CCLE: 24600
# GDSC: 33507

# number of biomarker associations with pval < 0.05
# gCSI: 408
# CCLE: 1161
# GDSC: 4930

# number of biomarker associations with FDR < 0.05
# gCSI: 2
# CCLE: 0
# GDSC: 788


############################################################
# Subset for significant associations
############################################################

# function to subset for significant associations
subset_dr <- function(dr, type = "pval") {

    if (type == "pval") {
        dr <- dr[which(dr$pval < 0.05),] 
        return(dr)
    } else {
        dr <- dr[which(dr$FDR < 0.05),] 
        return(dr)
    }
}

# subset binarized drug response
gcsi_pval_b <- subset_dr(gcsi_bin_dr)
gdsc_pval_b <- subset_dr(gdsc_bin_dr)
ccle_pval_b <- subset_dr(ccle_bin_dr)

gcsi_fdr_b <- subset_dr(gcsi_bin_dr, type = "fdr")
gdsc_fdr_b <- subset_dr(gdsc_bin_dr, type = "fdr")
ccle_fdr_b <- subset_dr(ccle_bin_dr, type = "fdr")


############################################################
# Upset plot of overlapping biomarkers
############################################################

# function to create upset plot
plot_upset <- function(comb_mat, set_order, filename) {
    png(paste0("../results/figures/figure9/upset_dr_", filename, ".png"), width = 6, height = 4, res = 600, units = "in")
    print({UpSet(comb_mat, set_order = set_order,
        top_annotation = upset_top_annotation(comb_mat, add_numbers = TRUE),
        comb_order = order(-comb_size(comb_mat))) 
    })
    dev.off()
}

# create list object of transcripts for upset plot
toPlot_bin <- make_comb_mat(list(
            gCSI = gcsi_pval_b$pair,
            CCLE = ccle_pval_b$pair,
            GDSC = gdsc_pval_b$pair))

# plot upset plots
plot_upset(toPlot_bin, set_order = c("gCSI", "CCLE", "GDSC"), filename = "bin0")

############################################################
# Compare Stat of Overlapping P-Val Sig Biomarkers
############################################################

# function to create dataframe of overlapping p-value significant biomarkers
plot_overlapping <- function(gcsi_pval, ccle_pval, gdsc_pval) {

    # identify p-value significant biomarkers in >1 PSet
    common_biomarker <- c(intersect(gcsi_pval$pair, gdsc_pval$pair),
                        intersect(ccle_pval$pair, gdsc_pval$pair),
                        intersect(gcsi_pval$pair, ccle_pval$pair))

    # save results from individual pset associations
    gcsi_common <- gcsi_pval[gcsi_pval$pair %in% common_biomarker,]  
    ccle_common <- ccle_pval[ccle_pval$pair %in% common_biomarker,]    
    gdsc_common <- gdsc_pval[gdsc_pval$pair %in% common_biomarker,]    

    # add pset labels
    gcsi_common$PSet <- "gCSI"
    ccle_common$PSet <- "CCLE"
    gdsc_common$PSet <- "GDSC"

    # create dataframe for plotting
    toPlot <- rbind(gcsi_common, ccle_common, gdsc_common)
    toPlot$pair <- gsub("_", "\n", toPlot$pair)

    return(toPlot)
}

# get overlapping p-value significant biomarkers
toPlot_bin <- plot_overlapping(gcsi_pval_b, ccle_pval_b, gdsc_pval_b)

# plot overlapping biomarkers
png("../results/figures/figure9/common_bin0_pval_biomarkers.png", width = 17, height = 5, res = 600, units = "in")
ggplot(toPlot_bin, aes(x = PSet, y = W, fill = pval)) + geom_bar(stat="identity", color = "black") +
    facet_nested(~ factor(pair), scales = "free_x") +
    labs(fill = "P-Value", y = "Wilcoxon Rank Sum Test Statistic", x = "PSet") + 
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5), axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()



############################################################
# Stats of P-value Significant Biomarkers in other PSets
############################################################

# function to create dataframe of p-value significant biomarkers
get_biomarkers <- function(gcsi_pval, ccle_pval, gdsc_pval, gcsi_all, ccle_all, gdsc_all) {

    # add pset labels
    gcsi_all$PSet <- "gCSI"
    ccle_all$PSet <- "CCLE"
    gdsc_all$PSet <- "GDSC"

    # keep only top 10 associations from each P-Set
    if (colnames(gcsi_pval)[3] == "W") {
        gcsi_pval <- gcsi_pval[order(gcsi_pval$W, decreasing = T),][1:10,]
        ccle_pval <- ccle_pval[order(ccle_pval$W, decreasing = T),][1:10,]
        gdsc_pval <- gdsc_pval[order(gdsc_pval$W, decreasing = T),][1:10,]
    } else {
        gcsi_pval <- gcsi_pval[order(gcsi_pval$estimate, decreasing = T),][1:10,]
        ccle_pval <- ccle_pval[order(ccle_pval$estimate, decreasing = T),][1:10,]
        gdsc_pval <- gdsc_pval[order(gdsc_pval$estimate, decreasing = T),][1:10,]
    }

    pval_biomarkers <- unique(c(gcsi_pval$pair, ccle_pval$pair, gdsc_pval$pair))

    # create dataframe of selected 30 associations
    toPlot <- rbind(gcsi_all[gcsi_all$pair %in% pval_biomarkers,],
                    ccle_all[ccle_all$pair %in% pval_biomarkers,],
                    gdsc_all[gdsc_all$pair %in% pval_biomarkers,])

    # keep order of circRNAs
    toPlot$pair <- factor(toPlot$pair, levels = rev(pval_biomarkers))
    toPlot$PSet <- factor(toPlot$PSet, levels = c("gCSI", "CCLE", "GDSC"))

    return(toPlot)
}

toPlot_bin <- get_biomarkers(gcsi_pval_b, ccle_pval_b, gdsc_pval_b, gcsi_bin_dr, ccle_bin_dr, gdsc_bin_dr)

# plot overlapping biomarkers
png("../results/figures/figure9/top10_bin0_pval_biomarkers.png", width = 7, height = 5, res = 600, units = "in")
ggplot(toPlot_bin, aes(x = PSet, y = pair, fill = W)) + geom_tile() +
    geom_hline(yintercept = c(10.5, 20.5), linetype = "dashed") +
    geom_text(aes(label = ifelse(pval < 0.05, "*", ""))) +
    labs(fill = "Wilcoxon Rank\nSum Test Statistic", y = "Drug-CircRNA Pair", x = "PSet") + 
    theme_classic() + scale_fill_gradient2(low = 'white', mid = '#E1E7DF', high = '#878E76') +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5))
dev.off()


############################################################
# Pair-wise correlation plots across PSets
############################################################

# function to create correlation plot
corr_pset <- function(pset1, pset2, type, filename) {

    # get common pairs
    common <- intersect(pset1$pair, pset2$pair)
    pset1 <- pset1[pset1$pair %in% common,]
    pset2 <- pset2[pset2$pair %in% common,]

    # rename metric for plotting
    metric = "Linear Regression\nEffect Size"
    if (type == "bin") {
        colnames(pset1)[colnames(pset1) == "W"] <- colnames(pset2)[colnames(pset2) == "W"] <- "estimate"
        metric = "Wilcoxon Rank Sum\nTest Statistic"
    } 

    # order dataframes
    pset1 <- pset1[order(pset1$pair),]
    pset2 <- pset2[order(pset2$pair),]

    # create dataframe to plot
    toPlot <- data.frame(pair = pset1$pair, 
                         estimate1 = pset1$estimate, estimate2 = pset2$estimate,
                         pval1 = pset1$pval, pval2 = pset2$pval)
    toPlot <- na.omit(toPlot)

    # create pset labels
    p1 <- str_split_1(filename, pattern = "_")[1]
    p2 <- str_split_1(filename, pattern = "_")[2]

    # create pvalue significance label
    toPlot$pval_sig <- ifelse(toPlot$pval1 < 0.05, 
                            ifelse(toPlot$pval2 < 0.05, "Both PSets", paste(p1, "Only")),
                        ifelse(toPlot$pval2 < 0.05, paste(p2, "Only"), "Neither PSet"))
    toPlot$pval_sig <- factor(toPlot$pval_sig, levels = c("Both PSets", paste(p1, "Only"), paste(p2, "Only"), "Neither PSet"))

    # set palette for plotting
    pal = c("#FFE573", "#62929E", "#CB807D", "grey")

    # plot scatter plot
    png(paste0(filename, ".png"), width = 6, height = 4, res = 600, units = "in")
    print({ggplot() + 
        geom_point(data = toPlot, aes(x = estimate1, y = estimate2, color = pval_sig), shape = 16) +
        geom_point(data = toPlot[toPlot$pval_sig == paste(p1, "Only"),], aes(x = estimate1, y = estimate2), size = 2, shape = 21, fill = pal[2]) +
        geom_point(data = toPlot[toPlot$pval_sig == paste(p2, "Only"),], aes(x = estimate1, y = estimate2), size = 2, shape = 21, fill = pal[3]) +
        geom_point(data = toPlot[toPlot$pval_sig == "Both PSets", ], aes(x = estimate1, y = estimate2), size = 2.5, shape = 21, fill = pal[1]) +
        scale_color_manual("P-Value < 0.05", values = pal) +
        theme_classic() + guides(color = guide_legend(override.aes = list(size = 3))) +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
        labs(x = paste(p1, metric), y = paste(p2, metric))})
    dev.off()
}


# binarized dr
gcsi_ccle_b <- corr_pset(gcsi_bin_dr, ccle_bin_dr, "bin", "gCSI_CCLE_bin0")
gcsi_gdsc_b <- corr_pset(gcsi_bin_dr, gdsc_bin_dr, "bin", "gCSI_GDSC_bin0")
ccle_gdsc_b <- corr_pset(ccle_bin_dr, gdsc_bin_dr, "bin", "CCLE_GDSC_bin0")

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


############################################################
# Specify analysis
############################################################

# analysis = c("circ", "GE")  
# circ for circRNA counts, GE for circ mapped to GE
analysis = "circ"

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
gcsi_df <- fread(paste0(path, "gcsi_counts.tsv"), data.table = F)
ccle_df <- fread(paste0(path, "ccle_counts.tsv"), data.table = F)
gdsc_df <- fread(paste0(path, "gdsc_counts.tsv"), data.table = F)

# formating
gcsi_df[is.na(gcsi_df)] <- 0
ccle_df[is.na(ccle_df)] <- 0
gdsc_df[is.na(gdsc_df)] <- 0

rownames(gcsi_df) <- gcsi_df$V1
rownames(ccle_df) <- ccle_df$V1
rownames(gdsc_df) <- gdsc_df$V1

gcsi_df$V1 <- ccle_df$V1 <- gdsc_df$V1 <- NULL      

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
# Filter expression dataframes
############################################################

# remove features with exp in less than threshold samples
dim(gcsi_df)        # 672 1766
dim(ccle_df)        # 1019 4437
dim(gdsc_df)        # 443 1491

# THRESHOLDS FOR circ (all samples):
# thres = 10, gcsi_df: 602, ccle_df: 1517, gdsc_df: 614 
# thres = 15, gcsi_df: 371, ccle_df: 975, gdsc_df: 444
# thres = 20, gcsi_df: 268, ccle_df: 716, gdsc_df: 354

dim(gcsi_df[,colnames(gcsi_df)[colSums(gcsi_df != 0) >= thres], drop = FALSE])
dim(ccle_df[,colnames(ccle_df)[colSums(ccle_df != 0) >= thres], drop = FALSE])
dim(gdsc_df[,colnames(gdsc_df)[colSums(gdsc_df != 0) >= thres], drop = FALSE])

gcsi_df <- gcsi_df[,colnames(gcsi_df)[colSums(gcsi_df != 0) >= thres], drop = FALSE]
ccle_df <- ccle_df[,colnames(ccle_df)[colSums(ccle_df != 0) >= thres], drop = FALSE]
gdsc_df <- gdsc_df[,colnames(gdsc_df)[colSums(gdsc_df != 0) >= thres], drop = FALSE]

# zero to NA for downstream analysis
gcsi_df[gcsi_df == 0] <- NA
ccle_df[ccle_df == 0] <- NA
gdsc_df[gdsc_df == 0] <- NA

save(gcsi_df, ccle_df, gdsc_df, file = "../results/data/biomarker_analysis_circ.RData")


############################################################
# Binarize transcript expression by zero expression
############################################################

# function to compute drug response associations
binary_dr <- function(counts_df, drug_df) {

    features <- colnames(counts_df)

    # create data frame to hold results
    combinations <- expand.grid(Drug = rownames(drug_df), Feature = features)
    combinations$num_samples <- combinations$num_high <- combinations$diff <- combinations$pval <-  combinations$W <- NA
    #combinations$drug <- combinations$feature <- NA
    
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

                # difference in average AAC of two groups
                AAC_diff <- mean(high[!is.na(high)]) - mean(low[!is.na(low)])

                # save results
                combinations$W[row] <- res$statistic
                combinations$pval[row] <- res$p.value
                combinations$diff[row] <- AAC_diff
                combinations$num_samples[row] <- nrow(subset)
                combinations$num_high[row] <- length(high_exp)
                #combinations$feature[row] <- feature
                #combinations$drug[row] <- rownames(drug_df)[j]

            } else {
                # save results
                combinations$W[row] <- NA
                combinations$pval[row] <- NA
                combinations$diff[row] <- NA
                combinations$num_samples[row] <- nrow(subset)
                combinations$num_high[row] <- NA
                #combinations$feature[row] <- feature
                #combinations$drug[row] <- rownames(drug_df)[j]
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
ccle_bin_dr <- binary_dr(ccle_df, ctrp_sen)
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

load(paste0(dr_out, "gcsi_bin0.RData"))
load(paste0(dr_out, "ccle_bin0.RData"))
load(paste0(dr_out, "gdsc_bin0.RData"))

############################################################
# FDR normalize per drug
############################################################

gcsi_bin_dr$FDR_drug <- ave(gcsi_bin_dr$pval, gcsi_bin_dr$Drug, FUN = function(p) p.adjust(p, method = "BH"))
ccle_bin_dr$FDR_drug <- ave(ccle_bin_dr$pval, ccle_bin_dr$Drug, FUN = function(p) p.adjust(p, method = "BH"))
gdsc_bin_dr$FDR_drug <- ave(gdsc_bin_dr$pval, gdsc_bin_dr$Drug, FUN = function(p) p.adjust(p, method = "BH"))


############################################################
# Subset for significant associations
############################################################

# function to subset for significant associations
subset_dr <- function(dr) {
    dr <- dr[which(dr$FDR_drug < 0.1),] 
    return(dr)
}

gcsi_fdr_b <- subset_dr(gcsi_bin_dr)
gdsc_fdr_b <- subset_dr(gdsc_bin_dr)
ccle_fdr_b <- subset_dr(ccle_bin_dr)

############################################################
# Save Table 4
############################################################

gcsi_fdr_b$dataset <- "gCSI"
ccle_fdr_b$dataset <- "CCLE"
gdsc_fdr_b$dataset <- "GDSC2"

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

# function to create upset plot
plot_upset <- function(comb_mat, set_order, filename) {
    #png(paste0("../results/figures/figure9/upset_dr_", filename, ".png"), width = 8, height = 3, res = 600, units = "in")
    png(paste0("upset_dr_", filename, ".png"), width = 8, height = 3, res = 600, units = "in")
    print({UpSet(comb_mat, set_order = set_order,
        top_annotation = upset_top_annotation(comb_mat, add_numbers = TRUE),
        comb_order = order(-comb_size(comb_mat)),
        right_annotation = upset_right_annotation(comb_mat, add_numbers = TRUE)) 
    })
    dev.off()
}

# create list object of transcripts for upset plot
toPlot_bin <- make_comb_mat(list(
            gCSI = gcsi_fdr_b$pair,
            CCLE = ccle_fdr_b$pair,
            GDSC = gdsc_fdr_b$pair))

# plot upset plots
plot_upset(toPlot_bin, set_order = c("gCSI", "CCLE", "GDSC"), filename = "bin0")


############################################################
# Stats of P-value Significant Biomarkers in other PSets
############################################################

# save top 10 associations from each PSet
gcsi_top <- gcsi_fdr_b[order(abs(gcsi_fdr_b$diff), decreasing = T),][1:20,]
ccle_top <- ccle_fdr_b[order(abs(ccle_fdr_b$diff), decreasing = T),][1:20,]
gdsc_top <- gdsc_fdr_b[order(abs(gdsc_fdr_b$diff), decreasing = T),][1:20,]

top_biomarkers <- unique(c(gcsi_top$pair, ccle_top$pair, gdsc_top$pair))

# function to create dataframe for plotting
format_bin <- function(bin_df, counts_df, drug_df) {

    # create dataframe for plotting
    df <- data.frame(Feature = gsub(".*_", "", top_biomarkers),
                     Drug = gsub("_.*", "", top_biomarkers), 
                     Pair = top_biomarkers,
                     PSet = bin_df$PSet[1], 
                     Status = NA, 
                     Diff = NA)

    # subset results
    subset <- bin_df[bin_df$pair %in% top_biomarkers,]

    for (i in seq_along(df$Pair)) {
        drug = df$Drug[i]
        feature = df$Feature[i]
        pair = df$Pair[i]

        if (!drug %in% rownames(drug_df)) { # drug no in PSet
            df$Status[i] <- NA
        } else if (!feature %in% colnames(counts_df)) { # circRNA not detected
            df$Status[i] <- "ND"
        } else if (!pair %in% subset$pair) { # not enough for association
            df$Diff[i] <- 0
            df$Status[i] <- "NE"
        } else {

            df$Diff[i] <- subset[subset$pair == pair,]$diff   # P-Val > 0.05
            df$Status[i] <- ""

            # check p-value and FDR
            if (subset[subset$pair == pair,]$FDR_drug < 0.1) { #FDR < 0.1
                df$Status[i] <- "*"
            }
            if (subset[subset$pair == pair,]$FDR_drug < 0.05) { #FDR < 0.05
                df$Status[i] <- "**"
            }
        }
    }
    return(df)
}

# add pset labels
gcsi_bin_dr$PSet <- "gCSI"
ccle_bin_dr$PSet <- "CCLE"
gdsc_bin_dr$PSet <- "GDSC"

toPlot <- rbind(format_bin(gcsi_bin_dr, gcsi_df, gcsi_sen), 
                format_bin(ccle_bin_dr, ccle_df, ctrp_sen), 
                format_bin(gdsc_bin_dr, gdsc_df, gdsc_sen))
toPlot$PSet <- factor(toPlot$PSet, levels = c("gCSI", "CCLE", "GDSC"))

# order by drug per PSet
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
png("../results/figures/figure9/top10_bin0_biomarkers.png", width = 6, height = 10, res = 600, units = "in")
ggplot(toPlot, aes(x = PSet, y = Pair, fill = Diff)) + 
    geom_tile() +
    geom_hline(yintercept = c(20.5, 40.5), linetype = "dashed") +
    geom_text(aes(label = Status), size = 3) +
    scale_y_discrete(labels = setNames(toPlot$label, toPlot$Pair)) +
    labs(fill = "Δ Mean AAC \n(Expression vs\nNon-expression)", y = "Drug-CircRNA Pair", x = "PSet") + 
    theme_classic() + 
    scale_fill_gradient2(low = "#9D3737", mid = "white", high = "#3670A0", midpoint = 0, na.value = "grey") +
    theme(
        panel.border = element_rect(color = "black", fill = NA, size = 0.5), 
        legend.title = element_text(size = 9))
dev.off()

############################################################
# Visualize Vorinostat_chr8.61680967.61684188
############################################################

circ <- "chr8.61680967.61684188"
drug <- "Vorinostat"

# helper function to label high vs low circ exp
label_exp <- function(counts_df, drug_df, dataset) {

    # binarize transcript expression by median
    subset <- counts_df[,which(colnames(counts_df) == circ), drop = FALSE]
    high_exp <- rownames(subset[complete.cases(subset), , drop = FALSE])
    low_exp <- rownames(subset)[-which(rownames(subset) %in% high_exp)]

    # get AAC for each expression group
    high = drug_df[rownames(drug_df) == drug,colnames(drug_df) %in% high_exp] |> as.numeric()
    low = drug_df[rownames(drug_df) == drug,colnames(drug_df) %in% low_exp] |> as.numeric()

    # return dataframe
    df <- data.frame(AAC = c(high, low),
                     Label = c(rep("Expression", length(high)), rep("Non-Expression", length(low))))
    df$Dataset <- dataset
    return(df)
}

toPlot <- rbind(label_exp(gcsi_df, gcsi_sen, "gCSI"),
                label_exp(ccle_df, ctrp_sen, "CCLE"),
                label_exp(gdsc_df, gdsc_sen, "GDSC"))

# format dataframe
toPlot$Dataset <- factor(toPlot$Dataset, levels = c("gCSI", "CCLE", "GDSC"))

# plot AAC by circ expression
png("../results/figures/figure9/example_biomarker.png", width = 7, height = 5.5, res = 600, units = "in")
ggplot(toPlot, aes(x = Dataset, y = AAC, fill = Label)) + 
    geom_boxplot() +
    theme_classic() + 
    scale_fill_manual(values = c("#586994", "#937D64")) +
    geom_signif(
        y_position = c(0.47, 0.49, 0.54), 
        xmin = c(0.8, 1.8, 2.8), xmax = c(1.2, 2.2, 3.2),
        annotation = c("*", "**", "**"),
        tip_length = 0.01,
        textsize = 4) +
    ylim(c(-0.01, 0.6)) +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
          legend.key.size = unit(0.5, 'cm'),
          legend.title = element_text(size = 9)) +
    labs(y = "Vorinostat Response (AAC)", fill = "Expression of\nchr8.61680967.61684188")
dev.off()

############################################################
# Pair-wise correlation plots across PSets
############################################################

# function to create correlation plot
corr_pset <- function(pset1, pset2, filename) {

    # get common pairs
    common <- intersect(pset1$pair, pset2$pair)
    pset1 <- pset1[pset1$pair %in% common,]
    pset2 <- pset2[pset2$pair %in% common,]

    # print number of common pairs
    print(paste("Number of common pairs:", length(common)))

    # order dataframes
    pset1 <- pset1[order(pset1$pair),]
    pset2 <- pset2[order(pset2$pair),]

    # create dataframe to plot
    toPlot <- data.frame(pair = pset1$pair, 
                         diff1 = pset1$diff, diff2 = pset2$diff,
                         fdr1 = pset1$FDR_drug, fdr2 = pset2$FDR_drug)
    toPlot <- na.omit(toPlot)

    # compute spearman correlation
    res <- cor.test(toPlot$diff1, toPlot$diff2, method = "spearman", alternative = "two.sided")
    print(res)

    # compute spearman correlation for those with FDR < 1 in at least one dataset
    to_keep <- c(toPlot$pair[toPlot$fdr1 < 0.1], toPlot$pair[toPlot$fdr2 < 0.1]) |> unique()
    res <- cor.test(toPlot$diff1[toPlot$pair %in% to_keep], toPlot$diff2[toPlot$pair %in% to_keep], method = "spearman", alternative = "two.sided")
    print(res)

    # create pset labels
    p1 <- gsub("../results/figures/figure9/", "", str_split_1(filename, pattern = "_")[1])
    p2 <- str_split_1(filename, pattern = "_")[2]

    # create fdr significance label
    toPlot$fdr_sig <- ifelse(toPlot$fdr1 < 0.05, 
                            ifelse(toPlot$fdr2 < 0.05, "Both Datasets", paste(p1, "Only")),
                        ifelse(toPlot$fdr2 < 0.05, paste(p2, "Only"), "Neither Dataset"))
    toPlot$fdr_sig <- factor(toPlot$fdr_sig, levels = c("Both Datasets", paste(p1, "Only"), paste(p2, "Only"), "Neither Dataset"))

    # set palette for plotting
    group_names <- c("Both Datasets", paste(p1, "Only"), paste(p2, "Only"), "Neither Dataset")
    col <- c("#FCD0A1", "#63535B", "#53917E", "grey")
    pal <- setNames(col, group_names)
    
    # plot scatter plot
    png(paste0(filename, ".png"), width = 6, height = 4, res = 600, units = "in")
    print({
        ggplot(toPlot, aes(x = diff1, y = diff2)) + 
        geom_point(data = toPlot, aes(x = diff1, y = diff2, color = fdr_sig), shape = 16) +
        geom_point(data = toPlot[toPlot$fdr_sig == paste(p1, "Only"),], aes(x = diff1, y = diff2), size = 2, shape = 21, fill = pal[2]) +
        geom_point(data = toPlot[toPlot$fdr_sig == paste(p2, "Only"),], aes(x = diff1, y = diff2), size = 2, shape = 21, fill = pal[3]) +
        geom_point(data = toPlot[toPlot$fdr_sig == "Both Datasets", ], aes(x = diff1, y = diff2), size = 2.5, shape = 21, fill = pal[1]) +
        geom_text_repel(data = toPlot[toPlot$pair == paste(drug, circ, sep = "_"), ],
                        aes(x = diff1, y = diff2, label = paste(drug, circ, sep = "\n")),
                        force_pull = 0.1, 
                        nudge_y = -0.19,
                        nudge_x = -0.11,
                        box.padding = 0.4) +
        geom_vline(xintercept = 0, linetype = "dashed") +
        geom_hline(yintercept = 0, linetype = "dashed") +
        scale_color_manual("FDR < 0.1", values = pal) +
        theme_classic() + guides(color = guide_legend(override.aes = list(size = 3))) +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
        labs(x = paste("Δ Mean AAC in", p1), y = paste("Δ Mean AAC in", p2))
        })
    dev.off()
}


# binarized dr
corr_pset(gcsi_bin_dr, ccle_bin_dr, "../results/figures/figure9/gCSI_CCLE_bin0")
corr_pset(gcsi_bin_dr, gdsc_bin_dr, "../results/figures/figure9/gCSI_GDSC_bin0")
corr_pset(ccle_bin_dr, gdsc_bin_dr, "../results/figures/figure9/CCLE_GDSC_bin0")


############################################################
# Plot distribution of 0 expression
############################################################

# function to plot the distribution of zero expression cross cell lines
plot_zeros <- function(bin, filename) {

    # format dataframe for plotting
    bin <- bin[,colnames(bin) %in% c("Feature", "num_high", "num_samples")] |> unique()
    bin <- bin[order(bin$num_high, decreasing = T),]
    toPlot <- melt(bin)

    # get circRNAs order by proportion of >0 expression
    circs <- bin$Feature
    toPlot$Feature <- factor(toPlot$Feature, levels = unique(circs))

    png(paste0("../results/figures/figure9/dist_zeros_", filename, ".png"), width = 4, height = 5, res = 600, units = "in")
    print({ggplot(toPlot, aes(y = Feature, x = value, fill = variable)) + 
        geom_bar(position="dodge", stat="identity") +
        scale_fill_manual("Expression", values = c("#60A090", "#DDDDDD"), label = c(">0", "0")) +
        theme(axis.text.y=element_blank(),
              axis.ticks.y=element_blank()) +
        labs(x = "Number of Cell Lines", y = "CircRNA Transcript")
    })
    dev.off()
}

plot_zeros(gcsi_bin_dr, "gcsi_bin0")
plot_zeros(ccle_bin_dr, "ccle_bin0")
plot_zeros(gdsc_bin_dr, "gdsc_bin0")
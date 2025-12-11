#' Function to compute drug response associations
#' of binarized transcript expression
#' 
#' @param counts_df dataframe. Counts matrix
#' @param drug_df dataframe. Drug sensitivity matrix
#' @param label str. {Pipeline}_{PSet} (only for MYC)
#' 
binary_dr <- function(counts_df, drug_df, label = NA) {

    features <- colnames(counts_df)

    # create data frame to hold results
    combinations <- expand.grid(Drug = rownames(drug_df), Feature = features)
    combinations$num_samples <- combinations$num_high <- combinations$diff <- combinations$pval <-  combinations$W <- NA
    combinations$drug <- combinations$feature <- NA # sanity check
    
    # initiate row count
    row <- 1
    print(paste("Number of features:", length(features)))

    for (i in seq_along(features)) {
        
        feature <- features[i]

        # binarize transcript expression by median (0)
        subset <- counts_df[,i, drop = FALSE]
        if (is.na(label)) subset[subset == 0] <- NA

        high_exp <- rownames(subset[complete.cases(subset), , drop = FALSE])
        low_exp <- rownames(subset)[-which(rownames(subset) %in% high_exp)]

        # loop through each drug

        for (j in seq_along(rownames(drug_df))) {
            
            high <- drug_df[j,colnames(drug_df) %in% high_exp] |> as.numeric()
            low <- drug_df[j,colnames(drug_df) %in% low_exp] |> as.numeric()

            if (length(high[!is.na(high)]) > 1 && length(low[!is.na(low)]) > 1) {

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
    if (is.na(label)) combinations$FDR_drug <- ave(combinations$pval, combinations$Drug, FUN = function(p) p.adjust(p, method = "BH"))

    # format dataframe for plotting 
    combinations <- combinations[order(combinations$W, decreasing = TRUE),]
    combinations$rank <- seq_len(nrow(combinations))
    combinations$pair <- paste(combinations$Drug, combinations$Feature, sep = "_")
    if (!is.na(label)) combinations$label <- label
    return(combinations)
}

#' Function to format binarized drug response associations
#' across datasets for plotting
#' 
#' Annotates transcripts and drugs based on expression / screen
#' 
#' @param bin_df dataframe. Drug response associations
#' @param counts_df data.frame. Expression matrix
#' @param drug_df data.frame. Drug sensitivity matrix
#' 
format_bin <- function(bin_df, counts_df, drug_df) {

    # create dataframe for plotting
    df <- data.frame(Feature = gsub(".*_", "", top_biomarkers),
                     Drug = gsub("_.*", "", top_biomarkers), 
                     Pair = top_biomarkers,
                     dataset = bin_df$dataset[1], 
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

#' Function to create dataframe for plotting MYC top biomarkers
#' Adapted from format_bin() above
#' 
formatMYC <- function(drug_df, label) {

    # subset bin results
    subset <- toPlot[toPlot$label == label,]

    # create dataframe for plotting
    df <- data.frame(Drug = drugs, Label = label, Status = NA, Diff = NA)

    for (i in seq_along(drugs)) {
        drug = df$Drug[i]

        if (!drug %in% rownames(drug_df)) { # drug no in PSet
            df$Status[i] <- NA
        } else if (!drug %in% subset$Drug) {
            df$Diff[i] <- 0
            df$Status[i] <- "NE"
        } else {

            df$Diff[i] <- subset[subset$Drug == drug,]$diff   # P-Val > 0.05
            df$Status[i] <- ""

            # check p-value and FDR
            if (subset[subset$Drug == drug,]$pval < 0.05) { #P-Val < 0.05
                df$Status[i] <- "*"
            }
            if (subset[subset$Drug == drug,]$FDR < 0.1) { #FDR < 0.05
                df$Status[i] <- "**"
            }
        }
    }
    return(df)
}


#' Function to label high vs low circ exp
#' for Vorinostat_chr8.61680967.61684188
#' 
label_exp <- function(counts_df, drug_df, dataset) {

    circ <- "chr8.61680967.61684188"
    drug <- "Vorinostat"

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


#' Function to create correlation plot
#' 
plot_corr_pset <- function(pset1, pset2, filename) {

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
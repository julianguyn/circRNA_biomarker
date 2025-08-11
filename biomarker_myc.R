# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggrepel)
})


############################################################
# Load in circRNA expression data
############################################################

path = "../data/processed_cellline/all_samples/"

# load circRNA expression data
ciri_gcsi <- fread(paste0(path, "CIRI2/ciri_gcsi_counts.tsv"), data.table = F)
ciri_gdsc <- fread(paste0(path, "CIRI2/ciri_gdsc_counts.tsv"), data.table = F)
ciri_ccle <- fread(paste0(path, "CIRI2/ciri_ccle_counts.tsv"), data.table = F)

circ_gcsi <- fread(paste0(path, "CIRCexplorer2/circ_gcsi_counts.tsv"), data.table = F)
circ_gdsc <- fread(paste0(path, "CIRCexplorer2/circ_gdsc_counts.tsv"), data.table = F)
circ_ccle <- fread(paste0(path, "CIRCexplorer2/circ_ccle_counts.tsv"), data.table = F)

cfnd_gcsi <- fread(paste0(path, "circRNA_finder/cfnd_gcsi_counts.tsv"), data.table = F)
cfnd_gdsc <- fread(paste0(path, "circRNA_finder/cfnd_gdsc_counts.tsv"), data.table = F)
cfnd_ccle <- fread(paste0(path, "circRNA_finder/cfnd_ccle_counts.tsv"), data.table = F)

fcrc_gcsi <- fread(paste0(path, "find_circ/fcrc_gcsi_counts.tsv"), data.table = F)
fcrc_gdsc <- fread(paste0(path, "find_circ/fcrc_gdsc_counts.tsv"), data.table = F)
fcrc_ccle <- fread(paste0(path, "find_circ/fcrc_ccle_counts.tsv"), data.table = F)


############################################################
# Identify circRNA mapping to MYC
############################################################

# save all circRNAs on chr8
all_circ <- c(colnames(ciri_gcsi), colnames(ciri_gdsc), colnames(ciri_ccle),
              colnames(circ_gcsi), colnames(circ_gdsc), colnames(circ_ccle),
              colnames(cfnd_gcsi), colnames(cfnd_gdsc), colnames(cfnd_ccle),
              colnames(fcrc_gcsi), colnames(fcrc_gdsc), colnames(fcrc_ccle))
chr8 <- gsub("chr8\\.", "", all_circ[grep("chr8", all_circ)])
split_vals <- strsplit(chr8, "\\.")
df <- data.frame(
    circ = paste0("chr8.", chr8),
    start = as.numeric(sapply(split_vals, `[`, 1)),
    end = as.numeric(sapply(split_vals, `[`, 2)),
    myc = 0
)

# get MYC bases (MYC gene: 8:127735434-127742951)
myc_start = 127735434
myc_end = 127742951
myc <- myc_start:myc_end

# save number of overlapping bases per circRNA on chr8
for (i in seq_along(df$circ)) {
    circ <- df$start[i]:df$end[i]
    if (sum(circ %in% myc) > 0) {
        df$myc[i] <- sum(circ %in% myc)
    }
}

# identify circRNAs that overlapp 100% with the MYC gene 
df$total <- df$end - df$start + 1
df$prop <- df$myc / df$total * 100

# plot proportion distribution
png("../results/figures/figure9/myc/circ_proportion_overlap.png", width = 7.5, height = 6, res = 600, units = "in")
ggplot(df, aes(x = prop)) + 
    geom_histogram(color="black", fill="#C57B57") + 
    geom_text(stat = "bin", aes(label = after_stat(count)), vjust = -1) +
    theme_classic() + 
    labs(x = "Proportion Overlap of CircRNA Transcript and MYC", y = "Count")
dev.off()

############################################################
# Keep circRNA overlapping 100% to MYC
############################################################

df <- df[df$prop == 100,] #n=508

# define function to filter circ exp dataframes
filter_myc <- function(circ) {

    rownames(circ) <- circ$sample
    circ <- circ[,colnames(circ) %in% df$circ]
    print(ncol(circ))
    circ <- data.frame(MYC = rowMeans(circ))

    return(circ)
}

ciri_gcsi <- filter_myc(ciri_gcsi) #8
ciri_gdsc <- filter_myc(ciri_gdsc) #2
ciri_ccle <- filter_myc(ciri_ccle) #9

circ_gcsi <- filter_myc(circ_gcsi) #85
circ_gdsc <- filter_myc(circ_gdsc) #3
circ_ccle <- filter_myc(circ_ccle) #19

cfnd_gcsi <- filter_myc(cfnd_gcsi) #116
cfnd_gdsc <- filter_myc(cfnd_gdsc) #3
cfnd_ccle <- filter_myc(cfnd_ccle) #155

fcrc_gcsi <- filter_myc(fcrc_gcsi) #41
fcrc_gdsc <- filter_myc(fcrc_gdsc) #11
fcrc_ccle <- filter_myc(fcrc_ccle) #56

# save dataframes
save(ciri_gcsi, ciri_gdsc, ciri_ccle,
     circ_gcsi, circ_gdsc, circ_ccle,
     cfnd_gcsi, cfnd_gdsc, cfnd_ccle,
     fcrc_gcsi, fcrc_gdsc, fcrc_ccle,
     file = "../results/data/circ_myc.RData")


############################################################
# Load in drug response data and subset
############################################################

# load in drug sensitivity from PSets
load("../data/temp/sensitivity_data.RData")

############################################################
# Binarize transcript expression by zero expression
############################################################

# function to compute drug response associations
binary_dr <- function(counts_df, drug_df, label) {

    features <- colnames(counts_df)

    # create data frame to hold results
    combinations <- expand.grid(Drug = rownames(drug_df), Feature = features)
    combinations$num_samples <- combinations$num_high <- combinations$diff <- combinations$pval <-  combinations$W <- NA
    combinations$drug <- combinations$feature <- NA
    
    # initiate row count
    row = 1
    print(paste("Number of features:", length(features)))

    for (i in seq_along(features)) {
        
        print(i)
        feature <- features[i]

        # binarize transcript expression by median (0)
        subset <- counts_df[,i, drop = FALSE]
        subset[subset == 0] <- NA

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
                combinations$feature[row] <- feature
                combinations$drug[row] <- rownames(drug_df)[j]

            } else {
                # save results
                combinations$W[row] <- NA
                combinations$pval[row] <- NA
                combinations$diff[row] <- NA
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
    combinations$label <- label
    return(combinations)
}


ciri_gcsi_bin <- binary_dr(ciri_gcsi, gcsi_sen, "CIRI_gCSI")
ciri_gdsc_bin <- binary_dr(ciri_gdsc, gdsc_sen, "CIRI_GDSC")
ciri_ccle_bin <- binary_dr(ciri_ccle, ccle_sen, "CIRI_CCLE")

circ_gcsi_bin <- binary_dr(circ_gcsi, gcsi_sen, "CIRC_gCSI")
circ_gdsc_bin <- binary_dr(circ_gdsc, gdsc_sen, "CIRC_GDSC")
circ_ccle_bin <- binary_dr(circ_ccle, ccle_sen, "CIRC_CCLE")

cfnd_gcsi_bin <- binary_dr(cfnd_gcsi, gcsi_sen, "CFND_gCSI")
cfnd_gdsc_bin <- binary_dr(cfnd_gdsc, gdsc_sen, "CFND_GDSC")
cfnd_ccle_bin <- binary_dr(cfnd_ccle, ccle_sen, "CFND_CCLE")

fcrc_gcsi_bin <- binary_dr(fcrc_gcsi, gcsi_sen, "FCRC_gCSI")
fcrc_gdsc_bin <- binary_dr(fcrc_gdsc, gdsc_sen, "FCRC_GDSC")
fcrc_ccle_bin <- binary_dr(fcrc_ccle, ccle_sen, "FCRC_CCLE")

# save dataframes
save(ciri_gcsi_bin, ciri_gdsc_bin, ciri_ccle_bin,
     circ_gcsi_bin, circ_gdsc_bin, circ_ccle_bin,
     cfnd_gcsi_bin, cfnd_gdsc_bin, cfnd_ccle_bin,
     fcrc_gcsi_bin, fcrc_gdsc_bin, fcrc_ccle_bin,
     file = "../results/data/circ_myc_avg_bindr.RData")
#     file = "../results/data/circ_myc_bindr.RData")



############################################################
# Volcano Plots
############################################################

set.seed(101)

plot_volcano <- function(bin_dr) {

    # get x limits
    x <- max(c(abs(min(bin_dr$diff)), max(bin_dr$diff)))

    # plot
    p <- ggplot() +
        geom_point(data = bin_dr, aes(x = diff, y = -log(FDR)), color = "gray") +
        geom_point(data = bin_dr[!bin_dr$to_label == "", ],
                    aes(x = diff, y = -log(FDR), color = label)) +
        geom_text_repel(data = bin_dr[!bin_dr$to_label == "", ],
                        aes(x = diff, y = -log(FDR), label = to_label),
                        force_pull = 0.5, max.overlaps = Inf, force = 5,
                        box.padding = 0.4) +
        scale_color_manual(
            values = c("CIRC_gCSI" = "#8181DA", "FCRC_gCSI" = "#271D80", "CFND_CCLE" = "#B23A48", "FCRC_CCLE" = "#6F1E27", "FCRC_GDSC" = "#496F5D"), 
            labels = c("CIRC_gCSI" = "gCSI (CIRCexplorer2)", "FCRC_gCSI" = "gCSI (find_circ)", "CFND_CCLE" = "CCLE (circRNA_finder)", "FCRC_CCLE" = "CCLE (find_circ)", "FCRC_GDSC" = "GDSC2 (find_circ)"),
            name = "PSet (Pipeline)") +
        theme_classic() + xlim(-x, x) +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                legend.key.size = unit(0.5, 'cm')) +
        labs(x = "Difference in Average AAC", y = "-log(FDR)")
    return(p)

}

toPlot <- rbind(ciri_gcsi_bin, ciri_gdsc_bin, ciri_ccle_bin,
                circ_gcsi_bin, circ_gdsc_bin, circ_ccle_bin,
                cfnd_gcsi_bin, cfnd_gdsc_bin, cfnd_ccle_bin,
                fcrc_gcsi_bin, fcrc_gdsc_bin, fcrc_ccle_bin)

toPlot <- toPlot[!is.na(toPlot$diff),]

# label overlapping associations
toPlot$temp <- paste(toPlot$pair, toPlot$label, sep = "_")
to_keep <- toPlot[toPlot$FDR < 0.05,]$temp
toPlot$to_label <- ifelse(toPlot$temp %in% to_keep, as.character(toPlot$Drug), "")

# keep only associations with magnitude difference > 0.03
toPlot[abs(toPlot$diff) < 0.03,]$to_label <- ""

# label by significance
toPlot$Label <- ifelse(toPlot$FDR < 0.05, "FDR\nSignificant", "Not FDR\nSignificant")
toPlot$Label <- factor(toPlot$Label, levels = c("FDR\nSignificant", "Not FDR\nSignificant"))
toPlot$label <- factor(toPlot$label, levels = c("CIRI_gCSI", "CIRC_gCSI", "CFND_gCSI", "FCRC_gCSI",
                                                "CIRI_CCLE", "CIRC_CCLE", "CFND_CCLE", "FCRC_CCLE",
                                                "CIRI_GDSC", "CIRC_GDSC", "CFND_GDSC", "FCRC_GDSC"))


#png("../results/figures/figure9/myc/volcano.png", width = 6, height = 4, res = 600, units = "in")
png("../results/figures/figure9/myc/myc_avg_volcano.png", width = 7, height = 5, res = 600, units = "in")
plot_volcano(toPlot)
dev.off()


############################################################
# Stats of MYC drug associations in other PSets/Pipelines
############################################################

# save significant MYC drug associations
drugs <- toPlot[!toPlot$to_label == "",]$Drug |> as.character() |> unique()

# function to create dataframe for plotting
formatMYC <- function(drug_df, label) {

    # subset bin results
    subset <- toPlot[toPlot$label == label,]

    # create dataframe for plotting
    df <- data.frame(Drug = drugs, Label = label, Status = NA, W = NA)

    for (i in seq_along(drugs)) {
        drug = df$Drug[i]

        if (!drug %in% rownames(drug_df)) {
            df$Status[i] <- ""
        } else if (!drug %in% subset$Drug) {
            df$W[i] <- 0
            df$Status[i] <- "NA"
        } else {

            df$W[i] <- subset[subset$Drug == drug,]$W   # P-Val > 0.05
            df$Status[i] <- ""

            # check p-value and FDR
            if (subset[subset$Drug == drug,]$pval < 0.05) { #P-Val < 0.05
                df$Status[i] <- "*"
            }
            if (subset[subset$Drug == drug,]$FDR < 0.05) { #FDR < 0.05
                df$Status[i] <- "**"
            }
        }
    }
    return(df)
}

sig_drug <- rbind(formatMYC(gcsi_sen, "CIRI_gCSI"), formatMYC(ccle_sen, "CIRI_CCLE"), formatMYC(gdsc_sen, "CIRI_GDSC"),
                  formatMYC(gcsi_sen, "CIRC_gCSI"), formatMYC(ccle_sen, "CIRC_CCLE"), formatMYC(gdsc_sen, "CIRC_GDSC"),
                  formatMYC(gcsi_sen, "CFND_gCSI"), formatMYC(ccle_sen, "CFND_CCLE"), formatMYC(gdsc_sen, "CFND_GDSC"),
                  formatMYC(gcsi_sen, "FCRC_gCSI"), formatMYC(ccle_sen, "FCRC_CCLE"), formatMYC(gdsc_sen, "FCRC_GDSC"))
sig_drug$PSet <- gsub(".*_", "", sig_drug$Label)
sig_drug$Pipeline <- rep(c("CIRI2", "CIRCexplorer2", "circRNA_finder", "find_circ"), each = 36)
sig_drug$Pipeline <- factor(sig_drug$Pipeline, levels = c("CIRI2", "CIRCexplorer2", "circRNA_finder", "find_circ"))
sig_drug$PSet <- factor(sig_drug$PSet, levels = c("gCSI", "CCLE", "GDSC"))


# plot overlapping biomarkers
png("myc_overlap.png", width = 6, height = 7, res = 600, units = "in")
ggplot(sig_drug, aes(x = PSet, y = Drug, fill = W)) +
  geom_tile() +
  geom_text(aes(label = Status), size = 2) +
  facet_grid(. ~ Pipeline, scales = "free_x", space = "free_x") +
  labs(fill = "Wilcoxon Rank\nSum Test Statistic", y = "Drug", x = "PSet (Pipeline)") +
  theme_classic() +
  scale_fill_gradient2(low = 'white', mid = '#E1E7DF', high = '#878E76', na.value="#A6A6A6")+
  theme(
    strip.background = element_rect(fill = "#f0f0f0"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.title = element_text(size = 8)
)
dev.off()
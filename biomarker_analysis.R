# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(stats)
})


############################################################
# Specify analysis
############################################################

# analysis = c("circ", "GE")  
# circ for circRNA counts, GE for circ mapped to GE

if (analysis == "circ") {
    path <- "../data/processed_cellline/merged_common_samples/"
    out <- "circ"
}
if (analysis == "GE") {
    path <- "../data/processed_cellline/mergedGE_common_samples/"
    out <- "GE"
}


############################################################
# Load in circRNA expression data
############################################################

# load circRNA expression data
gcsi_df <- fread(paste0(path, "gcsi_counts.tsv"), data.table = F)
ccle_df <- fread(paste0(path, "ccle_counts.tsv"), data.table = F)
gdsc_df <- fread(paste0(path, "gdsc_counts.tsv"), data.table = F)


############################################################
# Load in drug response data and subset
############################################################

# load in drug sensitivity from PSets
load("../data/temp/sensitivity_data.RData")

# drugs of interest
drugs <- intersect(intersect(rownames(gcsi_sen), rownames(ctrp_sen)), rownames(gdsc_sen))
# "Bortezomib", "Crizotinib", "Docetaxel", "Erlotinib", "Pictilisib", "Gemcitabine", "Lapatinib", "Entinostat", "Paclitaxel", "Sirolimus", "Vorinostat" 

# common samples
samples <- rownames(gcsi_df)

# keep common samples and drugs of interest
gcsi_sen <- gcsi_sen[rownames(gcsi_sen) %in% drugs, colnames(gcsi_sen) %in% samples]
ctrp_sen <- ctrp_sen[rownames(ctrp_sen) %in% drugs, colnames(ctrp_sen) %in% samples]
gdsc_sen <- gdsc_sen[rownames(gdsc_sen) %in% drugs, colnames(gdsc_sen) %in% samples]


############################################################
# Filter expression dataframes
############################################################

# remove features with exp in less than threshold samples
thres = 10

gcsi_df <- gcsi_df[,colnames(gcsi_df)[colSums(gcsi_df != 0) >= thres], drop = FALSE]
ccle_df <- ccle_df[,colnames(ccle_df)[colSums(ccle_df != 0) >= thres], drop = FALSE]
gdsc_df <- gdsc_df[,colnames(gdsc_df)[colSums(gdsc_df != 0) >= thres], drop = FALSE]

# zero to NA for downstream analysis
gcsi_df[gcsi_df == 0] <- NA
ccle_df[ccle_df == 0] <- NA
gdsc_df[gdsc_df == 0] <- NA


############################################################
# Get median expression for each transcript
############################################################

# keep only samples with drug response
gcsi_df <- gcsi_df[rownames(gcsi_df) %in% colnames(gcsi_sen),]
ccle_df <- ccle_df[rownames(ccle_df) %in% colnames(ccle_sen),]
gdsc_df <- gdsc_df[rownames(gdsc_df) %in% colnames(gdsc_sen),]

# compute median transcript expression
gcsi_median <- apply(gcsi_df, 2, median, na.rm = TRUE)
ccle_median <- apply(ccle_df, 2, median, na.rm = TRUE)
gdsc_median <- apply(gdsc_df, 2, median, na.rm = TRUE)


############################################################
# Binarize transcript expression by median
############################################################

# function to compute drug response associations
binary_dr <- function(counts_df, median, drug_df) {

    features <- colnames(counts_df)

    # create data frame to hold results
    combinations <- expand.grid(Drug = rownames(drug_df), Feature = features)
    combinations$num_samples <- combinations$pval <-  combinations$W <- NA
    combinations$drug <- combinations$feature <- NA
    
    # initiate row count
    row = 1

    for (i in seq_along(features)) {

        feature <- features[i]

        # binarize transcript expression by median
        subset <- counts_df[,i, drop = FALSE]
        subset <- subset[complete.cases(subset), , drop = FALSE]

        high_exp <- rownames(subset)[as.vector(subset >= median[[feature]])]
        low_exp <- rownames(subset)[-which(rownames(subset) %in% high_exp)]

        # loop through each drug

        for (j in seq_along(rownames(drug_df))) {
            
            high = drug_df[j,colnames(drug_df) %in% high_exp] |> as.numeric()
            low = drug_df[j,colnames(drug_df) %in% low_exp] |> as.numeric()

            # wilcoxon rank sum test
            res <- wilcox.test(high, low, alternative = "two.sided", exact = FALSE)

            # save results
            combinations$W[row] <- res$statistic
            combinations$pval[row] <- res$p.value
            combinations$num_samples[row] <- nrow(subset)
            combinations$feature[row] <- feature
            combinations$drug[row] <- rownames(drug_df)[j]

            row <- row + 1
        }
    }
    combinations$FDR <- p.adjust(combinations$pval, method = "BH")
    return(combinations)
}

gcsi_bin_dr <- binary_dr(gcsi_df, gcsi_median, gcsi_sen)
ccle_bin_dr <- binary_dr(ccle_df, ccle_median, ccle_sen)
gdsc_bin_dr <- binary_dr(gdsc_df, gdsc_median, gdsc_sen)


############################################################
# Save drug response associations
############################################################


############################################################
# Linear model controling for cancer type
############################################################
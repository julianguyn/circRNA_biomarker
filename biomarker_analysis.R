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
    path <- "../data/processed_cellline/merged_all_samples/"
    thres = 15
    dr_out <- "../results/data/bin_dr/circ_"
}
if (analysis == "GE") {                                                 
    path <- "../data/processed_cellline/mergedGE_common_samples/"
    thres = 10       
    dr_out <- "../results/data/bin_dr/GE_"                               
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
dim(gcsi_df)        # 545 30443     # ALL samples:  38858
dim(ccle_df)        # 524 36684                     62271
dim(gdsc_df)        # 141 5322                      11105

# THRESHOLDS FOR circ (intersected samples - in at least 2 PSets):
# thres = 10, gcsi_df: 480, ccle_df: 729, gdsc_df: 236
# thres = 15, gcsi_df: 292, ccle_df: 454, gdsc_df: 153
# thres = 20, gcsi_df: 215, ccle_df: 326, gdsc_df: 112

# THRESHOLDS FOR circ (all samples):
# thres = 10, gcsi_df: 656, ccle_df: , gdsc_df: 655 
# thres = 15, gcsi_df: 392, ccle_df: 1025, gdsc_df: 459
# thres = 20, gcsi_df: , ccle_df: , gdsc_df: 

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


############################################################
# Get median expression for each transcript
############################################################

# keep only samples with drug response
gcsi_df <- gcsi_df[rownames(gcsi_df) %in% colnames(gcsi_sen),]      # 334 392
ccle_df <- ccle_df[rownames(ccle_df) %in% colnames(ccle_sen),]      # 475 1025
gdsc_df <- gdsc_df[rownames(gdsc_df) %in% colnames(gdsc_sen),]      # 314 459

# compute median transcript expression
gcsi_median <- apply(gcsi_df, 2, median, na.rm = TRUE)
ccle_median <- apply(ccle_df, 2, median, na.rm = TRUE)
gdsc_median <- apply(gdsc_df, 2, median, na.rm = TRUE)

save(gcsi_df, ccle_df, gdsc_df, gcsi_median, ccle_median, gdsc_median, file = "../results/data/biomarker_analysis_circ.RData")

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
    print(paste("Number of features:", length(features)))

    for (i in seq_along(features)) {
        
        print(i)
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

            if (length(high[!is.na(high)]) > 1 & length(low[!is.na(low)]) > 1) {

                # wilcoxon rank sum test
                res <- wilcox.test(high, low, alternative = "two.sided", exact = FALSE)

                # save results
                combinations$W[row] <- res$statistic
                combinations$pval[row] <- res$p.value
                combinations$num_samples[row] <- nrow(subset)
                combinations$feature[row] <- feature
                combinations$drug[row] <- rownames(drug_df)[j]

            } else {
                # save results
                combinations$W[row] <- NA
                combinations$pval[row] <- NA
                combinations$num_samples[row] <- nrow(subset)
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

gcsi_bin_dr <- binary_dr(gcsi_df, gcsi_median, gcsi_sen)
ccle_bin_dr <- binary_dr(ccle_df, ccle_median, ccle_sen)
gdsc_bin_dr <- binary_dr(gdsc_df, gdsc_median, gdsc_sen)


############################################################
# Save drug response associations
############################################################

save(gcsi_bin_dr, file = paste0(dr_out, "gcsi.RData"))
save(ccle_bin_dr, file = paste0(dr_out, "ccle.RData"))
save(gdsc_bin_dr, file = paste0(dr_out, "gdsc.RData"))



############################################################
# Linear model controling for cancer type
############################################################


############################################################
# Waterfall plots for individual PSets
############################################################


ggplot(gcsi_bin_dr, aes(x = rank, y = W)) +
    geom_col() + #scale_x_continuous(limits = c(-0.5, 0.5), labels = function(x) x + 0.5) +
    scale_x_discrete(breaks = gcsi_bin_dr$rank, labels = gcsi_bin_dr$Drug) +
    scale_fill_manual(values = pal) +
    theme_classic() + geom_vline(xintercept = 0) + 
    theme(legend.text = element_text(size = 13),
          legend.title = element_text(size = 16),
          axis.text.x = element_text(size = 13),  
          axis.text.y = element_text(size = 13),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16)) +
    labs(y = "Drug", title = "", x = "Concordance Index (CI)", fill = "Signature") 
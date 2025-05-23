# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(stats)
    library(ggplot2)
    library(ComplexHeatmap)
    library(ggh4x)
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
# Quantify (TO REMOVE)
############################################################

load("circ_gdsc.RData")
load("circ_ccle.RData")
load("circ_gcsi.RData")

# number of biomarker associations
# gCSI: 5488
# CCLE: 24600
# GDSC: 33507

# number of biomarker associations with pval < 0.05
# gCSI: 200
# CCLE: 755
# GDSC: 1960

# number of biomarker associations with FDR < 0.05
# gCSI: 0
# CCLE: 0
# GDSC: 9


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

gcsi_pval <- subset_dr(gcsi_bin_dr)
gdsc_pval <- subset_dr(gdsc_bin_dr)
ccle_pval <- subset_dr(ccle_bin_dr)

gcsi_fdr <- subset_dr(gcsi_bin_dr, type = "fdr")
gdsc_fdr <- subset_dr(gdsc_bin_dr, type = "fdr")
ccle_fdr <- subset_dr(ccle_bin_dr, type = "fdr")


############################################################
# Waterfall plot of significant biomarkers
############################################################

plot_waterfall <- function(dr, scale, title) {

    dr <- dr[order(dr$rank),]
    dr$rank <- 1:nrow(dr)

    if (scale == "raw") {
        p <- ggplot(dr, aes(x = rank, y = W)) + geom_col(fill = "#899DA4") + 
        theme_classic() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
        labs(y = "Wilcoxon Rank Sum Test Statistic", x = "Biomarker Association", title = title) 
    }
    if (scale == "log2") {
        p <- ggplot(dr, aes(x = rank, y = log2(W))) + geom_col(fill = "#899DA4") + 
        theme_classic() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
        labs(y = "Wilcoxon Rank Sum Test Statistic", x = "Biomarker Association", title = title) 
    }
    return(p)
}

png("../results/figures/figure6/waterfall_gcsi_pval.png", width = 6, height = 4, res = 600, units = "in")
plot_waterfall(gcsi_pval, scale = "raw", title = "gCSI P-Val Significant Associations")
dev.off()

png("../results/figures/figure6/", width = 6, height = 4, res = 600, units = "in")
plot_waterfall(ccle_pval, scale = "raw", title = "CCLE P-Val Significant Associations")
dev.off()

png("../results/figures/figure6/waterfall_gdsc_pval.png", width = 6, height = 4, res = 600, units = "in")
plot_waterfall(gdsc_pval, scale = "raw", title = "GDSC P-Val Significant Associations")
dev.off()


############################################################
# Upset plot of overlapping biomarkers
############################################################

# function to create upset plot
plot_upset <- function(comb_mat, set_order) {
    p <- UpSet(comb_mat, set_order = set_order,
        top_annotation = upset_top_annotation(comb_mat, add_numbers = TRUE),
        comb_order = order(-comb_size(comb_mat)))
    return(p)
}

# create list object of transcripts for upset plot
toPlot <- make_comb_mat(list(
            gCSI = gcsi_pval$pair,
            CCLE = ccle_pval$pair,
            GDSC = gdsc_pval$pair))

# plot upset plots
png("../results/figures/figure6/bin_dr_upset.png", width = 6, height = 4, res = 600, units = "in")
plot_upset(toPlot, set_order = c("gCSI", "CCLE", "GDSC"))
dev.off()


############################################################
# Wilcoxon results of overlapping biomarkers
############################################################

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

# plot common biomarkers
png("../results/figures/figure6/common_pval_biomarkers.png", width = 17, height = 5, res = 600, units = "in")
ggplot(toPlot, aes(x = PSet, y = W, fill = pval)) + geom_bar(stat="identity", color = "black") +
    facet_nested(~ factor(pair), scales = "free_x") +
    labs(fill = "P-Value", y = "Wilcoxon Rank Sum Test Statistic", x = "PSet") + 
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5), axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()


############################################################
# Linear model controling for cancer type
############################################################

runLM <- function(exp, sen) {

    features = colnames(exp)

    # create data frame to hold results
    combinations <- expand.grid(Drug = rownames(sen), Feature = features)
    combinations$num_samples <- combinations$r2 <- combinations$se <- combinations$pval <-  combinations$estimate <- NA
    combinations$drug <- combinations$feature <- NA

    # initiate row count
    row = 1
    print(paste("Number of features:", length(features)))

    # loop through each feature
    for (i in seq_along(features)) {
        
        print(i)
        # extract feature vector
        feature <- features[i]
        exp_vector <- exp[,i,drop=F] 
        exp_vector <- exp_vector[complete.cases(exp_vector),,drop=F]

        # loop through each drug
        for (j in seq_along(rownames(sen))) {

            # get common cell lines
            common_ccls <- intersect(rownames(exp_vector), colnames(sen))

            # extract drug sen vector and match ccls
            drug_vector <- sen[j,colnames(sen) %in% common_ccls] |> t() |> as.data.frame()
            drug_vector <- drug_vector[complete.cases(drug_vector),,drop=F]

            if (nrow(drug_vector) < 2) {

                combinations$estimate[row] <- NA
                combinations$pval[row] <- NA
                combinations$se[row] <- NA
                combinations$r2[row] <- NA
                combinations$num_samples[row] <- 0
                combinations$feature[row] <- feature
                combinations$drug[row] <- colnames(drug_vector)

            } else {
                exp_vector_subset <- exp_vector[match(rownames(drug_vector), rownames(exp_vector)),,drop=F]
            
                # run linear regression
                circ_exp <- exp_vector_subset[,1]
                drug_vec <- drug_vector[,1]
                
                res <- summary(lm(drug_vec ~ circ_exp))

                # save results
                combinations$estimate[row] <- as.numeric(res$coefficients["circ_exp","Estimate"])
                combinations$pval[row] <- as.numeric(res$coefficients["circ_exp","Pr(>|t|)"])
                combinations$se[row] <- as.numeric(res$coefficients["circ_exp","Std. Error"])
                combinations$r2[row] <- res$'r.squared'
                combinations$num_samples[row] <- length(common_ccls)
                combinations$feature[row] <- feature
                combinations$drug[row] <- colnames(drug_vector)
            }

            row <- row + 1
        }
    }
    combinations$FDR <- p.adjust(combinations$pval, method = "BH")

    # format dataframe for plotting 
    combinations <- combinations[order(combinations$estimate, decreasing = T),]
    combinations$rank <- seq_len(nrow(combinations))
    combinations$pair <- paste(combinations$Drug, combinations$Feature, sep = "_")
    return(combinations)
}

gcsi_lm_dr <- runLM(gcsi_df, gcsi_sen)
ccle_lm_dr <- runLM(ccle_df, ccle_sen)
gdsc_lm_dr <- runLM(gdsc_df, gdsc_sen)


############################################################
# Save drug response associations
############################################################

save(gcsi_lm_dr, file = paste0(lm_out, "gcsi_lm.RData"))
save(ccle_lm_dr, file = paste0(lm_out, "ccle_lm.RData"))
save(gdsc_lm_dr, file = paste0(lm_out, "gdsc_lm.RData"))

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
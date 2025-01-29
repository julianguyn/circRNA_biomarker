# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggpubr)
    library(reshape2)
    library(matrixStats)
})

options(stringsAsFactors = FALSE)
set.seed(123)

############################################################
# Load in data 
############################################################

# load circRNA expression data
path <- "../data/processed_cellline/GE_common_samples/"

ciri_gcsi <- fread(paste0(path, "CIRI2/ciri_gcsi_counts.tsv"), data.table = F)                          #ncol: 2502
ciri_gdsc <- fread(paste0(path, "CIRI2/ciri_gdsc_counts.tsv"), data.table = F)                          #ncol: 3039
ciri_ccle <- fread(paste0(path, "CIRI2/ciri_ccle_counts.tsv"), data.table = F)                          #ncol: 2076

circ_gcsi <- fread(paste0(path, "CIRCexplorer2/circ_gcsi_counts.tsv"), data.table = F)                  #ncol: 3658
circ_gdsc <- fread(paste0(path, "CIRCexplorer2/circ_gdsc_counts.tsv"), data.table = F)                  #ncol: 3004
circ_ccle <- fread(paste0(path, "CIRCexplorer2/circ_ccle_counts.tsv"), data.table = F)                  #ncol: 2680

cfnd_gcsi <- fread(paste0(path, "circRNA_finder/cfnd_gcsi_counts.tsv"), data.table = F)                 #ncol: 9459
cfnd_gdsc <- fread(paste0(path, "circRNA_finder/cfnd_gdsc_counts.tsv"), data.table = F)                 #ncol: 9684
cfnd_ccle <- fread(paste0(path, "circRNA_finder/cfnd_ccle_counts.tsv"), data.table = F)                 #ncol: 3596

fcrc_gcsi <- fread(paste0(path, "find_circ/fcrc_gcsi_counts.tsv"), data.table = F)                      #ncol: 29898
fcrc_gdsc <- fread(paste0(path, "find_circ/fcrc_gdsc_counts.tsv"), data.table = F)                      #ncol: 30660
fcrc_ccle <- fread(paste0(path, "find_circ/fcrc_ccle_counts.tsv"), data.table = F)                      #ncol: 30728   


############################################################
# Filter low exp transcripts
############################################################

ciri_gcsi <- ciri_gcsi[,-which(colnames(ciri_gcsi) %in% names(which(colSums(ciri_gcsi == 0) > 45)))]    #ncol: 409
ciri_ccle <- ciri_ccle[,-which(colnames(ciri_ccle) %in% names(which(colSums(ciri_ccle == 0) > 45)))]    #ncol: 480
ciri_gdsc <- ciri_gdsc[,-which(colnames(ciri_gdsc) %in% names(which(colSums(ciri_gdsc == 0) > 45)))]    #ncol: 300

circ_gcsi <- circ_gcsi[,-which(colnames(circ_gcsi) %in% names(which(colSums(circ_gcsi == 0) > 45)))]    #ncol: 679
circ_ccle <- circ_ccle[,-which(colnames(circ_ccle) %in% names(which(colSums(circ_ccle == 0) > 45)))]    #ncol: 480
circ_gdsc <- circ_gdsc[,-which(colnames(circ_gdsc) %in% names(which(colSums(circ_gdsc == 0) > 45)))]    #ncol: 460

cfnd_gcsi <- cfnd_gcsi[,-which(colnames(cfnd_gcsi) %in% names(which(colSums(cfnd_gcsi == 0) > 45)))]    #ncol: 3798
cfnd_ccle <- cfnd_ccle[,-which(colnames(cfnd_ccle) %in% names(which(colSums(cfnd_ccle == 0) > 45)))]    #ncol: 3669
cfnd_gdsc <- cfnd_gdsc[,-which(colnames(cfnd_gdsc) %in% names(which(colSums(cfnd_gdsc == 0) > 45)))]    #ncol: 586

fcrc_gcsi <- fcrc_gcsi[,-which(colnames(fcrc_gcsi) %in% names(which(colSums(fcrc_gcsi == 0) > 45)))]    #ncol: 28390
fcrc_ccle <- fcrc_ccle[,-which(colnames(fcrc_ccle) %in% names(which(colSums(fcrc_ccle == 0) > 45)))]    #ncol: 29738
fcrc_gdsc <- fcrc_gdsc[,-which(colnames(fcrc_gdsc) %in% names(which(colSums(fcrc_gdsc == 0) > 45)))]    #ncol: 30029


############################################################
# Keep transcripts common across all psets
############################################################

# function to order circRNA dataframes and keep only transcripts found in all psets for each pipeline
subset_df <- function(df, common_transcripts) {

    # keep circRNA transcripts found in all psets for each pipeline
    df <- df[,which(colnames(df) %in% common_transcripts)]
    df <- df[order(df$sample),order(colnames(df))]

    # remove sample names
    rownames(df) <- df$sample
    df$sample <- NULL

    return(df)
}

# get common gene transcripts found in all psets for each pipeline
ciri_common <- intersect(intersect(colnames(ciri_gcsi), colnames(ciri_ccle)), colnames(ciri_gdsc))      #n: 124
circ_common <- intersect(intersect(colnames(circ_gcsi), colnames(circ_ccle)), colnames(circ_gdsc))      #n: 124
cfnd_common <- intersect(intersect(colnames(cfnd_gcsi), colnames(cfnd_ccle)), colnames(cfnd_gdsc))      #n: 311
fcrc_common <- intersect(intersect(colnames(fcrc_gcsi), colnames(fcrc_ccle)), colnames(fcrc_gdsc))      #n: 28185

# filter for common transcripts
ciri_gcsi <- subset_df(ciri_gcsi, ciri_common)
ciri_ccle <- subset_df(ciri_ccle, ciri_common)
ciri_gdsc <- subset_df(ciri_gdsc, ciri_common)

circ_gcsi <- subset_df(circ_gcsi, circ_common)
circ_ccle <- subset_df(circ_ccle, circ_common)
circ_gdsc <- subset_df(circ_gdsc, circ_common)

cfnd_gcsi <- subset_df(cfnd_gcsi, cfnd_common)
cfnd_ccle <- subset_df(cfnd_ccle, cfnd_common)
cfnd_gdsc <- subset_df(cfnd_gdsc, cfnd_common)

fcrc_gcsi <- subset_df(fcrc_gcsi, fcrc_common)
fcrc_ccle <- subset_df(fcrc_ccle, fcrc_common)
fcrc_gdsc <- subset_df(fcrc_gdsc, fcrc_common)


############################################################
# Compute pairwise Spearman corr from datasets 
############################################################

# function to compute pairwise spearman correlations
compute_spearman <- function(gcsi_df, ccle_df, gdsc_df, random = FALSE, iter = 1) {

    # INPUTS:
    #   gcsi_df, ccle_df, gdsc_df: PSet-specific dataframes from the subset_df() function
    #   random: TRUE for random sampling of sample names, FALSE otherwise
    #   iter: number of iterations to be performed (for random sampling and SI computation)

    # initialize dataframe to store results
    correlations <- data.frame(matrix(nrow=0, ncol=3))

    # loop through for number of iterations
    for (i in 1:iter) {
    
        if (random == TRUE) {
            # shuffle cell line names
            rownames(gcsi_df) <- sample(rownames(gcsi_df))
            rownames(ccle_df) <- sample(rownames(ccle_df))
            rownames(gdsc_df) <- sample(rownames(gdsc_df))
            # order dataframe
            gcsi_df <- gcsi_df[order(rownames(gcsi_df)),]
            ccle_df <- ccle_df[order(rownames(ccle_df)),]
            gdsc_df <- gdsc_df[order(rownames(gdsc_df)),]
        }
        
        # loop through each common transcript
        for (i in 1:ncol(gcsi_df)) {
            # compute correlations of transcript expression for pairs of psets
            gcsi_ccle_spearman <- suppressWarnings(cor(x = as.numeric(gcsi_df[, i]), y = as.numeric(ccle_df[, i]), method = "spearman")) #gCSI vs CCLE
            gcsi_gdsc_spearman <- suppressWarnings(cor(x = as.numeric(gcsi_df[, i]), y = as.numeric(gdsc_df[, i]), method = "spearman")) #gCSI vs GDSC
            gdsc_ccle_spearman <- suppressWarnings(cor(x = as.numeric(gdsc_df[, i]), y = as.numeric(ccle_df[, i]), method = "spearman")) #GDSC vs CCLE
            # combine results
            correlations <- rbind(correlations, c(gcsi_ccle_spearman, gcsi_gdsc_spearman, gdsc_ccle_spearman))
        }
    } 
    colnames(correlations) <- c("gCSI/CCLE", "gCSI/GDSC2", "GDSC2/CCLE")
    return(correlations)
}

# compute spearman correlations
ciri_stability <- compute_spearman(ciri_gcsi, ciri_ccle, ciri_gdsc)
circ_stability <- compute_spearman(circ_gcsi, circ_ccle, circ_gdsc)
cfnd_stability <- compute_spearman(cfnd_gcsi, cfnd_ccle, cfnd_gdsc)
fcrc_stability <- compute_spearman(fcrc_gcsi, fcrc_ccle, fcrc_gdsc)

# compute spearman correlations after random shuffling
ciri_stability_random <- compute_spearman(ciri_gcsi, ciri_ccle, ciri_gdsc, random = TRUE, iter = 100)
circ_stability_random <- compute_spearman(circ_gcsi, circ_ccle, circ_gdsc, random = TRUE, iter = 100)
cfnd_stability_random <- compute_spearman(cfnd_gcsi, cfnd_ccle, cfnd_gdsc, random = TRUE, iter = 100)
fcrc_stability_random <- compute_spearman(fcrc_gcsi, fcrc_ccle, fcrc_gdsc, random = TRUE, iter = 100)

save(ciri_stability, circ_stability, cfnd_stability, fcrc_stability, 
     file = "../results/data/temp/circ_ge_stability.RData")
save(ciri_stability_random, circ_stability_random, cfnd_stability_random, fcrc_stability_random, 
     file = "../results/data/temp/circ_ge_stability_random.RData")


############################################################
# Format stability index matrices for plotting
############################################################

# function to format stability dataframes
format_df <- function(df, label, random = "NonRandom") {
    colnames(df) <- c("gCSI/CCLE", "gCSI/GDSC2", "GDSC2/CCLE")
    toPlot <- melt(df)
    colnames(toPlot) <- c("PSet", "Stability")
    toPlot$label <- label
    toPlot$random <- random
    return(toPlot)
}

ciri_stability <- format_df(ciri_stability, "CIRI2", "NonRandom")
circ_stability <- format_df(circ_stability, "CIRCexplorer2", "NonRandom")
cfnd_stability <- format_df(cfnd_stability, "circRNA_finder", "NonRandom")
fcrc_stability <- format_df(fcrc_stability, "find_circ", "NonRandom")

ciri_stability_random <- format_df(ciri_stability_random, "CIRI2", "Random")
circ_stability_random <- format_df(circ_stability_random, "CIRCexplorer2", "Random")
cfnd_stability_random <- format_df(cfnd_stability_random, "circRNA_finder", "Random")
fcrc_stability_random <- format_df(fcrc_stability_random, "find_circ", "Random")

# load in gene expression and isoform stability
load("../results/data/temp/gene_isoform_stability.RData")

transcript_stability <- format_df(transcript_stability[,c("gcsi_ccle_spearman", "gcsi_gdsc_spearman", "gdsc_ccle_spearman")], "Isoforms")
gene_stability <- format_df(gene_stability[,c("gcsi_ccle_spearman", "gcsi_gdsc_spearman", "gdsc_ccle_spearman")], "Gene Expression")
gene_stability[is.na(gene_stability)] <- 0

############################################################
# Figure 5: Plot stability index distribution 
############################################################

# merge nonrandom results for plotting
toPlot <- rbind(gene_stability, transcript_stability, 
                ciri_stability, circ_stability, 
                cfnd_stability, fcrc_stability)
toPlot$label <- factor(toPlot$label, levels = c("Gene Expression", "Isoforms", "CIRI2", "CIRCexplorer2", "circRNA_finder", "find_circ"))

png("../results/figures/figure5/stability_nonrandom_ge.png", width=300, height=150, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = label, y = Stability)) + 
    geom_violin(aes(fill = label), alpha = 0.8) + geom_boxplot(width=0.1, alpha = 0.3) +
    facet_grid(factor(PSet)~.) +
    theme_classic() + labs(x = "", fill = "", y = "Stability Index") +
    scale_fill_manual(values = c("#23022E", "#611C35", "#839788", "#BFD7EA", "#BA9790", "#D5BC8A")) +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3), legend.key.size = unit(0.7, 'cm')) +
    geom_hline(yintercept = 0, linetype = "dotted")
dev.off()


# merge random results for plotting
toPlot <- rbind(ciri_stability, circ_stability, cfnd_stability, fcrc_stability,
                ciri_stability_random, circ_stability_random, cfnd_stability_random, fcrc_stability_random)
toPlot <- melt(toPlot)
toPlot$label <- factor(toPlot$label, levels = c("CIRI2", "CIRCexplorer2", "circRNA_finder", "find_circ"))

png("../results/figures/figure5/stability_random_ge.png", width=150, height=150, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = PSet, y = value, fill = random)) + 
    geom_boxplot() + facet_grid(label~.) + theme_classic() + 
    labs(x = "", fill = "", y = "Stability Index") + scale_fill_manual(values = c("#839788", "gray")) +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3), legend.key.size = unit(0.7, 'cm')) +
    geom_hline(yintercept = 0, linetype = "dotted")
dev.off()

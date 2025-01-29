# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggpubr)
    library(ggvenn)
    library(reshape2)
    library(matrixStats)
})

options(stringsAsFactors = FALSE)
set.seed(123)

# TODO: format counts dataframes and run randomization

############################################################
# Load in data 
############################################################

# load circRNA expression data
#path <- "../data/processed_cellline/GE_common_samples/"

#ciri_polyA <- fread(paste0(path, "CIRI2/ciri_gcsi_counts.tsv"), data.table = F)                          #ncol: 2502
#ciri_ribo0 <- fread(paste0(path, "CIRI2/ciri_gdsc_counts.tsv"), data.table = F)                          #ncol: 3039

#circ_polyA <- fread(paste0(path, "CIRCexplorer2/circ_gcsi_counts.tsv"), data.table = F)                  #ncol: 3658
#circ_ribo0 <- fread(paste0(path, "CIRCexplorer2/circ_gdsc_counts.tsv"), data.table = F)                  #ncol: 3004

#cfnd_polyA <- fread(paste0(path, "circRNA_finder/cfnd_gcsi_counts.tsv"), data.table = F)                 #ncol: 9459
#cfnd_ribo0 <- fread(paste0(path, "circRNA_finder/cfnd_gdsc_counts.tsv"), data.table = F)                 #ncol: 9684

#fcrc_polyA <- fread(paste0(path, "find_circ/fcrc_gcsi_counts.tsv"), data.table = F)                      #ncol: 29898
#fcrc_ribo0 <- fread(paste0(path, "find_circ/fcrc_gdsc_counts.tsv"), data.table = F)                      #ncol: 30660


load("../data/processed_lung/circ_lung_expression.RData")      # from circ_lung.R

####################################################
#  circRNA Quantification Comparisons 
####################################################


# set palette for plotting
pal = c("#839788", "#BFD7EA", "#BA9790", "#D5BC8A")


# create list object of transcripts
polyA_transcripts <- list(CIRI2 = colnames(ciri_polyA), CIRCexplorer2 = colnames(circ_polyA), circRNA_finder = colnames(cfnd_polyA), find_circ = colnames(fcrc_polyA))
ribo0_transcripts <- list(CIRI2 = colnames(ciri_ribo0), CIRCexplorer2 = colnames(circ_ribo0), circRNA_finder = colnames(cfnd_ribo0), find_circ = colnames(fcrc_ribo0))

# plot venn diagram
p1 <- ggvenn(polyA_transcripts, 
        fill_color = pal, stroke_size = 0.5, set_name_size = 4) + 
        theme(plot.title = element_text(hjust = 0.5, size = 15)) + labs(title = "poly(A)-Selection\n")
p2 <- ggvenn(ribo0_transcripts, 
        fill_color = pal, stroke_size = 0.5, set_name_size = 4) + 
        theme(plot.title = element_text(hjust = 0.5, size = 15)) + labs(title = "RiboZero\n")

png("../results/figures/figure7/venndiagram_per_protocol.png", width=250, height=150, units='mm', res = 600, pointsize=80)
ggarrange(p1, p2, ncol = 2, nrow = 1, common.legend = FALSE)
dev.off()


############################################################
# Filter low exp transcripts
############################################################

ciri_polyA <- ciri_polyA[,-which(colnames(ciri_polyA) %in% names(which(colSums(ciri_polyA == 0) > 48)))]    #ncol: 405
ciri_ribo0 <- ciri_ribo0[,-which(colnames(ciri_ribo0) %in% names(which(colSums(ciri_ribo0 == 0) > 48)))]    #ncol: 1103

circ_polyA <- circ_polyA[,-which(colnames(circ_polyA) %in% names(which(colSums(circ_polyA == 0) > 48)))]    #ncol: 7
circ_ribo0 <- circ_ribo0[,-which(colnames(circ_ribo0) %in% names(which(colSums(circ_ribo0 == 0) > 48)))]    #ncol: 750

cfnd_polyA <- cfnd_polyA[,-which(colnames(cfnd_polyA) %in% names(which(colSums(cfnd_polyA == 0) > 48)))]    #ncol: 357
cfnd_ribo0 <- cfnd_ribo0[,-which(colnames(cfnd_ribo0) %in% names(which(colSums(cfnd_ribo0 == 0) > 48)))]    #ncol: 853

fcrc_polyA <- fcrc_polyA[,-which(colnames(fcrc_polyA) %in% names(which(colSums(fcrc_polyA == 0) > 48)))]    #ncol: 742
fcrc_ribo0 <- fcrc_ribo0[,-which(colnames(fcrc_ribo0) %in% names(which(colSums(fcrc_ribo0 == 0) > 48)))]    #ncol: 476

# REMOVED CIRC DUE TO DISCREPANCY

############################################################
# Keep transcripts common across all pipelines
############################################################

# function to order circRNA dataframes and keep only transcripts found in both samples for each pipeline
subset_df <- function(df, common_transcripts) {

    # keep circRNA transcripts found in all psets for each pipeline
    df <- df[,which(colnames(df) %in% common_transcripts)]
    df <- df[,order(colnames(df))]

    return(df)
}

# get common gene transcripts found in all psets for each pipeline
polyA_common <- intersect(intersect(colnames(ciri_polyA), colnames(cfnd_polyA)), colnames(fcrc_polyA))     #n: 13
ribo0_common <- intersect(intersect(colnames(ciri_ribo0), colnames(cfnd_ribo0)), colnames(fcrc_ribo0))     #n: 687

# filter for common transcripts
ciri_polyA <- subset_df(ciri_polyA, polyA_common)
cfnd_polyA <- subset_df(cfnd_polyA, polyA_common)
fcrc_polyA <- subset_df(fcrc_polyA, polyA_common)

ciri_ribo0 <- subset_df(ciri_ribo0, ribo0_common)
cfnd_ribo0 <- subset_df(cfnd_ribo0, ribo0_common)
fcrc_ribo0 <- subset_df(fcrc_ribo0, ribo0_common)


############################################################
# Compute pairwise Spearman corr from datasets 
############################################################

# function to compute pairwise spearman correlations
compute_spearman <- function(
    ciri_df, cfnd_df, fcrc_df, 
    random = FALSE,         #   random: TRUE for random sampling of sample names, FALSE otherwise
    iter = 1                #   iter: number of iterations to be performed (for random sampling and SI computation)
) {

    # initialize dataframe to store results
    correlations <- data.frame(matrix(nrow=0, ncol=3))

    # loop through for number of iterations
    for (i in 1:iter) {
    
        if (random == TRUE) {
            # shuffle cell line names
            rownames(ciri_df) <- sample(rownames(ciri_df))
            rownames(cfnd_df) <- sample(rownames(cfnd_df))
            rownames(fcrc_df) <- sample(rownames(fcrc_df))
        }

        # order dataframe
        ciri_df <- ciri_df[order(rownames(ciri_df)),]
        cfnd_df <- cfnd_df[order(rownames(cfnd_df)),]
        fcrc_df <- fcrc_df[order(rownames(fcrc_df)),]
        
        # loop through each common transcript
        for (i in 1:ncol(ciri_df)) {
            ciri_cfnd <- suppressWarnings(cor(x = as.numeric(ciri_df[, i]), y = as.numeric(cfnd_df[, i]), method = "spearman")) 
            ciri_fcrc <- suppressWarnings(cor(x = as.numeric(ciri_df[, i]), y = as.numeric(fcrc_df[, i]), method = "spearman")) 
            cfnd_fcrc <- suppressWarnings(cor(x = as.numeric(cfnd_df[, i]), y = as.numeric(fcrc_df[, i]), method = "spearman")) 

            correlations <- rbind(correlations, c(ciri_cfnd, ciri_fcrc, cfnd_fcrc))
            rownames(correlations)[i] <- colnames(ciri_df)[i]
        }
    } 
    colnames(correlations) <- c("CIRI/CFND", "CIRI/FCRC", "CFND/FCRC")
    return(correlations)
}

# compute spearman correlations
polyA_stability <- compute_spearman(ciri_polyA, cfnd_polyA, fcrc_polyA)
ribo0_stability <- compute_spearman(ciri_ribo0, cfnd_ribo0, fcrc_ribo0)

save(polyA_stability, ribo0_stability, file = "../results/data/temp/circ_lung_stability.RData")

# compute spearman correlations after random shuffling
polyA_stability_random <- compute_spearman(ciri_polyA, circ_polyA, cfnd_polyA, fcrc_polyA, random = TRUE, iter = 100)
ribo0_stability_random <- compute_spearman(ciri_ribo0, circ_ribo0, cfnd_ribo0, fcrc_ribo0, random = TRUE, iter = 100)

save(polyA_stability_random, ribo0_stability_random, file = "../results/data/temp/circ_lung_stability_random.RData")


############################################################
# Format stability index matrices for plotting
############################################################

# function to format stability dataframes
format_df <- function(df, label, random = "NonRandom") {
    toPlot <- reshape2::melt(df)
    colnames(toPlot) <- c("Pipeline", "Stability")
    toPlot$label <- label
    toPlot$random <- random
    return(toPlot)
}

polyA_stability <- format_df(polyA_stability, "polyA", "NonRandom")
ribo0_stability <- format_df(ribo0_stability, "ribo0", "NonRandom")

polyA_stability_random <- format_df(polyA_stability_random, "polyA", "Random")
ribo0_stability_random <- format_df(ribo0_stability_random, "ribo0", "Random")


############################################################
# Plot stability index distribution 
############################################################

# merge nonrandom results for plotting
toPlot <- rbind(polyA_stability, ribo0_stability)
toPlot$label <- factor(toPlot$label, levels = c("polyA", "ribo0"))

png("../results/figures/figure7/stability_nonrandom_lung.png", width=100, height=150, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = label, y = Stability)) + 
    geom_violin(aes(fill = label), alpha = 0.8) + geom_boxplot(width=0.1, alpha = 0.3) +
    facet_grid(factor(Pipeline)~.) +
    theme_classic() + labs(x = "", fill = "", y = "Stability Index") +
    scale_fill_manual(values = c("#23022E", "#611C35", "#839788", "#BFD7EA", "#BA9790", "#D5BC8A")) +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3), legend.key.size = unit(0.7, 'cm')) +
    geom_hline(yintercept = 0, linetype = "dotted")
dev.off()


# merge random results for plotting
toPlot <- rbind(polyA_stability, ribo0_stability, polyA_stability_random, ribo0_stability_random)
#toPlot <- melt(toPlot)
toPlot$label <- factor(toPlot$label, levels = c("polyA", "ribo0"))

png("../results/figures/figure7/stability_random_lung.png", width=150, height=150, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = PSet, y = value, fill = random)) + 
    geom_boxplot() + facet_grid(label~.) + theme_classic() + 
    labs(x = "", fill = "", y = "Stability Index") + scale_fill_manual(values = c("#839788", "gray")) +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3), legend.key.size = unit(0.7, 'cm')) +
    geom_hline(yintercept = 0, linetype = "dotted")
dev.off()

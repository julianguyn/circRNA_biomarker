# Script to correlate circRNA expression across biological replicates (alternative to UMAP)

# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(reshape2)
    library(ggplot2)
})


############################################################
# Load in processed expression data
############################################################

# load processed expression data
load("../results/data/umapdf.RData")


############################################################
# Pairwise correlations of biological reps
############################################################

# function to parse processed dataframe and perform pairwise correlations
corr_reps <- function(df, cell_line_labels) {
    
    ## current df structure: 144 rows (48 cell lines x 3 PSets) x unique cell lines (columns)

    # split dataframe into three psets
    p1 <- df[1:48,]
    p2 <- df[49:96,]
    p3 <- df[97:144,]

    # add cell line labels
    rownames(p1) <- cell_line_labels[1:48]
    rownames(p2) <- cell_line_labels[49:96]
    rownames(p3) <- cell_line_labels[97:144]

    # order cell lines
    p1 <- p1[order(rownames(p1)),]
    p2 <- p2[order(rownames(p2)),]
    p3 <- p3[order(rownames(p3)),]

    # initiate dataframe for storing results
    correlations <- data.frame(matrix(nrow=0, ncol=4))
    colnames(correlations) <- c("Cells", "gCSI_CCLE", "gCSI_GDSC", "GDSC_CCLE")

    # loop through each common transcript
    for (i in 1:nrow(p1)) {

        cell <- rownames(p1)[i]

        # compute correlations of transcript expression for pairs of psets
        p1_p2 <- suppressWarnings(cor(x = as.numeric(p1[i,]), y = as.numeric(p2[i,]), method = "spearman")) #gCSI vs CCLE
        p1_p3 <- suppressWarnings(cor(x = as.numeric(p1[i,]), y = as.numeric(p3[i,]), method = "spearman")) #gCSI vs GDSC
        p2_p3 <- suppressWarnings(cor(x = as.numeric(p2[i,]), y = as.numeric(p3[i,]), method = "spearman")) #GDSC vs CCLE
    
        # combine results
        correlations <- rbind(correlations, 
                            data.frame(Cells = cell, gCSI_CCLE = p1_p2, gCSI_GDSC = p1_p3, GDSC_CCLE = p2_p3))
    }

    return(correlations)
}

gexpr_corr <- corr_reps(gexpr_df, cell_line_gexpr)
isoform_corr <- corr_reps(isoform_df, cell_line_isoforms)
ciri_corr <- corr_reps(ciri_df, cell_line_ciri)
circ_corr <- corr_reps(circ_df, cell_line_circ)
cfnd_corr <- corr_reps(cfnd_df, cell_line_cfnd)
fcrc_corr <- corr_reps(fcrc_df, cell_line_fcrc)

save(ciri_corr, circ_corr, cfnd_corr, fcrc_corr, gexpr_corr, isoform_corr, 
     file = "../results/data/corr_expr.RData")


############################################################
# Format corrleation matrices for plotting
############################################################

# function to format correlation dataframes
format_df <- function(df, label) {
  colnames(df) <- c("Cells", "gCSI/CCLE", "gCSI/GDSC2", "GDSC2/CCLE")
  toPlot <- melt(df)
  colnames(toPlot) <- c("Cells", "PSet", "Correlation")
  toPlot$label <- label
  return(toPlot)
}

gexpr_corr <- format_df(gexpr_corr, "Gene Expression")
isoform_corr <- format_df(isoform_corr, "Isoforms")
ciri_corr <- format_df(ciri_corr, "CIRI2")
circ_corr <- format_df(circ_corr, "CIRCexplorer2")
cfnd_corr <- format_df(cfnd_corr, "circRNA_finder")
fcrc_corr <- format_df(fcrc_corr, "find_circ")

# merge results for plotting
toPlot <- rbind(gexpr_corr, isoform_corr, 
                ciri_corr, circ_corr, 
                cfnd_corr, fcrc_corr)
toPlot[is.na(toPlot)] <- 0
toPlot$label <- factor(toPlot$label, levels = c("Gene Expression", "Isoforms", "CIRI2", "CIRCexplorer2", "circRNA_finder", "find_circ"))


############################################################
# Plot correlation of biological reps
############################################################

# violin plots
png("../results/figures/figure2/correlate_reps.png", width=200, height=75, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = Correlation, y = label)) + 
    geom_violin(aes(fill = label, alpha = PSet), scale = "width", width = 1.1, position = position_dodge(width = 0.8)) + 
    theme_classic() + labs(x = "Spearman Correlation", fill = "", alpha = "", y = "") +
    scale_fill_manual(values = c("#23022E", "#611C35", "#839788", "#BFD7EA", "#BA9790", "#D5BC8A")) +
    scale_alpha_manual(values = c(0.2, 0.5, 0.9)) +
    scale_y_discrete(limits = c("find_circ", "circRNA_finder", "CIRCexplorer2", "CIRI2", "Isoforms", "Gene Expression")) + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3), 
          legend.key.size = unit(0.5, 'cm')) 
dev.off()

# density plots
png("../results/figures/figure2/corr_dist_density.png", width=150, height=50, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = Correlation)) + geom_density(aes(fill = label), alpha = 0.4, size = 0.5) + 
        theme_classic() +  
        scale_fill_manual(values = c("#23022E", "#611C35", "#839788", "#BFD7EA", "#BA9790", "#D5BC8A")) +
        labs(x = "Spearman Correlation", y = "Density", fill = "") +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
            legend.key.size = unit(0.4, 'cm'))
dev.off()
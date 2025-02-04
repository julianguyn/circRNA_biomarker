# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(plyr)
    library(dplyr)
    library(reshape2)
    library(tidyverse)
    library(umap)
    library(ggplot2)
    library(ggpubr)
})

options(stringsAsFactors = FALSE)
set.seed(101)

############################################################
# Load in data 
############################################################

load("../data/processed_lung/circ_lung_expression.RData")      # from circ_lung.R


############################################################
# Filter matrices by transcripts
############################################################

# filter transcripts with low expression across samples
ciri_polyA <- ciri_polyA[,-which(colnames(ciri_polyA) %in% names(which(colSums(ciri_polyA == 0) > 48)))]    #ncol: 405
ciri_ribo0 <- ciri_ribo0[,-which(colnames(ciri_ribo0) %in% names(which(colSums(ciri_ribo0 == 0) > 48)))]    #ncol: 1103

circ_polyA <- circ_polyA[,-which(colnames(circ_polyA) %in% names(which(colSums(circ_polyA == 0) > 48)))]    #ncol: 7
circ_ribo0 <- circ_ribo0[,-which(colnames(circ_ribo0) %in% names(which(colSums(circ_ribo0 == 0) > 48)))]    #ncol: 750

cfnd_polyA <- cfnd_polyA[,-which(colnames(cfnd_polyA) %in% names(which(colSums(cfnd_polyA == 0) > 48)))]    #ncol: 357
cfnd_ribo0 <- cfnd_ribo0[,-which(colnames(cfnd_ribo0) %in% names(which(colSums(cfnd_ribo0 == 0) > 48)))]    #ncol: 853

fcrc_polyA <- fcrc_polyA[,-which(colnames(fcrc_polyA) %in% names(which(colSums(fcrc_polyA == 0) > 48)))]    #ncol: 742
fcrc_ribo0 <- fcrc_ribo0[,-which(colnames(fcrc_ribo0) %in% names(which(colSums(fcrc_ribo0 == 0) > 48)))]    #ncol: 476

# get common transcripts
transcripts <- data.frame(c(colnames(ciri_polyA), colnames(ciri_ribo0),
                            colnames(circ_polyA), colnames(circ_ribo0),
                            colnames(cfnd_polyA), colnames(cfnd_ribo0),
                            colnames(fcrc_polyA), colnames(fcrc_ribo0)))
colnames(transcripts) <- "transcriptID"
transcript_counts <- transcripts %>% count(transcriptID)

# remove transcripts only in one pipeline_protocol object (keep any in 2 or more)
transcript_counts <- transcript_counts[-which(transcript_counts$n == 1),] 
common_transcripts <- transcript_counts$transcriptID                                                    #n: 2256


############################################################
# Create merged dataframes 
############################################################

# REMOVE CIRCexplorer2

# function to create merged dataframe for each pipeline
mergePSet <- function(ciri_df, cfnd_df, fcrc_df) {
  # keep only common transcripts
  ciri_df <- ciri_df[,which(colnames(ciri_df) %in% common_transcripts)]
  cfnd_df <- cfnd_df[,which(colnames(cfnd_df) %in% common_transcripts)]
  fcrc_df <- fcrc_df[,which(colnames(fcrc_df) %in% common_transcripts)]

  # merge into one dataframe
  df <- rbind.fill(ciri_df, cfnd_df, fcrc_df)
  df[is.na(df)] <- 0
  df[] <- lapply(df, as.double)

  return(df)
}

polyA_df <- mergePSet(ciri_polyA, cfnd_polyA, fcrc_polyA)
ribo0_df <- mergePSet(ciri_ribo0, cfnd_ribo0, fcrc_ribo0)


############################################################
# Spearman corr of transcript exp across biological reps
############################################################

# function to parse processed dataframe and perform pairwise correlations
corr_reps <- function(exp_df) {

    # parse datasets        # already ordered by tumour label
    ciri_df <- exp_df[1:51,]
    cfnd_df <- exp_df[52:102,]
    fcrc_df <- exp_df[103:153,]
    
    # initiate dataframe for storing results
    correlations <- data.frame(matrix(nrow=1, ncol=3))
    colnames(correlations) <- c("CIRI_CNFD", "CIRI_FCRC", "CFND_FCRC")

    # loop through each common transcript
    for (i in 1:nrow(ciri_df)) {

        # compute correlations of transcript expression for pairs of psets
        ciri_cfnd <- suppressWarnings(cor(x = as.numeric(ciri_df[i,]), y = as.numeric(cfnd_df[i,]), method = "spearman")) #gCSI vs CCLE
        ciri_fcrc <- suppressWarnings(cor(x = as.numeric(ciri_df[i,]), y = as.numeric(fcrc_df[i,]), method = "spearman")) #gCSI vs GDSC
        cfnd_fcrc <- suppressWarnings(cor(x = as.numeric(cfnd_df[i,]), y = as.numeric(fcrc_df[i,]), method = "spearman")) #GDSC vs CCLE
    
        # combine results
        correlations <- rbind(correlations, c(ciri_cfnd, ciri_fcrc, cfnd_fcrc))
    }

    correlations <- correlations[complete.cases(correlations), ]
    return(correlations)
}


# compute correlations 
polyA_corr <- corr_reps(polyA_df)
ribo0_corr <- corr_reps(ribo0_df)


############################################################
# Plot Spearman corr
############################################################

# set palette for plotting
pal = c("#23022E", "#611C35", "#839788", "#BFD7EA", "#BA9790", "#D5BC8A")

# function to format correlation dataframes
format_df <- function(df, label) {
    toPlot <- reshape2::melt(df)
    colnames(toPlot) <- c("Pipeline", "Correlation")
    toPlot$label <- label
    return(toPlot)
}

polyA_corr <- format_df(polyA_corr, "poly(A)-selection")
ribo0_corr <- format_df(ribo0_corr, "RiboZero")

# merge results for plotting
toPlot <- rbind(polyA_corr, ribo0_corr)
toPlot[is.na(toPlot)] <- 0
toPlot$label <- factor(toPlot$label, levels = c("poly(A)-selection", "RiboZero"))

png("../results/figures/figure7/correlate_reps_lungs.png", width=200, height=75, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = Correlation, y = label)) + 
    geom_violin(aes(fill = label, alpha = Pipeline), 
                    scale = "width", width = 1.1, 
                    position = position_dodge(width = 0.8)) + 
    theme_classic() + 
    labs(x = "Spearman Correlation", fill = "", alpha = "", y = "") +
    scale_fill_manual(values = pal) + 
    scale_alpha_manual(values = c(0.2, 0.5, 0.9)) +
    scale_y_discrete(limits = c("RiboZero", "poly(A)-selection")) + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3), legend.key.size = unit(0.5, 'cm')) 
dev.off()

png("../results/figures/figure5/corr_dist_density_ge.png", width=150, height=50, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = Correlation)) + 
    geom_density(aes(fill = label), alpha = 0.4, size = 0.5) + 
    theme_classic() +  
    scale_fill_manual(values = pal) + 
    labs(x = "Spearman Correlation", y = "Density", fill = "") +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3), legend.key.size = unit(0.4, 'cm'))
dev.off()



############################################################
# Create UMAP projections
############################################################

# function to create cell line umap projections
umap_fn <- function(expr_df, cell_line) {
    umap_df <- umap(expr_df)
    umap_df <- as.data.frame(umap_df$layout)
    colnames(umap_df) <- c("UMAP1", "UMAP2")
    umap_df$cell_line <- cell_line$Sample
    umap_df$dataset <- cell_line$PSet
    return(umap_df)
}

ciri_umap <- umap_fn(ciri_df, cell_line_ciri) 
circ_umap <- umap_fn(circ_df, cell_line_circ) 
cfnd_umap <- umap_fn(cfnd_df, cell_line_cfnd) 
fcrc_umap <- umap_fn(fcrc_df, cell_line_fcrc) 

# load in gene expression and isoform df
load("../results/data/umapdf.RData")

cell_line_gexpr <- data.frame(Sample = cell_line_gexpr, PSet = rep(c("gCSI", "CCLE", "GDSC"), each = 48))
cell_line_isoforms <- data.frame(Sample = cell_line_isoforms, PSet = rep(c("gCSI", "CCLE", "GDSC"), each = 48))

gexpr_umap <- umap_fn(gexpr_df, cell_line_gexpr) 
isoform_umap <- umap_fn(isoform_df, cell_line_isoforms)

############################################################
# Plot UMAP projections
############################################################

# function to plot umap
plot_umap <- function(umap_df, title) {
  p <- ggplot(data = umap_df, aes(x = UMAP1, y = UMAP2, group = cell_line)) + 
          geom_line(show.legend = F) + 
          geom_point(aes(color = dataset, shape = dataset), size = 3) + 
          scale_color_manual(guide = guide_legend(reverse = FALSE, title = "Dataset"), 
                            labels=c("CCLE", "gCSI", "GDSC2"), 
                            values = c("#392C57", "#4CC5AB", "#3670A0")) + 
          guides(shape = guide_legend(ncol = 1)) + #make legend one column instead of 3
          theme_classic() +
          theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
                text = element_text(size = 15), 
                legend.key.size = unit(0.7, 'cm'),
                plot.title = element_text(hjust = 0.5, size = 18), 
                axis.text.x = element_text(size=15, vjust = 0.5), 
                axis.text.y = element_text(size=15)) +
          labs(title = title) +
          scale_shape_discrete(name = "Dataset", labels = c("CCLE", "gCSI", "GDSC2"))
  return(p)
}

p1 <- plot_umap(gexpr_umap, "Gene Expression")
p2 <- plot_umap(isoform_umap, "Isoform Expression")
p3 <- plot_umap(ciri_umap, "CIRI2 circRNA Expression")
p4 <- plot_umap(circ_umap, "CIRCexplorer2 circRNA Expression")
p5 <- plot_umap(cfnd_umap, "circRNA_finder circRNA Expression")
p6 <- plot_umap(fcrc_umap, "find_circ circRNA Expression")


png("../results/figures/figure5/umaps_ge.png", width=400, height=225, units='mm', res = 600, pointsize=80)
ggarrange(p1, p3, p5, p2, p4, p6,
          ncol = 3, nrow = 2,
          common.legend = TRUE,
          legend = "right")
dev.off()

############################################################
# Plot UMAP projection distances
############################################################

# function to compute euclidean distances
compute_dist <- function(umap_df, label) {

    # compute euclidean distance across replicates
    distances <- umap_df %>%
      group_by(cell_line) %>%
      summarize(euclidean_dist = dist(cbind(UMAP1, UMAP2)))

    distances$label <- label
    return(distances)
}

gexpr_dist <- compute_dist(gexpr_umap, "Gene Expression") |> suppressWarnings()
isoform_dist <- compute_dist(isoform_umap, "Isoforms") |> suppressWarnings()
ciri_dist <- compute_dist(ciri_umap, "CIRI2") |> suppressWarnings()
circ_dist <- compute_dist(circ_umap, "CIRCexplorer2") |> suppressWarnings()
cfnd_dist <- compute_dist(cfnd_umap, "circRNA_finder") |> suppressWarnings()
fcrc_dist <- compute_dist(fcrc_umap, "find_circ") |> suppressWarnings()

# format dataframe for plotting
toPlot <- rbind(gexpr_dist, isoform_dist, ciri_dist, circ_dist, cfnd_dist, fcrc_dist)
toPlot$label <- factor(toPlot$label, levels = c("Gene Expression", "Isoforms", "CIRI2", "CIRCexplorer2", "circRNA_finder", "find_circ"))

# plot density plot
png("../results/figures/figure5/umap_dist_density.png", width=150, height=100, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = euclidean_dist)) + geom_density(aes(fill = label), alpha = 0.4, size = 0.5) + 
        theme_classic() +  
        scale_fill_manual(values = c("#23022E", "#611C35", "#839788", "#BFD7EA", "#BA9790", "#D5BC8A")) +
        labs(x = "Euclidean Distance of UMAP Points", y = "Density") +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
            legend.key.size = unit(0.4, 'cm'))
dev.off()

# plot violin plots
png("../results/figures/figure5/umap_dist_boxplot.png", width=150, height=100, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = label, y = euclidean_dist)) + 
    geom_violin(aes(fill = label), alpha = 0.8) + geom_boxplot(width=0.1, alpha = 0.4) +
    theme_classic() + labs(x = "", fill = "", y = "Euclidean Distance of UMAP Points") +
    scale_fill_manual(values = c("#23022E", "#611C35", "#839788", "#BFD7EA", "#BA9790", "#D5BC8A")) +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3), legend.position = "none") 
dev.off()

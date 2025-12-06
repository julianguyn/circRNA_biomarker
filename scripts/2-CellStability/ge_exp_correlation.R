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

# save labels for cell lines
cell_line_ciri <- data.frame(Sample = c(ciri_gcsi$sample, ciri_ccle$sample, ciri_gdsc$sample), PSet = rep(c("gCSI", "CCLE", "GDSC"), each = 48))
cell_line_circ <- data.frame(Sample = c(circ_gcsi$sample, circ_ccle$sample, circ_gdsc$sample), PSet = rep(c("gCSI", "CCLE", "GDSC"), each = 48))
cell_line_cfnd <- data.frame(Sample = c(cfnd_gcsi$sample, cfnd_ccle$sample, cfnd_gdsc$sample), PSet = rep(c("gCSI", "CCLE", "GDSC"), each = 48))
cell_line_fcrc <- data.frame(Sample = c(fcrc_gcsi$sample, fcrc_ccle$sample, fcrc_gdsc$sample), PSet = rep(c("gCSI", "CCLE", "GDSC"), each = 48))

# remove cell line labels
ciri_gcsi$sample <- ciri_ccle$sample <- ciri_gdsc$sample <- NULL
circ_gcsi$sample <- circ_ccle$sample <- circ_gdsc$sample <- NULL
cfnd_gcsi$sample <- cfnd_ccle$sample <- cfnd_gdsc$sample <- NULL
fcrc_gcsi$sample <- fcrc_ccle$sample <- fcrc_gdsc$sample <- NULL


############################################################
# Filter matrices by transcripts
############################################################

# filter transcripts with low expression across samples
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

# get common transcripts
transcripts <- data.frame(c(colnames(ciri_gcsi), colnames(ciri_ccle), colnames(ciri_gdsc),
                            colnames(circ_gcsi), colnames(circ_ccle), colnames(circ_gdsc),
                            colnames(cfnd_gcsi), colnames(cfnd_ccle), colnames(cfnd_gdsc),
                            colnames(fcrc_gcsi), colnames(fcrc_ccle), colnames(fcrc_gdsc)))
colnames(transcripts) <- "transcriptID"
transcript_counts <- transcripts %>% count(transcriptID)

# remove transcripts only in one method_dataset object (keep any in 2 or more)
transcript_counts <- transcript_counts[-which(transcript_counts$n == 1),] 
common_transcripts <- transcript_counts$transcriptID                                                    #n: 29764


############################################################
# Create merged dataframes 
############################################################

# function to create merged dataframe for each pipeline
mergePSet <- function(gcsi_df, ccle_df, gdsc_df) {
  # keep only common transcripts
  gcsi_filtered <- gcsi_df[,which(colnames(gcsi_df) %in% common_transcripts)]
  ccle_filtered <- ccle_df[,which(colnames(ccle_df) %in% common_transcripts)]
  gdsc_filtered <- gdsc_df[,which(colnames(gdsc_df) %in% common_transcripts)]

  # merge into one dataframe
  df <- rbind.fill(gcsi_filtered, ccle_filtered, gdsc_filtered)
  df[is.na(df)] <- 0
  df[] <- lapply(df, as.double)

  return(df)
}

ciri_df <- mergePSet(ciri_gcsi, ciri_ccle, ciri_gdsc)
circ_df <- mergePSet(circ_gcsi, circ_ccle, circ_gdsc)
cfnd_df <- mergePSet(cfnd_gcsi, cfnd_ccle, cfnd_gdsc)
fcrc_df <- mergePSet(fcrc_gcsi, fcrc_ccle, fcrc_gdsc)

# save all dataframes
save(ciri_df, circ_df, cfnd_df, fcrc_df, cell_line_ciri, cell_line_circ, cell_line_cfnd, cell_line_fcrc,
     file="../results/data/ge_exp_corr_df.RData")

############################################################
# Spearman corr of transcript exp across biological reps
############################################################

# function to parse processed dataframe and perform pairwise correlations
corr_reps <- function(exp_df, samples) {

    # parse datasets
    gcsi_df <- exp_df[1:48,]
    ccle_df <- exp_df[49:96,]
    gdsc_df <- exp_df[97:144,]
    
    # add cell line labels
    rownames(gcsi_df) <- samples$Sample[samples$PSet == "gCSI"]
    rownames(ccle_df) <- samples$Sample[samples$PSet == "CCLE"]
    rownames(gdsc_df) <- samples$Sample[samples$PSet == "GDSC"]

    # order cell lines
    gcsi_df <- gcsi_df[order(rownames(gcsi_df)),]
    ccle_df <- ccle_df[order(rownames(ccle_df)),]
    gdsc_df <- gdsc_df[order(rownames(gdsc_df)),]

    # initiate dataframe for storing results
    correlations <- data.frame(matrix(nrow=0, ncol=4))
    colnames(correlations) <- c("Cells", "gCSI_CCLE", "gCSI_GDSC", "GDSC_CCLE")

    # loop through each common transcript
    for (i in 1:nrow(gcsi_df)) {

        cell <- rownames(gcsi_df)[i]

        # compute correlations of transcript expression for pairs of psets
        gcsi_ccle <- suppressWarnings(cor(x = as.numeric(gcsi_df[i,]), y = as.numeric(ccle_df[i,]), method = "spearman")) #gCSI vs CCLE
        gcsi_gdsc <- suppressWarnings(cor(x = as.numeric(gcsi_df[i,]), y = as.numeric(gdsc_df[i,]), method = "spearman")) #gCSI vs GDSC
        ccle_gdsc <- suppressWarnings(cor(x = as.numeric(ccle_df[i,]), y = as.numeric(gdsc_df[i,]), method = "spearman")) #GDSC vs CCLE
    
        # combine results
        correlations <- rbind(correlations, 
                            data.frame(Cells = cell, gCSI_CCLE = gcsi_ccle, gCSI_GDSC = gcsi_gdsc, GDSC_CCLE = ccle_gdsc))
    }

    return(correlations)
}

# load in gene expression and isoform correlations
load("../results/data/corr_expr.RData")

# compute correlations 
ciri_corr <- corr_reps(ciri_df, cell_line_ciri)
circ_corr <- corr_reps(circ_df, cell_line_circ)
cfnd_corr <- corr_reps(cfnd_df, cell_line_cfnd)
fcrc_corr <- corr_reps(fcrc_df, cell_line_fcrc)

save(gexpr_corr, isoform_corr, ciri_corr, circ_corr, cfnd_corr, fcrc_corr, file = "../results/data/ge_corr_expr.RData")


############################################################
# Compute average spearman correlations (Table 3)
############################################################

colMeans(gexpr_corr[,2:4])
colMeans(isoform_corr[,2:4])
colMeans(ciri_corr[,2:4], na.rm = TRUE)
colMeans(circ_corr[,2:4], na.rm = TRUE)
colMeans(cfnd_corr[,2:4], na.rm = TRUE)
colMeans(fcrc_corr[,2:4], na.rm = TRUE)

############################################################
# Plot Spearman corr
############################################################

# set palette for plotting
pal = c("#23022E", "#611C35", "#839788", "#BFD7EA", "#BA9790", "#D5BC8A")

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

png("../results/figures/figure5/correlate_reps_ge.png", width=200, height=75, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = Correlation, y = label)) + 
    geom_violin(aes(fill = label, alpha = PSet), 
                    scale = "width", width = 1.1, 
                    position = position_dodge(width = 0.8)) + 
    theme_classic() + 
    labs(x = "Spearman Correlation", fill = "", alpha = "", y = "") +
    scale_fill_manual(values = pal) + 
    scale_alpha_manual(values = c(0.2, 0.5, 0.9)) +
    scale_y_discrete(limits = c("find_circ", "circRNA_finder", "CIRCexplorer2", "CIRI2", "Isoforms", "Gene Expression")) + 
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

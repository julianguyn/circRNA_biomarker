suppressMessages(library(data.table))
suppressMessages(library(plyr))
suppressMessages(library(tidyverse))
suppressMessages(library(umap))  
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(dplyr))
set.seed(101)


###########################################
######### Load in Expression Data #########
###########################################


# ========== Process isoform and gene expression matrices for UMAP projection ========== #

# load in isoform and gene expression data
load("../results/data/gene_expression.RData")
load("../results/data/isoform_expression.RData")

# save labels for cell lines
cell_line_gexpr <- c(colnames(expr_gcsi_p), colnames(expr_ccle_p), colnames(expr_gdsc_p))
cell_line_isoforms <- c(colnames(expr_gcsi_i), colnames(expr_ccle_i), colnames(expr_gdsc_i))

# format gene expression data
expr_gcsi_p <- as.data.frame(t(expr_gcsi_p))
expr_ccle_p <- as.data.frame(t(expr_ccle_p))
expr_gdsc_p <- as.data.frame(t(expr_gdsc_p))

rownames(expr_gcsi_p) <- NULL
rownames(expr_ccle_p) <- NULL
rownames(expr_gdsc_p) <- NULL


# format isoform expression data
expr_gcsi_i <- as.data.frame(t(expr_gcsi_i))
expr_ccle_i <- as.data.frame(t(expr_ccle_i))
expr_gdsc_i <- as.data.frame(t(expr_gdsc_i))

rownames(expr_gcsi_i) <- NULL
rownames(expr_ccle_i) <- NULL
rownames(expr_gdsc_i) <- NULL


# merge all expression values from each data set 
gexpr_df <- rbind(expr_gcsi_p, expr_ccle_p, expr_gdsc_p)
isoform_df <- rbind(expr_gcsi_i, expr_ccle_i, expr_gdsc_i)


# ========== Process cell line circRNA matrices for UMAP projection ========== #

# load circRNA expression data
ciri_gcsi_sub <- fread("../data/processed_cellline/common_samples/CIRI2/ciri_gcsi_counts.tsv", data.table = F)
ciri_gdsc_sub <- fread("../data/processed_cellline/common_samples/CIRI2/ciri_gdsc_counts.tsv", data.table = F)
ciri_ccle_sub <- fread("../data/processed_cellline/common_samples/CIRI2/ciri_ccle_counts.tsv", data.table = F)

circ_gcsi_sub <- fread("../data/processed_cellline/common_samples/CIRCexplorer2/circ_gcsi_counts.tsv", data.table = F)
circ_gdsc_sub <- fread("../data/processed_cellline/common_samples/CIRCexplorer2/circ_gdsc_counts.tsv", data.table = F)
circ_ccle_sub <- fread("../data/processed_cellline/common_samples/CIRCexplorer2/circ_ccle_counts.tsv", data.table = F)

cfnd_gcsi_sub <- fread("../data/processed_cellline/common_samples/circRNA_finder/cfnd_gcsi_counts.tsv", data.table = F)
cfnd_gdsc_sub <- fread("../data/processed_cellline/common_samples/circRNA_finder/cfnd_gdsc_counts.tsv", data.table = F)
cfnd_ccle_sub <- fread("../data/processed_cellline/common_samples/circRNA_finder/cfnd_ccle_counts.tsv", data.table = F)

fcrc_gcsi_sub <- fread("../data/processed_cellline/common_samples/find_circ/fcrc_gcsi_counts.tsv", data.table = F)
fcrc_gdsc_sub <- fread("../data/processed_cellline/common_samples/find_circ/fcrc_gdsc_counts.tsv", data.table = F)
fcrc_ccle_sub <- fread("../data/processed_cellline/common_samples/find_circ/fcrc_ccle_counts.tsv", data.table = F)


# save labels for cell lines
cell_line_ciri <- c(ciri_gcsi_sub$sample, ciri_ccle_sub$sample, ciri_gdsc_sub$sample)
cell_line_circ <- c(circ_gcsi_sub$sample, circ_ccle_sub$sample, circ_gdsc_sub$sample)
cell_line_cfnd <- c(cfnd_gcsi_sub$sample, cfnd_ccle_sub$sample, cfnd_gdsc_sub$sample)
cell_line_fcrc <- c(fcrc_gcsi_sub$sample, fcrc_ccle_sub$sample, fcrc_gdsc_sub$sample)


# remove cell line labels
ciri_gcsi_sub$sample <- ciri_ccle_sub$sample <- ciri_gdsc_sub$sample <- NULL
circ_gcsi_sub$sample <- circ_ccle_sub$sample <- circ_gdsc_sub$sample <- NULL
cfnd_gcsi_sub$sample <- cfnd_ccle_sub$sample <- cfnd_gdsc_sub$sample <- NULL
fcrc_gcsi_sub$sample <- fcrc_ccle_sub$sample <- fcrc_gdsc_sub$sample <- NULL


# filter circRNA transcripts with low detection rates
# distribution: table(colSums(ciri_gcsi_sub == 0)), shows the number of 0 in each column

ciri_gcsi_filtered <- ciri_gcsi_sub[,-which(colnames(ciri_gcsi_sub) %in% names(which(colSums(ciri_gcsi_sub == 0) > 45)))]
ciri_ccle_filtered <- ciri_ccle_sub[,-which(colnames(ciri_ccle_sub) %in% names(which(colSums(ciri_ccle_sub == 0) > 45)))]
ciri_gdsc_filtered <- ciri_gdsc_sub[,-which(colnames(ciri_gdsc_sub) %in% names(which(colSums(ciri_gdsc_sub == 0) > 45)))]

print(dim(ciri_gcsi_filtered))
print(dim(ciri_ccle_filtered))
print(dim(ciri_gdsc_filtered))

circ_gcsi_filtered <- circ_gcsi_sub[,-which(colnames(circ_gcsi_sub) %in% names(which(colSums(circ_gcsi_sub == 0) > 45)))]
circ_ccle_filtered <- circ_ccle_sub[,-which(colnames(circ_ccle_sub) %in% names(which(colSums(circ_ccle_sub == 0) > 45)))]
circ_gdsc_filtered <- circ_gdsc_sub[,-which(colnames(circ_gdsc_sub) %in% names(which(colSums(circ_gdsc_sub == 0) > 45)))]

print(dim(circ_gcsi_filtered))
print(dim(circ_ccle_filtered))
print(dim(circ_gdsc_filtered))

cfnd_gcsi_filtered <- cfnd_gcsi_sub[,-which(colnames(cfnd_gcsi_sub) %in% names(which(colSums(cfnd_gcsi_sub == 0) > 45)))]
cfnd_ccle_filtered <- cfnd_ccle_sub[,-which(colnames(cfnd_ccle_sub) %in% names(which(colSums(cfnd_ccle_sub == 0) > 45)))]
cfnd_gdsc_filtered <- cfnd_gdsc_sub[,-which(colnames(cfnd_gdsc_sub) %in% names(which(colSums(cfnd_gdsc_sub == 0) > 45)))]

print(dim(cfnd_gcsi_filtered))
print(dim(cfnd_ccle_filtered))
print(dim(cfnd_gdsc_filtered))

fcrc_gcsi_filtered <- fcrc_gcsi_sub[,-which(colnames(fcrc_gcsi_sub) %in% names(which(colSums(fcrc_gcsi_sub == 0) > 45)))]
fcrc_ccle_filtered <- fcrc_ccle_sub[,-which(colnames(fcrc_ccle_sub) %in% names(which(colSums(fcrc_ccle_sub == 0) > 45)))]
fcrc_gdsc_filtered <- fcrc_gdsc_sub[,-which(colnames(fcrc_gdsc_sub) %in% names(which(colSums(fcrc_gdsc_sub == 0) > 45)))]

print(dim(fcrc_gcsi_filtered))
print(dim(fcrc_ccle_filtered))
print(dim(fcrc_gdsc_filtered))


# get common circRNA transcripts
transcripts <- data.frame(c(colnames(ciri_gcsi_filtered), colnames(ciri_ccle_filtered), colnames(ciri_gdsc_filtered),
                            colnames(circ_gcsi_filtered), colnames(circ_ccle_filtered), colnames(circ_gdsc_filtered),
                            colnames(cfnd_gcsi_filtered), colnames(cfnd_ccle_filtered), colnames(cfnd_gdsc_filtered),
                            colnames(fcrc_gcsi_filtered), colnames(fcrc_ccle_filtered), colnames(fcrc_gdsc_filtered)))
colnames(transcripts) <- "transcriptID"
transcript_counts <- transcripts %>% count(transcriptID)

# remove circRNA transcripts that are only in one method_dataset object (keep any in 2 or more)
transcript_counts <- transcript_counts[-which(transcript_counts$n == 1),] 
common_transcripts <- transcript_counts$transcriptID

print(length(common_transcripts))

# create merged dataframe for each pipeline
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

ciri_df <- mergePSet(ciri_gcsi_filtered, ciri_ccle_filtered, ciri_gdsc_filtered)
circ_df <- mergePSet(circ_gcsi_filtered, circ_ccle_filtered, circ_gdsc_filtered)
cfnd_df <- mergePSet(cfnd_gcsi_filtered, cfnd_ccle_filtered, cfnd_gdsc_filtered)
fcrc_df <- mergePSet(fcrc_gcsi_filtered, fcrc_ccle_filtered, fcrc_gdsc_filtered)


# save all dataframes
save(gexpr_df, isoform_df, ciri_df, circ_df, cfnd_df, fcrc_df,
     cell_line_gexpr, cell_line_isoforms, cell_line_ciri, cell_line_circ, cell_line_cfnd, cell_line_fcrc,
     file="../results/data/umapdf.RData")



###########################################
######### Create UMAP Projections #########
###########################################


# ========== Create UMAP projections for cell line data ========== #

# function to create cell line umap projections
umap_fn <- function(expr_df, cell_line) {
  # umap
  umap_df <- umap(expr_df)

  # format umap df
  umap_df <- as.data.frame(umap_df$layout)
  colnames(umap_df) <- c("UMAP1", "UMAP2")
  umap_df$cell_line <- cell_line
  umap_df$dataset <- c(rep("gCSI", 48), rep("CCLE", 48), rep("GDSC", 48))

  return(umap_df)
}

# create umap projections
gexpr_umap <- umap_fn(gexpr_df, cell_line_gexpr) #gene expression
isoform_umap <- umap_fn(isoform_df, cell_line_isoforms) #isoforms
ciri_umap <- umap_fn(ciri_df, cell_line_ciri) #CIRI2 circRNA
circ_umap <- umap_fn(circ_df, cell_line_circ) #CIRCexplorer2 circRNA
cfnd_umap <- umap_fn(cfnd_df, cell_line_cfnd) #CIRI2 circRNA
fcrc_umap <- umap_fn(fcrc_df, cell_line_fcrc) #CIRCexplorer2 circRNA



##############################################
######### Visualize UMAP Projections #########
##############################################

# ========== Plot UMAP projections for cell line data ========== #

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


png("../results/figures/figure2/umaps.png", width=400, height=225, units='mm', res = 600, pointsize=80)
ggarrange(p1, p3, p5, p2, p4, p6,
          ncol = 3, nrow = 2,
          common.legend = TRUE,
          legend = "right")
dev.off()


# ========== Plot distances between replicates on UMAP projection ========== #

# function to compute euclidean distances
compute_dist <- function(umap_df, label) {

    # compute euclidean distance across replicates
    distances <- umap_df %>%
      group_by(cell_line) %>%
      summarize(euclidean_dist = dist(cbind(UMAP1, UMAP2)))

    distances$label <- label
    return(distances)
}

suppressWarnings(gexpr_dist <- compute_dist(gexpr_umap, "Gene Expression"))
suppressWarnings(isoform_dist <- compute_dist(isoform_umap, "Isoforms"))
suppressWarnings(ciri_dist <- compute_dist(ciri_umap, "CIRI2"))
suppressWarnings(circ_dist <- compute_dist(circ_umap, "CIRCexplorer2"))
suppressWarnings(cfnd_dist <- compute_dist(cfnd_umap, "circRNA_finder"))
suppressWarnings(fcrc_dist <- compute_dist(fcrc_umap, "find_circ"))

# format dataframe for plotting
toPlot <- rbind(gexpr_dist, isoform_dist, ciri_dist, circ_dist, cfnd_dist, fcrc_dist)
toPlot$label <- factor(toPlot$label, levels = c("Gene Expression", "Isoforms", "CIRI2", "CIRCexplorer2", "circRNA_finder", "find_circ"))

# plot density plot
png("../results/figures/figure2/umap_dist_density.png", width=150, height=100, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = euclidean_dist)) + geom_density(aes(fill = label), alpha = 0.4, size = 0.5) + 
        theme_classic() +  
        scale_fill_manual(values = c("#23022E", "#611C35", "#839788", "#BFD7EA", "#BA9790", "#D5BC8A")) +
        labs(x = "Euclidean Distance of UMAP Points", y = "Density") +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
            legend.key.size = unit(0.4, 'cm'))
dev.off()

# plot violin plots
png("../results/figures/figure2/umap_dist_boxplot.png", width=150, height=100, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = label, y = euclidean_dist)) + 
    geom_violin(aes(fill = label), alpha = 0.8) + geom_boxplot(width=0.1, alpha = 0.4) +
    theme_classic() + labs(x = "", fill = "", y = "Euclidean Distance of UMAP Points") +
    scale_fill_manual(values = c("#23022E", "#611C35", "#839788", "#BFD7EA", "#BA9790", "#D5BC8A")) +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3), legend.position = "none") 
dev.off()


############################################################################################################

# ========== Process tumour circRNA matrices for UMAP projection ========== #

# load in lung circRNA expression data 
"../results/data/circ_lung_expression.RData"

rownames(lung_ribozero_ciri) <- NULL
rownames(lung_polyA_ciri) <- NULL
rownames(lung_ribozero_circ) <- NULL
rownames(lung_polyA_circ) <- NULL

lung_ribozero_ciri <- as.data.frame(lung_ribozero_ciri)
lung_polyA_ciri <- as.data.frame(lung_polyA_ciri)
lung_ribozero_circ <- as.data.frame(lung_ribozero_circ)
lung_polyA_circ <- as.data.frame(lung_polyA_circ)

lung_df <- rbind.fill(lung_ribozero_ciri, lung_polyA_ciri, lung_ribozero_circ, lung_polyA_circ)
lung_df[is.na(lung_df)] <- 0
lung_df[] <- lapply(lung_df, as.double)


# save all dataframes
save(gexpr_df, isoform_df, ciri_df, circ_df, lung_df,
     cell_line_gexpr, cell_line_isoforms, cell_line_ciri, cell_line_circ,
     file="../results/umapdf.RData")


# ========== Create UMAP projection for lung tumor data ========== #

# create lung umap projection
lung_umap <- umap(lung_df)

# format data frame
lung_umap <- as.data.frame(lung_umap$layout)
colnames(lung_umap) <- c("UMAP1", "UMAP2")
lung_umap$tumourID <- c(rep(paste0("tumor",1:51), 4))
lung_umap$Selection <- c(rep("poly(A)", 51), rep("RiboZero", 51), rep("poly(A)", 51), rep("RiboZero", 51))
lung_umap$Method <- c(rep("CIRI2", 102), rep("CIRCexplorer2", 102))


# ========== Plot UMAP projections for lung tumor data ========== #

png("../results/supplementary_figure7.png", width=150, height=125, units='mm', res = 600, pointsize=80)
ggplot(data = lung_umap, aes(x = UMAP1, y = UMAP2)) + 
  geom_point(aes(color = Selection, shape = Method), size = 3) + 
  guides(shape = guide_legend(ncol = 1), color = guide_legend(override.aes=list(shape=15, size = 8))) + #make legend one column instead of 3, change shape to square
  scale_color_manual(guide = guide_legend(reverse = FALSE), labels=c("poly(A)", "RiboZero"), values = c("#4CC5AB", "#392C57")) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "Lung Samples circRNA Expression", x = "UMAP1 \n", y = "\nUMAP2") 
dev.off()




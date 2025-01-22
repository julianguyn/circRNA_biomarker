# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(matrixStats)
    library(reshape2)
    library(ggplot2)
    library(ggpubr)
})


############################################################
# Load in circRNA expression data
############################################################

# load circRNA expression data
path = "../data/processed_cellline/common_samples/" 

ciri_gcsi <- fread(paste0(path, "CIRI2/ciri_gcsi_counts.tsv", data.table = F))
ciri_gdsc <- fread(paste0(path, "CIRI2/ciri_gdsc_counts.tsv", data.table = F))
ciri_ccle <- fread(paste0(path, "CIRI2/ciri_ccle_counts.tsv", data.table = F))

circ_gcsi <- fread(paste0(path, "CIRCexplorer2/circ_gcsi_counts.tsv", data.table = F))
circ_gdsc <- fread(paste0(path, "CIRCexplorer2/circ_gdsc_counts.tsv", data.table = F))
circ_ccle <- fread(paste0(path, "CIRCexplorer2/circ_ccle_counts.tsv", data.table = F))

cfnd_gcsi <- fread(paste0(path, "circRNA_finder/cfnd_gcsi_counts.tsv", data.table = F))
cfnd_gdsc <- fread(paste0(path, "circRNA_finder/cfnd_gdsc_counts.tsv", data.table = F))
cfnd_ccle <- fread(paste0(path, "circRNA_finder/cfnd_ccle_counts.tsv", data.table = F))

fcrc_gcsi <- fread(paste0(path, "find_circ/fcrc_gcsi_counts.tsv", data.table = F))
fcrc_gdsc <- fread(paste0(path, "find_circ/fcrc_gdsc_counts.tsv", data.table = F))
fcrc_ccle <- fread(paste0(path, "find_circ/fcrc_ccle_counts.tsv", data.table = F))


############################################################
# Filter circRNA transcripts with low detection rates
############################################################
# distribution: table(colSums(ciri_gcsi_sub == 0)), shows the number of 0 in each column

# function to filter low detection transcripts
filter_circ <- function(df) {
    df <- df[,-which(colnames(df) %in% names(which(colSums(df == 0) > 45)))]
    print(dim(df))
    return(df)
}

ciri_gcsi <- filter_circ(ciri_gcsi)
ciri_ccle <- filter_circ(ciri_ccle)
ciri_gdsc <- filter_circ(ciri_gdsc)

circ_gcsi <- filter_circ(circ_gcsi)
circ_ccle <- filter_circ(circ_ccle)
circ_gdsc <- filter_circ(circ_gdsc)

cfnd_gcsi <- filter_circ(cfnd_gcsi)
cfnd_ccle <- filter_circ(cfnd_ccle)
cfnd_gdsc <- filter_circ(cfnd_gdsc)

fcrc_gcsi <- filter_circ(fcrc_gcsi)
fcrc_ccle <- filter_circ(fcrc_ccle)
fcrc_gdsc <- filter_circ(fcrc_gdsc)


############################################################
# Get common transcripts
############################################################

# get common circRNA transcripts found in all psets for each pipeline
ciri_common <- intersect(intersect(colnames(ciri_gcsi), colnames(ciri_ccle)), colnames(ciri_gdsc))
circ_common <- intersect(intersect(colnames(circ_gcsi), colnames(circ_ccle)), colnames(circ_gdsc))
cfnd_common <- intersect(intersect(colnames(cfnd_gcsi), colnames(cfnd_ccle)), colnames(cfnd_gdsc))
fcrc_common <- intersect(intersect(colnames(fcrc_gcsi), colnames(fcrc_ccle)), colnames(fcrc_gdsc))

print(length(ciri_common))
print(length(circ_common))
print(length(cfnd_common))
print(length(fcrc_common))


############################################################
# Keep only common transcripts
############################################################

# function to order circRNA dataframes and keep only transcripts found in all psets for each pipeline
subset_df <- function(df, common_transcripts) {
    df <- df[,which(colnames(df) %in% common_transcripts)]
    df <- df[order(df$sample),order(colnames(df))]
    rownames(df) <- df$sample
    df$sample <- NULL
    return(df)
}

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

# save subsetted dataframes for featuer importance
save(ciri_gcsi, ciri_ccle, ciri_gdsc, circ_gcsi, circ_ccle, circ_gdsc,
     cfnd_gcsi, cfnd_ccle, cfnd_gdsc, fcrc_gcsi, fcrc_ccle, fcrc_gdsc,
     file = "../results/data/temp/circ_stability_subsetdf.RData")



############################################################
# Compute pairwise spearman correlations
############################################################

# function to compute pairwise spearman correlations
compute_spearman <- function(
    gcsi_df, ccle_df, gdsc_df,      # PSet-specific dataframes from the subset_df() function
    random = FALSE,                 # random: TRUE for random sampling of sample names, FALSE otherwise
    iter = 1                        # iter: number of iterations (for random sampling)
) {
    # initialize dataframe to store results
    correlations <- data.frame(matrix(nrow=0, ncol=3))
    for (i in 1:iter) {

        if (random == TRUE) {
            # randomly shuffle cell line names
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
     file = "../results/data/temp/circ_stability.RData")

save(ciri_stability_random, circ_stability_random, cfnd_stability_random, fcrc_stability_random, 
     file = "../results/data/temp/circ_stability_random.RData")



############################################################
# Load isoform and gene expression data
############################################################

load("../results/data/isoform_expression.RData")
load("../results/data/gene_expression.RData")


############################################################
# Compute isoform stability
############################################################

# function to compute spearman correlation for each dataset pair (gCSI/CCLE; gCSI/GDSC; GDSC/CCLE)                                 
compute_spearman  <- function(x){
  i <- x
  gcsi_ccle_spearman <- suppressWarnings(cor(x=as.numeric(expr_gcsi_i[i,]), y=as.numeric(expr_ccle_i[i,]), method="spearman"))
  gcsi_gdsc_spearman <- suppressWarnings(cor(x=as.numeric(expr_gcsi_i[i,]), y=as.numeric(expr_gdsc_i[i,]), method="spearman"))
  gdsc_ccle_spearman <- suppressWarnings(cor(x=as.numeric(expr_gdsc_i[i,]), y=as.numeric(expr_ccle_i[i,]), method="spearman"))
  
  combined <- c(gcsi_ccle_spearman, gcsi_gdsc_spearman, gdsc_ccle_spearman)
  return(combined)
}                                         

transcripts <- rownames(expr_gcsi_i)                                     

spearman_compute <- sapply(transcripts, compute_spearman)
transcript_stability <- as.data.frame(t(as.data.frame(spearman_compute)))
colnames(transcript_stability) <- c("gcsi_ccle_spearman", "gcsi_gdsc_spearman", "gdsc_ccle_spearman")  


############################################################
# Compute gene expression stability
############################################################

# function to compute spearman correlation for each dataset pair (gCSI/CCLE; gCSI/GDSC; GDSC/CCLE)                                 
compute_spearman  <- function(x){
  i <- x
  gcsi_ccle_spearman <- suppressWarnings(cor(x=as.numeric(expr_gcsi_p[i,]), y=as.numeric(expr_ccle_p[i,]), method="spearman"))
  gcsi_gdsc_spearman <- suppressWarnings(cor(x=as.numeric(expr_gcsi_p[i,]), y=as.numeric(expr_gdsc_p[i,]), method="spearman"))
  gdsc_ccle_spearman <- suppressWarnings(cor(x=as.numeric(expr_gdsc_p[i,]), y=as.numeric(expr_ccle_p[i,]), method="spearman"))
  
  combined <- c(gcsi_ccle_spearman, gcsi_gdsc_spearman, gdsc_ccle_spearman)
  return(combined)
}                                         

genes <- rownames(expr_gcsi_p)                                     

spearman_compute <- sapply(genes, compute_spearman)
gene_stability <- as.data.frame(t(as.data.frame(spearman_compute)))
colnames(gene_stability) <- c("gcsi_ccle_spearman", "gcsi_gdsc_spearman", "gdsc_ccle_spearman")  


############################################################
# Compute median expression for transcripts
############################################################

# format matrices
gcsi_matrix <- as.matrix(expr_gcsi_i)                                 
ccle_matrix <- as.matrix(expr_ccle_i)  
gdsc_matrix <- as.matrix(expr_gdsc_i)

gcsi_matrix <- as.matrix(expr_gcsi_p)                                 
ccle_matrix <- as.matrix(expr_ccle_p)  
gdsc_matrix <- as.matrix(expr_gdsc_p)


# compute median expression for each transcript across biological replicates for each dataset
gcsi_median <- as.numeric(rowMedians(gcsi_matrix))                                 
ccle_median <- as.numeric(rowMedians(ccle_matrix))                                  
gdsc_median <- as.numeric(rowMedians(gdsc_matrix))        

gcsi_median <- as.numeric(rowMedians(gcsi_matrix))                                 
ccle_median <- as.numeric(rowMedians(ccle_matrix))                                  
gdsc_median <- as.numeric(rowMedians(gdsc_matrix)) 


# compile isoforms into one data frame                                 
transcript_stability$transcript_id <- transcripts
rownames(transcript_stability) <- transcript_stability$transcript_id
transcript_stability$gcsi_median <- gcsi_median
transcript_stability$ccle_median <- ccle_median
transcript_stability$gdsc_median <- gdsc_median

#remove transcripts that have median expression of 0 across all datasets
transcript_stability <- transcript_stability[-which(transcript_stability$gcsi_median == 0 & transcript_stability$ccle_median == 0 & transcript_stability$gdsc_median == 0),]  
           
# compile gene expression into one data frame                                 
gene_stability$gene_id <- genes
rownames(gene_stability) <- gene_stability$gene_id
gene_stability$gcsi_median <- gcsi_median
gene_stability$ccle_median <- ccle_median
gene_stability$gdsc_median <- gdsc_median

save(gene_stability, transcript_stability, 
     file = "../results/data/temp/gene_isoform_stability.RData")

############################################################
# Plot comparison of SI indices
############################################################

# function to format stability dataframes
format_df <- function(df, label) {
  colnames(df) <- c("gCSI/CCLE", "gCSI/GDSC2", "GDSC2/CCLE")
  toPlot <- melt(df)
  colnames(toPlot) <- c("PSet", "Stability")
  toPlot$label <- label
  return(toPlot)
}

ciri_stability <- format_df(ciri_stability, "CIRI2")
circ_stability <- format_df(circ_stability, "CIRCexplorer2")
cfnd_stability <- format_df(cfnd_stability, "circRNA_finder")
fcrc_stability <- format_df(fcrc_stability, "find_circ")

transcript_stability <- format_df(transcript_stability[,c("gcsi_ccle_spearman", "gcsi_gdsc_spearman", "gdsc_ccle_spearman")], "Isoforms")
gene_stability <- format_df(gene_stability[,c("gcsi_ccle_spearman", "gcsi_gdsc_spearman", "gdsc_ccle_spearman")], "Gene Expression")
gene_stability[is.na(gene_stability)] <- 0


# merge results for plotting
toPlot <- rbind(gene_stability, transcript_stability, 
                ciri_stability, circ_stability, 
                cfnd_stability, fcrc_stability)

toPlot$label <- factor(toPlot$label, levels = c("Gene Expression", "Isoforms", "CIRI2", "CIRCexplorer2", "circRNA_finder", "find_circ"))


png("../results/figures/figure2/stability.png", width=300, height=150, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = label, y = Stability)) + 
    geom_violin(aes(fill = label), alpha = 0.8) + geom_boxplot(width=0.1, alpha = 0.3) +
    facet_grid(factor(PSet)~.) +
    theme_classic() + labs(x = "", fill = "", y = "Stability Index") +
    scale_fill_manual(values = c("#23022E", "#611C35", "#839788", "#BFD7EA", "#BA9790", "#D5BC8A")) +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3), legend.key.size = unit(0.7, 'cm')) +
    geom_hline(yintercept = 0, linetype = "dotted")
dev.off()


############################################################
# Plot comparison of SI indices with random shuffling
############################################################

# label dataframes
ciri_stability$label <- "CIRI2-NonRandom"
circ_stability$label <- "CIRCexplorer2-NonRandom"
cfnd_stability$label <- "circRNA_finder-NonRandom"
fcrc_stability$label <- "find_circ-NonRandom"

ciri_stability_random$label <- "CIRI2-Random"
circ_stability_random$label <- "CIRCexplorer2-Random"
cfnd_stability_random$label <- "circRNA_finder-Random"
fcrc_stability_random$label <- "find_circ-Random"

# format dataframe for plotting
toPlot <- rbind(ciri_stability, circ_stability, cfnd_stability, fcrc_stability,
                ciri_stability_random, circ_stability_random, cfnd_stability_random, fcrc_stability_random)
toPlot <- melt(toPlot)
toPlot$Label <- gsub(".*-", "", toPlot$label)
toPlot$PSet <- gsub("-.*", "", toPlot$label)
toPlot$PSet <- factor(toPlot$PSet, levels = c("CIRI2", "CIRCexplorer2", "circRNA_finder", "find_circ"))

# plot
png("../results/figures/figure2/stability_randomized.png", width=300, height=150, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = variable, y = value, fill = Label)) + 
    #geom_violin(aes(fill = Label), alpha = 0.8) + 
    geom_boxplot() + facet_grid(PSet~.) + theme_classic() + 
    labs(x = "", fill = "", y = "Stability Index") + scale_fill_manual(values = c("#839788", "gray")) +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3), legend.key.size = unit(0.7, 'cm')) +
    geom_hline(yintercept = 0, linetype = "dotted")
dev.off()


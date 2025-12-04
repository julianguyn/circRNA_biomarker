# revised to include transcripts present in at least 2/3 datasets instead of all 3

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

ciri_gcsi <- fread(paste0(path, "CIRI2/ciri_gcsi_counts.tsv"), data.table = F)
ciri_gdsc <- fread(paste0(path, "CIRI2/ciri_gdsc_counts.tsv"), data.table = F)
ciri_ccle <- fread(paste0(path, "CIRI2/ciri_ccle_counts.tsv"), data.table = F)

circ_gcsi <- fread(paste0(path, "CIRCexplorer2/circ_gcsi_counts.tsv"), data.table = F)
circ_gdsc <- fread(paste0(path, "CIRCexplorer2/circ_gdsc_counts.tsv"), data.table = F)
circ_ccle <- fread(paste0(path, "CIRCexplorer2/circ_ccle_counts.tsv"), data.table = F)

cfnd_gcsi <- fread(paste0(path, "circRNA_finder/cfnd_gcsi_counts.tsv"), data.table = F)
cfnd_gdsc <- fread(paste0(path, "circRNA_finder/cfnd_gdsc_counts.tsv"), data.table = F)
cfnd_ccle <- fread(paste0(path, "circRNA_finder/cfnd_ccle_counts.tsv"), data.table = F)

fcrc_gcsi <- fread(paste0(path, "find_circ/fcrc_gcsi_counts.tsv"), data.table = F)
fcrc_gdsc <- fread(paste0(path, "find_circ/fcrc_gdsc_counts.tsv"), data.table = F)
fcrc_ccle <- fread(paste0(path, "find_circ/fcrc_ccle_counts.tsv"), data.table = F)



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

# function to keep only transcripts in 2/3 of the psets 
get_common <- function(gcsi, ccle, gdsc) {
    all_transcripts <- c(colnames(gcsi), colnames(ccle), colnames(gdsc))
    freq <- table(all_transcripts)
    to_keep <- names(freq[freq >= 2])
    return(to_keep)
}

# get common circRNA transcripts found in all psets for each pipeline
ciri_common <- get_common(ciri_gcsi, ciri_ccle, ciri_gdsc)
circ_common <- get_common(circ_gcsi, circ_ccle, circ_gdsc)
cfnd_common <- get_common(cfnd_gcsi, cfnd_ccle, cfnd_gdsc)
fcrc_common <- get_common(fcrc_gcsi, fcrc_ccle, fcrc_gdsc)

print(length(ciri_common)) #138
print(length(circ_common)) #206
print(length(cfnd_common)) #607
print(length(fcrc_common)) #1569


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
     file = "../results/data/temp/circ_stability_subsetdf2.RData")



############################################################
# Compute pairwise spearman correlations
############################################################

# helper function to compute pairwise spearman correlations
pairwise_spearman <- function(pset1, pset2, label) {

    # get common transcripts
    common <- intersect(colnames(pset1), colnames(pset2))

    # initialize dataframe to store results
    correlations <- data.frame(transcript = common, si = NA)

    # keep common transcripts
    pset1 <- pset1[,which(colnames(pset1) %in% common)]
    pset2 <- pset2[,which(colnames(pset2) %in% common)]

    for (i in 1:length(common)) {
        # compute spearman
        sp <- suppressWarnings(cor(x = as.numeric(pset1[, i]), y = as.numeric(pset2[, i]), method = "spearman"))
        correlations$si[i] <- sp
    }
    correlations$pair <- label
    return(correlations)
}

# function to compute pairwise spearman correlations
compute_spearman <- function(
    gcsi_df, ccle_df, gdsc_df,      # PSet-specific dataframes from the subset_df() function
    random = FALSE,                 # random: TRUE for random sampling of sample names, FALSE otherwise
    iter = 1                        # iter: number of iterations (for random sampling)
) {

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

        # compute correlations of transcript expression for pairs of psets
        gcsi_ccle_spearman <- pairwise_spearman(gcsi_df, ccle_df, "gCSI/CCLE") 
        gcsi_gdsc_spearman <- pairwise_spearman(gcsi_df, gdsc_df, "gCSI/GDSC2") 
        gdsc_ccle_spearman <- pairwise_spearman(gdsc_df, ccle_df, "GDSC2/CCLE")
        
        # combine results
        correlations <- rbind(gcsi_ccle_spearman, gcsi_gdsc_spearman, gdsc_ccle_spearman)
        
    }
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
     file = "../results/data/temp/circ_stability2.RData")

save(ciri_stability_random, circ_stability_random, cfnd_stability_random, fcrc_stability_random, 
     file = "../results/data/temp/circ_stability_random2.RData")

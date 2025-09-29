# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(reshape2)
    library(ggplot2)
})


############################################################
# Specify analysis
############################################################

# analysis = c("circ", "GE")  
# circ for circRNA counts, GE for circ mapped to GE

if (analysis == "circ") {
    path <- "../data/processed_cellline/all_samples/" 
    out <- "all_circ"
    df <- "../data/processed_cellline/merged_all_samples/"
}
if (analysis == "GE") {
    path <- "../data/processed_cellline/GE_common_samples/"
    out <- "GE"
    df <- "../data/processed_cellline/mergedGE_common_samples/"
}


############################################################
# Load in circRNA expression data
############################################################

print("Reading in dataframes")

# load circRNA expression data
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
# Format circRNA dataframes
############################################################

print("Formating dataframes")

# function to return 
format_df <- function(circ_counts) {
    rownames(circ_counts) <- circ_counts$sample
    circ_counts$sample <- NULL
    circ_counts[is.na(circ_counts)] <- 0  

    return(circ_counts)
}

ciri_gcsi <- format_df(ciri_gcsi)
ciri_gdsc <- format_df(ciri_gdsc)
ciri_ccle <- format_df(ciri_ccle)

circ_gcsi <- format_df(circ_gcsi)
circ_gdsc <- format_df(circ_gdsc)
circ_ccle <- format_df(circ_ccle)

cfnd_gcsi <- format_df(cfnd_gcsi)
cfnd_gdsc <- format_df(cfnd_gdsc)
cfnd_ccle <- format_df(cfnd_ccle)

fcrc_gcsi <- format_df(fcrc_gcsi)
fcrc_gdsc <- format_df(fcrc_gdsc)
fcrc_ccle <- format_df(fcrc_ccle)


############################################################
# Filter low expressed transcripts
############################################################

print("Filtering low expressed transcripts")

# function to count and remove low expressed transcripts
filter_transcripts <- function(df) {
    counts <- colSums(df != 0)

    to_print <- data.frame(Category = c("Total Transcripts", ">5 CCLs", ">10 CCLs", ">15 CCLs"),
                           Number = c(length(counts), length(counts[counts > 5]), length(counts[counts > 10]), length(counts[counts > 15])))
    print(to_print)

    to_keep <- names(counts[counts > 5])
    df <- df[,which(colnames(df) %in% to_keep)]
    return(df)
}

ciri_gcsi <- filter_transcripts(ciri_gcsi)
ciri_gdsc <- filter_transcripts(ciri_gdsc)
ciri_ccle <- filter_transcripts(ciri_ccle)

circ_gcsi <- filter_transcripts(circ_gcsi)
circ_gdsc <- filter_transcripts(circ_gdsc)
circ_ccle <- filter_transcripts(circ_ccle)

cfnd_gcsi <- filter_transcripts(cfnd_gcsi)
cfnd_gdsc <- filter_transcripts(cfnd_gdsc)
cfnd_ccle <- filter_transcripts(cfnd_ccle)

fcrc_gcsi <- filter_transcripts(fcrc_gcsi)
fcrc_gdsc <- filter_transcripts(fcrc_gdsc)
fcrc_ccle <- filter_transcripts(fcrc_ccle)

############################################################
# Plot distribution of circRNA expression per pipeline
############################################################

print("Plotting expression distribution")

# function to vectorize each dataframe
vectorize <- function(df, pset, pipeline) {
    values <- as.vector(df[df != 0])
    df <- data.frame(Exp = values, 
                     PSet = rep(pset, length(values)),
                     Pipeline = rep(pipeline, length(values)))
    return(df)
}

toPlot <- rbind(vectorize(ciri_gcsi, "gCSI", "CIRI2"),
                vectorize(ciri_gdsc, "GDSC2", "CIRI2"),
                vectorize(ciri_ccle, "CCLE", "CIRI2"),
                vectorize(circ_gcsi, "gCSI", "CIRCexplorer2"),
                vectorize(circ_gdsc, "GDSC2", "CIRCexplorer2"),
                vectorize(circ_ccle, "CCLE", "CIRCexplorer2"),
                vectorize(cfnd_gcsi, "gCSI", "circRNA_finder"),
                vectorize(cfnd_gdsc, "GDSC2", "circRNA_finder"),
                vectorize(cfnd_ccle, "CCLE", "circRNA_finder"),
                vectorize(fcrc_gcsi, "gCSI", "find_circ"),
                vectorize(fcrc_gdsc, "GDSC2", "find_circ"),
                vectorize(fcrc_ccle, "CCLE", "find_circ"))

png("../results/figures/figure6/exp_distribution.png", width=150, height=75, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = PSet, y = Exp, fill = Pipeline)) + 
    geom_violin(trim = FALSE, scale = "width", adjust = 1.0) +
    theme_classic() + theme(legend.key.size = unit(0.5, 'cm'))
dev.off()

png("../results/figures/figure6/exp_distribution_log10.png", width=150, height=75, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = PSet, y = Exp, fill = Pipeline)) + 
    geom_violin(trim = FALSE, scale = "width", adjust = 1.0) +
    scale_y_log10() +
    theme_classic() + theme(legend.key.size = unit(0.5, 'cm'))
dev.off()


############################################################
# Robust normalization of circRNA expression
############################################################

print("Robust normalization")

# function to count and remove low expressed transcripts
robust_norm <- function(df) {
    for (circ in colnames(df)) {
        exp <- df[[circ]]
        non_zero <- exp != 0
        
        # median and IQR on non-zero values only
        med <- median(exp[non_zero])
        iqr <- IQR(exp[non_zero])
        
        # avoid division by zero in case IQR is 0
        if (iqr == 0) {
            df[non_zero, circ] <- 0  # or NA if preferred
            print(paste("IQR issue for"), circ)
        } else {
            df[non_zero, circ] <- (exp[non_zero] - med) / iqr
        }
    }
    return(df)
}


ciri_gcsi <- robust_norm(ciri_gcsi)
ciri_gdsc <- robust_norm(ciri_gdsc)
ciri_ccle <- robust_norm(ciri_ccle)

circ_gcsi <- robust_norm(circ_gcsi)
circ_gdsc <- robust_norm(circ_gdsc)
circ_ccle <- robust_norm(circ_ccle)

cfnd_gcsi <- robust_norm(cfnd_gcsi)
cfnd_gdsc <- robust_norm(cfnd_gdsc)
cfnd_ccle <- robust_norm(cfnd_ccle)

fcrc_gcsi <- robust_norm(fcrc_gcsi)
fcrc_gdsc <- robust_norm(fcrc_gdsc)
fcrc_ccle <- robust_norm(fcrc_ccle)


############################################################
# Plot distribution of circRNA expression per pipeline after normalization
############################################################

print("Plotting expression distribution")

# function to vectorize each dataframe
toPlot <- rbind(vectorize(ciri_gcsi, "gCSI", "CIRI2"),
                vectorize(ciri_gdsc, "GDSC2", "CIRI2"),
                vectorize(ciri_ccle, "CCLE", "CIRI2"),
                vectorize(circ_gcsi, "gCSI", "CIRCexplorer2"),
                vectorize(circ_gdsc, "GDSC2", "CIRCexplorer2"),
                vectorize(circ_ccle, "CCLE", "CIRCexplorer2"),
                vectorize(cfnd_gcsi, "gCSI", "circRNA_finder"),
                vectorize(cfnd_gdsc, "GDSC2", "circRNA_finder"),
                vectorize(cfnd_ccle, "CCLE", "circRNA_finder"),
                vectorize(fcrc_gcsi, "gCSI", "find_circ"),
                vectorize(fcrc_gdsc, "GDSC2", "find_circ"),
                vectorize(fcrc_ccle, "CCLE", "find_circ"))

png("../results/figures/figure6/exp_distribution_afternorm.png", width=150, height=75, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = PSet, y = Exp, fill = Pipeline)) + 
    geom_violin(trim = FALSE, scale = "width", adjust = 1.0) +
    theme_classic() + theme(legend.key.size = unit(0.5, 'cm'))
dev.off()

png("../results/figures/figure6/exp_distribution_log10_afternorm.png", width=150, height=75, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = PSet, y = Exp, fill = Pipeline)) + 
    geom_violin(trim = FALSE, scale = "width", adjust = 1.0) +
    scale_y_log10() +
    theme_classic() + theme(legend.key.size = unit(0.5, 'cm'))
dev.off()


############################################################
# Combine counts across pipelines
############################################################

print("Combining pipelines")

# function to combine exp counts for each sample across pipelines
combine_pipelines <- function(ciri_df, circ_df, cfnd_df, fcrc_df) {

    samples <- rownames(ciri_df)

    # initiate list to store results
    list_res <- vector("list", length(samples))

    # split by sample
    ciri_split <- split(ciri_df, rownames(ciri_df))
    circ_split <- split(circ_df, rownames(circ_df))
    cfnd_split <- split(cfnd_df, rownames(cfnd_df))
    fcrc_split <- split(fcrc_df, rownames(fcrc_df))

    for (i in seq_along(samples)) {

        print(i)
        sample <- samples[i]
        ciri_sample <- ciri_split[[sample]]
        circ_sample <- circ_split[[sample]]
        cfnd_sample <- cfnd_split[[sample]]
        fcrc_sample <- fcrc_split[[sample]]

        # remove transcripts detected in only one pipeline
        transcript_count <- c(colnames(ciri_sample), colnames(circ_sample), colnames(cfnd_sample), colnames(fcrc_sample)) |> table()
        keep <- transcript_count[transcript_count > 1] |> names()

        ciri_sample <- ciri_sample[,colnames(ciri_sample) %in% keep]
        circ_sample <- circ_sample[,colnames(circ_sample) %in% keep]
        cfnd_sample <- cfnd_sample[,colnames(cfnd_sample) %in% keep]
        fcrc_sample <- fcrc_sample[,colnames(fcrc_sample) %in% keep]

        merged <- rbindlist(list(ciri_sample, circ_sample, cfnd_sample, fcrc_sample), fill = TRUE) |> as.data.frame()
        merged[is.na(merged)] <- 0

        # keep only transcripts that have expression in at least 2 pipelines
        merged <- merged[,colnames(merged)[colSums(merged != 0) >= 2], drop = FALSE]

        # compute average of only non-zero pipelines
        list_res[[i]] <- sapply(merged, function(x) mean(x[x != 0], na.rm = TRUE))
    }

    # get sample names of samples with gene expression
    names(list_res) <- samples 
    list_res_dt <- lapply(list_res, function(x) as.data.table(as.list(x)))
    samples <- names(list_res_dt)[sapply(list_res_dt, function(x) !(is.null(x) || (is.data.table(x) && nrow(x) == 0)))]

    # bind
    pset_df <- rbindlist(list_res_dt, fill = TRUE) |> as.data.frame()
    pset_df[is.na(pset_df)] <- 0 
    rownames(pset_df) <- samples

    return(pset_df)
}

gcsi_df <- combine_pipelines(ciri_gcsi, circ_gcsi, cfnd_gcsi, fcrc_gcsi)
ccle_df <- combine_pipelines(ciri_ccle, circ_ccle, cfnd_ccle, fcrc_ccle)
gdsc_df <- combine_pipelines(ciri_gdsc, circ_gdsc, cfnd_gdsc, fcrc_gdsc)


############################################################
# Save dataframes
############################################################

print("Saving dataframes")

write.table(gcsi_df, file = paste0(df, "gcsi_counts.tsv"), quote = F, sep = "\t", col.names = T, row.names = T)
write.table(ccle_df, file = paste0(df, "ccle_counts.tsv"), quote = F, sep = "\t", col.names = T, row.names = T)
write.table(gdsc_df, file = paste0(df, "gdsc_counts.tsv"), quote = F, sep = "\t", col.names = T, row.names = T)

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
# Filter samples
############################################################

# load in cells to keep
load("../data/temp/all_intersect_cells.RData")

# function to filter circ expressiond dataframes 
filter_circ <- function(circ_counts) {
    circ_counts <- circ_counts[circ_counts$sample %in% to_keep,]
    return(circ_counts)
}

ciri_gcsi <- filter_circ(ciri_gcsi)
ciri_gdsc <- filter_circ(ciri_gdsc)
ciri_ccle <- filter_circ(ciri_ccle)

circ_gcsi <- filter_circ(circ_gcsi)
circ_gdsc <- filter_circ(circ_gdsc)
circ_ccle <- filter_circ(circ_ccle)

cfnd_gcsi <- filter_circ(cfnd_gcsi)
cfnd_gdsc <- filter_circ(cfnd_gdsc)
cfnd_ccle <- filter_circ(cfnd_ccle)

fcrc_gcsi <- filter_circ(fcrc_gcsi)
fcrc_gdsc <- filter_circ(fcrc_gdsc)
fcrc_ccle <- filter_circ(fcrc_ccle)


############################################################
# Format circRNA dataframes
############################################################

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
# Count distribution of 0 exp across cell lines
############################################################

# function to count number of 0 expression across each row
count_zero <- function(circ_counts) {
    df <- data.frame(Sample = rownames(circ_counts),
                     CountZero = rowSums(circ_counts == 0))
    df$Proportion <- df$CountZero / ncol(circ_counts) * 100
    df <- melt(df)
    return(df)
}

ciri_gcsi_cz <- count_zero(ciri_gcsi)
ciri_gdsc_cz <- count_zero(ciri_gdsc)
ciri_ccle_cz <- count_zero(ciri_ccle)

circ_gcsi_cz <- count_zero(circ_gcsi)
circ_gdsc_cz <- count_zero(circ_gdsc)
circ_ccle_cz <- count_zero(circ_ccle)

cfnd_gcsi_cz <- count_zero(cfnd_gcsi)
cfnd_gdsc_cz <- count_zero(cfnd_gdsc)
cfnd_ccle_cz <- count_zero(cfnd_ccle)

fcrc_gcsi_cz <- count_zero(fcrc_gcsi)
fcrc_gdsc_cz <- count_zero(fcrc_gdsc)
fcrc_ccle_cz <- count_zero(fcrc_ccle)



############################################################
# Format dataframes for plotting
############################################################

# function to combine df from each pset for per pipeline
merge_psets <- function(gcsi_df, ccle_df, gdsc_df, label) {

    gcsi_df$PSet <- "gCSI"
    ccle_df$PSet <- "CCLE"
    gdsc_df$PSet <- "GDSC"
    toPlot <- rbind(gcsi_df, ccle_df, gdsc_df)
    toPlot$PSet <- factor(toPlot$PSet, levels = c("gCSI", "CCLE", "GDSC"))
    toPlot$label <- label
    toPlot <- toPlot[toPlot$variable == "Proportion",]

    return(toPlot)

}

toPlot <- rbind(merge_psets(ciri_gcsi_cz, ciri_ccle_cz, ciri_gdsc_cz, "CIRI2"),
                merge_psets(circ_gcsi_cz, circ_ccle_cz, circ_gdsc_cz, "CIRCexplorer2"),
                merge_psets(cfnd_gcsi_cz, cfnd_ccle_cz, cfnd_gdsc_cz, "circRNA_finder"),
                merge_psets(fcrc_gcsi_cz, fcrc_ccle_cz, fcrc_gdsc_cz, "find_circ"))

toPlot$label <- factor(toPlot$label, levels = c("CIRI2", "CIRCexplorer2", "circRNA_finder", "find_circ"))


############################################################
# Plot distribution of zeros
############################################################

# set palette for plotting
pal = c("#839788", "#BFD7EA", "#BA9790", "#D5BC8A")

png(paste0("../results/figures/figure6/dist_zeros_", out, ".png"), width=150, height=75, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = value, y = label)) + 
    geom_violin(aes(fill = label, alpha = PSet), 
                    scale = "width", width = 1.1, 
                    position = position_dodge(width = 0.8)) + 
    theme_classic() + xlim(0, 100) +
    labs(x = "Proportion Zero Expression in Cell Lines", fill = "", alpha = "", y = "") +
    scale_fill_manual(values = pal) + 
    scale_alpha_manual(values = c(0.2, 0.5, 0.9)) +
    scale_y_discrete(limits = c("find_circ", "circRNA_finder", "CIRCexplorer2", "CIRI2")) + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3), legend.key.size = unit(0.5, 'cm')) 
dev.off()


############################################################
# Combine counts across pipelines
############################################################

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

        sample <- samples[i]
        ciri_sample <- ciri_split[[sample]]
        circ_sample <- circ_split[[sample]]
        cfnd_sample <- cfnd_split[[sample]]
        fcrc_sample <- fcrc_split[[sample]]

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


#keep = intersect(intersect(intersect(colnames(ciri_gcsi), colnames(circ_gcsi)),colnames(cfnd_gcsi)), colnames(fcrc_gcsi))
#keep = keep[1:10]
#ciri_df = ciri_gcsi[,colnames(ciri_gcsi) %in% keep]
#circ_df = circ_gcsi[,colnames(circ_gcsi) %in% keep]
#cfnd_df = cfnd_gcsi[,colnames(cfnd_gcsi) %in% keep]
#fcrc_df = fcrc_gcsi[,colnames(fcrc_gcsi) %in% c(keep, "A1BG", "A1BG-AS1", "A1CF", "A2M")]

gcsi_df <- combine_pipelines(ciri_gcsi, circ_gcsi, cfnd_gcsi, fcrc_gcsi)
ccle_df <- combine_pipelines(ciri_ccle, circ_ccle, cfnd_ccle, fcrc_ccle)
gdsc_df <- combine_pipelines(ciri_gdsc, circ_gdsc, cfnd_gdsc, fcrc_gdsc)


############################################################
# Save dataframes
############################################################

write.table(gcsi_df, file = paste0(df, "gcsi_counts.tsv"), quote = F, sep = "\t", col.names = T, row.names = T)
write.table(ccle_df, file = paste0(df, "ccle_counts.tsv"), quote = F, sep = "\t", col.names = T, row.names = T)
write.table(gdsc_df, file = paste0(df, "gdsc_counts.tsv"), quote = F, sep = "\t", col.names = T, row.names = T)

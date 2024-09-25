suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(ggplotify))
suppressMessages(library(matrixStats))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(RColorBrewer))


# ========== Plot distribution of cancer types across common cell lines for each PSet pair ========== #


# load in metadata
load("../results/data/tissue-metadata.RData")

# figure out tissue overlap
gcsi_ccle <- data.frame(sample = intersect(gcsi$cellid, ccle$cellid),
                        gcsi_tissue = gcsi$tissue[match(intersect(gcsi$cellid, ccle$cellid), gcsi$cellid)],
                        ccle_tissue = ccle$tissue[match(intersect(gcsi$cellid, ccle$cellid), ccle$cellid)])
gcsi_gdsc <- data.frame(sample = intersect(gcsi$cellid, gdsc$cellid),
                        gcsi_tissue = gcsi$tissue[match(intersect(gcsi$cellid, gdsc$cellid), gcsi$cellid)],
                        gdsc_tissue = gdsc$tissue[match(intersect(gcsi$cellid, gdsc$cellid), gdsc$cellid)])
ccle_gdsc <- data.frame(sample = intersect(ccle$cellid, gdsc$cellid),
                        ccle_tissue = ccle$tissue[match(intersect(ccle$cellid, gdsc$cellid), ccle$cellid)],
                        gdsc_tissue = gdsc$tissue[match(intersect(ccle$cellid, gdsc$cellid), gdsc$cellid)])


# function to plot mismatches
plot_mismatches <- function(mismatch, label1, label2) {

    colnames(mismatch) <- c("sample", "pset1", "pset2")

    mismatch$pair <- paste0(mismatch$pset1, "-", mismatch$pset2)
    mismatch <- mismatch[order(mismatch$pair),]
    mismatch$rank <- c(1:nrow(mismatch))

    p1 <- ggplot(mismatch, aes(x = rank, y = 1, fill = pset1)) +
        geom_tile(color = "white") +  scale_fill_brewer(palette = "Paired") + 
        labs(fill = paste0(label1, " Tissue")) + theme_void() 
    p2 <- ggplot(mismatch, aes(x = rank, y = 1, fill = pset2)) +
        geom_tile(color = "white") +  labs(fill = paste0(label2, " Tissue")) + theme_void() 

    l1 <- as.ggplot(get_legend(p1))
    l2 <- as.ggplot(get_legend(p2))

    p1 <- p1 + guides(fill="none")
    p2 <- p2 + guides(fill="none")

    png("temp.png", width=150, height=50, units='mm', res = 600, pointsize=80)
    ggarrange(p1, p2, nrow = 2)
    dev.off()
    
}

# plot mismatches in tissue groupings
#mismatch <- gcsi_ccle[-which(gcsi_ccle$gcsi_tissue == gcsi_ccle$ccle_tissue),]
#mismatch <- gcsi_gdsc[-which(gcsi_gdsc$gcsi_tissue == gcsi_gdsc$gdsc_tissue),]
#mismatch <- ccle_gdsc[-which(ccle_gdsc$gdsc_tissue == ccle_gdsc$ccle_tissue),]

# TODO: properly call function and output plots

# plot the matching
#match <- gcsi_ccle[which(gcsi_ccle$gcsi_tissue == gcsi_ccle$ccle_tissue),]
#match <- gcsi_gdsc[which(gcsi_gdsc$gcsi_tissue == gcsi_gdsc$gdsc_tissue),]
#match <- ccle_gdsc[which(ccle_gdsc$gdsc_tissue == ccle_gdsc$ccle_tissue),]

# TODO: turn the following into a function:
to_create <- function(match) {
    colnames(match) <- c("sample", "pset1", "pset2")

    match_count <- match %>% group_by(pset1) %>% summarise(Count = n()) 

    ggplot(match_count, aes(x = "", y = Count, fill = pset1)) +
    geom_bar(stat="identity", width = 1, color = "black") +
    geom_text(aes(label = Count), position = position_stack(vjust = 0.5)) +
    #scale_fill_manual(values = pal_primary_count) +
    coord_polar("y", start = 0) + theme_void() + labs(fill = "Tissue")

}


# create consensus for each pair
gcsi_ccle$tissue <- gcsi_ccle$gcsi_tissue
gcsi_gdsc$tissue <- gcsi_gdsc$gcsi_tissue
ccle_gdsc$tissue <- ccle_gdsc$ccle_tissue
ccle_gdsc$tissue[ccle_gdsc$tissue %in% c("Autonomic Ganglia", "Central Nervous System")] <- "Brain"
ccle_gdsc$tissue[ccle_gdsc$tissue %in% c("Oesophagus")] <- "Head-Neck"

# keep only tissue groups with at least 4 cell lines
keep <- names(table(gcsi_ccle$tissue)[table(gcsi_ccle$tissue) > 4])
gcsi_ccle <- gcsi_ccle[gcsi_ccle$tissue %in% keep,]

keep <- names(table(gcsi_gdsc$tissue)[table(gcsi_gdsc$tissue) > 4])
gcsi_gdsc <- gcsi_gdsc[gcsi_gdsc$tissue %in% keep,]

keep <- names(table(ccle_gdsc$tissue)[table(ccle_gdsc$tissue) > 4])
ccle_gdsc <- ccle_gdsc[ccle_gdsc$tissue %in% keep,]


# ========== Compute pairwise spearman correlations from circRNA datasets ========== #

# load circRNA expression data
ciri_gcsi <- fread("../data/processed_cellline/all_samples/CIRI2/ciri_gcsi_counts.tsv", data.table = F)
ciri_gdsc <- fread("../data/processed_cellline/all_samples/CIRI2/ciri_gdsc_counts.tsv", data.table = F)
ciri_ccle <- fread("../data/processed_cellline/all_samples/CIRI2/ciri_ccle_counts.tsv", data.table = F)

circ_gcsi <- fread("../data/processed_cellline/all_samples/CIRCexplorer2/circ_gcsi_counts.tsv", data.table = F)
circ_gdsc <- fread("../data/processed_cellline/all_samples/CIRCexplorer2/circ_gdsc_counts.tsv", data.table = F)
circ_ccle <- fread("../data/processed_cellline/all_samples/CIRCexplorer2/circ_ccle_counts.tsv", data.table = F)

cfnd_gcsi <- fread("../data/processed_cellline/all_samples/circRNA_finder/cfnd_gcsi_counts.tsv", data.table = F)
cfnd_gdsc <- fread("../data/processed_cellline/all_samples/circRNA_finder/cfnd_gdsc_counts.tsv", data.table = F)
cfnd_ccle <- fread("../data/processed_cellline/all_samples/circRNA_finder/cfnd_ccle_counts.tsv", data.table = F)

fcrc_gcsi <- fread("../data/processed_cellline/all_samples/find_circ/fcrc_gcsi_counts.tsv", data.table = F)
fcrc_gdsc <- fread("../data/processed_cellline/all_samples/find_circ/fcrc_gdsc_counts.tsv", data.table = F)
fcrc_ccle <- fread("../data/processed_cellline/all_samples/find_circ/fcrc_ccle_counts.tsv", data.table = F)


# function to compute tissue-specific pairwise spearman correlations 
compute_spearman <- function(df1, df2, pset_label, pipeline_label, meta) {

    # INPUTS:
    #   df1: first circRNA dataframe
    #   df2: second circRNA dataframe
    #   pset_label: options are "gCSI/CCLE", "gCSI/GDSC2", or "GDSC2/CCLE"
    #   pipeline_label: options are "CIRI2", "CIRCexplorer2", "circRNA_finder", or "find_circ"
    #   meta: metadata with tissue column

    # initialize dataframe to store results
    correlations <- data.frame(matrix(nrow=0, ncol=2))

    print(paste("--------Starting", pset_label, "-", pipeline_label))
    print(unique(meta$tissue))


    # loop through each unique tissue
    for (tissue in unique(meta$tissue)) {

        print(paste("--------Tissue:", tissue))

        # initialize vector to store results
        tissue_res <- c()

        # keep only cell lines in tissue of interest
        cells_to_keep <- meta$sample[meta$tissue == tissue]
        print(paste("Cells:", length(cells_to_keep)))
        df1_sub <- df1[df1$sample %in% cells_to_keep,]
        df2_sub <- df2[df2$sample %in% cells_to_keep,]

        # filter circRNA transcripts with low detection rates
        min_num_cells = length(cells_to_keep) - 2
        df1_sub <- df1_sub[,-which(colnames(df1_sub) %in% names(which(colSums(df1_sub == 0) > min_num_cells)))]
        df2_sub <- df2_sub[,-which(colnames(df2_sub) %in% names(which(colSums(df2_sub == 0) > min_num_cells)))]

        # keep only common circRNA transcripts
        common_transcripts <- intersect(colnames(df1_sub), colnames(df2_sub))
        print(paste("Transcripts:", length(common_transcripts)))

        # only compute spearman if there is at least 3 transcripts
        if (length(common_transcripts) > 2) {

            df1_sub <- df1_sub[,colnames(df1_sub) %in% common_transcripts]
            df2_sub <- df2_sub[,colnames(df2_sub) %in% common_transcripts]

            # order dataframe
            df1_sub <- df1_sub[order(df1_sub$sample),order(colnames(df1_sub))]
            df2_sub <- df2_sub[order(df2_sub$sample),order(colnames(df2_sub))]

            # remove sample names
            rownames(df1_sub) <- df1_sub$sample
            rownames(df2_sub) <- df2_sub$sample
            df1_sub$sample <- NULL
            df2_sub$sample <- NULL

            # loop through each common transcript
            for (i in 1:ncol(df1_sub)) {

                # compute correlations of transcript expression for pairs of psets
                spearman_results <- suppressWarnings(cor(x = as.numeric(df1_sub[, i]), y = as.numeric(df2_sub[, i]), method = "spearman"))
                
                # combine results
                tissue_res <- c(tissue_res, spearman_results)
            }
            
            tissue_res <- data.frame(Spearman = tissue_res, Tissue = tissue)
            correlations <- rbind(correlations, tissue_res)

        } else { print(paste("-- SKIPPING", tissue))}

    }

    if (nrow(correlations) > 0) {
        correlations$PSet_Pair <- pset_label
        correlations$Pipeline <- pipeline_label
    }

    return(correlations)
}


# compute spearman correlations
gcsi_ccle_ciri <- compute_spearman(ciri_gcsi, ciri_ccle, "gCSI/CCLE", "CIRI2", gcsi_ccle)
gcsi_ccle_circ <- compute_spearman(circ_gcsi, circ_ccle, "gCSI/CCLE", "CIRCexplorer2", gcsi_ccle)
gcsi_ccle_cfnd <- compute_spearman(cfnd_gcsi, cfnd_ccle, "gCSI/CCLE", "circRNA_finder", gcsi_ccle)
gcsi_ccle_fcrc <- compute_spearman(fcrc_gcsi, fcrc_ccle, "gCSI/CCLE", "find_circ", gcsi_ccle)

gcsi_gdsc_ciri <- compute_spearman(ciri_gcsi, ciri_gdsc, "gCSI/GDSC2", "CIRI2", gcsi_gdsc)
gcsi_gdsc_circ <- compute_spearman(circ_gcsi, circ_gdsc, "gCSI/GDSC2", "CIRCexplorer2", gcsi_gdsc)
gcsi_gdsc_cfnd <- compute_spearman(cfnd_gcsi, cfnd_gdsc, "gCSI/GDSC2", "circRNA_finder", gcsi_gdsc)
gcsi_gdsc_fcrc <- compute_spearman(fcrc_gcsi, fcrc_gdsc, "gCSI/GDSC2", "find_circ", gcsi_gdsc)

ccle_gdsc_ciri <- compute_spearman(ciri_ccle, ciri_gdsc, "CCLE/GDSC2", "CIRI2", ccle_gdsc)
ccle_gdsc_circ <- compute_spearman(circ_ccle, circ_gdsc, "CCLE/GDSC2", "CIRCexplorer2", ccle_gdsc)
ccle_gdsc_cfnd <- compute_spearman(cfnd_ccle, cfnd_gdsc, "CCLE/GDSC2", "circRNA_finder", ccle_gdsc)
ccle_gdsc_fcrc <- compute_spearman(fcrc_ccle, fcrc_gdsc, "CCLE/GDSC2", "find_circ", ccle_gdsc)

# merge results for each pipeline
ciri_stability <- rbind(gcsi_ccle_ciri, gcsi_gdsc_ciri, ccle_gdsc_ciri)
circ_stability <- rbind(gcsi_ccle_circ, gcsi_gdsc_circ, ccle_gdsc_circ)
cfnd_stability <- rbind(gcsi_ccle_cfnd, gcsi_gdsc_cfnd, ccle_gdsc_cfnd)
fcrc_stability <- rbind(gcsi_ccle_fcrc, gcsi_gdsc_fcrc, ccle_gdsc_fcrc)

#save(ciri_stability, circ_stability, cfnd_stability, fcrc_stability, file = "../results/data/si_distribution_tissue.RData")

#######################
### Plot Stability ####
#######################

# format dataframe for plotting
toPlot <- rbind(ciri_stability, circ_stability, cfnd_stability, fcrc_stability)
toPlot$label <- paste(toPlot$Pipeline, "-", toPlot$PSet_Pair)

for (group in unique(toPlot$label)) {

    # subset dataframe to group
    subset_df <- toPlot[toPlot$label == group,]

    # get number of transcripts per tissue type to include in plot
    counts <- as.data.frame(subset_df %>% group_by(Tissue) %>% summarise(Count = n()))
    subset_df$Count <- counts[match(subset_df$Tissue, counts$Tissue),]$Count
    subset_df$Tissue <- paste0(subset_df$Tissue, "\n", subset_df$Count)

    # get parameters for plotting
    filename <- paste0("../results/figures/figure2/tissue/", gsub(" - ", "_", gsub("/", "-", group)), ".png")
    w <- ifelse(subset_df$PSet_Pair[1] == "gCSI/CCLE", 300, 100)

    png(filename, width=w, height=100, units='mm', res = 600, pointsize=80)
    print(ggplot(subset_df, aes(x = Tissue, y = Spearman)) + ylim(-1, 1) +
        geom_violin(alpha = 0.8, fill = "#5B4B49") + geom_boxplot(width=0.1, alpha = 0.3) +
        theme_classic() + labs(x = "", y = "Stability Index") +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
            plot.title = element_text(hjust = 0.5)) + ggtitle(subset_df$label[1]) +
        geom_hline(yintercept = 0, linetype = "dotted"))
    dev.off()

}









# ========== Plot distribution of cancer types across 48 common cell lines ========== #

# TODO: fix code below

# mapping cell line names across metadata and analysis
mapping <- c("DU 145" = "DU145", "HCT-15" = "HCT 15", "HT-115" = "HT115", 
             "KARPAS-620" = "Karpas-620", "KU812" = "Ku812", "LS 180" = "LS180",
             "SJCRH30" = "Rh30", "SW 1116" = "SW1116", "SW 1463" = "SW1463", 
             "SW 1573" = "SW1573", "SW 403" = "SW403", "SW 48" = "SW48", 
             "SW 620" = "SW620", "SW 780" = "SW780", "SW 837" = "SW837", "SW 948" = "SW948")

for (i in 1:length(metadata$Cell_line)) {
    cell = metadata$Cell_line[i]
    if (cell %in% names(mapping)) {metadata$Cell_line[i] <- unname(mapping[cell])}
}

# keep only cell lines and variables of interest
metadata <- as.data.frame(metadata[metadata$Cell_line %in% intersected_rnacells,colnames(metadata) %in% c("Cell_line", "Tissue_supergroup", "Primary_Tissue")])
write.csv(metadata, file = "intersected_cells_metadata.csv", quote = F, col.names = T, row.names = F)



metadata <- read.csv("intersected_cells_metadata.csv")


# count number of cell lines in each group of Primary_Tissue
primary_tissue_count <- metadata %>% group_by(Tissue_supergroup, Primary_Tissue) %>% summarise(Count = n()) 
primary_tissue_count <- as.data.frame(primary_tissue_count)
primary_tissue_count$Primary_Tissue <- factor(primary_tissue_count$Primary_Tissue, levels = primary_tissue_count$Primary_Tissue)

# count number of cell lines in each group of Tissue_supergroup
tissue_count <- metadata %>% group_by(Tissue_supergroup) %>% summarise(Count = n()) 
tissue_count$hsize <- 4

# set colour palettes for plotting
pal_tissue_count <- brewer.pal(n = 12, name = "Set3")
pal_primary_count <- c(pal_tissue_count[1:3], "#DD5142", "#FFA89E", 
                       pal_tissue_count[5], "#DD9441", "#FFC98C", "#90BF40", "#D5FB91", 
                       pal_tissue_count[8:12])

# save plots
p1 <- ggplot(primary_tissue_count, aes(x = "", y = Count, fill = Primary_Tissue)) +
  geom_bar(stat="identity", width = 1, color = "black") +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = pal_primary_count) +
  coord_polar("y", start = 0) + theme_void() + labs(fill = "Primary Tissue")

p2 <- ggplot(tissue_count, aes(x = hsize, y = Count, fill = Tissue_supergroup)) +
  geom_bar(stat="identity", width = 1, color = "black") + 
  coord_polar("y", start = 0) + theme_void() +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = pal_tissue_count) + xlim(c(0.2, 4.5)) + labs(fill = "Tissue Supergroup")

ggsave(p1, filename = "primary.png", width = 5.5, height = 4.5, units = "in", bg = "transparent")
ggsave(p2, filename = "tissue.png", width = 7, height = 6, units = "in", bg = "transparent")



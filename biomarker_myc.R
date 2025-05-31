# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ComplexHeatmap)
    library(ggh4x)
    library(ggVennDiagram)
    library(stringr)
})


############################################################
# Load in circRNA gene data
############################################################

# save dataframes
path = "../data/processed_cellline/GE_all_samples/"

ciri_gcsi_ge <- fread(paste0(path, "CIRI2/ciri_gcsi_counts.tsv"), data.table = F)
ciri_gdsc_ge <- fread(paste0(path, "CIRI2/ciri_gdsc_counts.tsv"), data.table = F)
ciri_ccle_ge <- fread(paste0(path, "CIRI2/ciri_ccle_counts.tsv"), data.table = F)

circ_gcsi_ge <- fread(paste0(path, "CIRCexplorer2/circ_gcsi_counts.tsv"), data.table = F)
circ_gdsc_ge <- fread(paste0(path, "CIRCexplorer2/circ_gdsc_counts.tsv"), data.table = F)
circ_ccle_ge <- fread(paste0(path, "CIRCexplorer2/circ_ccle_counts.tsv"), data.table = F)

cfnd_gcsi_ge <- fread(paste0(path, "circRNA_finder/cfnd_gcsi_counts.tsv"), data.table = F)
cfnd_gdsc_ge <- fread(paste0(path, "circRNA_finder/cfnd_gdsc_counts.tsv"), data.table = F)
cfnd_ccle_ge <- fread(paste0(path, "circRNA_finder/cfnd_ccle_counts.tsv"), data.table = F)

fcrc_gcsi_ge <- fread(paste0(path, "find_circ/fcrc_gcsi_counts.tsv"), data.table = F)
fcrc_gdsc_ge <- fread(paste0(path, "find_circ/fcrc_gdsc_counts.tsv"), data.table = F)
fcrc_ccle_ge <- fread(paste0(path, "find_circ/fcrc_ccle_counts.tsv"), data.table = F)


############################################################
# Load in drug response data and subset
############################################################

# load in drug sensitivity from PSets
load("../data/temp/sensitivity_data.RData")

# belinostat (CTRP)
belinostat <- ctrp_sen[rownames(ctrp_sen) == "Belinostat",,drop=F]
belinostat <- belinostat[,colSums(is.na(belinostat))<nrow(belinostat)]


############################################################
# Keep only circRNAs mapped to MYC gene
############################################################

keep_myc <- function(counts) {
    if ("MYC" %in% colnames(counts)) {
        counts <- counts[,colnames(counts) %in% c("sample", "MYC")]
        counts[counts == 0] <- NA
        rownames(counts) <- counts$sample
        counts$sample <- NULL
        return(counts)
    }
    return(0)
}

ciri_gcsi_ge <- keep_myc(ciri_gcsi_ge)
ciri_gdsc_ge <- keep_myc(ciri_gdsc_ge)
ciri_ccle_ge <- keep_myc(ciri_ccle_ge)

circ_gcsi_ge <- keep_myc(circ_gcsi_ge)
circ_gdsc_ge <- keep_myc(circ_gdsc_ge)
circ_ccle_ge <- keep_myc(circ_ccle_ge)

cfnd_gcsi_ge <- keep_myc(cfnd_gcsi_ge)
cfnd_gdsc_ge <- keep_myc(cfnd_gdsc_ge)
cfnd_ccle_ge <- keep_myc(cfnd_ccle_ge)

fcrc_gcsi_ge <- keep_myc(fcrc_gcsi_ge)
fcrc_gdsc_ge <- keep_myc(fcrc_gdsc_ge)
fcrc_ccle_ge <- keep_myc(fcrc_ccle_ge)


############################################################
# Check number of 0 expression
############################################################

nrow(ciri_ccle_ge[ciri_ccle_ge$MYC == 0,]) #504/524
nrow(circ_ccle_ge[circ_ccle_ge$MYC == 0,]) #519/524
nrow(cfnd_ccle_ge[cfnd_ccle_ge$MYC == 0,]) #489/524
nrow(fcrc_ccle_ge[fcrc_ccle_ge$MYC == 0,]) #69/524


############################################################
# Binarize transcript expression by median
############################################################

# function to compute drug response associations
binary_dr <- function(counts_df, drug_df, pset) {

    # create data frame to hold results
    combinations <- expand.grid(Drug = rownames(drug_df), Feature = colnames(counts_df))
    combinations$num_samples <- combinations$pval <-  combinations$W <- NA

    # compute median
    median <- median(counts_df$MYC, na.rm = T)
    
    # initiate row count
    row = 1
        
    # binarize transcript expression by median
    subset <- counts_df[complete.cases(counts_df),,drop=F]

    high_exp <- rownames(subset)[as.vector(subset >= median)]
    low_exp <- rownames(subset)[-which(rownames(subset) %in% high_exp)]

    # loop through each drug

    for (j in seq_along(rownames(drug_df))) {
        
        high = drug_df[j,colnames(drug_df) %in% high_exp] |> as.numeric()
        low = drug_df[j,colnames(drug_df) %in% low_exp] |> as.numeric()

        if (length(high[!is.na(high)]) > 1 & length(low[!is.na(low)]) > 1) {

            # wilcoxon rank sum test
            res <- wilcox.test(high, low, alternative = "two.sided", exact = FALSE)

            # save results
            combinations$W[row] <- res$statistic
            combinations$pval[row] <- res$p.value
            combinations$num_samples[row] <- nrow(subset)

        } else {
            # save results
            combinations$W[row] <- NA
            combinations$pval[row] <- NA
            combinations$num_samples[row] <- nrow(subset)
        }

        row <- row + 1
        
    }
    combinations$pset <- pset
    combinations$FDR <- p.adjust(combinations$pval, method = "BH")

    # format dataframe for plotting 
    combinations <- combinations[order(combinations$W, decreasing = T),]
    combinations$rank <- seq_len(nrow(combinations))
    combinations$pair <- paste(combinations$Drug, combinations$Feature, sep = "_")
    return(combinations)
}

belinostat_res <- rbind(binary_dr(ciri_ccle_ge, belinostat, "CIRI2"),
                        binary_dr(circ_ccle_ge, belinostat, "CIRCexplorer2"),
                        binary_dr(cfnd_ccle_ge, belinostat, "circRNA_finder"),
                        binary_dr(fcrc_ccle_ge, belinostat, "find_circ"))


############################################################
# Plot binary drug response association
############################################################

# function to create dataframe for binary drug response plot
binarize_myc <- function(counts_df, drug_df) {

    # compute median
    median <- median(counts_df$MYC, na.rm = T)
    
    # binarize transcript expression by median
    subset <- counts_df[complete.cases(counts_df),,drop=F]
    high_exp <- rownames(subset)[as.vector(subset >= median)]
    low_exp <- rownames(subset)[-which(rownames(subset) %in% high_exp)]

    # create dataframe
    toPlot <- rbind(t(drug_df[j,colnames(drug_df) %in% high_exp]), t(drug_df[j,colnames(drug_df) %in% low_exp])) |> as.data.frame()
    toPlot$label <- ifelse(rownames(toPlot) %in% high_exp, "high", "low")
    rownames(toPlot) <- NULL

    return(toPlot)
}

# function to plot box plot
plot_box <- function(toPlot, drugname, filename) {
    colnames(toPlot)[1] <- "Drug"
    png(paste0(filename, ".png"), width = 4, height = 4, res = 600, units = "in")
    print({
        ggplot(toPlot, aes(x = label, y = Drug, fill = label)) + geom_boxplot() +
            geom_jitter(shape=16, position=position_jitter(0.2)) +
            scale_fill_manual("MYC Expression", values = c("#B23A48", "#3E92CC")) + theme_classic() +
            labs(x = "MYC Expression", y = paste(drugname, "Response (AAC)"))
    })
    dev.off()
}

plot_box(binarize_myc(ciri_ccle_ge, belinostat), "Belinostat", "CIRI_CCLE_Belinostat")
plot_box(binarize_myc(circ_ccle_ge, belinostat), "Belinostat", "CIRC_CCLE_Belinostat")
plot_box(binarize_myc(cfnd_ccle_ge, belinostat), "Belinostat", "CFND_CCLE_Belinostat")
plot_box(binarize_myc(fcrc_ccle_ge, belinostat), "Belinostat", "FCRC_CCLE_Belinostat")

belinostat_plot <- rbind(binarize_myc(ciri_ccle_ge, belinostat),
                         binarize_myc(circ_ccle_ge, belinostat),
                         binarize_myc(cfnd_ccle_ge, belinostat),
                         binarize_myc(fcrc_ccle_ge, belinostat))

plot_box(belinostat_plot, "Belinostat", "CCLE_Belinostat")


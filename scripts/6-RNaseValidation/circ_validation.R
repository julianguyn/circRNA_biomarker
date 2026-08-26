# load packages
suppressPackageStartupMessages({
    library(ggplot2)
    library(reshape2)
    library(ggpubr)
    library(ggsignif)
    library(patchwork)
})

options(stringsAsFactors = F)

source("scripts/utils/circ_validation_processing.R")
source("scripts/utils/palettes.R")

############################################################
# Summarize circRNA counts
############################################################

# CIRI2
ciri_ribo0 <- summarizeCIRI(dir_path = "data/processed_cellline/CIRI2/hansen/result") #both Ribo-zero + RNAse-R
ciri_polyA <- summarizeCIRI(dir_path = "data/processed_cellline/CIRI2/hansen_match/result") #poly-A equivalent to Hansen data

# CIRCexplorer2
circ_ribo0 <- summarizeCIRCexplorer(dir_path = "data/processed_cellline/CIRCexplorer2/hansen/annotate") #both Ribo-zero + RNAse-R
circ_polyA <- summarizeCIRCexplorer(dir_path = "data/processed_cellline/CIRCexplorer2/hansen_match/annotate") #poly-A equivalent to Hansen data

############################################################
# Summarize RNaseR validated circRNAs
############################################################

# using matched Hansen RNAse-R enriched samples

# CIRI - 22Rv1 (CCLE), LNCaP(CCLE), PC3(CCLE/gCSI)
ciri_polyA_val <- filterCIRI(
  nonRNAseR_dir = "data/processed_cellline/CIRI2/hansen_match/result", 
  RNAseR_dir = "data/processed_cellline/CIRI2/hansen/result/",
  suff = "\\.tsv$"
)

# CIRI - 22Rv1, LNCaP, and PC3 Hansen Ribo-Zero samples
ciri_ribo0_val <- filterCIRI(
  nonRNAseR_dir = "data/processed_cellline/CIRI2/hansen/result", 
  RNAseR_dir = "data/processed_cellline/CIRI2/hansen/result/",
  suff = "*neg.tsv$"
)

# CIRCexplorer2 - 22Rv1 (CCLE), LNCaP(CCLE), PC3(CCLE/gCSI)
circ_polyA_val <- filterCIRC(
  nonRNAseR_dir = "data/processed_cellline/CIRCexplorer2/hansen_match/annotate", 
  RNAseR_dir = "data/processed_cellline/CIRCexplorer2/hansen/annotate",
  suff = "\\.txt$"
)

# CIRCexplorer2 - 22Rv1, LNCaP, and PC3 Hansen Ribo-Zero samples
circ_ribo0_val <- filterCIRC(
  nonRNAseR_dir = "data/processed_cellline/CIRCexplorer2/hansen/annotate", 
  RNAseR_dir = "data/processed_cellline/CIRCexplorer2/hansen/annotate",
  suff = "*neg*"
)

############################################################
# Compile counts for plotting
############################################################

compile_counts <- function(ribo0, polyA, ribo0_val, polyA_val) {
    ribo0 <- ribo0[which(ribo0$sample %in% c("22Rv1neg", "LNCaPneg", "PC3neg")),]
    polyA <- polyA[which(polyA$sample %in% c("22Rv1ccle", "LNCaPccle", "PC3ccle")),]
    polyA_val <- polyA_val[which(polyA_val$sample %in% c("22Rv1ccle", "LNCaPccle", "PC3ccle")),]

    compiled <- rbind(ribo0, polyA, ribo0_val, polyA_val)
    compiled$count <- as.numeric(compiled$count)

    compiled$sample <- c(rep(c("22Rv1", "LNCaP", "PC3"), 4))
    compiled$method <- c(rep(c(rep("ribo0", 3), rep("polyA", 3)), 2))
    compiled$validation <- c(rep("RNase R-", 6), rep("RNase R+", 6))
    
    return(compiled)
}

ciri <- compile_counts(ciri_ribo0, ciri_polyA, ciri_ribo0_val, ciri_polyA_val)
circ <- compile_counts(circ_ribo0, circ_polyA, circ_ribo0_val, circ_polyA_val)

############################################################
# Compute fold changes
############################################################

# fold difference calculations
# (new value / original value) - 1

calculate_fc <- function(df) {
    for (group in unique(df$validation)) {
        message(paste("FC for", group))
        for (cell in unique(df$sample)) {
            ribo0 <- df$count[which(df$validation == group & df$sample == cell & df$method == "ribo0")]
            polyA <- df$count[which(df$validation == group & df$sample == cell & df$method == "polyA")]
            fc <- ribo0 / polyA - 1
            message(paste("---- Cell:", cell, "- FC:", fc))
        }
    }
}
calculate_fc(ciri)
calculate_fc(circ)

############################################################
# Plot counts
############################################################

plot_combined_counts <- function(df, pipeline) {
    p <- ggplot(df, aes(x=sample, y=log2(count), fill=method)) + 
    geom_bar(stat="identity", width = 0.7, position = position_dodge(), color = "black") + 
    facet_grid(~validation) + 
    scale_y_continuous(limits = c(0, 14), expand=c(0,0)) +
    scale_fill_manual("Method", values=protocol_pal, limits=c("polyA", "ribo0"), labels=c("poly(A)-selected", "rRNA-depleted")) + 
    theme_bw() +
    theme(strip.background = element_rect(fill = "white", color = "black"),
            legend.key.size = unit(0.5, 'cm'),
            plot.title = element_text(hjust = 0.5), 
            axis.text.x = element_text(vjust = 1, hjust = 1, angle = 45)) +
    labs(x = "Cell line \n", y = "Log2 Normalized Counts", title = pipeline)
    return(p)
}

p1 <- plot_combined_counts(ciri, "CIRI2")
p2 <- plot_combined_counts(circ, "CIRCexplorer2")

p <- p1 + p2 + plot_layout(guides = "collect")
filename <- "results/figures/suppfig5/rnaseR_validation.png"
ggsave(filename, p, width = 7, height = 3)


############################################################
# Old code
############################################################

# annotations for geom_signif
#ciri_annotations <- data.frame(validation = c("RNase R-", "RNase R+"), method = c("poly-A", "RiboMinus"),
#                                      y = c(5297, 4291, 11079, 9700, 10913, 9542), start = c(0.85, 0.85, 1.85, 1.85, 2.85, 2.85),
#                                      end = c(1.15, 1.15, 2.15, 2.15, 3.15, 3.15), labels = c("", "", "", "", "", ""))
#ciri_labels <- data.frame(label = c("", "****", "*", "***", "*", "****"),
#                                 cyl = c("RNase R-", "RNase R+"), validation = c("RNase R-", "RNase R+"), method = c("poly-A", "RiboMinus"))

#circ_annotations <- data.frame(validation = c("RNase R-", "RNase R+"), method = c("poly-A", "RiboMinus"),
#                                      y = c(9534, 5838, 13158, 10423, 14522, 11011), start = c(0.85, 0.85, 1.85, 1.85, 2.85, 2.85),
#                                      end = c(1.15, 1.15, 2.15, 2.15, 3.15, 3.15), labels = c("", "", "", "", "", ""))
#circ_labels <- data.frame(label = c("", "**", "", "***", "", "**"),
#                                 cyl = c("RNase R-", "RNase R+"), validation = c("RNase R-", "RNase R+"), method = c("poly-A", "RiboMinus"))


# plot circRNA counts
p1 <- ggplot(ciri, aes(x=sample, y=count, fill=method)) + 
  geom_bar(stat="identity", size = 0.5, width = 0.7, position = position_dodge(), color = "black") + 
  geom_signif(data = ciri_annotations,
              aes(y_position = y, 
                  xmin = start,
                  xmax = end, 
                  annotations = labels),
              tip_length = c(0.05, 0.01),
              manual = TRUE) +
  geom_text(data = ciri_labels,
            mapping = aes(x = c(0.9, 0.75, 1.9, 1.8, 2.9, 2.75),
                          y = c(5297, 4291, 11079, 9700, 10913, 9542),
                          label = label),
            hjust = -0.1, vjust = -1) +
  facet_grid(~validation) + 
  scale_y_continuous(limits = c(0, 15500), expand=c(0,0)) +
  scale_fill_manual("Method",values=c("#8B7B96", "#71A2B6"), limits=c("poly-A", "RiboMinus"), labels=c("poly(A)", "RiboMinus")) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        text = element_text(size = 17), 
        legend.key.size = unit(0.4, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5, angle = 45), 
        axis.text.y = element_text(size=15), 
        legend.text = element_text(size=15)) +
  labs(x = "Cell line \n", y = "No. of circRNAs", title = "CIRI2 circRNA Counts")

p2 <- ggplot(hansen_circ, aes(x=sample, y=count, fill=method)) + 
  geom_bar(stat="identity", size = 0.5, width = 0.7, position = position_dodge(), color = "black") + 
  geom_signif(data = hansen_circ_annotations,
              aes(y_position = y, 
                  xmin = start,
                  xmax = end, 
                  annotations = labels),
              tip_length = c(0.05, 0.01),
              manual = TRUE) +
  geom_text(data = hansen_circ_labels,
            mapping = aes(x = c(0.9, 0.85, 1.9, 1.8, 2.9, 2.85),
                          y = c(9534, 5838, 13158, 10423, 14522, 11011),
                          label = label),
            hjust = -0.1, vjust = -1) +
  facet_grid(~validation) + 
  scale_y_continuous(limits = c(0, 15500), expand=c(0,0)) +
  scale_fill_manual("Method",values=c("#8B7B96", "#71A2B6"), limits=c("poly-A", "RiboMinus"), labels=c("poly(A)", "RiboMinus")) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        text = element_text(size = 17), 
        legend.key.size = unit(0.4, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5, angle = 45), 
        axis.text.y = element_text(size=15), 
        legend.text = element_text(size=15)) +
  labs(x = "Cell line \n", y = "No. of circRNAs", title = "CIRCexplorer2 circRNA Counts")


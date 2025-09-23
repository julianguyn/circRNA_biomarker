# Create analysis-ready count matrices for circRNA, isoforms, and gene expression #

# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggpubr)
    library(ggvenn)
    library(stringr)
    library(ComplexHeatmap)
})

options(stringsAsFactors = FALSE)


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
### circRNA Quantification Comparisons 
############################################################

# set palette for plotting
pal = c("#51C7AD", "#392C57", "#3670A0")

# remove sample names column for downstream analysis
ciri_gcsi <- ciri_gcsi[,-which(colnames(ciri_gcsi) %in% c("sample"))]
ciri_ccle <- ciri_ccle[,-which(colnames(ciri_ccle) %in% c("sample"))]
ciri_gdsc <- ciri_gdsc[,-which(colnames(ciri_gdsc) %in% c("sample"))]

circ_gcsi <- circ_gcsi[,-which(colnames(circ_gcsi) %in% c("sample"))]
circ_ccle <- circ_ccle[,-which(colnames(circ_ccle) %in% c("sample"))]
circ_gdsc <- circ_gdsc[,-which(colnames(circ_gdsc) %in% c("sample"))]

cfnd_gcsi <- cfnd_gcsi[,-which(colnames(cfnd_gcsi) %in% c("sample"))]
cfnd_ccle <- cfnd_ccle[,-which(colnames(cfnd_ccle) %in% c("sample"))]
cfnd_gdsc <- cfnd_gdsc[,-which(colnames(cfnd_gdsc) %in% c("sample"))]

fcrc_gcsi <- fcrc_gcsi[,-which(colnames(fcrc_gcsi) %in% c("sample"))]
fcrc_ccle <- fcrc_ccle[,-which(colnames(fcrc_ccle) %in% c("sample"))]
fcrc_gdsc <- fcrc_gdsc[,-which(colnames(fcrc_gdsc) %in% c("sample"))]

# ========== Unique Transcript Detection: Dataset Comparison per Pipeline ========== #

# create list object of transcripts
ciri_transcripts <- list(gCSI = colnames(ciri_gcsi), CCLE = colnames(ciri_ccle), GDSC2 = colnames(ciri_gdsc))
circ_transcripts <- list(gCSI = colnames(circ_gcsi), CCLE = colnames(circ_ccle), GDSC2 = colnames(circ_gdsc))
cfnd_transcripts <- list(gCSI = colnames(cfnd_gcsi), CCLE = colnames(cfnd_ccle), GDSC2 = colnames(cfnd_gdsc))
fcrc_transcripts <- list(gCSI = colnames(fcrc_gcsi), CCLE = colnames(fcrc_ccle), GDSC2 = colnames(fcrc_gdsc))

# plot venn diagram
p1 <- ggvenn(ciri_transcripts, 
        fill_color = pal, stroke_size = 0.5, set_name_size = 4) + 
        theme(plot.title = element_text(hjust = 0.5, size = 15)) + labs(title = "CIRI2")
p2 <- ggvenn(circ_transcripts, 
        fill_color = pal, stroke_size = 0.5, set_name_size = 4) + 
        theme(plot.title = element_text(hjust = 0.5, size = 15)) + labs(title = "CIRCexplorer2")
p3 <- ggvenn(cfnd_transcripts, 
        fill_color = pal, stroke_size = 0.5, set_name_size = 4) + 
        theme(plot.title = element_text(hjust = 0.5, size = 15)) + labs(title = "circRNA_finder")
p4 <- ggvenn(fcrc_transcripts, 
        fill_color = pal, stroke_size = 0.5, set_name_size = 4) + 
        theme(plot.title = element_text(hjust = 0.5, size = 15)) + labs(title = "find_circ")

#png("../results/figures/figure1/venndiagram_per_pipeline.png", width=500, height=150, units='mm', res = 600, pointsize=80)
png("../results/figures/figure1/common_venndiagram_per_pipeline.png", width=500, height=150, units='mm', res = 600, pointsize=80)
ggarrange(p1, p2, p3, p4, ncol = 4, nrow = 1, common.legend = FALSE)
dev.off()

# ========== Proportion plots ========== #

# function to count proportions
count_prop <- function(transcripts, label) {
    
    gcsi <- transcripts$gCSI
    ccle <- transcripts$CCLE
    gdsc <- transcripts$GDSC2

    # get all transcripts
    all_transcripts <- unique(c(gcsi, ccle, gdsc))

    # label transcript detection per dataset
    transcript_df <- data.frame(
        transcript = all_transcripts,
        gcsi = all_transcripts %in% gcsi,
        ccle = all_transcripts %in% ccle,
        gdsc = all_transcripts %in% gdsc
    )

    # label transcript detection
    transcript_df$category <- with(transcript_df, 
        ifelse(gcsi & !ccle & !gdsc, "gCSI only",
        ifelse(!gcsi & ccle & !gdsc, "CCLE only",
        ifelse(!gcsi & !ccle & gdsc, "GDSC2 only",
        ifelse(gcsi & ccle & !gdsc, "gCSI & CCLE",
        ifelse(gcsi & !ccle & gdsc, "gCSI & GDSC2",
        ifelse(!gcsi & ccle & gdsc, "CCLE & GDSC2",
        "All datasets"))))))
    )

    # create dataframe for plotting
    toPlot <- table(transcript_df$category) |> as.data.frame()
    toPlot$Var1 <- factor(toPlot$Var1, levels = c("All datasets",
                                                  "gCSI & CCLE",
                                                  "gCSI & GDSC2",
                                                  "CCLE & GDSC2",
                                                  "gCSI only",
                                                  "CCLE only",
                                                  "GDSC2 only"))
    toPlot$Prop <- round(toPlot$Freq / sum(toPlot$Freq) * 100, digits = 2)
    toPlot$label <- paste0(toPlot$Freq, " (", toPlot$Prop, "%)")
    toPlot$pipeline <- label

    return(toPlot)
}

toPlot <- rbind(count_prop(ciri_transcripts, "CIRI2"),
                count_prop(circ_transcripts, "CIRCexplorer2"),
                count_prop(cfnd_transcripts, "circRNA_finder"),
                count_prop(fcrc_transcripts, "find_circ"))
toPlot$pipeline <- factor(toPlot$pipeline, 
                          levels = c("CIRI2", "CIRCexplorer2", 
                                     "circRNA_finder", "find_circ"))

pal = c("#987A82", "#DACCAB", "#C78B76", "#9D3737", "#51C7AD", "#392C57", "#3670A0")

# plot
png("../results/figures/figure1/proportion_pipelines.png", width=200, height=150, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(fill = Var1, y = Freq, x = pipeline)) + 
  geom_bar(position = "fill", stat = "identity", color = "black") +
  geom_text(aes(label = ifelse(Var1 %in% c("gCSI only",
                                           "CCLE only",
                                           "GDSC2 only"), label, "")), 
            position = position_fill(vjust = 0.5)) +
  theme_classic() + theme(legend.key.size = unit(0.5, 'cm')) +
  scale_fill_manual(values = pal) +
  labs(fill = "Category", x = "Pipeline", y = "Proportion of Unique Transcripts")
dev.off()

# ========== Unique Transcript Detection: Pipeline Comparison ========== #

# create list object of transcripts for venn diagram
all_comparison <- list(CIRI2 = c(colnames(ciri_gcsi), colnames(ciri_ccle), colnames(ciri_gdsc)), 
                       CIRCexplorer2 = c(colnames(circ_gcsi), colnames(circ_ccle), colnames(circ_gdsc)), 
                       circRNA_finder = c(colnames(cfnd_gcsi), colnames(cfnd_ccle), colnames(cfnd_gdsc)),
                       find_circ = c(colnames(fcrc_gcsi), colnames(fcrc_ccle), colnames(fcrc_gdsc)))

# plot venn diagram
#png("../results/figures/figure1/venndiagram.png", width=200, height=150, units='mm', res = 600, pointsize=80)
png("../results/figures/figure1/common_venndiagram.png", width=200, height=150, units='mm', res = 600, pointsize=80)
ggvenn(all_comparison, 
        fill_color = c("#E8F6B1", "#A5DBB7", "#2088BC", "#26479D"),
        stroke_size = 0.5, set_name_size = 4)
dev.off()


# create list object of transcripts for upset plot
set.seed(123)
toPlot <- make_comb_mat(list(
            gCSI_CIRI2 = colnames(ciri_gcsi), CCLE_CIRI2 = colnames(ciri_ccle), GDSC_CIRI2 = colnames(ciri_gdsc),
            gCSI_CIRCexplorer2 = colnames(circ_gcsi), CCLE_CIRCexplorer2 = colnames(circ_ccle), GDSC_CIRCexplorer2 = colnames(circ_gdsc),
            gCSI_circRNA_finder = colnames(cfnd_gcsi), CCLE_circRNA_finder = colnames(cfnd_ccle), GDSC_circRNA_finder = colnames(cfnd_gdsc),
            gCSI_find_circ = colnames(fcrc_gcsi), CCLE_find_circ = colnames(fcrc_ccle), GDSC_find_circ = colnames(fcrc_gdsc)))

# remove combinations of less than 100 pairs
toPlot <- toPlot[comb_size(toPlot) >= 100]

# upset plot
#pdf("../results/figures/figure1/upset_by_pipeline.pdf", width=10, height=5)
pdf("../results/figures/figure1/common_upset_by_pipeline.pdf", width=10, height=5)
UpSet(toPlot, set_order = c("gCSI_CIRI2", "CCLE_CIRI2", "GDSC_CIRI2", "gCSI_CIRCexplorer2", "CCLE_CIRCexplorer2", "GDSC_CIRCexplorer2",
                            "gCSI_circRNA_finder", "CCLE_circRNA_finder", "GDSC_circRNA_finder", "gCSI_find_circ", "CCLE_find_circ", "GDSC_find_circ"),
        top_annotation = upset_top_annotation(toPlot, add_numbers = TRUE),
        comb_order = order(comb_size(toPlot)))
dev.off()

#pdf("../results/figures/figure1/upset_by_pset.pdf", width=10, height=5)
#UpSet(toPlot, set_order = c("gCSI_CIRI2", "gCSI_CIRCexplorer2", "gCSI_circRNA_finder", "gCSI_find_circ",
#                            "CCLE_CIRCexplorer2", "CCLE_CIRI2", "CCLE_circRNA_finder", "CCLE_find_circ",
#                            "GDSC_CIRI2", "GDSC_CIRCexplorer2", "GDSC_circRNA_finder", "GDSC_find_circ"),
#        top_annotation = upset_top_annotation(toPlot, add_numbers = TRUE),
#        comb_order = order(comb_size(toPlot)))
#dev.off()


# ========== Transcript Quantification: Dataset Comparison per Pipeline ========== #

# TODO: remove
ciri_gdsc[is.na(ciri_gdsc)] <- 0
circ_gdsc[is.na(circ_gdsc)] <- 0
cfnd_gdsc[is.na(cfnd_gdsc)] <- 0
fcrc_gdsc[is.na(fcrc_gdsc)] <- 0

# create data frame of counts for plotting
df <- data.frame(Count = c(sum(ciri_gcsi), sum(ciri_ccle), sum(ciri_gdsc),
                           sum(circ_gcsi), sum(circ_ccle), sum(circ_gdsc), 
                           sum(cfnd_gcsi), sum(cfnd_ccle), sum(cfnd_gdsc), 
                           sum(fcrc_gcsi), sum(fcrc_ccle), sum(fcrc_gdsc)),
                PSet = c(rep(c("gCSI", "CCLE", "GDSC2"), 4)),
                Pipeline = c(rep("CIRI2", 3), rep("CIRCexplorer2", 3), rep("circRNA_finder", 3), rep("find_circ", 3)))
df$Pipeline <- factor(df$Pipeline, levels = c("CIRI2", "CIRCexplorer2", "circRNA_finder", "find_circ"))
df$PSet <- factor(df$PSet, levels = c("gCSI", "CCLE", "GDSC2"))

# plot bar plot of counts
#png("../results/figures/figure1/counts.png", width=200, height=150, units='mm', res = 600, pointsize=80)
#ggplot(df, aes(x = Pipeline, y = Count, fill = PSet)) + geom_bar(stat="identity", position = "dodge", color = "black") +
#  scale_fill_manual(values=pal, limits=c("gCSI", "CCLE", "GDSC2")) + 
#  scale_y_continuous(limits = c(0, 240000), expand=c(0,0))  + theme_classic() + 
#  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
#        legend.key.size = unit(0.4, 'cm'))
#dev.off()

png("../results/figures/figure1/common_counts.png", width=150, height=100, units='mm', res = 600, pointsize=80)
ggplot(df, aes(x = Pipeline, y = log2(Count), fill = PSet)) + geom_bar(stat="identity", position = "dodge", color = "black") +
  scale_fill_manual(values=pal, limits=c("gCSI", "CCLE", "GDSC2")) + 
  scale_y_continuous(limits = c(0, 13), expand=c(0,0))  + theme_classic() + 
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        legend.key.size = unit(0.4, 'cm')) + labs(y = "Log2 Normalized Counts")
dev.off()


# following code is not used in figures:
# ========== Distribution of Transcript Expression ========== #

suppressWarnings(ciri_gcsi <- fread("../data/processed_cellline/all_samples/CIRI2/ciri_gcsi_counts.tsv", data.table = F))
suppressWarnings(ciri_gdsc <- fread("../data/processed_cellline/all_samples/CIRI2/ciri_gdsc_counts.tsv", data.table = F))
suppressWarnings(ciri_ccle <- fread("../data/processed_cellline/all_samples/CIRI2/ciri_ccle_counts.tsv", data.table = F))

suppressWarnings(circ_gcsi <- fread("../data/processed_cellline/all_samples/CIRCexplorer2/circ_gcsi_counts.tsv", data.table = F))
suppressWarnings(circ_gdsc <- fread("../data/processed_cellline/all_samples/CIRCexplorer2/circ_gdsc_counts.tsv", data.table = F))
suppressWarnings(circ_ccle <- fread("../data/processed_cellline/all_samples/CIRCexplorer2/circ_ccle_counts.tsv", data.table = F))

suppressWarnings(cfnd_gcsi <- fread("../data/processed_cellline/all_samples/circRNA_finder/cfnd_gcsi_counts.tsv", data.table = F))
suppressWarnings(cfnd_gdsc <- fread("../data/processed_cellline/all_samples/circRNA_finder/cfnd_gdsc_counts.tsv", data.table = F))
suppressWarnings(cfnd_ccle <- fread("../data/processed_cellline/all_samples/circRNA_finder/cfnd_ccle_counts.tsv", data.table = F))

suppressWarnings(fcrc_gcsi <- fread("../data/processed_cellline/all_samples/find_circ/fcrc_gcsi_counts.tsv", data.table = F))
suppressWarnings(fcrc_gdsc <- fread("../data/processed_cellline/all_samples/find_circ/fcrc_gdsc_counts.tsv", data.table = F))
suppressWarnings(fcrc_ccle <- fread("../data/processed_cellline/all_samples/find_circ/fcrc_ccle_counts.tsv", data.table = F))


# remove sample names column for downstream analysis
ciri_gcsi <- ciri_gcsi[,-which(colnames(ciri_gcsi) %in% c("sample"))]
ciri_ccle <- ciri_ccle[,-which(colnames(ciri_ccle) %in% c("sample"))]
ciri_gdsc <- ciri_gdsc[,-which(colnames(ciri_gdsc) %in% c("sample"))]

circ_gcsi <- circ_gcsi[,-which(colnames(circ_gcsi) %in% c("sample"))]
circ_ccle <- circ_ccle[,-which(colnames(circ_ccle) %in% c("sample"))]
circ_gdsc <- circ_gdsc[,-which(colnames(circ_gdsc) %in% c("sample"))]

cfnd_gcsi <- cfnd_gcsi[,-which(colnames(cfnd_gcsi) %in% c("sample"))]
cfnd_ccle <- cfnd_ccle[,-which(colnames(cfnd_ccle) %in% c("sample"))]
cfnd_gdsc <- cfnd_gdsc[,-which(colnames(cfnd_gdsc) %in% c("sample"))]

fcrc_gcsi <- fcrc_gcsi[,-which(colnames(fcrc_gcsi) %in% c("sample"))]
fcrc_ccle <- fcrc_ccle[,-which(colnames(fcrc_ccle) %in% c("sample"))]
fcrc_gdsc <- fcrc_gdsc[,-which(colnames(fcrc_gdsc) %in% c("sample"))]


# create data frame of expression values
print("making expression")
Expression = c(as.vector(unlist(ciri_gcsi)), as.vector(unlist(ciri_ccle)), as.vector(unlist(ciri_gdsc)),
                                    as.vector(unlist(circ_gcsi)), as.vector(unlist(circ_ccle)), as.vector(unlist(circ_gdsc)),
                                    as.vector(unlist(cfnd_gcsi)), as.vector(unlist(cfnd_ccle)), as.vector(unlist(cfnd_gdsc)),
                                    as.vector(unlist(fcrc_gcsi)), as.vector(unlist(fcrc_ccle)), as.vector(unlist(fcrc_gdsc)))
print("making label")
Label = c(rep("gCSI-CIRI2", prod(dim(ciri_gcsi))), rep("CCLE-CIRI2", prod(dim(ciri_ccle))), rep("GDSC2-CIRI2", prod(dim(ciri_gdsc))),
                               rep("gCSI-CIRCexplorer2", prod(dim(circ_gcsi))), rep("CCLE-CIRCexplorer2", prod(dim(circ_ccle))), rep("GDSC2-CIRCexplorer2", prod(dim(circ_gdsc))),
                               rep("gCSI-circRNA_finder", prod(dim(cfnd_gcsi))), rep("CCLE-circRNA_finder", prod(dim(cfnd_ccle))), rep("GDSC2-circRNA_finder", prod(dim(cfnd_gdsc))),
                               rep("gCSI-find_circ", prod(dim(fcrc_gcsi))), rep("CCLE-find_circ", prod(dim(fcrc_ccle))), rep("GDSC2-find_circ", prod(dim(fcrc_gdsc))))
                
toPlot <- data.frame(Expression, Label)
dim(toPlot)
toPlot$Pipeline <- gsub(".*-", "", toPlot$Label)
toPlot$PSet <- gsub("-.*", "", toPlot$Label)
toPlot$Pipeline <- factor(toPlot$Pipeline, levels = c("CIRI2", "CIRCexplorer2", "circRNA_finder", "find_circ"))
toPlot$PSet <- factor(toPlot$PSet, levels = c("gCSI", "CCLE", "GDSC2"))


max = max(toPlot$Expression)

png("../results/figures/figure1/count_dist_all.png", width=150, height=100, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = Expression)) + geom_density(aes(color = Pipeline), alpha = 0.5, size = 1.5) + 
        theme_classic() + xlim(0, max+5) + 
        labs(x = "Normalized circRNA Count", y = "Density") +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
            legend.key.size = unit(0.4, 'cm'))
dev.off()

png("../results/figures/figure1/count_dist_pipeline.png", width=250, height=100, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = Expression)) + geom_density(aes(color = Pipeline), alpha = 0.5, size = 1.5) + 
        theme_classic() + xlim(0, max+5) + facet_grid(.~Pipeline) +
        labs(x = "Normalized circRNA Count", y = "Density") +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
            legend.key.size = unit(0.4, 'cm'))
dev.off()

png("../results/figures/figure1/count_dist_pset.png", width=250, height=100, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = Expression)) + geom_density(aes(color = PSet), alpha = 0.5, size = 1.5) + 
        theme_classic() + xlim(0, max+5) + facet_grid(.~PSet) +
        labs(x = "Normalized circRNA Count", y = "Density") +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
            legend.key.size = unit(0.4, 'cm'))
dev.off()

# COMMON:
png("../results/figures/figure1/common_count_dist_all.png", width=150, height=100, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = Expression)) + geom_density(aes(color = Pipeline), alpha = 0.5, size = 1.5) + 
        theme_classic() + xlim(0, max+5) + 
        labs(x = "Normalized circRNA Count", y = "Density") +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
            legend.key.size = unit(0.4, 'cm'))
dev.off()

png("../results/figures/figure1/common_count_dist_pipeline.png", width=250, height=100, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = Expression)) + geom_density(aes(color = Pipeline), alpha = 0.5, size = 1.5) + 
        theme_classic() + xlim(0, max+5) + facet_grid(.~Pipeline) +
        labs(x = "Normalized circRNA Count", y = "Density") +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
            legend.key.size = unit(0.4, 'cm'))
dev.off()

png("../results/figures/figure1/common_count_dist_pset.png", width=250, height=100, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = Expression)) + geom_density(aes(color = PSet), alpha = 0.5, size = 1.5) + 
        theme_classic() + xlim(0, max+5) + facet_grid(.~PSet) +
        labs(x = "Normalized circRNA Count", y = "Density") +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
            legend.key.size = unit(0.4, 'cm'))
dev.off()

# compute distribution of 0 values in normalized circRNA counts
library(dplyr)
toPlot$Group <- toPlot$Expression == 0
toPlot_counts <- toPlot %>% group_by(Pipeline, Group) %>% summarise(Count = n()) 
toPlot_counts

toPlot <- toPlot[toPlot$Expression != 0,]
save(toPlot, toPlot_counts, file = "temp.RData")


toPlot_counts <- as.data.frame(toPlot_counts)
# pie chart
p1 <- ggplot(toPlot_counts[toPlot_counts$Pipeline == "CIRI2",], aes(x = "", y = Count, fill = Group)) +
    geom_bar(stat="identity", width = 1, color = "black") +
    theme(legend.key.size = unit(0.4, 'cm'), plot.title = element_text(hjust = 0.5, size = 18)) +
    geom_text(aes(label = Count), position = position_stack(vjust = 0.5)) +
    scale_fill_manual(values = c("#839788", "gray")) +
    coord_polar("y", start = 0) + theme_void() + labs(fill = "Expression = 0") + ggtitle("CIRI2")
p2 <- ggplot(toPlot_counts[toPlot_counts$Pipeline == "CIRCexplorer2",], aes(x = "", y = Count, fill = Group)) +
    geom_bar(stat="identity", width = 1, color = "black") +
    theme(legend.key.size = unit(0.4, 'cm'), plot.title = element_text(hjust = 0.5, size = 18)) +
    geom_text(aes(label = Count), position = position_stack(vjust = 0.5)) +
    scale_fill_manual(values = c("#839788", "gray")) +
    coord_polar("y", start = 0) + theme_void() + labs(fill = "Expression = 0") + ggtitle("CIRCexplorer2")
p3 <- ggplot(toPlot_counts[toPlot_counts$Pipeline == "circRNA_finder",], aes(x = "", y = Count, fill = Group)) +
    geom_bar(stat="identity", width = 1, color = "black") +
    theme(legend.key.size = unit(0.4, 'cm'), plot.title = element_text(hjust = 0.5, size = 18)) +
    geom_text(aes(label = Count), position = position_stack(vjust = 0.5)) +
    scale_fill_manual(values = c("#839788", "gray")) +
    coord_polar("y", start = 0) + theme_void() + labs(fill = "Expression = 0") + ggtitle("circRNA_finder")
p4 <- ggplot(toPlot_counts[toPlot_counts$Pipeline == "find_circ",], aes(x = "", y = Count, fill = Group)) +
    geom_bar(stat="identity", width = 1, color = "black") +
    theme(legend.key.size = unit(0.4, 'cm'), plot.title = element_text(hjust = 0.5, size = 18)) +
    geom_text(aes(label = Count), position = position_stack(vjust = 0.5)) +
    scale_fill_manual(values = c("#839788", "gray")) +
    coord_polar("y", start = 0) + theme_void() + labs(fill = "Expression = 0") + ggtitle("find_circ")
png("count_zero_pipeline.png", width=100, height=100, units='mm', res = 600, pointsize=80)
ggarrange(p1, p2, p3, p4,
          ncol = 2, nrow = 2,
          common.legend = TRUE,
          legend = "right")
dev.off()

png("../results/figures/figure1/common_count_zero_pipeline.png", width=250, height=100, units='mm', res = 600, pointsize=80)
ggplot(toPlot_counts, aes(x = "", y = Count, fill = Label)) +
    geom_bar(stat="identity", width = 1, color = "black") +
    facet_grid(.~Pipeline) + theme(legend.key.size = unit(0.4, 'cm')) +
    geom_text(aes(label = Count), position = position_stack(vjust = 0.5)) +
    coord_polar("y", start = 0) + theme_void() + labs(fill = "Tissue")
dev.off()


toPlot_counts <- toPlot %>% group_by(PSet, Group) %>% summarise(Count = n()) 
toPlot_counts

# pie chart
png("../results/figures/figure1/count_zero_pset.png", width=250, height=100, units='mm', res = 600, pointsize=80)
ggplot(toPlot_counts, aes(x = "", y = Count, fill = Label)) +
    geom_bar(stat="identity", width = 1, color = "black") +
    facet_grid(.~PSet) + theme(legend.key.size = unit(0.4, 'cm')) +
    geom_text(aes(label = Count), position = position_stack(vjust = 0.5)) +
    coord_polar("y", start = 0) + theme_void() + labs(fill = "Tissue")
dev.off()

png("../results/figures/figure1/common_count_zero_pset.png", width=250, height=100, units='mm', res = 600, pointsize=80)
ggplot(toPlot_counts, aes(x = "", y = Count, fill = Label)) +
    geom_bar(stat="identity", width = 1, color = "black") +
    facet_grid(.~PSet) + theme(legend.key.size = unit(0.4, 'cm')) +
    geom_text(aes(label = Count), position = position_stack(vjust = 0.5)) +
    coord_polar("y", start = 0) + theme_void() + labs(fill = "Tissue")
dev.off()


# remove all 0s
toPlot <- toPlot[toPlot$Expression != 0,]

max = max(toPlot$Expression)

png("count_dist_common_nozero.png", width=200, height=150, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = Pipeline, y = log2(Expression))) + 
    geom_violin(aes(fill = Pipeline), alpha = 0.8) + geom_boxplot(width=0.1, alpha = 0.3) +
    theme_classic() + labs(x = "", fill = "", y = "log2 Normalized circRNA Count Value") +
    scale_fill_manual(values = c("#839788", "#BFD7EA", "#BA9790", "#D5BC8A")) +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3), legend.key.size = unit(0.7, 'cm')) +
    coord_flip()
dev.off()

ggplot(toPlot, aes(x = log2(Expression))) + geom_density(aes(fill = Pipeline), alpha = 0.5, size = 1) + 
        theme_classic() + facet_grid(factor(Pipeline)~.) +
        scale_fill_manual(values = c("#839788", "#BFD7EA", "#BA9790", "#D5BC8A")) +
        labs(x = "log2 Normalized circRNA Count", y = "Density") +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
            legend.key.size = unit(0.4, 'cm'))

png("count_dist_common_nozero.png", width=200, height=150, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = log2(Expression))) + geom_density(aes(fill = Pipeline), alpha = 0.5) + 
        theme_classic() + 
        scale_fill_manual(values = c("#839788", "#BFD7EA", "#BA9790", "#D5BC8A")) +
        labs(x = "log2 Normalized circRNA Count", y = "Density") +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
            legend.key.size = unit(0.4, 'cm'))
dev.off()

png("../results/figures/figure1/count_dist_all_nozero.png", width=150, height=100, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = log2(Expression))) + geom_density(aes(color = Pipeline), alpha = 0.5, size = 1) + 
        theme_classic() + xlim(0, max+5) + 
        labs(x = "log2 Normalized circRNA Count", y = "Density") +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
            legend.key.size = unit(0.4, 'cm'))
dev.off()

png("../results/figures/figure1/count_dist_pipeline_nozero.png", width=250, height=100, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = Expression)) + geom_density(aes(color = Pipeline), alpha = 0.5, size = 1.5) + 
        theme_classic() + xlim(0, max+5) + facet_grid(.~Pipeline) +
        labs(x = "Normalized circRNA Count", y = "Density") +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
            legend.key.size = unit(0.4, 'cm'))
dev.off()

png("../results/figures/figure1/count_dist_pset_nozero.png", width=250, height=100, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = Expression)) + geom_density(aes(color = PSet), alpha = 0.5, size = 1.5) + 
        theme_classic() + xlim(0, max+5) + facet_grid(.~PSet) +
        labs(x = "Normalized circRNA Count", y = "Density") +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
            legend.key.size = unit(0.4, 'cm'))
dev.off()


#COMMON:
png("../results/figures/figure1/common_count_dist_all_nozero.png", width=150, height=100, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = Expression)) + geom_density(aes(color = Pipeline), alpha = 0.5, size = 1.5) + 
        theme_classic() + xlim(0, max+5) + 
        labs(x = "Normalized circRNA Count", y = "Density") +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
            legend.key.size = unit(0.4, 'cm'))
dev.off()

png("../results/figures/figure1/common_count_dist_pipeline_nozero.png", width=250, height=100, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = Expression)) + geom_density(aes(color = Pipeline), alpha = 0.5, size = 1.5) + 
        theme_classic() + xlim(0, max+5) + facet_grid(.~Pipeline) +
        labs(x = "Normalized circRNA Count", y = "Density") +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
            legend.key.size = unit(0.4, 'cm'))
dev.off()

png("../results/figures/figure1/common_count_dist_pset_nozero.png", width=250, height=100, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = Expression)) + geom_density(aes(color = PSet), alpha = 0.5, size = 1.5) + 
        theme_classic() + xlim(0, max+5) + facet_grid(.~PSet) +
        labs(x = "Normalized circRNA Count", y = "Density") +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
            legend.key.size = unit(0.4, 'cm'))
dev.off()


##TODO: Add code for lung, isoform, and gene expression processing, normalization